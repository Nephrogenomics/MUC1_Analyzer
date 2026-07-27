#!/usr/bin/env python3
"""dupc_dispatch — auto-route a MUC1 dupC / frameshift call by input regime (locus depth).

Two regimes need two DIFFERENT callers, and running the wrong one is a documented failure mode either way:
  - DEEP PCR / amplicon (thousands of reads spanning the VNTR) → the **per-allele** caller (`pcr_dupc` /
    `pcr_variant`). The pooled statistical test over-calls here: the PCR-depleted mutant dilutes to the ~5 %
    ONT homopolymer floor, which is "significant" (p≈1e-45) in NEGATIVES too at amplicon depth.
  - WGS / adaptive-sampling / urine (~10–40×) → the **statistical** caller (`dupc_caller`, which auto-routes
    the anchor: genomic for hg38-phased, sequence for unphased/T2T/urine). The per-allele effect-size call has
    no power at that depth; the PoN/posterior model is what detects there. ⚠ WGS/AS dupC is coverage-limited
    (~half the carriers at full AS depth; collapses below) — this dispatch picks the RIGHT caller, it does not
    beat the coverage floor.

Routing is by the read count spanning the locus — a clean separator (AS/WGS ≈ 10–115 reads vs PCR thousands;
our AS CRAMs max ~115, our PCR amplicons min ~458).

    python3 -m muc1_analyzer.dupc_dispatch -b sample.bam --region chr1:155188000-155192000 [--ref REF] \
        [--variant insG] [--pon pon_dupc.json] [--amplicon-min 200]

⚠ pysam (compute). `pick_regime` is pure + unit-tested.
"""
from __future__ import annotations
import argparse
import json
import sys
from typing import Optional


def pick_regime(n_reads: int, *, amplicon_min: int = 200) -> str:
    """'pcr' (deep amplicon → per-allele caller) vs 'wgs' (low-cov → statistical caller). Pure."""
    return "pcr" if n_reads >= amplicon_min else "wgs"


def _count_locus_reads(bam: str, chrom: str, start: int, end: int, ref: Optional[str] = None) -> int:
    import pysam
    mode = "rc" if str(bam).endswith(".cram") else "rb"
    kw = {"reference_filename": ref} if (mode == "rc" and ref) else {}
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        return sum(1 for r in af.fetch(chrom, start, end)
                   if not (r.is_secondary or r.is_supplementary or r.is_unmapped))


def dispatch_dupc(bam: str, *, chrom: str, start: int, end: int, ref: Optional[str] = None,
                  amplicon_min: int = 200, mut_ctract: int = 8, variant: Optional[str] = None,
                  pon=None, fast: bool = False) -> dict:
    """Count locus reads → pick the regime → run the matching caller. Returns the routed result + the regime.
    `fast` (WGS/AS only): use the positional path (`scan_dupc_positional`) instead of the context+PoN scan —
    ~1500× faster (~0.2 s vs ~5-7 min/sample) with IDENTICAL specificity (0 FP on 233 MAT+1000G negatives).
    Recommended at scale; the context+PoN scan only adds marginal per-context-null sensitivity."""
    n = _count_locus_reads(bam, chrom, start, end, ref)
    regime = pick_regime(n, amplicon_min=amplicon_min)
    if regime == "pcr":
        if variant:
            from .pcr_variant import pcr_variant
            res = pcr_variant(bam, chrom=chrom, start=start, end=end, ref=ref, variant=variant)
            caller = "pcr_variant"
        else:
            from .pcr_dupc import pcr_dupc
            res = pcr_dupc(bam, chrom=chrom, start=start, end=end, ref=ref, mut_ctract=mut_ctract)
            caller = "pcr_dupc"
    elif fast:
        from .detectors.vntr_dupc import scan_dupc_positional
        d = scan_dupc_positional(bam, genome_ref=ref, mut_len=mut_ctract)   # genomic positional, PoN-free
        res = {"called": bool(d.get("called")), "tier": "positional", "best": d.get("best"),
               "n_span": d.get("n_span")}
        caller = "positional(fast)"
    else:
        from .dupc_caller import call_dupc
        pon_d = json.load(open(pon)) if (pon and isinstance(pon, str)) else pon
        # Let call_dupc AUTO-ROUTE the anchor: genomic (hg38 phased, e.g. a homozygous-length sample) vs
        # sequence (unphased / T2T / urine, keyed on HP presence). Do NOT force seq_anchor — that pushes an
        # hg38 CRAM to the seq path (wrong + errors). Pass the region so the fetch is LOCUS-scoped (region=None
        # reads to EOF → "truncated file" on some CRAMs).
        probe = "unitlen" if variant in ("del8_27",) else "ctract"    # del8_27 = body-deletion unit-length probe
        res = call_dupc(bam, pon=pon_d, genome_ref=ref, mut_len=mut_ctract,
                        region=f"{chrom}:{start}-{end}", probe=probe)
        caller = f"dupc_caller(auto-anchor,{probe})"
    return {"regime": regime, "caller": caller, "n_locus_reads": n, "amplicon_min": amplicon_min,
            "region": f"{chrom}:{start}-{end}", "result": res}


def dispatch_dupc_vntr(vntr_bam: str, vntr_ref: str, *, lengths=None, mut_len: int = 8,
                       wt_len: int = 7) -> dict:
    """dupC verdict IN VNTR SPACE — the route for an input with no genomic locus (fastq / uBAM).

    Every other route reads the ORIGINAL genomic reads at the hg38 dupC locus, so a fastq got no reliable
    verdict at all: the pipeline fell back to the consensus token, which is documented depth-fragile. The
    run-length-shift caller needs no genomic anchor and no PoN — it calibrates on the sample's own C-runs —
    so it closes exactly that hole.

    Scaffolded PER ALLELE (from the arbiter's called `lengths` when available): a heterozygote's long
    allele is where the dupC usually sits, and a sample scaffolded on its modal contig never looks at it.
    Validated before wiring: specificity 1.000 on 92 genomic and on 15 tandem non-carriers, sensitivity 2/2.

    Returns the `dispatch_dupc` shape, plus `power_by_allele` so the report can say whether a NEGATIVE is
    adequately powered rather than printing a bare "not called"."""
    from .allele_scaffold import allele_contigs
    from .detectors import runlen_shift as R
    from .dupc_power import min_reads as _min_reads

    contigs = allele_contigs(vntr_bam, lengths=lengths)
    per_allele, best, power = {}, None, {}
    for contig in contigs:
        try:
            r = R.call_bam(vntr_bam, contig, vntr_ref, mut_len=mut_len, wt_len=wt_len)
        except Exception as e:
            per_allele[contig] = {"error": f"{type(e).__name__}: {e}"[:80]}
            continue
        per_allele[contig] = r
        cand = r.get("candidate") or {}
        n_cover = cand.get("n") or 0
        need = _min_reads(0.3, target=0.90) or 12      # reads to reach 90 % power at d=0.3
        power[contig] = {"n_cover": n_cover, "adequately_powered": bool(n_cover >= need),
                         "reads_for_power0.90": {"d=0.3": need}}
        if r.get("called") and (best is None or cand.get("p_bonf", 1) < best[1].get("p_bonf", 1)):
            best = (contig, cand)

    if best:
        contig, cand = best
        interp, called = "CONFIRMED", True
        variant = {"repeat": cand.get("position"), "label": "59dupC" if mut_len == 8 else f"mut_len{mut_len}",
                   "k_ge_mut": cand.get("k_ge_mut"), "n": cand.get("n"), "p_bonf": cand.get("p_bonf")}
    else:
        contig, variant, called = None, None, False
        tested = [c for c, r in per_allele.items() if (r.get("n_tested") or 0)]
        if not tested:
            interp = "INDETERMINATE"          # nothing measurable — never a clean negative
        elif all(power.get(c, {}).get("adequately_powered") for c in tested):
            interp = "NEGATIVE"
        else:
            interp = "NEGATIVE_BORDERLINE" if any(
                power.get(c, {}).get("n_cover", 0) >= 6 for c in tested) else "INDETERMINATE"

    return {"regime": "vntr", "caller": "runlen_shift(per-allele)", "region": None,
            "scaffolds": contigs,
            "result": {"assessed": bool(contigs), "called": called, "interpretation": interp,
                       "method": "runlen_shift", "carrier_contig": contig, "variant": variant,
                       "power_by_allele": power, "per_allele": per_allele,
                       "reason": None if contigs else "no scaffold contig could be chosen"}}


def _parse_region(s: str):
    chrom, rng = s.split(":")
    a, b = rng.replace(",", "").split("-")
    return chrom, int(a), int(b)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="muc1_analyzer.dupc_dispatch",
                                 description="Auto-route a MUC1 dupC/frameshift call by input regime (depth)")
    ap.add_argument("-b", "--bam", required=True)
    ap.add_argument("--region", required=True, help="e.g. chr1:155188000-155192000")
    ap.add_argument("--ref", default=None, help="reference FASTA (required for a CRAM)")
    ap.add_argument("--variant", default=None, help="PCR regime: call this dictionary variant (insG/dupG/…) via pcr_variant")
    ap.add_argument("--pon", default=None, help="WGS regime: pon_dupc.json for the statistical caller")
    ap.add_argument("--mut-ctract", type=int, default=8, help="8 = dupC (default), 5 = delCC")
    ap.add_argument("--amplicon-min", type=int, default=200, help="reads spanning the locus ≥ this ⇒ PCR regime")
    ap.add_argument("--fast", action="store_true",
                    help="RECOMMENDED for WGS/AS: positional path — ~1500× faster (~0.2 s vs ~5-7 min/sample) "
                         "with IDENTICAL specificity (0 FP on 233 MAT+1000G). Without it you get the ~7-min "
                         "context+PoN scan (only a marginal per-context-null sensitivity gain).")
    a = ap.parse_args(argv)
    chrom, start, end = _parse_region(a.region)
    out = dispatch_dupc(a.bam, chrom=chrom, start=start, end=end, ref=a.ref, amplicon_min=a.amplicon_min,
                        mut_ctract=a.mut_ctract, variant=a.variant, pon=a.pon, fast=a.fast)
    if out["regime"] == "wgs" and not a.fast:
        sys.stderr.write("[dupc_dispatch] WGS/AS context+PoN scan (~5-7 min). Add --fast for the ~0.2 s "
                         "positional path (same 0-FP specificity).\n")
    print(json.dumps(out, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
