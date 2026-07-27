#!/usr/bin/env python3
"""pcr_dupc — per-allele dupC caller for HIGH-DEPTH PCR / amplicon MUC1 data.

Why a separate mode (validated 2026-07-14 on the Madritsch ENA PRJEB92208 amplicons): on deep amplicons
the pooled statistical caller (`dupc_caller`) FAILS. PCR depletes the long (often mutant) allele, so its
dupC dilutes to the ONT 7C→8C homopolymer error floor (~5%); and at thousands of reads that floor is
"significant" (p ~1e-45) in NEGATIVES too → a p-value is useless here. The fix, validated: MP1 long-allele
dupC fraction jumps 0.073 (pool) → 0.275 (long allele only) while HG001 stays 0.053 (floor).

So the PCR mode:
  1. **separates reads by array length** (per-read unit count → 2 alleles), which also de-biases the PCR
     amplification skew (the long allele is under-amplified but still present);
  2. **scans the dupC PER allele** (de-diluted from the pool);
  3. **calls by EFFECT SIZE** above the homopolymer floor (not by p-value), gated by a **coverage/focality**
     check — the mutant index must be covered by a substantial share of the allele's reads (rejects the
     array-EDGE small-number 8C spikes that false-positived the HG001 reference negative);
  4. **tests BOTH alleles** — the mutant may sit on the short OR the long one (Madritsch: MP1/2/4 long, MP3 short);
  5. requires **per-strand support** to reject strand-dependent homopolymer artifacts;
  6. reports a FOUR-way status — positive / negative / **undetermined** / low_coverage — so a sample whose
     internal null is corrupted (BON: the WT/control allele's own 8C rate is itself at calling level, from a
     strand-asymmetric homopolymer artifact) surfaces as UNDETERMINED, not a silent False. BON's *mutant*
     allele is indistinguishable from the HG002 negative (both ~0.10, rev-only, tot ~20-37) → the honest
     discriminator is the CONTROL allele's reliability, not the candidate's strand.

The low-coverage WGS/urine statistical mode (`dupc_caller`) is UNCHANGED — this is the deep-amplicon regime
only (dispatch by context). Reuses the C-tract machinery (`_CTX` / `gene_oriented` / `best_8c_index`).
Length here comes from the read unit-count (NOT the false-SHORT contig method, NOT the flank arbiter which
needs flanks amplicons lack).

⚠ pysam (compute). The cores (`split_by_length`, `call_by_effect`, `call_from_reads`) are pysam-free + tested.

  python3 -m muc1_analyzer.pcr_dupc -b sample.muc1win.bam --region muc1win:2328605-2330538   # dupC (8C)
  python3 -m muc1_analyzer.pcr_dupc -b sample.muc1win.bam --region … --mut-ctract 5           # delCC (5C)
"""
from __future__ import annotations
import argparse
import collections
import json
import statistics as st
import sys
from typing import Optional

from .allele_balance import best_8c_index


# ── Pure cores (pysam-free, unit-tested) ──────────────────────────────────────

def split_by_length(nunits, *, min_sep: int = 5, min_reads: int = 6) -> dict:
    """1D 2-means on per-read unit counts → allele split threshold. `het=False` ⇒ one length allele
    (homozygous by length; can't de-dilute). Returns {het, threshold, medians:[...], n}."""
    vals = sorted(int(x) for x in nunits)
    n = len(vals)
    if n < min_reads:
        return {"het": False, "threshold": None, "medians": [vals[n // 2]] if vals else [], "n": n}
    best = None
    for i in range(2, n - 1):
        lo, hi = vals[:i], vals[i:]
        ml, mh = sum(lo) / len(lo), sum(hi) / len(hi)
        w = sum((x - ml) ** 2 for x in lo) + sum((x - mh) ** 2 for x in hi)
        if best is None or w < best[0]:
            best = (w, i)
    i = best[1]
    smed, lmed = int(st.median(vals[:i])), int(st.median(vals[i:]))
    if lmed - smed >= min_sep:
        return {"het": True, "threshold": (smed + lmed) / 2.0, "medians": [smed, lmed], "n": n}
    return {"het": False, "threshold": None, "medians": [int(st.median(vals))], "n": n}


def call_by_effect(best: Optional[dict], *, floor: float = 0.10, min_tot: int = 20) -> dict:
    """Effect-size dupC call at an allele's best C-tract index (`best` from `best_8c_index`). Pure.
    Called iff coverage ≥ min_tot AND frac ≥ floor (the mutant fraction must clear the homopolymer
    error floor — a p-value would flag the ~5% ONT error as 'significant' at amplicon depth)."""
    if not best or best.get("tot", 0) < min_tot:
        return {"called": False, "reason": "low_cov", "best": best}
    return {"called": bool(best["frac"] >= floor), "reason": None, "frac": best["frac"],
            "over_floor": round(best["frac"] - floor, 3), "index": best["index"],
            "n8C": best["n8C"], "tot": best["tot"], "best": best}


def _idx8c(reads, mut_ctract: int = 8, min_tot: int = 8):
    """From reads = [(n_units, strand, [ctract_len,…]),…] build {index: Counter{ctract_len}} (both strands)
    → best C-tract index, plus that index's per-strand mutant fraction. Pure."""
    both = collections.defaultdict(collections.Counter)
    per_strand = {"fwd": collections.defaultdict(collections.Counter),
                  "rev": collections.defaultdict(collections.Counter)}
    for _n, s, cts in reads:
        for i, c in enumerate(cts):
            both[i][c] += 1
            if s in per_strand:
                per_strand[s][i][c] += 1
    best = best_8c_index(both, mut_ctract, min_tot)
    strand_frac = None
    if best:
        idx = best["index"]
        strand_frac = {}
        for s in ("fwd", "rev"):
            c = per_strand[s].get(idx, collections.Counter())
            tot = sum(c.values())
            strand_frac[s] = round(c.get(mut_ctract, 0) / tot, 3) if tot else None
    return best, strand_frac


def call_from_reads(reads, *, mut_ctract: int = 8, floor: float = 0.10, min_abs: float = 0.07,
                    min_ratio: float = 1.4, min_tot: int = 20, min_sep: int = 5,
                    min_strand_frac: float = 0.04, min_cand: int = 5, min_cov_frac: float = 0.03) -> dict:
    """Full per-allele dupC logic from scanned reads (no pysam) — testable end-to-end.
    reads = [(n_units, strand, [ctract lengths per unit]), …].

    Decision — INTERNAL NULL when both alleles are measurable: the mutant allele's dupC fraction must
    exceed the OTHER allele's (the sample's own WT homopolymer rate) by a RATIO (`min_ratio`) and clear
    `min_abs`. A ratio (not an absolute delta) scales with depletion, so a heavily under-amplified mutant
    still calls (Madritsch MP4 0.086/0.05=1.7×; one of our samples 0.083/0.057=1.5× — both depleted ~9-21×), while
    a negative's two ~equal alleles give ~1.0×. When only one allele is measurable (or length-homozygous),
    fall back to the absolute `floor`. A per-strand guard rejects one-strand homopolymer artifacts.

    STATUS — four outcomes, not two, so neither a corrupted null nor a coverage limit is silently reported
    as a negative (each carries `flags` for transparency):
      * "positive"      — dupC called.
      * "negative"      — both length-allele bins measurable, the WT/control allele is at the homopolymer
                          floor (< min_abs, trustworthy null), and no mutant → a CLEAN negative (HG001-4).
      * "undetermined"  — no call, but the internal NULL is corrupted: the WT/control (lower) allele's own 8C
                          rate reaches calling level (≥ min_abs), from a sample-wide strand-asymmetric
                          homopolymer artifact → a real dupC would be indistinguishable from this inflated
                          null → we cannot trust a negative. This is BON (WT allele 0.085 rev-only; its mutant
                          long allele — 0.108, rev-only, tot 37 — is identical to the HG002 NEGATIVE's long
                          allele, so the mutant is uncallable; only the control's unreliability is diagnostic).
                          Hand the case to the deep WGS/urine mode (`dupc_caller`); do NOT fine-tune to call it.
      * "low_coverage"  — no call and a real length-allele bin (≥ min_cand reads) is too thin to evaluate
                          (best index < min_tot ⇒ best=None) → coverage-limited, also UNDETERMINED-in-spirit."""
    split = split_by_length([n for n, _, _ in reads], min_sep=min_sep)
    if split["het"]:
        thr = split["threshold"]
        bins = {"short": [r for r in reads if r[0] <= thr], "long": [r for r in reads if r[0] > thr]}
    else:
        bins = {"single": reads}   # length-homozygous → can't de-dilute; absolute-floor fallback
    per = {}
    for name, rs in bins.items():
        # COVERAGE/FOCALITY gate — a dupC index must be covered by a SUBSTANTIAL share of the allele's reads,
        # not an array-EDGE index reached only by a thin tail (length-split noise). eff = max(absolute floor,
        # min_cov_frac × depth): HG001's false call came entirely from edge indices (best i23 = 6/65 of a
        # 7831-read allele = 0.8% coverage, small-number noise); its well-covered indices sit at the 0.047
        # homopolymer floor. A real dupC (MP2/MP4) is at a deep within-array index (tot 600-840, ≥15% of the
        # allele, bi-strand). 0.03 sits 3.75× above HG001 (0.008) and 4.8× below MP4 (0.146) — wide margins.
        eff_min_tot = max(min_tot, int(min_cov_frac * len(rs)))
        best, sf = _idx8c(rs, mut_ctract, min_tot=eff_min_tot)
        ok_cov = bool(best and best.get("tot", 0) >= eff_min_tot)
        frac = best["frac"] if (best and ok_cov) else None
        # coverage status of THIS allele: measurable / a real-but-thin candidate we can't rule out / noise
        cov = "ok" if ok_cov else ("low" if len(rs) >= min_cand else "none")
        strand_ok = True
        if frac is not None:
            vals = [v for v in (sf or {}).values() if v is not None]
            strand_ok = len(vals) >= 2 and min(vals) >= min_strand_frac
        per[name] = {"n_reads": len(rs), "best": best, "frac": frac, "min_tot_eff": eff_min_tot,
                     "strand_frac": sf, "strand_ok": strand_ok, "cov": cov}
    out = {"n_reads": len(reads), "length_split": split, "mut_ctract": mut_ctract, "alleles": per}
    out.update(_decide(per, floor=floor, min_abs=min_abs, min_ratio=min_ratio))
    return out


def _decide(per, *, floor: float = 0.10, min_abs: float = 0.07, min_ratio: float = 1.4) -> dict:
    """Variant-AGNOSTIC decision core, shared by the dupC (C-tract) and the generic (segment) callers.
    `per` = {allele_name: {frac, strand_ok, cov, …}} where `frac` is the mutant-FEATURE fraction at the
    allele's best index (8C for dupC, or any dictionary variant for the generic caller). Returns
    {called, mut_allele, criterion, status, flags}. The internal-null / focality / strand logic is identical
    whatever the feature — only how `frac` is computed differs (that's what makes this generalize)."""
    out = {"called": False, "mut_allele": None, "criterion": None, "status": "negative", "flags": []}
    measurable = {n: p for n, p in per.items() if p["frac"] is not None}
    null_corrupt = False                                       # internal null (WT allele) untrustworthy?
    if len(measurable) == 2:                                   # internal null (mutant vs the other allele)
        hi = max(measurable, key=lambda n: measurable[n]["frac"])
        lo = min(measurable, key=lambda n: measurable[n]["frac"])
        hf, lf = measurable[hi]["frac"], measurable[lo]["frac"]
        # RATIO, not absolute delta: the mutant allele's variant rate must exceed the WT allele's by a factor
        # (scales with depletion — a heavily under-amplified mutant still shows the same ratio). max(lf,·)
        # keeps the denominator off zero.
        if hf >= min_abs and hf >= min_ratio * max(lf, 0.02) and per[hi]["strand_ok"]:
            out.update(called=True, mut_allele=hi, criterion=f"internal_null({hf:.3f}/{lf:.3f}={hf/max(lf,1e-9):.2f}x)")
        # NULL RELIABILITY — the internal null IS the WT (lower) allele. A "negative" is only trustworthy if
        # that control sits at the error floor. When the control's OWN variant rate reaches calling level
        # (>= min_abs), a sample-wide (strand-asymmetric) basecalling artifact has inflated it → the null is
        # corrupted and a drowned real variant is indistinguishable from a negative. This is exactly what
        # separates BON (WT allele 0.085, rev-only, mutant only 1.27× above it → report UNDETERMINED) from
        # HG002 (WT allele clean 0.062 → its NEGATIVE is trustworthy) even though their mutant alleles are
        # identical (both ~0.10, rev-only, tot ~20-37 = the ONT strand homopolymer artifact, uncallable).
        if lf >= min_abs:
            out["flags"].append(f"control_hot({lf:.3f}>={min_abs})")
            null_corrupt = True
        if not per[lo]["strand_ok"]:                          # WT allele itself strand-imbalanced (transparency)
            out["flags"].append("control_strand_imbalance")
        if not per[hi]["strand_ok"] and hf >= min_abs:        # elevated mutant blocked ONLY by strand (HG002-long)
            out["flags"].append("candidate_strand_imbalance")
    else:                                                      # one allele measurable → absolute floor
        for name, p in per.items():
            if p["frac"] is not None and p["frac"] >= floor and p["strand_ok"]:
                out.update(called=True, mut_allele=name, criterion=f"abs_floor({p['frac']:.3f}>={floor})")
                break
    # Status: a call → positive; else a corrupted internal null (untrustworthy negative = BON) → undetermined;
    # else a real-but-thin allele we couldn't evaluate → low_coverage; else a CLEAN negative (HG001-4).
    if out["called"]:
        out["status"] = "positive"
    elif null_corrupt:
        out["status"] = "undetermined"
    elif any(p["cov"] == "low" for p in per.values()):
        out["status"] = "low_coverage"
    else:
        out["status"] = "negative"
    return out


# ── pysam reader ──────────────────────────────────────────────────────────────

def _scan_reads(bam: str, chrom: str, start: int, end: int, ref: Optional[str] = None):
    """Per read spanning [start,end): (n_units, strand, [ctract_len per unit]). Reuses the C-tract regex."""
    import pysam
    from .detectors.vntr_dupc import _CTX, gene_oriented
    mode = "rc" if str(bam).endswith(".cram") else "rb"
    kw = {"reference_filename": ref} if (mode == "rc" and ref) else {}
    out = []
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        for r in af.fetch(chrom, start, end):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or r.query_sequence is None:
                continue
            matches = list(_CTX.finditer(gene_oriented(r.query_sequence.upper())))
            if not matches:
                continue
            cts = [len(m.group(1)) for m in matches]
            out.append((len(matches), "rev" if r.is_reverse else "fwd", cts))
    return out


def pcr_dupc(bam: str, *, chrom: str, start: int, end: int, ref: Optional[str] = None,
             mut_ctract: int = 8, floor: float = 0.10, min_tot: int = 20, min_sep: int = 5,
             min_cov_frac: float = 0.03) -> dict:
    """Per-allele dupC on a deep amplicon BAM. Separates alleles by unit count, scans dupC per allele,
    calls by effect size above the homopolymer floor with per-strand support and the coverage/focality gate."""
    reads = _scan_reads(bam, chrom, start, end, ref)
    res = call_from_reads(reads, mut_ctract=mut_ctract, floor=floor, min_tot=min_tot, min_sep=min_sep,
                          min_cov_frac=min_cov_frac)
    return {"region": f"{chrom}:{start}-{end}", **res}


def _parse_region(s: str):
    chrom, rng = s.split(":")
    a, b = rng.replace(",", "").split("-")
    return chrom, int(a), int(b)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="muc1_analyzer.pcr_dupc",
                                 description="Per-allele dupC caller for high-depth PCR/amplicon MUC1 data")
    ap.add_argument("-b", "--bam", required=True)
    ap.add_argument("--region", required=True, help="e.g. muc1win:2328605-2330538 (the VNTR span)")
    ap.add_argument("--ref", default=None, help="reference FASTA (required for a CRAM)")
    ap.add_argument("--mut-ctract", type=int, default=8, help="mutant C-tract length: 8=dupC (default), 5=delCC")
    ap.add_argument("--floor", type=float, default=0.10, help="min per-allele frac to call (homopolymer floor)")
    ap.add_argument("--min-tot", type=int, default=20, help="min reads at the best index per allele")
    ap.add_argument("--min-sep", type=int, default=5, help="min unit-count gap to split the two alleles")
    ap.add_argument("--min-cov-frac", type=float, default=0.03,
                    help="coverage/focality gate: the mutant index must have tot >= this fraction of the "
                         "allele depth (rejects array-edge small-number 8C spikes; 0 = off)")
    a = ap.parse_args(argv)
    chrom, start, end = _parse_region(a.region)
    out = pcr_dupc(a.bam, chrom=chrom, start=start, end=end, ref=a.ref, mut_ctract=a.mut_ctract,
                   floor=a.floor, min_tot=a.min_tot, min_sep=a.min_sep, min_cov_frac=a.min_cov_frac)
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
