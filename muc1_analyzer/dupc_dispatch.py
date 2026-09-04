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
import re
from collections import Counter

#: contig name -> copy number, for the replication annotation
CONTIG_LEN_RE = re.compile(r"MUC1_VNTR_(\d+)repeats")
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
                  pon=None, fast: bool = False, f_null="adaptive") -> dict:
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
        d = scan_dupc_positional(bam, genome_ref=ref, mut_len=mut_ctract,
                                 f_null=f_null)                            # genomic positional, PoN-free
        # The null DIAGNOSIS travels with the verdict. `scan_dupc_positional` always reports it (even in
        # fixed mode), and this branch used to drop it — so a `run` output could not say which null it
        # was judged against, which is exactly what `adaptive` made a question worth asking. Measured on
        # the release gate: the f_null column came back None on all 21 samples.
        # ⚠ NO repeat number is derivable here, and the output must SAY so rather than leave a silent
        # None downstream. `best.index` is the ordinal of a C-tract in the hg38 locus, where the VNTR is
        # COLLAPSED — there are no array units to count (measured on the gate: index 2-3 on carriers whose
        # arrays hold 44-70 repeats). A carrier detected by this route alone therefore has a verdict and
        # no onset axis. To get the position, run the sample from its FASTQ: the VNTR-space route returns
        # a scaffold `ref_pos`, which converts (cf. `_unit_index`). That is the nominal route for an
        # amplicon anyway — `docs/benchmark_comparison.md` §8.
        res = {"called": bool(d.get("called")), "tier": "positional", "best": d.get("best"),
               "n_span": d.get("n_span"),
               "variant": ({"repeat": None, "repeat_index_space": "unavailable_genomic",
                            "ctract_index": (d.get("best") or {}).get("index"),
                            "label": "59dupC" if mut_ctract == 8 else f"mut_len{mut_ctract}",
                            "reason": "the hg38 locus is VNTR-collapsed: no array unit to number; "
                                      "re-run from the FASTQ for the repeat position"}
                           if d.get("called") else None),
               "f_null_used": d.get("f_null_used"), "f_null_mode": d.get("f_null_mode"),
               "f_null_observed": d.get("f_null_observed"), "rho": d.get("rho"),
               "noise_flags": d.get("noise_flags")}
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
                        region=f"{chrom}:{start}-{end}", probe=probe, f_null=f_null)
        caller = f"dupc_caller(auto-anchor,{probe})"
    return {"regime": regime, "caller": caller, "n_locus_reads": n, "amplicon_min": amplicon_min,
            "region": f"{chrom}:{start}-{end}", "result": res}


def _unit_index(vntr_ref: str, contig: str, ref_pos):
    """Scaffold C-run coordinate → 1-based ARRAY UNIT index, reading the 5' flank off the reference
    itself (the flank length is recorded nowhere the caller can read). None on any failure — an absent
    repeat number is honest, a C-run rank published as a repeat number is not."""
    if ref_pos is None or not vntr_ref or not contig:
        return None
    try:
        import pysam
        from .build_vntr_ref import MOTIF_A
        from .config import FLANK_TO_PHYSICAL_OFFSET
        from .detectors.runlen_shift import flank5_len, unit_index_from_ref_pos
        with pysam.FastaFile(vntr_ref) as fa:
            seq = fa.fetch(contig)
        return unit_index_from_ref_pos(ref_pos, flank5_len(seq, MOTIF_A),
                                       offset=FLANK_TO_PHYSICAL_OFFSET)
    except Exception:
        return None


#: C-tract length → the variant it means. The run-length caller is parameterised by the MUTANT tract
#: length, so every all-C family is the same machinery with a different target: 7C is wild type, 8C is
#: 59dupC, 5C is 58_59delCC.
#: ⚠⚠ 5C IS OFF BY DEFAULT — CLOSED on 2026-08-04, on 182 samples, and re-opening it needs a different
#: STATISTIC rather than a different threshold.
#:
#: The contraction statistic has exactly one detectable feature on this terrain and it is not a variant.
#: 113 calls: **111 at `ref_pos` 4954, 2 at 4894** (the same C-run one 60 bp unit earlier), and **zero
#: anywhere else** — across 119 negatives, 58 carriers and 5 unknowns alike. The only powered delCC
#: carrier, fires at 4954 too, so masking the site would buy specificity 1.000 and ZERO true positives:
#: the probe goes silent, it does not become correct.
#:
#: The reference is not at fault — 4954 reads `...CTCCACCG|CCCCCCC|A|GCCCACGG...`, an ordinary unit-ending
#: C-run identical to the ~45 others. It is the FIRST canonical C-run of the array, which is where
#: minimap2 places the indel when a read's copy number differs from its contig, and an indel truncates the
#: walked run DOWNWARD: invisible to an expansion target (`x >= 8` = 0, which is why the 8C path never met
#: this) and universal to a contraction one. CONFIRMED by construction: scaffold a sample on contigs that
#: match its measured allele lengths and the site disappears — there is no copy-number mismatch left to
#: place an indel.
#:
#: And with the artifact gone there is nothing under it. The only powered delCC carrier scaffolded
#: on its clinical 75|77: NOT CALLED, best position k=1/23 (frac 0.043, p_bonf 1.0) on the 75 and k=0
#: everywhere on the 77. Meanwhile two NEGATIVES, both faithfully scaffolded, call at ranks 51 and 59.
#: Sensitivity 0/1 with false positives — the probe is worse than silent.
#:
#: Earlier readings of this failure were wrong and are kept only as warnings. "1 FP / 0 TP on 30 amplicons,
#: signal under the noise floor" (2026-08-02) came from an inequality running backwards. The retention
#: floor was refuted, and so was the flank — but the ALIGNMENT explanation above, which `f_null` seemed to
#: refute on 2026-08-04, was right: that refutation read a table built from replays that had silently
#: re-picked their own scaffolds. See `docs/MUC1_log.md`.
CTRACT_FAMILIES = {8: "59dupC", 5: "58_59delCC", 6: "58_59delC", 11: "dupCCCC"}


def scan_segments(vntr_bam: str, *, max_reads: int = 2000) -> dict:
    """Per-read inter-anchor segment types — the families the C-run caller CANNOT see. Impure.

    `58_59insG` and `52dupG` insert a G INSIDE the C-tract, which BREAKS the run, so a `(C+)` run-length
    caller misses them by construction (proven 2026-07-05). `classify_segment` reads the segment whole and
    types it by length AND composition, which is why it sees them.

    ⚠ Returns COUNTS, not a verdict. The run-length path carries a measured specificity (1.000 on 92
    genomic negatives); this one does not yet, so it feeds `frameshift_alarms` — visible, not deciding —
    exactly as the attested-variant rule of 2026-07-29 requires. That design is what kept it from doing
    damage when it was measured.

    ⚠⚠ MEASURED AND IT DOES NOT DISCRIMINATE (2026-08-02, 30 raw amplicons). The fractions are the same
    in carriers and in negatives — delC 2.6-6.5 %, dupC 1.0-4.4 % everywhere — and neither `dupG` nor
    `insG` appears at all, including in the insG carrier (2.7 %/1.4 %) and the dupG carrier (
    4.9 %/2.9 %). It is reading ONT homopolymer error, not variants. Kept as a QC number; it is NOT a
    detector on this substrate and must not be presented as one.
    """
    import pysam
    from .detectors.vntr_segment import segment_runs
    tally, n_reads = Counter(), 0
    with pysam.AlignmentFile(vntr_bam, "rb") as af:
        for r in af.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.query_sequence:
                continue
            n_reads += 1
            if n_reads > max_reads:
                break
            for c in segment_runs(r.query_sequence.upper()):
                tally[c.get("type")] += 1
    wt = tally.get("WT", 0)
    fs = {k: v for k, v in tally.items() if k not in ("WT", "other") and v}
    total = sum(tally.values()) or 1
    return {"n_reads": n_reads, "n_segments": total, "wt": wt,
            "frameshift_types": dict(sorted(fs.items(), key=lambda kv: -kv[1])),
            "fractions": {k: round(v / total, 4) for k, v in
                          sorted(fs.items(), key=lambda kv: -kv[1])}}


#: Gates for the opt-in second detector, measured 2026-08-07 as its zero-false-positive point on the
#: 180-PCR lot (43 carrier / 116 negative subjects). ⚠ CHOSEN ON THE COHORT THEY WERE SCORED ON — this is
#: why the probe is opt-in and not a default. Price on independent negatives before that changes.
PROBE_GATES = {"min_abs": 0.12, "floor": 0.15}


def dispatch_dupc_vntr(vntr_bam: str, vntr_ref: str, *, lengths=None, mut_len: int = 8,
                       wt_len: int = 7, extra_mut_lens=(), segments: bool = True,
                       variant_probe: str = None, probe_gates: dict = None,
                       replicate_check: bool = True) -> dict:
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
    per_allele, best, power, best_len = {}, None, {}, mut_len
    # EVERY all-C family, not just dupC: same caller, same statistics, one parameter apart.
    for target in (mut_len, *[m for m in (extra_mut_lens or ()) if m != mut_len]):
        for contig in contigs:
            key = contig if target == mut_len else f"{contig}@{target}C"
            try:
                r = R.call_bam(vntr_bam, contig, vntr_ref, mut_len=target, wt_len=wt_len)
            except Exception as e:
                per_allele[key] = {"error": f"{type(e).__name__}: {e}"[:80]}
                continue
            per_allele[key] = r
            cand = r.get("candidate") or {}
            n_cover = cand.get("n") or 0
            need = _min_reads(0.3, target=0.90) or 12      # reads to reach 90 % power at d=0.3
            if target == mut_len:                          # power is reported on the primary target
                power[contig] = {"n_cover": n_cover, "adequately_powered": bool(n_cover >= need),
                                 "reads_for_power0.90": {"d=0.3": need}}
            if r.get("called") and (best is None or cand.get("p_bonf", 1) < best[1].get("p_bonf", 1)):
                best, best_len = (contig, cand), target

    seg = None
    if segments:
        try:
            seg = scan_segments(vntr_bam)
        except Exception as e:
            seg = {"error": f"{type(e).__name__}: {e}"[:80]}

    # ── the SECOND detector, opt-in ───────────────────────────────────────────────────────────────
    # `runlen_shift` tests a C-tract run-length shift over ~81 positions under Bonferroni; `pcr_variant`
    # types each array unit by its LITERAL sequence and searches for the best repeat index. Measured
    # 2026-08-07, with BOTH arms scaffolded by the same code: the probe recovers ONE carrier the
    # run-length caller misses, at ZERO added false positives — carrier axis 0.61 -> 0.64 on the 151
    # cleaned samples and 0.50 -> 0.57 on the 30 raw ones.
    # ⚠ An earlier "+2 / 0.814" was RETRACTED: it differenced a STORED verdict against a freshly
    # recomputed one on a DIFFERENT scaffold, and the extra carrier turned out to be the scaffold fix,
    # not the probe. On a correct scaffold the two detectors largely agree, which is why the marginal
    # value is one carrier and not more.
    # ⚠ OPT-IN, and it must stay opt-in until `PROBE_GATES` is priced on independent negatives: those
    # values were chosen on the cohort they were scored on (`decisions.md` 2026-07-28).
    probe_hit = None
    if variant_probe:
        _g = {**PROBE_GATES, **(probe_gates or {})}
        try:
            # ⚠ `_scan_segments` + `call_variant_from_reads`, NOT the `pcr_variant()` wrapper: the wrapper
            # does not accept `min_abs`, so passing the measured gates through it raises TypeError. The
            # first integration did exactly that, the except below swallowed it, and the probe silently
            # never called — the union looked wired and measured nothing (2026-08-07). This is the same
            # code path the measurement itself ran through, which is the only way the shipped behaviour
            # and the priced behaviour are the same thing.
            from .pcr_variant import _scan_segments, call_variant_from_reads
            import pysam as _ps
            with _ps.AlignmentFile(vntr_bam, "rb") as _fh:
                _len = {c: _fh.get_reference_length(c) for c in contigs}
            for c in contigs:
                reads = _scan_segments(vntr_bam, c, 1, int(_len[c]), None)
                r = call_variant_from_reads(reads, variant=variant_probe, **_g)
                if str(r.get("status", "")).lower().startswith("pos"):
                    probe_hit = {"contig": c, "variant": variant_probe, "status": r.get("status"),
                                 "alleles": r.get("alleles"), "gates": _g}
                    break
        except Exception as e:
            # LOUD. A swallowed probe failure is indistinguishable from a probe that found nothing, and
            # that is precisely how the first integration passed a full 152-sample benchmark measuring
            # the old pipeline.
            probe_hit = {"error": f"{type(e).__name__}: {e}"[:120], "gates": _g}
            print(f"[dupc_dispatch] ⚠ variant probe FAILED: {type(e).__name__}: {e}", file=sys.stderr)

    # POWER IS NOT ONLY A NEGATIVE'S BUSINESS. Until 2026-08-07 `power_by_allele` gated the NEGATIVE
    # branch alone, so a call came out "CONFIRMED" with nothing said about the other allele — even when
    # that allele carried 10 reads and was never testable. FER: called on its 43-copy scaffold while its
    # 82-copy allele had 11 reads and no array unit reached the coverage floor. Whatever the call means,
    # "this subject is a heterozygous carrier and nothing else" is NOT what was measured, and a positive
    # that hides an untested allele is the same silence as a negative that hides one.
    underpowered = sorted(c for c, r in per_allele.items()
                          if (r.get("n_tested") or 0) and not power.get(c, {}).get("adequately_powered"))
    # AN ALLELE NEVER SCAFFOLDED IS NOT AN ALLELE THAT PASSED. `underpowered` can only speak about alleles
    # that were TESTED, so a subject scaffolded on ONE contig came out `fully_powered: True` with half the
    # genotype missing — measured on the VNTRtools trio 2026-08-07, where both parents were reported
    # NEGATIVE on a single scaffold (18 and 26 reads) while the index carried two alleles. A dominant
    # variant on the unexamined allele would have been invisible, and nothing in the output said so.
    # A single scaffold is legitimate for a length-HOMOZYGOTE, so this reports rather than judges.
    n_scaffolds = len([c for c in (contigs or []) if "@" not in str(c)])
    single_allele = n_scaffolds < 2

    # SCAFFOLD REPLICATION — an ANNOTATION, never a gate (measured 2026-08-07).
    # A real variant is carried by the READS and survives a change of reference length; a padding
    # artifact is made by the PLACEMENT (aligning an N-copy allele to a shorter contig forces the
    # surplus into gaps, and a gap inside a C-run is a run-length shift) and evaporates. On the cleaned
    # arm 24 of 27 called carriers replicate on at least one neighbouring scaffold; the lone false
    # positive, FER, replicates on none.
    # ⚠ NOT a filter: requiring one replication would drop 3 carriers to remove 1 false positive — the
    # same bad trade as the proximity scaffold rule. It is reported so a reader can weigh a call, and
    # deliberately does not change one.
    replication = None
    if best and replicate_check:
        try:
            from .allele_scaffold import contig_counts as _cc
            _m = CONTIG_LEN_RE.search(best[0])
            if _m:
                _L0, _counts = int(_m.group(1)), _cc(vntr_bam)
                _rep = _tried = 0
                for _d in (-2, -1, 1, 2):
                    _nm = f"MUC1_VNTR_{_L0 + _d}repeats"
                    if _counts.get(_nm, 0) <= 0:
                        continue
                    _tried += 1
                    _rep += bool(R.call_bam(vntr_bam, _nm, vntr_ref,
                                            mut_len=best_len, wt_len=wt_len).get("called"))
                replication = {"replicated": _rep, "neighbours_tested": _tried,
                               "note": ("call did NOT survive a change of scaffold length — the "
                                        "signature of a padding artifact, not of a variant"
                                        if _tried and not _rep else None)}
        except Exception as e:
            replication = {"error": f"{type(e).__name__}: {e}"[:80]}

    if best:
        contig, cand = best
        interp, called = "CONFIRMED", True
        # `position` is the RANK among the scaffold's C-runs, not the repeat number — units without a
        # C-run are skipped. Publishing the rank as `repeat` fed the ONSET AXIS a number from the wrong
        # space: measured on one family, the same inherited variant on the same 45-repeat allele read
        # 11 here and 16 through the consensus (onset 0.7556 vs 0.6444). Convert; if it cannot be
        # converted, emit NO repeat rather than a number in the wrong unit.
        _unit_idx = _unit_index(vntr_ref, contig, cand.get("ref_pos"))
        variant = {"repeat": _unit_idx, "repeat_ctract_rank": cand.get("position"),
                   "repeat_index_space": "array_unit" if _unit_idx is not None else "unconverted",
                   "ref_pos": cand.get("ref_pos"),
                   "label": CTRACT_FAMILIES.get(best_len, f"mut_len{best_len}"),
                   "ctract_len": best_len,
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

    # PROVENANCE, never a silent merge. A union call that cannot say which detector made it is
    # unauditable — and the two do not carry the same evidence: `runlen_shift` yields a scaffold
    # `ref_pos` (hence a repeat number and the onset axis), the probe yields an array INDEX and no p.
    called_by = "runlen_shift" if best else None
    if probe_hit and not probe_hit.get("error"):
        called_by = "both" if best else "pcr_variant"
        if not best:
            called, interp = True, "CONFIRMED_BY_PROBE"
            contig = probe_hit["contig"]
            _a = max((v for v in (probe_hit.get("alleles") or {}).values()
                      if v.get("frac") is not None), key=lambda v: v["frac"], default={})
            variant = {"repeat": None, "repeat_index_space": "array_unit_probe",
                       "repeat_probe_index": (_a.get("best") or {}).get("index"),
                       "label": CTRACT_FAMILIES.get(mut_len, variant_probe),
                       "frac": _a.get("frac"), "n": _a.get("n_reads"),
                       "reason": "called by the segment-typing probe, which yields an array index and "
                                 "no p-value; runlen_shift did not call"}

    return {"regime": "vntr", "caller": "runlen_shift(per-allele)"
            + (f"+pcr_variant({variant_probe})" if variant_probe else ""), "region": None,
            "scaffolds": contigs,
            "result": {"assessed": bool(contigs), "called": called, "interpretation": interp,
                       "method": "runlen_shift", "called_by": called_by, "probe": probe_hit,
                       "carrier_contig": contig, "variant": variant,
                       # Which alleles could NOT be assessed, whatever the verdict. On a CALL this is the
                       # difference between "heterozygous carrier" and "carrier on one allele, the other
                       # never tested"; on a NEGATIVE it is what `interpretation` already encodes.
                       "underpowered_alleles": underpowered,
                       "alleles_scaffolded": n_scaffolds,
                       "single_allele_examined": single_allele,
                       # BOTH conditions: every tested allele powered AND both alleles actually tested.
                       "fully_powered": (not underpowered) and not single_allele,
                       "coverage_note": ("only ONE allele was scaffolded — either the subject is "
                                         "length-homozygous or the second allele was not resolved; a "
                                         "variant on an unexamined allele cannot be excluded"
                                         if single_allele else None),
                       "scaffold_replication": replication,
                       "power_by_allele": power, "per_allele": per_allele,
                       "ctract_families_tested": sorted({mut_len, *(extra_mut_lens or ())}),
                       "segment_scan": seg,
                       # The G-insertions the C-run caller cannot see, as ALARMS: visible, not deciding,
                       # until this path has its own measured specificity.
                       "frameshift_alarms": sorted((seg or {}).get("frameshift_types", {}) or {},
                                                   key=lambda k: -(seg or {}).get(
                                                       "frameshift_types", {}).get(k, 0)),
                       "reason": None if contigs else "no scaffold contig could be chosen"}}


def _parse_region(s: str):
    chrom, rng = s.split(":")
    a, b = rng.replace(",", "").split("-")
    return chrom, int(a), int(b)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="muc1_analyzer.dupc_dispatch",
                                 description="Auto-route a MUC1 dupC/frameshift call by input regime (depth)")
    ap.add_argument("-b", "--bam", default=None, help="genomic bam/cram (not needed with --vntr-bam)")
    ap.add_argument("--region", default=None, help="e.g. chr1:155188000-155192000 (genomic route)")
    # VNTR-space REPLAY: `run` already wrote a scaffolded `.vntr.bam`, and `prepare` is the expensive
    # step. Re-measuring a different C-tract target on an existing one costs the C-run scan alone, so a
    # target can be priced on the cohort without re-aligning anything.
    ap.add_argument("--vntr-bam", default=None,
                    help="VNTR-space route on an ALREADY scaffolded bam (skips prepare); needs --vntr-ref")
    ap.add_argument("--vntr-ref", default=None, help="the 150-contig VNTR reference FASTA")
    ap.add_argument("--extra-ctract", default="", metavar="5,11",
                    help="additional C-tract targets to test alongside --mut-ctract (comma-separated)")
    # The original run scaffolds on the ARBITER's lengths. A replay that lets `allele_contigs` choose on
    # its own can land on a different contig, and on a carrier the scaffold decides whether the mutant
    # allele is looked at at all — so a replay meant to reproduce a verdict must be given the same lengths.
    ap.add_argument("--lengths", default=None, metavar="40,70",
                    help="the arbiter's called allele lengths, as the run used them (comma-separated)")
    ap.add_argument("--ref", default=None, help="reference FASTA (required for a CRAM)")
    ap.add_argument("--variant", default=None, help="PCR regime: call this dictionary variant (insG/dupG/…) via pcr_variant")
    ap.add_argument("--variant-probe", dest="variant_probe", default=None, metavar="dupC",
                    help="VNTR route: ALSO run pcr_variant for this family and UNION the verdicts. "
                         "Measured +1 carrier at 0 added false positives, both arms scaffolded alike "
                         "(carrier axis 0.61 -> 0.64 cleaned, 0.50 -> 0.57 raw); opt-in because its "
                         "gates were chosen on that same cohort")
    ap.add_argument("--pon", default=None, help="WGS regime: pon_dupc.json for the statistical caller")
    ap.add_argument("--mut-ctract", type=int, default=8, help="8 = dupC (default), 5 = delCC")
    ap.add_argument("--amplicon-min", type=int, default=200, help="reads spanning the locus ≥ this ⇒ PCR regime")
    ap.add_argument("--fast", action="store_true",
                    help="RECOMMENDED for WGS/AS: positional path — ~1500× faster (~0.2 s vs ~5-7 min/sample) "
                         "with IDENTICAL specificity (0 FP on 233 MAT+1000G). Without it you get the ~7-min "
                         "context+PoN scan (only a marginal per-context-null sensitivity gain).")
    a = ap.parse_args(argv)
    extra = tuple(int(x) for x in a.extra_ctract.split(",") if x.strip())
    if a.vntr_bam:
        if not a.vntr_ref:
            ap.error("--vntr-bam needs --vntr-ref (the 150-contig VNTR reference)")
        lens = [int(x) for x in a.lengths.split(",") if x.strip()] if a.lengths else None
        out = dispatch_dupc_vntr(a.vntr_bam, a.vntr_ref, lengths=lens,
                                 mut_len=a.mut_ctract, extra_mut_lens=extra,
                                 variant_probe=a.variant_probe)
        print(json.dumps(out, indent=2, default=str))
        return 0
    if not (a.bam and a.region):
        ap.error("the genomic route needs -b/--bam and --region (or use --vntr-bam/--vntr-ref)")
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


def merge_dupc_position(genomic: dict, vntr: dict) -> dict:
    """Give a GENOMIC dupC verdict the array POSITION only the VNTR-space route can produce. Pure.

    Why both routes on the same aligned input: `run` already builds the VNTR-space BAM at `prepare`, so
    the position was never missing from the DATA — only from the ROUTING. An aligned CRAM sent its dupC to
    the genomic path, where the hg38 array is collapsed and no repeat number exists, and the VNTR path that
    does yield one was never asked (measured 2026-07-29 on the cohort gate: P2 and P6 called carriers with
    `onset_index=None`).

    Rules, in order:
      · the GENOMIC verdict decides `called` — it is the route whose specificity was measured on the
        independent negatives; this function never changes a verdict;
      · the position is taken ONLY when both routes called. A position from a route that says "not a
        carrier" describes nothing, and silently borrowing it would manufacture agreement;
      · disagreement is RECORDED (`position_source`, `routes_agree`) rather than smoothed over.
    """
    out = dict(genomic or {})
    res = dict(out.get("result") or {})
    vres = (vntr or {}).get("result") or {}
    g_called, v_called = bool(res.get("called")), bool(vres.get("called"))
    res["vntr_route"] = {"called": v_called, "interpretation": vres.get("interpretation"),
                         "carrier_contig": vres.get("carrier_contig"), "variant": vres.get("variant")}
    res["routes_agree"] = (g_called == v_called)
    var = res.get("variant") or {}
    has_pos = var.get("repeat") is not None
    if g_called and v_called and not has_pos:
        vvar = vres.get("variant") or {}
        if vvar.get("repeat") is not None:
            merged = dict(var)
            merged.update({k: vvar.get(k) for k in ("repeat", "repeat_index_space",
                                                    "repeat_ctract_rank", "ref_pos")})
            merged["position_source"] = "vntr_route"
            merged.pop("reason", None)
            res["variant"] = merged
    elif g_called and not v_called:
        res.setdefault("variant", {})
        if isinstance(res["variant"], dict):
            res["variant"]["position_source"] = "unavailable_routes_disagree"
    out["result"] = res
    return out
