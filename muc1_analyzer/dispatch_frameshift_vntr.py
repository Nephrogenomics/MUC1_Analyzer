#!/usr/bin/env python3
"""dispatch_frameshift_vntr.py — VNTR-space frameshift verdict via the VALIDATED detector chain.

Drop-in alternative to ``dupc_dispatch.dispatch_dupc_vntr`` for the PCR-FASTQ route. Instead of the
run-length-shift caller (``detectors.runlen_shift`` — 7 false positives on the LR-PCR cohort), it runs
the validated chain

    allele_lengths.two_alleles  ->  frameshift_vntr.call     (four-verdict; Se/Sp = 1.000 on the 171-sample cohort)

and REMAPS the result into the exact ``dispatch_dupc_vntr`` shape, so every downstream consumer works
unchanged: ``run_pcr`` (interpretation / carrier length), ``caller`` via ``--dupc-json`` (variant repeat
+ label, named authoritatively by the consensus), and ``score`` via ``carrier_contig``.

The wiring reproduces ``eval_frameshift_cohort.scan_lengths`` VERBATIM — the bench that measured the
100 % Se/Sp — so the shipped behaviour and the priced behaviour are the same thing. Naming is NOT done
here: the detector only delivers the VERDICT + the carrier repeat (already in array-unit space, so no
``_unit_index`` conversion is needed); ``caller`` names the variant from the reconstructed consensus.
"""
from __future__ import annotations

import os
import sys

# The detector chain lives at the REPOSITORY ROOT and imports its siblings by top-level name
# (``import vntr_raw_length``, ``from consensus_arbiter import ...``, ``import frameshift_vntr``). Make
# the repo root importable from inside the package, mirroring what ``dispatch._delegate_fromfastq`` does
# with ``os.chdir(here)`` for MUC1_Analyzer_fromfastq.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import allele_lengths as AL      # noqa: E402  (repository-root module)
import frameshift_vntr as FS     # noqa: E402  (repository-root module)


# Verdict order and positivity — identical to frameshift_vntr / eval_frameshift_cohort.
VERDICT_RANK = {"CALLED": 4, "POSITIVE_LOW_DEPTH": 3, "NEEDS_IGV": 2,
                "NEEDS_IGV_FOR_ZYGOSITY": 1, "NEG": 0}
POSITIVE_VERDICTS = ("CALLED", "POSITIVE_LOW_DEPTH")

# Verdict -> the `interpretation` string the rest of the pipeline reads (dispatch_dupc_vntr vocabulary).
_INTERP = {
    "CALLED": "CONFIRMED",
    "POSITIVE_LOW_DEPTH": "CONFIRMED_LOW_DEPTH",
    "NEEDS_IGV": "NEEDS_IGV",
    "NEEDS_IGV_FOR_ZYGOSITY": "NEEDS_IGV_FOR_ZYGOSITY",
    "NEG": "NEGATIVE",
}


def _carrier_reads(bam, ref, contig, rep, net, min_mapq=0, ont_slack=4):
    """Read names carrying the called net indel at array unit `rep`, over the ONT overcall window — the
    SAME window frameshift_vntr.call aggregates on. Re-derived from frameshift_vntr WITHOUT modifying it
    (calls its public vntr_bounds / per_repeat_shift), so the validated detector file stays byte-for-byte
    unchanged. These reads let caller.py phase a size-homozygous heterozygous carrier BY THE VARIANT."""
    import pysam
    try:
        fa = pysam.FastaFile(ref)
        n_rep = FS.contig_n_repeats(contig)
        if n_rep is None:
            return []
        vstart, vend = FS.vntr_bounds(fa, contig, n_rep)
        if vstart is None:
            return []
        _shift, _cover, _ins, carriers = FS.per_repeat_shift(bam, contig, vstart, vend, min_mapq=min_mapq)
    except Exception:
        return []
    if net > 0:
        win = range(net, net + ont_slack + 1)           # insertion: absorb larger overcalls
    else:
        win = range(net - ont_slack, net + 1)           # deletion: absorb larger undercounts
    names = set()
    for d in win:
        names |= carriers.get((rep, d), set())
    return sorted(names)


def _scan_lengths(bam, ref, lengths, min_mapq, zygosity_uncertain=False, force_pooled=False):
    """Strongest call across the allele contig(s) — VERBATIM from eval_frameshift_cohort.scan_lengths.

    ``same_length_alleles=True`` when a SINGLE length is returned (the two alleles share one contig, so
    phasing / n//2 applies). ``zygosity_uncertain`` is forwarded only in that single-length case, so the
    detector can flag NEEDS_IGV_FOR_ZYGOSITY when the n//2 bascule rests on < zygosity_depth reads.

    Returns the strongest call by verdict rank then VAF, or None. The VERDICT is taken at the SAMPLE
    level (``res['verdict']`` — which already applies the per-sample zygosity cap), while the position and
    variant come from the top call (``res['calls'][0]``); this is exactly what the bench does.
    """
    # ``same_len`` gates the halved(n//2) per-allele VAF, valid ONLY for a CONFIRMED size-homozygous single
    # length (both alleles same length -> a carrier is ~50 % of reads -> halve). ``force_pooled`` disables it
    # when the single length is NOT trustworthy — e.g. it came from the cassette fallback, which cannot
    # resolve a Delta=1 length-heterozygote (VNTRtools' 44|45 mother collapsed to 45, and the halved
    # doubling turned the ~0.34 background into a spurious 0.68 POSITIVE). Pooled is the safe estimate there.
    same_len = (len(lengths) == 1) and not force_pooled
    best = None
    for L in lengths:
        contig = f"MUC1_VNTR_{L}repeats"
        try:
            res = FS.call(bam, ref, contig, same_length_alleles=same_len,
                          zygosity_uncertain=(zygosity_uncertain and same_len), min_mapq=min_mapq)
        except Exception:
            continue
        if res.get("error") or not res.get("calls"):
            continue
        top = res["calls"][0]                          # calls are sorted by verdict rank then VAF
        verdict = res.get("verdict", "NEG")
        cand = {"contig": contig, "length": L, "repeat": top["repeat"], "net_bp": top["net_bp"],
                "variant": top.get("variant"), "vaf": top.get("vaf"), "vaf_mode": top.get("vaf_mode"),
                "raw_frac": top.get("raw_frac"), "kind": top.get("kind"), "verdict": verdict}
        key = (VERDICT_RANK.get(verdict, 0), cand["vaf"] or 0)
        if best is None or key > (VERDICT_RANK.get(best["verdict"], 0), best["vaf"] or 0):
            best = cand
    return best


def _consensus_frameshift_result(vntr_bam, vntr_ref, lengths, note, min_mapq=0, low_conf=False):
    """FOREIGN-amplicon + PacBio path (VNTRtools & co): decode the frameshift from a CASSETTE-binned dense
    per-allele CONSENSUS (``vntr_scaffold.scaffold_both_alleles`` with the cassette ``copies_fn``), the way
    ``run`` handles a PacBio smear — instead of the ONT-calibrated STATISTICAL ``frameshift_vntr.call``,
    which mis-fires on the HiFi C-tract +1 background. A consensus needs the indel in > 50 % of an allele's
    reads, so the ~0.34 background is rejected by construction. Reached ONLY under the strict foreign+PacBio
    gate, so our native ONT LR-PCR is never routed here.

    ``low_conf`` (the cassette arbiter's ``long_low_conf``) marks that a candidate allele's length rests on
    too few spanning reads to trust a dense consensus. When it is set AND no carrier consensus was found, the
    allele can be neither confirmed nor excluded, so the verdict is the PRUDENT ``NEEDS_IGV_FOR_ZYGOSITY``
    rather than a confident ``NEG`` — a confident negative on an under-powered allele is exactly the
    dangerous false negative. A found carrier still wins (``CALLED``); ``low_conf`` only guards the else."""
    import tempfile
    from . import vntr_scaffold as _SC
    import vntr_raw_length as _VRL
    _cass = lambda seq, maxmm, offset: _VRL.read_copies_cassette(seq, maxmm)
    contigs = [f"MUC1_VNTR_{int(L)}repeats" for L in lengths]
    try:
        sc = _SC.scaffold_both_alleles(
            vntr_bam, vntr_ref, workdir=tempfile.mkdtemp(prefix="muc1_hifi_"),
            offset=4, alleles=[int(L) for L in lengths], pacbio=True, min_mq=min_mapq,
            threads=4, copies_fn=_cass, snp_phase=True, min_snps=3)
    except Exception as e:
        sc = {"available": False, "note": f"scaffold error: {type(e).__name__}: {e}", "per_allele": []}
    per = sc.get("per_allele", [])
    # SNP phasing may have split a length-collapsed cluster into two alleles -> reflect its lengths/contigs.
    eff_lengths = sc.get("alleles") or list(lengths)
    contigs = [f"MUC1_VNTR_{int(L)}repeats" for L in eff_lengths]
    phased_by_snp = bool(sc.get("phased_by_snp"))
    carriers = [p for p in per if p.get("indel")]          # allele whose dense consensus shows a frameshift
    method_note = (note + " | consensus-scaffold(cassette,PacBio)").strip(" |")
    if carriers:
        c0 = max(carriers, key=lambda p: p.get("n_aligned", 0) or 0)
        verdict, called, status, interp = "CALLED", True, "positive", "CONFIRMED"
        carrier_contig = c0.get("contig")
        variant = {"repeat": None, "repeat_index_space": "array_unit", "label": c0.get("indel"),
                   "net_bp": None, "kind": None, "vaf": None, "vaf_mode": "consensus",
                   "raw_frac": None, "verdict": verdict,
                   "n_carrier_reads": c0.get("n_aligned"), "carrier_reads": []}
    elif low_conf:
        # No carrier consensus, but a candidate allele is under-powered (long_low_conf): we can neither
        # confirm nor exclude a frameshift on it -> PRUDENT verdict, never a confident NEG.
        verdict, called, status, interp = "NEEDS_IGV_FOR_ZYGOSITY", False, "needs_igv", "INCONCLUSIVE"
        carrier_contig, variant = None, None
    else:
        verdict, called, status, interp = "NEG", False, "negative", "NEGATIVE"
        carrier_contig, variant = None, None
    underpowered = list(contigs) if (low_conf and not carriers) else []
    cov_note = ("candidate allele under-powered (cassette long_low_conf): dense consensus not "
                "trustworthy — frameshift neither confirmed nor excluded, review in IGV"
                if (low_conf and not carriers) else None)
    return {
        "regime": "vntr", "caller": "frameshift_vntr(consensus-scaffold/PacBio-foreign)", "region": None,
        "scaffolds": contigs,
        "result": {
            "assessed": bool(per), "called": called, "status": status, "interpretation": interp,
            "verdict": verdict, "method": "consensus_scaffold", "called_by": "vntr_scaffold" if carriers else None,
            "carrier_contig": carrier_contig, "variant": variant,
            "alleles_lengths": [int(L) for L in eff_lengths], "alleles_scaffolded": len(contigs),
            "single_allele_examined": (len(eff_lengths) == 1), "coverage_note": cov_note,
            "underpowered_alleles": underpowered, "scaffold_replication": None,
            "length_detail": {"note": method_note + (" | snp-phased" if phased_by_snp else ""),
                              "arbitration": None},
            "reason": (cov_note if (low_conf and not carriers)
                       else (None if per else "cassette-scaffold found no spanning reads")),
        },
    }


def dispatch_frameshift_vntr(vntr_bam, vntr_ref, *, lengths=None, min_mapq=0,
                             variant_probe=None, **_ignored):
    """Frameshift/dupC verdict IN VNTR SPACE via the validated detector — ``dispatch_dupc_vntr`` shape.

    Parameters
    ----------
    vntr_bam, vntr_ref : the VNTR-aligned BAM (MUC1 profile) and the multi-contig VNTR reference.
    lengths : the two allele lengths if already measured; otherwise they are measured here with
        ``allele_lengths.two_alleles`` (the validated caller), which also yields the picker note that
        drives the zygosity flag.
    variant_probe : accepted for signature parity with ``dispatch_dupc_vntr`` and IGNORED — this detector
        already covers every frameshift family (dupC, 60dupA, 58_59insG, 58_59delCC, del8_27, ...).

    Returns the same top-level dict as ``dispatch_dupc_vntr`` (``regime`` / ``caller`` / ``scaffolds`` /
    ``result``), with ``result`` carrying ``called``, ``status``, ``interpretation``, ``carrier_contig``,
    ``variant`` (``repeat`` in array-unit space, plus ``label``), ``alleles_lengths`` and the coverage
    annotations the rest of the pipeline expects.
    """
    # 1. allele lengths (validated flank caller) + the picker note that drives the zygosity flag.
    detail = {}
    cassette_note = None
    from_cassette = False          # True when the length came from the cassette fallback (foreign amplicon)
    arb_kind = arb_platform = None
    arb_platform_source = None      # how platform was decided: "read-name" (trusted) vs "bq-fallback" (weak)
    arb_low_conf = False           # cassette arbiter's long_low_conf: candidate allele under-powered
    if lengths is None:
        lengths, detail = AL.two_alleles(vntr_bam, vntr_ref, return_detail=True)
        lengths = lengths or []
        # 1b. CASSETTE FALLBACK — foreign amplicon (our LR-PCR flanks absent). The validated flank caller
        # (``two_alleles``) anchors on our AL/AH flanks and returns NOTHING on a different PCR design (e.g.
        # VNTRtools' short-flank amplicon); the amplicon-agnostic cassette arbiter (motif 1->9,
        # ``vntr_raw_length.auto_length``) still measures the allele length there. STRICT fallback: it fires
        # ONLY when the flank caller yielded no length, so it is a NO-OP on our cohort (whose reads always
        # carry the flanks -> two_alleles returns lengths). The frameshift scan below then runs on the
        # cassette-derived contigs; ``FS.call`` anchors on the REFERENCE contig's structure, not on the
        # reads' flanks, so it works on a foreign design once the right contigs are supplied. The arbiter's
        # own low-confidence long allele (``long_low_conf``) is passed through as-is — the scan's depth
        # handling routes a thin minor allele to NEEDS_IGV/NEG, not to a false CALLED.
        if not lengths:
            try:
                import vntr_raw_length as _VRL
                try:
                    from .config import FLANK_TO_PHYSICAL_OFFSET as _OFF
                except Exception:
                    _OFF = 4
                # offset only affects the flank method; the cassette method (the foreign-amplicon path
                # this fallback targets) ignores it. Passed for parity with the rest of the pipeline.
                _arb = _VRL.auto_length(vntr_bam, chrom=None, offset=_OFF)
                if _arb and _arb.get("available") and _arb.get("alleles"):
                    lengths = [int(L) for L in _arb["alleles"]]
                    from_cassette = True
                    arb_kind = _arb.get("kind")
                    arb_platform = _arb.get("platform")
                    arb_platform_source = _arb.get("platform_source")
                    arb_low_conf = bool(_arb.get("long_low_conf"))
                    cassette_note = ("cassette-fallback "
                                     f"kind={_arb.get('kind')} platform={_arb.get('platform')} "
                                     f"flank_bp={_arb.get('flank_bp')} alleles={_arb.get('alleles')} "
                                     f"counts={_arb.get('counts')} long_low_conf={_arb.get('long_low_conf')}")
            except Exception as _e:                       # never let the fallback sink the call
                cassette_note = f"cassette-fallback-error {type(_e).__name__}: {_e}"
    note = (detail.get("note") or "")
    if cassette_note:
        note = (note + " | " + cassette_note).strip(" |")

    # ── FOREIGN amplicon + PacBio → CONSENSUS-scaffold path (not the ONT statistical detector) ──
    # STRICT gate: cassette fired AND foreign short-PCR AND PacBio HiFi CONFIRMED BY THE READ NAME (a `/ccs`
    # vote). A weak bq-fallback platform is NOT trusted — a processed ONT BAM can read median Q >= 30 and be
    # mis-called PacBio; with stripped read names we DEFAULT to ONT and keep the classic FS.call thresholds,
    # warning loudly. Our native ONT LR-PCR (flanks present -> fallback never fires) is NEVER routed here, so
    # the validated 100 %/100 % cohort performance is untouched.
    _pacbio_confirmed = (arb_platform == "PacBio HiFi" and arb_platform_source == "read-name")
    if lengths and from_cassette and arb_kind == "foreign/short-PCR" and _pacbio_confirmed:
        return _consensus_frameshift_result(vntr_bam, vntr_ref, lengths, note,
                                            min_mapq=min_mapq, low_conf=arb_low_conf)
    if from_cassette and arb_platform_source != "read-name":
        # Names uninformative (stripped) -> platform came from the weak bq-fallback (or none). We do NOT
        # trust it; a confidently-named ONT foreign amplicon (source == read-name, platform ONT) is silent.
        print("[dispatch_frameshift_vntr] foreign amplicon (no AL/AH anchors) and read names are "
              f"uninformative (platform={arb_platform!r} via source={arb_platform_source!r}) → DEFAULTING "
              "to ONT and the classic FS.call thresholds. If these ARE PacBio HiFi reads, supply read names "
              "ending in '/ccs' so the consensus (>50 %) path is used.", file=sys.stderr)

    # UNCERTAIN zygosity: a single length was returned because the 2nd peak was rejected ON CONTRAST
    # (a "tail bump" that might be a real allele — JOT), as opposed to a clean single peak or a
    # consensus-verified Delta<=1 fusion. Detected exactly as the validation bench does.
    zygosity_uncertain = (len(lengths) == 1) and ("contrast" in note.lower())

    contigs = [f"MUC1_VNTR_{int(L)}repeats" for L in lengths]

    # 2. strongest verdict across the allele contig(s) — VERBATIM eval recipe.
    best = _scan_lengths(vntr_bam, vntr_ref, lengths, min_mapq,
                         zygosity_uncertain=zygosity_uncertain,
                         force_pooled=from_cassette) if lengths else None

    # 3. remap into the dispatch_dupc_vntr shape.
    # EXACTLY one allele examined (length-homozygous, or the 2nd allele unresolved). Zero lengths is a
    # different case — nothing was measurable — handled as INDETERMINATE below, not as "single allele".
    single_allele = (len(lengths) == 1)
    if best:
        verdict = best["verdict"]
        called = verdict in POSITIVE_VERDICTS
        interp = _INTERP.get(verdict, "NEGATIVE")
        status = ("positive" if called
                  else "needs_igv" if verdict.startswith("NEEDS_IGV") else "negative")
        carrier_contig = best["contig"]
        # Read names carrying the variant — so caller.py can phase a size-homozygous heterozygous carrier
        # BY THE VARIANT (no SNP separates the alleles) and emit two haplotypes (one mutant, one wild-type).
        carrier_reads = _carrier_reads(vntr_bam, vntr_ref, best["contig"], best["repeat"],
                                       best["net_bp"], min_mapq=min_mapq)
        variant = {"repeat": best["repeat"], "repeat_index_space": "array_unit",
                   "label": best["variant"], "net_bp": best["net_bp"], "kind": best["kind"],
                   "vaf": best["vaf"], "vaf_mode": best["vaf_mode"], "raw_frac": best["raw_frac"],
                   "verdict": verdict, "n_carrier_reads": len(carrier_reads),
                   "carrier_reads": carrier_reads}
    else:
        # No call anywhere. INDETERMINATE only when NO length could be measured (nothing was scanned);
        # otherwise the contigs were scanned and came back clean -> a real NEGATIVE.
        verdict, called = "NEG", False
        interp = "INDETERMINATE" if not lengths else "NEGATIVE"
        status = "indeterminate" if not lengths else "negative"
        carrier_contig, variant = None, None

    return {
        "regime": "vntr",
        "caller": "frameshift_vntr(per-allele)",
        "region": None,
        "scaffolds": contigs,
        "result": {
            "assessed": bool(contigs),
            "called": called,
            "status": status,
            "interpretation": interp,
            "verdict": verdict,
            "method": "frameshift_vntr",
            "called_by": "frameshift_vntr" if best else None,
            "carrier_contig": carrier_contig,
            "variant": variant,
            "alleles_lengths": [int(L) for L in lengths],
            "alleles_scaffolded": len(contigs),
            # A single scaffold is legitimate for a length-homozygote, so this reports rather than judges.
            "single_allele_examined": single_allele,
            "coverage_note": ("only ONE allele length was resolved — either the subject is "
                              "length-homozygous or the second allele was not measured; a variant on an "
                              "unexamined allele cannot be excluded" if single_allele else None),
            # This detector routes low depth to POSITIVE_LOW_DEPTH / NEEDS_IGV rather than an
            # "underpowered" list, and its own VNTR-border + background guards stand in for the
            # run-length caller's scaffold-replication annotation. Kept present (empty / None) so a
            # consumer written for dispatch_dupc_vntr never KeyErrors.
            "underpowered_alleles": [],
            "scaffold_replication": None,
            "length_detail": {"note": note, "arbitration": detail.get("arbitration")},
            "reason": None if contigs else "no allele length could be measured",
        },
    }


# Optional CLI: a one-sample smoke test of the adapter alone (no prepare / caller / score).
def main(argv=None) -> int:
    import argparse
    import json
    ap = argparse.ArgumentParser(description="Frameshift VNTR verdict (validated chain), dupc_dispatch shape.")
    ap.add_argument("-b", "--bam", required=True, help="VNTR-aligned BAM (MUC1 profile)")
    ap.add_argument("-r", "--ref", required=True, help="multi-contig VNTR reference (...Kirby.fa)")
    ap.add_argument("--lengths", default=None, help="comma-separated allele lengths (else measured here)")
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)
    lengths = [int(x) for x in a.lengths.split(",")] if a.lengths else None
    dv = dispatch_frameshift_vntr(a.bam, a.ref, lengths=lengths, min_mapq=a.min_mapq)
    out = json.dumps(dv, indent=2, default=str)
    if a.json:
        with open(a.json, "w") as fh:
            fh.write(out)
    _r = dv["result"]
    print(out)
    return 0 if _r["assessed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
