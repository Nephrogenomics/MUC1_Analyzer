#!/usr/bin/env python3
"""run_pcr — the PCR-FASTQ clinical pipeline: prepare (MUC1 profile) -> frameshift detector -> call -> score.

Same overall shape as ``python -m muc1_analyzer run`` on an amplicon, but the detection step uses the
VALIDATED frameshift chain (``dispatch_frameshift_vntr``: allele_lengths.two_alleles -> frameshift_vntr.call,
Se/Sp = 1.000 on the 171-sample cohort) instead of ``runlen_shift``. The naming, scoring and report of the
existing pipeline are REUSED UNCHANGED — the variant NAME is produced authoritatively by ``caller`` from the
reconstructed consensus (motif dictionary), not by the detector. ``run`` and its detector are NOT modified.

This is Phase 1 of the integration: ``run_pcr`` is a standalone entry, callable directly, so it can be
validated without touching the production routing. Phase 2 (only after validation) makes the ``run`` fork
delegate a PCR FASTQ to ``run_pcr``, symmetrically to how it delegates adaptive-sampling to
``MUC1_Analyzer_fromfastq``.

Usage
-----
    python -m muc1_analyzer.run_pcr -i reads.fastq.gz -r MUC1_fakedVNTR1to150revcomplKirby.fa -s SAMPLE -o out/

    # already VNTR-aligned BAM (e.g. an existing xavier_bams/*.bam) -> prepare bypasses alignment,
    # which isolates the DETECTOR-in-pipeline from the alignment step during validation:
    python -m muc1_analyzer.run_pcr -i SAMPLE.xavier.bam -r MUC1_fakedVNTR1to150revcomplKirby.fa -s SAMPLE -o out/
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys


_READ_EXTS = {"gz", "bz2", "fastq", "fq", "fasta", "fa", "bam", "cram", "ubam", "sam"}


def _sample_from_input(path):
    """Derive a sample name from the INPUT filename when ``-s`` is not given, so outputs are not all
    named ``sample.*`` (a real sample-mix-up risk in a clinical run). Strips trailing read extensions
    (``.fastq.gz``, ``.fq``, ``.bam`` …); a non-extension token such as ``.muc1`` is kept, so it stays
    unique per sample. Falls back to ``"sample"`` if nothing usable remains."""
    base = os.path.basename(str(path))
    parts = base.split(".")
    while len(parts) > 1 and parts[-1].lower() in _READ_EXTS:
        parts.pop()
    name = ".".join(parts).strip()
    return name or "sample"


def _cmd_run_pcr(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="muc1_analyzer run_pcr",
        description="PCR-FASTQ pipeline: prepare (MUC1 profile) -> frameshift detector (validated) -> call "
                    "(consensus + naming) -> two-axis MUC1_Score. Standalone; does not modify `run`.")
    ap.add_argument("-i", "--input", required=True,
                    help="LR-PCR amplicon FASTQ(.gz) / uBAM, or an already VNTR-aligned BAM (prepare bypasses)")
    ap.add_argument("-r", "--vntr-ref", dest="vntr_ref", required=True,
                    help="multi-contig VNTR reference FASTA (…Kirby.fa)")
    ap.add_argument("-o", "--outdir", default="muc1_out", help="output directory")
    ap.add_argument("-s", "--sample", default=None,
                    help="sample name for the output files; default: derived from the input filename "
                         "(so outputs are not all named sample.*)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--align-profile", dest="align_profile", default="muc1",
                    choices=["setB", "muc1", "muc1_scoring", "xavier"],
                    help="LR-PCR minimap2 parameter set for `prepare`. Default `muc1` — the profile the "
                         "frameshift detector was VALIDATED on. Change only to price an "
                         "alternative against the clinical truth table.")
    ap.add_argument("--min-mapq", type=int, default=0,
                    help="min mapping quality for the frameshift detector (validated at 0)")
    ap.add_argument("--pcr-cap", type=int, default=None,
                    help="subsample a very deep amplicon to this many reads before alignment (prevents the "
                         "multi-contig minimap2 stall; preserves the allele ratio)")
    ap.add_argument("--mut-variant", default=None,
                    help="validated mutation name for the ONSET axis (e.g. 59dupC) if the caller's top-N missed it")
    ap.add_argument("--rs4072037-mut", dest="rs4072037_mut", choices=["C", "T"], default=None,
                    help="observed rs4072037 base of the MUTANT allele (severity axis); otherwise genotyped "
                         "T2T-natively on the VNTR-ref BAM")
    args = ap.parse_args(argv)
    if not args.sample:
        args.sample = _sample_from_input(args.input)

    os.makedirs(args.outdir, exist_ok=True)
    base = os.path.join(args.outdir, args.sample)

    # ── 1. prepare — any input -> a BAM aligned on the VNTR reference (MUC1 PCR profile) ────────────
    # An input that is ALREADY on the VNTR target is returned as-is (no realignment), so passing an
    # existing xavier_bams/*.bam validates the DETECTOR-in-pipeline in isolation from the alignment.
    from .prepare import prepare
    print("[run_pcr] step 1/4 — prepare (align on VNTR ref)", file=sys.stderr)
    bam = prepare(args.input, args.vntr_ref, base + ".vntr.bam",
                  threads=args.threads, pcr=True, no_pcr_autodetect=True,
                  recover_unmapped=False, sample=args.sample, workdir=args.outdir,
                  pcr_cap=args.pcr_cap, align_profile=args.align_profile)

    # ── 2. detect — validated frameshift chain, mapped to the dispatch_dupc_vntr shape ───────────────
    print("[run_pcr] step 2/4 — frameshift detection (validated chain)", file=sys.stderr)
    from .dispatch_frameshift_vntr import dispatch_frameshift_vntr
    dv = dispatch_frameshift_vntr(bam, args.vntr_ref, min_mapq=args.min_mapq)
    dupc_json = base + ".dupc.json"
    with open(dupc_json, "w") as f:
        json.dump(dv, f, indent=2, default=str)
    _r = dv.get("result", {}) or {}
    lengths = _r.get("alleles_lengths") or []
    _v = _r.get("variant") or {}
    print(f"[run_pcr] alleles={lengths}  verdict={_r.get('verdict')} ({_r.get('interpretation')})  "
          f"carrier={_r.get('carrier_contig')}  variant={_v.get('label')}@rep{_v.get('repeat')}  "
          f"vaf={_v.get('vaf')}", file=sys.stderr)
    if _r.get("single_allele_examined"):
        print("[run_pcr] ⚠ only ONE allele length resolved — a variant on an unexamined allele cannot be "
              "excluded", file=sys.stderr)

    # ── 3. call — consensus reconstruction + motif NAMING (authoritative variant name) + PDF/FASTA/JSON ─
    # The name shown for the carrier is `caller`'s consensus motif (KNOWN_REPEATS), NOT the detector's
    # label. `--dupc-json` carries the verdict + carrier repeat so the congruence layer can cross-check.
    # NB: no `--pcr` here on purpose — that caps the caller at max-depth 200, and a subsampled consensus
    # is the one place delCC / insG / del8_27 stop being nameable. Alignment already used the PCR profile.
    print("[run_pcr] step 3/4 — call (consensus + naming)", file=sys.stderr)
    from . import caller
    call_json = base + ".analyzer.json"
    consensus_fa = base + ".consensus.fa"
    call_argv = ["-b", bam, "-r", args.vntr_ref, "-s", args.sample,
                 "--min-mq", "0", "--threads", str(args.threads),
                 "--json", call_json, "--pdf", base + ".pdf", "--dupc-json", dupc_json,
                 "--fasta-consensus", consensus_fa]
    # Reconstruct the haplotypes AT the detector's measured lengths, so the long (PCR-depleted) allele —
    # where delCC / del8_27 / insG sit — actually gets a consensus.
    if lengths:
        call_argv += ["--alleles", ",".join(str(int(x)) for x in lengths)]
    if args.mut_variant:
        call_argv += ["--mut-variant", args.mut_variant]
    rc = caller.main(call_argv)
    if rc:
        print("[run_pcr] the caller failed — stopping before score", file=sys.stderr)
        return rc or 1

    # Index the consensus FASTA (one sequence per resolved haplotype length — a single one on a
    # length-homozygote, two on a length-heterozygote). `--fasta-consensus` writes the FASTA
    # but not its .fai; downstream tools (and IGV) want the index.
    try:
        import pysam
        if os.path.exists(consensus_fa) and os.path.getsize(consensus_fa) > 0:
            pysam.faidx(consensus_fa)
            print(f"[run_pcr] consensus FASTA → {consensus_fa} (+ .fai)", file=sys.stderr)
    except Exception as e:
        print(f"[run_pcr] consensus FASTA index skipped ({e})", file=sys.stderr)

    # ── 3b. ribbon — allele-reconstruction figure as a DEFAULT output ────────────────────────────────
    # The motif-array picture of the two reconstructed haplotypes (variant highlighted only on a POSITIVE
    # verdict, per the clinical rule). Never allowed to sink the run.
    from . import vntr_ribbon
    ribbon_png = base + ".ribbon.png"
    if vntr_ribbon.render_safe(call_json, ribbon_png, verdict=_r.get("verdict")):
        print(f"[run_pcr] ribbon figure → {ribbon_png}", file=sys.stderr)
    else:
        print("[run_pcr] ribbon figure skipped (no haplotypes or matplotlib unavailable)", file=sys.stderr)

    # ── 4. score — two-axis MUC1_Score ───────────────────────────────────────────────────────────────
    print("[run_pcr] step 4/4 — score", file=sys.stderr)
    score_argv = ["--analyzer-json", call_json, "-s", args.sample, "-o", base + ".score.json"]
    if lengths:
        score_argv += ["--arbiter-alleles", ",".join(str(int(x)) for x in lengths)]
    # Which allele is the MUTANT comes from the detector's carrier contig — a measurement — rather than
    # from the caller's long/short ordering, which is a deduction.
    _cc = re.search(r"MUC1_VNTR_(\d+)repeats", str(_r.get("carrier_contig") or ""))
    if _cc:
        score_argv += ["--carrier-len", _cc.group(1)]
    # rs4072037 severity axis: T2T-native on the VNTR-ref BAM (fixed 5' flank offset) — no GRCh38 needed
    # from a FASTQ. A no-call simply leaves the severity axis unset.
    if args.rs4072037_mut:
        score_argv += ["--rs4072037-mut", args.rs4072037_mut]
    else:
        try:
            from .detectors.splice_snp import genotype_rs4072037_vntr_ref
            gt = genotype_rs4072037_vntr_ref(bam, args.vntr_ref)
            dosage = gt.get("dosage") if gt else None
            if dosage is not None:
                score_argv += ["--snp-genotype", str(dosage)]
                print(f"[run_pcr] rs4072037: {gt.get('genotype') or gt.get('gt_label')} "
                      f"depth={gt.get('depth')}", file=sys.stderr)
            else:
                print("[run_pcr] rs4072037 not callable here → severity axis unset "
                      "(pass --rs4072037-mut to set it)", file=sys.stderr)
        except Exception as e:
            print(f"[run_pcr] rs4072037 genotyping skipped ({e})", file=sys.stderr)

    from .cli import main as score_main
    rc = score_main(score_argv)

    # ── report: re-render the PDF now that the score exists (same single renderer, one extra second) ──
    if rc == 0:
        try:
            res = json.load(open(base + ".score.json"))
            print(f"[run_pcr] {args.sample}: carrier={res.get('carrier')} "
                  f"frameshift_tail={res.get('fs_tail')} → {res.get('fs_category')} "
                  f"(rs4072037 {res.get('splice_base')}, info only)", file=sys.stderr)
        except Exception:
            res = None
        if res:
            try:
                caller.render_pdf_from_json(call_json, base + ".pdf", score=res)
                print(f"[run_pcr] report re-rendered with the MUC1 Score → {base}.pdf", file=sys.stderr)
            except Exception as e:
                print(f"[run_pcr] score block not added to the PDF ({e}) — the report from `call` stands",
                      file=sys.stderr)
        print(f"\n[run_pcr] done → {base}.score.json  (calls: {call_json}, dupc: {dupc_json}, reads: {bam})",
              file=sys.stderr)
    return rc


def main(argv=None) -> int:
    return _cmd_run_pcr(argv)


if __name__ == "__main__":
    raise SystemExit(main())
