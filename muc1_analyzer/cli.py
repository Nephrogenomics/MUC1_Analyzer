"""MUC1_Score orchestrator CLI.

Phase 1 (`--cram`/`--bam`, or `--fastq` aligned on the fly): automatic detectors from the
patient alignment → DEL (MUC1-SV8533) + rs4072037 + phased VNTR + dupC gate → MUC1_Score, with
**contig-naming auto-detection** (`chr1` vs `1`). Outputs: `--table` (TSV) and `--report` (JSON).

Phase 0 (`--analyzer-json`): assembles a score from already-available fragments (a MUC1_Analyzer
JSON for the VNTR/dupC + DEL/SNP passed explicitly).

Usage (Phase 1):
    python -m muc1_analyzer --cram patient.cram --ref GRCh38.fa --vntr-ref vntr_ref.fa \
        --pon pon_dupc.json -s PATIENT --table PATIENT.tsv --report PATIENT.json

Usage (Phase 0):
    python -m muc1_analyzer --analyzer-json PATIENT.json --del-genotype 1 --del-phase cis_mut \
        --snp-genotype 1 -s PATIENT -o PATIENT.score.json
"""
from __future__ import annotations
import argparse
import json
import os
import subprocess
import sys

from .detectors import vntr
# NB: `score.py` (LEGACY composite, DEL) is imported LAZILY, only in the legacy path (build_params /
# --legacy-composite). The canonical TWO-AXIS path (clinical_call / analyzer-JSON) does NOT touch it, so the
# public build can drop score.py + regulatory_del + batch_score_phased without breaking the two-axis CLI.

# chr1 lengths (GRCh38 / GRCh37) — to recover the contig by size when the name is atypical.
_CHR1_LENGTHS = (248_956_422, 249_250_621)


def pick_chr1_contig(sqs) -> str | None:
    """Pick the chr1 contig name among a BAM/CRAM's @SQ entries (list of (name, length)).
    Prefers the explicit name `chr1` then `1`; falls back to the chr1 LENGTH. Pure (testable)."""
    names = {n for n, _ in sqs}
    for cand in ("chr1", "1", "CM000663.2", "NC_000001.11"):
        if cand in names:
            return cand
    for n, ln in sqs:
        if ln in _CHR1_LENGTHS:
            return n
    return None


def _detect_chrom(aln: str, ref: str | None) -> str | None:
    import pysam
    mode = "rc" if str(aln).endswith(".cram") else "rb"
    kw = {"reference_filename": ref} if (mode == "rc" and ref) else {}
    with pysam.AlignmentFile(aln, mode, **kw) as af:
        return pick_chr1_contig(list(zip(af.references, af.lengths)))


def minimap_preset(pacbio: bool = False) -> str:
    """minimap2 preset for raw-read alignment: `map-hifi` for PacBio HiFi, else `map-ont`.
    Pure (testable). MUC1_Analyzer is tech-agnostic — only the aligner preset differs by platform
    (cf. docs/MUC1_log.md 2026-07-13: analyzer ran directly on PacBio HiFi, structure decoded)."""
    return "map-hifi" if pacbio else "map-ont"


def _align_fastq(fastq: str, ref: str, out_bam: str, threads: int = 8, pacbio: bool = False) -> str:
    """Align raw reads (ONT or PacBio HiFi) to `ref` with minimap2 → sorted+indexed BAM. NOT phased
    (no HP tags) → DEL/SNP/VNTR phasing will be unavailable; genotypes are still read."""
    preset = minimap_preset(pacbio)
    print(f"[fastq] minimap2 {preset} {fastq} → {out_bam} (NOT phased: limited phasing)", file=sys.stderr)
    p1 = subprocess.Popen(["minimap2", "-ax", preset, "-t", str(threads), ref, fastq],
                          stdout=subprocess.PIPE)
    subprocess.run(["samtools", "sort", "-@", str(threads), "-o", out_bam, "-"],
                   stdin=p1.stdout, check=True)
    p1.wait()
    subprocess.run(["samtools", "index", out_bam], check=True)
    return out_bam


# ── Phase 0 (explicit fragments) ───────────────────────────────────────────────
def build_params(args):
    from .score import MUC1Params                        # LEGACY composite only (lazy)
    frag = vntr.from_analyzer_json(args.analyzer_json) if args.analyzer_json else {}
    mutation = frag.get("mutation_present")
    if args.mutation is not None:
        mutation = (args.mutation == "yes")
    mv = frag.get("mutation_variant") or {}
    mut_pos = args.mut_position if getattr(args, "mut_position", None) is not None else mv.get("repeat_index")
    return MUC1Params(
        mutation_present=mutation,
        vntr_len_mut=frag.get("vntr_len_mut"), vntr_len_healthy=frag.get("vntr_len_healthy"),
        mut_position=mut_pos,                    # -> combined prognostic component (falls back to ratio if None)
        del_genotype=args.del_genotype, del_phase=args.del_phase,
        snp_genotype=args.snp_genotype, snp_phase=args.snp_phase,
        notes=frag.get("notes", {}))


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="muc1_score",
        description="MUC1_Score — composite ADTKD-MUC1 score from a patient CRAM/BAM (Phase 1) "
                    "or from fragments (Phase 0). See docs/ROADMAP.")
    p.add_argument("-s", "--sample", default="", help="sample name")
    # ── Phase 1: alignment input ──
    p.add_argument("--cram", "--bam", dest="aln", help="aligned patient CRAM/BAM (ideally haplotagged)")
    p.add_argument("--fastq", help="raw reads (ONT or PacBio HiFi) → aligned on the fly to --ref (NOT-phased fallback)")
    p.add_argument("--pacbio", action="store_true",
                   help="PacBio HiFi reads: align --fastq with minimap2 map-hifi (default map-ont, ONT)")
    p.add_argument("--ref", "--genome-ref", dest="ref", help="genome FASTA (GRCh38; required for CRAM and fetch)")
    p.add_argument("--genome", choices=["hg38", "t2t"], default="hg38",
                   help="coordinate system of the input alignment (default hg38; DEL/rs4072037 = GRCh38)")
    p.add_argument("--vntr-ref", dest="vntr_ref", help="multi-contig VNTR reference (VNTR/dupC; else baseline DEL+SNP score)")
    p.add_argument("--pon", help="dupC PoN table (pon_dupc.json)")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--table", help="[Phase 1] TSV output (one line: score + components)")
    p.add_argument("--report", help="[Phase 1] full JSON output")
    # ── Phase 0: explicit fragments ──
    p.add_argument("-o", "--output", default="-", help="[Phase 0] score JSON (- = stdout)")
    p.add_argument("--analyzer-json", help="[Phase 0] MUC1_Analyzer JSON (VNTR/dupC)")
    p.add_argument("--arbiter-alleles", default=None,
                   help="authoritative allele lengths from the alignment-free arbiter, e.g. '45,77'. They "
                        "REPLACE the caller's haplotype lengths in mut_len/healthy_len/onset_index, which "
                        "the caller can get wrong by losing a PCR-depleted allele.")
    p.add_argument("--carrier-len", type=int, default=None,
                   help="length of the allele the dupC caller found the variant on — a MEASUREMENT of "
                        "which arbiter allele is the mutant. Without it the mutant is inferred from the "
                        "caller's long/short ordering and the report says so.")
    p.add_argument("--mut-position", type=int, default=None,
                   help="[Phase 0] position (repeat) of the mutation in the array -> combined "
                        "prognostic component ratio+position (default: repeat_index from the JSON, else ratio only)")
    p.add_argument("--del-genotype", type=int, choices=[0, 1, 2], default=None,
                   help="[LEGACY — mechanistic, weight 0, NOT in the two-axis score] short-VNTR-array dosage")
    p.add_argument("--del-phase", choices=["cis_mut", "trans", "cis_healthy"], default=None,
                   help="[LEGACY — mechanistic, NOT prognostic]")
    p.add_argument("--snp-genotype", type=int, choices=[0, 1, 2], default=None,
                   help="rs4072037 dosage (0/1/2). In the two-axis score it drives the SEVERITY axis "
                        "(with --clinical-call); in the legacy composite it is mechanistic (weight 0).")
    p.add_argument("--snp-phase", choices=["cis_mut", "trans", "cis_healthy"], default=None)
    p.add_argument("--mutation", choices=["yes", "no"], default=None,
                   help="force dupC presence (otherwise derived from the analyzer JSON)")
    # ── canonical TWO-AXIS score (public) : onset ← frameshift position, severity ← rs4072037 splice ──
    p.add_argument("--clinical-call", dest="clinical_call", default=None,
                   help="[TWO-AXIS] visually/clinically finalized call, e.g. "
                        "'64 repeats | 80 repeats (del8_27 @ repeat 37)' → canonical two-axis MUC1_Score "
                        "(no DEL). For variants finalized by eye (del8_27) or coverage-floored het long-VNTR.")
    p.add_argument("--rs4072037-mut", dest="rs4072037_mut", choices=["C", "T"], default=None,
                   help="[TWO-AXIS] observed rs4072037 base of the MUTANT allele (visual review / dRNA-seq); "
                        "highest precedence for the severity axis, else --snp-genotype dosage + cis rule.")
    p.add_argument("--w-onset", type=float, default=2.0,
                   help="[TWO-AXIS] onset-axis weight (direction-of-effect; UNCALIBRATED, n small)")
    p.add_argument("--w-splice", type=float, default=1.0,
                   help="[TWO-AXIS] severity-axis weight (direction-of-effect; UNCALIBRATED)")
    p.add_argument("--legacy-composite", action="store_true",
                   help="[dev] use the OLD additive ratio composite (score.py) instead of the canonical "
                        "two-axis (DEL/rs4072037 weight 0 there). For reproducibility only.")
    return p.parse_args(argv)


def _run_clinical_call(args) -> int:
    """Canonical TWO-AXIS MUC1_Score from a visual/clinical call string (no DEL, no methylation)."""
    from .clinical_call import score_call
    res = score_call(args.clinical_call, rs4072037_mut=args.rs4072037_mut,
                     dosage=args.snp_genotype, w_onset=args.w_onset, w_splice=args.w_splice)
    res["sample"] = args.sample
    res["score_model"] = "two-axis (onset = frameshift position ; severity = rs4072037 splice of the mutant allele)"
    out = json.dumps(res, indent=2, ensure_ascii=False)
    if args.output == "-":
        print(out)
    else:
        with open(args.output, "w") as fh:
            fh.write(out + "\n")
        print(f"[OK] two-axis score → {args.output}", file=sys.stderr)
    if res.get("carrier"):
        print(f"[MUC1_Score two-axis] {args.sample or '?'}: "
              f"ONSET onset_index={res.get('onset_index')} (score {res.get('onset_score')}) · "
              f"SEVERITY splice={res.get('splice_base')} (score {res.get('severity_score')}, "
              f"{res.get('splice_source')}) · mut={res.get('mut_len')}/healthy={res.get('healthy_len')} "
              f"ratio={res.get('ratio')}", file=sys.stderr)
    else:
        print(f"[MUC1_Score two-axis] {args.sample or '?'}: non-carrier "
              f"({res.get('short')}/{res.get('long')})", file=sys.stderr)
    return 0


def _run_phase1(args) -> int:
    aln = args.aln
    if not aln and args.fastq:
        if not args.ref:
            print("[ERROR] --ref required to align --fastq", file=sys.stderr); return 2
        aln = _align_fastq(args.fastq, args.ref, (args.sample or "sample") + ".mm2.bam",
                           args.threads, pacbio=args.pacbio)
    if not aln:
        print("[ERROR] provide --cram/--bam or --fastq", file=sys.stderr); return 2
    if str(aln).endswith(".cram") and not args.ref:
        print("[ERROR] --ref (GRCh38 FASTA) required for a CRAM", file=sys.stderr); return 2
    try:                                                    # one-shot detect+score = dev/advanced (legacy engine)
        from .batch_score_phased import score_one, _COLS
    except ImportError:
        print("[ERROR] the --cram/--bam/--fastq path (one-shot detection+score) is not included in this "
              "build. Use `MUC1_Analyzer.py` (detection) then `muc1_score --analyzer-json`, or "
              "`muc1_score --clinical-call`.", file=sys.stderr); return 2
    if args.genome == "t2t":
        print("[WARNING] --genome t2t: DEL MUC1-SV8533 & rs4072037 are annotated on GRCh38; "
              "the input must be GRCh38-aligned for these components.", file=sys.stderr)

    chrom = _detect_chrom(aln, args.ref)
    print(f"[INFO] chr1 contig detected: {chrom or '(config default)'}", file=sys.stderr)
    pon = json.load(open(args.pon)) if args.pon else None
    if not args.vntr_ref:
        print("[INFO] no --vntr-ref → baseline score (DEL + rs4072037; VNTR/dupC ignored)", file=sys.stderr)

    row = score_one(aln, args.ref, args.vntr_ref, pon=pon, threads=args.threads, chrom=chrom)
    if args.sample:
        row["sample"] = args.sample

    # canonical TWO-AXIS, attached to the detection row (onset ← position, severity ← rs4072037 splice)
    from .clinical_call import score_from_fields
    ta = score_from_fields(row.get("vntr_len_mut"), row.get("vntr_len_healthy"), row.get("mut_position"),
                           rs4072037_mut=args.rs4072037_mut, dosage=row.get("snp_genotype"),
                           w_onset=args.w_onset, w_splice=args.w_splice)
    row.update({k: ta[k] for k in ("mut_is_long", "ratio", "onset_index", "onset_score",
                                   "splice_base", "severity_score", "splice_source")})
    row["score_model"] = "two-axis"

    if args.table:
        with open(args.table, "w") as fh:
            fh.write("\t".join(_COLS) + "\n")
            fh.write("\t".join("" if row.get(c) is None else str(row.get(c, "")) for c in _COLS) + "\n")
        print(f"[OK] table → {args.table}", file=sys.stderr)
    if args.report:
        with open(args.report, "w") as fh:
            json.dump(row, fh, indent=2, ensure_ascii=False)
        print(f"[OK] report → {args.report}", file=sys.stderr)
    if not args.table and not args.report:
        print(json.dumps(row, indent=2, ensure_ascii=False))
    print(f"[MUC1_Score two-axis] {row.get('sample')}: ONSET onset_index={row.get('onset_index')} "
          f"(score {row.get('onset_score')}) · SEVERITY splice={row.get('splice_base')} "
          f"(score {row.get('severity_score')}, {row.get('splice_source')}) · "
          f"ratio={row.get('ratio')} · mutation={row.get('mutation_present')}", file=sys.stderr)
    return 0


def main(argv=None):
    args = parse_args(argv)
    if args.clinical_call:                                  # canonical TWO-AXIS (manual/visual call)
        return _run_clinical_call(args)
    if args.aln or args.fastq:                              # Phase 1 : from the alignment
        return _run_phase1(args)

    # Phase 0 : explicit fragments
    if not args.analyzer_json and args.del_genotype is None and args.snp_genotype is None:
        print("[ERROR] nothing to score: provide --cram/--bam/--fastq (Phase 1) or "
              "--analyzer-json/--del-*/--snp-* (Phase 0). See --help.", file=sys.stderr)
        return 2

    # canonical TWO-AXIS from the analyzer JSON (default); legacy composite only on request / no JSON
    if args.analyzer_json and not args.legacy_composite:
        from .clinical_call import apply_arbiter_lengths, score_from_fields
        frag = vntr.from_analyzer_json(args.analyzer_json)
        # The ARBITER is authoritative on length; the caller's haplotype lengths are a by-product of
        # consensus reconstruction and can lose a depleted allele entirely (measured: arbiter 45/77,
        # caller 43/44). The variant POSITION still comes from the caller — that is what it is for.
        if getattr(args, "arbiter_alleles", None):
            _arb = [int(x) for x in str(args.arbiter_alleles).replace(",", " ").split()]
            frag = apply_arbiter_lengths(frag, _arb, carrier_len=getattr(args, "carrier_len", None))
            print(f"[score] lengths from the {frag.get('length_source')} "
                  f"→ mut={frag.get('vntr_len_mut')} healthy={frag.get('vntr_len_healthy')} "
                  f"(caller said {frag.get('caller_lengths')})", file=sys.stderr)
        mv = frag.get("mutation_variant") or {}
        pos = args.mut_position if getattr(args, "mut_position", None) is not None else mv.get("repeat_index")
        res = score_from_fields(frag.get("vntr_len_mut"), frag.get("vntr_len_healthy"), pos,
                                rs4072037_mut=args.rs4072037_mut, dosage=args.snp_genotype,
                                w_onset=args.w_onset, w_splice=args.w_splice)
        res.update(sample=args.sample, carrier=frag.get("mutation_present"), variant=mv.get("name"),
                   score_model="two-axis (onset = frameshift position ; severity = rs4072037 splice)")
        out = json.dumps(res, indent=2, ensure_ascii=False)
        (print(out) if args.output == "-" else open(args.output, "w").write(out + "\n"))
        print(f"[MUC1_Score two-axis] {args.sample or '?'}: onset_index={res.get('onset_index')} · "
              f"severity splice={res.get('splice_base')} ({res.get('splice_source')}) · "
              f"ratio={res.get('ratio')} (mut={res.get('mut_len')}/healthy={res.get('healthy_len')})",
              file=sys.stderr)
        return 0

    from .score import compute_score                     # LEGACY composite only (lazy)
    params = build_params(args)
    result = compute_score(params)
    result["sample"] = args.sample
    result["params"] = {k: v for k, v in params.__dict__.items() if k != "notes"}
    result["score_model"] = "legacy_composite (DEL/rs4072037 weight 0)"
    out = json.dumps(result, indent=2, ensure_ascii=False)
    if args.output == "-":
        print(out)
    else:
        with open(args.output, "w") as fh:
            fh.write(out + "\n")
        print(f"[OK] score → {args.output}  (legacy MUC1_Score = {result.get('score')})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
