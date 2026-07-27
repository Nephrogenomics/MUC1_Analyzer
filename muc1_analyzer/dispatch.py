"""`python -m muc1_analyzer <subcommand>` — the single entry point of the tool.

    run       full pipeline: prepare → call → score        (the beginner one-liner)
    prepare   any input → BAM aligned on the VNTR reference (extract + align)
    call      VNTR haplotype caller (length, motifs, frameshifts incl. dupC)
    score     two-axis MUC1_Score (onset ← frameshift position, severity ← rs4072037 splice)
    bundle    IGV review bundle

Sub-parsers are NOT declared here: each subcommand owns its own argparse and is imported
LAZILY, so `muc1_analyzer score --clinical-call` (pure JSON, no genomics deps) never pays
for pysam, and a public build that drops a research module still runs the rest.

Backward compatibility: a flat legacy invocation (`python -m muc1_analyzer --clinical-call …`,
as used by existing sbatch launchers on the cluster) is routed to `score` with a notice
rather than failing, since those launchers are pulled and run by hand on the cluster.
"""
from __future__ import annotations

import sys

SUBCOMMANDS = ("run", "prepare", "call", "score", "bundle")

USAGE = """usage: python -m muc1_analyzer <subcommand> [options]

MUC1_Analyzer — MUC1 VNTR analysis for ADTKD-MUC1 (long reads).

subcommands:
  run       full pipeline (prepare -> call -> score) — start here
  prepare   any input (fastq / uBAM / hg38 / T2T) -> BAM on the VNTR reference
  call      VNTR haplotype caller: length, motifs, frameshifts (incl. 59dupC)
  score     two-axis MUC1_Score (onset = frameshift position, severity = rs4072037 splice)
  bundle    IGV review bundle

  python -m muc1_analyzer <subcommand> --help   for the options of one subcommand

quickstart:
  python -m muc1_analyzer run -i reads.fastq.gz -r MUC1_fakedVNTR1to150revcomplKirby.fa \\
      -s SAMPLE -o out/
"""


def signature_from_reads(path: str, *, cap: int = 2000, maxmm: int = 9,
                         max_scan: int = 400_000) -> dict:
    """`amplicon_signature`'s measurements on RAW READS — fastq or uBAM, no alignment. Impure.

    `read_cassette_span` works on the sequence, so the whole fork can be decided before a single read is
    aligned. Returns {n, n_scanned, al_frac, ah_frac, tight5, tight3}.

    BOTH ends are bounded, and the second one matters: `cap` stops once enough MUC1-like reads are found,
    `max_scan` stops after that many reads have been LOOKED at. On an adaptive-sampling uBAM (20 GB, MUC1 a
    tiny fraction of it) a cap on found-reads alone would stream most of the file just to answer a routing
    question — the very cost this pre-alignment fork exists to avoid."""
    import statistics as st

    import pysam
    from vntr_raw_length import AL, AH, fuzzy_find, rc, read_cassette_span

    f5, f3, n, have_al, have_ah, scanned = [], [], 0, 0, 0, 0
    lengths = []

    def _consume(seq):
        nonlocal n, have_al, have_ah, scanned
        scanned += 1
        lengths.append(len(seq))
        # A read too short to hold flanks + a tandem cannot carry the signature, and an adaptive-sampling
        # file is mostly REJECTED fragments — testing them is the whole cost for none of the information.
        if len(seq) < 3000:
            return
        span = read_cassette_span(seq, maxmm)
        if span is None:
            return
        n += 1
        f5.append(span[1])
        f3.append(span[2])
        rcs = rc(seq)
        have_al += bool(fuzzy_find(AL, seq, maxmm) or fuzzy_find(AL, rcs, maxmm))
        have_ah += bool(fuzzy_find(AH, seq, maxmm) or fuzzy_find(AH, rcs, maxmm))

    p = str(path).lower()
    if p.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz", ".fasta", ".fa", ".fa.gz")):
        with pysam.FastxFile(path) as fh:
            for e in fh:
                _consume(e.sequence.upper())
                if n >= cap or scanned >= max_scan:
                    break
    else:
        with pysam.AlignmentFile(path, check_sq=False) as af:
            for r in af.fetch(until_eof=True):
                if r.is_secondary or r.is_supplementary or r.query_sequence is None:
                    continue
                _consume(r.query_sequence.upper())
                if n >= cap or scanned >= max_scan:
                    break
    # Read-LENGTH tightness, over EVERY read scanned — computable without finding a single MUC1 read.
    # An LR-PCR is a tight band around its product; adaptive sampling is dominated by short rejected
    # fragments and long on-target reads, i.e. broad. This is what answers "is this an amplicon" when the
    # cassette scan comes back empty, which on a 20 GB AS uBAM is the normal outcome.
    len_tight = 0.0
    if lengths:
        lm = st.median(lengths)
        len_tight = sum(1 for x in lengths if abs(x - lm) <= 0.25 * lm) / len(lengths)
    base = {"n_scanned": scanned, "len_tight": round(len_tight, 3),
            "median_len": int(st.median(lengths)) if lengths else 0}
    if not n:
        return {"n": 0, "al_frac": 0.0, "ah_frac": 0.0, "tight5": 0.0, "tight3": 0.0, **base}
    f5m, f3m = st.median(f5), st.median(f3)
    return {"n": n, "al_frac": have_al / n, "ah_frac": have_ah / n,
            "tight5": sum(1 for x in f5 if abs(x - f5m) <= 200) / n,
            "tight3": sum(1 for x in f3 if abs(x - f3m) <= 200) / n, **base}


def _delegate_fromfastq(args, route: str) -> int:
    """Hand an AS/WGS read set to `MUC1_Analyzer_fromfastq.py` end to end. Returns its rc, or None when
    the tool is unavailable (caller then falls back to our own pipeline).

    That script owns the AS terrain: a ruler-anchored streaming length (tens of millions of reads without
    OOM), the indel-burden scaffold sweep, read recovery, and SNP phasing for alleles the length cannot
    separate. Our pipeline owns the amplicon terrain — the frozen chaining flags and the cassette fallback.
    Neither reimplements the other; the fork picks the one whose assumptions the sample actually meets.

    The report is the same either way: its `write_merged_pdf` delegates to `caller.write_pdf_report`."""
    import importlib.util
    import os
    import sys as _sys
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script = os.path.join(here, "MUC1_Analyzer_fromfastq.py")
    if not os.path.exists(script):
        return None
    try:
        spec = importlib.util.spec_from_file_location("muc1_fromfastq", script)
        mod = importlib.util.module_from_spec(spec)
        _cwd = os.getcwd()
        os.chdir(here)                     # it imports `vntr_raw_length` from the repository root
        try:
            spec.loader.exec_module(mod)
        finally:
            os.chdir(_cwd)
    except Exception as e:
        print(f"[run] fromfastq unavailable ({e}) — staying on our pipeline", file=_sys.stderr)
        return None
    argv = ["-b", args.input, "-r", args.vntr_ref, "-s", args.sample, "-o", args.outdir,
            "-t", str(args.threads)]
    if args.pacbio:
        argv.append("--pacbio")
    if args.ref38:
        argv += ["--ref-cram", args.ref38]
    if args.rs4072037_mut:
        argv += ["--rs4072037", args.rs4072037_mut]
    print(f"[run] route {route} → delegating to MUC1_Analyzer_fromfastq "
          f"(ruler length + scaffold sweep + phasing)", file=_sys.stderr)
    return mod.main(argv)


def _looks_unaligned(path: str) -> bool:
    """True for a uBAM — a .bam whose header declares no reference. Impure, never raises.

    A uBAM carries the .bam extension but has no genomic locus to read a dupC from, so it must route like
    a fastq. Deciding on the extension alone sent it down the genomic path, where it silently found
    nothing."""
    try:
        import pysam
        with pysam.AlignmentFile(path, check_sq=False) as af:
            return not af.header.to_dict().get("SQ")
    except Exception:
        return False


def _cmd_run(argv) -> int:
    """prepare → call → score, chained on one sample."""
    import argparse
    import json
    import os
    import re

    ap = argparse.ArgumentParser(
        prog="muc1_analyzer run",
        description="Full pipeline on one sample: prepare (extract+align) → call (VNTR "
                    "haplotypes + frameshift) → score (two-axis MUC1_Score).")
    ap.add_argument("-i", "--input", required=True,
                    help="fastq(.gz) | uBAM | aligned BAM/CRAM (hg38 or T2T)")
    ap.add_argument("-r", "--vntr-ref", dest="vntr_ref", required=True,
                    help="multi-contig VNTR reference FASTA")
    ap.add_argument("-o", "--outdir", default="muc1_out", help="output directory")
    ap.add_argument("-s", "--sample", default="sample")
    ap.add_argument("--ref", dest="ref38", default=None, help="genome FASTA (required for a CRAM input)")
    ap.add_argument("--ont", action="store_true", help="Oxford Nanopore reads (default)")
    ap.add_argument("--pacbio", action="store_true", help="PacBio HiFi reads")
    ap.add_argument("--pcr", action="store_true",
                    help="LR-PCR amplicon (sets the caller PCR defaults: depth cap + frozen chaining flags)")
    ap.add_argument("--no-pcr-autodetect", action="store_true",
                    help="disable prepare's read-count LR-PCR auto-detection")
    ap.add_argument("--pcr-cap", type=int, default=None,
                    help="--pcr: subsample a deep amplicon to this many reads before minimap2 (prevents the "
                         "multi-contig stall; preserves the allele ratio). Default from config.")
    ap.add_argument("--recover-unmapped", action="store_true",
                    help="(DEFAULT for aligned inputs) recover the unplaced VNTR reads (thin-locus rescue, ~35x)")
    ap.add_argument("--no-recover-unmapped", action="store_true",
                    help="disable the default unplaced-read recovery on an aligned input")
    ap.add_argument("--copy-offset", type=int, default=None,
                    help="flank→physical copy offset (default 4 ALWAYS: the AL/AH anchors sit ~4 units inside the "
                         "array, so a flank measurement is 4 short — a property of the METHOD, not of the assay)")
    ap.add_argument("--snp-chrom", default="chr1", help="rs4072037 contig in the genomic input")
    ap.add_argument("--snp-pos", type=int, default=None, help="rs4072037 1-based position (default GRCh38)")
    ap.add_argument("--mut-variant", default=None,
                    help="validated mutation for the ONSET axis (e.g. 59dupC) when the top-N missed it")
    ap.add_argument("--mut-repeat", type=int, default=None, help="repeat position of --mut-variant")
    ap.add_argument("--mut-allele-len", type=int, default=None, help="copies of the mutant allele")
    ap.add_argument("--pon", default=None, help="dupC panel-of-normals JSON")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--input-type", default="auto",
                    choices=["auto", "pcr-ont", "pcr-pacbio", "as-ont", "targeted-pacbio",
                             "wgs-ont", "wgs-pacbio"],
                    help="what this read set IS. Declaring it beats any heuristic and also sets the "
                         "aligner preset: pcr-* runs our amplicon path (frozen chaining flags + cassette "
                         "fallback), every non-amplicon type is delegated to MUC1_Analyzer_fromfastq (no primers, "
                         "streaming search over the whole file). 'auto' infers it from the reads and "
                         "STOPS rather than guessing when it cannot.")
    ap.add_argument("--no-fromfastq", action="store_true",
                    help="never delegate an AS/WGS read set to MUC1_Analyzer_fromfastq — run everything "
                         "on this pipeline (the fork is otherwise decided from the read sequences)")
    ap.add_argument("--rs4072037-mut", dest="rs4072037_mut", choices=["C", "T"], default=None,
                    help="observed rs4072037 base of the MUTANT allele (severity axis); "
                         "otherwise genotyped from the alignment when possible")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    base = os.path.join(args.outdir, args.sample)

    # ── 0. FORK, decided on the raw reads — before a single one is aligned ────────────────────────
    # Aligning tens of millions of adaptive-sampling reads only to discover they belonged to the other
    # pipeline is exactly what this must avoid, so the decision is made from the read SEQUENCES.
    from .config import FLANK_TO_PHYSICAL_OFFSET
    from . import length_congruence as LCG
    _reads_only = not str(args.input).lower().endswith((".bam", ".cram")) or _looks_unaligned(args.input)
    _declared, _dec_pacbio = LCG.route_from_declared(args.input_type)
    if _declared:
        args.pacbio = args.pacbio or _dec_pacbio          # the declaration also sets the aligner preset
        print(f"[run] input declared as {args.input_type} → route {_declared}"
              f"{' (PacBio preset)' if _dec_pacbio else ''}", file=sys.stderr)
    if _reads_only and not args.pcr and not args.no_fromfastq:
        _r0 = _declared
        if _r0 is None:
            try:
                _sig = signature_from_reads(args.input)
                _r0 = LCG.route_from_signature(**_sig)
                print(f"[run] input signature: {_sig['n']} MUC1 reads sampled, "
                      f"AL {_sig['al_frac']:.0%} / AH {_sig['ah_frac']:.0%}, "
                      f"flank tightness {_sig['tight5']:.0%}/{_sig['tight3']:.0%} → route {_r0}",
                      file=sys.stderr)
            except Exception as e:
                print(f"[run] could not characterise the reads ({e})", file=sys.stderr)
                _r0 = "unknown"
        if _r0 == "unknown":
            # Read statistics cannot settle this: flank tightness and median read length were both tried
            # on real data and both misclassify (see length_congruence). Guessing here would silently run
            # the wrong pipeline on a clinical sample, so ask instead of defaulting.
            print("[run] ERROR: cannot tell what this read set is — no MUC1 read was identifiable in the "
                  "sampled reads.\n"
                  "      Declare it with --input-type: pcr-ont | pcr-pacbio | as-ont | "
                  "targeted-pacbio | wgs-ont | wgs-pacbio\n"
                  "      (pcr-* stays on this pipeline; the non-amplicon types are delegated to the streaming "
                  "FASTQ pipeline.)", file=sys.stderr)
            return 2
        if _r0 == "as":
            rc = _delegate_fromfastq(args, _r0)
            if rc is not None:
                return rc
            print("[run] delegation unavailable — continuing on our pipeline", file=sys.stderr)
        else:
            print("[run] amplicon input → our PCR path (frozen chaining flags + cassette fallback, "
                  "which a ruler built for AS does not carry)", file=sys.stderr)

    # ── 1. prepare ───────────────────────────────────────────────────────────
    from .prepare import prepare
    print("[run] step 1/3 — prepare", file=sys.stderr)
    bam = prepare(args.input, args.vntr_ref, base + ".vntr.bam", ref38=args.ref38,
                  threads=args.threads, pacbio=args.pacbio, pcr=args.pcr,
                  no_pcr_autodetect=args.no_pcr_autodetect,
                  recover_unmapped=not args.no_recover_unmapped, sample=args.sample, workdir=args.outdir,
                  pcr_cap=args.pcr_cap)

    # ── 1a. ROUTE + cross-check the length with a SECOND, independent method ────────────────────────
    # The fork is `is_amplicon` x `our flanks present`, not simply AS vs PCR:
    #   · amplicon           -> our PCR path owns it (frozen chaining flags -z/-r, without which a read
    #                           spanning a long tandem does not chain, and the cassette fallback);
    #   · AS/WGS with flanks -> the terrain the burden sweep was written for -> cross-check the length;
    #   · no flanks at all   -> only our cassette can measure it; the sweep would return nothing.
    # The arbiter is authoritative in every case. The sweep is a control, never a vote.
    arbiter = length_cong = None
    try:
        import vntr_raw_length as _V
        arbiter = _V.auto_length(bam, chrom=None, offset=FLANK_TO_PHYSICAL_OFFSET)
    except Exception as e:
        print(f"[run] arbiter length unavailable ({e})", file=sys.stderr)
    from . import length_congruence as LCG
    _route = LCG.route_for(arbiter)
    print(f"[run] route: {_route} (amplicon={(arbiter or {}).get('is_amplicon')}, "
          f"method={(arbiter or {}).get('method')})", file=sys.stderr)
    if LCG.cross_check_applies(arbiter):
        try:
            from .allele_scaffold import sweep_lengths
            _alleles = arbiter.get("alleles") or []
            length_cong = LCG.compare(_alleles, sweep_lengths(bam, _alleles))
            _msg = LCG.render(length_cong)
            if _msg:
                print("[run] " + _msg.rstrip(), file=sys.stderr)
        except Exception as e:
            print(f"[run] length cross-check skipped ({e})", file=sys.stderr)

    # ── 1b. reliable dupC (depth-routed), auto-chained so the user never runs dupc_dispatch by hand ──
    # Amplicon → per-allele pcr_dupc ; AS/WGS → statistical/positional. Reads the ORIGINAL genomic reads at
    # the hg38 dupC locus (the royal-road input = an hg38 haplotagged CRAM). A fastq/uBAM has no genomic
    # locus → skipped (the consensus X-59dupC token still appears, flagged depth-fragile).
    dupc_json = None
    # ROUTING BY INPUT: an aligned BAM/CRAM has a genomic dupC locus to read; a fastq/uBAM does NOT, and
    # used to get no reliable verdict at all — only the consensus token, which is documented depth-fragile.
    # Those inputs now go through VNTR space on the bam `prepare` just produced, per allele. Same question,
    # different anchor.
    _is_aligned = str(args.input).lower().endswith((".bam", ".cram")) and not _looks_unaligned(args.input)
    if not _is_aligned:
        try:
            from .dupc_dispatch import dispatch_dupc_vntr
            # reuse the arbiter computed at step 1a — measuring the same lengths twice would be a second
            # chance to disagree with ourselves
            _lengths = (arbiter or {}).get("alleles") if (arbiter or {}).get("available") else None
            dv = dispatch_dupc_vntr(bam, args.vntr_ref, lengths=_lengths)
            dupc_json = base + ".dupc.json"
            with open(dupc_json, "w") as f:
                json.dump(dv, f, indent=2, default=str)
            _r = dv.get("result", {}) or {}
            print(f"[run] dupC (vntr/runlen_shift, scaffolds {dv.get('scaffolds')}): "
                  f"{_r.get('interpretation')}", file=sys.stderr)
        except Exception as e:
            print(f"[run] VNTR-space dupC call skipped ({e})", file=sys.stderr)
    if _is_aligned:
        try:
            from .dupc_dispatch import dispatch_dupc
            dv = dispatch_dupc(args.input, chrom=args.snp_chrom, start=155188000, end=155192000,
                               ref=args.ref38, fast=True, pon=args.pon)
            dupc_json = base + ".dupc.json"
            with open(dupc_json, "w") as f:
                json.dump(dv, f, indent=2, default=str)
            _r = dv.get("result", {}) or {}
            print(f"[run] dupC ({dv.get('regime')}/{dv.get('caller')}): "
                  f"{_r.get('status') or ('called' if _r.get('called') else 'not called')}", file=sys.stderr)
        except Exception as e:
            print(f"[run] reliable dupC auto-call skipped ({e})", file=sys.stderr)

    # ── 1c. clinical coordinates on a CONFIRMED positive (DETECT with dupc → REPORT here) ──────────────
    # pcr_report gives the Kirby REPEAT index + the carrier-allele decomposition, which the fast detector
    # cannot (its index numbers C-tracts, not array units). Only run on a positive: it is a reporting layer,
    # documented to over-call if used as a detector.
    report_json = None
    # pcr_report reads the ORIGINAL genomic reads at hg38 coordinates to give the Kirby repeat index. A
    # fastq/uBAM has no genomic locus, so it cannot run there at all — attempting it produced a pysam
    # "file has no sequences defined" that reads like a crash rather than a limitation. The position is
    # not lost: the VNTR-space caller reports its own index, and the consensus its motif index.
    if dupc_json and not _is_aligned:
        print("[run] clinical hg38 coordinates (pcr_report) need an aligned BAM/CRAM — not available from "
              "a read file; the variant position comes from the VNTR-space caller instead", file=sys.stderr)
    elif dupc_json:
        try:
            _r = json.load(open(dupc_json)).get("result", {}) or {}
            if _r.get("status") == "positive" or _r.get("called"):
                from .pcr_report import pcr_report
                rep = pcr_report(args.input, chrom=args.snp_chrom, start=155188000, end=155192000,
                                 ref=args.ref38, jobs=args.threads)
                report_json = base + ".pcr_report.json"
                with open(report_json, "w") as f:
                    json.dump(rep, f, indent=2, default=str)
                print(f"[run] clinical coordinates: {rep.get('report', '(none)')}", file=sys.stderr)
        except Exception as e:
            print(f"[run] pcr_report skipped ({e})", file=sys.stderr)

    # ── 2. call ──────────────────────────────────────────────────────────────
    print("[run] step 2/3 — call", file=sys.stderr)
    from . import caller
    call_json = base + ".analyzer.json"
    call_argv = ["-b", bam, "-r", args.vntr_ref, "-s", args.sample,
                 "--min-mq", "0", "--threads", str(args.threads), "--json", call_json,
                 "--pdf", base + ".pdf"]
    if dupc_json:
        call_argv += ["--dupc-json", dupc_json]
    if report_json:
        call_argv += ["--report-json", report_json]
    if args.pcr:
        call_argv.append("--pcr")
    # Two-axis add-on: feed the ORIGINAL genomic reads (a BAM/CRAM) so the report gets the arbiter
    # length + rs4072037/hap. The arbiter is alignment-free (works on a uBAM too); rs4072037 needs the
    # aligned SNP contig and degrades gracefully if absent. A fastq input has no genomic BAM → skipped.
    if str(args.input).lower().endswith((".bam", ".cram")):
        offset = args.copy_offset if args.copy_offset is not None else FLANK_TO_PHYSICAL_OFFSET
        call_argv += ["--genomic-bam", args.input, "--copy-offset", str(offset),
                      "--snp-chrom", args.snp_chrom]
        if args.ref38:
            call_argv += ["--genome-ref", args.ref38]
        if args.snp_pos is not None:
            call_argv += ["--snp-pos", str(args.snp_pos)]
        if args.mut_variant:
            call_argv += ["--mut-variant", args.mut_variant]
            if args.mut_repeat is not None:
                call_argv += ["--mut-repeat", str(args.mut_repeat)]
            if args.mut_allele_len is not None:
                call_argv += ["--mut-allele-len", str(args.mut_allele_len)]
    rc = caller.main(call_argv)
    if rc:
        print("[run] the caller failed — stopping before score", file=sys.stderr)
        return rc or 1

    # ── 3. score ─────────────────────────────────────────────────────────────
    print("[run] step 3/3 — score", file=sys.stderr)
    score_argv = ["--analyzer-json", call_json, "-s", args.sample, "-o", base + ".score.json"]
    # THE ARBITER IS AUTHORITATIVE on length. The caller's haplotype lengths feed mut_len/healthy_len/
    # onset_index today, and it can lose a PCR-depleted allele outright (measured: arbiter 45/77 from two
    # entry points, caller 43/44). Which allele is the MUTANT comes from the dupC caller's carrier contig
    # — a measurement — rather than from the caller's long/short ordering, which is a deduction.
    if arbiter and arbiter.get("available") and arbiter.get("alleles"):
        score_argv += ["--arbiter-alleles", ",".join(str(a) for a in arbiter["alleles"])]
        _carrier_len = None
        try:
            _cc = ((json.load(open(dupc_json)) if dupc_json else {}) or {}).get("result", {}) or {}
            _m = re.search(r"MUC1_VNTR_(\d+)repeats", str(_cc.get("carrier_contig") or ""))
            if _m:
                _carrier_len = int(_m.group(1))
        except Exception:
            _carrier_len = None
        if _carrier_len is not None:
            score_argv += ["--carrier-len", str(_carrier_len)]
    if args.rs4072037_mut:
        score_argv += ["--rs4072037-mut", args.rs4072037_mut]
    else:
        # T2T-native rs4072037: genotypeable on the VNTR-ref BAM itself (fixed 5' flank offset),
        # so `run` needs no GRCh38. A no-call simply leaves the severity axis unset.
        from .detectors.splice_snp import genotype_rs4072037_vntr_ref, genotype_snp
        gt = genotype_rs4072037_vntr_ref(bam, args.vntr_ref)
        dosage = gt.get("dosage") if gt else None
        # Fallback (LM 2026-07-25): VNTR-native starved by a multi-contig smear → read rs4072037 off the
        # original genomic hg38 BAM/CRAM at chr1:155,192,276 (only when the input is aligned).
        if dosage is None and str(args.input).lower().endswith((".bam", ".cram")):
            try:
                g = genotype_snp(args.input, args.ref38, chrom=args.snp_chrom, pos1=args.snp_pos)
                if g.get("snp_genotype") is not None:
                    dosage, gt = g["snp_genotype"], g
                    print(f"[run] rs4072037 VNTR-native undetermined → genomic fallback: "
                          f"{g.get('gt_label')} depth={g.get('depth')}", file=sys.stderr)
            except Exception as e:
                print(f"[run] rs4072037 genomic fallback failed: {e}", file=sys.stderr)
        if dosage is not None:
            score_argv += ["--snp-genotype", str(dosage)]
            print(f"[run] rs4072037: {gt.get('genotype') or gt.get('gt_label')} "
                  f"depth={gt.get('depth')}", file=sys.stderr)
        else:
            print("[run] rs4072037 not callable here → severity axis unset "
                  "(pass --rs4072037-mut to set it)", file=sys.stderr)
    from .cli import main as score_main
    rc = score_main(score_argv)

    # ── the check that matters clinically: arbiter vs the lengths the SCORE rode on ───────────────
    # `clinical_call` turns the caller's haplotype lengths into mut_len / healthy_len / onset_index. When
    # the caller misses an allele the arbiter measured — observed on a real carrier: arbiter 45/77, caller
    # 43/44 — those clinical numbers are computed on alleles that do not exist. Say so, loudly.
    caller_cong = None
    if arbiter and arbiter.get("available"):
        try:
            _haps = (json.load(open(call_json)) or {}).get("haplotypes") or []
            caller_cong = LCG.compare_caller(arbiter.get("alleles"), _haps)
            if caller_cong.get("message"):
                print("[run] " + caller_cong["message"], file=sys.stderr)
        except Exception as e:
            print(f"[run] arbiter-vs-caller length check skipped ({e})", file=sys.stderr)

    # The cross-check belongs in the machine-readable output too, not only on stderr: a length the two
    # methods disagreed about must stay flagged wherever the number is consumed.
    if length_cong is not None or caller_cong is not None:
        try:
            _sj = base + ".score.json"
            _d = json.load(open(_sj))
            if length_cong is not None:
                _d["length_congruence"] = length_cong
            if caller_cong is not None:
                _d["arbiter_vs_caller_length"] = caller_cong
            with open(_sj, "w") as f:
                json.dump(_d, f, indent=2, default=str)
        except Exception as e:
            print(f"[run] could not attach the length cross-check ({e})", file=sys.stderr)

    if rc == 0:
        print(f"\n[run] done → {base}.score.json  (calls: {call_json}, reads: {bam})", file=sys.stderr)
        try:
            res = json.load(open(base + ".score.json"))
            print(f"[run] {args.sample}: carrier={res.get('carrier')} "
                  f"onset_index={res.get('onset_index')} splice={res.get('splice_base')}", file=sys.stderr)
        except Exception:
            pass
    return rc


def _cmd_prepare(argv) -> int:
    from .prepare import main
    return main(argv)


def _cmd_call(argv) -> int:
    from . import caller
    return caller.main(argv)


def _cmd_score(argv) -> int:
    from .cli import main
    return main(argv)


def _cmd_bundle(argv) -> int:
    from .review_bundle import main
    return main(argv)


# Subcommand -> handler NAME, resolved at call time (not bound at import): the handler is then
# a single overridable seam, which keeps routing testable without touching the real pipelines.
_DISPATCH = {"run": "_cmd_run", "prepare": "_cmd_prepare", "call": "_cmd_call",
             "score": "_cmd_score", "bundle": "_cmd_bundle"}


def _handler(sub):
    return globals()[_DISPATCH[sub]]


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    if not argv or argv[0] in ("-h", "--help", "help"):
        print(USAGE)
        return 0
    if argv[0] in ("-V", "--version"):
        from . import __version__
        print(f"MUC1_Analyzer {__version__}")
        return 0

    sub = argv[0]
    if sub in _DISPATCH:
        return _handler(sub)(argv[1:])

    if sub.startswith("-"):                                 # legacy flat call (existing launchers)
        print("[muc1_analyzer] NOTE: flat options without a subcommand are the legacy form; "
              "routing to `score`. Prefer: python -m muc1_analyzer score ...", file=sys.stderr)
        return _cmd_score(argv)

    print(f"[muc1_analyzer] unknown subcommand: {sub}\n", file=sys.stderr)
    print(USAGE, file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
