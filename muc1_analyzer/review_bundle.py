"""Per-haplotype VISUAL REVIEW bundle + feature capture for a future classifier.

Goal (Ilias/Xavier spec): for each patient, produce what is needed to ADJUDICATE BY EYE a
suspected frameshift, and — at the same time — accumulate a labelled training set.

Per patient, for EACH haplotype (HP1/HP2):
  - `hpN.bam` (+ .bai)  : the HP reads RE-ALIGNED onto the SINGLE best contig (the winning VNTR
                          length) → MAPQ 60, single-contig header, reads spanning the array → clean
                          IGV view. (⚠ do NOT subset the multi-contig: flank fragments at MAPQ 0,
                          which IGV hides → empty view.)
  - `hpN_ref.fa` (+.fai): that single contig, extracted from the multi-contig ref → ad hoc IGV reference.
  - `hpN.consensus.fa`  : the haplotype consensus (MUC1_Analyzer `--fasta-consensus`).
  - `hpN.report.pdf`    : the native visual report (MUC1_Analyzer `--pdf`).
  - `hpN.units.tsv`     : motif-by-motif decomposition (from the analyzer JSON).
Plus, at the patient level:
  - `igv_session.batch` : IGV batch script that loads ref + BAM of each HP (one click → view).
  - `ALARM.txt`         : 1 line — suspected frameshift? HP? variant? position? in-frame? tier.
  - append to `review_features.tsv` (shared path): the FEATURE VECTOR + an EMPTY `verdict`
                          column, to be filled in review (real/artifact/unsure) = the LABEL.

⚠ The `build_bundle` orchestrator needs samtools/minimap2 + data (on an HPC cluster); the pure bricks
(contig extraction, IGV script, feature line) are tested locally. Does NOT import pysam at
module level — stays usable without data.
"""
from __future__ import annotations
import csv
import os
import shutil
import subprocess
import sys

def _sid(path: str) -> str:
    """Patient identifier from the CRAM path (same convention as batch_score_phased)."""
    b = os.path.basename(path)
    for suf in (".haplotagged.cram", ".cram"):
        if b.endswith(suf):
            b = b[: -len(suf)]
    return b.replace(".5mc.sup.unaligned", "").replace(".5mC.sup.unaligned", "")


_FEATURE_COLS = [
    "sample", "hp", "is_mut_hp", "n_reads", "vntr_len", "vntr_len_other",
    "has_mut", "variant", "delta_bp", "repeat_index", "in_frame", "array_wide_artifact",
    "dupc_called", "dupc_tier", "dupc_posterior", "n8C", "tot", "frac_8C", "f_null",
    "p_bonf", "p_paired", "hp_asymmetry",
    "verdict", "verdict_by", "verdict_date",   # <- to fill in review: the classifier LABEL
]


def best_contig_name(length) -> str | None:
    """Ad hoc contig name for a VNTR length (multi-contig ref convention)."""
    if length in (None, "", 0):
        return None
    return f"MUC1_VNTR_{int(length)}repeats"


def extract_contig_ref(vntr_ref: str, contig: str, out_fa: str) -> str:
    """Write `contig` alone from the multi-contig ref into `out_fa` (+ .fai). Ad hoc IGV ref."""
    import pysam
    fa = pysam.FastaFile(vntr_ref)
    try:
        seq = fa.fetch(contig)
    finally:
        fa.close()
    with open(out_fa, "w") as fh:
        fh.write(f">{contig}\n")
        for i in range(0, len(seq), 70):
            fh.write(seq[i:i + 70] + "\n")
    pysam.faidx(out_fa)
    return out_fa


def extract_contig_gff3(gff3: str, contig: str, out_gff3: str):
    """Write the GFF3 features of `contig` only (VNTR_repeat_N + exons) into `out_gff3` → IGV track
    annotating the repeats. Pure (filters by col1). Returns the path, or None if absent/empty."""
    if not gff3 or not os.path.exists(gff3):
        return None
    kept = [ln for ln in open(gff3)
            if not ln.startswith("#") and ln.split("\t", 1)[0] == contig]
    if not kept:
        return None
    with open(out_gff3, "w") as fh:
        fh.write("##gff-version 3\n")
        fh.writelines(kept)
    return out_gff3


def igv_batch(sample: str, tracks: list[dict], out_path: str) -> str:
    """Generate an IGV batch script (one per HP: ad hoc genome + BAM). `tracks` = [{hp,ref,bam}].

    IGV loads only one genome at a time: we emit one block per HP (genome → load → snapshot),
    runnable via `igv.sh -b igv_session.batch`. Paths RELATIVE to the bundle folder.
    """
    lines = ["new", "maxPanelHeight 2000"]
    for t in tracks:
        base = os.path.dirname(out_path)
        ref = os.path.relpath(t["ref"], base)
        bam = os.path.relpath(t["bam"], base)
        lines += [f"genome {ref}", f"load {bam}"]
        if t.get("gff3"):                                    # repeats/exons annotation track
            lines.append(f"load {os.path.relpath(t['gff3'], base)}")
        lines += ["snapshotDirectory .", f"snapshot {sample}_hp{t['hp']}.png", "new"]
    text = "\n".join(lines) + "\n"
    with open(out_path, "w") as fh:
        fh.write(text)
    return text


def _fs_variant(hp_data: dict) -> dict:
    """Extract the representative frameshift variant of an HP (from the analyzer JSON)."""
    mv = hp_data.get("mutation_variant") or {}
    return {"variant": mv.get("motif"), "delta_bp": mv.get("delta_bp"),
            "repeat_index": mv.get("repeat_index"),
            "in_frame": (mv.get("delta_bp") is not None and mv.get("delta_bp") % 3 == 0)}


def feature_row(sample: str, hp: str, hp_data: dict, other_len, *, is_mut_hp: bool,
                dupc: dict = None) -> dict:
    """Build the feature line (`verdict` column left EMPTY = label to collect)."""
    fs = _fs_variant(hp_data)
    best = (dupc or {}).get("best") or {}
    # inter-HP asymmetry of the 8C (signature of an allele-specific event)
    tot, n8 = best.get("tot"), best.get("n8C")
    toth, n8h = best.get("tot_healthy"), best.get("n8C_healthy")
    asym = None
    if tot and n8 is not None and toth and n8h is not None and toth > 0:
        r_other = n8h / toth
        asym = round((n8 / tot) / r_other, 2) if r_other > 0 else None
    row = {c: "" for c in _FEATURE_COLS}
    row.update({
        "sample": sample, "hp": hp, "is_mut_hp": is_mut_hp,
        "n_reads": hp_data.get("n_reads"), "vntr_len": hp_data.get("length"),
        "vntr_len_other": other_len, "has_mut": hp_data.get("has_mut"),
        "variant": fs["variant"], "delta_bp": fs["delta_bp"], "repeat_index": fs["repeat_index"],
        "in_frame": fs["in_frame"], "array_wide_artifact": bool(hp_data.get("array_wide")),
        "dupc_called": bool(dupc and dupc.get("called")), "dupc_tier": (dupc or {}).get("tier"),
        "dupc_posterior": best.get("posterior"), "n8C": n8, "tot": tot,
        "frac_8C": best.get("frac_8C"), "f_null": best.get("f"),
        "p_bonf": best.get("p_bonf"), "p_paired": best.get("p_paired"), "hp_asymmetry": asym,
    })
    return row


def append_features(path: str, rows: list[dict]) -> None:
    """Append rows to the shared feature TSV (header written if the file is new)."""
    new = not os.path.exists(path) or os.path.getsize(path) == 0
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=_FEATURE_COLS, delimiter="\t", extrasaction="ignore")
        if new:
            w.writeheader()
        for r in rows:
            w.writerow({k: ("" if v is None else v) for k, v in r.items()})


def alarm_line(sample: str, phased: dict, dupc: dict = None) -> str:
    """1 ALARM line: suspected frameshift? HP? variant? position? in-frame? tier."""
    per_hp = (phased.get("notes") or {}).get("per_hp") or {}
    hits = []
    for hp, d in per_hp.items():
        if d.get("has_mut"):
            fs = _fs_variant(d)
            tag = "IN-FRAME(benign)" if fs["in_frame"] else "FRAMESHIFT"
            hits.append(f"HP{hp}:{tag} {fs['variant']} @rep{fs['repeat_index']}")
    tier = (dupc or {}).get("tier") or "none"
    flag = "SUSPECT" if any("FRAMESHIFT" in h for h in hits) else "clear"
    detail = " ; ".join(hits) if hits else "no frameshift"
    return f"{sample}\t{flag}\tdupC_tier={tier}\t{detail}"


def build_bundle(cram: str, vntr_ref: str, ref: str, outdir: str, *, sample: str = None,
                 pon: dict = None, features_tsv: str = None, threads: int = 8,
                 analyzer: str = None, gff3: str = None) -> dict:
    """Orchestrator: phased (keep=True) → per HP {ad hoc BAM, ref, consensus, PDF, units} + IGV +
    alarm + features. Best-effort per HP. Requires samtools/minimap2 + data (on an HPC cluster)."""
    from .detectors import vntr as V
    sample = sample or _sid(cram)
    if gff3 is None:                                         # default: the .gff3 next to the VNTR ref
        cand = os.path.splitext(vntr_ref)[0] + ".gff3"
        gff3 = cand if os.path.exists(cand) else None
    os.makedirs(outdir, exist_ok=True)
    work = os.path.join(outdir, "work")
    analyzer = analyzer or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "caller.py")

    phased = V.detect_analyzer_phased(cram, vntr_ref, sample=sample, genome_ref=ref,
                                      threads=threads, keep=True, workdir=work)
    dupc = None
    try:
        from .dupc_caller import call_dupc
        dupc = call_dupc(cram, genome_ref=ref, pon=pon)
    except Exception as e:
        print(f"[warn] dupC {sample}: {e}", file=sys.stderr)

    per_hp = (phased.get("notes") or {}).get("per_hp") or {}
    len_by_hp = {hp: per_hp.get(hp, {}).get("length") for hp in ("1", "2")}
    tracks, rows = [], []
    for hp in ("1", "2"):
        d = per_hp.get(hp) or {}
        contig = best_contig_name(d.get("length"))
        realigned = os.path.join(work, f"{sample}_hp{hp}.realigned.bam")
        if contig and os.path.exists(realigned):
            try:
                hp_bam = os.path.join(outdir, f"hp{hp}.bam")
                hp_ref = extract_contig_ref(vntr_ref, contig, os.path.join(outdir, f"hp{hp}_ref.fa"))
                # RE-ALIGN all the HP reads onto the SINGLE contig `hp_ref` (instead of subsetting
                # the multi-contig). The subset kept only flank fragments at MAPQ 0 (IGV hides them)
                # → empty view; here single-contig header + MAPQ 60 + reads spanning the array
                # present → clean IGV view (cf. MUC1_log 2026-07-09).
                reads_fq = os.path.join(work, f"{sample}_hp{hp}.reads.fq")
                with open(reads_fq, "w") as fq:
                    subprocess.run(["samtools", "fastq", "-F", "0x900", realigned],
                                   check=True, stdout=fq, stderr=subprocess.DEVNULL)
                mm = subprocess.run(["minimap2", "-ax", "map-ont", "--secondary=no",
                                     "-t", str(threads), hp_ref, reads_fq],
                                    check=True, capture_output=True)
                subprocess.run(["samtools", "sort", "-o", hp_bam, "-"], input=mm.stdout, check=True)
                subprocess.run(["samtools", "index", hp_bam], check=True)
                subprocess.run([sys.executable, analyzer, "-b", hp_bam, "-r", hp_ref,
                                "-s", f"{sample}_hp{hp}", "--min-mq", "0", "--top-n", "1",
                                "--pdf", os.path.join(outdir, f"hp{hp}.report.pdf"),
                                "--fasta-consensus", os.path.join(outdir, f"hp{hp}.consensus.fa"),
                                "--json", os.path.join(outdir, f"hp{hp}.units.json")], check=True)
                hp_gff3 = extract_contig_gff3(gff3, contig, os.path.join(outdir, f"hp{hp}.gff3"))
                track = {"hp": hp, "ref": hp_ref, "bam": hp_bam}
                if hp_gff3:
                    track["gff3"] = hp_gff3
                tracks.append(track)
            except (subprocess.CalledProcessError, OSError) as e:
                print(f"[warn] bundle {sample} hp{hp}: {e}", file=sys.stderr)
        other = "2" if hp == "1" else "1"
        rows.append(feature_row(sample, hp, d, len_by_hp.get(other),
                                is_mut_hp=(phased.get("hp_mut") == hp), dupc=dupc))

    if tracks:
        igv_batch(sample, tracks, os.path.join(outdir, "igv_session.batch"))
    with open(os.path.join(outdir, "ALARM.txt"), "w") as fh:
        fh.write(alarm_line(sample, phased, dupc) + "\n")
    if features_tsv:
        append_features(features_tsv, rows)
    return {"sample": sample, "outdir": outdir, "tracks": len(tracks),
            "alarm": alarm_line(sample, phased, dupc)}


def main(argv=None):
    import argparse
    import glob as globmod
    import json
    ap = argparse.ArgumentParser(prog="muc1_analyzer.review_bundle",
                                 description="Per-haplotype IGV review bundle + feature capture")
    ap.add_argument("--glob", help="pattern of patient CRAM/BAM files (or --list)")
    ap.add_argument("--list", dest="list_file", help="file of CRAM paths (takes precedence over --glob)")
    ap.add_argument("--ref", required=True, help="GRCh38 genome (for CRAM + genotyping)")
    ap.add_argument("--vntr-ref", required=True, help="multi-contig VNTR reference")
    ap.add_argument("--gff3", default=None,
                    help="multi-contig GFF3 (repeats/exons); default = <vntr-ref without ext>.gff3 if it exists")
    ap.add_argument("--outdir", required=True, help="root folder of the bundles (one subfolder per patient)")
    ap.add_argument("--pon", default=None, help="dupC PoN table (for the alarm/features)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--skip-done", action="store_true",
                    help="skip a patient if its <sid>/ALARM.txt already exists (resume after timeout)")
    args = ap.parse_args(argv)
    pon = json.load(open(args.pon)) if args.pon else None
    if args.list_file:
        with open(args.list_file) as fh:
            paths = [ln.strip() for ln in fh if ln.strip()]
    elif args.glob:
        paths = sorted(globmod.glob(args.glob))
    else:
        print("[ERROR] provide --glob OR --list", file=sys.stderr)
        return 1
    if not paths:
        print(f"[ERROR] no CRAM: {args.list_file or args.glob}", file=sys.stderr)
        return 1
    os.makedirs(args.outdir, exist_ok=True)
    features_tsv = os.path.join(args.outdir, "review_features.tsv")
    for i, cram in enumerate(paths, 1):
        sid = _sid(cram)
        if args.skip_done and os.path.exists(os.path.join(args.outdir, sid, "ALARM.txt")):
            print(f"  [{i}/{len(paths)}] {sid}: already done", file=sys.stderr)
            continue
        print(f"  [{i}/{len(paths)}] {sid}…", file=sys.stderr, flush=True)
        try:
            r = build_bundle(cram, args.vntr_ref, args.ref, os.path.join(args.outdir, sid),
                             sample=sid, pon=pon, features_tsv=features_tsv, threads=args.threads,
                             gff3=args.gff3)
            print(f"      {r['tracks']} HP · {r['alarm']}", file=sys.stderr)
        except Exception as e:
            print(f"      ERROR {sid}: {e}", file=sys.stderr)
    print(f"\n[OK] bundles -> {args.outdir}  | features -> {features_tsv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
