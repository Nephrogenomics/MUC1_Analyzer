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
    # `.vntr.bam` before `.bam`: the VNTR route's input is `<sample>.vntr.bam`, and leaving either suffix
    # on renames the sample — which is not cosmetic. `arbiter_lengths_of` then looks for
    # `<sample>.bam.score.json`, misses it, and the scaffolds silently fall back to the read modes, i.e.
    # the exact failure this route was built to avoid (measured on one carrier, 2026-08-07).
    for suf in (".haplotagged.cram", ".cram", ".vntr.bam", ".bam"):
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
        # WITHOUT a locus IGV opens the whole contig — ~11.9 kb for a 73-repeat array, where the thing
        # to adjudicate is a 360 bp block. A bundle that produces the right files and does not say where
        # to look has not saved the reviewer anything (measured on one carrier, 2026-08-07).
        if t.get("locus"):
            lines.append(f"goto {t['locus']}")
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


def missing_tools() -> list:
    """External tools the bundle cannot run without. Pure-ish (queries PATH).

    ⚠ On a typical HPC conda environment neither is on the base PATH — `module load samtools minimap2`
    (or the site equivalent) is required first. Without this check the failure is silent and looks like a
    scientific result: `build_allele_view` raises, every allele is caught by the per-allele guard, and the
    run ends with no tracks, no IGV session and a bundle folder that merely looks empty (measured on a
    cluster login node, 2026-08-07)."""
    return [t for t in ("samtools", "minimap2") if shutil.which(t) is None]


def build_allele_view(label: str, contig: str, reads_bam: str, vntr_ref: str, outdir: str, *,
                      sample: str, threads: int = 8, analyzer: str = None, gff3: str = None,
                      region: str = None) -> "dict | None":
    """ONE allele's reviewable view: ad hoc ref + single-contig BAM + consensus + PDF + units + gff3.

    This is the bundle. Everything above it only decides WHERE the per-allele read sets come from —
    HP tags on a haplotagged CRAM, or the allele scaffolds of a `.vntr.bam`. Factored out on 2026-08-07
    rather than copied for the PCR route: the IGV view, the consensus and the feature vector are the same
    artefacts whatever separated the alleles, and two copies would drift.

    `region` subsets `reads_bam` to one contig before extraction — the PCR route's `.vntr.bam` holds all
    150 contigs, exactly like the HP route's multi-contig source.

    ⚠ The reads are RE-ALIGNED onto the single contig, never subsetted in place. A subset keeps only
    flank fragments at MAPQ 0, which IGV hides → an empty view (MUC1_log 2026-07-09). That trap is
    identical on both routes, which is the clearest sign they are one procedure."""
    analyzer = analyzer or os.path.join(os.path.dirname(os.path.abspath(__file__)), "caller.py")
    work = os.path.join(outdir, "work")
    os.makedirs(work, exist_ok=True)
    out_bam = os.path.join(outdir, f"{label}.bam")
    out_ref = extract_contig_ref(vntr_ref, contig, os.path.join(outdir, f"{label}_ref.fa"))
    reads_fq = os.path.join(work, f"{sample}_{label}.reads.fq")
    view = ["samtools", "view", "-b", reads_bam] + ([region] if region else [])
    with open(reads_fq, "w") as fq:
        if region:
            sub = subprocess.run(view, check=True, capture_output=True)
            subprocess.run(["samtools", "fastq", "-F", "0x900", "-"], input=sub.stdout,
                           check=True, stdout=fq, stderr=subprocess.DEVNULL)
        else:
            subprocess.run(["samtools", "fastq", "-F", "0x900", reads_bam],
                           check=True, stdout=fq, stderr=subprocess.DEVNULL)
    mm = subprocess.run(["minimap2", "-ax", "map-ont", "--secondary=no", "-t", str(threads),
                         out_ref, reads_fq], check=True, capture_output=True)
    subprocess.run(["samtools", "sort", "-o", out_bam, "-"], input=mm.stdout, check=True)
    subprocess.run(["samtools", "index", out_bam], check=True)
    subprocess.run([sys.executable, analyzer, "-b", out_bam, "-r", out_ref,
                    "-s", f"{sample}_{label}", "--min-mq", "0", "--top-n", "1",
                    "--pdf", os.path.join(outdir, f"{label}.report.pdf"),
                    "--fasta-consensus", os.path.join(outdir, f"{label}.consensus.fa"),
                    "--json", os.path.join(outdir, f"{label}.units.json")], check=True)
    track = {"hp": label, "ref": out_ref, "bam": out_bam}
    g = extract_contig_gff3(gff3, contig, os.path.join(outdir, f"{label}.gff3")) if gff3 else None
    if g:
        track["gff3"] = g
    return track


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
                tracks.append(build_allele_view(f"hp{hp}", contig, realigned, vntr_ref, outdir,
                                                sample=sample, threads=threads, analyzer=analyzer,
                                                gff3=gff3))
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


#: C-tract target -> (motif label, base-pair change). The run-length caller is parameterised by the
#: mutant tract length, so the family it means is a lookup, not a second detector.
_TARGET_MOTIF = {8: ("59dupC", 1), 5: ("58_59delCC", -2), 6: ("58_59delC", -1), 11: ("dupCCCC", 4)}

#: C-tract target -> the name `pcr_variant`/`variant_profile` know it by (the segment-typing vocabulary).
_TARGET_PROBE = {8: "dupC", 5: "delCC", 6: "delC", 11: "dupCCCC"}


def array_origin(units_json: str) -> "int | None":
    """Scaffold coordinate of array unit 0, from the caller's OWN units export. Impure.

    ⚠ NOT `flank5_len`. That finds the first EXACT `MOTIF_A`, and the caller's decomposition starts
    earlier, on units (`1-2-3-4-5`) that are not that motif — measured on a carrier's 73-repeat contig,
    `flank5_len` = 4578 against a unit-0 at 4304. The gap is 274 bp = 4.57 units, so the two frames are
    out of PHASE, not merely offset by whole repeats, and a locus built on the wrong one lands ~5 units
    away. Taking the origin from the same artefact the reviewer reads makes the frames agree by
    construction."""
    import json as _json
    try:
        with open(units_json) as fh:
            haps = (_json.load(fh) or {}).get("haplotypes") or []
        motifs = (haps[0] or {}).get("motifs") if haps else None
        return motifs[0].get("pos") if motifs else None
    except Exception:
        return None


def block_locus(units_json: str, contig: str, block: dict, *, unit_bp: int = 60, pad: int = 120):
    """`contig:start-end` for a run of array units, or None. Impure (reads the caller's units export).

    Returns None on any failure: sending a reviewer to a computed-wrong locus is worse than sending them
    to the whole contig, because they would believe the frame."""
    if not block or not block.get("len"):
        return None
    origin = array_origin(units_json)
    if origin is None:
        return None
    start = int(origin) + int(block["start"]) * unit_bp - pad
    end = int(origin) + (int(block["end"]) + 1) * unit_bp + pad
    return f"{contig}:{max(1, start)}-{end}"


def pcr_allele_data(result: dict, *, length=None, repeat_index=None, mut_len: int = 8) -> dict:
    """A PCR scaffold's runlen_shift result -> the per-allele shape `feature_row`/`alarm_line` expect. Pure.

    The two routes separate alleles differently (HP tags vs scaffolds) and agree on everything after, so
    the adapter lives here and the reporting code stays single-source. Fields the PCR route genuinely
    cannot supply stay None — an absent value is honest, a fabricated one is not."""
    cand = (result or {}).get("candidate") or {}
    motif, delta = _TARGET_MOTIF.get(mut_len, (f"mut_len{mut_len}", None))
    return {"length": length, "n_reads": cand.get("n"),
            "has_mut": bool((result or {}).get("called")),
            "array_wide": bool((result or {}).get("array_wide")),
            "mutation_variant": {"motif": motif if (result or {}).get("called") else None,
                                 "delta_bp": delta if (result or {}).get("called") else None,
                                 "repeat_index": repeat_index}}


def _write_profile(vntr_bam: str, contig: str, outdir: str, label: str, mut_len: int) -> "str | None":
    """Write `<label>.profile.txt` (variant fraction repeat by repeat) and return the IGV locus. Impure.

    Runs AFTER `build_allele_view`, because the locus frame comes from that step's `<label>.units.json`
    — the grid the reviewer actually sees — and not from an independently recomputed flank."""
    import pysam
    from .pcr_variant import _scan_segments
    from .variant_profile import focality, index_profile, longest_block, render
    probe = _TARGET_PROBE.get(mut_len)
    if not probe:
        return None
    with pysam.AlignmentFile(vntr_bam, "rb") as fh:
        length = fh.get_reference_length(contig)
    prof = index_profile(_scan_segments(vntr_bam, contig, 1, int(length), None), probe)
    blk = longest_block(prof)
    with open(os.path.join(outdir, f"{label}.profile.txt"), "w") as fh:
        fh.write(f"{contig}  probe={probe}\n")
        fh.write(render(prof, focality(prof), target=probe) + "\n")
    return block_locus(os.path.join(outdir, f"{label}.units.json"), contig, blk)


def build_bundle_pcr(vntr_bam: str, vntr_ref: str, outdir: str, *, sample: str = None,
                     lengths=None, dupc_json: str = None, features_tsv: str = None,
                     threads: int = 8, analyzer: str = None, gff3: str = None,
                     mut_len: int = 8) -> dict:
    """The SAME bundle, sourced from an already-scaffolded `.vntr.bam` instead of a haplotagged CRAM.

    There is no second bundle. `run` on a FASTQ has already separated the alleles onto per-length
    contigs, so the step the HP route spends `detect_analyzer_phased` on is simply already done — this
    route is shorter, not different. Everything downstream is `build_allele_view`, shared verbatim.

    ⚠ Scaffolds come from `allele_contigs(lengths=arbiter)`, never from the workbook: a run scaffolded on
    a false-SHORT allele pair looks at the wrong contigs and never examines the mutant allele
    (MUC1_log 2026-08-07, VAL/SIM)."""
    import json as _json
    from .allele_scaffold import CONTIG_LEN, allele_contigs
    from .dupc_dispatch import _unit_index

    sample = sample or _sid(vntr_bam)
    if gff3 is None:
        cand = os.path.splitext(vntr_ref)[0] + ".gff3"
        gff3 = cand if os.path.exists(cand) else None
    os.makedirs(outdir, exist_ok=True)

    doc = {}
    if dupc_json and os.path.exists(dupc_json):
        with open(dupc_json) as fh:
            doc = _json.load(fh)
    per_allele = ((doc.get("result") or doc).get("per_allele")) or {}

    contigs = allele_contigs(vntr_bam, lengths=lengths)
    tracks, rows, per_label = [], [], {}
    lens = {}
    for c in contigs:
        m = CONTIG_LEN.search(c)
        lens[c] = int(m.group(1)) if m else None
    for c in contigs:
        label = f"allele{lens[c]}" if lens[c] else c
        r = per_allele.get(c) or {}
        idx = _unit_index(vntr_ref, c, ((r.get("candidate") or {}).get("ref_pos")))
        per_label[label] = pcr_allele_data(r, length=lens[c], repeat_index=idx, mut_len=mut_len)
        # The per-repeat profile, written next to the view. A CONSENSUS cannot show a minority variant —
        # That carrier's delCC is ~9 % of reads and its 73 units all decode exact — so without this the bundle
        # produces the right files and still does not point at anything (2026-08-07).
        try:
            t = build_allele_view(label, c, vntr_bam, vntr_ref, outdir, sample=sample,
                                  threads=threads, analyzer=analyzer, gff3=gff3, region=c)
            try:                                   # AFTER the view: the locus frame comes from its JSON
                locus = _write_profile(vntr_bam, c, outdir, label, mut_len)
                if locus:
                    t["locus"] = locus
            except Exception as e:
                print(f"[warn] profile {sample} {label}: {e}", file=sys.stderr)
            tracks.append(t)
        except (subprocess.CalledProcessError, OSError) as e:
            print(f"[warn] bundle {sample} {label}: {e}", file=sys.stderr)
    for label, d in per_label.items():
        other = [v["length"] for k, v in per_label.items() if k != label]
        rows.append(feature_row(sample, label, d, other[0] if other else None,
                                is_mut_hp=bool(d.get("has_mut"))))

    if tracks:
        igv_batch(sample, tracks, os.path.join(outdir, "igv_session.batch"))
    # `alarm_line` reads {"notes": {"per_hp": …}} — handing it the same shape keeps ONE alarm format
    # across both routes rather than a second one a reader would have to learn.
    alarm = alarm_line(sample, {"notes": {"per_hp": per_label}})
    # The reason this bundle was opened, when there is one. A call that survived no change of scaffold
    # length is the case the reviewer is here to settle, and the ALARM file is the first thing read.
    _rep = ((doc.get("result") or {}) if doc else {}).get("scaffold_replication") or {}
    if _rep.get("neighbours_tested") and not _rep.get("replicated"):
        alarm += (f"\t⚠ did not replicate on any of {_rep['neighbours_tested']} neighbouring "
                  f"scaffold(s) — signature of a PLACEMENT artifact, not a variant")
    with open(os.path.join(outdir, "ALARM.txt"), "w") as fh:
        fh.write(alarm + "\n")
    if features_tsv:
        append_features(features_tsv, rows)
    return {"sample": sample, "outdir": outdir, "tracks": len(tracks), "alarm": alarm}


def main(argv=None):
    import argparse
    import glob as globmod
    import json
    ap = argparse.ArgumentParser(prog="muc1_analyzer.review_bundle",
                                 description="Per-haplotype IGV review bundle + feature capture")
    ap.add_argument("--glob", help="pattern of patient CRAM/BAM files (or --list)")
    ap.add_argument("--list", dest="list_file", help="file of CRAM paths (takes precedence over --glob)")
    # PCR ROUTE: the alleles are already on their own contigs, so no CRAM, no --ref, no phasing step.
    ap.add_argument("--vntr-bam", nargs="+", default=None,
                    help="PCR route: `<sample>.vntr.bam` file(s) `run` already wrote (no --ref needed)")
    ap.add_argument("--lengths", default=None, metavar="43,96",
                    help="PCR route: the arbiter's allele lengths (single --vntr-bam only)")
    ap.add_argument("--mut-len", type=int, default=8, help="PCR route: 8 = dupC, 5 = delCC")
    ap.add_argument("--ref", default=None, help="GRCh38 genome (CRAM route: required)")
    ap.add_argument("--vntr-ref", required=True, help="multi-contig VNTR reference")
    ap.add_argument("--gff3", default=None,
                    help="multi-contig GFF3 (repeats/exons); default = <vntr-ref without ext>.gff3 if it exists")
    ap.add_argument("--outdir", required=True, help="root folder of the bundles (one subfolder per patient)")
    ap.add_argument("--pon", default=None, help="dupC PoN table (for the alarm/features)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--skip-done", action="store_true",
                    help="skip a patient if its <sid>/ALARM.txt already exists (resume after timeout)")
    args = ap.parse_args(argv)
    if args.vntr_bam:
        miss = missing_tools()
        if miss:
            print(f"[ERROR] not on PATH: {', '.join(miss)} — run `module load samtools minimap2` first. "
                  "Refusing to start: without them every allele fails silently and the bundle folder "
                  "just looks empty.", file=sys.stderr)
            return 1
        lens = [int(x) for x in args.lengths.split(",")] if args.lengths else None
        if lens and len(args.vntr_bam) > 1:
            print("[ERROR] --lengths applies to ONE --vntr-bam; with several, let the arbiter of each "
                  "sample speak via its own score.json", file=sys.stderr)
            return 1
        from .miss_report import arbiter_lengths_of
        rc = 0
        for bam in args.vntr_bam:
            sid = _sid(bam)
            out = os.path.join(args.outdir, sid)
            if args.skip_done and os.path.exists(os.path.join(out, "ALARM.txt")):
                print(f"[bundle] {sid}: done, kept", file=sys.stderr)
                continue
            l = lens or arbiter_lengths_of(os.path.dirname(bam), sid)
            try:
                r = build_bundle_pcr(bam, args.vntr_ref, out, sample=sid, lengths=l,
                                     dupc_json=os.path.join(os.path.dirname(bam), f"{sid}.dupc.json"),
                                     features_tsv=os.path.join(args.outdir, "review_features.tsv"),
                                     threads=args.threads, gff3=args.gff3, mut_len=args.mut_len)
                print(f"[bundle] {sid}: {r['tracks']} allele(s) · {r['alarm']}", file=sys.stderr)
            except Exception as e:                     # one sample must not lose the batch
                print(f"[warn] bundle {sid}: {e}", file=sys.stderr)
                rc = 1
        return rc
    if not args.ref:
        ap.error("--ref is required on the CRAM route (or use --vntr-bam)")
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
