"""`muc1_analyzer prepare` — universal input → MUC1 VNTR-ready BAM.

Hides samtools/minimap2 behind one step. Ports the semantics of the validated
in-house alignment recipe into the package so a third-party user never has to know the
recipe. Two things it must not lose:

  * **phase + methylation tags** — `samtools fastq -T MM,ML,HP,PS` + `minimap2 -y` (the reason we
    do not just pipe plain fastq): HP/PS keep a haplotagged input's PHASE alive across the realign
    to the VNTR reference (so per-HP work needs no re-phasing), MM/ML keep 5mC/5hmC for paper 2.
  * **reads that hg38 hides** — an hg38-aligned input is RESTRICTED to the MUC1 window
    before re-extraction, but a uBAM/fastq is taken WHOLE: hg38 cannot represent the long
    VNTR, so its alignment omits reads that the VNTR reference does recover.

Input kind is auto-detected (fastq / uBAM / aligned hg38 / aligned T2T / already on the
VNTR ref) and an input already aligned to the target is BYPASSED, not realigned.

The decision layer (`detect_input_kind`, `classify_reference`, `is_already_target`,
`extraction_region`) is pure and unit-tested; only the plumbing shells out.
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys

from .config import PCR_ALIGN_CAP, PCR_AUTODETECT_MIN_READS, PCR_MINIMAP_EXTRA

# ── reference fingerprints ────────────────────────────────────────────────────
# chr1 lengths that identify the coordinate system of an aligned input.
_HG38_CHR1 = (248_956_422, 249_250_621)      # GRCh38 / GRCh37
_T2T_CHR1 = (248_387_328,)                   # CHM13v2.0
_T2T_NAMES = ("NC_060925.1", "CP068277.2")

# MUC1 windows. hg38 is deliberately WIDER than config.GRCh38["LOCUS"]: prepare wants every
# read touching the locus (incl. soft-clipped VNTR reads), not just the scoring window.
# The chr1 coordinates are fixed; the CONTIG NAME is resolved from the header at run time
# (`chr1` vs `1` vs RefSeq) — a `1`-named GRCh38 BAM (e.g. a noalt build) would otherwise extract
# 0 reads from a hardcoded `chr1:…` region. MUC1_HG38_WINDOW is the chr1-named default/reference.
MUC1_HG38_COORDS = "155180000-155200000"
MUC1_HG38_WINDOW = f"chr1:{MUC1_HG38_COORDS}"

# Tags carried through the fastq round-trip into the realigned BAM (`samtools fastq -T` + `minimap2 -y`).
# HP,PS keep an upstream HAPLOTAGGED input's phase alive across the realignment to the VNTR reference
# (verified: samtools fastq -T HP,PS then minimap2 -y round-trips both) → per-HP rs4072037 / per-allele
# work needs no re-phasing. MM,ML keep 5mC/5hmC methylation (paper-2 use; free to carry). Reads with none
# of these tags simply carry none — safe for fastq/uBAM/unphased inputs.
FASTQ_CARRY_TAGS = "MM,ML,HP,PS"
MUC1_T2T_WINDOW = "muc1win:2327000-2331500"
# Full-genome T2T-CHM13 (chr1 / NC_060925.1): the MUC1 VNTR sits at ~154.32–154.34 Mb
# (config.T2T["VNTR"] = 154,324,904–154,335,104). Like the hg38 window, deliberately WIDER (pad ~5 kb)
# to keep soft-clipped VNTR reads; the CONTIG NAME is resolved from the header at run time. Without this
# a full-genome T2T BAM/CRAM falls through extraction_region → None → every read kept → bloated BAM
# (identical failure to a whole-run fastq: the input is aligned but never restricted to the locus).
MUC1_T2T_COORDS = "154320000-154340000"
KIRBY_CONTIG_PREFIX = "MUC1_VNTR_"

FASTQ_EXTS = (".fastq", ".fastq.gz", ".fq", ".fq.gz")


# ── pure decision layer ───────────────────────────────────────────────────────
def detect_input_kind(path: str, sqs: list | None) -> str:
    """'fastq' | 'ubam' | 'aligned'. `sqs` = [(contig, length)] from the header (None for fastq)."""
    if str(path).endswith(FASTQ_EXTS):
        return "fastq"
    return "aligned" if sqs else "ubam"


def classify_reference(sqs: list | None) -> str:
    """Which coordinate system an aligned input sits in: 'kirby' | 't2t_win' | 'hg38' | 't2t' | 'unknown'."""
    if not sqs:
        return "unknown"
    names = {n for n, _ in sqs}
    if any(n.startswith(KIRBY_CONTIG_PREFIX) for n in names):
        return "kirby"
    if "muc1win" in names:
        return "t2t_win"
    if names & set(_T2T_NAMES):
        return "t2t"
    for n, ln in sqs:
        if ln in _T2T_CHR1:
            return "t2t"
        if ln in _HG38_CHR1:
            return "hg38"
    return "unknown"


def autopcr_decision(*, amplicon_scale: bool, explicit_pcr: bool, disabled: bool) -> bool:
    """Should prepare auto-enable PCR mode? Only when the locus is amplicon-scale AND the user neither
    set --pcr already nor disabled autodetection. Pure → unit-tested."""
    return amplicon_scale and not explicit_pcr and not disabled


def vntr_aligned_count(fq: str, vntr_ref: str, *, threads: int = 8, pacbio: bool = False) -> int:
    """How many reads of `fq` actually MAP to the VNTR reference (default preset, primary only).

    The amplicon-vs-AS/WGS discriminator must count REAL MUC1-VNTR reads, not the raw locus window — the
    20 kb hg38 window also holds flanking-gene / off-target reads, so a busy AS window could otherwise be
    mis-read as a deep amplicon. This is a cheap detection pass on the window (a few seconds) before the
    final alignment. Runs only when auto-detection is actually in play (not --pacbio / --pcr / disabled)."""
    _require_tools("minimap2", "samtools")
    preset = minimap_preset(pacbio)
    p1 = subprocess.Popen(["minimap2", "-ax", preset, "--secondary=no", "-t", str(threads), vntr_ref, fq],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    out = subprocess.run(["samtools", "view", "-c", "-F", "0x904", "-"],
                         stdin=p1.stdout, stdout=subprocess.PIPE, text=True).stdout
    p1.wait()
    try:
        return int(out.strip() or 0)
    except ValueError:
        return 0


def _reads_at_least(fq: str, k: int) -> bool:
    """True if the fastq holds >= k reads. Early-exits at k*4 lines so a huge amplicon (or a big fastq)
    is not fully scanned — we only need the threshold verdict, not the exact count. Handles .gz."""
    import gzip
    opener = gzip.open if str(fq).endswith(".gz") else open
    lines = 0
    with opener(fq, "rt") as fh:
        for _ in fh:
            lines += 1
            if lines >= k * 4:
                return True
    return False


def is_already_target(sqs: list | None, target: str) -> bool:
    """True when the input is ALREADY aligned to the requested target → realigning is wasted work."""
    ref = classify_reference(sqs)
    return (target == "kirby" and ref == "kirby") or (target == "t2t" and ref == "t2t_win")


def chr1_contig_name(sqs: list | None) -> str | None:
    """The actual chr1 contig name in this header — `chr1`, `1`, or a RefSeq/GenBank id.
    Resolved by name first, then by the GRCh38 chr1 LENGTH so a `1`-named BAM still matches."""
    if not sqs:
        return None
    names = {n for n, _ in sqs}
    for cand in ("chr1", "1", "CM000663.2", "NC_000001.11"):
        if cand in names:
            return cand
    for n, ln in sqs:
        if ln in _HG38_CHR1:
            return n
    return None


def t2t_chr1_contig_name(sqs: list | None) -> str | None:
    """The chr1 contig name of a full-genome T2T-CHM13 header — a RefSeq/GenBank id (NC_060925.1 /
    CP068277.2) or `chr1`/`1`. The RefSeq/GenBank ids are tried FIRST so a header that mixes an hg38
    `chr1` with a T2T RefSeq contig still resolves to the T2T one; fall back to the CHM13 chr1 LENGTH."""
    if not sqs:
        return None
    names = {n for n, _ in sqs}
    for cand in ("NC_060925.1", "CP068277.2", "chr1", "1"):
        if cand in names:
            return cand
    for n, ln in sqs:
        if ln in _T2T_CHR1:
            return n
    return None


def extraction_region(sqs: list | None) -> str | None:
    """Region to restrict an aligned input to before re-extracting reads.

    None means "take every read": for a uBAM/fastq there is nothing to restrict, and for an
    unknown reference guessing a window would silently drop reads. For hg38 / full-genome T2T the
    window uses the header's REAL chr1 name (not a hardcoded `chr1`) so a `1`-named BAM does not
    extract 0 reads. A full-genome T2T input (ref == 't2t') MUST be windowed here too: otherwise it
    falls through to None and every read is kept — the same bloat as a whole-run fastq.
    """
    ref = classify_reference(sqs)
    if ref == "hg38":
        return f"{chr1_contig_name(sqs) or 'chr1'}:{MUC1_HG38_COORDS}"
    if ref == "t2t":
        return f"{t2t_chr1_contig_name(sqs) or 'chr1'}:{MUC1_T2T_COORDS}"
    if ref == "t2t_win":
        return MUC1_T2T_WINDOW
    return None


def minimap_preset(pacbio: bool = False) -> str:
    """`map-hifi` for PacBio HiFi, else `map-ont`. The only platform-dependent knob."""
    return "map-hifi" if pacbio else "map-ont"


def minimap_argv(ref: str, fq: str, *, threads: int = 8, pacbio: bool = False, pcr: bool = False) -> list:
    """The minimap2 command line (pure → unit-testable). `-y` carries MM/ML/HP/PS tags across the realign;
    `--secondary=no` for the ~150 near-identical VNTR contigs. `pcr=True` appends the FROZEN LR-PCR
    chaining flags (`config.PCR_MINIMAP_EXTRA` = `-z 600,200 -r 2000,20000`, set B) that let a VNTR-spanning
    amplicon read chain across the long tandem — the default bandwidth fragments/loses it on a repeat."""
    preset = minimap_preset(pacbio)
    extra = list(PCR_MINIMAP_EXTRA) if pcr else []
    return ["minimap2", "-ax", preset, "-y", "--secondary=no", *extra, "-t", str(threads), ref, fq]


# ── plumbing ──────────────────────────────────────────────────────────────────
def _header_sqs(path: str, ref: str | None = None) -> list | None:
    """[(contig, length)] from a BAM/CRAM header; None for fastq or an unreadable header."""
    if str(path).endswith(FASTQ_EXTS):
        return None
    import pysam
    mode = "rc" if str(path).endswith(".cram") else "rb"
    kw = {"reference_filename": ref} if (mode == "rc" and ref) else {}
    try:
        with pysam.AlignmentFile(path, mode, **kw) as af:
            return list(zip(af.references, af.lengths))
    except Exception:
        return None


def _require_tools(*tools: str) -> None:
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        raise SystemExit(
            f"[prepare] missing external tool(s): {', '.join(missing)}. "
            "Install them (conda env: `conda env create -f environment.yml`) or "
            "`module load samtools minimap2` on an HPC.")


def _samtools_view(inp: str, ref38: str | None):
    v = ["samtools", "view", "-b"]
    if str(inp).endswith(".cram") and ref38:
        v += ["-T", ref38]
    return v


def clean_fastq_header(header: str) -> str:
    """Keep a fastq header's read id + ONLY the whitelisted SAM tags (MM/ML/HP/PS), dropping any other
    comment field. PURE / unit-tested. A RAW fastq — especially PacBio HiFi — often carries comment fields
    that are NOT valid SAM tags (`ccs`, `np:i:..`, a `zm:4:..` with a digit type, bare coordinates).
    `minimap2 -y` copies the comment verbatim as SAM aux fields, and samtools then ABORTS the whole pipe
    (`[E::aux_parse] unrecognized type '4'` / `Incomplete aux field`). Stripping to the whitelist makes -y
    safe on any fastq while still preserving ONT methylation/phase (dorado writes MM/ML/HP/PS here)."""
    parts = header.split()
    if not parts:
        return header
    kept = [parts[0]] + [t for t in parts[1:] if t[:3] in ("MM:", "ML:", "HP:", "PS:")]
    return " ".join(kept)


def sanitize_fastq_comments(inp: str, out_fq: str) -> str:
    """Stream a (optionally gzipped) raw fastq to `out_fq`, cleaning every header via clean_fastq_header so
    the downstream `minimap2 -y` cannot be broken by a non-SAM-tag comment. Returns out_fq."""
    opn = gzip.open if str(inp).endswith(".gz") else open
    with opn(inp, "rt") as fi, open(out_fq, "w") as fo:
        for i, line in enumerate(fi):
            fo.write((clean_fastq_header(line.rstrip("\n")) + "\n") if i % 4 == 0 else line)
    return out_fq


_SAM_AUX = re.compile(r"^[A-Za-z][A-Za-z0-9]:[AifcCsSIZHB]:")   # a valid SAM aux field TAG:TYPE:…


def fastq_comment_is_clean(inp: str, k: int = 20) -> bool:
    """True when the first k reads carry NO comment or ONLY valid SAM-tag comments → `minimap2 -y` is safe
    and NO sanitize rewrite is needed (the fast path). ONT dorado (`MM:Z:`/`ML:B:`) and tag-free fastqs stay
    on the fast path; a PacBio comment (`ccs`, `rq:0.998`, `zm:4:`) trips this → sanitize. Peeking ~20 reads
    is O(1); it avoids a full-file copy on every ONT run (a perf regression the blanket sanitize introduced)."""
    opn = gzip.open if str(inp).endswith(".gz") else open
    try:
        with opn(inp, "rt") as fh:
            seen = 0
            for i, line in enumerate(fh):
                if i % 4 != 0:
                    continue
                for tok in line.split()[1:]:
                    if not _SAM_AUX.match(tok):
                        return False
                seen += 1
                if seen >= k:
                    break
    except Exception:
        return False                                           # unreadable → be safe, sanitize
    return True


def to_fastq(inp: str, out_fq: str, *, region: str | None = None, ref38: str | None = None) -> str:
    """Reads → fastq, PRESERVING MM/ML methylation tags. `region` restricts an aligned input."""
    _require_tools("samtools")
    if str(inp).endswith(FASTQ_EXTS):
        # A raw fastq is not round-tripped through samtools (which would clean the comment), so an
        # arbitrary comment would reach `minimap2 -y` and abort samtools. Only rewrite when a peek shows
        # a non-SAM-tag comment (PacBio); a clean/tag-free ONT fastq is used AS-IS (no full-file copy).
        return inp if fastq_comment_is_clean(inp) else sanitize_fastq_comments(inp, out_fq)
    view = _samtools_view(inp, ref38) + [inp] + ([region] if region else [])
    with open(out_fq, "wb") as fh:
        if region:
            p1 = subprocess.Popen(view, stdout=subprocess.PIPE)
            subprocess.run(["samtools", "fastq", "-T", FASTQ_CARRY_TAGS, "-"],
                           stdin=p1.stdout, stdout=fh, check=True)
            p1.wait()
        else:
            subprocess.run(["samtools", "fastq", "-T", FASTQ_CARRY_TAGS, inp], stdout=fh, check=True)
    if os.path.getsize(out_fq) == 0:
        raise SystemExit(f"[prepare] no reads extracted from {inp}"
                         + (f" in {region}" if region else "")
                         + " — wrong region, or the input has no MUC1 coverage.")
    return out_fq


def append_unmapped_fastq(inp: str, out_fq: str, *, ref38: str | None = None) -> str:
    """APPEND the unplaced reads (`RNAME=*`) of an aligned input to an existing fastq. On hg38 the VNTR is
    collapsed, so VNTR-heavy reads cannot be placed and are flagged unmapped — the window fetch misses them
    (validated: 264 window reads → +9060 recovered). The VNTR reference filters non-MUC1 at the realign
    (≈98 % of the unplaced pool drops, ~0 contamination). Needs the input indexed (the `*` fetch seeks the
    unplaced block). Kept SEPARATE from the window so the auto-PCR read count keys off the window alone."""
    _require_tools("samtools")
    with open(out_fq, "ab") as fh:
        p = subprocess.Popen(_samtools_view(inp, ref38) + [inp, "*"], stdout=subprocess.PIPE)
        subprocess.run(["samtools", "fastq", "-T", FASTQ_CARRY_TAGS, "-"],
                       stdin=p.stdout, stdout=fh, check=True)
        p.wait()
    return out_fq


def align(fq: str, ref: str, out_bam: str, *, threads: int = 8, pacbio: bool = False,
          pcr: bool = False) -> str:
    """fastq → sorted+indexed BAM on `ref`, keeping methylation tags (`-y`), mapped reads only.
    `pcr=True` applies the frozen LR-PCR chaining flags (see `minimap_argv`)."""
    _require_tools("minimap2", "samtools")
    argv = minimap_argv(ref, fq, threads=threads, pacbio=pacbio, pcr=pcr)
    print(f"[prepare] {' '.join(argv[:-2])} → {out_bam}", file=sys.stderr)
    # --secondary=no: the VNTR reference is ~150 NEAR-IDENTICAL length-contigs, so a read scores almost
    # equally on all of them → default minimap2 emits ~150 SECONDARY alignments/read, bloating the BAM
    # (observed 2.3 GB) with records NOTHING downstream reads (count/consensus/arbiter all exclude
    # secondary). Suppress them at the source; also drop any secondary/supplementary at `view` (-F 0x904).
    p1 = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(["samtools", "view", "-b", "-F", "0x904", "-"],
                          stdin=p1.stdout, stdout=subprocess.PIPE)
    subprocess.run(["samtools", "sort", "-@", "4", "-o", out_bam, "-"], stdin=p2.stdout, check=True)
    p1.wait(); p2.wait()
    subprocess.run(["samtools", "index", out_bam], check=True)
    return out_bam


def _cap_fastq(fq: str, cap: int, work: str, tag: str = "pcr") -> str:
    """Subsample a (plain) fastq to at most `cap` reads → bound minimap2 on a DEEP amplicon vs the ~150-contig
    reference (chaining × ~150 near-identical contigs stalls the aligner). Reservoir sample (fixed seed →
    reproducible) so the allele RATIO is preserved and the PCR-depleted long allele survives at a
    representative fraction. Returns `fq` unchanged when it already holds ≤ cap reads."""
    import random
    with open(fq) as fh:
        n = sum(1 for _ in fh) // 4
    if n <= cap:
        return fq
    rng = random.Random(1234)
    keep, idx, rec = [], 0, []
    with open(fq) as fh:
        for line in fh:
            rec.append(line)
            if len(rec) == 4:
                if len(keep) < cap:
                    keep.append(rec)
                else:
                    j = rng.randint(0, idx)
                    if j < cap:
                        keep[j] = rec
                idx += 1
                rec = []
    out = os.path.join(work, f"{tag}.capped.fq")
    with open(out, "w") as w:
        for r in keep:
            w.writelines(r)
    print(f"[prepare] --pcr: capped {n} → {cap} reads before alignment (a deep amplicon vs the multi-contig "
          f"reference stalls minimap2; the allele ratio is preserved). Raise/lower with --pcr-cap.", file=sys.stderr)
    return out


def prepare(inp: str, vntr_ref: str, out_bam: str, *, ref38: str | None = None,
            threads: int = 8, pacbio: bool = False, pcr: bool = False,
            no_pcr_autodetect: bool = False, recover_unmapped: bool = True,
            workdir: str | None = None, sample: str | None = None, pcr_cap: int | None = None) -> str:
    """Any MUC1 input → a BAM aligned on the VNTR reference, ready for `call`.

    Returns the path to the ready BAM (`inp` itself when it is already on the target).
    """
    sqs = _header_sqs(inp, ref38)
    kind = detect_input_kind(inp, sqs)
    refkind = classify_reference(sqs)

    if is_already_target(sqs, "kirby"):
        print(f"[prepare] input is ALREADY aligned to the VNTR reference ({refkind}) → bypass",
              file=sys.stderr)
        return inp

    region = extraction_region(sqs)
    print(f"[prepare] input={kind} reference={refkind} "
          f"extract={region or 'ALL reads (uBAM/fastq/unknown → nothing dropped)'}", file=sys.stderr)
    # RECOVER-UNMAPPED: DEFAULT-ON for any ALIGNED input (region present) — hg38 COLLAPSES the VNTR so
    # VNTR-heavy reads are flagged unmapped and the window fetch misses them; recovering them lifts a thin
    # locus by ~35× (validated: 264 → 9324 reads). The VNTR reference filters non-MUC1 at the realign. Same
    # philosophy as auto-PCR (optimise per input type). Disable with --no-recover-unmapped.
    # NOTE: on a full-genome T2T input the VNTR is faithfully REPRESENTED, so its reads already map inside
    # the window and recovery is low-yield (harmless — the VNTR ref still filters — but adds the whole
    # unplaced-pool pass); drop --no-recover-unmapped to skip it on a deep T2T WGS.
    recover = recover_unmapped and region is not None
    if recover:
        print("[prepare] recover-unmapped ON (aligned input): appending the unplaced (hg38-unmappable VNTR) "
              "reads; the VNTR reference filters non-MUC1 at the realign. Disable with --no-recover-unmapped. "
              "⚠ on a deep WGS the unplaced pool can be millions of reads (a few extra minutes).", file=sys.stderr)
    elif kind == "aligned" and refkind == "hg38":
        print("[prepare] NOTE: recover-unmapped is OFF — hg38 collapses the VNTR, so the window alone misses "
              "the VNTR-heavy (unmapped) reads. Drop --no-recover-unmapped to pull them back.", file=sys.stderr)

    work = workdir or os.path.dirname(os.path.abspath(out_bam)) or "."
    os.makedirs(work, exist_ok=True)
    fq = os.path.join(work, f"{sample or 'sample'}.muc1.fq")
    fq = to_fastq(inp, fq, region=region, ref38=ref38)     # window (or ALL for uBAM/fastq)

    # AUTO-amplicon: apply the frozen chaining flags when the sample is amplicon-scale, even if the user
    # forgot --pcr. Count the reads that actually ALIGN to the VNTR reference (NOT the raw 20 kb window,
    # which also holds flanking/off-target reads → a busy AS window would be mis-read as PCR). Detect on the
    # WINDOW, before recovery. Only run the detection align when auto-detection is genuinely in play.
    cap = pcr_cap if pcr_cap is not None else PCR_ALIGN_CAP
    if not pacbio and not pcr and not no_pcr_autodetect:
        # detect on a CAPPED sample so a deep amplicon doesn't stall the detection align too (≥ threshold in
        # the sample ⇒ amplicon-scale; the cap ≫ threshold so the decision holds).
        n_vntr = vntr_aligned_count(_cap_fastq(fq, cap, work, (sample or "s") + ".detect"),
                                    vntr_ref, threads=threads, pacbio=pacbio)
        if autopcr_decision(amplicon_scale=n_vntr >= PCR_AUTODETECT_MIN_READS,
                            explicit_pcr=False, disabled=False):
            pcr = True
            print(f"[prepare] ⚠ AUTO-amplicon: {n_vntr} reads align to the MUC1 VNTR (≥ "
                  f"{PCR_AUTODETECT_MIN_READS}) → amplicon-scale → applying the frozen chaining flags "
                  f"(set B). Add --pcr to also cap the `call` depth, or --no-pcr-autodetect to disable.",
                  file=sys.stderr)

    if recover:                                            # append the hg38-unmappable VNTR reads AFTER detection
        append_unmapped_fastq(inp, fq, ref38=ref38)
    if pcr:                                                # deep amplicon → bound minimap2 (preserves allele ratio)
        fq = _cap_fastq(fq, cap, work, sample or "pcr")
    out = align(fq, vntr_ref, out_bam, threads=threads, pacbio=pacbio, pcr=pcr)
    report_amplicon(out, is_pcr=pcr, pacbio=pacbio)
    return out


def report_amplicon(out_bam: str, *, is_pcr: bool, pacbio: bool) -> None:
    """Self-configuring log line: characterise the amplicon from the aligned reads — kind (native LR-PCR vs
    foreign/short PCR, from our AL/AH flank anchors), primer flank sizes, a platform guess (median base
    quality), and which VNTR-length method fits. Best-effort; never raises into the pipeline."""
    try:
        import vntr_raw_length as V
        sig = V.amplicon_signature(out_bam)
        if not sig.get("n_spanning"):
            return
        plat = "PacBio HiFi" if pacbio else sig.get("platform", "unknown")   # --pacbio wins over the Q heuristic
        lenmethod = ("flank arbiter (native)" if sig["kind"] == "native_LR-PCR"
                     else "cassette 1→9 (amplicon-agnostic)")
        prod = ""
        if sig.get("product_bp"):
            prod = " — product " + "/".join(f"~{round(p/1000,1)}kb" for p in sig["product_bp"]) + " (flanks+VNTR)"
        tag = "PCR amplicon" if is_pcr else "MUC1 reads"
        print(f"[prepare] {tag} detected: {sig['kind']}{prod} — flanks 5'≈{sig['flank5_med']}/"
              f"3'≈{sig['flank3_med']} bp (AL {sig['has_AL_frac']*100:.0f}%/AH {sig['has_AH_frac']*100:.0f}%), "
              f"platform≈{plat}. VNTR length → {lenmethod}.", file=sys.stderr)
    except Exception:
        pass


def sample_read_names(inp: str, ref38: str | None = None, k: int = 200) -> list:
    """Up to k read names from the input (fastq(.gz) / uBAM / BAM / CRAM) — enough to guess the platform."""
    names = []
    if str(inp).endswith(FASTQ_EXTS):
        opn = gzip.open if str(inp).endswith(".gz") else open
        with opn(inp, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0 and line:
                    names.append(line[1:].split()[0] if line[0] == "@" else line.split()[0])
                if len(names) >= k:
                    break
        return names
    try:
        import pysam
        mode = "rc" if str(inp).endswith(".cram") else "rb"
        kw = {"reference_filename": ref38} if (mode == "rc" and ref38) else {}
        with pysam.AlignmentFile(inp, mode, **kw) as af:
            for r in af.fetch(until_eof=True):
                if r.query_name:
                    names.append(r.query_name)
                if len(names) >= k:
                    break
    except Exception:
        pass
    return names


def detect_platform(inp: str, ref38: str | None = None, k: int = 200) -> str | None:
    """'PacBio HiFi' | 'ONT' | None from the input read names (majority vote). Read names are written by the
    instrument (ONT = a UUID, PacBio CCS = a movie/zmw/ccs name), so this works for PCR, AS and WGS alike —
    it identifies the SEQUENCER, letting prepare pick minimap2 map-ont vs map-hifi when the user forgot to."""
    import vntr_raw_length as V
    from collections import Counter
    votes = Counter(v for v in (V.platform_from_name(n) for n in sample_read_names(inp, ref38, k)) if v)
    return votes.most_common(1)[0][0] if votes else None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="muc1_analyzer prepare",
        description="Extract MUC1 reads from any input and align them to the VNTR reference "
                    "(methylation tags preserved). Auto-detects fastq / uBAM / hg38 / T2T and "
                    "bypasses an input already aligned to the VNTR reference.")
    ap.add_argument("-i", "--input", required=True,
                    help="fastq(.gz) | uBAM (unaligned) | aligned BAM/CRAM (hg38 or T2T)")
    ap.add_argument("-r", "--vntr-ref", dest="vntr_ref", required=True,
                    help="multi-contig VNTR reference FASTA (MUC1_fakedVNTR1to150revcomplKirby.fa)")
    ap.add_argument("-o", "--out", required=True, help="output BAM (sorted + indexed)")
    ap.add_argument("-s", "--sample", default=None, help="sample name (temp-file prefix)")
    ap.add_argument("--ref", dest="ref38", default=None,
                    help="genome FASTA — REQUIRED to read a CRAM input")
    ap.add_argument("--ont", action="store_true", help="Oxford Nanopore reads (default)")
    ap.add_argument("--pacbio", action="store_true", help="PacBio HiFi reads (minimap2 map-hifi)")
    ap.add_argument("--pcr", action="store_true",
                    help="LR-PCR amplicon input: apply the frozen chaining flags (-z 600,200 -r 2000,20000) "
                         "so a VNTR-spanning amplicon read chains across the tandem (ONT; config set B). "
                         "Auto-enabled when the locus read count is amplicon-scale unless --no-pcr-autodetect.")
    ap.add_argument("--no-pcr-autodetect", action="store_true",
                    help="disable the read-count LR-PCR auto-detection (keep default map-ont chaining)")
    ap.add_argument("--recover-unmapped", action="store_true",
                    help="(DEFAULT for aligned inputs) also pull the UNPLACED reads (VNTR-heavy reads hg38/T2T "
                         "could not place) → the VNTR ref filters them → recovers a thin locus (~35x).")
    ap.add_argument("--no-recover-unmapped", action="store_true",
                    help="disable the default unplaced-read recovery on an aligned input")
    ap.add_argument("--pcr-cap", type=int, default=None,
                    help=f"--pcr: subsample a deep amplicon to this many reads before minimap2 "
                         f"(default {PCR_ALIGN_CAP}; preserves the allele ratio; prevents the multi-contig stall)")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args(argv)

    if args.ont and args.pacbio:
        print("[prepare] --ont and --pacbio are mutually exclusive", file=sys.stderr)
        return 2
    if args.pcr and args.pacbio:
        print("[prepare] --pcr (ONT LR-PCR, frozen set B) and --pacbio are mutually exclusive",
              file=sys.stderr)
        return 2
    if str(args.input).endswith(".cram") and not args.ref38:
        print("[prepare] --ref (genome FASTA) is required to read a CRAM", file=sys.stderr)
        return 2

    # Auto-pick the sequencer preset when the user gave NO platform/mode flag: read names identify PacBio
    # (…/ccs) vs ONT (UUID) for PCR / AS / WGS alike, so the right minimap2 preset is used even if forgotten.
    if not args.ont and not args.pacbio and not args.pcr:
        plat = detect_platform(args.input, args.ref38)
        if plat == "PacBio HiFi":
            args.pacbio = True
            print("[prepare] auto-detected PacBio HiFi (read names …/ccs) → minimap2 map-hifi "
                  "(pass --ont/--pacbio to override).", file=sys.stderr)
        elif plat == "ONT":
            print("[prepare] auto-detected ONT (UUID read names) → minimap2 map-ont.", file=sys.stderr)

    out = prepare(args.input, args.vntr_ref, args.out, ref38=args.ref38,
                  threads=args.threads, pacbio=args.pacbio, pcr=args.pcr,
                  no_pcr_autodetect=args.no_pcr_autodetect,
                  recover_unmapped=not args.no_recover_unmapped, sample=args.sample,
                  pcr_cap=args.pcr_cap)
    print(f"[prepare] ready → {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
