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

from .config import (PCR_ALIGN_CAP, PCR_AUTODETECT_MIN_READS, PCR_MINIMAP_EXTRA,
                     PCR_MINIMAP_PROFILES, resolve_pcr_profile)

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
    rc = p1.wait()
    # Same trap as `align`: a KILLED aligner counts 0, and 0 reads reads as "not amplicon-scale" — the run
    # would then pick the wrong chaining profile for a reason nobody could see. Fail here instead; the real
    # alignment is about to die anyway, and an early message names the cause.
    if rc:
        raise SystemExit(f"[prepare] the amplicon-detection alignment FAILED (minimap2 rc={rc}) — a read "
                         f"count of 0 would be indistinguishable from a non-amplicon sample. rc=-9 means "
                         f"the aligner was KILLED (memory, or a shared front-end's cgroup).")
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


def minimap_argv(ref: str, fq: str, *, threads: int = 8, pacbio: bool = False, pcr: bool = False,
                 profile: str = "setB") -> list:
    """The minimap2 command line (pure → unit-testable). `-y` carries MM/ML/HP/PS tags across the realign;
    `--secondary=no` for the ~150 near-identical VNTR contigs. `pcr=True` appends the LR-PCR flags of
    `profile` (default `setB` = the FROZEN `-z 600,200 -r 2000,20000`, which lets a VNTR-spanning amplicon
    read chain across the long tandem — the default bandwidth fragments/loses it on a repeat).

    `profile` selects a candidate parameter set from `config.PCR_MINIMAP_PROFILES` so an alternative can be
    PRICED on data rather than argued about. It changes nothing unless asked for by name, and it is a no-op
    outside the LR-PCR route. A profile may also override the platform preset (the MUC1 set runs `lr:hq`)."""
    _canonical, prof = resolve_pcr_profile(profile)   # RAISES on an unknown name (no silent setB fallback)
    preset = (prof.get("preset") if pcr else None) or minimap_preset(pacbio)
    extra = list(prof.get("extra") or PCR_MINIMAP_EXTRA) if pcr else []
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


def _preset_platform_mismatch(preset, detected):
    """Is `preset` committing to the OTHER platform than the reads' names indicate? Pure (testable).

    Only the platform-committing presets can mismatch: `map-ont` on PacBio reads (the exact
    HiFi-aligned-as-nanopore class) or `map-hifi` on ONT reads. `lr:hq` (the muc1 set) is
    high-quality-agnostic, and an absent signal (None) is never a mismatch."""
    return ((preset == "map-ont" and detected == "PacBio HiFi") or
            (preset == "map-hifi" and detected == "ONT"))


def align(fq: str, ref: str, out_bam: str, *, threads: int = 8, pacbio: bool = False,
          pcr: bool = False, profile: str = "setB") -> str:
    """fastq → sorted+indexed BAM on `ref`, keeping methylation tags (`-y`), mapped reads only.
    `pcr=True` applies the LR-PCR flags of `profile` (default: the frozen set B — see `minimap_argv`)."""
    _require_tools("minimap2", "samtools")
    argv = minimap_argv(ref, fq, threads=threads, pacbio=pacbio, pcr=pcr, profile=profile)
    print(f"[prepare] {' '.join(argv[:-2])} → {out_bam}", file=sys.stderr)
    # ── platform × preset guard (non-fatal) — catch a HiFi library aligned with a nanopore preset (or
    # the reverse) even when the profile/flags were chosen deliberately. Read names identify the
    # sequencer; if they contradict a platform-committing preset, warn and PROCEED (never change what
    # was asked). lr:hq (muc1) is agnostic → not checked; stripped names → no signal → silent.
    _preset = argv[argv.index("-ax") + 1] if "-ax" in argv else None
    if _preset in ("map-ont", "map-hifi"):
        try:
            _seen = detect_platform(fq)
        except Exception:
            _seen = None
        if _preset_platform_mismatch(_preset, _seen):
            _hint = "pass --pacbio" if _seen == "PacBio HiFi" else "drop --pacbio"
            print(f"[prepare] ⚠ read names look {_seen} but aligning with -ax {_preset} — {_hint} "
                  f"(or --align-profile muc1) for the matching preset. Proceeding as requested.",
                  file=sys.stderr)
    # --secondary=no: the VNTR reference is ~150 NEAR-IDENTICAL length-contigs, so a read scores almost
    # equally on all of them → default minimap2 emits ~150 SECONDARY alignments/read, bloating the BAM
    # (observed 2.3 GB) with records NOTHING downstream reads (count/consensus/arbiter all exclude
    # secondary). Suppress them at the source; also drop any secondary/supplementary at `view` (-F 0x904).
    #
    # ⚠ A DEAD ALIGNER MUST NOT READ AS A SAMPLE WITHOUT COVERAGE. minimap2's stderr used to go to
    # DEVNULL and `p1.wait()`'s status was discarded, so nothing here could fail: `samtools sort` on an
    # empty stream exits 0 and writes a valid, EMPTY BAM. Measured 2026-09-03 on a login node — the
    # cgroup killed minimap2 mid-run and the pipeline carried on to a silent `carrier=False`. The
    # empty-input guard in `to_fastq` sits one step earlier and does not cover this.
    # stderr goes to a FILE, not a pipe: minimap2 is chatty and a pipe nobody drains would deadlock it.
    err = out_bam + ".minimap2.err"
    with open(err, "wb") as efh:
        p1 = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=efh)
        p2 = subprocess.Popen(["samtools", "view", "-b", "-F", "0x904", "-"],
                              stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()   # so a dying p2 SIGPIPEs the aligner instead of leaving it to fill a pipe
        subprocess.run(["samtools", "sort", "-@", "4", "-o", out_bam, "-"], stdin=p2.stdout, check=True)
        p2.stdout.close()
        rc1, rc2 = p1.wait(), p2.wait()
    if rc1 or rc2:
        raise SystemExit(f"[prepare] ALIGNMENT FAILED (minimap2 rc={rc1}, samtools view rc={rc2}) — "
                         f"{out_bam} holds whatever was written before the failure and MUST NOT be called "
                         f"on. A negative rc is a signal, not an error code: rc=-9 means the aligner was "
                         f"KILLED, out of memory or by a shared front-end's cgroup — rerun on a compute "
                         f"node.\n{_tail(err)}")
    try:
        os.remove(err)      # kept only when it has something to say
    except OSError:
        pass
    subprocess.run(["samtools", "index", out_bam], check=True)
    return out_bam


def _tail(path: str, n: int = 8) -> str:
    """Last `n` lines of a file, best-effort — for putting a failed tool's own words in our message."""
    try:
        with open(path, errors="replace") as fh:
            return "".join(fh.readlines()[-n:]).strip()
    except OSError:
        return ""


def aligned_read_count(bam: str) -> int:
    """Primary aligned records in `bam` (-F 0x904), or -1 when it cannot be counted."""
    try:
        out = subprocess.run(["samtools", "view", "-c", "-F", "0x904", bam],
                             stdout=subprocess.PIPE, text=True, check=True).stdout
        return int(out.strip() or 0)
    except (subprocess.CalledProcessError, OSError, ValueError):
        return -1


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
            workdir: str | None = None, sample: str | None = None, pcr_cap: int | None = None,
            two_stage: bool = False, align_profile: str = "setB") -> str:
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
    out = None
    if two_stage:
        out, info = align_two_stage(fq, vntr_ref, out_bam, threads=threads, pacbio=pacbio, pcr=pcr,
                                    profile=align_profile,
                                    work=work, sample=sample or "s")
        if out is None:
            print(f"[prepare] two-stage align unavailable ({info}) → full {os.path.basename(vntr_ref)}",
                  file=sys.stderr)
        else:
            print(f"[prepare] two-stage align: ruler → alleles {info['alleles']} "
                  f"({info['method']}) → {info['contigs_aligned']} contig alignments instead of the full "
                  f"reference", file=sys.stderr)
    if out is None:
        out = align(fq, vntr_ref, out_bam, threads=threads, pacbio=pacbio, pcr=pcr,
                    profile=align_profile)
    # Belt to the return-code brace: refuse to hand `call` a BAM with nothing in it. Reads WERE extracted
    # (`to_fastq` guarantees it), so zero primary alignments here means the aligner produced nothing usable
    # — not a biological result. Checked on the DELIVERED bam only: `align_two_stage`'s internal ruler pass
    # is allowed to come back empty, that is exactly how it decides to fall back to the full reference.
    n = aligned_read_count(out)
    if n == 0:
        raise SystemExit(f"[prepare] {out} contains ZERO aligned reads while {fq} was non-empty — the "
                         f"alignment produced nothing. Calling on this file would report a NON-CARRIER "
                         f"for a sample that was never measured. Check memory/tooling, then rerun.")
    report_amplicon(out, is_pcr=pcr, pacbio=pacbio)
    return out


def subset_ref(vntr_ref: str, contigs, out_fa: str) -> str | None:
    """Write a mini-reference holding only `contigs`. Impure, returns None when none could be extracted.

    The point of the two-stage align: minimap2's cost scales with the reference, and 150 near-identical
    VNTR contigs is 150x the work of the two that actually matter. Everything downstream keeps working on
    a subset — `find_rs4072037_offset` is a CONSTANT offset shared by every contig, and the length arbiter
    reads sequences, not contig names."""
    import pysam
    fa = pysam.FastaFile(vntr_ref)
    have = set(fa.references)
    keep = [c for c in contigs if c in have]
    if not keep:
        return None
    with open(out_fa, "w") as fh:
        for c in keep:
            seq = fa.fetch(c)
            fh.write(f">{c}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")
    pysam.faidx(out_fa)
    return out_fa


def align_two_stage(fq: str, vntr_ref: str, out_bam: str, *, threads: int = 8, pacbio: bool = False,
                    pcr: bool = False, work: str = ".", sample: str = "s",
                    ruler: str = "MUC1_VNTR_150repeats", window: int = 5,
                    profile: str = "setB"):
    """Align to ONE ruler contig to measure, then realign to the TWO allele contigs. Impure.

    Measured motivation: aligning an LR-PCR amplicon against the full 150-contig reference took 43 min for
    one sample, while the streaming AS pipeline did 20 GB in 25 min — because it aligns to a single contig.
    The reference, not the read count, is the cost.

    Everything the single-stage path provides is preserved: the same preset (ONT/PacBio), the same frozen
    LR-PCR chaining flags, `-y` tag carry-through, and the amplicon signature / strand balance, which are
    computed from read SEQUENCES and flags rather than from the number of contigs.

    Returns (bam, info) or (None, reason) when the measurement stage cannot call alleles — the caller then
    falls back to the full reference rather than guessing."""
    import vntr_raw_length as V
    from .config import FLANK_TO_PHYSICAL_OFFSET

    ruler_fa = subset_ref(vntr_ref, [ruler], os.path.join(work, f"{sample}.ruler.fa"))
    if ruler_fa is None:
        return None, f"ruler contig {ruler} absent from the reference"
    ruler_bam = os.path.join(work, f"{sample}.ruler.bam")
    align(fq, ruler_fa, ruler_bam, threads=threads, pacbio=pacbio, pcr=pcr, profile=profile)

    arb = V.auto_length(ruler_bam, chrom=None, offset=FLANK_TO_PHYSICAL_OFFSET)
    alleles = arb.get("alleles") if arb.get("available") else None
    if not alleles:
        return None, "the ruler stage could not call allele lengths"

    from .allele_scaffold import CONTIG_LEN
    import pysam
    have = {}
    for c in pysam.FastaFile(vntr_ref).references:
        m = CONTIG_LEN.search(c)
        if m:
            have[int(m.group(1))] = c
    # A WINDOW around each allele, not a single contig. Measured: forcing every read onto exactly two
    # contigs sent the reads that belong at intermediate lengths onto the nearest allele — the 45-contig
    # took 1378 reads instead of ~500, its consensus mixed both alleles, and the dupC insertion fraction
    # fell under the calling threshold. The variant was lost while the LENGTHS were right. Reads must be
    # free to settle at their own length; the speed comes from dropping 150 contigs to ~15, not to 2.
    #
    # ⚠ This window must be >= `allele_scaffold.contigs_for_lengths`' own window (5). They model the SAME
    # thing — the smear of ONT length calls across neighbouring contigs — and the scaffold step picks the
    # read MODE within its window. A narrower window here silently deletes the very contig that step would
    # have chosen: measured on SER, the allele's 139 reads sit 4 copies from the called length, which a
    # window of 3 would have excluded from the mini reference altogether.
    want = []
    for a in alleles:
        near = min(have, key=lambda k: (abs(k - int(a)), k))
        for L in sorted(k for k in have if abs(k - near) <= window):
            if have[L] not in want:
                want.append(have[L])
    mini = subset_ref(vntr_ref, want, os.path.join(work, f"{sample}.alleles.fa"))
    if mini is None:
        return None, f"none of {want} is in the reference"
    align(fq, mini, out_bam, threads=threads, pacbio=pacbio, pcr=pcr, profile=profile)
    return out_bam, {"ruler": ruler, "alleles": alleles, "scaffolds": want,
                     "contigs_aligned": 1 + len(want), "method": arb.get("method")}


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
    ap.add_argument("--two-stage", action="store_true",
                    help="align to ONE ruler contig to measure, then realign to the TWO allele contigs, "
                        "instead of the full 150-contig reference. minimap2's cost scales with the "
                        "reference, not the read count: 43 min measured on one amplicon against 150 "
                        "contigs. Everything else is unchanged (preset, frozen PCR flags, tag carry-"
                        "through, amplicon signature, strand balance).")
    ap.add_argument("--align-profile", dest="align_profile", default="setB",
                    choices=sorted(PCR_MINIMAP_PROFILES) + ["xavier"],  # xavier = deprecated alias → muc1
                    help="LR-PCR minimap2 parameter set (default setB = the shipped frozen flags). "
                         "`muc1` = the MUC1 production set VERBATIM, `muc1_scoring` = the MUC1 scoring with OUR "
                         "chaining, so the anti-slide bias and the -z/-r bandwidth can be priced separately. "
                         "No effect outside the LR-PCR route.")
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
                  pcr_cap=args.pcr_cap, two_stage=getattr(args, 'two_stage', False),
                  align_profile=getattr(args, 'align_profile', 'setB'))
    print(f"[prepare] ready → {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
