"""Entry-point input policy (v1.0.1): MUC1_Analyzer reads FASTQ or unaligned BAM (uBAM), nothing else.

An ALIGNED BAM/CRAM used to be accepted and re-extracted. Both external user reports of 2026-11 came from
that path: the extraction from a GRCh38 CRAM pulled every unplaced read of the file (81 235 rejected
adaptive-sampling fragments buried 98 real locus reads; a 569 MB fastq on a WGS CRAM), and the result
depended on how the CRAM had been made. The entry points (`run`, `run_pcr`, `prepare`, `from-fastq`)
now refuse an aligned input with the command that turns it into a locus FASTQ. `call` keeps taking the
VNTR-reference BAM it is designed for; it is an internal step, not an entry point."""
from __future__ import annotations

import os

FASTQ_EXTS = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
ALIGNMENT_EXTS = (".bam", ".ubam", ".cram", ".sam")

EXTRACT_HINT = (
    "Give the MUC1 reads as FASTQ or as an unaligned BAM (uBAM). From an aligned BAM/CRAM:\n"
    "  samtools view -b -F 0x904 [-T genome.fa] INPUT.cram chr1:155180000-155200000 \\\n"
    "    | samtools fastq -T MM,ML - > SAMPLE.MUC1.fq\n"
    "(chr1:155180000-155200000 is the GRCh38 MUC1 window; -T is needed for a CRAM; -T MM,ML keeps the "
    "methylation tags.) Then declare the data type, e.g. --input-type as-ont or wgs-ont.")


def classify_input(path: str) -> str:
    """'fastq' | 'ubam' | 'aligned' | 'unreadable' | 'unsupported'. Opens only the header of a BAM/CRAM.

    A uBAM is a .bam whose header declares no reference sequence (@SQ). Any @SQ line means the reads were
    aligned to something (a genome, or the VNTR reference itself)."""
    p = str(path).lower()
    if p.endswith(FASTQ_EXTS):
        return "fastq"
    # Any other extension is SNIFFED, not judged by its name: a uBAM is often saved as `.ubam`, and
    # refusing it for its suffix rejected a valid input (2026-11). pysam recognises BAM/CRAM/SAM by content.
    try:
        import pysam
        with pysam.AlignmentFile(str(path), check_sq=False) as af:
            return "aligned" if af.header.to_dict().get("SQ") else "ubam"
    except Exception:
        return "unreadable" if p.endswith(ALIGNMENT_EXTS) or not os.path.exists(str(path)) else "unsupported"


def input_error(path: str, prog: str = "muc1_analyzer", allow_aligned: bool = False) -> str | None:
    """None when `path` is an accepted input, else the error message to print (then exit 2)."""
    kind = classify_input(path)
    if kind in ("fastq", "ubam"):
        return None
    if kind == "aligned":
        if allow_aligned:
            return None
        return (f"[{prog}] ERROR: {path} is an ALIGNED BAM/CRAM. Since v1.0.1 the input must be FASTQ or an "
                f"unaligned BAM.\n" + EXTRACT_HINT)
    if kind == "unreadable":
        return (f"[{prog}] ERROR: cannot read the header of {path} (missing file, or a CRAM that needs its "
                f"reference). " + EXTRACT_HINT)
    return (f"[{prog}] ERROR: unsupported input {path}: expected FASTQ (.fastq/.fq, optionally .gz) or an "
            f"unaligned BAM.\n"
            + EXTRACT_HINT)
