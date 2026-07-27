#!/usr/bin/env python3
"""
BAM Cleaner — filter a BAM's reads to keep only those that actually cover a
target region (e.g. the MUC1 VNTR), with a reasonable soft-clip and a minimal
overlap. Serves as an upstream step for MUC1_Analyzer.

(Reconstructed version: the original file had been corrupted by stray line
breaks that broke the Python tokens.)

Usage:
    python3 Bam_cleaner.py input.bam -o output.bam \
        -region NC_060925.1:154324904-154335104 \
        --softclip-threshold 10000 --min-span 5000
"""
import argparse
import sys
from pathlib import Path

import pysam


def parse_region(region_str):
    """'chr:start-end' -> (chrom, start, end)."""
    try:
        chrom, coords = region_str.split(":")
        start, end = coords.split("-")
        return chrom, int(start), int(end)
    except Exception as e:
        raise ValueError(
            f"Invalid region format: {region_str}. "
            f"Use 'chr:start-end'."
        ) from e


def calculate_overlap(read_start, read_end, region_start, region_end):
    """Read/region overlap length (0 if disjoint)."""
    if read_start is None or read_end is None:
        return 0
    return max(0, min(read_end, region_end) - max(read_start, region_start))


def get_softclip_length(read):
    """Total soft-clip length (start + end; CIGAR operation 4)."""
    total = 0
    cig = read.cigartuples
    if cig:
        if cig[0][0] == 4:
            total += cig[0][1]
        if cig[-1][0] == 4:
            total += cig[-1][1]
    return total


def classify_read(read, chrom, region_start, region_end, softclip_threshold,
                  min_span, min_mq, keep_supplementary):
    """Return None if the read is kept, otherwise a rejection reason (str)."""
    if read.is_unmapped:
        return "unmapped"
    if read.is_secondary:
        return "secondary"
    if read.is_supplementary and not keep_supplementary:
        return "supplementary"
    if read.is_duplicate:
        return "duplicate"
    if read.mapping_quality < min_mq:
        return "low_mapq"
    if read.reference_name != chrom:
        return "wrong_chrom"
    if calculate_overlap(read.reference_start, read.reference_end,
                         region_start, region_end) < min_span:
        return "low_span"
    if get_softclip_length(read) > softclip_threshold:
        return "high_softclip"
    return None


def should_keep_read(read, chrom, region_start, region_end,
                     softclip_threshold, min_span, min_mq=0, keep_supplementary=False):
    return classify_read(read, chrom, region_start, region_end, softclip_threshold,
                         min_span, min_mq, keep_supplementary) is None


def _open_input(path, reference):
    """Open a BAM or CRAM (CRAM -> reference_filename required)."""
    mode = "rc" if str(path).endswith(".cram") else "rb"
    kw = {"reference_filename": reference} if (mode == "rc" and reference) else {}
    return pysam.AlignmentFile(path, mode, **kw)


def clean_bam(input_bam, output_bam, region_str, softclip_threshold, min_span,
              min_mq=0, keep_supplementary=False, reference=None):
    from collections import Counter
    chrom, region_start, region_end = parse_region(region_str)
    print(f"File          : {input_bam}")
    print(f"Target region : {chrom}:{region_start}-{region_end}")
    print(f"Thresholds    : soft-clip≤{softclip_threshold}  span≥{min_span}  MAPQ≥{min_mq}"
          f"  supplementary={'kept' if keep_supplementary else 'dropped'}")

    bam_in = _open_input(input_bam, reference)
    # CRAM output if the extension asks for it (else BAM)
    out_mode = "wc" if str(output_bam).endswith(".cram") else "wb"
    out_kw = {"reference_filename": reference} if (out_mode == "wc" and reference) else {}
    bam_out = pysam.AlignmentFile(output_bam, out_mode, template=bam_in, **out_kw)

    total = kept = 0
    dropped = Counter()
    try:
        try:
            iterator = bam_in.fetch(chrom, region_start, region_end)
        except (ValueError, KeyError):
            iterator = bam_in.fetch(until_eof=True)
        for read in iterator:
            total += 1
            reason = classify_read(read, chrom, region_start, region_end,
                                   softclip_threshold, min_span, min_mq, keep_supplementary)
            if reason is None:
                bam_out.write(read)
                kept += 1
            else:
                dropped[reason] += 1
            if total % 10000 == 0:
                print(f"Processed: {total} reads, kept: {kept}", end="\r")
    finally:
        bam_in.close()
        bam_out.close()

    print(f"\nDone. Total reads: {total}; kept: {kept} "
          f"({(kept/total*100 if total else 0):.2f} %)")
    if dropped:
        print("Rejections by reason: " + "  ".join(f"{k}={v}" for k, v in dropped.most_common()))
    print(f"Output: {output_bam}")
    # output index (handy for MUC1_Analyzer)
    try:
        pysam.index(str(output_bam))
        print(f"Index : {output_bam}.bai")
    except Exception as e:  # noqa
        print(f"(index not created: {e})", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(
        description="Clean a BAM to keep only the reads covering a target "
                    "region (controlled soft-clip + overlap).")
    ap.add_argument("input_bam", help="input BAM or CRAM")
    ap.add_argument("-o", "--output", required=True, help="output BAM/CRAM (.cram -> CRAM)")
    ap.add_argument("-region", "--region", required=True,
                    help="target region 'chr:start-end'")
    ap.add_argument("-r", "--reference", default=None,
                    help="reference FASTA (required for a CRAM)")
    ap.add_argument("--softclip-threshold", type=int, default=10000,
                    help="max allowed soft-clip (default: 10000)")
    ap.add_argument("--min-span", type=int, default=5000,
                    help="minimal overlap with the region (default: 5000)")
    ap.add_argument("--min-mq", type=int, default=0,
                    help="minimal mapping quality (default: 0 = no filter)")
    ap.add_argument("--keep-supplementary", action="store_true",
                    help="keep supplementary alignments (default: dropped)")
    args = ap.parse_args()

    if not Path(args.input_bam).exists():
        sys.exit(f"Error: {args.input_bam} not found")
    if str(args.input_bam).endswith(".cram") and not args.reference:
        sys.exit("Error: a CRAM requires --reference <FASTA>")
    clean_bam(args.input_bam, args.output, args.region,
              args.softclip_threshold, args.min_span, min_mq=args.min_mq,
              keep_supplementary=args.keep_supplementary, reference=args.reference)


if __name__ == "__main__":
    main()
