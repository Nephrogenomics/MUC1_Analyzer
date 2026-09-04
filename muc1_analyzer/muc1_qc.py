"""Per-subject QC — MUC1 coverage OUTSIDE the VNTR + protocol (WGS/AS) + chemistry, to STRATIFY PoN/validation.

The VNTR is unreliable (collapsed, GC-rich) → we measure depth over **unique non-VNTR** MUC1 windows
(5′UTR/exon1, CDS exon6-7, intron6 flank) = an honest proxy of the subject's coverage. A control called
"negative" at LOW coverage is **non-informative** (not enough reads for a false positive) → to be
distinguished from well-covered ones. An **off-target** window (gene desert) gives the on/off ratio:
**adaptive sampling** (bed_CKD) enriches MUC1 → ratio ≫ 1 ; **WGS** → ratio ≈ 1.
The **chemistry** (R9/R10) is read from the file name (not the BAM). GRCh38 "chr" coordinates.
"""
from __future__ import annotations
import re

import pysam

# UNIQUE MUC1 windows outside the VNTR (chr1, GRCh38) — cf. ASM Note / config
WINDOWS = {
    "utr5_exon1": ("chr1", 155192300, 155192900),
    "cds_exon6_7": ("chr1", 155186500, 155187500),
    "intron6_flank": ("chr1", 155188000, 155188800),
}
OFFTARGET = ("chr1", 50_000_000, 50_010_000)      # gene desert ≈ off-panel proxy (AS detection)


def chemistry_from_name(name: str) -> str:
    """Chemistry inferred from the file name: 'R10' / 'R9' / '?'."""
    if re.search(r"R10|LSK114", name, re.I):
        return "R10"
    if re.search(r"R9|LSK110|LSK109", name, re.I):
        return "R9"
    return "?"


def _mean_depth(af, chrom, start, end) -> float:
    cov = af.count_coverage(chrom, start, end)     # 4 arrays (A,C,G,T), base-quality-filtered depth
    n = end - start
    return sum(sum(c) for c in cov) / n if n else 0.0


def qc(bam: str, *, genome_ref: str = None, name: str = None) -> dict:
    """Mean depth per non-VNTR window + off-target + ratio + flags (protocol, coverage)."""
    mode = "rc" if bam.endswith(".cram") else "rb"
    kw = {"reference_filename": genome_ref} if (mode == "rc" and genome_ref) else {}
    depths = {}
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        for label, (c, s, e) in WINDOWS.items():
            depths[label] = round(_mean_depth(af, c, s, e), 1)
        off = round(_mean_depth(af, *OFFTARGET), 2)
    muc1_mean = round(sum(depths.values()) / len(depths), 1)
    ratio = round(muc1_mean / off, 1) if off > 0 else None
    return {"name": name or bam, "chemistry": chemistry_from_name(name or bam),
            "depth_by_window": depths, "muc1_nonvntr_depth": muc1_mean,
            "offtarget_depth": off, "on_off_ratio": ratio,
            "protocol_guess": ("AS" if (ratio and ratio >= 5) else "WGS"),
            "cov_flag": ("LOW" if muc1_mean < 10 else "OK")}


if __name__ == "__main__":
    import argparse
    import json
    ap = argparse.ArgumentParser(prog="muc1_analyzer.muc1_qc",
                                 description="MUC1 non-VNTR coverage QC + protocol/chemistry (stratification)")
    ap.add_argument("-b", "--bam", required=True)
    ap.add_argument("--genome-ref", default=None)
    ap.add_argument("--name", default=None, help="file name (for the chemistry; default = path)")
    ap.add_argument("--tsv", action="store_true", help="one TSV line (to aggregate a cohort)")
    args = ap.parse_args()
    r = qc(args.bam, genome_ref=args.genome_ref, name=args.name or args.bam)
    if args.tsv:
        print("\t".join(str(x) for x in [
            r["name"].split("/")[-1], r["chemistry"], r["muc1_nonvntr_depth"],
            r["offtarget_depth"], r["on_off_ratio"], r["protocol_guess"], r["cov_flag"]]))
    else:
        import json as _j
        print(_j.dumps(r, indent=2, ensure_ascii=False))
