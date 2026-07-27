"""Build the multi-contig MUC1 VNTR reference (MUC1_Analyzer's `-r`).

Starts from the real MUC1 locus sequence in a genome (T2T-CHM13 recommended),
detects the VNTR array (tandem repeats of the 60 bp unit), then emits one contig
per possible length: `MUC1_VNTR_Nrepeats` = 5′ flank + N×unit + 3′ flank, N=1..150.

Orientation: MUC1 is on the − strand. The `MOTIF_A` unit is in gene/coding orientation
(the one MUC1_Analyzer expects). Use --revcomp when the region extracted from the genome is
in forward orientation (T2T/GRCh38 case) to bring it back into gene orientation.

Usage (on a host where the T2T genome is available):
    python -m muc1_analyzer.build_vntr_ref \
        --genome /path/to/chm13v2.0.fa \
        --region NC_060925.1:154324904-154335104 --revcomp \
        -o /path/to/vntr_ref.fa
"""
from __future__ import annotations
import argparse
import sys

import pysam

# Canonical VNTR unit (60 bp, gene/coding orientation) = MUC1_Analyzer's 'A' motif.
MOTIF_A = "GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCGCCCGCA"

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def find_array(seq: str, unit: str, max_mismatch: int) -> tuple:
    """Locate the longest tandem array of `unit` (step = len(unit), Hamming tolerance).

    Returns (n_repeats, array_start, array_end). (0,0,0) if nothing found.
    """
    L = len(unit)
    n = len(seq)
    best = (0, 0, 0)
    s = 0
    while s <= n - L:
        if hamming(seq[s:s + L], unit) <= max_mismatch:
            k, reps = s, 0
            while k + L <= n and hamming(seq[k:k + L], unit) <= max_mismatch:
                reps += 1
                k += L
            if reps > best[0]:
                best = (reps, s, k)
            s = k                      # skip past the array just found
        else:
            s += 1
    return best


def build(genome: str, chrom: str, start: int, end: int, *, unit: str = MOTIF_A,
          do_revcomp: bool = False, nmin: int = 1, nmax: int = 150,
          max_mismatch: int = 6, array_start: int = None, array_end: int = None,
          flank_len: int = None, out: str = "vntr_ref.fa") -> dict:
    fa = pysam.FastaFile(genome)
    seq = fa.fetch(chrom, start, end).upper()
    fa.close()
    if do_revcomp:
        seq = revcomp(seq)

    if array_start is not None and array_end is not None:
        a0, a1 = array_start, array_end
        reps = (a1 - a0) // len(unit)
    else:
        reps, a0, a1 = find_array(seq, unit, max_mismatch)
        if reps == 0:
            raise RuntimeError(
                "VNTR array not found: check --region, --revcomp, or --max-mismatch "
                "(the unit must appear in tandem in the extracted sequence).")

    flank5, flank3 = seq[:a0], seq[a1:]
    # short flanks: large minimap2 speed gain (the identical flanks × N contigs
    # are the main source of chaining blow-up). ~200 bp are enough to anchor.
    if flank_len:
        flank5 = flank5[-flank_len:]
        flank3 = flank3[:flank_len]
    with open(out, "w") as fh:
        for N in range(nmin, nmax + 1):
            contig = flank5 + unit * N + flank3
            fh.write(f">MUC1_VNTR_{N}repeats\n")
            for i in range(0, len(contig), 70):
                fh.write(contig[i:i + 70] + "\n")
    pysam.faidx(out)

    info = {"detected_repeats": reps, "array_start": a0, "array_end": a1,
            "flank5_len": len(flank5), "flank3_len": len(flank3),
            "extracted_len": len(seq), "contigs": nmax - nmin + 1, "out": out}
    return info


def main(argv=None):
    ap = argparse.ArgumentParser(prog="muc1_analyzer.build_vntr_ref",
                                 description="Build the multi-contig MUC1 VNTR reference")
    ap.add_argument("--genome", required=True, help="genome FASTA (T2T-CHM13 recommended)")
    ap.add_argument("--region", required=True, help="chrom:start-end (1-based) of the MUC1 locus")
    ap.add_argument("--revcomp", action="store_true", help="reverse-complement (gene strand)")
    ap.add_argument("--unit", default=MOTIF_A, help="60 bp VNTR unit (default motif A)")
    ap.add_argument("--nmin", type=int, default=1)
    ap.add_argument("--nmax", type=int, default=150)
    ap.add_argument("--max-mismatch", type=int, default=6, help="Hamming tolerance for array detection")
    ap.add_argument("--array-start", type=int, default=None, help="override: array start (0-based, post-revcomp)")
    ap.add_argument("--array-end", type=int, default=None)
    ap.add_argument("--flank-len", type=int, default=None,
                    help="trim the flanks to N bp on each side (minimap2 speed; e.g. 200)")
    ap.add_argument("-o", "--output", default="vntr_ref.fa")
    args = ap.parse_args(argv)

    chrom, coords = args.region.split(":")
    s, e = coords.replace(",", "").split("-")
    start, end = int(s) - 1, int(e)                # 1-based -> 0-based half-open


    info = build(args.genome, chrom, start, end, unit=args.unit, do_revcomp=args.revcomp,
                 nmin=args.nmin, nmax=args.nmax, max_mismatch=args.max_mismatch,
                 array_start=args.array_start, array_end=args.array_end,
                 flank_len=args.flank_len, out=args.output)
    print(f"[OK] {info['out']}  ({info['contigs']} contigs)", file=sys.stderr)
    print(f"  VNTR array detected: {info['detected_repeats']} repeats "
          f"(pos {info['array_start']}-{info['array_end']} in the extracted sequence)", file=sys.stderr)
    print(f"  5′ flank={info['flank5_len']} bp  3′ flank={info['flank3_len']} bp  "
          f"(total extract {info['extracted_len']} bp)", file=sys.stderr)
    print("  ⚠ check: 5′ flank consistent with the GFF3 (VNTR ~4.3 kb from the start in the 7-repeats ref).",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
