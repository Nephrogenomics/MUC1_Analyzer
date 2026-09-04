#!/usr/bin/env python3
"""variant_profile — the variant fraction REPEAT BY REPEAT, so a borderline call can be adjudicated.

WHY. `pcr_variant` returns one number: the best index and its fraction. On 2026-08-07 a delCC carrier (a delCC
carrier) came back `frac=0.091` against a `floor=0.10` — a miss by nine thousandths — and that single
number cannot say whether it is a real variant just under the line or ONT homopolymer background that
happened to peak somewhere. Those two have completely different remedies and the same summary statistic.

The discriminator is SHAPE, not height. A real frameshift sits at ONE array unit: a focal spike above a
flat baseline. Homopolymer error is spread: every unit carries a few percent and the "best" index is just
the luckiest draw. `scan_segments` was measured on 2026-08-02 as non-discriminating in POOLED counts for
exactly this reason — pooling destroys the only feature that separates them.

⚠ THIS REPORTS, IT DOES NOT DECIDE. `focality` is a descriptive ratio (peak over median), not a calibrated
test: it has no null, no specificity, and must never gate a call. It exists so a human can look at the
shape — and so the visual review bundle has something to draw.

⚠ The index is the ordinal of the unit ALONG THE SCAFFOLD, which is not the clinical repeat number unless
the scaffold length equals the true allele length. Pass `--expect` to mark where the clinic says the
variant is, and read the distance as a sanity check, never as agreement.

    python3 -m muc1_analyzer.variant_profile -b pcr_lot180_raw_2s/SAMPLE.vntr.bam \\
        --contig MUC1_VNTR_73repeats --variant delCC --expect 50
"""
from __future__ import annotations
import argparse
import collections
import statistics
import sys


def index_profile(reads, target: str, *, min_tot: int = 20) -> list:
    """[{index, tot, n, frac, other}] per array unit, over units covered >= min_tot. Pure.

    reads = [(n_units, strand, [seg_type per unit]), …] — the shape `pcr_variant._scan_segments` returns."""
    idx = collections.defaultdict(collections.Counter)
    for _n, _strand, classes in reads:
        for i, t in enumerate(classes):
            idx[i][t] += 1
    out = []
    for i in sorted(idx):
        tot = sum(idx[i].values())
        if tot < min_tot:
            continue
        out.append({"index": i, "tot": tot, "n": idx[i][target],
                    "frac": idx[i][target] / tot, "other": idx[i]["other"] / tot})
    return out


#: Absolute level a unit must reach to count as part of a block. Set from the NEGATIVE controls, not from
#: a carrier: 8 lot negatives peak at 0.016-0.052 delCC with a MEDIAN of 0.000, so the background is not
#: "a few percent everywhere" — it is zero with isolated one-unit accidents. 0.035 sits above the ordinary
#: accident and below the 0.039-0.091 of that carrier's block. ⚠ PROVISIONAL: 8 negatives is a first look, not a
#: licence. Price it on `lists/neg_validation.txt` before it gates anything.
BLOCK_LEVEL = 0.035


def longest_block(profile: list, *, level: float = BLOCK_LEVEL) -> dict:
    """Longest run of CONSECUTIVE array units at or above `level`. Pure.

    This is the discriminator `focality` missed. Height alone cannot separate a real variant from noise
    here, because pooling two alleles that differ by a couple of copies puts the SAME physical unit at
    different indices in different reads — it smears one spike into a short block. Noise does the
    opposite: it hits isolated units. So the question is not "how high" but "how many in a row".

    ⚠ Contiguity is measured over the units PRESENT in `profile`. Units dropped for low coverage break a
    run rather than silently bridging it — a gap you cannot see is not a gap you may assume is filled."""
    best = {"len": 0, "start": None, "end": None, "max": None}
    run, start, peak, prev_idx = 0, None, 0.0, None
    for d in profile:
        contiguous = prev_idx is None or d["index"] == prev_idx + 1
        if d["frac"] >= level and (contiguous or run == 0):
            if run == 0:
                start, peak = d["index"], 0.0
            run += 1
            peak = max(peak, d["frac"])
        elif d["frac"] >= level:                       # above level but the index jumped — a new run
            run, start, peak = 1, d["index"], d["frac"]
        else:
            run = 0
        if run > best["len"]:
            best = {"len": run, "start": start, "end": d["index"], "max": round(peak, 4)}
        prev_idx = d["index"]
    return best


def focality(profile: list) -> dict:
    """Peak over baseline: {peak_index, peak_frac, median_frac, ratio, n_above_half_peak}. Pure.

    `ratio` is descriptive. A focal variant gives a large ratio with n_above_half_peak == 1; a flat
    homopolymer background gives a ratio near 1 with many units near the peak. Neither is a p-value."""
    if not profile:
        return {"peak_index": None, "peak_frac": None, "median_frac": None,
                "ratio": None, "n_above_half_peak": 0}
    peak = max(profile, key=lambda d: d["frac"])
    med = statistics.median([d["frac"] for d in profile])
    half = peak["frac"] / 2 if peak["frac"] else 0
    return {"peak_index": peak["index"], "peak_frac": round(peak["frac"], 4),
            "median_frac": round(med, 4),
            "ratio": round(peak["frac"] / med, 2) if med else None,
            "n_above_half_peak": sum(1 for d in profile if d["frac"] >= half > 0)}


def render(profile: list, foc: dict, *, target: str, expect=None, width: int = 40) -> str:
    """The profile as text bars — the shape is the point, so it must be visible without a plot. Pure."""
    if not profile:
        return "no array unit reached the coverage floor — nothing to profile"
    top = max(d["frac"] for d in profile) or 1.0
    L = [f"{'idx':>4} {'tot':>5} {'n':>4} {target+' frac':>11}  {'other':>6}",
         "-" * (28 + width)]
    for d in profile:
        bar = "#" * int(round(d["frac"] / top * width))
        mark = "  <= clinic" if (expect is not None and d["index"] == expect) else ""
        L.append(f"{d['index']:>4} {d['tot']:>5} {d['n']:>4} {d['frac']:>11.3f}  "
                 f"{d['other']:>6.2f} {bar}{mark}")
    blk = longest_block(profile)
    L += ["", f"peak at index {foc['peak_index']} = {foc['peak_frac']}, median {foc['median_frac']}, "
              f"peak/median {foc['ratio']}, units within half the peak: {foc['n_above_half_peak']}",
          f"longest block >= {BLOCK_LEVEL}: {blk['len']} unit(s)"
          + (f" at {blk['start']}-{blk['end']}, max {blk['max']}" if blk["len"] else "")]
    # The verdict reads the BLOCK, not the height. Height alone mislabelled that carrier as background on
    # 2026-08-07: its 0.091 looked like one lucky unit until the negatives showed their peaks are
    # single-unit accidents at 0.016-0.052 and the carrier's is six units wide.
    if blk["len"] >= 3:
        L.append(f"SHAPE: BLOCK — {blk['len']} consecutive units above {BLOCK_LEVEL}. Noise here hits "
                 "ISOLATED units; a variant smeared by allele pooling occupies adjacent ones. This is "
                 "the carrier-like shape.")
    elif blk["len"]:
        L.append(f"SHAPE: ISOLATED — {blk['len']} unit(s) above {BLOCK_LEVEL}, no run. That is what the "
                 "lot negatives look like (single accidents, peaks 0.016-0.052).")
    else:
        L.append(f"SHAPE: FLAT — nothing reaches {BLOCK_LEVEL}. No candidate to adjudicate.")
    L.append("⚠ descriptive only — no null, no specificity, gates nothing; BLOCK_LEVEL is priced on 8 "
             "negatives and is provisional")
    if expect is not None and foc["peak_index"] is not None:
        L.append(f"peak is {abs(foc['peak_index'] - expect)} unit(s) from the clinical repeat {expect} "
                 f"(scaffold indices are not clinical repeat numbers — a sanity check, not agreement)")
    return "\n".join(L)


BATCH_COLS = ("sample", "status", "contig", "n_reads", "n_units", "peak_index", "peak_frac",
              "median_frac", "block_len", "block_start", "block_end", "block_max")


def _batch(a) -> int:
    """One row per sample, so the block statistic gets PRICED rather than admired on one carrier."""
    import os
    from .miss_report import arbiter_lengths_of, find_bam
    from .miss_report import read_truth
    from .pcr_truth import is_carrier, is_negative

    rows = []
    for t in read_truth(a.truth):
        st = "pos" if is_carrier(t.get("status")) else ("neg" if is_negative(t.get("status")) else None)
        if st is None:                                 # `unknown` is neither arm — never a silent negative
            continue
        bam = find_bam(t["sample"], a.outdir or [])
        if not bam:
            continue
        lens = arbiter_lengths_of(os.path.dirname(bam), t["sample"])
        try:
            for contig, prof, reads in _profiles(bam, a, lens):
                foc, blk = focality(prof), longest_block(prof)
                rows.append({"sample": t["sample"], "status": st, "contig": contig, "n_reads": reads,
                             "n_units": len(prof), "peak_index": foc["peak_index"],
                             "peak_frac": foc["peak_frac"], "median_frac": foc["median_frac"],
                             "block_len": blk["len"], "block_start": blk["start"],
                             "block_end": blk["end"], "block_max": blk["max"]})
        except Exception as e:                         # one unreadable bam must not lose the cohort
            print(f"[warn] {t['sample']}: {e}", file=sys.stderr)
    if a.out:
        with open(a.out, "w") as fh:
            fh.write("\t".join(BATCH_COLS) + "\n")
            for r in rows:
                fh.write("\t".join("" if r.get(c) is None else str(r.get(c, "")) for c in BATCH_COLS) + "\n")
        print(f"[profile] {len(rows)} rows → {a.out}", file=sys.stderr)
    pos = [r for r in rows if r["status"] == "pos"]
    neg = [r for r in rows if r["status"] == "neg"]
    print(f"{'block_len':>10}{'carriers':>10}{'negatives':>11}")
    print("  a block of N consecutive units above "
          f"{BLOCK_LEVEL}. Read DOWN the negative column: the length at which it empties is the")
    print("  only thing that licenses a block rule, and this is the lot's own negatives, not the "
          "independent set.")
    for n in range(1, 8):
        print(f"{'>= ' + str(n):>10}{sum(1 for r in pos if r['block_len'] >= n):>10}"
              f"{sum(1 for r in neg if r['block_len'] >= n):>11}")
    return 0


def _profiles(bam: str, a, lens):
    """[(contig, profile, n_reads)] over the chosen scaffolds. Impure."""
    import pysam
    from .allele_scaffold import allele_contigs
    from .pcr_variant import _scan_segments
    contigs = [a.contig] if a.contig else allele_contigs(bam, lengths=lens)
    out = []
    with pysam.AlignmentFile(bam, "rb") as fh:
        lengths = {c: fh.get_reference_length(c) for c in contigs}
    for c in contigs:
        reads = _scan_segments(bam, c, 1, int(lengths[c]), None)
        out.append((c, index_profile(reads, a.variant, min_tot=a.min_tot), len(reads)))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="muc1_analyzer.variant_profile",
                                 description="Per-repeat variant fraction along one allele scaffold")
    ap.add_argument("-b", "--bam", default=None, help="the `.vntr.bam` a run already wrote")
    ap.add_argument("--contig", default=None, help="scaffold contig (default: the allele scaffolds)")
    ap.add_argument("--variant", default="delCC")
    ap.add_argument("--expect", type=int, default=None, help="clinical repeat number, marked in the output")
    ap.add_argument("--min-tot", type=int, default=20)
    ap.add_argument("--lengths", default=None, metavar="43,96", help="arbiter lengths, to pick scaffolds")
    # BATCH: one row per sample so the block statistic can be PRICED on the negatives instead of admired
    # on one carrier. `-b` is not required in this mode.
    ap.add_argument("--truth", default=None, help="batch mode: profile every sample in this truth table")
    ap.add_argument("--outdir", nargs="+", default=None, help="batch mode: dirs holding <sample>.vntr.bam")
    ap.add_argument("--out", default=None, help="batch mode: TSV (patient data — keep on /scratch)")
    a = ap.parse_args(argv)
    if a.truth:
        return _batch(a)
    if not a.bam:
        ap.error("-b/--bam is required (or use --truth/--outdir for batch mode)")

    lens = [int(x) for x in a.lengths.split(",")] if a.lengths else None
    got = _profiles(a.bam, a, lens)                    # ONE implementation, shared with batch mode
    if not got:
        print("no scaffold contig — pass --contig", file=sys.stderr)
        return 1
    for contig, prof, n_reads in got:
        print(f"\n═══ {contig}  ({n_reads} reads, {a.variant}) ═══")
        print(render(prof, focality(prof), target=a.variant, expect=a.expect))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
