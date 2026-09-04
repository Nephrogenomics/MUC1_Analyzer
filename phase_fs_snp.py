#!/usr/bin/env python3
"""phase_fs_snp — DIRECT single-molecule phasing of rs4072037 (the splice SNP) with the VNTR allele.

The severity axis of the two-axis model rests on the MUTANT allele's rs4072037 status (C = MUC1-TR, VNTR
retained → severe ; T = MUC1-Y, VNTR spliced out → protected). Until now that phase was INFERRED from a cis
rule (rs4072037-T is cis with the short VNTR ~94 %). Long-read LR-PCR lets us OBSERVE it: rs4072037
(chr1:155,192,276, ~370 bp 5′ of the tandem) sits inside the amplicon, and full-length reads span it AND the
whole VNTR — so each read carries the SNP base AND its allele's copy number. Grouping reads by SNP base gives
each allele's VNTR length directly (chimeras are shorter → the true length is the tallest peak per base).

Then: the mutant allele = the length that carries the frameshift (from the clinical call / pcr_dupc). Its
observed rs4072037 base = the OBSERVED splice status — replacing the inference for the severity axis, and giving
a per-read phasing figure ("why long-read"). Verified on one carrier/P3 (2026-07-17): C↔78 (n=96), T↔44 (n=953);
mutant=44 → T → protected → matches the Normal phenotype and confirms the cis rule.

    python3 phase_fs_snp.py --bam pcr_screen/sample.chr1.bam --mut-len 44 --name SAMPLE
"""
from __future__ import annotations
import argparse
import collections

import pysam

from vntr_raw_length import read_length, GAP

SNP_POS = 155_192_276          # rs4072037, GRCh38 (1-based). Ref C / alt T on this (minus) strand.
SNP_REF, SNP_ALT = "C", "T"    # C = MUC1-TR (severe) · T = MUC1-Y (protected)


def collect_reads(bam, chrom="chr1", snp_pos=SNP_POS, offset=4, maxmm=9, ref=None):
    """[(snp_base, copies), ...] over reads that cover the SNP AND span the whole VNTR. Needs pysam."""
    mode = "rc" if str(bam).endswith(".cram") else "rb"
    kw = {"reference_filename": ref} if (mode == "rc" and ref) else {}
    pairs = []
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        for col in af.pileup(chrom, snp_pos - 1, snp_pos, truncate=True, min_base_quality=0):
            for pr in col.pileups:
                if pr.is_del or pr.query_position is None:
                    continue
                r = pr.alignment
                if r.is_secondary or r.is_supplementary or r.query_sequence is None:
                    continue
                base = r.query_sequence[pr.query_position].upper()
                if base not in (SNP_REF, SNP_ALT):
                    continue
                g = read_length(r.query_sequence.upper(), maxmm)   # only reads that also span the VNTR
                if g is None:
                    continue
                pairs.append((base, round((g - GAP) / 60) + offset))
    return pairs


def _allele_length(cs, sep=12, frac=0.05, min_reads=6, smooth=2):
    """PURE: one SNP base comes from ONE allele. Its real length = the HIGHEST-copy SUBSTANTIAL peak
    (≥ max(min_reads, frac·n) reads); everything shorter is chimeric truncation, and a tiny higher bump is
    noise/SNP-miscall contamination. This threads the needle that neither 'mode' (fails on a depleted long
    allele: P4/P11 C → 27 instead of 77) nor 'max peak' (fails on the short allele: a small 57 bump beats 44)
    can. None if no substantial peak."""
    if not cs:
        return None
    h = collections.Counter(cs)
    lo, hi, total = min(cs), max(cs), len(cs)
    floor = max(min_reads, round(frac * total))

    def sc(c):
        return sum(h.get(c + d, 0) for d in range(-smooth, smooth + 1))

    w = max(2, sep // 2)
    peaks = [c for c in range(lo, hi + 1)
             if sc(c) >= floor and sc(c) >= max(sc(c + d) for d in range(-w, w + 1))]
    return max(peaks) if peaks else None


def phase_snp_to_length(pairs, sep=12, min_reads=8, frac=0.05):
    """PURE: [(snp_base, copies)] -> per-base allele length via `_allele_length`. Returns {base:{length,n}},
    plus `resolved` (both bases present → distinct lengths ≥ sep apart)."""
    by = collections.defaultdict(list)
    for b, c in pairs:
        if b and c is not None:
            by[b].append(c)
    phase = {}
    for b, cs in by.items():
        if len(cs) < min_reads:
            continue
        phase[b] = {"length": _allele_length(cs, sep=sep, frac=frac), "n": len(cs)}
    lengths = [v["length"] for v in phase.values() if v["length"] is not None]
    resolved = (len(phase) == 2 and len(lengths) == 2 and abs(lengths[0] - lengths[1]) >= sep)
    return {"phase": phase, "resolved": resolved}


def mutant_snp(phase_result, mut_len, tol=6):
    """PURE: given the phasing and the MUTANT allele length, return the mutant allele's observed SNP base
    (= splice status). Picks the base whose phased length is closest to `mut_len` (within tol). None if ambiguous."""
    cand = [(b, v["length"]) for b, v in phase_result["phase"].items() if v["length"] is not None]
    hits = [(abs(l - mut_len), b) for b, l in cand if abs(l - mut_len) <= tol]
    if not hits:
        return None
    hits.sort()
    if len(hits) >= 2 and hits[0][0] == hits[1][0]:
        return None                      # equidistant → can't assign (length-homozygous)
    return hits[0][1]


_SPLICE = {"C": ("MUC1-TR", "VNTR retained → severe"), "T": ("MUC1-Y", "VNTR spliced out → protected")}


def render(pairs, out_prefix, name, mut_len=None):
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    res = phase_snp_to_length(pairs)
    fig, ax = plt.subplots(figsize=(9, 4.6))
    col = {"C": "#b71c1c", "T": "#1565c0"}
    for base in ("C", "T"):
        cs = [c for b, c in pairs if b == base]
        if not cs:
            continue
        h = collections.Counter(c for c in cs if 0 < c < 130)
        xs = sorted(h)
        ax.bar([x for x in xs], [h[x] for x in xs], width=0.9, color=col[base], alpha=0.65,
               label=f"rs4072037-{base} ({_SPLICE[base][0]}, n={len(cs)})")
        L = res["phase"].get(base, {}).get("length")
        if L is not None:
            ax.axvline(L, color=col[base], lw=2, ls="--")
            ax.text(L, ax.get_ylim()[1] * 0.9, f" {L}", color=col[base], fontweight="bold")
    if mut_len is not None:
        mb = mutant_snp(res, mut_len)
        tag = f"mutant allele = {mut_len} copies → rs4072037-{mb} ({_SPLICE.get(mb, ('?', '?'))[1]})" if mb else \
              f"mutant allele = {mut_len} copies → unresolved"
        ax.set_title(f"{name} — direct single-molecule phasing of rs4072037\n{tag}",
                     fontsize=11, fontweight="bold")
    else:
        ax.set_title(f"{name} — direct single-molecule phasing of rs4072037", fontsize=11, fontweight="bold")
    ax.set_xlabel("VNTR copies (per read that also spans rs4072037)")
    ax.set_ylabel("reads")
    ax.legend(fontsize=8.5)
    ax.grid(True, axis="y", alpha=0.2)
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{out_prefix}.{ext}", dpi=200, bbox_inches="tight")
    print(f"wrote {out_prefix}.png / .pdf")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="phase_fs_snp")
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", default=None)
    ap.add_argument("--name", default=None)
    ap.add_argument("--chrom", default="chr1")
    ap.add_argument("--snp-pos", type=int, default=SNP_POS)
    ap.add_argument("--offset", type=int, default=4, help="LR-PCR copy convention (+4)")
    ap.add_argument("--maxmm", type=int, default=9)
    ap.add_argument("--mut-len", type=int, default=None, help="mutant allele copies → report its observed SNP base")
    ap.add_argument("--fig", default=None)
    a = ap.parse_args(argv)
    name = a.name or a.bam.split("/")[-1]
    pairs = collect_reads(a.bam, a.chrom, a.snp_pos, a.offset, a.maxmm, a.ref)
    res = phase_snp_to_length(pairs)
    print(f"=== {name} — direct rs4072037 phasing ({len(pairs)} reads span SNP+VNTR) ===")
    for base in ("C", "T"):
        v = res["phase"].get(base)
        if v:
            print(f"  rs4072037-{base} ({_SPLICE[base][0]:8}): allele = {v['length']} copies  (n={v['n']})")
    print(f"  resolved (two distinct alleles): {res['resolved']}")
    if a.mut_len is not None:
        mb = mutant_snp(res, a.mut_len)
        if mb:
            print(f"  -> mutant allele ({a.mut_len} copies) carries rs4072037-{mb} = {_SPLICE[mb][0]} "
                  f"({_SPLICE[mb][1]})  [OBSERVED, not inferred]")
        else:
            print(f"  -> mutant allele ({a.mut_len} copies): unresolved (length-homozygous or ambiguous)")
    if a.fig:
        render(pairs, a.fig, name, a.mut_len)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
