#!/usr/bin/env python3
"""allele_scaffold — pick ONE scaffold contig per allele on the multi-contig VNTR reference.

A `.vntr.bam` is aligned to all 150 `MUC1_VNTR_Nrepeats` contigs, so a near-identical tandem smears its
reads across neighbouring lengths. Anything that must read a per-allele signal (the run-length dupC test,
a per-allele consensus) first has to decide which contig IS each allele.

Taking the single modal contig is wrong on a heterozygote: LR-PCR amplifies the short allele better, so the
mode is the short one every time — measured on the cohort, scaffold 24 repeats for true alleles 59/79, 11
for 12/43, 13 for 11/75. The long allele, which is where the dupC usually sits, was never examined.

Two strategies, in order of preference:
  · `contigs_for_lengths` — the arbiter has ALREADY called the allele lengths, on the physical scale these
    contigs are named in. Use them. This is the only thing that recovers an allele the PCR suppressed
    (a real long allele sat at ~1.7 % of the mapped reads on a cohort carrier).
  · `pick_allele_contigs` — fallback when no length call exists: read modes, with the distance-dependent
    scatter rule (a neighbouring second peak must reach 0.5 of the dominant, a distant one 0.25, with an
    absolute read floor as the only guard against tiny distant noise).
"""
from __future__ import annotations

import re

CONTIG_LEN = re.compile(r"MUC1_VNTR_(\d+)repeats")


def contig_counts(bam: str) -> dict:
    """Mapped reads per contig, from the index. Impure."""
    import pysam
    with pysam.AlignmentFile(bam) as af:
        return {x.contig: x.mapped for x in af.get_index_statistics() if x.mapped > 0}


def modal_contig(counts: dict):
    """The single busiest contig. Pure. Deterministic on a tie."""
    if not counts:
        return None
    return max(counts.items(), key=lambda kv: (kv[1], kv[0]))[0]


def contigs_for_lengths(counts: dict, lengths) -> list:
    """Scaffold contigs for KNOWN allele lengths — nearest read-carrying contig per length. Pure."""
    avail = {}
    for name, n in counts.items():
        m = CONTIG_LEN.search(name)
        if m and n > 0:
            avail[int(m.group(1))] = name
    if not avail:
        return []
    out = []
    for L in lengths or []:
        nearest = min(avail, key=lambda k: (abs(k - int(L)), k))
        if avail[nearest] not in out:
            out.append(avail[nearest])
    return sorted(out, key=lambda nm: int(CONTIG_LEN.search(nm).group(1)))


def pick_allele_contigs(counts: dict, *, window: int = 2, min_sep: int = 8,
                        ratio_near: float = 0.5, ratio_far: float = 0.25,
                        min_reads: int = 5, max_alleles: int = 2) -> list:
    """Up to `max_alleles` scaffold contigs, one per read MODE. Pure. Fallback — prefer the arbiter.

    A mode is scored over a +/-`window` window (stutter lands on neighbouring lengths). Whether a second
    mode is an allele or scatter uses the distance-dependent rule: within `min_sep` it must reach
    `ratio_near` of the dominant, beyond it only `ratio_far`; the guard against tiny distant noise is the
    absolute `min_reads` floor, never the ratio.

    ⚠ Not a depletion rescue: an allele the PCR pushed to ~1.7 % of reads fails `ratio_far` too."""
    lens = {}
    for name, n in counts.items():
        m = CONTIG_LEN.search(name)
        if m:
            lens.setdefault(int(m.group(1)), [0, name])[0] += n
    if not lens:
        return [modal_contig(counts)] if counts else []
    picked, masked, dominant = [], set(), None
    for _ in range(max_alleles):
        best, best_key = None, None
        for L in sorted(lens):
            if L in masked:
                continue
            w = sum(v[0] for k, v in lens.items() if abs(k - L) <= window)
            key = (w, lens[L][0])          # ties go to the contig actually carrying the reads
            if best_key is None or key > best_key:
                best, best_key = L, key
        if best is None:
            break
        if picked:
            thr = ratio_near if abs(best - dominant[0]) < min_sep else ratio_far
            if best_key[0] < min_reads or best_key[0] < thr * dominant[1]:
                break
        else:
            dominant = (best, best_key[0])
        picked.append((best, lens[best][1]))
        masked |= {k for k in lens if abs(k - best) <= 2 * window}
    return [name for _, name in sorted(picked)]


def indel_burden(bam: str, contig: str) -> tuple:
    """(total indel bp, reads) charged to the reads aligned on `contig`. Impure.

    The reads themselves say which contig they fit: on a scaffold too short the extra units pile up as
    insertions, on one too long as deletions. The minimum over a window is the length that needs the least
    explaining — an ALIGNMENT-based estimate, independent of the alignment-free flank count."""
    import pysam
    total = n = 0
    with pysam.AlignmentFile(bam) as af:
        for r in af.fetch(contig):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or not r.cigartuples:
                continue
            total += sum(ln for op, ln in r.cigartuples if op in (1, 2))     # I and D
            n += 1
    return total, n


def burden_sweep(bam: str, *, center: int, window: int = 6, min_reads: int = 3) -> dict:
    """Scan +/-`window` contigs around `center` and keep the lowest per-read indel burden. Impure.

    Returns {best, burden_per_read, table, n_reads} or {} when nothing carries enough reads. Per-READ
    burden, not total: a contig with more reads would otherwise always look worse."""
    counts = contig_counts(bam)
    have = {}
    for name, n in counts.items():
        m = CONTIG_LEN.search(name)
        if m:
            have[int(m.group(1))] = name
    table = []
    for L in sorted(k for k in have if abs(k - center) <= window):
        total, n = indel_burden(bam, have[L])
        if n >= min_reads:
            table.append({"length": L, "burden_per_read": round(total / n, 1), "n_reads": n})
    if not table:
        return {}
    best = min(table, key=lambda t: (t["burden_per_read"], abs(t["length"] - center)))
    return {"best": best["length"], "burden_per_read": best["burden_per_read"],
            "n_reads": best["n_reads"], "table": table}


def sweep_lengths(bam: str, centers, *, window: int = 6) -> list:
    """Burden-swept length for each allele centre. Impure, never raises. [] when it cannot measure."""
    out = []
    for c in centers or []:
        try:
            r = burden_sweep(bam, center=int(c), window=window)
        except Exception:
            r = {}
        if r:
            out.append(r["best"])
    return sorted(out)


def allele_contigs(bam: str, *, lengths=None) -> list:
    """The scaffold contig of each allele: from `lengths` when known, else from the read modes. Impure."""
    counts = contig_counts(bam)
    if lengths:
        got = contigs_for_lengths(counts, lengths)
        if got:
            return got
    return pick_allele_contigs(counts)
