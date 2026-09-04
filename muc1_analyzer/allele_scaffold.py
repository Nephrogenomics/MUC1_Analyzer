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


def contigs_for_lengths(counts: dict, lengths, *, window: int = 5) -> list:
    """Scaffold contigs for KNOWN allele lengths — the read MODE near each length, not the nearest name.

    Picking `min(abs(k - L))` over any contig carrying at least one read was wrong on a tandem, and the
    cost is measured. ONT length calls SMEAR across neighbouring contigs, so a target of 82 finds 78
    holding the allele's 139 reads while 79-82 each catch a handful — and the old rule took 81 with 24
    reads because it was one copy closer. On SER that turned `k=13/98, p_bonf 1.7e-05` into
    `k=5/24, p_bonf 0.045`: a 2700x loss of significance, from choosing a name over the data.

    The smear widens with array length, which is why the sensitivity fell with it: 0.80 on mutant alleles
    under 50 copies against 0.41 at 70+.

    Within `window` copies of the target the READ COUNT decides; outside it, the old nearest-with-reads
    rule still applies. A contig already taken by the first allele is excluded from the second, so two
    alleles that sit close together still get two scaffolds instead of collapsing onto one mode. Pure."""
    avail = {}
    for name, n in counts.items():
        m = CONTIG_LEN.search(name)
        if m and n > 0:
            avail[int(m.group(1))] = (name, n)
    if not avail:
        return []
    out, seen_targets = [], []
    for L in lengths or []:
        L = int(L)
        if L in seen_targets:
            continue            # length-homozygous: ONE allele, so one scaffold — never the same twice
        seen_targets.append(L)
        free = [k for k in avail if avail[k][0] not in out]
        if not free:
            break
        near = [k for k in free if abs(k - L) <= window]
        pool = near or free
        best = max(pool, key=lambda k: (avail[k][1], -abs(k - L), -k))
        out.append(avail[best][0])
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


#: The arbiter must have measured BOTH alleles for its lengths to constrain the scaffold. One length is
#: not a measurement — it comes with verdict INSUFFICIENT — and treating it as a constraint trusts exactly
#: what the arbiter declared unreliable.
MIN_LENGTHS_TO_CONSTRAIN = 2


def usable_lengths(lengths) -> list:
    """The arbiter lengths that may constrain the scaffold — [] unless there are at least two. Pure.

    ⚠ MEASURED 2026-08-07 on the 152-sample cleaned arm. A SINGLE arbiter length still narrowed the pick
    to +/-5 around it, and on that window `contigs_for_lengths` can take a THINNER contig than the read
    mode: one subject's arbiter said [84], the constrained pick took 79repeats, the read mode was 81repeats, and
    the caller found nothing on the thinner one. Dropping the constraint in that case recovered 2 carriers
    (two subjects) AND removed the cohort's only false positive — better on both axes, which a threshold move
    never is. The two-length constraint itself stays: it is what stopped the 5C debacle, where reads from
    both alleles piled onto one contig and halved the mutant fraction."""
    vals = [int(x) for x in (lengths or []) if isinstance(x, (int, float)) and not isinstance(x, bool)]
    return vals if len(vals) >= MIN_LENGTHS_TO_CONSTRAIN else []


def allele_contigs(bam: str, *, lengths=None) -> list:
    """The scaffold contig of each allele: from `lengths` when the arbiter measured BOTH, else from the
    read modes. Impure."""
    counts = contig_counts(bam)
    usable = usable_lengths(lengths)
    if usable:
        got = contigs_for_lengths(counts, usable)
        if got:
            return got
    return pick_allele_contigs(counts)
