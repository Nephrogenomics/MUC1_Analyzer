#!/usr/bin/env python3
"""allele_lengths.py — the two allele lengths to hand the dupC dispatch, chosen dropout-aware.

WHY THIS MODULE EXISTS. On LR-PCR the dispatch is only as good as the two lengths it is told to scaffold on.
Left to a plain read-mode arbiter, the LONG allele — depleted by preferential PCR dropout — loses to its own
SHORTER slippage by-products, which are more amplified. Measured on a carrier (truth 77 wt | 82 mut): the modal
arbiter returned [44, 77], scaffolded the variant onto the 44 by-products, and reported the carrier at the
wrong length; feeding [77, 82] instead put the dupC on the 82 contig at frac 0.741 (vs 0.16). So the length
CHOICE, not the alignment or the scaffold pick, was the failing link — and this module fixes exactly it.

It measures each read's VNTR copy number ALIGNMENT-FREE (flank AL→AH, +offset), then picks the two allele
lengths by a score that is height × local-contrast × a graded length weight f(L)=L/dominant. The length
weight is the dropout correction: between two peaks longer than the dominant one (e.g. a real long allele vs
a slippage plateau), the LONGER survives even at lower depth, because slippage products are always shorter
than their parent and a long peak reaching that length despite dropout is the more likely true allele. This
is the same picker validated inside rescaffold on [50,82], [44,106] and [68,77]; here it is a
standalone so the direct amplicon→dispatch path (no binning, VAF preserved) can call it.

    from allele_lengths import two_alleles
    lengths = two_alleles("SAMPLE.xavier.bam")            # e.g. [77, 82]
    dispatch_dupc_vntr(bam, ref, lengths=lengths, variant_probe="dupC")
"""
from collections import Counter

import pysam
import vntr_raw_length as V

DEFAULT_OFFSET = 4        # flank AL→AH sits ~4 units inside the array; +4 recovers the physical count
DEFAULT_MAXMM = 9
DEFAULT_SEP = 1           # do NOT fuse alleles that differ by ≥2 copies; Δ≤1 is arbitrated downstream
DEFAULT_MIN_READS = 5     # a real allele can measure as few as 5 flank-spanning reads (JOT 77: 5 reads,
                          # truth-confirmed) — at 6 it was dropped, making a false length-homozygote → n/2
                          # → a spurious del-1 FP. 5 recovers it; measured impact re-checked on the cohort.

ARTEFACT_LENGTHS = (44, 70)   # recurrent PCR/slippage by-product lengths → small score malus (not a blacklist)
ARTEFACT_MALUS = 0.9


def load_reads(bam):
    """{name: seq} primary reads only, each once. Input alignment ignored — only the sequence is used, so a
    vntr BAM, a xavier BAM, a uBAM or a CRAM all work identically."""
    reads = {}
    mode = "rc" if str(bam).endswith(".cram") else "rb"
    with pysam.AlignmentFile(bam, mode, check_sq=False) as af:
        for r in af.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.query_sequence is None:
                continue
            reads.setdefault(r.query_name, r.query_sequence.upper())
    return reads


def measure_lengths(reads, *, offset=DEFAULT_OFFSET, maxmm=DEFAULT_MAXMM):
    """{name: copies} over reads whose flanks span AL→AH; others get no length."""
    out = {}
    for name, seq in reads.items():
        g = V.read_length(seq, maxmm)
        if g is not None:
            out[name] = round((g - V.GAP) / 60) + offset
    return out


def _peaks(hist, *, min_reads=DEFAULT_MIN_READS, sep=5, nbr_frac=0.40):
    """Transparent peak-pick: local maxima ≥ min_reads, thinned by `sep`, tallest first. Catches a minority
    LONG allele that a smoothing/fractional arbiter drops.

    A plain local-maximum test (n ≥ both immediate neighbours) STRUCTURALLY loses one of two true alleles
    one copy apart: the slightly smaller one sits on the taller's slope and is never a candidate (a
    43|44 pair: 620 vs 623 reads — the 620-read allele is dropped for a 3-read difference). So a length also
    qualifies if it is a SUBSTANTIAL Δ1 neighbour — ≥ nbr_frac of its taller immediate neighbour AND ≥
    min_reads — which admits a real adjacent allele (620/623 = 99 %; 129/148 = 87 %) while still
    rejecting a slippage shoulder (a few % of its parent). Whether two kept Δ1 peaks are truly two alleles
    is decided downstream (consensus / divergent-dupC arbitration), not here."""
    if not hist:
        return []
    def is_local_max(L, n):
        return n >= max(hist.get(L + d, 0) for d in range(-1, 2))
    def is_substantial_neighbour(L, n):
        # a Δ1 neighbour of a TALLER peak, carrying ≥ nbr_frac of it → keep as a distinct candidate
        for d in (-1, 1):
            hi = hist.get(L + d, 0)
            if hi > n and n >= nbr_frac * hi:
                return True
        return False
    cand = [(L, n) for L, n in hist.items()
            if n >= min_reads and (is_local_max(L, n) or is_substantial_neighbour(L, n))]
    cand.sort(key=lambda x: (-x[1], x[0]))
    kept = []
    for L, n in cand:
        if all(abs(L - k[0]) >= sep for k in kept):
            kept.append((L, n))
    return sorted(kept)


def _contrast(hist, L, halfwin=3):
    """Peak height / (mean neighbours within ±halfwin + 1). +1 floor so an isolated micro-spike (no
    neighbours) does not score infinite. A real allele is a sharp spike; a slippage plateau ≈ 1."""
    neigh = [hist.get(L + d, 0) for d in range(-halfwin, halfwin + 1) if d != 0]
    m = (sum(neigh) / len(neigh)) if neigh else 0
    return hist[L] / (m + 1.0)


def pick_two(lengths, *, sep=DEFAULT_SEP, min_reads=DEFAULT_MIN_READS, min_second_contrast=3.0):
    """The two allele lengths from per-read copy numbers.

    Every local-maxima peak is scored by height × local-contrast × f(L), f(L)=10^(L/L_maj − 1) where L_maj
    is the MAJORITY (most-read) peak (so the length weight lifts a peak LONGER than the majority — a real
    allele depleted by long-allele dropout — and penalises a peak SHORTER than it — a slippage by-product,
    always shorter than its parent). The two highest-scoring peaks are the alleles. We do NOT assume the
    most-abundant peak is a real allele: on a carrier whose long mutant allele slippage-decays, the MOST
    abundant peak can be the by-product pile (peak 44 is slippage of the 82 mutant; the true
    alleles are 77 and 82, neither of them the mode). Scoring all peaks and taking the top two avoids
    anchoring on a spurious dominant.

    Returns (alleles, detail). One peak → single length (homozygous / same-length is the SNP-phasing case).

    KNOWN LIMIT (length alone cannot fix it). When the MOST abundant peak is itself a slippage pile of a
    long mutant allele (peak 44, n=180, is by-product of the 82 mutant; true alleles 77 & 82), it
    is taller AND sharper than the real alleles, so no length-shape score demotes it — it wrongly takes an
    allele slot, and the depleted 82 loses the second slot to 77. Distinguishing "44 = by-product of 82"
    from "44 = a real short allele" is a SEQUENCE question (the 44 reads carry the 82 haplotype), not a
    length one. That is the step-B SNP-phasing case; this picker returns [44, 77] there and is not expected
    to be right. It IS right when the dominant peak is a true allele (all three cases above pass)."""
    hist = Counter(lengths)
    pk = _peaks(hist, min_reads=min_reads, sep=sep)
    if not pk:
        return [], {"peaks": []}
    if len(pk) == 1:
        L = pk[0][0]
        return [L], {"peaks": pk, "chosen": [L], "note": "single peak"}
    L_maj = max(pk, key=lambda x: x[1])[0]             # majority peak = length pivot for f(L)
    scored = []
    for L, n in pk:
        c = _contrast(hist, L)
        fL = 10.0 ** (L / L_maj - 1.0)                 # exponential dropout/slippage correction
        malus = ARTEFACT_MALUS if L in ARTEFACT_LENGTHS else 1.0   # 44/70 = recurrent by-products → ×0.9
        scored.append((L, n, round(c, 2), round(fL, 3), round(c * n * fL * malus, 1)))
    scored.sort(key=lambda x: -x[4])                   # best score first
    # take the top two, but enforce `sep` between them (a single allele can raise two adjacent maxima)
    chosen = [scored[0]]
    for row in scored[1:]:
        if all(abs(row[0] - c[0]) >= sep for c in chosen):
            chosen.append(row)
        if len(chosen) == 2:
            break
    # GUARD: a second peak is a real allele only if it is a SHARP peak, not a bump in the slippage/chimera
    # tail a big allele drags behind it. Use CONTRAST (height / local neighbourhood), NOT read fraction:
    # dropout leaves a real long allele legitimately low in READS (106@50 vs 44@1509 = 3 %; and
    # 82@12 vs 50@218 = 5 %), so a read-fraction floor would reject exactly those true long alleles. But a
    # real allele — however depleted — is a SHARP spike (106 contrast 21, 82 contrast 12, and
    # 82 contrast ~7), whereas a tail bump is flat (29 contrast 1.9, 32 contrast 1.8). A contrast
    # floor keeps depleted long alleles and rejects the noise tail — the discrimination read fraction can't
    # make. This applies at ALL Δ: a slippage/chimera bump can sit several copies from its parent (29 vs
    # 44 = Δ15, and 32 vs 44 = Δ12), so the guard MUST NOT be restricted to Δ≤1.
    if len(chosen) == 2:
        second_contrast = min(chosen[0][2], chosen[1][2])
        if second_contrast < min_second_contrast:
            keep = max(chosen, key=lambda r: r[4])     # keep the higher-scoring of the two → homozygous
            return [keep[0]], {"peaks": pk, "baseline": L_maj, "scored": scored, "chosen": [keep[0]],
                               "note": f"2nd peak contrast {second_contrast} < {min_second_contrast} "
                                       f"(tail bump, not a sharp allele) → length-homozygous"}
    alleles = sorted(r[0] for r in chosen)
    return alleles, {"peaks": pk, "baseline": L_maj, "scored": scored, "chosen": alleles}

def two_alleles(bam, ref=None, *, offset=DEFAULT_OFFSET, maxmm=DEFAULT_MAXMM,
                sep=DEFAULT_SEP, min_reads=DEFAULT_MIN_READS, min_second_contrast=3.0,
                min_diff=3, nbr_frac=0.40, dupc_div=0.40, outdir=None, return_detail=False):
    """The two allele lengths to hand the dispatch, straight from a BAM.

    With sep=1, pick_two no longer fuses alleles that differ by ≥2 copies, so a close true pair (77|80,
    76|78) is kept. Two candidate lengths are then resolved by their gap Δ:
      · Δ ≥ 2  → two distinct alleles (a single allele's flank-length does not scatter by 2);
      · Δ ≤ 1  → ambiguous — compare the two groups' CONSENSUS at the base level (consensus_arbiter):
                 ≥ min_diff substitutions ⇒ two alleles; else ⇒ one allele (homozygous / measurement
                 scatter). Needs `ref`; if ref is None the pair is kept as-is (prudent, no silent merge).
    When pick_two returns a SINGLE allele, a Δ1 neighbour carrying ≥ nbr_frac of its reads AND a dupC
    fraction diverging by ≥ dupc_div is rescued as a second allele (the 43|44 case, invisible to the
    local-maximum peak test — see below).
    Returns [] if nothing spans; [L] for a single length; [L1,L2] for two. `return_detail=True` adds the trace."""
    reads_qual = load_reads(bam)                       # {name: (seq, qual)} or {name: seq}
    meas = measure_lengths(reads_qual, offset=offset, maxmm=maxmm)   # {name: copies}
    lengths = list(meas.values())
    alleles, detail = pick_two(lengths, sep=sep, min_reads=min_reads, min_second_contrast=min_second_contrast)
    detail["n_reads"] = len(reads_qual)
    detail["n_spanning"] = len(lengths)

    # Δ rule + consensus arbitration for a close pair
    if len(alleles) == 2:
        L1, L2 = sorted(alleles)
        delta = L2 - L1
        detail["delta"] = delta
        if delta >= 2:
            detail["arbitration"] = f"Δ={delta} ≥ 2 → two alleles by length"
        elif ref is not None:
            # build {name: seq} and the two length-groups, then compare consensuses
            seqs = {nm: (r[0] if isinstance(r, tuple) else r) for nm, r in reads_qual.items()}
            names_a = [nm for nm, L in meas.items() if L == L1]
            names_b = [nm for nm, L in meas.items() if L == L2]
            try:
                from consensus_arbiter import same_or_two_alleles
                two, cdet = same_or_two_alleles(
                    seqs, names_a, names_b,
                    f"MUC1_VNTR_{L1}repeats", f"MUC1_VNTR_{L2}repeats", ref,
                    min_diff=min_diff, outdir=outdir, return_detail=True)
                detail["consensus_cmp"] = cdet
                if two:
                    detail["arbitration"] = (f"Δ={delta}, consensus differ "
                                             f"({cdet.get('n_mismatch')} subs ≥ {min_diff}) → two alleles")
                else:
                    keep = L1 if meas and list(meas.values()).count(L1) >= list(meas.values()).count(L2) else L2
                    alleles = [keep]
                    detail["chosen"] = [keep]
                    detail["arbitration"] = (f"Δ={delta}, consensus ~identical "
                                             f"({cdet.get('n_mismatch')} subs < {min_diff}) → one allele")
            except Exception as e:
                detail["arbitration"] = f"Δ={delta}, consensus compare failed ({str(e)[:40]}) → kept both (prudent)"
        else:
            detail["arbitration"] = f"Δ={delta} ≤ 1 but no ref given → kept both (prudent, no silent merge)"
    # Neighbour rescue by DIVERGENT frameshift fraction. The local-maximum test in _peaks needs a length to
    # dominate BOTH immediate neighbours, so a true second allele a few copies from a taller one (43|44:
    # 148 vs 129; 76|78: the 78 carrier rejected on contrast 2.96<3.0) can be structurally invisible or
    # rejected as a "tail bump". So, when a SINGLE allele L was returned, we look at each neighbour L±1..±3
    # that carries ≥ nbr_frac of L's reads, and ask the FUNCTIONAL question the height/contrast test cannot:
    # do the two length-groups carry the frameshift at very different rates? A slippage shoulder has the SAME
    # dupC fraction as its parent (same DNA, mis-replicated); a true carrier allele diverges sharply (
    # 0.00 vs 0.76; 76 0.115 vs 78 0.763). |Δfrac| ≥ dupc_div ⇒ two alleles. Extended to Δ≤3 (the 76|78 case is Δ2)
    # with a STRICTER divergence floor (0.40 vs the 0.30 that was used for Δ1) because the further apart the
    # lengths, the stronger the functional signal we require to justify splitting. Population floor (≥40 % of
    # the dominant peak) applies at every Δ. This can only ADD a second allele when a real frameshift signal
    # separates the groups, so on a variant-negative sample (both fractions ~0) it never fires — it cannot
    # turn a negative positive, only recover the true carrier allele and place the variant on it.
    if ref is not None and len(alleles) == 1:
        L = alleles[0]
        n_L = list(meas.values()).count(L)
        best_nbr = None
        for d in (-3, -2, -1, 1, 2, 3):
            Ln = L + d
            n_n = list(meas.values()).count(Ln)
            if n_n >= max(nbr_frac * n_L, min_reads):
                fL = _contig_dupc_fraction(bam, ref, L)
                fN = _contig_dupc_fraction(bam, ref, Ln)
                if fL is not None and fN is not None and abs(fL - fN) >= dupc_div:
                    cand = (Ln, n_n, round(abs(fL - fN), 3), round(fL, 3), round(fN, 3))
                    if best_nbr is None or n_n > best_nbr[1]:
                        best_nbr = cand
        if best_nbr:
            Ln = best_nbr[0]
            alleles = sorted([L, Ln])
            detail["chosen"] = alleles
            detail["nbr_rescue"] = {"kept": L, "added": Ln, "delta": abs(Ln - L), "n_added": best_nbr[1],
                                    "dupc_div": best_nbr[2], "frac_L": best_nbr[3], "frac_nbr": best_nbr[4]}
            detail["arbitration"] = (f"Δ{abs(Ln - L)} rescue: neighbour {Ln} ({best_nbr[1]} reads) has "
                                     f"divergent frameshift fraction (|Δ|={best_nbr[2]} ≥ {dupc_div}) → two alleles")

    return (alleles, detail) if return_detail else alleles


def _contig_dupc_fraction(bam, ref, L):
    """Raw pooled dupC fraction on the length-L contig (best called/strongest repeat), or None. Used only
    by the Δ1 rescue to compare a length-group's variant load to its neighbour's."""
    try:
        import frameshift_vntr as _FS
        res = _FS.call(bam, ref, f"MUC1_VNTR_{L}repeats", same_length_alleles=False)
    except Exception:
        return None
    if res.get("error") or not res.get("calls"):
        return 0.0                                   # scanned fine, no frameshift → fraction 0 (a clean allele)
    return max((c.get("raw_frac") or 0.0) for c in res["calls"])


if __name__ == "__main__":
    import argparse
    import sys
    ap = argparse.ArgumentParser(description="Print the two allele lengths for the dispatch from a BAM.")
    ap.add_argument("-b", "--bam", required=True)
    ap.add_argument("--offset", type=int, default=DEFAULT_OFFSET)
    ap.add_argument("--sep", type=int, default=DEFAULT_SEP)
    ap.add_argument("--min-reads", type=int, default=DEFAULT_MIN_READS)
    a = ap.parse_args()
    alleles, detail = two_alleles(a.bam, offset=a.offset, sep=a.sep,
                                  min_reads=a.min_reads, return_detail=True)
    print(f"alleles: {alleles}", file=sys.stderr)
    print(f"  {detail['n_spanning']}/{detail['n_reads']} reads span the VNTR", file=sys.stderr)
    if detail.get("peaks"):
        print(f"  peaks (L:n): {[(L, n) for L, n in detail['peaks']]}", file=sys.stderr)
    if detail.get("second_candidates"):
        print("  second-allele candidates (L, n, contrast, f(L), score):", file=sys.stderr)
        for row in detail["second_candidates"][:5]:
            print(f"    {row}", file=sys.stderr)
    print(",".join(str(x) for x in alleles))     # stdout: machine-readable
