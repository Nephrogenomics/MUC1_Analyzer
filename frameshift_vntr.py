#!/usr/bin/env python3
"""frameshift_vntr.py — frameshift caller by MOTIF-LENGTH SHIFT, bounded to the VNTR.

Why this replaces the run-length detector
-----------------------------------------
The shipped dupC detector (`runlen_shift`) measures the LENGTH OF A C-RUN at reference positions and
calls a read a carrier when that run is long. Two failure modes, both seen on the real cohort:

  1. It confuses a C-run lengthened by a true INSERTION (dupC = motif 60→61 bp, a frameshift) with a
     C-run lengthened by anything else (ONT homopolymer noise, a SNP, a local mis-alignment). Its
     denominator filter `rl >= wt_len - keep_slack` also drops the reads whose run is normal, so on a
     negative the denominator collapses to a handful of reads and the fraction is spuriously high
     (measured: a clean negative called at frac 0.80 on n=15 of 170 covering reads).
  2. It does NOT bound the scan to the VNTR, so an insertion artefact sitting in the 3' flank (outside
     the tandem) is counted (measured: a flank insertion at "repeat 73" of a 44-repeat array).

This detector instead asks the biologically correct question: **is a repeat unit 61 bp instead of 60
(insertion) or 59 bp (deletion)?** — i.e. the NET indel (inserted − deleted) inside each repeat, read
by read, from the CIGAR, and ONLY between the VNTR boundaries. That is:

  * SNP-immune by construction (a substitution changes no length);
  * flank-immune (the scan is clamped to the tandem, delimited by the border motifs);
  * general to ALL frameshifts, not just +1 C (58_59insG, 60dupA, 58_59delCC, longer indels …);
  * honest in its VAF: denominator = ALL reads covering the repeat. When the two alleles are the same
    length and were NOT physically separated (a length-homozygote carrier, sequence-identical save the
    variant), both alleles pile on this contig, so the pooled fraction k/n ≈ 0.5 for a real carrier; the
    per-allele VAF is recovered by halving the DENOMINATOR (n/2), i.e. VAF = 2·k/n — see
    `--per-allele-homozygote`. (This DOUBLES the reported VAF; it does not halve it.)

Once a carrier repeat is localised, its consensus 61/59-bp sequence is matched against the motif
dictionary (`KNOWN_REPEATS`, which carries indel motifs like `X-59dupC`, `B-58_59insG`,
`X-58_59delCC`) to NAME the exact variant.

VNTR delimitation
-----------------
The VNTR starts with border motifs 1–5 and ends with 6–9 (distinct sequences from the standard body).
The start is found by the exact sequence of motif 1; the array then runs `n_repeats × 60 bp`. Using
the standard-body anchor alone MISSES motifs 1–5 and shifts every repeat index by 5 — so we anchor on
motif 1.

This is a PROTOTYPE for evaluation. It does not touch the production `runlen_shift`. Validate on the
labelled cohort (carriers must be recalled at the truth repeat index; negatives, incl. the known
run-length false positives, must be rejected) before any production use.
"""
import argparse
import collections
import json
import os
import sys

import pysam

# Border motifs (1–5 open the VNTR, 6–9 close it) — exact 60 bp sequences from the lab motif dictionary.
BORDER_MOTIFS = {
    "1": "AAGGAGACTTCGGCTACCCAGAGAAGTTCAGTGCCCAGCTCTACTGAGAAGAATGCTGTG",
    "2": "AGTATGACCAGCAGCGTACTCTCCAGCCACAGCCCCGGTTCAGGCTCCTCCACCACTCAG",
    "3": "GGACAGGATGTCACTCTGGCCCCGGCCACGGAACCAGCTTCAGGTTCAGCTGCCACCTGG",
    "4": "GGACAGGATGTCACCTCGGTCCCAGTCACCAGGCCAGCCCTGGGCTCCACCACCCCGCCA",
    "5": "GCCCACGATGTCACCTCAGCCCCGGACAACAAGCCAGCCCCGGGCTCCACCGCCCCCCCA",
}
MOTIF = 60


# ── motif dictionary (loaded from the repo's caller if importable, else a minimal built-in) ───────────
def load_known_repeats():
    """Return {sequence: name}. Prefer the lab dictionary from the haplotype caller so indel motifs
    (X-59dupC, B-58_59insG, X-58_59delCC, …) are available for naming; fall back to a tiny built-in."""
    for modname in ("vntr_haplotype_caller_v5_6_hmz_newmotifs", "vntr_haplotype_caller"):
        try:
            mod = __import__(modname)
            kr = getattr(mod, "KNOWN_REPEATS", None)
            if kr:
                return dict(kr)
        except Exception:
            pass
    # minimal fallback: the standard body motif + the common dupC indel motifs
    return {
        "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA": "X",
        "GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCCA": "X-59dupC",
    }


# ── VNTR boundaries ───────────────────────────────────────────────────────────
def vntr_bounds(fa, contig, n_repeats):
    """(start, end) of the VNTR on `contig`. Start = exact position of border motif 1; the array spans
    n_repeats × 60 bp from there. Returns (None, None) if motif 1 is not found."""
    seq = fa.fetch(contig).upper()
    start = seq.find(BORDER_MOTIFS["1"])
    if start < 0:
        return None, None
    return start, start + n_repeats * MOTIF


def contig_n_repeats(contig):
    import re
    m = re.search(r"MUC1_VNTR_(\d+)repeats", contig)
    return int(m.group(1)) if m else None


# ── per-repeat net indel, read by read ────────────────────────────────────────
def per_repeat_shift(bam, contig, vstart, vend, min_mapq=0):
    """Walk every primary read's CIGAR; accumulate the NET indel (inserted − deleted) per repeat inside
    [vstart, vend). Returns:
        shift[rep]        = Counter{net_bp: n_reads}   (only reads whose net for that repeat != 0)
        cover[rep]        = n_reads covering that repeat
        ins_seq[rep][net] = list of inserted sequences (to name the variant later)
        carriers[(rep,net)] = set of read names carrying that net indel at that repeat (for phasing)
    A read 'covers' a repeat if it has at least one aligned (M/=/X) base within that repeat."""
    shift = collections.defaultdict(collections.Counter)
    cover = collections.Counter()
    ins_seq = collections.defaultdict(lambda: collections.defaultdict(list))
    carriers = collections.defaultdict(set)
    with pysam.AlignmentFile(bam, "rb") as af:
        for r in af.fetch(contig):
            if r.is_secondary or r.is_supplementary or r.cigartuples is None:
                continue
            if r.mapping_quality < min_mapq:
                continue
            refpos = r.reference_start
            qpos = 0
            seq = r.query_sequence or ""
            net = collections.Counter()       # net indel per repeat for THIS read
            covered = set()
            ins_here = collections.defaultdict(list)
            for op, length in r.cigartuples:
                if op in (0, 7, 8):                      # M/=/X
                    for i in range(length):
                        p = refpos + i
                        if vstart <= p < vend:
                            covered.add((p - vstart) // MOTIF + 1)
                    refpos += length
                    qpos += length
                elif op == 1:                            # insertion
                    if vstart <= refpos < vend:
                        rep = (refpos - vstart) // MOTIF + 1
                        net[rep] += length
                        ins_here[rep].append(seq[qpos:qpos + length])
                    qpos += length
                elif op == 2:                            # deletion
                    if vstart <= refpos < vend:
                        rep = (refpos - vstart) // MOTIF + 1
                        net[rep] -= length
                    refpos += length
                elif op == 4:                            # soft clip
                    qpos += length
            for rep in covered:
                cover[rep] += 1
            for rep, nb in net.items():
                if nb != 0:
                    shift[rep][nb] += 1
                    carriers[(rep, nb)].add(r.query_name)
                    for s in ins_here.get(rep, []):
                        ins_seq[rep][nb].append(s)
    return shift, cover, ins_seq, carriers


def name_variant(fa, contig, vstart, rep, net, known):
    """Build the observed repeat sequence with the net indel applied and match it (and its length class)
    against KNOWN_REPEATS. Returns a name like 'X-59dupC' when a length-matching indel motif is found,
    else a generic descriptor ('ins+1'/'del-2'…). Best-effort — naming never gates the call."""
    want_len = MOTIF + net
    # candidate indel motifs of the right length
    cands = {seq: nm for seq, nm in known.items() if len(seq) == want_len and "-" in nm}
    if not cands:
        return f"{'ins' if net > 0 else 'del'}{net:+d}"
    # the reference repeat sequence (60 bp) to compare the stem against
    ref_rep = fa.fetch(contig, vstart + (rep - 1) * MOTIF, vstart + rep * MOTIF).upper()
    # score each candidate by identity of its non-indel stem to the reference repeat
    def stem_score(seq):
        return sum(1 for a, b in zip(seq[:MOTIF], ref_rep) if a == b)
    best = max(cands, key=stem_score)
    return cands[best]


def _aggregate_overcalls(net_counts, ont_slack=4, max_size=4):
    """ONT miscounts homopolymers: a true +1 dupC is basecalled +1, but also +2/+3/+4 (and sometimes
    -1). Counting only the exact dominant shift undercounts a real frameshift (measured on a carrier: +1=33%
    but +2/+3 add ~8% more, all the same dupC). So we take the DOMINANT non-zero shift (within
    +-max_size — larger nets are STRUCTURAL indels, a different event, e.g. +18 bp = X-33_34ins, which
    belong to run's motif-dictionary consensus, not this point-frameshift screen) as the variant and
    ABSORB same-sign neighbours within `ont_slack` bp of it.

    Returns (variant_net, k_aggregated): variant_net is the dominant small shift (defines the variant to
    name), k_aggregated sums all reads whose net is same-sign and within [variant_net, variant_net+-slack]."""
    nz = {d: n for d, n in net_counts.items() if d != 0 and abs(d) <= max_size}
    if not nz:
        return 0, 0
    variant_net = max(nz, key=lambda d: nz[d])          # dominant non-zero shift within max_size
    if variant_net > 0:
        lo, hi = variant_net, variant_net + ont_slack   # insertion: absorb larger overcalls
    else:
        lo, hi = variant_net - ont_slack, variant_net   # deletion: absorb larger undercounts
    all_nz = {d: n for d, n in net_counts.items() if d != 0}   # aggregation window may reach beyond max_size
    k_agg = sum(n for d, n in all_nz.items() if lo <= d <= hi and (d > 0) == (variant_net > 0))
    return variant_net, k_agg


def _background_insertion_rate(shift, cover, min_reads):
    """Empirical background: the median per-repeat +1-insertion fraction across the VNTR. ONT makes
    ~1 bp indel errors in the C-runs at a low rate on EVERY repeat; a true frameshift is an OUTLIER
    above this. Returns the median fraction over repeats with enough coverage (excludes the strongest
    peak so a real carrier does not inflate its own null)."""
    fracs = []
    for rep, c in cover.items():
        if c < min_reads:
            continue
        net_counts = shift.get(rep)
        f = (net_counts.get(1, 0) / c) if net_counts else 0.0   # +1 insertion fraction
        fracs.append(f)
    if not fracs:
        return 0.0
    fracs.sort()
    # drop the top one (a carrier) before taking the median, so the null is the true background
    core = fracs[:-1] if len(fracs) > 2 else fracs
    return core[len(core) // 2]


def _binom_sf(k, n, p):
    """P(X >= k) for X~Binom(n,p). Exact for small n; for large n (thousands of reads) the exact sum
    overflows / is slow, so use a normal approximation there. (The exact path was raising 'int too large
    to convert' on a deep contig — ELE.)"""
    if p <= 0:
        return 1.0 if k <= 0 else 0.0
    if p >= 1:
        return 1.0
    if k <= 0:
        return 1.0
    if n > 1000:
        # normal approximation with continuity correction
        import math
        mu = n * p
        sd = math.sqrt(n * p * (1 - p))
        if sd == 0:
            return 1.0 if k <= mu else 0.0
        z = (k - 0.5 - mu) / sd
        return 0.5 * math.erfc(z / math.sqrt(2))
    from math import comb
    return sum(comb(n, i) * p**i * (1 - p)**(n - i) for i in range(k, n + 1))


def _phase_denominator(bam, contig, vstart, vend, rep, net, ref_len, carrier_read_test,
                       pooled_frac):
    """For same-length alleles on one contig, try to split the two haplotypes by SNP and count the
    variant only against the CARRIER haplotype. Returns (mode, k, n):
      · ('phased', k_carrier, n_carrier)  ONLY if phasing REALLY separated the alleles — i.e. the
        carrier haplotype is strongly enriched for the variant AND the other haplotype is strongly
        depleted. Per-allele VAF, no ÷2.
      · ('halved', None, None)   if no usable phasing (too few SNPs, OR the split did NOT separate
        carriers from non-carriers — the sequence-homozygous case, where find_phasing_snps
        picks up noise positions that do not distinguish the two identical alleles). Caller uses n//2.
      · ('pooled', None, None)   on error → caller keeps the raw pooled fraction.

    The separation check is the whole point: on a real length-homozygote-but-heterozygous carrier the
    two alleles are identical in sequence except the dupC, so ANY SNP-based split is spurious and both
    'haplotypes' end up with the SAME pooled carrier fraction (~0.33 there). We detect exactly that and
    refuse to trust the phasing."""
    try:
        from muc1_analyzer.caller import find_phasing_snps, phase_reads
    except Exception:
        return ("pooled", None, None)
    try:
        snps = find_phasing_snps(bam, contig, ref_len)
        if not snps or len(snps) < 3:
            return ("halved", None, None)
        names_a, names_b = phase_reads(bam, contig, ref_len, snps)
        set_a, set_b = set(names_a), set(names_b)
        if not set_a or not set_b:
            return ("halved", None, None)
        k_a = sum(1 for nm in set_a if carrier_read_test(nm))
        k_b = sum(1 for nm in set_b if carrier_read_test(nm))
        fa_frac = k_a / len(set_a)
        fb_frac = k_b / len(set_b)
        hi = max(fa_frac, fb_frac)
        lo = min(fa_frac, fb_frac)
        # REAL separation: carrier haplotype clearly enriched, other clearly depleted. If both
        # haplotypes carry the variant at ~the pooled rate, the split did not separate the alleles
        # (sequence-homozygous) → use n//2 instead.
        if hi >= 0.60 and lo <= 0.20 and (hi - lo) >= 0.40:
            if fa_frac >= fb_frac:
                return ("phased", k_a, len(set_a))
            return ("phased", k_b, len(set_b))
        return ("halved", None, None)
    except Exception:
        return ("pooled", None, None)


# ── calling ───────────────────────────────────────────────────────────────────
def call(bam, ref, contig, *, floor=0.20, called_at=0.35,
         del_floor=0.35, del_called_at=0.50,
         bg_mult=3.0, alpha=0.001, min_cover_igv=8, min_mapq=0, same_length_alleles=False,
         zygosity_uncertain=False, zygosity_depth=10,
         ont_slack=4, max_size=4):
    """Scan one contig for POINT frameshift-carrying repeats (|net| <= max_size; structural indels are
    run's motif-dictionary job, not this screen).

    Verdict on the corrected VAF (thresholds stricter for deletions — more ONT homopolymer artefact):
        VAF >= call-threshold, coverage OK, above background   -> CALLED
        VAF >= call-threshold but low coverage OR not above bg -> POSITIVE_LOW_DEPTH (a flagged POSITIVE,
                                                                  counts as detected, needs depth confirm)
        floor <= VAF < call-threshold, OR not above background -> NEEDS_IGV (grey zone, prudent NEGATIVE)
        VAF < floor                                            -> not reported (NEG)
    (insertion: floor / called_at ; deletion: del_floor / del_called_at)

    VAF / denominator: different-length alleles -> pooled k/n is already per-allele; same-length alleles
    -> phase per-allele when SNPs truly separate the haplotypes, else n//2 (sequence-homozygous).
    No min-read rejection: allele lengths are trusted upstream, so a low-coverage repeat on a real allele
    is genuine signal — it is never dropped, only routed to POSITIVE_LOW_DEPTH / NEEDS_IGV."""
    fa = pysam.FastaFile(ref)
    n_rep = contig_n_repeats(contig)
    if n_rep is None:
        return {"contig": contig, "error": "cannot parse repeat count from contig name"}
    vstart, vend = vntr_bounds(fa, contig, n_rep)
    if vstart is None:
        return {"contig": contig, "error": "border motif 1 not found — cannot delimit VNTR"}
    ref_len = fa.get_reference_length(contig)
    known = load_known_repeats()
    shift, cover, ins_seq, carriers = per_repeat_shift(bam, contig, vstart, vend, min_mapq=min_mapq)

    f_null = _background_insertion_rate(shift, cover, min_cover_igv)

    calls = []
    for rep in sorted(shift):
        c = cover.get(rep, 0)
        if c == 0:
            continue
        # dominant SMALL shift + ONT overcall aggregation (a +1 dupC also shows as +2/+3/+4). Structural
        # indels (|net| > max_size, e.g. +18) are ignored here — they belong to run's motif-dict consensus.
        net, k = _aggregate_overcalls(shift[rep], ont_slack=ont_slack, max_size=max_size)
        if net == 0 or k == 0:
            continue

        raw_frac = k / c
        low_cover = c < min_cover_igv          # too few reads to trust a CALLED — flag for confirmation

        # significance vs the empirical background. At low coverage the binomial has no power, so we do
        # NOT require it (that would silently drop a real 6/6 signal); we route to POSITIVE_LOW_DEPTH.
        p_bg = _binom_sf(k, c, max(f_null, 1e-6))
        above_bg = (p_bg <= alpha) and (raw_frac >= bg_mult * max(f_null, 1e-6))

        # corrected VAF (per-allele / n//2) for same-length alleles. Carrier reads = those whose net is
        # in the same overcall window as the variant (so phasing sees all the true-carrier reads).
        if net > 0:
            win = range(net, net + ont_slack + 1)
        else:
            win = range(net - ont_slack, net + 1)
        carrier_names = set()
        for d in win:
            carrier_names |= carriers.get((rep, d), set())

        vaf_mode = "pooled"
        k_use, n_use = k, c
        if same_length_alleles:
            mode, kc, nc = _phase_denominator(
                bam, contig, vstart, vend, rep, net, ref_len,
                carrier_read_test=lambda nm: nm in carrier_names, pooled_frac=raw_frac)
            if mode == "phased" and nc:
                k_use, n_use, vaf_mode = kc, nc, "phased"
            elif mode == "halved":
                k_use, n_use, vaf_mode = k, max(1, c // 2), "halved(n//2)"
        vaf = k_use / n_use

        # thresholds are stricter for deletions (more ONT artefact)
        is_del = net < 0
        this_floor = del_floor if is_del else floor
        this_call = del_called_at if is_del else called_at

        # below the absolute floor → NEG (not reported).
        if vaf < this_floor:
            continue
        # Three distinct low-confidence outcomes, which must NOT be conflated:
        #   · NEEDS_IGV_FOR_ZYGOSITY (NEGATIVE): the VAF reaches the call threshold ONLY because n//2 was
        #     applied (raw pooled fraction is below threshold), that n//2 came from a length-homozygote
        #     ASSUMPTION whose second peak was rejected on CONTRAST alone (zygosity_uncertain — it may be a
        #     real second allele, as in JOT 77|80), AND coverage is < zygosity_depth so the n//2 rests on
        #     too few reads to trust. A het would give half this VAF and be NEG, so we do not call it — we
        #     flag the zygosity for IGV. (that case is spared: same n//2 bascule, but coverage ≫ zygosity_depth.)
        #   · POSITIVE_LOW_DEPTH (POSITIVE): VAF at the CALLED level but coverage < min_cover_igv OR not
        #     above background — a real variant, flagged for depth confirmation; counts as detected.
        #   · NEEDS_IGV (NEGATIVE): grey zone (floor ≤ VAF < call-threshold) — variant itself uncertain.
        bascule_par_n2 = (vaf_mode == "halved(n//2)") and (raw_frac < this_call)
        if (vaf >= this_call and zygosity_uncertain and bascule_par_n2 and c < zygosity_depth):
            verdict = "NEEDS_IGV_FOR_ZYGOSITY"
        elif vaf >= this_call and (low_cover or not above_bg):
            verdict = "POSITIVE_LOW_DEPTH"
        elif low_cover or not above_bg:
            verdict = "NEEDS_IGV"
        elif vaf >= this_call:
            verdict = "CALLED"
        else:
            verdict = "NEEDS_IGV"
        # Border-motif guard: a real ADTKD frameshift sits in the BODY of the tandem. The opening border
        # motifs (repeats 1–5) and the closing border motifs (the last 4, n_rep−3..n_rep) are distinct
        # sequences where a called indel is an edge/alignment artefact, not a carrier signal — measured:
        # the only repeat-1 calls in the cohort were false positives, every true carrier was at rep ≥ 11.
        # Such a call is not a clean negative either (something is there), so it is flagged NEEDS_IGV, not
        # dropped — a prudent negative for the Sp objective, with a visual-review flag.
        in_border = (rep <= 5) or (rep >= n_rep - 3)
        if in_border and verdict in ("CALLED", "POSITIVE_LOW_DEPTH"):
            verdict = "NEEDS_IGV"
        calls.append({
            "repeat": rep, "net_bp": net, "k": k_use, "n_cover": n_use,
            "vaf": round(vaf, 3), "vaf_mode": vaf_mode, "verdict": verdict,
            "k_pooled": k, "n_pooled": c, "raw_frac": round(raw_frac, 3),
            "low_cover": low_cover, "in_border": in_border,
            "f_null": round(f_null, 4), "p_vs_bg": p_bg,
            "variant": name_variant(fa, contig, vstart, rep, net, known),
            "kind": "insertion" if net > 0 else "deletion",
        })
    # ranking for picking the strongest call. NEEDS_IGV_FOR_ZYGOSITY and NEEDS_IGV are both NEGATIVE
    # outcomes (ranked above NEG only so a flagged repeat is chosen over a silent one for the IGV note).
    order = {"CALLED": 4, "POSITIVE_LOW_DEPTH": 3, "NEEDS_IGV": 2, "NEEDS_IGV_FOR_ZYGOSITY": 1, "NEG": 0}
    # sort by verdict rank first, then VAF — so calls[0] is the strongest VERDICT, not merely the highest
    # VAF (a border call downgraded to NEEDS_IGV must not shadow a real CALLED in the body).
    calls.sort(key=lambda d: (-order.get(d["verdict"], 0), -d["vaf"]))
    best_verdict = calls[0]["verdict"] if calls else "NEG"
    # PER-SAMPLE zygosity cap: the n//2 divides the denominator of EVERY position on this contig, so if any
    # repeat is NEEDS_IGV_FOR_ZYGOSITY (zygosity uncertain + low depth + bascule-by-n//2), the whole scan
    # rests on that unreliable n//2 — no call from it may be POSITIVE. Cap the sample verdict at
    # NEEDS_IGV_FOR_ZYGOSITY (a prudent NEGATIVE), overriding any CALLED/POSITIVE_LOW_DEPTH here.
    if any(cl["verdict"] == "NEEDS_IGV_FOR_ZYGOSITY" for cl in calls) and \
       best_verdict in ("CALLED", "POSITIVE_LOW_DEPTH"):
        best_verdict = "NEEDS_IGV_FOR_ZYGOSITY"
    # both CALLED and POSITIVE_LOW_DEPTH are POSITIVE calls (the latter flagged for depth confirmation)
    is_positive = best_verdict in ("CALLED", "POSITIVE_LOW_DEPTH")
    return {"contig": contig, "vntr_start": vstart, "vntr_end": vend, "n_repeats": n_rep,
            "same_length_alleles": same_length_alleles, "f_null": round(f_null, 4),
            "called": is_positive, "verdict": best_verdict,
            "carrier_repeat": calls[0]["repeat"] if calls else None,
            "calls": calls}


def main(argv=None):
    ap = argparse.ArgumentParser(description="Frameshift caller by motif-length shift, bounded to the VNTR.")
    ap.add_argument("-b", "--bam", required=True, help="VNTR-aligned BAM (e.g. the muc1 align profile)")
    ap.add_argument("-r", "--ref", required=True, help="multi-contig VNTR reference (…Kirby.fa)")
    ap.add_argument("-c", "--contig", action="append", required=True,
                    help="contig to scan, e.g. MUC1_VNTR_44repeats (repeatable)")
    ap.add_argument("-s", "--sample", default="sample")
    ap.add_argument("--floor", type=float, default=0.20, help="insertion NEEDS_IGV floor")
    ap.add_argument("--called-at", type=float, default=0.35, help="insertion CALLED threshold")
    ap.add_argument("--del-floor", type=float, default=0.35, help="deletion NEEDS_IGV floor")
    ap.add_argument("--del-called-at", type=float, default=0.50, help="deletion CALLED threshold")
    ap.add_argument("--bg-mult", type=float, default=3.0, help="raw fraction must exceed bg_mult x background")
    ap.add_argument("--alpha", type=float, default=0.001, help="binomial significance vs background")
    ap.add_argument("--min-cover-igv", type=int, default=5, help="coverage below this -> POSITIVE_LOW_DEPTH/NEEDS_IGV")
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--same-length", action="store_true",
                    help="the two alleles are the same length (one contig): apply per-allele phasing / n//2")
    ap.add_argument("--ont-slack", type=int, default=4, help="ONT overcall aggregation window (bp)")
    ap.add_argument("--max-size", type=int, default=4, help="max |net| treated as a point frameshift")
    ap.add_argument("--json", default=None, help="write full result JSON here")
    a = ap.parse_args(argv)

    result = {"sample": a.sample, "contigs": []}
    for contig in a.contig:
        res = call(a.bam, a.ref, contig, floor=a.floor, called_at=a.called_at,
                   del_floor=a.del_floor, del_called_at=a.del_called_at,
                   bg_mult=a.bg_mult, alpha=a.alpha, min_cover_igv=a.min_cover_igv,
                   min_mapq=a.min_mapq, same_length_alleles=a.same_length,
                   ont_slack=a.ont_slack, max_size=a.max_size)
        result["contigs"].append(res)
        if res.get("error"):
            print(f"[frameshift] {contig}: {res['error']}", file=sys.stderr)
        elif res["calls"]:
            for cll in res["calls"]:
                print(f"[frameshift] {a.sample} {contig}: repeat {cll['repeat']} "
                      f"{cll['variant']} ({cll['net_bp']:+d}bp) [{cll['verdict']}] "
                      f"VAF={cll['vaf']} ({cll['vaf_mode']}, {cll['k']}/{cll['n_cover']}; "
                      f"pooled {cll['raw_frac']})", file=sys.stderr)
        else:
            print(f"[frameshift] {a.sample} {contig}: no frameshift in VNTR", file=sys.stderr)
    print(f"[frameshift] {a.sample}: verdict = "
          + ", ".join(f"{r['contig'].split('_')[-1]}:{r.get('verdict','ERR')}" for r in result["contigs"]),
          file=sys.stderr)
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(result, fh, indent=2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
