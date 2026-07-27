"""S2 — targeted dupC scan, POSITION-AGNOSTIC (per read, phased). No alignment, no consensus.

The canonical dupC (27dupC = 59dupC) = the **terminal C-tract** of a unit goes from **7C (WT) -> 8C
(dupC)**. Gene-oriented context (cf `MUC1_Analyzer.KNOWN_REPEATS`, units B/X/D/5/aA... ending in
`...GGGCTCCACCGC[7C]AGCCCAC...`) :

        GGGCTCCACCG (C+) AGCCCAC        # captures the terminal C-tract, 7=WT, 8=dupC

We scan EACH read (`query_sequence`, reference-oriented, + its reverse-complement to handle the
strand — MUC1 is on the − strand) for this context, record the **length of the C-tract** at each occurrence,
aggregated by HP. The dupC is DEEP in the array -> local detection, wherever it is. **Paired null** :
compare the **mutant** allele to the patient's **internal HEALTHY** allele (same context, same run/chemistry) ->
an excess of 8C on the mutant vs the healthy = signal, calibrated by the ONT homopolymer noise of the healthy.
Robust at low depth; does not replace S3 (Bayesian) but supplies its counts.
"""
from __future__ import annotations
import collections
import re

import pysam

from ..config import GRCh38

# Gene-oriented context of the terminal C-tract of dupC-able units. (C+) = C-tract (WT 7, dupC 8).
_CTX = re.compile(r"GGGCTCCACCG(C+)AGCCCAC")
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def ctract_runs(seq: str) -> list:
    """Lengths of the terminal C-tracts found in `seq` (both orientations)."""
    return [len(m.group(1)) for m in _CTX.finditer(seq)] + \
           [len(m.group(1)) for m in _CTX.finditer(_rc(seq))]


def _n_ctx(s: str) -> int:
    return len(_CTX.findall(s))


def gene_oriented(seq: str) -> str:
    """Return the sequence in GENE orientation (the one where the dupC context matches)."""
    return seq if _n_ctx(seq) >= _n_ctx(_rc(seq)) else _rc(seq)


# ── FLANK anchors (sequence) for the ALIGNMENT-FREE positional scan ──────────────────────────
# Mirror of `vntr_raw_length.py` (same unique T2T just-outside-tandem flanks). Purpose: run the dupC
# positional scan on ANY bam (T2T 'muc1win', hg38, or UNPHASED) by anchoring on the flanks
# IN the read's sequence — without depending on the hg38 genomic anchors (`_VNTR_ANCHORS`) nor on HP tags.
# Essential for URINE (unphasable -> the HP-dependent caller-context is off-limits).
_AL = "GCTCTTGCTGGCTGGGGTTGTGGTAGCCCTGGCAGAGGTGCCGTTGTGCACCAGAGTAGAA"   # below the tandem (ref+)
_AH = "GCTGGCCTGGTGACTGGGACCGAGGTGACATCCTGTCCCCAGGTGGCAGCTGAACCTGAAG"   # above the tandem (ref+)


def _fuzzy_find(probe: str, seq: str, maxmm: int, k: int = 11):
    """Best position (min Hamming) of `probe` in `seq` via exact k-mer seeds + verification.
    None if the best match exceeds `maxmm`. (Identical to vntr_raw_length.fuzzy_find.)"""
    L = len(probe)
    best = None
    seen = set()
    for so in range(0, L - k + 1, k):
        seed = probe[so:so + k]
        start = 0
        while True:
            j = seq.find(seed, start)
            if j < 0:
                break
            start = j + 1
            p = j - so                      # implied probe start
            if p < 0 or p + L > len(seq) or p in seen:
                continue
            seen.add(p)
            mm = sum(a != b for a, b in zip(probe, seq[p:p + L]))
            if mm <= maxmm and (best is None or mm < best[1]):
                best = (p, mm)
    return best


def array_between_flanks(seq: str, maxmm: int = 9):
    """Sub-sequence of the tandem between the two unique flanks, in GENE orientation (ready for `_CTX`).
    Tries both orientations of the read (− strand); None if the read does not cover both flanks."""
    for s in (seq, _rc(seq)):
        l = _fuzzy_find(_AL, s, maxmm)
        h = _fuzzy_find(_AH, s, maxmm)
        if l and h:
            a, b = l[0] + len(_AL), h[0]
            if b - a > 0:
                return gene_oriented(s[a:b])
    return None


# GENE-5' flank (= rc(AH)) : AH is "above" the tandem in ref+ = HIGH coord; MUC1 is on the − strand
# so high coord = 5' of the gene (cf. memory muc1-t2t-geometry-minus-strand). In gene orientation, this
# 5' flank therefore appears as rc(AH), JUST before the 1st unit. The 59dupC is 5'-proximal (start of
# the array) -> anchoring HERE lets us index the C-tracts with PARTIAL reads (5' flank + first units
# only, without requiring the low flank) -> far more coverage than a complete flank-to-flank span.
_GENE5 = _rc(_AH)


def ctracts_from_gene5(seq: str, maxmm: int = 9):
    """C-tracts indexed from the GENE-5' flank (`_GENE5`) downstream, PARTIAL reads allowed (only the
    5' flank is required). Index 0 = 1st C-tract after the 5' flank. None if the 5' flank is absent.
    Consistent index frame between reads (same anchor) -> a heterozygous dupC = 1 enriched position."""
    g = gene_oriented(seq)
    a = _fuzzy_find(_GENE5, g, maxmm)
    if a is None:
        return None
    arr = g[a[0] + len(_GENE5):]
    return [len(m.group(1)) for m in _CTX.finditer(arr)]


_GENE3 = _rc(_AL)                                            # GENE-3' flank (rc(AL), bottom of the tandem in ref+)


def ctracts_from_gene3(seq: str, maxmm: int = 9):
    """C-tracts indexed from the GENE-3' flank (`_GENE3`) UPSTREAM (reversed index = distance to the
    3' flank). Reads partial on the 3' side allowed (only the 3' flank required) — complementary to gene5 :
    together they cover reads truncated on BOTH sides. None if the 3' flank is absent."""
    g = gene_oriented(seq)
    a = _fuzzy_find(_GENE3, g, maxmm)
    if a is None:
        return None
    arr = g[:a[0]]                                           # everything preceding the 3' flank
    return [len(m.group(1)) for m in _CTX.finditer(arr)][::-1]   # index from the 3' flank


def _runs_pair(seq: str, maxmm: int):
    a = array_between_flanks(seq, maxmm)
    return None if a is None else [len(m.group(1)) for m in _CTX.finditer(a)]


_ANCHOR_FN = {"pair": _runs_pair, "gene5": ctracts_from_gene5, "gene3": ctracts_from_gene3}


# ── UNIT-LENGTH probe (del8_27 & motif-body deletions) ─────────────────────────────────────────
# del8_27 removes ~20 bp from the BODY of the unit but LEAVES the terminal C-tract anchor (`_CTX`)
# intact -> the unit goes from ~60 to ~40 bp. We measure the gap between consecutive `_CTX` starts (=
# unit length) BY POSITION; a heterozygous del = ~one position where half the reads have a SHORT unit.
# Same flank-sequence anchoring (gene5/gene3) as the C-tract -> unphased / T2T / urine OK.
def _unit_gaps(arr: str) -> list:
    """Unit lengths = gaps between consecutive `_CTX` starts in the gene-oriented array."""
    starts = [m.start() for m in _CTX.finditer(arr)]
    return [starts[i + 1] - starts[i] for i in range(len(starts) - 1)]


def unitlens_from_gene5(seq: str, maxmm: int = 9):
    g = gene_oriented(seq)
    a = _fuzzy_find(_GENE5, g, maxmm)
    return None if a is None else _unit_gaps(g[a[0] + len(_GENE5):])


def unitlens_from_gene3(seq: str, maxmm: int = 9):
    g = gene_oriented(seq)
    a = _fuzzy_find(_GENE3, g, maxmm)
    return None if a is None else _unit_gaps(g[:a[0]])[::-1]   # index from the 3' flank


_UNITLEN_FN = {"gene5": unitlens_from_gene5, "gene3": unitlens_from_gene3}


def positional_profile_unitlen_from_seqs(read_seqs, *, anchor: str = "gene5", short_max: int = 50,
                                         min_reads: int = 8):
    """PURE (testable) : per position, fraction of reads whose UNIT is SHORT (<= `short_max`, WT ~= 60 bp,
    del8_27 ~= 40 bp). Returns (profile, n_span) in the same format as the C-tract -> reuses `_decide_positional`."""
    fn = _UNITLEN_FN[anchor]
    pos = collections.defaultdict(collections.Counter)
    n_span = 0
    for seq in read_seqs:
        if not seq:
            continue
        gaps = fn(seq)
        if gaps is None:
            continue
        for i, gap in enumerate(gaps):
            pos[i]["short" if gap <= short_max else "ok"] += 1
        n_span += 1
    profile = []
    for i in sorted(pos):
        tot = pos[i]["short"] + pos[i]["ok"]
        if tot >= min_reads:
            profile.append({"index": i, "tot": tot, "n_mut": pos[i]["short"],
                            "frac": round(pos[i]["short"] / tot, 3)})
    return profile, n_span


def scan_del_positional_seq(bam: str, *, region: str = None, genome_ref: str = None, maxmm: int = 9,
                            short_max: int = 50, min_reads: int = 8, f_null: float = 0.03,
                            alpha: float = 0.001, min_tot: int = 15, min_frac: float = 0.15,
                            min_mut: int = 4, anchor: str = "both") -> dict:
    """Positional UNIT-LENGTH caller (del8_27 / body deletions) — flank-sequence anchoring, unphased/T2T/urine
    OK. `f_null` = expected fraction of falsely short units on WT (ONT noise + missed anchors). `min_frac`
    slightly higher than the C-tract (0.15). Union gene5/gene3, best frame."""
    mode = "rc" if str(bam).endswith(".cram") else "rb"
    kw = {"reference_filename": genome_ref} if (mode == "rc" and genome_ref) else {}
    seqs = []
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        for r in _fetch_region(af, region):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.query_sequence:
                continue
            seqs.append(r.query_sequence)
    frames = ("gene5", "gene3") if anchor == "both" else (anchor,)
    nf = len(frames)
    per_frame = {}
    for fr in frames:
        profile, n_span = positional_profile_unitlen_from_seqs(seqs, anchor=fr, short_max=short_max,
                                                               min_reads=min_reads)
        dec = _decide_positional(profile, f_null=f_null, alpha=alpha / nf, min_tot=min_tot,
                                 min_frac=min_frac, min_mut=min_mut)
        dec.update({"n_span": n_span, "profile": profile, "frame": fr})
        per_frame[fr] = dec
    best = min(per_frame.values(),
               key=lambda d: (not d["called"], (d.get("best") or {}).get("p_bonf", 1.0)))
    return {"best": best.get("best"), "called": best["called"], "n_span": best["n_span"],
            "profile": best["profile"], "probe": "unitlen", "short_max": short_max,
            "anchor": anchor, "frame": best["frame"],
            "frames": {fr: {"n_span": r["n_span"], "called": r["called"], "best": r.get("best")}
                       for fr, r in per_frame.items()}}


def scan_dupc_context(bam: str, *, region: str = None, genome_ref: str = None,
                      hp_mut: str = None, ctx_units: int = 3, min_support: int = 2,
                      mut_len: int = 8, wt_len: int = 7) -> dict:
    """S2-v2 : group C-tracts by local CONTEXT (neighboring units decoded via match_motifs).

    The true variant recurs at ONE context (same neighboring units) with a C-tract of length `mut_len`
    (dupC = 8C ; delCC = 5C); the ONT errors are scattered across varied contexts. We compare, per
    context, mutant vs internal healthy. Warning: the probe reads the **LITERAL C-tract** of the read (regex
    `_CTX`), NOT the `match_motifs` decoding -> the -2 decoding artifact (which drowned the generic caller)
    does not appear here. The keys `n8C`/`n7C`/`hp_8C` remain (compat) but actually count `mut_len`/`wt_len`
    (see `mut_len`)."""
    from .vntr_scaffold import _import_analyzer
    M = _import_analyzer()
    loc = GRCh38["LOCUS"]
    region = region or f"{loc.chrom}:{loc.start}-{loc.end}"
    chrom, coords = region.split(":")
    start, end = (int(x) for x in coords.replace(",", "").split("-"))
    mode = "rc" if bam.endswith(".cram") else "rb"
    kw = {"reference_filename": genome_ref} if (mode == "rc" and genome_ref) else {}

    # per_hp[hp][context] = Counter{ctract_len: n}
    per_hp = {"1": collections.defaultdict(collections.Counter),
              "2": collections.defaultdict(collections.Counter)}
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        for rd in af.fetch(chrom, start, end):
            if rd.is_secondary or rd.is_supplementary or rd.is_unmapped or not rd.query_sequence:
                continue
            hp = str(rd.get_tag("HP")) if rd.has_tag("HP") else None
            if hp not in ("1", "2"):
                continue
            g = gene_oriented(rd.query_sequence)
            for m in _CTX.finditer(g):
                clen = len(m.group(1))
                win = g[max(0, m.start() - 60 * ctx_units): m.start() + 11]   # preceding units
                matched, _ = M.match_motifs(win)
                ctx = "-".join(x["name"] for x in matched[-ctx_units:]) if matched else "?"
                per_hp[hp][ctx][clen] += 1

    # candidates : context where one HP has C-tracts of `mut_len` (support >= min) and the other HP mostly `wt_len`
    cand = []
    for ctx in set(per_hp["1"]) | set(per_hp["2"]):
        n8 = {hp: per_hp[hp][ctx].get(mut_len, 0) for hp in ("1", "2")}
        n7 = {hp: per_hp[hp][ctx].get(wt_len, 0) for hp in ("1", "2")}
        for hp in ("1", "2"):
            other = "2" if hp == "1" else "1"
            if n8[hp] >= min_support and n8[hp] > n8[other]:
                cand.append({"context": ctx, "hp_8C": hp, "n8C": n8, "n7C": n7,
                             "enrichment": n8[hp] - n8[other]})
    cand.sort(key=lambda c: (-c["n8C"][c["hp_8C"]], -c["enrichment"]))
    return {"hp_mut_expected": hp_mut, "mut_len": mut_len, "wt_len": wt_len,
            "n_contexts": {hp: len(per_hp[hp]) for hp in ("1", "2")},
            "dupC_candidates": cand[:8]}


def _decide_positional(profile: list, *, f_null: float = 0.02, alpha: float = 0.001,
                       min_tot: int = 15, min_frac: float = 0.10, min_mut: int = 4) -> dict:
    """PURE decision (testable) from a {index, tot, n_mut, frac} profile per position.

    A true heterozygous dupC = ONE position where the 8C fraction SIGNIFICANTLY exceeds the homopolymer
    noise `f_null`, **after correction for the number of positions tested** (we take the max over
    ~40 positions -> Bonferroni mandatory, otherwise a negative/artifact concentrates by chance : e.g. 4/36,
    4/12). Minimal support per position (`min_tot`) -> kills low-coverage fluctuations
    (e.g. tot=12). One-tailed binomial test (reuses `dupc_power._tail`)."""
    from ..dupc_power import _tail
    tested = [p for p in profile if p["tot"] >= min_tot]
    n_tests = max(1, len(tested))
    scored = []
    for p in tested:
        pv = _tail(p["tot"], p["n_mut"], f_null)
        scored.append({**p, "p": pv, "p_bonf": min(1.0, pv * n_tests)})
    best = min(scored, key=lambda p: (p["p_bonf"], -p["frac"]), default=None)
    called = bool(best and best["p_bonf"] <= alpha
                  and best["frac"] >= min_frac and best["n_mut"] >= min_mut)
    return {"best": best, "called": called, "n_tests": n_tests}


def scan_dupc_positional(bam: str, *, genome_ref: str = None, mut_len: int = 8, margin: int = 60,
                         min_reads: int = 8, f_null: float = 0.02, alpha: float = 0.001,
                         min_tot: int = 15, min_frac: float = 0.10, min_mut: int = 4) -> dict:
    """HP-AGNOSTIC dupC fallback, POSITION-ANCHORED — for phasing failure (length-homozygous).

    whatshap does not separate two alleles of the SAME length (44|44) -> it dumps ~90% of the reads into one
    HP and starves the other (garbage len); the dupC becomes a DROWNED sub-population and the per-HP methods
    (paired caller, analyzer consensus) all miss it. Here we IGNORE the HP : we anchor each read on the unique
    flanks (`l_anchor`/`r_anchor`, like `measure_vntr_lengths`), enumerate the C-tracts IN ORDER (index =
    repeat position, consistent across reads because gene-oriented), and aggregate the fraction of `mut_len`
    (8C) BY POSITION. A heterozygous dupC = ~one position with a high fraction; the artifact/negative stays
    flat. `_decide_positional` decides."""
    from .vntr import _qpos_at_ref, _VNTR_ANCHORS
    a = _VNTR_ANCHORS
    chrom, l_anchor, r_anchor = a["chrom"], a["l_anchor"], a["r_anchor"]
    mode = "rc" if str(bam).endswith(".cram") else "rb"
    kw = {"reference_filename": genome_ref} if (mode == "rc" and genome_ref) else {}
    pos = collections.defaultdict(collections.Counter)      # index -> Counter{ctract_len}
    n_span = 0
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        for r in af.fetch(chrom, l_anchor - 200, r_anchor + 200):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or not r.query_sequence:
                continue
            if r.reference_start > l_anchor - margin or (r.reference_end or 0) < r_anchor + margin:
                continue                                    # must cover both anchors (full span)
            qL, qR = _qpos_at_ref(r, l_anchor), _qpos_at_ref(r, r_anchor)
            if qL is None or qR is None or not (0 < qR - qL <= 12000):
                continue
            g = gene_oriented(r.query_sequence[qL:qR])      # array only, consistent gene orientation
            for i, m in enumerate(_CTX.finditer(g)):
                pos[i][len(m.group(1))] += 1
            n_span += 1
    profile = []
    for i in sorted(pos):
        tot = sum(pos[i].values())
        if tot >= min_reads:
            profile.append({"index": i, "tot": tot, "n_mut": pos[i].get(mut_len, 0),
                            "frac": round(pos[i].get(mut_len, 0) / tot, 3)})
    out = _decide_positional(profile, f_null=f_null, alpha=alpha, min_tot=min_tot,
                             min_frac=min_frac, min_mut=min_mut)
    out.update({"n_span": n_span, "profile": profile, "mut_len": mut_len})
    return out


def _fetch_region(af, region: str = None):
    """Iterate the reads of `af` over `region` ('contig:start-end'); otherwise the WHOLE bam (until_eof)."""
    if not region:
        return af.fetch(until_eof=True)
    chrom, coords = region.split(":")
    start, end = (int(x) for x in coords.replace(",", "").split("-"))
    return af.fetch(chrom, start, end)


def positional_profile_from_seqs(read_seqs, *, mut_len: int = 8, maxmm: int = 9,
                                 min_reads: int = 8, anchor: str = "pair"):
    """PURE (testable without BAM) : from read sequences, anchor by FLANK-SEQUENCE, enumerate the
    C-tracts BY POSITION (repeat index, gene-oriented), and return (profile, n_span). A heterozygous
    dupC = ~one position with a high `mut_len` fraction.

    `anchor` : `pair` = complete flank-to-flank span (SPECIFIC but STARVED, ~25 reads/6.6M); `gene5` =
    GENE-5' flank alone (partial 5' reads allowed); `gene3` = GENE-3' flank alone (partial 3' reads). The
    gene5/gene3 frames are complementary -> cover reads truncated on both sides (cf. `anchor=both`
    in `scan_dupc_positional_seq`)."""
    runs_fn = _ANCHOR_FN[anchor]
    pos = collections.defaultdict(collections.Counter)
    n_span = 0
    for seq in read_seqs:
        if not seq:
            continue
        runs = runs_fn(seq, maxmm)
        if runs is None:
            continue
        for i, clen in enumerate(runs):
            pos[i][clen] += 1
        n_span += 1
    profile = []
    for i in sorted(pos):
        tot = sum(pos[i].values())
        if tot >= min_reads:
            profile.append({"index": i, "tot": tot, "n_mut": pos[i].get(mut_len, 0),
                            "frac": round(pos[i].get(mut_len, 0) / tot, 3)})
    return profile, n_span


def scan_dupc_positional_seq(bam: str, *, region: str = None, genome_ref: str = None,
                             mut_len: int = 8, maxmm: int = 9, min_reads: int = 8,
                             f_null: float = 0.02, alpha: float = 0.001, min_tot: int = 15,
                             min_frac: float = 0.10, min_mut: int = 4,
                             anchor: str = "both") -> dict:
    """Positional SEQUENCE-ANCHORED (ALIGNMENT-FREE) dupC fallback — T2T-native variant of
    `scan_dupc_positional`. Works on any bam (T2T 'muc1win', hg38, unphased) because it anchors on the
    flanks IN the read (like vntr_raw_length), without hg38 genomic anchors nor HP tags. It is the ONLY
    dupC path for URINE (unphasable).

    `anchor` : `both` (default) = union of the gene-5' AND gene-3' frames -> MAX coverage (partial reads on
    both sides); we decide each frame and keep the best (called first, then min p_bonf, with
    correction x n_frames to avoid inflating the FP). `gene5`/`gene3`/`pair` = a single frame.
    Warning: specificity is VALIDATED on a true negative (a MUC1-neg sample), not only in theory."""
    mode = "rc" if str(bam).endswith(".cram") else "rb"
    kw = {"reference_filename": genome_ref} if (mode == "rc" and genome_ref) else {}
    seqs = []
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        for r in _fetch_region(af, region):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.query_sequence:
                continue
            seqs.append(r.query_sequence)
    frames = ("gene5", "gene3") if anchor == "both" else (anchor,)
    nf = len(frames)
    per_frame = {}
    for fr in frames:
        profile, n_span = positional_profile_from_seqs(seqs, mut_len=mut_len, maxmm=maxmm,
                                                       min_reads=min_reads, anchor=fr)
        # multi-frame correction : tighten alpha by the nb of frames (Bonferroni over the frames too)
        dec = _decide_positional(profile, f_null=f_null, alpha=alpha / nf, min_tot=min_tot,
                                 min_frac=min_frac, min_mut=min_mut)
        dec.update({"n_span": n_span, "profile": profile, "frame": fr})
        per_frame[fr] = dec
    best = min(per_frame.values(),
               key=lambda d: (not d["called"], (d.get("best") or {}).get("p_bonf", 1.0)))
    out = {"best": best.get("best"), "called": best["called"], "n_span": best["n_span"],
           "profile": best["profile"], "mut_len": mut_len, "anchor": anchor, "frame": best["frame"],
           "frames": {fr: {"n_span": r["n_span"], "called": r["called"], "best": r.get("best")}
                      for fr, r in per_frame.items()}}
    return out


def scan_ctract_hp(bam: str, *, region: str = None, genome_ref: str = None, hp: str = None) -> dict:
    """Histogram of C-tract lengths over the reads (optionally filtered by HP). Per primary read."""
    loc = GRCh38["LOCUS"]
    region = region or f"{loc.chrom}:{loc.start}-{loc.end}"
    chrom, coords = region.split(":")
    start, end = (int(x) for x in coords.replace(",", "").split("-"))
    mode = "rc" if bam.endswith(".cram") else "rb"
    kw = {"reference_filename": genome_ref} if (mode == "rc" and genome_ref) else {}

    hist = collections.Counter()
    n_reads = n_with = 0
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        for rd in af.fetch(chrom, start, end):
            if rd.is_secondary or rd.is_supplementary or rd.is_unmapped:
                continue
            if hp in ("1", "2") and not (rd.has_tag("HP") and str(rd.get_tag("HP")) == hp):
                continue
            seq = rd.query_sequence
            if not seq:
                continue
            n_reads += 1
            runs = ctract_runs(seq)
            if runs:
                n_with += 1
                hist.update(runs)
    total = sum(hist.values())
    return {"hp": hp or "all", "n_reads": n_reads, "n_reads_with_ctract": n_with,
            "n_ctracts": total, "ctract_hist": dict(sorted(hist.items())),
            "n_7C": hist.get(7, 0), "n_8C": hist.get(8, 0),
            "frac_8C": round(hist.get(8, 0) / total, 3) if total else None}


def scan_dupc(bam: str, *, region: str = None, genome_ref: str = None, hp_mut: str = None) -> dict:
    """Scan both HP + all, and compare (paired null). `hp_mut` (if known) marks the mutant allele."""
    per_hp = {hp: scan_ctract_hp(bam, region=region, genome_ref=genome_ref, hp=hp)
              for hp in ("1", "2", None)}
    out = {"hp1": per_hp["1"], "hp2": per_hp["2"], "all": per_hp[None]}
    # excess of 8C of one HP vs the other (paired null) : the mutant should show more 8C than the healthy
    f1, f2 = per_hp["1"]["frac_8C"], per_hp["2"]["frac_8C"]
    if f1 is not None and f2 is not None:
        out["frac_8C_hp1_vs_hp2"] = (f1, f2)
        out["hp_more_8C"] = "1" if f1 > f2 else ("2" if f2 > f1 else "tie")
    out["hp_mut_expected"] = hp_mut
    return out


if __name__ == "__main__":   # e.g. a sample : see if the mutant allele has an excess of 8C vs the healthy
    import argparse
    import json
    ap = argparse.ArgumentParser(prog="muc1_analyzer.detectors.vntr_dupc",
                                 description="S2 : dupC scan (C-tract 7C/8C) per read, phased")
    ap.add_argument("-b", "--bam", required=True)
    ap.add_argument("--genome-ref", default=None)
    ap.add_argument("--hp-mut", choices=["1", "2"], default=None, help="expected mutant HP (annotation)")
    ap.add_argument("--context", action="store_true",
                    help="S2-v2 : group C-tracts by local context (concentrates the true dupC)")
    ap.add_argument("--ctx-units", type=int, default=3,
                    help="nb of neighboring units decoded for the context (higher = RARER context, "
                         "un-dilutes a dupC in a common context; e.g. 5). With --context.")
    ap.add_argument("--min-support", type=int, default=2,
                    help="minimal 8C on a context to retain a candidate. With --context.")
    args = ap.parse_args()
    if args.context:
        res = scan_dupc_context(args.bam, genome_ref=args.genome_ref, hp_mut=args.hp_mut,
                                ctx_units=args.ctx_units, min_support=args.min_support)
    else:
        res = scan_dupc(args.bam, genome_ref=args.genome_ref, hp_mut=args.hp_mut)
    print(json.dumps(res, indent=2, ensure_ascii=False))
