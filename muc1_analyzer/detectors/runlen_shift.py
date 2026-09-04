#!/usr/bin/env python3
"""runlen_shift — homopolymer frameshift caller by RUN-LENGTH SHIFT, calibrated in-sample.

Extracted from `MUC1_Analyzer_fromfastq.py` (branch `muc1-analyzer-fromfastq`) so the statistic can be
tested and put through the specificity harness independently of that script's FASTQ pipeline.
Now WIRED into `dupc_dispatch` (per-allele, VNTR space) and measured: specificity **1.000 on 92 independent
non-carriers**, against the 0.66 our first frequentist dupC caller scored on that same set before it was
hardened. The caveat that remains is the TERRAIN, not the statistic — see "the terrain caveat" below.

## The idea worth keeping

ONT cannot count a homopolymer: a TRUE 8C run is basecalled 6, 7, 8 or 9. Counting reads showing *exactly*
`mut_len` C's therefore keeps only ~1/3 of a carrier's reads — our own "3/9" puzzle — and throws away the
9's, which are over-called 8's. The right statistic is the fraction of reads whose C-tract is **>= mut_len**.

## The null

Compared against an IN-SAMPLE empirical null: how this very sample's reads render its OTHER C-runs of the
same reference length. That adapts to the patient's basecaller and needs no PoN — the real usability win
for anyone who does not have our 1000G panel.

Two corrections applied to the extracted version:

* **Leave-one-out.** The null is computed EXCLUDING the position under test. Pooling it in (the original)
  lets a true carrier's own mutant reads inflate the null it is tested against — conservative, so it costs
  sensitivity rather than specificity, but it is still wrong and it is free to fix.
* **The null must rest on enough runs.** With a single other C-run the "empirical null" is one position's
  noise; `min_null_obs` refuses to call rather than test against a null it cannot estimate.

## The pooled null — the risk, and where it actually sits

The null is POOLED over the other runs, i.e. it assumes the ONT error rate is homogeneous across contexts.
Our PoN work says `f` varies by context, so this is the method's load-bearing assumption. It has now been
measured rather than argued.

**MEASURED** (`spec_runlen_shift.py`, 92 independent non-carriers, GRCh38 window, 34 positions):
specificity **1.000 (0 FP)**, and the per-position null spans **f = 0.0081 → 0.0286, median 0.0223 —
a 4x spread, all of it low**.

**Simulated at that measured heterogeneity**: 0.2 % false calls at depth 20 and 0 % at 40/80/150.
The pooled null is therefore defensible here. The failure mode is real but needs a contrast the data
does not show — injecting one hot position:

    f_hot   x median   FP at depth 40
    0.06      2.7x          0 %
    0.10      4.5x        0.5 %
    0.15      6.7x        5.5 %
    0.25     11.2x       55.2 %

So the method breaks only once a context is ~5x noisier than the median; the noisiest position observed is
1.3x. That is a wide margin, and an earlier version of this docstring quoted the 0.25 row as if it were the
expected regime — it is 11x anything measured.

**The terrain caveat, now with numbers.** The VNTR is COLLAPSED on GRCh38 (`C{7}A` = 0 matches in that
window), so those 34 runs are not the tandem contexts. Measured in VNTR space on scaffold contigs
(`sens_runlen_shift.py`, job 1916470):

* **Sensitivity 2/2** on known carriers, p_bonf 9.5e-07 and 1.6e-16, both at the same run index —
  the `>= mut_len` statistic works on the real terrain.
* But **the tandem is noisier**: f_null 0.030-0.070 there versus 0.022 genome-side, i.e. up to ~3x.
* And the margin on negatives collapses. In the genomic sweep the best non-carrier sat at p_bonf 0.11;
  in the tandem one non-carrier reached **p_bonf 2.2e-03 with frac 0.242** — a factor 2.2 from the
  alpha=0.001 gate and 0.008 from the then-0.25 frac gate. Two gates each held by a hair, on a negative.
  (Since 2026-07-28 `min_frac` is 0.20, so that negative is held by ALPHA alone — see the note on
  DEFAULTS for what the move was measured to cost, and the condition to revert it.)

So the genomic **specificity of 1.000 does NOT transfer to the tandem**, and it should not be quoted as
if it did. A tandem-space specificity needs a proper set of non-carrier scaffold BAMs (the cohort's PCR
negatives, put through `prepare`), not the two available here.

Where a context-resolved PoN exists (`dupc_pon.f_for_context`) it remains the better null; this caller is
for when none does.
"""
from __future__ import annotations

import collections

# `min_frac` 0.25 -> 0.20 on 2026-07-28, PROVISIONALLY. What it costs, measured, not reasoned:
#   · the 92 INDEPENDENT genomic non-carriers: 0 false positive at 0.20 (83 re-measured directly; the
#     other 9-10 could not be re-read but sit at frac 0.038-0.111, far under the gate, with p_bonf
#     0.149-1.0 — provably immune). Specificity stays 1.000.
#   · the cohort's 15 tandem negatives: 0 false positive at 0.20.
#   · what it BUYS: a real carrier rejected at 0.25 with p_bonf 6.7e-07 and k=10, because its frac was
#     0.227. The non-carrier that motivated the 0.25 gate sat HIGHER (frac 0.242), so the fraction never
#     separated those two — the p-value does, by 3300x, and alpha still blocks that negative at 0.20
#     (2.2e-03 > 1e-03). The defence is unchanged; only its redundant second line moved.
# ⚠ PROVISIONAL, and the reason is a real gap: the independent negatives are GENOMIC, while the tandem is
#   noisier (f_null 0.030-0.070 against 0.022). The tandem negatives that say 0 FP here are the cohort's
#   own, so they are NOT independent of the question. REVERT to 0.25 unless fresh PCR negatives confirm
#   it on the tandem terrain. Decision LM 2026-07-28; see docs/MUC1_log.md and decisions.md.
# `min_cover` is `dupc_power.min_reads(d=0.3, target=0.90)` = 12 — the caller ALREADY computes it and
# already publishes `adequately_powered` from it, but until 2026-08-04 it used it only to qualify a
# NEGATIVE. That asymmetry was the defect: a scaffold too thin to believe a negative is too thin to
# believe a positive. Measured on the cohort — a subject whose clinical alleles are 43|44 got a phantom
# 68-copy scaffold from the arbiter and was CALLED there on k=4 of 8 reads (p_bonf 8.9e-04 against an
# alpha of 1e-03), the only false positive in 106 negatives. And it costs nothing: 0 of the 23 called
# carriers rests on a candidate under 12 reads. Not a patch fitted to one case — the same number,
# applied to both verdicts. `tests/test_runlen_shift.py` pins it against `dupc_power`.
# `alpha` 0.001 -> 0.005 on 2026-08-04, on LM's decision, and priced on the POOLED homogeneous arms
# (`pcr_lot180_t2t_clean_2s` + `pcr_lot180_raw_2s`, 177 subjects: 58 carriers / 119 negatives, of which
# 13 are patients disjoint from the arm the threshold was chosen on):
#     alpha 0.001  ->  29/58 = 0.500, specificity 1.000 (119/119)
#     alpha 0.005  ->  31/58 = 0.534, specificity 1.000 (119/119)
#     alpha 0.01   ->  31/58 = 0.534, specificity 1.000 (119/119)
# The reason to move is the PLATEAU rather than the +2: 0.005 and 0.01 give identical counts, so the
# threshold is not perched on a cliff and would not flip on one more negative. It also stays stringent —
# Bonferroni is already applied, and 0.005 is 10x tighter than the conventional 0.05. The first false
# positive appears somewhere between 0.01 and 0.05 (alpha 0.05 cost one negative on the T2T-clean arm).
# ⚠ The genuinely INDEPENDENT slice is only 13 subjects: 0/13 does not exclude an FP rate near 20 % at the
# top of the interval. Revisit with a larger independent tandem negative set. See [[muc1-arbiter-no-tuning]].
# `min_frac` 0.20 -> 0.10 on 2026-08-04, LM's decision, and the reason is a REGIME, not a better score.
# Traced on the long mutant alleles (>= 70 copies, 63 % of the cohort, and the whole of the sensitivity
# gap): their mutant fraction is structurally compressed to 0.10-0.23 — ONT drops a true 8C to 6C/7C, and
# the scaffold pool is diluted — while on short alleles it reaches 0.24-0.83. So at 0.20 this gate stops
# being a noise floor and starts cutting into signal. Forcing the clinical lengths on the missed long
# carriers showed what it costs:
#     OUE  contig 63  k=14/94  frac 0.149  p_bonf 3.7e-07   <- rejected on effect size alone
#     SER  contig 78  k= 9/71  frac 0.127  p_bonf 2.6e-03
# Seven orders of magnitude past alpha, refused because 0.149 < 0.20. In that regime `p_bonf` is what
# discriminates and the fraction is not informative.
# MEASURED PRICE at the SHIPPED alpha of 0.005, pooled homogeneous arms (58 carriers / 119 negatives):
# 31/58 -> **33/58 = 0.569**, and specificity 1.000 -> 0.992 — ONE false positive, which survives
# `min_cover`. (34/58 = 0.586 is the alpha 0.01 row of the same sweep and was quoted here by mistake when
# the change landed; it is not what ships.) A real trade, taken deliberately rather than discovered.
# ⚠ Re-price the moment a larger independent tandem negative set exists; the independent slice is 13.
# ⚠ And note what the all-families figure hides: all 33 true positives are dupC. The route tests the 8C
# target only, so the 13 non-dupC carriers are 0/13 — a different gap, with a different fix (family
# wiring), not something a threshold reaches.
DEFAULTS = dict(mut_len=8, wt_len=7, alpha=0.005, min_k=4, min_frac=0.10,
                min_reads=6, min_null_obs=30, min_cover=12)


def is_expansion(mut_len: int, wt_len: int) -> bool:
    """True when the target is LONGER than WT (dupC 8 vs 7) — i.e. the `>=` statistic applies. Pure."""
    return mut_len >= wt_len


def k_shifted(vals, mut_len: int, wt_len: int) -> int:
    """Reads rendered at least as extreme as the mutant, IN THE MUTANT'S OWN DIRECTION. Pure.

    An expansion and a contraction are not the same test with a different number. `>= mut_len` is right
    for a dupC (it keeps the 9's, which are over-called 8's) and MEANINGLESS for a delCC: every WT read
    (7C) satisfies `>= 5`, so `frac` pins to ~1.0 at every position. That does NOT yield p≈1 and a silent
    non-call — at amplicon depth the null sits near 0.95 (reads rendering 4C elsewhere) and 0.95**n over
    thousands of reads is overwhelming, so the test ranks positions by how FEW short alignment artifacts
    they carry and calls the winner. Measured 2026-08-02 on the 30 raw amplicons: 1 FP, 0 TP — zero
    sensitivity by construction (a carrier's 5C reads are still `>= 5`) plus a live false-positive rate."""
    return sum(x >= mut_len for x in vals) if is_expansion(mut_len, wt_len) \
        else sum(x <= mut_len for x in vals)


def null_fraction(per_position: dict, *, exclude=None, mut_len: int = 8, wt_len: int = 7,
                  floor: float = 1e-4) -> float:
    """Empirical P(read renders a WT run as far as the mutant), pooled over every position but `exclude`.

    Direction follows `is_expansion(mut_len, wt_len)`; see `k_shifted`. Pure.

    Jeffreys-smoothed, NOT a raw ratio against a constant floor. Observing zero errors in N reads does not
    mean the error rate is the floor; it means it is below roughly 1/N, and a null pinned at 1e-4 turns any
    2-read blip into overwhelming significance. Measured on a real non-carrier: f_null collapsed to the
    1e-4 floor and the sample scored p_bonf 7.2e-06 — 140x past the alpha gate — surviving only because a
    support gate rejected it. (k+0.5)/(n+1) keeps the null tied to how much evidence there actually is."""
    obs = [x for pos, vals in per_position.items() if pos != exclude for x in vals]
    if not obs:
        return floor
    k = k_shifted(obs, mut_len, wt_len)
    return max((k + 0.5) / (len(obs) + 1), floor)


def verdict(cand: dict, **kw) -> tuple:
    """(called, reason) for ONE candidate position under a gate setting. Pure.

    Factored out so a gate sweep can re-score an already-written `scored` list offline — the JSON holds
    k/n/frac/p/p_bonf per position, so trading sensitivity against specificity costs no compute at all.
    It must stay the ONLY expression of these gates: a sweep that re-implemented them would price a
    threshold the caller does not actually apply, which is worse than not measuring."""
    p = {**DEFAULTS, **kw}
    if (cand.get("n_null_obs") or 0) < p["min_null_obs"]:
        return False, (f"null rests on {cand.get('n_null_obs')} observations "
                       f"(< min_null_obs={p['min_null_obs']}) — not estimable")
    if (cand.get("n") or 0) < p["min_cover"]:
        return False, (f"candidate rests on {cand.get('n')} reads (< min_cover={p['min_cover']}, the "
                       f"caller's own 90 %-power requirement) — too thin to believe either way")
    ok = bool(cand["p_bonf"] <= p["alpha"] and cand["k_ge_mut"] >= p["min_k"]
              and cand["frac"] >= p["min_frac"])
    return ok, None if ok else "no position passed p_bonf/k/frac gates"


def call_from_runlengths(per_position: dict, **kw) -> dict:
    """per_position = {position_id: [observed C-run length per read]} at reference-WT runs. Pure.

    Returns {f_null, n_tested, candidate, called, scored, reason}. `candidate` is the most significant
    position: one-sided binomial on (k reads >= mut_len out of n), Bonferroni over the tested positions."""
    from ..dupc_power import _tail

    p = {**DEFAULTS, **kw}
    tested = sorted(pos for pos, v in per_position.items() if len(v) >= p["min_reads"])
    if not tested:
        return {"called": False, "n_tested": 0, "candidate": None, "scored": [],
                "reason": f"no position reached min_reads={p['min_reads']}"}
    scored = []
    for pos in tested:
        v = per_position[pos]
        k = k_shifted(v, p["mut_len"], p["wt_len"])
        # leave-one-out: the position under test must not calibrate the null it is tested against
        n_null_obs = sum(len(vals) for q, vals in per_position.items() if q != pos)
        f_null = null_fraction(per_position, exclude=pos, mut_len=p["mut_len"], wt_len=p["wt_len"])
        pv = _tail(len(v), k, f_null)
        scored.append({"position": pos, "k_ge_mut": k, "n": len(v), "frac": round(k / len(v), 3),
                       "f_null": round(f_null, 4), "n_null_obs": n_null_obs,
                       "p": pv, "p_bonf": min(1.0, pv * len(tested))})
    scored.sort(key=lambda d: (d["p_bonf"], -d["k_ge_mut"]))
    cand = scored[0]
    called, reason = verdict(cand, **p)
    return {"called": called, "n_tested": len(tested), "candidate": cand, "scored": scored,
            "f_null": cand["f_null"], "reason": reason,
            "direction": "expansion" if is_expansion(p["mut_len"], p["wt_len"]) else "contraction",
            "mut_len": p["mut_len"], "wt_len": p["wt_len"]}


# ── IO layer (pysam) ──────────────────────────────────────────────────────────

def measure_runlengths(bam: str, contig: str, ref: str, *, wt_len: int = 7, keep_slack: int = 3,
                       start: int = None, end: int = None, ref_cram: str = None,
                       mut_len: int = None) -> dict:
    """Per reference WT C-run, the run length each read renders there. Impure (reads the bam + fasta).

    Reference runs are found as `C{wt_len}A` on `contig`; each read's own run is walked in ITS sequence from
    the aligned position, so a read carrying an insertion is measured at its true length, not the reference's.
    Reads whose run is not captured (`< wt_len - keep_slack`) are alignment artifacts and dropped.

    ⚠ `mut_len` is accepted and DELIBERATELY does not move the floor. It was lowered for the contraction
    target on 2026-08-03, reasoning that `wt_len - keep_slack` = 4 clips the tail a `<= 5` test wants (a 5C
    read drifting to 3C). Measured the next morning on 5 amplicons: **5/5 called delCC at frac 0.89-0.96,
    including two 59dupC carriers** — a heterozygous delCC gives ~0.50, so those were artifacts. The floor
    is not only about capturing the mutant tail: in a tandem, a read whose copy number differs from the
    contig carries a large indel, and at that junction the coordinate mapping breaks and the walked run is
    truncated to 0-2 bases. Those truncations are indistinguishable from a contraction and outnumber it.
    At 4 the floor keeps 4 and 5 (mutant, and mutant with a -1 error) and discards the junk, which is what
    a 5C target needs. Any future change here must be re-measured on carriers of ANOTHER family: they are
    the control that caught this one.

    `start`/`end` scope the scan to a window — REQUIRED in genomic space (scanning a whole chromosome for
    C-runs is neither tractable nor meaningful); omit them on a small faked VNTR contig. `ref_cram` is the
    genome FASTA a CRAM needs to decode, per the repo-wide CRAM rule."""
    import re
    import pysam

    fa = pysam.FastaFile(ref)
    off = int(start or 0)
    refseq = fa.fetch(contig, start, end).upper() if (start is not None or end is not None) \
        else fa.fetch(contig).upper()
    # MUC1 is on the MINUS strand. The faked VNTR reference is reverse-complemented, so the run reads
    # `C{wt}A` there — but on a PLUS-strand genome (GRCh38/T2T) the very same run reads `TG{wt}`. Searching
    # only for the C form found ZERO positions in the genomic window and the sweep reported a vacuous
    # specificity of 1.000 on 92 subjects. Detect the orientation from the sequence instead of assuming it.
    c_hits = [m.start() + off for m in re.finditer("C{%d}A" % wt_len, refseq)]
    g_hits = [m.start() + 1 + off for m in re.finditer("TG{%d}" % wt_len, refseq)]
    run_base, runpos = ("C", c_hits) if len(c_hits) >= len(g_hits) else ("G", g_hits)
    per = collections.defaultdict(list)
    mode = "rc" if str(bam).endswith(".cram") else "rb"
    with pysam.AlignmentFile(bam, mode, reference_filename=ref_cram) as af:
        for r in af.fetch(contig, start, end):
            if r.is_unmapped or not r.query_sequence:
                continue
            r2q = {rp: qp for qp, rp in r.get_aligned_pairs() if rp is not None and qp is not None}
            s = r.query_sequence
            for idx, rs in enumerate(runpos):
                q = r2q.get(rs)
                if q is None:
                    continue
                i = q
                while i < len(s) and s[i] == run_base:
                    i += 1
                j = q
                while j > 0 and s[j - 1] == run_base:
                    j -= 1
                rl = i - j
                if rl >= wt_len - keep_slack:
                    per[idx].append(rl)
    return {"per_position": dict(per), "ref_positions": runpos, "contig": contig,
            "run_base": run_base, "n_ref_runs": len(runpos),
            "n_reads_measured": sum(len(v) for v in per.values())}


def call_bam(bam: str, contig: str, ref: str, *, start: int = None, end: int = None,
             ref_cram: str = None, **kw) -> dict:
    """measure + call, on one scaffold contig (or one genomic window). Impure."""
    p = {**DEFAULTS, **kw}
    meas = measure_runlengths(bam, contig, ref, wt_len=p["wt_len"], mut_len=p["mut_len"],
                              start=start, end=end, ref_cram=ref_cram)
    res = call_from_runlengths(meas["per_position"], **kw)
    res["contig"] = contig
    res["method"] = "runlen_shift"
    # Carried so a caller can tell "measured, found nothing" from "measured NOTHING" — the two look
    # identical in `called: False` and only the second invalidates a specificity sweep.
    res["run_base"] = meas["run_base"]
    res["n_ref_runs"] = meas["n_ref_runs"]
    res["n_reads_measured"] = meas["n_reads_measured"]
    if res.get("candidate") is not None:
        res["candidate"]["ref_pos"] = meas["ref_positions"][res["candidate"]["position"]]
    return res


# ── C-run index → ARRAY UNIT index ────────────────────────────────────────────────────────────────
# `candidate.position` is the RANK of the C-run among the scaffold's reference C-runs. It is NOT the
# repeat number: units that carry no C-run are skipped, so the two spaces drift apart. Measured on the
# One family (2026-07-29): the same inherited variant, on the same 45-repeat allele, came out at repeat
# 11 through this caller and repeat 16 through the consensus — 0.7556 vs 0.6444 on the onset axis, which
# is the deliverable. `candidate.ref_pos` already carries the scaffold coordinate, so the conversion is
# arithmetic, not a heuristic.
UNIT_BP = 60


def _revcomp(s: str) -> str:
    return (s or "").translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def flank5_len(contig_seq: str, unit: str) -> "int | None":
    """Length of the 5' flank of a `MUC1_VNTR_Nrepeats` contig = offset of the first exact unit. Pure.

    The reference is built as `flank5 + N x unit + flank3` (`build_vntr_ref`), and the flank lengths are
    not recorded anywhere the caller can read — but the unit sequence is known, so the boundary is found
    by looking rather than by being told.

    ⚠ BOTH ORIENTATIONS are searched. The shipped reference is
    `MUC1_fakedVNTR1to150revcomplKirby.fa` — reverse-complemented — while `MOTIF_A` is in gene/coding
    orientation, so a single-orientation search returns -1 on the very reference the pipeline uses. Same
    trap as the 2026-07-26 sweep that scanned the C form against a plus-strand window and reported a
    vacuous specificity. Index order is unaffected: it is counted along the CONTIG, which is the frame the
    consensus numbers its motifs in too — the two witnesses must agree, and that is what makes them agree.
    Returns None if neither orientation occurs (wrong reference)."""
    seq = (contig_seq or "").upper()
    for u in ((unit or "").upper(), _revcomp(unit or "").upper()):
        if not u:
            continue
        i = seq.find(u)
        if i >= 0:
            return i
    return array_start_by_periodicity(seq)


def array_start_by_periodicity(contig_seq: str, unit_bp: int = UNIT_BP):
    """Offset of the first position from which the sequence is `unit_bp`-periodic. Pure.

    The reference-agnostic fallback, and on the SHIPPED reference the only one that works:
    `MUC1_fakedVNTR1to150revcomplKirby.fa` repeats a DIFFERENT motif variant (`GGACACCAGG` where
    `MOTIF_A` has `GGAGAGCAGG`), in another phase, so neither `MOTIF_A` nor its reverse complement
    occurs at all (measured: find() = -1 both ways, array start
    4578, flanks 7500 bp for a 2700 bp array). Periodicity assumes nothing about WHICH unit was used —
    only that the array is a tandem repeat, which is what makes it an array."""
    seq = (contig_seq or "").upper()
    n = len(seq) - 2 * unit_bp
    if n <= 0:
        return None
    best_start = best_len = cur_start = cur_len = None
    for i in range(n):
        if seq[i:i + unit_bp] == seq[i + unit_bp:i + 2 * unit_bp]:
            cur_start = i if cur_len is None else cur_start
            cur_len = 1 if cur_len is None else cur_len + 1
            if best_len is None or cur_len > best_len:
                best_start, best_len = cur_start, cur_len
        else:
            cur_start = cur_len = None
    # LONGEST periodic stretch, not the first: a low-complexity flank is periodic too (60 identical
    # bases satisfy any period), and taking the first match would put the array start inside it. The
    # array is the longest tandem stretch in the contig by construction — that is what makes it the array.
    return best_start


def unit_index_from_ref_pos(ref_pos, flank5, unit_bp: int = UNIT_BP, offset: int = 0):
    """Scaffold coordinate of a C-run → 1-based PHYSICAL repeat number. Pure.

    `offset` puts the index in the same frame as the LENGTH it will be divided by. The arbiter's allele
    length is already physical (it applies `FLANK_TO_PHYSICAL_OFFSET`, the ~4 units the AL/AH anchors sit
    inside the array), so an onset computed as `1 - index/length` from a contig-relative index would mix
    two frames. Proved on one family 2026-07-29: mother and son share ref_pos 5254 — the same
    inherited variant — giving contig unit 12, while the mother's CONSENSUS, an independent witness of
    the same physical position, reads 16. The residual is the offset, not noise.

    Returns None when the inputs do not allow the conversion — an ABSENT repeat number is honest, a
    C-run rank published as a repeat number is not."""
    if ref_pos is None or flank5 is None or unit_bp <= 0:
        return None
    off = int(ref_pos) - int(flank5)
    return (off // unit_bp) + 1 + int(offset) if off >= 0 else None
