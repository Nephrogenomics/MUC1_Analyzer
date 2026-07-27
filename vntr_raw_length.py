#!/usr/bin/env python3
"""GROUND-TRUTH VNTR length, immune to alignment (no CIGAR, no collapse artefact).

Per read, find the UNIQUE flank anchor BELOW the tandem and the one ABOVE it directly in the
READ SEQUENCE (fuzzy, both orientations); the raw #bases between them ÷ 60 = copies. This
bypasses minimap2's tandem-collapse (which makes sv8533_length under-estimate) and the
multi-contig smear (which makes MUC1_Analyzer mis-call). See memory muc1-vntr-length-pitfall.

Anchors are T2T-CHM13 unique flanks just outside the tandem (chr1:154,328,120-154,330,561):
  aL = chr1:154,328,040-154,328,099 (below) ; aH = chr1:154,330,585-154,330,644 (above).
Between-anchor CHM13 distance = 2485 bp; minus the 44 bp non-tandem gap => /60 = copies.

  python3 vntr_raw_length.py --bam SANG.t2t.bam --chrom muc1win --start 2327000 --end 2331500
  python3 vntr_raw_length.py --bam FAB_AS.t2t.bam --chrom chr1 --start 154327000 --end 154331500
"""
import argparse, re, statistics as st, sys
import pysam

# Sequencing platform from the READ NAME — the instrument writes it, so it is the reliable signal (a
# processed BAM can carry misleadingly high base qualities: e.g. our ONT LR-PCR bams read median Q≈40).
# ONT = a UUID (8-4-4-4-12 hex); PacBio CCS/HiFi = a movie/zmw/ccs name (`m64…_…/…/ccs`).
_ONT_NAME = re.compile(r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$", re.I)
_PB_NAME = re.compile(r"^m\w+_\d+.*/\d+/ccs$|/ccs$")


def platform_from_name(name):
    """'ONT' | 'PacBio HiFi' | None from a read name (UUID → ONT, .../ccs → PacBio)."""
    if not name:
        return None
    if _ONT_NAME.match(name):
        return "ONT"
    if "/ccs" in name or _PB_NAME.match(name):
        return "PacBio HiFi"
    return None

AL = "GCTCTTGCTGGCTGGGGTTGTGGTAGCCCTGGCAGAGGTGCCGTTGTGCACCAGAGTAGAA"   # below tandem (ref+)
AH = "GCTGGCCTGGTGACTGGGACCGAGGTGACATCCTGTCCCCAGGTGGCAGCTGAACCTGAAG"   # above tandem (ref+)
GAP = 44            # non-tandem bp between the two anchors in CHM13 (2485 - 2441)
# Absolute read floor used by the (withdrawn) length-bias rescue — kept only as the reference point of
# `detectable_delta`, the honest reach statement.
_BIAS_MIN_RAW = 3

def rc(s): return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

def fuzzy_find(probe, seq, maxmm, k=11):
    """Best (min-Hamming) position of probe in seq via exact k-mer seeds + verify. None if > maxmm."""
    L = len(probe); best = None
    seen = set()
    for so in range(0, L - k + 1, k):
        seed = probe[so:so + k]
        start = 0
        while True:
            j = seq.find(seed, start)
            if j < 0: break
            start = j + 1
            p = j - so                      # implied probe start
            if p < 0 or p + L > len(seq) or p in seen: continue
            seen.add(p)
            mm = sum(a != b for a, b in zip(probe, seq[p:p + L]))
            if mm <= maxmm and (best is None or mm < best[1]):
                best = (p, mm)
    return best

def read_length(seq, maxmm):
    """Raw tandem bp between aL and aH; try both read orientations. None if not spanning."""
    for s in (seq, rc(seq)):
        l = fuzzy_find(AL, s, maxmm); h = fuzzy_find(AH, s, maxmm)
        if l and h:
            gap = h[0] - (l[0] + len(AL))    # bases between aL end and aH start
            if gap > 0:
                return gap
    return None


def copies_from_bam(bam, chrom="chr1", start=None, end=None, offset=0, maxmm=9, ref=None):
    """Per-read VNTR copies via flank-to-flank `read_length`, alignment-free. Returns [(copies, HP)]
    (HP = the whatshap tag if present, else None), each read counted once. Shared by the CLI and by
    the two-axis report so both compute copies identically."""
    mode = "rc" if str(bam).endswith(".cram") else "rb"
    kw = {"reference_filename": ref} if (mode == "rc" and ref) else {}
    data = []
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        it = af.fetch(chrom, start, end) if start is not None else af.fetch(until_eof=True)
        seen = set()
        for r in it:
            if r.is_secondary or r.is_supplementary or r.query_sequence is None:
                continue
            if r.query_name in seen:                # count each read once
                continue
            seen.add(r.query_name)
            g = read_length(r.query_sequence.upper(), maxmm)
            if g is not None:
                hp = r.get_tag("HP") if r.has_tag("HP") else None
                data.append((round((g - GAP) / 60) + offset, hp))
    return data


# ── Cassette length — AMPLICON-AGNOSTIC ───────────────────────────────────────
# AL/AH above are OUR LR-PCR amplicon's unique FLANKS; a DIFFERENT MUC1 PCR need not span them (a shorter
# amplicon whose 3' primer sits inside AL — e.g. VNTRtools: 81-93 % of reads carry AH but 0 % carry AL, so
# the flank arbiter returns nothing). But the VNTR itself ALWAYS begins with the Kirby '1' motif and ends
# with the '9' motif (every reconstruction reads 1-2-3-4-5-…-6-7-8-9), so anchoring on those INVARIANT
# cassette motifs measures the length on ANY MUC1 amplicon. These two sequences == KNOWN_REPEATS['1'] /
# ['9'] in muc1_analyzer/caller.py (kept literal here so this stays a zero-dependency standalone tool).
CASSETTE_5 = "AAGGAGACTTCGGCTACCCAGAGAAGTTCAGTGCCCAGCTCTACTGAGAAGAATGCTGTG"   # motif '1' — 5' tandem start
CASSETTE_3 = "GGCTCCACCGCCCCTCCAGTCCACAATGTCACCTCGGCCTCAGGCTCTGCATCAGGCTCA"   # motif '9' — 3' tandem end


def read_cassette_span(seq, maxmm=9):
    """(copies, flank5_bp, flank3_bp) for ONE read via the invariant cassette motifs ('1' → '9'), or None.
    `copies` = 60 bp units between the cassettes (amplicon-agnostic, per-read). `flank5_bp` / `flank3_bp` =
    the non-VNTR bp BEFORE motif '1' and AFTER motif '9' — the amplicon's PRIMER flanks, i.e. its signature
    (our native LR-PCR has long flanks both sides; a shorter/foreign PCR truncates one)."""
    for s in (seq, rc(seq)):
        a = fuzzy_find(CASSETTE_5, s, maxmm)
        b = fuzzy_find(CASSETTE_3, s, maxmm)
        if a and b and b[0] > a[0]:
            return round(((b[0] + 60) - a[0]) / 60), a[0], len(s) - (b[0] + 60)
    return None


def read_copies_cassette(seq, maxmm=9):
    """VNTR copy number of ONE read via the cassette motifs (see read_cassette_span). None if not spanning."""
    r = read_cassette_span(seq, maxmm)
    return r[0] if r else None


def amplicon_signature(bam, chrom="chr1", start=None, end=None, maxmm=9, ref=None, cap=5000):
    """Characterise the amplicon from the reads themselves — WHAT PCR is this? Over up to `cap` reads:
    the median 5'/3' primer flank (bp, from the cassette anchors) and the fraction carrying OUR LR-PCR flank
    anchors AL / AH. A read spanning AL+AH = our native design; AH-only (no AL) = a shorter/foreign PCR that
    truncates the 3' flank (e.g. VNTRtools). Returns a dict (n_spanning=0 → not amplicon-scale / no MUC1)."""
    mode = "rc" if str(bam).endswith(".cram") else "rb"
    kw = {"reference_filename": ref} if (mode == "rc" and ref) else {}
    f5, f3, copies, quals, n = [], [], [], [], 0
    have_al = have_ah = 0
    n_fwd = n_rev = 0          # STRAND BALANCE: a real variant shows on BOTH strands; an ONT homopolymer
                               # artefact often on one. If the library/PCR is strand-skewed, a variant living
                               # on the weak strand is under-powered EVEN at huge depth — so measure it here,
                               # at sample level, instead of only discovering it inside the variant caller.
    from collections import Counter
    plat_votes = Counter()
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        it = af.fetch(chrom, start, end) if start is not None else af.fetch(until_eof=True)
        seen = set()
        for r in it:
            if r.is_secondary or r.is_supplementary or r.query_sequence is None:
                continue
            if r.query_name in seen:
                continue
            seen.add(r.query_name)
            seq = r.query_sequence.upper()
            span = read_cassette_span(seq, maxmm)
            if span is None:
                continue
            n += 1
            n_rev += bool(r.is_reverse); n_fwd += not r.is_reverse
            copies.append(span[0]); f5.append(span[1]); f3.append(span[2])
            p = platform_from_name(r.query_name)
            if p:
                plat_votes[p] += 1
            if r.query_qualities is not None and len(r.query_qualities):
                quals.append(sum(r.query_qualities) / len(r.query_qualities))
            rcs = rc(seq)
            if fuzzy_find(AL, seq, maxmm) or fuzzy_find(AL, rcs, maxmm):
                have_al += 1
            if fuzzy_find(AH, seq, maxmm) or fuzzy_find(AH, rcs, maxmm):
                have_ah += 1
            if n >= cap:
                break
    if not n:
        return {"n_spanning": 0}
    al_frac, ah_frac = have_al / n, have_ah / n
    kind = ("native_LR-PCR" if al_frac > 0.5 and ah_frac > 0.5
            else "foreign/short-PCR" if ah_frac > 0.5 or al_frac > 0.5
            else "unknown")
    f5med, f3med = int(st.median(f5)), int(st.median(f3))
    flank_bp = f5med + f3med
    # PCR vs AS/WGS: a PCR has FIXED primers → the cassette-to-read-end flanks cluster tightly; adaptive
    # sampling / WGS reads start at random positions → the flanks are spread. High consistency ⇒ amplicon.
    tight5 = sum(1 for x in f5 if abs(x - f5med) <= 200) / n
    tight3 = sum(1 for x in f3 if abs(x - f3med) <= 200) / n
    is_amplicon = tight5 > 0.6 and tight3 > 0.6
    alleles = call_alleles(sorted(copies)).get("alleles", [])
    # full PCR product per allele = constant flanks + the VNTR itself (copies × 60 bp) — the primer-to-primer
    # band size (varies with allele length: a short allele is a shorter product; an unexpectedly short
    # product on a long-allele carrier is the long-allele DROPOUT signature).
    product_bp = [flank_bp + c * 60 for c in alleles]
    mbq = int(st.median(quals)) if quals else None
    # platform: the read-name vote is authoritative (see platform_from_name); base Q is a weak fallback ONLY
    # when names are uninformative — and even then conservative, since a processed BAM can inflate Q.
    if plat_votes:
        platform = plat_votes.most_common(1)[0][0]
    elif mbq is not None:
        platform = "PacBio HiFi" if mbq >= 30 else "ONT"
    else:
        platform = "unknown"
    return {"n_spanning": n, "n_fwd": n_fwd, "n_rev": n_rev,
            "strand_minor_frac": round(min(n_fwd, n_rev) / n, 3) if n else None,
            "flank5_med": f5med, "flank3_med": f3med, "amplicon_flank_bp": flank_bp,
            "has_AL_frac": round(al_frac, 2), "has_AH_frac": round(ah_frac, 2),
            "kind": kind, "is_amplicon": is_amplicon, "alleles": alleles, "product_bp": product_bp,
            "median_bq": mbq, "platform": platform}


def _alleles_by_hp(data, sep=8, min_reads=2):
    """HP-aware allele call: when the reads carry HP tags (haplotagged AS/WGS), ONE allele per haplotype is
    the trustworthy split — no chimera heuristic (which wrongly drops a real SHORT allele as 'truncation' on
    non-PCR data, e.g. a haplotagged AS with alleles 32/66 mis-called HOM ~66). Returns {alleles, counts, long_low_conf}
    or None when fewer than 2 haplotypes have support (→ caller falls back to the pooled peak-caller)."""
    per = []
    for hp in sorted({h for _, h in data if h is not None}):
        v = sorted(c for c, h in data if h == hp)
        if len(v) >= min_reads:
            per.append((int(st.median(v)), len(v)))
    if len(per) < 2:
        return None
    per.sort(key=lambda x: x[0])
    if per[-1][0] - per[0][0] < sep:                     # HP medians agree → homozygous length
        return {"alleles": [per[0][0]], "counts": [sum(c for _, c in per)], "long_low_conf": False}
    return {"alleles": [per[0][0], per[-1][0]], "counts": [per[0][1], per[-1][1]], "long_low_conf": False}


def calibrate_length_bias(samples):
    """Aggregate a per-amplicon-design length bias from SAMPLES WHOSE BOTH ALLELES ARE VISIBLE.

    `samples` = iterable of per-read copy lists (one per sample). Only clean length-hets inform the fit
    (a homozygous or single-peak sample says nothing about depletion). Returns
    {f, n_samples, per_sample:[…], spread} using the MEDIAN (robust to one bad amplicon) — the value to
    then APPLY to samples where the long allele is depleted below the peak floor. Calibrating on the
    visible cases and applying to the invisible ones is what makes this non-circular."""
    import math
    import statistics as _st
    ests = [e for e in (estimate_length_bias(c) for c in samples) if e]
    if not ests:
        return None
    fs = sorted(e["f"] for e in ests)
    # POOLED least-squares slope through the origin on log(ratio) = Δ·log f — the primary estimate.
    # A per-sample f = ratio**(1/Δ) has a variance that COLLAPSES as Δ grows (any ratio to the power 1/64
    # tends to 1), so a median over per-sample f over-weights small-Δ (noisy) samples and drags large-Δ ones
    # toward 1. The pooled fit uses every sample at its true leverage. The median is kept as a robustness
    # cross-check: a large gap between the two means the geometric per-unit model fits poorly.
    num = sum(e["delta"] * math.log(e["n_long"] / e["n_short"]) for e in ests)
    den = sum(e["delta"] ** 2 for e in ests)
    pooled = math.exp(num / den) if den else None
    # model check: does the estimate drift with Δ? (it should not, if depletion really compounds per unit)
    small = [e["f"] for e in ests if e["delta"] <= 20]
    large = [e["f"] for e in ests if e["delta"] >= 33]
    drift = (round(_st.mean(small), 4), round(_st.mean(large), 4)) if small and large else None
    return {"f": round(min(pooled, 1.0), 4) if pooled else round(_st.median(fs), 4),
            "f_pooled": round(pooled, 4) if pooled else None,
            "f_median": round(_st.median(fs), 4),
            "n_samples": len(fs), "per_sample": fs,
            "spread": (round(fs[0], 4), round(fs[-1], 4)), "drift_small_vs_large_delta": drift}


def auto_length(bam, chrom="chr1", start=None, end=None, offset=0, maxmm=9, ref=None):
    """AUTO alignment-free VNTR length on a BAM — the ONE call to get the length end-to-end. Uses the flank
    arbiter (AL→AH) when the reads carry our flanks (native LR-PCR / AS / WGS → calibrated numbers) and falls
    back to the cassette (motif 1→9) otherwise (a foreign/shorter PCR). When the reads are HAPLOTAGGED, the
    two alleles come from the HP split (robust); else from the pooled peak-caller. Returns a dict:
    {available, method, kind, n, alleles, counts, long_low_conf}. `available` is False when nothing spans."""
    sig = amplicon_signature(bam, chrom, start, end, maxmm, ref)
    native = sig.get("has_AL_frac", 0) > 0.5 and sig.get("has_AH_frac", 0) > 0.5
    if native:
        data = copies_from_bam(bam, chrom, start, end, offset, maxmm, ref)
        method = "flank AL→AH"
    else:
        data = cassette_copies_from_bam(bam, chrom, start, end, maxmm, ref)
        method = "cassette 1→9"
    extra = {"platform": sig.get("platform", "unknown"), "kind": sig.get("kind", "?"),
             "is_amplicon": sig.get("is_amplicon", False), "flank_bp": sig.get("amplicon_flank_bp"),
             "n_fwd": sig.get("n_fwd"), "n_rev": sig.get("n_rev"),
             "strand_minor_frac": sig.get("strand_minor_frac")}
    if not data:
        return {"available": False, "method": method, "n": 0, **extra}
    cops = [c for c, _ in data]
    res = _alleles_by_hp(data) or call_alleles(sorted(cops))   # HP split (haplotagged) → else pooled peaks
    # full PCR product per allele (bp) = fixed flanks + copies×60 (the gel band) — recomputed on the CALLED
    # alleles (sig.product_bp is on a capped sample); only meaningful when is_amplicon.
    product_bp = [extra["flank_bp"] + c * 60 for c in res["alleles"]] if extra["flank_bp"] else []
    return {"available": True, "method": method, "n": len(cops),
            "alleles": res["alleles"], "counts": res.get("counts", []),
            "long_low_conf": res.get("long_low_conf", False), "product_bp": product_bp, **extra}


def cassette_copies_from_bam(bam, chrom="chr1", start=None, end=None, maxmm=9, ref=None):
    """Per-read cassette copies over a region (or until_eof when start is None). Returns [(copies, HP), …]
    with each read counted once — SAME shape as copies_from_bam, so call_alleles / the two-axis report
    consume it identically."""
    mode = "rc" if str(bam).endswith(".cram") else "rb"
    kw = {"reference_filename": ref} if (mode == "rc" and ref) else {}
    data = []
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        it = af.fetch(chrom, start, end) if start is not None else af.fetch(until_eof=True)
        seen = set()
        for r in it:
            if r.is_secondary or r.is_supplementary or r.query_sequence is None:
                continue
            if r.query_name in seen:
                continue
            seen.add(r.query_name)
            c = read_copies_cassette(r.query_sequence.upper(), maxmm)
            if c is not None:
                hp = r.get_tag("HP") if r.has_tag("HP") else None
                data.append((c, hp))
    return data


def estimate_length_bias(cops, *, sep=8, min_frac=0.006, min_reads=6, smooth=2):
    """MEASURE the per-unit PCR length bias `f` on a length-HET sample, instead of assuming a constant.

    In a competitive LR-PCR the longer allele amplifies less, and the depletion compounds per extra 60 bp
    unit: n_long / n_short ≈ f**Δcopies (0 < f ≤ 1 ; f = 1 = no bias). Both alleles are equimolar in the
    genome, so the observed ratio measures f directly:  f = (n_long/n_short) ** (1/Δ).

    This is deliberately measured PER SAMPLE / PER AMPLICON DESIGN rather than hardcoded (a competitor uses
    a fixed 0.9 calibrated on its own PacBio chemistry — wrong for another amplicon/chemistry: our LR-PCR
    was optimised to LIMIT long-allele dropout, so our f should sit closer to 1). Aggregate the per-sample
    f over a cohort (median) to calibrate a design. Returns {f, delta, n_short, n_long} or None when the
    sample is not informative (not a clean het)."""
    res = call_alleles(cops, sep=sep, min_frac=min_frac, min_reads=min_reads, smooth=smooth)
    al, ct = res.get("alleles", []), res.get("counts", [])
    if len(al) != 2 or len(ct) != 2 or min(ct) <= 0:
        return None
    delta = al[1] - al[0]
    if delta <= 0:
        return None
    f = (ct[1] / ct[0]) ** (1.0 / delta)
    return {"f": round(min(f, 1.0), 4), "delta": delta, "n_short": ct[0], "n_long": ct[1]}


def detectable_delta(total, f, *, min_frac=0.006, min_reads=6):
    """How much LONGER than the dominant allele can a second allele be and still be SEEN at this depth?

    Under the two-sided model a real long allele carries ≈ n_dominant·f**Δ reads and is accepted at a
    fraction of that prediction, so the binding constraint is the ABSOLUTE read floor: the allele vanishes
    once its expected count drops below `_BIAS_MIN_RAW` → Δ ≤ log(_BIAS_MIN_RAW / n_dominant) / log(f).
    An honest statement of the assay's reach (a competitor reports the same idea as
    "longest detectable VNTR"). Returns an int, or None when f is uninformative (f ≥ 1 = no bias)."""
    import math
    if not f or f >= 1.0 or total <= 0:
        return None
    if total <= _BIAS_MIN_RAW:
        return 0
    return max(0, int(math.log(_BIAS_MIN_RAW / total) / math.log(f)))


def call_alleles(cops, sep=8, min_frac=0.006, min_reads=6, smooth=2):
    """PURE: call the 1–2 VNTR alleles from a pool of per-read copy numbers, robust to PCR chimeras.

    KEY biology of the noise (verified 2026-07-17 on our LR-PCR): chimeric/truncated amplicons are always
    SHORTER than the real alleles (both flanks kept, copies lost between) → they pile up at LOW copy numbers,
    never exceed the real alleles. So the two real alleles are the dominant peak (the short, best-amplified
    allele) plus the tallest well-separated peak of HIGHER copy — NOT the two tallest peaks (that is the
    2-means bug: it picked a chimera median). The long allele can be tiny (competitive PCR depletion scales
    with the length ratio: POV 50/95 → long allele only ~0.9 % of reads) → keep a low floor but flag it.

    Returns {alleles:[copies...], counts:[...], mode, n_peaks, long_low_conf}."""
    from collections import Counter
    cops = [c for c in cops if c is not None]
    if not cops:
        return {"alleles": [], "counts": [], "mode": None, "n_peaks": 0, "long_low_conf": False}
    h = Counter(cops)
    lo, hi = min(cops), max(cops)
    total = len(cops)
    floor = max(min_reads, round(min_frac * total))

    def sc(c):                                   # smoothed count in a ±`smooth` window
        return sum(h.get(c + d, 0) for d in range(-smooth, smooth + 1))

    w = max(2, sep // 2)
    peaks = [(c, sc(c)) for c in range(lo, hi + 1)
             if sc(c) >= floor and sc(c) >= max(sc(c + d) for d in range(-w, w + 1))]
    # non-max suppression: tallest first, drop anything within `sep`
    peaks.sort(key=lambda x: (-x[1], -x[0]))
    kept = []
    for c, s in peaks:
        if all(abs(c - k[0]) >= sep for k in kept):
            kept.append((c, s))
    if not kept:
        return {"alleles": [], "counts": [], "mode": None, "n_peaks": 0, "long_low_conf": False}

    dom = kept[0]                                # tallest peak = short (best-amplified) real allele
    highers = [p for p in kept if p[0] >= dom[0] + sep]   # real 2nd allele is LONGER (chimeras are shorter)
    chosen = [dom] + ([max(highers, key=lambda x: x[1])] if highers else [])
    chosen.sort()
    counts = [s for _, s in chosen]
    long_low = len(chosen) == 2 and chosen[1][1] < max(0.02 * total, 15)
    return {"alleles": [c for c, _ in chosen], "counts": counts, "mode": dom[0],
            "n_peaks": len(kept), "long_low_conf": long_low}

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", default=None)
    ap.add_argument("--name", default=None)
    ap.add_argument("--chrom", default="chr1")
    ap.add_argument("--start", type=int, default=None)
    ap.add_argument("--end", type=int, default=None)
    ap.add_argument("--maxmm", type=int, default=9, help="max mismatches per 60bp anchor (~15%)")
    ap.add_argument("--sep", type=int, default=8, help="min copies between alleles to call het (pooled)")
    ap.add_argument("--min-frac", type=float, default=0.006, help="peak floor as a fraction of spanning reads")
    ap.add_argument("--min-reads", type=int, default=6, help="absolute peak floor (reads)")
    ap.add_argument("--offset", type=int, default=0,
                    help="copies added to every FLANK measurement. The AL/AH anchors sit ~4 units INSIDE the "
                         "array, so the flank scheme reads 4 short: use 4 to report the PHYSICAL/clinical count "
                         "(the cassette scheme already does, and needs no offset). Default 0 here = raw flank.")
    ap.add_argument("--cassette", action="store_true",
                    help="FORCE the amplicon-agnostic cassette length (motifs 1→9). Default is AUTO: use our "
                         "flank anchors when the reads span them (native LR-PCR, calibrated numbers) and fall "
                         "back to the cassette otherwise (a foreign/shorter PCR). --offset is ignored here.")
    ap.add_argument("--flank", action="store_true",
                    help="FORCE the flank-to-flank arbiter (AL+AH) even if the reads look like a foreign PCR.")
    ap.add_argument("--calibrate-bias", nargs="+", metavar="BAM", default=None,
                    help="measure the length bias on these BAMs (length-HET samples only inform it) → the "
                         "median f to feed back via --length-bias")
    ap.add_argument("--amplicon", action="store_true",
                    help="report the detected amplicon signature (5'/3' primer flank bp, which of our flank "
                         "anchors AL/AH the reads carry → native LR-PCR vs foreign/short PCR) and exit.")
    a = ap.parse_args(argv)
    name = a.name or a.bam.split("/")[-1]
    if a.calibrate_bias:
        print("\n=== PCR length-bias calibration (per-design, measured — not imported) ===")
        per, samples = [], []
        for b in a.calibrate_bias:
            d = (copies_from_bam(b, a.chrom, a.start, a.end, a.offset, a.maxmm, a.ref) if not a.cassette
                 else cassette_copies_from_bam(b, a.chrom, a.start, a.end, a.maxmm, a.ref))
            cops = [c for c, _ in d]
            samples.append(cops)
            e = estimate_length_bias(cops)
            tag = (f"f={e['f']:.4f}  (Δ={e['delta']}, n {e['n_short']}/{e['n_long']})" if e
                   else "not a clean het → uninformative")
            print(f"  {b.split('/')[-1]:<32} {tag}")
            if e:
                per.append(e)
        cal = calibrate_length_bias(samples)
        if not cal:
            print("  → no informative sample (need length-HET amplicons)"); return 0
        print(f"\n  POOLED f = {cal['f']}  over {cal['n_samples']} het sample(s)   "
              f"[median cross-check {cal['f_median']}, spread {cal['spread']}]")
        if cal.get("drift_small_vs_large_delta"):
            s_, l_ = cal["drift_small_vs_large_delta"]
            print(f"  model check: mean f small-Δ {s_} vs large-Δ {l_} (a big gap ⇒ per-unit model fits poorly)")
        print(f"  → apply with:  --length-bias {cal['f']}")
        print("  (f close to 1 = little long-allele depletion; a competitor hardcodes 0.9 from its own "
              "PacBio chemistry — measure yours instead.)")
        return 0
    if a.amplicon:
        sig = amplicon_signature(a.bam, a.chrom, a.start, a.end, a.maxmm, a.ref)
        print(f"\n=== {name} === amplicon signature")
        if not sig.get("n_spanning"):
            print("  no cassette-spanning reads — not an MUC1 amplicon, or too shallow"); return 0
        print(f"  cassette-spanning reads : {sig['n_spanning']}")
        print(f"  primer flanks (non-VNTR): 5'={sig['flank5_med']} bp  3'={sig['flank3_med']} bp  "
              f"(≈ {sig['amplicon_flank_bp']} bp)")
        if sig.get("product_bp"):
            prod = "  ".join(f"~{p} bp ({c} cp)" for p, c in zip(sig['product_bp'], sig['alleles']))
            print(f"  full PCR product        : {prod}   [= flanks + VNTR, i.e. the gel band]")
        print(f"  our flank anchors       : AL in {sig['has_AL_frac']*100:.0f}% of reads, "
              f"AH in {sig['has_AH_frac']*100:.0f}%")
        bqs = f" (median base Q {sig['median_bq']})" if sig.get("median_bq") is not None else ""
        print(f"  platform (heuristic)    : {sig['platform']}{bqs}")
        print(f"  → amplicon kind         : {sig['kind']}")
        return 0
    def _flank():
        return copies_from_bam(a.bam, a.chrom, a.start, a.end, a.offset, a.maxmm, a.ref)

    def _cassette():
        return cassette_copies_from_bam(a.bam, a.chrom, a.start, a.end, a.maxmm, a.ref)

    if a.cassette:
        data, method = _cassette(), "cassette 1→9 (forced)"
    elif a.flank:
        data, method = _flank(), "flank-to-flank (forced)"
    else:
        # AUTO: our flank anchors when the reads carry them (native LR-PCR → keep the calibrated numbers,
        # no regression), else the amplicon-agnostic cassette (a foreign/shorter PCR). Zero-config.
        sig = amplicon_signature(a.bam, a.chrom, a.start, a.end, a.maxmm, a.ref)
        native = sig.get("has_AL_frac", 0) > 0.5 and sig.get("has_AH_frac", 0) > 0.5
        if native:
            data, method = _flank(), "flank-to-flank (auto: native LR-PCR)"
        else:
            data, method = _cassette(), f"cassette 1→9 (auto: {sig.get('kind', 'flanks absent')})"
    anchor_miss = "no reads span the length anchors (wrong locus / coordinates / too shallow)"
    print(f"\n=== {name} === RAW VNTR length ({method}, alignment-free)")
    if not data:
        print(f"  {anchor_miss}"); return 0
    cops = sorted(c for c, _ in data)
    print(f"  spanning reads: {len(cops)} | copies: min={cops[0]} med={int(st.median(cops))} max={cops[-1]}")
    hps = {h for _, h in data if h is not None}
    if hps:
        for hp in sorted(hps):
            v = sorted(c for c, h in data if h == hp)
            print(f"  HP{hp}: n={len(v):3} median={int(st.median(v))} copies (range {v[0]}-{v[-1]})")
    else:
        # peak-based caller (robust to PCR chimeras: the real alleles are the dominant peak + the tallest
        # HIGHER-copy peak; everything shorter is chimeric truncation). Replaces the old 2-means, which
        # mistook a chimera median for an allele (e.g. WOO 26&66 instead of 66&76).
        res = call_alleles(cops, sep=a.sep, min_frac=a.min_frac, min_reads=a.min_reads)
        al, ct = res["alleles"], res["counts"]
        if len(al) == 2:
            flag = "  [LONG allele near the PCR-chimera floor — low confidence]" if res["long_low_conf"] else ""
            print(f"  -> HET: {al[0]} (n={ct[0]}) & {al[1]} (n={ct[1]}) copies{flag}")
        elif len(al) == 1:
            print(f"  -> HOM: ~{al[0]} (n={ct[0]}) copies")
            # possible long-allele PCR DROPOUT: reads clearly longer than the single called allele that sit
            # below the peak floor. A shorter PCR that under-amplifies the long allele looks homozygous —
            # this shoulder is the tell (our LR-PCR was optimised to avoid it; a foreign/short PCR may not).
            longer = [c for c in cops if c >= al[0] + a.sep]
            if longer:
                floor = max(a.min_reads, round(a.min_frac * len(cops)))
                print(f"     ⚠ {len(longer)} read(s) longer than {al[0]} (up to {max(longer)} copies) below the "
                      f"peak floor ({floor}) — possible LONG-ALLELE DROPOUT / low-coverage mosaic, not a clean HOM")
        else:
            print(f"  -> no clear peak (n_peaks={res['n_peaks']})")
    # histogram (5-copy bins)
    from collections import Counter
    b = Counter((c // 5) * 5 for c in cops)
    print("  histogram (copies):")
    for k in sorted(b):
        print(f"    {k:3d}-{k+4:3d} | {'#'*b[k]} {b[k]}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
