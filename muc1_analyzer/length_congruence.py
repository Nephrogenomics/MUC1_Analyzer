#!/usr/bin/env python3
"""length_congruence — cross-check two INDEPENDENT VNTR length measurements.

Two methods now measure the same array by different physics:

  · the **arbiter** (`vntr_raw_length.auto_length`) — alignment-free, counts units between the AL/AH flanks
    (or the 1→9 cassette when the flanks are absent). This one is AUTHORITATIVE: it is the scale validated
    against 7 clinical truths, and every number in the report is on it.
  · the **burden sweep** (the `MUC1_Analyzer_fromfastq` idea) — alignment-based, scans the faked contigs
    around an estimate and keeps the one the reads fit with the lowest indel burden.

Agreement is evidence; disagreement is information. What this module refuses to do is average them.

## Why the tiers are not a single tolerance

A CONSTANT offset and a DISPERSION are different failures and must not share a threshold. The flank→physical
conversion is exactly **+4** (`config.FLANK_TO_PHYSICAL_OFFSET`): the AL/AH anchors sit ~4 units inside what
the cassette counts. So a discrepancy of about 4 is not a noisy measurement — it is the signature of that
offset being applied on one side and not the other, or applied twice. It is the single most likely bug in
the whole length stack, and a tolerance of "flag only beyond 4" would swallow it silently.

Hence a dedicated tier at |Δ| ≈ 4, louder than the tier above it.
"""
from __future__ import annotations

# |Δ| in copies, per allele
TOL_CONCORDANT = 1          # measurement noise on a tandem: the two agree
TOL_SPREAD = 3              # up to here, imprecision — report both, no alarm
OFFSET_LO, OFFSET_HI = 3.5, 4.5     # the offset-shaped window; see the module docstring

VERDICTS = ("CONCORDANT", "SPREAD", "OFFSET_SUSPECT", "DIVERGENT", "INSUFFICIENT")


MIN_READS_TO_ROUTE = 30     # below this the flank-tightness fractions are noise, not a signature


# What the user can DECLARE about the input. An explicit declaration always wins over any heuristic: the
# person who ran the assay knows what it is, and no read statistic is as reliable as that.
INPUT_TYPES = {
    "pcr-ont":     {"route": "pcr", "pacbio": False},
    "pcr-pacbio":  {"route": "pcr", "pacbio": True},
    # Adaptive sampling is NANOPORE-only (real-time rejection of a molecule); there is no PacBio
    # equivalent. PacBio enrichment is capture-based — 'targeted'. PacBio amplicons are `pcr-pacbio`.
    "as-ont":            {"route": "as", "pacbio": False},
    "targeted-pacbio":   {"route": "as", "pacbio": True},
    "wgs-ont":           {"route": "as", "pacbio": False},
    "wgs-pacbio":        {"route": "as", "pacbio": True},
}


def route_from_declared(input_type: str):
    """(route, is_pacbio) for a declared input type, or (None, None) for 'auto'/unknown. Pure.

    Every non-amplicon type routes the same way on purpose: adaptive sampling, PacBio capture and WGS
    share the one property the fork cares about — no primers, so the MUC1 reads have to be FOUND by a
    streaming search over the whole file. How they were enriched changes nothing downstream. What does
    change is the aligner preset, which the declaration also sets."""
    d = INPUT_TYPES.get(str(input_type or "").lower())
    return (d["route"], d["pacbio"]) if d else (None, None)


def route_from_signature(n: int, al_frac: float, ah_frac: float,
                         tight5: float, tight3: float, n_scanned: int = None,
                         len_tight: float = None, median_len: int = None) -> str:
    """Same fork, decided from the READ SEQUENCES alone — no alignment. Pure.

    This is what lets `run` fork BEFORE `prepare`: aligning tens of millions of adaptive-sampling reads
    just to discover they should have gone to the other pipeline is the one thing the fork must avoid.

    The measurements are `amplicon_signature`'s, computed on a capped sample of reads:
    `tight5`/`tight3` = the fraction whose cassette-to-read-end flank sits within 200 bp of the median. A
    PCR has FIXED primers so those cluster; adaptive sampling starts at random positions so they spread.
    `al_frac`/`ah_frac` = the fraction carrying our LR-PCR anchors, which decides flank vs cassette.

    Below `MIN_READS_TO_ROUTE` MUC1-like reads the anchor fractions are noise, and on a real 20 GB
    adaptive-sampling uBAM the cassette scan found ZERO in 400 000 reads — the normal outcome there, not
    a failure. So the fallback answers the question that does NOT need a MUC1 read.

    When no MUC1 read is identifiable this returns 'unknown' and the caller must ASK: read statistics
    cannot answer it. Two candidates were tried on real data and both failed.

    · flank TIGHTNESS: 63 % measured on a real LR-PCR against 69 % on an AS uBAM — the AS file is tighter
      than the amplicon, because a PCR on a length-heterozygote makes two product sizes plus truncated reads.
    · MEDIAN read length: 578 bp on that AS uBAM looked decisive against 6727 bp for the PCR, but 578 is a
      property of THAT FILE (it kept the rejected fragments), not of adaptive sampling. A filtered AS run,
      or any WGS — ONT or PacBio — has a median in the tens of kb and would be misrouted by that rule.

    So an inconclusive signature is reported as such, and the user declares the input with
    `--input-type` (`route_from_declared`). The person who ran the assay knows what it is."""
    if not n or n < MIN_READS_TO_ROUTE:
        return "unknown"
    if tight5 > 0.6 and tight3 > 0.6:
        return "pcr"
    return "as" if (al_frac > 0.5 and ah_frac > 0.5) else "cassette"


def route_for(arbiter) -> str:
    """Which pipeline owns this sample: 'pcr' | 'as' | 'cassette' | 'unknown'. Pure.

    The fork is `is_amplicon` x `our flanks present`, not simply AS vs PCR:

    · **pcr** — our PCR path owns it. It carries the frozen chaining flags (`-z 600,200 -r 2000,20000`),
      without which a read spanning a long tandem does not chain and is silently misclassified, and the
      cassette fallback. A ruler built for AS has neither.
    · **as** — adaptive sampling / WGS carrying our AL/AH flanks: the terrain the burden sweep was written
      for, so the length gets a second, independent measurement.
    · **cassette** — no AL/AH flanks (a foreign or shorter amplicon). Only the 1→9 cassette can measure it;
      a flank-anchored sweep would return nothing at all, so there is nothing to cross-check against.
    """
    if not arbiter or not arbiter.get("available"):
        return "unknown"
    if arbiter.get("is_amplicon"):
        return "pcr"
    return "as" if "flank" in str(arbiter.get("method", "")) else "cassette"


def cross_check_applies(arbiter) -> bool:
    """True when a second length method can meaningfully run. Pure."""
    return route_for(arbiter) == "as" and bool(arbiter.get("alleles"))


def _pair(a, b):
    """Pair the two allele lists by sorted length. Pure. Returns (pairs, ploidy_note|None)."""
    a = sorted(int(x) for x in (a or []))
    b = sorted(int(x) for x in (b or []))
    note = None
    if a and b and len(a) != len(b):
        note = (f"ploidy disagreement: the arbiter called {len(a)} allele(s), the sweep {len(b)} — "
                f"comparing the {min(len(a), len(b))} that can be paired")
    n = min(len(a), len(b))
    return list(zip(a[:n], b[:n])), note


def compare(arbiter_alleles, sweep_alleles, *, arbiter_name="arbiter (alignment-free)",
            sweep_name="burden sweep (alignment-based)") -> dict:
    """Compare two allele-length calls. Pure.

    Returns {verdict, max_delta, pairs, authoritative, alternative, notes, message}. `authoritative` is
    ALWAYS the arbiter's call — the sweep is a control, never a vote. When they disagree the report is
    expected to show both, which is what `alternative` is for."""
    pairs, ploidy = _pair(arbiter_alleles, sweep_alleles)
    notes = [ploidy] if ploidy else []
    if not pairs:
        return {"verdict": "INSUFFICIENT", "max_delta": None, "pairs": [],
                "authoritative": list(arbiter_alleles or []), "alternative": list(sweep_alleles or []),
                "notes": notes + ["only one method produced a length — nothing to cross-check"],
                "message": ""}

    deltas = [abs(x - y) for x, y in pairs]
    d = max(deltas)
    if d <= TOL_CONCORDANT:
        verdict = "CONCORDANT"
        msg = (f"length CONCORDANT across two independent methods "
               f"({'/'.join(str(x) for x, _ in pairs)} copies, max Δ {d})")
    elif OFFSET_LO <= d <= OFFSET_HI:
        # checked BEFORE the spread tier: 4 falls inside "> 3" too, and the offset reading is the useful one
        verdict = "OFFSET_SUSPECT"
        msg = (f"⚠ length discrepancy of {d} copies — THE SIZE OF THE FLANK→PHYSICAL OFFSET. This is the "
               f"signature of a conversion error (offset applied on one side only, or twice), not of "
               f"imprecision. Do not average: check which measurement carries the +4.")
    elif d <= TOL_SPREAD:
        verdict = "SPREAD"
        msg = (f"length methods differ by {d} copies — measurement spread; both reported, "
               f"the {arbiter_name} value stands")
    else:
        verdict = "DIVERGENT"
        msg = (f"⚠ length DIVERGENT between methods (max Δ {d} copies): "
               f"{arbiter_name} {[x for x, _ in pairs]} vs {sweep_name} {[y for _, y in pairs]} — "
               f"the length is not reliable on this sample")
    return {"verdict": verdict, "max_delta": d,
            "pairs": [{"arbiter": x, "sweep": y, "delta": abs(x - y)} for x, y in pairs],
            "authoritative": list(arbiter_alleles or []), "alternative": list(sweep_alleles or []),
            "notes": notes, "message": msg}


def compare_caller(arbiter_alleles, caller_haplotypes) -> dict:
    """Arbiter lengths vs the CALLER's haplotype lengths — the pair the score actually rides on. Pure.

    This is not a variant of `compare`: the sweep is a research control, whereas the caller's haplotype
    lengths are what `clinical_call` turns into `mut_len` / `healthy_len` / `onset_index`. A disagreement
    here is a wrong clinical number, not a methodological curiosity.

    Observed on a real carrier: the arbiter measured 45/77 (979 and 96 reads) from the FASTQ AND from the
    chr1 BAM, identically, while the caller returned 43 and 44 — it missed the PCR-depleted long allele
    entirely and split the dominant one into two neighbours. The score then reported
    `mut=44 / healthy=43`, onset_index 0.6364, on two alleles that do not exist.

    `caller_haplotypes` = the analyzer JSON's haplotype list (dicts with a `length`, or plain ints)."""
    def _len(h):
        return h.get("length") if isinstance(h, dict) else h
    caller = sorted(int(_len(h)) for h in (caller_haplotypes or []) if _len(h))
    arb = sorted(int(x) for x in (arbiter_alleles or []))
    if not arb or not caller:
        return {"verdict": "INSUFFICIENT", "arbiter": arb, "caller": caller, "message": "",
                "missing_allele": None}
    pairs, _ = _pair(arb, caller)
    deltas = [abs(x - y) for x, y in pairs]
    d = max(deltas, default=0)
    # ORDER MATTERS. A SYSTEMATIC shift of the offset size on every allele is a conversion error, and it
    # must be recognised before the missing-allele test — otherwise "45/77 vs 41/73" reads as two lost
    # alleles when it is one arithmetic mistake applied twice.
    systematic = (len(pairs) == len(arb) == len(caller) and deltas
                  and all(OFFSET_LO <= x <= OFFSET_HI for x in deltas))
    # An allele the caller never produced is the failure that matters: it is the one the variant may sit
    # on. The tolerance here is deliberately generous — this asks "did it find this allele at all", not
    # "did it size it exactly".
    missing = [] if systematic else [a for a in arb
                                     if not any(abs(a - c) <= OFFSET_HI for c in caller)]
    if systematic:
        verdict = "OFFSET_SUSPECT"
        msg = (f"⚠ caller lengths are shifted from the arbiter by ~{d} copies on EVERY allele — the size "
               f"of the flank→physical offset ({arb} vs {caller}): a conversion error, not imprecision.")
    elif missing:
        verdict = "ALLELE_MISSING"
        msg = (f"⚠ the caller did NOT recover the {missing} allele(s) the arbiter measured "
               f"(arbiter {arb} vs caller {caller}). The score's mut/healthy lengths and onset_index are "
               f"computed on the caller's haplotypes, so they are WRONG on this sample.")
    elif d <= TOL_CONCORDANT:
        verdict, msg = "CONCORDANT", ""
    elif d <= TOL_SPREAD:
        verdict = "SPREAD"
        msg = (f"caller lengths {caller} differ from the arbiter {arb} by up to {d} copies — measurement "
               f"spread on a tandem; the arbiter value stands")
    else:
        verdict = "DIVERGENT"
        msg = (f"⚠ caller lengths {caller} disagree with the arbiter {arb} (max Δ {d}); the arbiter is "
               f"authoritative, so the score's lengths should not be read as clinical.")
    return {"verdict": verdict, "arbiter": arb, "caller": caller, "max_delta": d,
            "missing_allele": missing or None, "message": msg}


def render(cong) -> str:
    """One report line, or '' when there was nothing to compare. Pure."""
    if not cong or cong["verdict"] == "INSUFFICIENT":
        return ""
    L = [cong["message"]]
    if cong["verdict"] != "CONCORDANT":
        L.append("    " + " · ".join(f"allele {p['arbiter']} vs {p['sweep']} (Δ{p['delta']})"
                                     for p in cong["pairs"]))
    L += ["    note: " + n for n in cong.get("notes", [])]
    return "\n".join(L) + "\n"


def show_both(cong) -> bool:
    """Should the report print BOTH measurements? Pure. Concordant → one number; otherwise both."""
    return bool(cong) and cong.get("verdict") in ("SPREAD", "OFFSET_SUSPECT", "DIVERGENT")
