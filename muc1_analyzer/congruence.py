#!/usr/bin/env python3
"""congruence — cross-check the INDEPENDENT layers that call the pathogenic MUC1 variant.

Why
---
`run` stacks three layers that answer overlapping questions about the SAME variant:
  · `dupc_dispatch` → `pcr_dupc` / statistical caller : is there a dupC, and on which allele BIN?
  · `pcr_report`                                      : which motif, at which array REPEAT, on which allele?
  · the haplotype caller's nomenclature               : the consensus decomposition of each allele.
Until now they were printed side by side and NEVER compared, so a disagreement was silent. On a real case
the three appeared to contradict each other on the carrier allele and it took a manual computation to see
which one was wrong — no clinician would do that. This module makes the comparison explicit: agreement is
stated, and a DIVERGENCE is surfaced rather than hidden.

The decisive calculation
------------------------
`infer_carrier_allele` reproduces the reasoning that resolved that case, so it happens automatically:
a heterozygous variant is carried by ONE allele, so the fraction of reads showing it should match THAT
allele's abundance. With carrier reads 74/952 = 7.8 % and alleles 34 cp (93.1 % of reads) / 70 cp (6.9 %),
the variant can only sit on the 70-copy allele. This is independent of any read BINNING — which matters,
because `pcr_dupc` bins by C-tract count and that does not track array length (it mis-split that sample
366/541 when the true split was 351/26, diluting the signal).

⚠ It is an INFERENCE, not proof: it assumes a heterozygous variant and detection proportional to allele
abundance. It is reported as such, with the observed numbers, never as a hard call.

Pure logic (no pysam / no I/O) → unit-tested.
"""
from __future__ import annotations

#: minimum share of a variant's reads that must decode EXACTLY for it to be considered a real call at all.
#: Below this the motif was only reached by approximate matching — the signature of a decoding artefact
#: (e.g. the 58_59delCC cascade: 122 reads but 6 % exact).
MIN_EXACT_FRAC = 0.5


def infer_carrier_allele(carrier_reads, total_reads, alleles, counts, *, tol_factor=2.5):
    """Which allele carries the variant, from the CARRIER READ FRACTION vs each allele's abundance.

    A het variant sits on one allele, so `carrier_reads / total_reads` should track that allele's share of
    the reads. Returns {carrier_frac, abundances, inferred, ratio, confident} — `inferred` is the allele
    whose abundance is closest (in log-ratio) to the observed carrier fraction, `confident` is False when
    the two alleles are too close to tell apart or the numbers are too thin to mean anything.
    """
    if not alleles or not counts or len(alleles) != len(counts) or not total_reads or carrier_reads is None:
        return None
    tot = sum(counts)
    if tot <= 0 or carrier_reads <= 0:
        return None
    frac = carrier_reads / total_reads
    ab = [c / tot for c in counts]
    # compare in RATIO space: being 2x off matters the same whether the abundance is 5 % or 50 %
    dev = [max(frac, a) / min(frac, a) if min(frac, a) > 0 else float("inf") for a in ab]
    best = min(range(len(ab)), key=lambda i: dev[i])
    others = [dev[i] for i in range(len(ab)) if i != best]
    # confident when the best allele is within tolerance AND clearly better than the alternative
    confident = dev[best] <= tol_factor and bool(others) and min(others) > dev[best] * 1.5
    return {"carrier_frac": round(frac, 4),
            "abundances": [round(a, 4) for a in ab],
            "inferred": alleles[best], "ratio_to_inferred": round(dev[best], 2),
            "confident": confident}


def _norm_variant(v):
    """'X:59dupC' / 'X-59dupC' / '59dupC' → '59dupC' so layers using different renderings compare equal."""
    if not v:
        return None
    s = str(v).strip()
    for sep in (":", "-"):
        if sep in s:
            s = s.split(sep, 1)[1] if len(s.split(sep, 1)) > 1 else s
    return s or None


def best_reported_variant(report):
    """The credible variant of a `pcr_report` result, or None. Pure.

    Credible = most EXACTLY-decoded reads, among those decoded credibly at all. Ranking by exact FRACTION
    alone was wrong: a small perfect call outranks a large real one (observed on 4 cohort samples, where a
    22/22 X:60dupA displaced the true 59dupC). Ranking by RAW read count is equally wrong (a 122-read
    6 %-exact decoding artefact would win). n_exact combines both: support AND decode quality; the
    fraction floor only discards artefacts.

    Extracted so the SCORE reads the same variant the congruence block does — one ranking rule, one place.
    Re-deriving it elsewhere is what once printed `X:60dupA` next to a CONCORDANT verdict."""
    variants = (report or {}).get("variants") or []
    credible = [v for v in variants
                if (v.get("n_exact") or 0) / max(v.get("n_reads") or 1, 1) >= MIN_EXACT_FRAC]
    return max(credible, key=lambda v: v.get("n_exact") or 0) if credible else None


def check(*, dupc=None, report=None, arbiter=None, hap_variant=None, hap_repeat=None):
    """Compare the layers on VARIANT / REPEAT / CARRIER ALLELE. Pure.

    `dupc`   = a dupc_dispatch result dict ; `report` = a pcr_report result dict ;
    `arbiter`= an auto_length dict (alleles/counts) ; `hap_*` = the haplotype caller's nomenclature call.
    Returns {verdict, checks:[…], allele_inference, notes}. `verdict` is CONCORDANT when every axis that
    could be compared agrees, DIVERGENT when at least one disagrees, INSUFFICIENT when nothing could be.
    """
    res = (dupc or {}).get("result", dupc) or {}
    best = best_reported_variant(report)

    checks, notes = [], []

    # ── variant identity ──
    d_pos = res.get("status") == "positive" or bool(res.get("called"))
    r_var = _norm_variant(best.get("kmoch") or best.get("motif")) if best else None
    h_var = _norm_variant(hap_variant)
    seen = {k: v for k, v in (("pcr_report", r_var), ("haplotype_caller", h_var)) if v}
    # A dupC NAMED by the consensus ranking while the DEDICATED dupC caller returned a POWERED negative is
    # a contradiction, not a silence. Measured on the truth-scored cohort: four non-dupC carriers (delCC,
    # dupG, insG, del8_27) each came back `X:59dupC` from the ranking with the dupC caller negative, and
    # this block called them CONCORDANT — because it only compared the sources that spoke, and read the
    # caller's negative as "did not speak". ONT manufactures C-run tokens; the caller is the authority on
    # dupC, and an underpowered/indeterminate verdict still counts as silence.
    d_neg = (res.get("assessed", True) and not d_pos
             and str(res.get("status") or "").lower() in ("negative", "not called", "false"))
    if d_pos or seen:
        agree = (len(set(seen.values())) <= 1) and (not d_pos or not seen or "dupC" in (r_var or h_var or ""))
        # EITHER naming source contradicts a powered negative, not just `pcr_report`. The guard below was
        # written for the consensus ranking and read only `r_var`, so the identical contradiction coming
        # from the haplotype caller printed CONCORDANT. (Not what happened on one adaptive-sampling run: there the
        # dupC caller had 6 reads, so its negative is silence and the single speaking layer is no
        # contradiction — which is why this needs `d_neg`, i.e. a negative WITH power, not any negative.)
        contradicted = bool(d_neg and any("dupc" in (v or "").lower() for v in seen.values()))
        if contradicted:
            agree = False
            notes.append("the consensus ranking names a dupC that the dedicated dupC caller REJECTED "
                         "with adequate power — on this terrain the ranking is not a detector (ONT "
                         "manufactures C-run tokens); trust the caller, and ask `pcr_variant` for the "
                         "variant family actually suspected")
        checks.append({"axis": "variant", "agree": bool(agree),
                       "sources": {"dupc_caller": ("positive" if d_pos else
                                                   "NEGATIVE (powered)" if d_neg else "not called"),
                                   **seen}})

    # ── repeat position ──
    # A repeat number the reads do not agree on cannot arbitrate anything: on truncated reads the reported
    # index is the most common STARTING OFFSET, so comparing it would manufacture both false agreements and
    # false divergences. Report it, refuse to vote on it.
    rep_unreliable = bool(best) and best.get("repeat_confident") is False
    reps = {k: v for k, v in (("pcr_report", best.get("repeat_index") if best else None),
                              ("haplotype_caller", hap_repeat)) if v is not None}
    if rep_unreliable and reps:
        checks.append({"axis": "repeat", "agree": None, "sources": reps,
                       "note": "position not comparable — carrier reads are truncated"})
        notes.append(f"repeat position from `pcr_report` is UNRELIABLE ("
                     f"{best.get('repeat_concentration', 0):.0%} of carrier reads agree, spread over "
                     f"{best.get('repeat_index_distinct', 0)} values): the reads do not span enough of the "
                     f"array to place the unit. The VARIANT stands; only its position is unresolved.")
    elif len(reps) >= 2:
        vals = list(reps.values())
        checks.append({"axis": "repeat", "agree": max(vals) - min(vals) <= 1, "sources": reps})
    elif reps:
        notes.append(f"repeat position available from one layer only ({list(reps)[0]})")

    # ── carrier allele ──
    inf = None
    if best and arbiter and arbiter.get("alleles"):
        inf = infer_carrier_allele(best.get("n_reads"), (report or {}).get("n_reads"),
                                   arbiter["alleles"], arbiter.get("counts") or [])
    allele_sources = {}
    if res.get("mut_allele"):
        allele_sources["dupc_caller(bin)"] = res["mut_allele"]
    if best and best.get("allele_units") is not None:
        allele_sources["pcr_report(decomposition)"] = best["allele_units"]
    if inf and inf.get("inferred") is not None:
        allele_sources["read-fraction inference"] = inf["inferred"]
    if allele_sources:
        # only the read-fraction inference is trustworthy on a depleted allele; the other two are known-weak
        # (C-tract binning does not track length; truncated reads decode short) → report, do not vote.
        checks.append({"axis": "carrier_allele", "agree": None, "sources": allele_sources,
                       "note": "sources use different definitions — trust the read-fraction inference"})
        notes.append("carrier-allele: `dupc_caller` reports its own C-tract BIN (not the array length) and "
                     "`pcr_report` the decomposition length of carrier reads (short-biased when reads are "
                     "truncated); the read-fraction inference is the one comparable to the true alleles.")

    comparable = [c for c in checks if c.get("agree") is not None]
    verdict = ("INSUFFICIENT" if not comparable
               else "CONCORDANT" if all(c["agree"] for c in comparable) else "DIVERGENT")
    action = None
    if any(c["axis"] == "variant" and c.get("agree") is False for c in checks):
        action = review_action(best, dupc_rejected=d_neg)
        notes.append(action["message"])
    return {"verdict": verdict, "checks": checks, "allele_inference": inf, "notes": notes,
            "best_variant": best, "action": action}


# ONT's failure mode is the HOMOPOLYMER: a C-run rendered one base long or short. So a contested C-tract
# token is artefact-prone, while an inserted G is a distinct base the basecaller does not invent — the
# reason `pcr_variant` calls dupG/insG at least as well as dupC (docs/clef_generique_nonDupC.md).
_C_TRACT_FAMILIES = ("dupc", "delcc", "delc", "dupcccc")


def artefact_prior(variant, *, exact_frac=None, dupc_rejected=False) -> tuple:
    """(prior, why) for a CONTESTED variant token — 'likely' | 'unlikely' | 'unknown'. Pure.

    Not a verdict, a triage: it decides whether the divergence is worth a human's time BEFORE sending
    anyone to the images. Three signals, all already measured:
      · the family — a C-tract token is what ONT manufactures, an inserted G is not;
      · the decode quality — `n_exact / n_reads` of the token;
      · the dedicated dupC caller having REJECTED it with power, which outranks a consensus token.
    """
    fam = str(_norm_variant(variant) or "").lower()
    if not fam:                       # no token at all — triaging it either way would be inventing
        return "unknown", "no variant token to assess"
    c_tract = any(f in fam for f in _C_TRACT_FAMILIES)
    if c_tract and dupc_rejected:
        return "likely", ("a C-tract token that the dedicated dupC caller rejected with power — the ONT "
                          "homopolymer artefact looks exactly like this")
    if exact_frac is not None and exact_frac < 0.7:
        return "likely", f"only {exact_frac:.0%} of the supporting reads decode it exactly"
    if not c_tract and (exact_frac is None or exact_frac >= 0.9):
        return "unlikely", ("an inserted/deleted G is a distinct base, not a homopolymer slip — the "
                            "basecaller does not manufacture it")
    return "unknown", "no signal strong enough to triage this token either way"


def review_action(best, *, dupc_rejected=False) -> dict:
    """What to DO about a contested variant type. Pure.

    The point of the program is to SIGNAL a frameshift, not to name one it cannot type. So when the layers
    disagree on WHICH variant it is: state that a frameshift is present, triage the artefact risk, and send
    the reader to the visual bundle — which exists for exactly this (`review_bundle`: per-HP BAM on a single
    contig, ad hoc reference, consensus, units, IGV session, ALARM line, and a `verdict` column to fill).
    Typing it by eye is the documented route for the variants no dedicated caller covers."""
    frac = None
    if best and best.get("n_reads"):
        frac = (best.get("n_exact") or 0) / best["n_reads"]
    prior, why = artefact_prior((best or {}).get("kmoch") or (best or {}).get("motif"),
                                exact_frac=frac, dupc_rejected=dupc_rejected)
    if prior == "likely":
        msg = (f"FRAMESHIFT SIGNAL, TYPE CONTESTED — and it is PROBABLY AN ARTEFACT ({why}). "
               f"Do not report a variant name from this. Build the review bundle and adjudicate by eye.")
    elif prior == "unlikely":
        msg = (f"FRAMESHIFT PRESENT, TYPE NOT SETTLED — unlikely to be an artefact ({why}). "
               f"Build the review bundle and type it by eye.")
    else:
        msg = ("FRAMESHIFT SIGNAL, TYPE NOT SETTLED ({}). Build the review bundle and adjudicate by "
               "eye.".format(why))
    return {"review": True, "artefact_prior": prior, "why": why, "message": msg,
            "how": "python -m muc1_analyzer bundle --glob <patient BAM/CRAM> --ref <GRCh38> "
                   "--vntr-ref <VNTR ref> --outdir <bundles>",
            "exact_frac": round(frac, 3) if frac is not None else None}


def render_lines(cong) -> list:
    """[(text, kind)] for a rendered report — `kind` in {head, agree, disagree, info}. Pure, no reportlab.

    The text report and the PDF must say the SAME thing: a divergence that only reaches the terminal is a
    divergence the clinician never sees. Both renderers are built from this one list."""
    if not cong or cong.get("verdict") == "INSUFFICIENT":
        return []
    out = [(f"LAYER CONGRUENCE — {cong['verdict']}",
            "disagree" if cong["verdict"] == "DIVERGENT" else "head")]
    for c in cong["checks"]:
        src = "  ·  ".join(f"{k}: {v}" for k, v in c["sources"].items())
        if c.get("agree") is None:
            out.append((f"{c['axis']}: {src}", "info"))
        else:
            out.append((f"{c['axis']} — {'agree' if c['agree'] else 'DISAGREE'}: {src}",
                        "agree" if c["agree"] else "disagree"))
    act = cong.get("action")
    if act:
        out.append((act["message"], "disagree"))
        out.append((f"→ {act['how']}", "info"))
    inf = cong.get("allele_inference")
    if inf:
        conf = "" if inf["confident"] else " [not confident — alleles too close or reads too thin]"
        out.append((f"carrier reads {inf['carrier_frac']:.1%} of the locus → likely the "
                    f"{inf['inferred']}-copy allele{conf}", "info"))
    return out


def render(cong) -> str:
    """Text block for the report. Empty when there was nothing to compare."""
    if not cong or cong.get("verdict") == "INSUFFICIENT":
        return ""
    bar = "=" * 68
    icon = {"CONCORDANT": "✓", "DIVERGENT": "⚠"}.get(cong["verdict"], "·")
    L = [bar, f"LAYER CONGRUENCE — {icon} {cong['verdict']}"]
    for c in cong["checks"]:
        src = "  ·  ".join(f"{k}: {v}" for k, v in c["sources"].items())
        mark = "" if c.get("agree") is None else (" ✓" if c["agree"] else "  ⚠ DISAGREE")
        L.append(f"  {c['axis']:<15}{mark}")
        L.append(f"      {src}")
    inf = cong.get("allele_inference")
    if inf:
        conf = "" if inf["confident"] else "  [not confident — alleles too close or reads too thin]"
        L.append(f"  carrier reads {inf['carrier_frac']:.1%} of the locus vs allele abundances "
                 f"{', '.join(f'{a:.1%}' for a in inf['abundances'])} → likely the {inf['inferred']}-copy "
                 f"allele{conf}")
    for n in cong.get("notes", []):
        L.append(f"  note: {n}")
    L.append(bar)
    return "\n".join(L) + "\n\n"
