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


def check(*, dupc=None, report=None, arbiter=None, hap_variant=None, hap_repeat=None):
    """Compare the layers on VARIANT / REPEAT / CARRIER ALLELE. Pure.

    `dupc`   = a dupc_dispatch result dict ; `report` = a pcr_report result dict ;
    `arbiter`= an auto_length dict (alleles/counts) ; `hap_*` = the haplotype caller's nomenclature call.
    Returns {verdict, checks:[…], allele_inference, notes}. `verdict` is CONCORDANT when every axis that
    could be compared agrees, DIVERGENT when at least one disagrees, INSUFFICIENT when nothing could be.
    """
    res = (dupc or {}).get("result", dupc) or {}
    variants = (report or {}).get("variants") or []
    # Credible variant = most EXACTLY-decoded reads, among those decoded credibly at all.
    # Ranking by exact FRACTION alone was wrong: a small perfect call outranks a large real one (observed on
    # 4 cohort samples, where a 22/22 X:60dupA displaced the true 59dupC and the congruence check flagged it).
    # Ranking by RAW read count is equally wrong (a 122-read 6 %-exact decoding artefact would win).
    # n_exact combines both: support AND decode quality. The fraction floor only discards artefacts.
    best = None
    if variants:
        credible = [v for v in variants
                    if (v.get("n_exact") or 0) / max(v.get("n_reads") or 1, 1) >= MIN_EXACT_FRAC]
        if credible:
            best = max(credible, key=lambda v: v.get("n_exact") or 0)

    checks, notes = [], []

    # ── variant identity ──
    d_pos = res.get("status") == "positive" or bool(res.get("called"))
    r_var = _norm_variant(best.get("kmoch") or best.get("motif")) if best else None
    h_var = _norm_variant(hap_variant)
    seen = {k: v for k, v in (("pcr_report", r_var), ("haplotype_caller", h_var)) if v}
    if d_pos or seen:
        agree = (len(set(seen.values())) <= 1) and (not d_pos or not seen or "dupC" in (r_var or h_var or ""))
        checks.append({"axis": "variant", "agree": bool(agree),
                       "sources": {"dupc_caller": "positive" if d_pos else "not called", **seen}})

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
    return {"verdict": verdict, "checks": checks, "allele_inference": inf, "notes": notes,
            "best_variant": best}


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
