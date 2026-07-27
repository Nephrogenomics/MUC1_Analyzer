#!/usr/bin/env python3
"""clinical_call — score a MUC1 carrier from a VISUALLY/clinically finalized call string.

Auto-detection (PCR/AS/WGS → dupC/frameshift caller) is the default path. But some variants are
finalized BY EYE (e.g. del8_27, invisible to the automated caller; or a coverage-floored het with a long
VNTR, as seen in some samples). For those, the human provides the call in the lab nomenclature, and this module feeds
the SAME two-axis model:

    "64 repeats | 80 repeats (del8_27 @ repeat 37)"
      → short = 64, mutant = the 80-repeat (LONG) allele, del8_27 at repeat 37
      → Axis ONSET   : onset_index = 1 - 37/80 = 0.54   (fraction translated as MUC1fs neoprotein)
      → Axis SEVERITY: rs4072037 splice of the mutant allele (visual/dRNAseq base > dosage/cis inference)
      → descriptor   : ratio = 80/64 = 1.25, mutant is the long allele

The parenthetical `(variant @ repeat N)` marks the MUTANT allele; `@` is canonical, `~` also accepted
(the historical `DATA-Patients-MUC1-vntr.xlsx` uses `~`). rs4072037 can likewise be given by eye
(`rs4072037_mut='C'|'T'`) or, one day, by dRNA-seq — highest precedence; else inferred from dosage + the
cis rule (T cis-short). Weights `w_onset`/`w_splice` are configurable, default = direction-of-effect
(NOT a calibrated model; a `calibrate()` hook is left for when the cohort grows).
"""
from __future__ import annotations
import re

# <N> repeats  [ ( <variant>  @|~  repeat <pos> ) ]   — the parenthetical marks the MUTANT allele
_ALLELE = re.compile(
    r"(\d+)\s*repeats\s*(?:\(\s*([\w.]+?)\s*[@~]\s*repeat\s*(\d+)\s*\))?",
    re.IGNORECASE,
)


def parse_clinical_call(s: str) -> dict:
    """PURE: parse '<a> repeats [(variant @ repeat N)] | <b> repeats [(variant @ repeat N)]'.
    Returns the structural layer; `carrier=False` when neither allele carries a variant annotation.
    Raises ValueError on no allele / a variant annotated on >1 allele (ADTKD-MUC1 is mono-allelic)."""
    alleles = []
    for part in str(s).split("|"):
        m = _ALLELE.search(part)
        if not m:
            continue
        alleles.append({"length": int(m.group(1)),
                        "variant": m.group(2),
                        "position": int(m.group(3)) if m.group(3) else None})
    if not alleles:
        raise ValueError(f"clinical_call: no '<N> repeats' allele in {s!r}")
    lens = sorted(a["length"] for a in alleles)
    muts = [a for a in alleles if a["variant"]]
    if len(muts) > 1:
        raise ValueError(f"clinical_call: a variant is annotated on >1 allele (mono-allelic expected): {s!r}")
    if not muts:
        return {"carrier": False, "short": lens[0], "long": lens[-1]}
    mut = muts[0]
    others = [a for a in alleles if a is not mut]
    healthy_len = others[0]["length"] if others else mut["length"]   # homozygous-length: healthy == mut length
    mut_len, pos = mut["length"], mut["position"]
    onset_index = (1 - pos / mut_len) if (pos is not None and mut_len) else None
    return {
        "carrier": True,
        "short": lens[0], "long": lens[-1],
        "mut_len": mut_len, "healthy_len": healthy_len,
        "mut_is_long": (mut_len > healthy_len) if mut_len != healthy_len else None,
        "variant": mut["variant"], "position": pos,
        "ratio": round(mut_len / healthy_len, 3) if healthy_len else None,
        "onset_index": round(onset_index, 4) if onset_index is not None else None,
    }


def apply_arbiter_lengths(frag: dict, arbiter_alleles, carrier_len=None) -> dict:
    """Replace the CALLER's allele lengths with the ARBITER's in a score fragment. Pure.

    The arbiter is authoritative — it is the alignment-free scale validated against clinical truths — and
    the caller's haplotype lengths are only a by-product of consensus reconstruction. On a real carrier the
    arbiter measured 45/77 identically from two entry points while the caller returned 43/44, losing the
    PCR-depleted long allele; `onset_index` was then computed on alleles that do not exist.

    WHICH arbiter allele is the mutant comes from `carrier_len` — the length of the allele the dupC caller
    actually found the variant on, a MEASUREMENT. Without it we fall back to preserving the caller's
    long/short ordering, which is a deduction and is recorded as such.

    Returns a new fragment; the caller's values are kept under `caller_lengths` so nothing is lost."""
    out = dict(frag or {})
    alleles = sorted(int(a) for a in (arbiter_alleles or []))
    if not alleles or not out.get("mutation_present"):
        return out
    old_mut, old_healthy = out.get("vntr_len_mut"), out.get("vntr_len_healthy")
    if len(alleles) == 1:
        mut = healthy = alleles[0]                       # length-homozygous: both alleles measure the same
        source = "arbiter (length-homozygous)"
    elif carrier_len is not None:
        mut = min(alleles, key=lambda a: abs(a - int(carrier_len)))
        healthy = [a for a in alleles if a != mut][0] if len(set(alleles)) > 1 else mut
        source = f"arbiter, carrier allele MEASURED at {carrier_len}"
    else:
        # no measurement of which allele carries it: keep the caller's ordering, and say so
        mut_is_long = out.get("mut_is_long")
        if mut_is_long is None and old_mut is not None and old_healthy is not None:
            mut_is_long = old_mut > old_healthy
        mut = max(alleles) if mut_is_long else min(alleles)
        healthy = min(alleles) if mut_is_long else max(alleles)
        source = "arbiter lengths, mutant allele INFERRED from the caller's ordering (not measured)"
    out["vntr_len_mut"], out["vntr_len_healthy"] = mut, healthy
    out["length_source"] = source
    out["caller_lengths"] = {"vntr_len_mut": old_mut, "vntr_len_healthy": old_healthy}
    return out


def _mut_T_from_dosage(dosage, mut_is_long, cis_anchor=None):
    """Mutant allele carries rs4072037-T? Reuses the cohort cis rule (T cis-short). None if unresolvable
    (het with length-homozygous alleles → no cis handle → caller must pass rs4072037_mut)."""
    from .haplotype_score import _mut_splice_is_T
    if str(dosage) == "1" and mut_is_long is None:
        return None                                  # length-hom het → cannot cis-infer
    return _mut_splice_is_T(dosage, cis_anchor, mut_is_long)


def score_from_fields(mut_len, healthy_len, position, *, rs4072037_mut: str | None = None,
                      dosage=None, cis_anchor=None, mut_is_long=None,
                      w_onset: float = 2.0, w_splice: float = 1.0) -> dict:
    """PURE two-axis core, shared by the manual (clinical-call string) and the AUTO (analyzer-JSON /
    score_one detection) paths. Axis A ONSET = 1 − position/mut_len ; Axis B SEVERITY = mutant-allele
    rs4072037 splice (visual/dRNA-seq base > dosage+cis inference). No DEL, no methylation."""
    if mut_is_long is None and mut_len is not None and healthy_len is not None and mut_len != healthy_len:
        mut_is_long = mut_len > healthy_len
    onset_index = (1 - position / mut_len) if (position is not None and mut_len) else None
    out = {
        "mut_len": mut_len, "healthy_len": healthy_len, "mut_is_long": mut_is_long, "position": position,
        "ratio": round(mut_len / healthy_len, 3) if (mut_len and healthy_len) else None,
        "onset_index": round(onset_index, 4) if onset_index is not None else None,
        "onset_score": round(w_onset * onset_index, 3) if onset_index is not None else None,
    }
    if rs4072037_mut in ("C", "T"):
        mut_T, src = (rs4072037_mut == "T"), "visual/dRNAseq"
    elif dosage is not None:
        mut_T, src = _mut_T_from_dosage(dosage, mut_is_long, cis_anchor), "inferred(dosage/cis)"
    else:
        mut_T, src = None, None
    out.update(mut_T=mut_T,
               splice_base=(None if mut_T is None else ("T" if mut_T else "C")),
               severity_score=(None if mut_T is None else (-w_splice if mut_T else +w_splice)),
               splice_source=src)                                   # T→protected(−), C→severe(+)
    return out


def score_call(s: str, *, rs4072037_mut: str | None = None, dosage=None, cis_anchor=None,
               w_onset: float = 2.0, w_splice: float = 1.0) -> dict:
    """Two-axis score from a clinical/visual call STRING. Parses the nomenclature, then delegates the
    scoring to `score_from_fields`. Provenance `position_source='clinical_visual'`."""
    c = parse_clinical_call(s)
    if not c["carrier"]:
        return dict(c, position_source="clinical_visual", onset_score=None, severity_score=None,
                    mut_T=None, splice_base=None, splice_source=None)
    f = score_from_fields(c["mut_len"], c["healthy_len"], c["position"], rs4072037_mut=rs4072037_mut,
                          dosage=dosage, cis_anchor=cis_anchor, mut_is_long=c["mut_is_long"],
                          w_onset=w_onset, w_splice=w_splice)
    return dict(c, **f, position_source="clinical_visual", variant=c["variant"])


# calibrate(cohort) — hook left UNWIRED on purpose: w_onset/w_splice stay at direction-of-effect defaults
# until the phenotyped cohort is large enough to fit them (n=12 now → premature). See ROADMAP_MUC1_score.md.
def calibrate(*_a, **_k):
    raise NotImplementedError("weight calibration deferred until the cohort grows (n=12 → underpowered)")
