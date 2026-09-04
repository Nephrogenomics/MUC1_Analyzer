#!/usr/bin/env python3
"""clinical_call — score a MUC1 carrier from a VISUALLY/clinically finalized call string.

Auto-detection (PCR/AS/WGS → dupC/frameshift caller) is the default path. But some variants are
finalized BY EYE (e.g. del8_27, invisible to the automated caller; or a coverage-floored het with a long
VNTR, as seen in some samples). For those, the human provides the call in the lab nomenclature, and this module feeds
the SAME single-axis model:

    "64 repeats | 80 repeats (del8_27 @ repeat 37)"
      → short = 64, mutant = the 80-repeat (LONG) allele, del8_27 at repeat 37
      → SCORE AXIS : frameshift tail = mut_len − position = 80 − 37 = 43 repeats downstream of the variant
                     → placed on the Gaussian cohort curve (fs_tail_score) → 4-category verdict
      → rs4072037  : splice base of the mutant allele, reported FOR INFORMATION only (NOT scored)
      → descriptor : ratio = 80/64 = 1.25, mutant is the long allele

The MUC1_Score depends ONLY on the frameshift tail length: a longer tail → more MUC1fs neoprotein → a
more severe phenotype (concordant with Vrbacká/Kmoch 2025). rs4072037 has no independent prognostic
value and is reported but not weighted; ONSET has been abandoned. See `fs_tail_score` / `config.FS_TAIL_CAL`.

The parenthetical `(variant @ repeat N)` marks the MUTANT allele; `@` is canonical, `~` also accepted
(the historical `DATA-Patients-MUC1-vntr.xlsx` uses `~`). rs4072037 can likewise be given by eye
(`rs4072037_mut='C'|'T'`) or, one day, by dRNA-seq — highest precedence; else inferred from dosage + the
cis rule (T cis-short). `w_onset`/`w_splice` are accepted for signature compatibility but no longer
affect the score.
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
                      dosage=None, cis_anchor=None, mut_is_long=None, position_confident=None,
                      w_onset: float = 2.0, w_splice: float = 1.0) -> dict:
    """PURE SINGLE-AXIS core, shared by the manual (clinical-call string) and the AUTO (analyzer-JSON /
    score_one detection) paths. Used by `run`, `run_pcr` and `from-fastq` alike.

    The MUC1_Score depends on ONE quantity: the FRAMESHIFT TAIL length = repeats downstream of the
    variant on the mutant allele = `mut_len − position`. A longer tail → more MUC1fs neoprotein → a more
    severe phenotype. The tail is placed on the Gaussian cohort curve (`fs_tail_score`, calibrated on
    sheet Index_v2, n=33) and mapped to a 4-category verdict (very mild / mild / severe / very severe).

    rs4072037 is NO LONGER part of the score — it has no independent prognostic value. It is still
    resolved (`mut_T`, `splice_base`) and returned FOR INFORMATION only; `severity_score` is kept as a
    back-compat key but is now always None. `onset_index`/`ratio` are likewise still emitted for
    back-compat and audit but no longer drive anything.

    `position_confident=False` (from `pcr_report`'s `repeat_confident`) does NOT withhold the score, but
    is recorded as a caveat (`fs_confident=False` + `fs_caveat`): the tail is only as reliable as the
    position it is measured from. `w_onset`/`w_splice` are accepted for signature compatibility and no
    longer affect the result."""
    from .fs_tail_score import fs_tail_score
    if mut_is_long is None and mut_len is not None and healthy_len is not None and mut_len != healthy_len:
        mut_is_long = mut_len > healthy_len
    onset_index = (1 - position / mut_len) if (position is not None and mut_len) else None
    # ── the score axis: frameshift tail = repeats downstream of the variant ──
    fs_tail = (mut_len - position) if (mut_len is not None and position is not None) else None
    fst = fs_tail_score(fs_tail, position_confident=position_confident)
    out = {
        "mut_len": mut_len, "healthy_len": healthy_len, "mut_is_long": mut_is_long, "position": position,
        "ratio": round(mut_len / healthy_len, 3) if (mut_len and healthy_len) else None,
        "onset_index": round(onset_index, 4) if onset_index is not None else None,
        "onset_score": None,                       # onset abandoned — kept as a key for back-compat only
        "position_confident": position_confident,
        "score_axis": "frameshift_tail",
    }
    out.update(fst)                                # fs_tail, fs_z, fs_percentile, fs_category, fs_label, …
    # rs4072037 — resolved and reported FOR INFORMATION ONLY (no severity contribution).
    if rs4072037_mut in ("C", "T"):
        mut_T, src = (rs4072037_mut == "T"), "visual/dRNAseq"
    elif dosage is not None:
        mut_T, src = _mut_T_from_dosage(dosage, mut_is_long, cis_anchor), "inferred(dosage/cis)"
        if str(dosage) in ("0", "2"):
            src = "homozygous_determined"           # not inferred — see `phasing_note`
    else:
        mut_T, src = None, None
    out.update(mut_T=mut_T,
               splice_base=(None if mut_T is None else ("T" if mut_T else "C")),
               severity_score=None,                 # rs4072037 NOT weighted any more (informational)
               splice_source=src,
               rs4072037_informational=True)
    out.update(phasing_note(dosage, mut_len, healthy_len))
    return out


def phasing_note(dosage, mut_len, healthy_len) -> dict:
    """Is phasing the variant against rs4072037 WORTH DOING, and if not, why? Pure.

    The clinical rule (LM, 2026-07-30), and it saves most of the work: rs4072037 modulates severity
    through the exon-2 splice, so what matters is the splice status OF THE MUTANT ALLELE. When the SNP is
    HOMOZYGOUS, that status is the same whichever allele carries the frameshift —
      · **TT** → the splice is promoted on BOTH alleles → protected either way;
      · **CC** → promoted on NEITHER → severe either way.
    Phasing then cannot change the answer, and reporting "unphased" as if something were missing is
    misleading: nothing is missing.

    The effort is worth it in exactly ONE configuration: **heterozygous rs4072037**. And it becomes hard
    only when the two VNTRs are also the SAME LENGTH, because length no longer separates the alleles.
    There the discriminator is not length at all — it is the motif composition and its order along reads
    that span the array END TO END and reach rs4072037: one such read carries the variant AND the SNP
    base, so it settles the phase by itself, with no statistical phasing.

    Returns `phasing_informative` (does phasing change the severity?), `phasing_resolved` (do we have the
    answer?) and a human-readable `phasing_note`."""
    d = str(dosage) if dosage is not None else None
    same_len = (mut_len is not None and healthy_len is not None and mut_len == healthy_len)
    if d == "2":
        return {"phasing_informative": False, "phasing_resolved": True,
                "phasing_note": "rs4072037 TT — the splice is promoted on BOTH alleles, so the mutant "
                                "allele is protected whichever one it is. Phasing is uninformative here, "
                                "not missing."}
    if d == "0":
        return {"phasing_informative": False, "phasing_resolved": True,
                "phasing_note": "rs4072037 CC — the splice is promoted on NEITHER allele, so the mutant "
                                "allele is severe whichever one it is. Phasing is uninformative here, "
                                "not missing."}
    if d == "1":
        if same_len:
            return {"phasing_informative": True, "phasing_resolved": False,
                    "phasing_note": "rs4072037 HETEROZYGOUS and the two VNTRs are the SAME LENGTH — the "
                                    "one configuration where phasing changes the severity AND length "
                                    "cannot separate the alleles. Resolve it on the READS: a read "
                                    "spanning the array end to end and reaching rs4072037 carries the "
                                    "variant and the SNP base together, so it settles the phase on its "
                                    "own. Motif composition and order along such reads is the "
                                    "discriminator, not length."}
        return {"phasing_informative": True, "phasing_resolved": True,
                "phasing_note": "rs4072037 heterozygous with alleles of different length — phasing "
                                "matters and length provides the handle (cis rule)."}
    return {"phasing_informative": None, "phasing_resolved": False,
            "phasing_note": "rs4072037 not genotyped — severity axis unset."}


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
