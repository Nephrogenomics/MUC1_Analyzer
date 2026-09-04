#!/usr/bin/env python3
"""fs_tail_score — the SINGLE-AXIS MUC1_Score.

The MUC1_Score depends on ONE quantity only: the frameshift TAIL length, i.e. the number of VNTR
repeats DOWNSTREAM of the pathogenic variant on the mutant allele:

    fs_tail = mut_len − position

A longer tail carries more of the mutant MUC1fs neoprotein, and that tracks a MORE SEVERE phenotype
(earlier ESRD), concordant with Vrbacká/Kmoch 2025 (fast progressors carry more frameshifted repeats).
The tail is placed on a Gaussian curve fitted to the phenotyped cohort (`config.FS_TAIL_CAL`, sheet
Index_v2, n=33, μ=36.94, σ=14.07, normality confirmed) and mapped to a percentile → one of 4 categories:

    percentile  < 0.10          →  very mild    "Genotype associated with very mild phenotype"
    0.10 ≤ pct  < 0.50          →  mild         "Genotype associated with mild phenotype"
    0.50 ≤ pct  ≤ 0.90          →  severe       "Genotype associated with severe phenotype"
    percentile  > 0.90          →  very severe  "Genotype associated with very severe phenotype"

Design notes:
  · rs4072037 is NOT part of this score (no independent prognostic value; reported for information only,
    see clinical_call.score_from_fields and the PDF renderer).
  · The old two-axis fields (onset_index, severity_score, ratio) are still emitted by score_from_fields
    for backward compatibility, but they no longer drive the verdict.
  · PURE and stdlib-only (math.erf) — importable anywhere (no pysam, no numpy/scipy), so a light JSON
    utility can score without pulling heavy deps.
"""
from __future__ import annotations
import math

from .config import FS_TAIL_CAL

# Category colours reused by the PDF renderer (green → red, short → long tail).
_VERY_MILD = ("very_mild", "Genotype associated with very mild phenotype", "#1B7A2F")
_MILD = ("mild", "Genotype associated with mild phenotype", "#4CAF50")
_SEVERE = ("severe", "Genotype associated with severe phenotype", "#E67E22")
_VERY_SEVERE = ("very_severe", "Genotype associated with very severe phenotype", "#C00000")


def _normal_cdf(x: float, mu: float, sd: float) -> float:
    """Φ((x−μ)/σ) via math.erf — no scipy/numpy dependency."""
    return 0.5 * (1.0 + math.erf((x - mu) / (sd * math.sqrt(2.0))))


def category_for_percentile(pct: float):
    """(key, label, color) for a Gaussian percentile. Boundaries: <0.10 very mild, <0.50 mild,
    ≤0.90 severe, >0.90 very severe. Symmetric around the median (0.50 falls in the severe side)."""
    if pct < 0.10:
        return _VERY_MILD
    if pct < 0.50:
        return _MILD
    if pct <= 0.90:
        return _SEVERE
    return _VERY_SEVERE


def fs_tail_score(fs_tail, *, cal: dict | None = None, position_confident=None) -> dict:
    """PURE single-axis score for a frameshift tail length (repeats).

    Returns a dict with the tail, its z-score and Gaussian percentile against the calibration cohort,
    and the 4-category verdict (`category` key + human `label`). `fs_tail=None` (non-carrier / no
    position) returns an all-None result with `category=None`.

    `position_confident=False` (forwarded from a caller that could not pin the repeat position — e.g.
    `pcr_report.repeat_confident`) does NOT suppress the score, but is recorded as a caveat: the tail
    length is only as trustworthy as the position it is measured from, so the report flags it."""
    cal = cal or FS_TAIL_CAL
    mu, sd, n = cal["mu"], cal["sd"], cal["n"]
    if fs_tail is None:
        return {"fs_tail": None, "fs_z": None, "fs_percentile": None,
                "fs_category": None, "fs_label": None, "fs_color": None,
                "fs_mu": mu, "fs_sd": sd, "fs_n": n,
                "fs_confident": (position_confident is not False)}
    z = (fs_tail - mu) / sd
    pct = _normal_cdf(fs_tail, mu, sd)
    key, label, color = category_for_percentile(pct)
    out = {
        "fs_tail": fs_tail,
        "fs_z": round(z, 3),
        "fs_percentile": round(pct, 4),
        "fs_category": key,
        "fs_label": label,
        "fs_color": color,
        "fs_mu": mu, "fs_sd": sd, "fs_n": n,
        "fs_confident": (position_confident is not False),
    }
    if position_confident is False:
        out["fs_caveat"] = ("frameshift POSITION not reliably measured (the carrier reads did not agree "
                            "on it) → the tail length, and therefore this category, is provisional")
    return out
