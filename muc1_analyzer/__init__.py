"""MUC1_Analyzer — long-read MUC1 VNTR analysis and the MUC1_Score for ADTKD-MUC1.

MUC1_Score (since 2026-08-31): ONE axis, the frameshift tail (repeats downstream of the variant on the
mutant allele) placed on the calibration cohort's Gaussian → percentile + four categories
(`fs_tail_score`). rs4072037 is informational only. The former two-axis fields are still emitted for
backward compatibility. The legacy additive composite (`score.py`) is dev-only and imported LAZILY.
"""
from .clinical_call import parse_clinical_call, score_call, score_from_fields

__all__ = ["parse_clinical_call", "score_call", "score_from_fields"]
__version__ = "1.0.1"


def __getattr__(name):
    # legacy composite kept reachable (dev/repro) without hard-importing score.py at package load
    if name in ("MUC1Params", "compute_score"):
        from . import score
        return getattr(score, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
