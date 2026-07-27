"""MUC1_Score — two-axis genetic prognostic score for ADTKD-MUC1.

Canonical model (manuscript M1): ONSET ← frameshift position, SEVERITY ← rs4072037 splice of the mutant
allele. See docs/ROADMAP_MUC1_score.md. The legacy additive composite (`score.py`, DEL/rs4072037 weight 0)
is dev-only and imported LAZILY — the public build ships the two-axis (no DEL, no methylation).
"""
from .clinical_call import parse_clinical_call, score_call, score_from_fields

__all__ = ["parse_clinical_call", "score_call", "score_from_fields"]
__version__ = "0.1.0"


def __getattr__(name):
    # legacy composite kept reachable (dev/repro) without hard-importing score.py at package load
    if name in ("MUC1Params", "compute_score"):
        from . import score
        return getattr(score, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
