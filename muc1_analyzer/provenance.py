"""Which MUC1_Analyzer produced this output — version + git commit, stamped into every result file.

Two user reports (2026-11) could not be tied to a code version from their outputs alone. Every JSON
and the PDF now carry {"tool", "version", "commit", "dirty"}. The commit is read from the git checkout
the package runs from; outside a checkout (a copied tree) it is "unknown" rather than a guess."""
from __future__ import annotations

import functools
import os
import subprocess

from . import __version__

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _git(*args: str) -> str | None:
    try:
        r = subprocess.run(["git", "-c", "safe.directory=*", "-C", _ROOT, *args],
                           capture_output=True, text=True, timeout=5)
    except (OSError, subprocess.SubprocessError):
        return None
    return r.stdout.strip() if r.returncode == 0 else None


@functools.lru_cache(maxsize=1)
def provenance() -> dict:
    commit = os.environ.get("MUC1_ANALYZER_COMMIT") or _git("rev-parse", "--short=7", "HEAD") or "unknown"
    dirty = None
    if commit != "unknown":
        st = _git("status", "--porcelain", "--untracked-files=no")
        dirty = bool(st) if st is not None else None
    return {"tool": "MUC1_Analyzer", "version": __version__, "commit": commit, "dirty": dirty}


def provenance_line() -> str:
    p = provenance()
    tail = "" if not p["dirty"] else " (local modifications)"
    return f"{p['tool']} v{p['version']} · commit {p['commit']}{tail}"
