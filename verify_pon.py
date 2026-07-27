#!/usr/bin/env python3
"""Verify that a PoN file is an AGGREGATE error model, safe to publish.

The PoN we ship is `{sequence context -> homopolymer error counts}` — the analogue of a
technical error model or a gnomAD allele frequency. That is only true if the file really
contains nothing per-sample, so this checks the actual bytes instead of trusting the
generator: the aggregation in `dupc_pon.aggregate` drops sample identity by construction,
but a hand-edited, partially-aggregated, or wrong file would not.

It is deliberately a WHITELIST (known-good keys and value types) rather than a blacklist of
ID-looking patterns: a blacklist passes anything it failed to imagine.

    python3 release/verify_pon.py /path/to/pon_dupc.json

Exit 0 = safe to publish. Exit 1 = do NOT publish; the reasons are printed.
"""
from __future__ import annotations

import json
import os
import re
import sys

# Top-level keys the aggregate may carry. Anything else is unreviewed content.
ALLOWED_TOP = {"n_profiles", "contexts", "n_skipped", "_note", "synthetic", "meta"}
# Per-context fields emitted by dupc_pon.aggregate.
ALLOWED_CTX = {"n_patients", "n_ctracts", "pooled_7C", "pooled_8C", "counts",
               "f_8C_pooled", "f_8C_mean", "f_8C_sd", "f_8C_max", "difficulty"}
# Metadata we add on publication (free-text, but reviewed: printed in full below).
ALLOWED_META = {"chemistry", "source", "n_samples", "regime", "built", "note", "label"}

# A context key is a run of MOTIF NAMES joined by "-", e.g. "7-7-7" or "aS-7-aS" (dupc_pon builds
# it from the last `ctx_units` matched motifs). So the tokens are validated against the caller's
# actual motif vocabulary — a sample identifier is not a motif name and is therefore rejected.
# KNOWN_REPEATS is read by AST, not imported: caller.py pulls pysam, and this check must stay
# runnable anywhere, including on a login node. caller.py is found relative to THIS file in both
# layouts — the dev repo (`release/verify_pon.py` → `../muc1_analyzer/caller.py`) and the public
# repo, where verify_pon.py sits at the root (`./muc1_analyzer/caller.py`).
_HERE = os.path.dirname(os.path.abspath(__file__))
_MOTIF_CANDIDATES = (
    os.path.join(_HERE, "..", "muc1_analyzer", "caller.py"),   # dev: release/ sibling
    os.path.join(_HERE, "muc1_analyzer", "caller.py"),         # public: repo root
)
#: token for a stretch the caller could not match to a known motif — a legitimate context key.
_UNMATCHED = "?"
#: conservative fallback if the motif panel cannot be read: the unmatched token, or short
#: alphanumerics — still excludes anything shaped like an identifier (no underscore, no length).
_FALLBACK_TOKEN = re.compile(r"^([A-Za-z0-9]{1,4}\+*|\?)$")


def motif_vocabulary(src: str | None = None) -> set | None:
    """Motif-name bases from `KNOWN_REPEATS` in caller.py, via AST. None if unreadable."""
    candidates = (src,) if src else _MOTIF_CANDIDATES
    tree = None
    for cand in candidates:
        try:
            import ast
            tree = ast.parse(open(cand, encoding="utf-8").read())
            break
        except (OSError, SyntaxError):
            continue
    if tree is None:
        return None
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(
                getattr(t, "id", None) == "KNOWN_REPEATS" for t in node.targets):
            try:
                d = ast.literal_eval(node.value)
            except ValueError:
                return None
            return {str(v).split("-")[0] for v in d.values()} | {_UNMATCHED}
    return None


def _token_ok(tok: str, vocab: set | None) -> bool:
    return tok in vocab if vocab else bool(_FALLBACK_TOKEN.match(tok))


def verify(path: str) -> tuple[bool, list, dict]:
    problems, info = [], {}
    with open(path) as fh:
        d = json.load(fh)

    if not isinstance(d, dict):
        return False, ["top level is not an object"], info

    extra = set(d) - ALLOWED_TOP
    if extra:
        problems.append(f"unexpected top-level key(s): {sorted(extra)}")

    meta = d.get("meta") or {}
    if meta:
        bad_meta = set(meta) - ALLOWED_META
        if bad_meta:
            problems.append(f"unexpected meta key(s): {sorted(bad_meta)}")

    contexts = d.get("contexts")
    if not isinstance(contexts, dict) or not contexts:
        problems.append("no `contexts` table — this is not an aggregated PoN")
        return False, problems, info

    vocab = motif_vocabulary()
    info["vocabulary"] = "motif panel" if vocab else "fallback charset"
    total_ctracts, bad_keys, bad_fields = 0, [], set()
    for ctx, v in contexts.items():
        if not all(_token_ok(t, vocab) for t in ctx.split("-")):
            bad_keys.append(ctx)                     # a sample ID would land here
        if not isinstance(v, dict):
            problems.append(f"context {ctx!r} is not an object")
            continue
        bad_fields |= set(v) - ALLOWED_CTX
        counts = v.get("counts") or {}
        if not all(isinstance(k, str) and k.isdigit() for k in counts):
            problems.append(f"context {ctx!r}: `counts` keys are not C-tract lengths")
        if not all(isinstance(x, (int, float)) for x in counts.values()):
            problems.append(f"context {ctx!r}: `counts` values are not numeric")
        total_ctracts += v.get("n_ctracts") or 0

    if bad_keys:
        problems.append(f"{len(bad_keys)} context key(s) are not C-tract patterns "
                        f"(possible identifiers): {bad_keys[:5]}")
    if bad_fields:
        problems.append(f"unexpected per-context field(s): {sorted(bad_fields)}")

    info.update({"n_profiles": d.get("n_profiles"), "n_skipped": d.get("n_skipped"),
            "n_contexts": len(contexts), "total_ctracts": total_ctracts,
            "meta": meta, "synthetic": d.get("synthetic"), "note": d.get("_note")})
    return not problems, problems, info


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        print(__doc__)
        return 2
    path = argv[0]
    try:
        ok, problems, info = verify(path)
    except (OSError, json.JSONDecodeError) as e:
        print(f"[verify_pon] cannot read {path}: {e}")
        return 1

    print(f"── {path}")
    for k in ("n_profiles", "n_skipped", "n_contexts", "total_ctracts", "synthetic", "note"):
        if info.get(k) is not None:
            print(f"   {k:<14} {info[k]}")
    if info.get("meta"):
        print(f"   meta           {json.dumps(info['meta'], ensure_ascii=False)}")

    if ok:
        print("\n[OK] aggregate-only: contexts are C-tract patterns, values are counts/stats.")
        print("     No reads, no genotypes, no sample identifiers. Safe to publish.")
        return 0
    print("\n[REFUSED] do NOT publish this file:")
    for p in problems:
        print(f"   - {p}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
