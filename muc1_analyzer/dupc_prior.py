"""Frameshift prior/posterior matrix ON THE VNTR UNIT — mutation hotspots × ONT difficulty.

Idea (duality): within the 60 bp unit, the **homopolymers** are BOTH (i) the mutation hotspots
(replication slippage → 27dupC in the 7C tract) AND (ii) the sites ONT reads worst. The
`MUC1_Analyzer.KNOWN_REPEATS` dictionary already encodes the DESCRIBED mutations (names
`X-59dupC`, `X-58_59delCC`, `X-56_59dupCCCC`…) → an empirical catalogue of hotspots with position.

For a reference motif (default 'X') we build a per-HOMOPOLYMER matrix:
  - `slippage_score` = structural mutation prior (grows with run length);
  - `n_described`    = number of described pathogenic variants falling on this run (empirical evidence);
  - `d_prior`        = expected ONT per-read detection rate (DECREASES with the mutant run length:
                       ONT misses long homopolymers) — this is the PRIOR, to be updated by the PoN/data;
  - `reads_needed`   = mutant reads required (via `dupc_power`) to call at 95% given `d_prior`.
`posterior(...)` performs the Bayesian update: mutation prior × binomial likelihood of the `k`
8C observed over `n` reads (mutant rate `d` vs PoN noise `f`) → posterior probability of mutation.

⚠ `d_prior(L)` is a structural PRIOR (monotonic heuristic): the real `d` varies by patient (0.16
in one sample … 0.50 in another at the same L) with read quality/phasing → to be recalibrated with `dupc_pon`.
Complements the PoN (empirical, per context) with a MECHANISTIC layer (per position in the unit).
"""
from __future__ import annotations
import re

from . import dupc_power


def homopolymer_runs(seq: str, min_len: int = 3) -> list:
    """Homopolymer runs (base, start 1-based, end, length) of `seq`, length >= `min_len`."""
    return [{"base": m.group(1), "start": m.start() + 1, "end": m.end(), "len": m.end() - m.start()}
            for m in re.finditer(r"(.)\1{%d,}" % (min_len - 1), seq)]


_VAR = re.compile(r"^[A-Za-z0-9+]+-(\d+)(?:_(\d+))?(dup|del|ins)([ACGT]+)$")


def described_mutations(known_repeats: dict) -> list:
    """Parse the variant names from the dictionary → [{name, pos1, pos2, kind, bases}]."""
    out = []
    for name in known_repeats.values():
        m = _VAR.match(name)
        if m:
            out.append({"name": name, "pos1": int(m.group(1)),
                        "pos2": int(m.group(2)) if m.group(2) else int(m.group(1)),
                        "kind": m.group(3), "bases": m.group(4)})
    return out


def slippage_prior(L: int) -> int:
    """Structural slippage mutation prior: grows with run length (∝ (L-2)²)."""
    return max(0, L - 2) ** 2


def d_prior(mut_len: int) -> float:
    """Expected ONT per-read detection rate by MUTANT homopolymer length (decreasing).

    Monotonic heuristic (prior): ONT ~reliable up to ~4-5, degrades beyond. Bounded [0.1, 0.95].
    Roughly anchored on our data (8C tract → observed d 0.16-0.50). To be recalibrated via `dupc_pon`.
    """
    return round(min(0.95, max(0.10, 0.95 - 0.11 * (mut_len - 4))), 3)


def hotspot_matrix(motif: str = "X", *, alpha: float = 0.05, target: float = 0.95) -> list:
    """Per-homopolymer prior matrix of motif `motif`: mutation × ONT difficulty × reads required."""
    from .detectors.vntr_scaffold import _import_analyzer
    kr = _import_analyzer().KNOWN_REPEATS
    seq = {v: k for k, v in kr.items()}[motif]
    described = described_mutations(kr)
    rows = []
    for run in homopolymer_runs(seq):
        # described variants falling on/adjacent to the run (dup at the end, del/ins within the window)
        near = [d for d in described if run["start"] <= d["pos2"] <= run["end"] + 1
                or run["start"] <= d["pos1"] <= run["end"] + 1]
        mut_len = run["len"] + 1                       # dup → homopolymer +1 (the most frequent)
        d = d_prior(mut_len)
        rows.append({
            "run": f"{run['len']}x{run['base']}", "pos": f"{run['start']}-{run['end']}",
            "len": run["len"], "slippage_score": slippage_prior(run["len"]),
            "n_described": len(near), "described": sorted({d_["name"] for d_ in near}),
            "d_prior": d, "reads_needed": dupc_power.min_reads(d, alpha=alpha, target=target)})
    rows.sort(key=lambda r: (-r["slippage_score"], -r["n_described"]))
    return rows


def posterior(prior_mut: float, k: int, n: int, *, d: float, f: float = 0.02) -> float:
    """Posterior probability of mutation at the hotspot: prior × binomial likelihood of the k 8C / n reads.

    P(mut | data) = [P(data|mut)·prior] / [P(data|mut)·prior + P(data|¬mut)·(1-prior)]
    with P(data|mut)=Binom(k;n,d) (mutant allele, detection rate d) and P(data|¬mut)=Binom(k;n,f) (noise).
    """
    from math import comb
    like_mut = comb(n, k) * d ** k * (1 - d) ** (n - k)
    like_not = comb(n, k) * f ** k * (1 - f) ** (n - k)
    num = like_mut * prior_mut
    den = num + like_not * (1 - prior_mut)
    return round(num / den, 4) if den else 0.0


if __name__ == "__main__":
    import argparse
    import json
    ap = argparse.ArgumentParser(prog="muc1_analyzer.dupc_prior",
                                 description="Frameshift hotspots × ONT difficulty matrix (VNTR unit)")
    ap.add_argument("--motif", default="X", help="reference motif (default X = wild-type)")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()
    mat = hotspot_matrix(args.motif)
    if args.json:
        print(json.dumps(mat, indent=2, ensure_ascii=False))
    else:
        print(f"Frameshift hotspots of motif {args.motif} (sorted by mutation prior):\n")
        print(f"{'run':>6} {'pos':>7} {'slip':>5} {'n_var':>5} {'d_prior':>7} {'reads95':>7}  described variants")
        for r in mat:
            print(f"{r['run']:>6} {r['pos']:>7} {r['slippage_score']:>5} {r['n_described']:>5} "
                  f"{r['d_prior']:>7} {str(r['reads_needed']):>7}  {', '.join(r['described']) or '—'}")
