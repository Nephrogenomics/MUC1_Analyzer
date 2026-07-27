"""ONT dupC detection-limit model — adaptive-sampling sizing.

The dupC is an 8C homopolymer that ONT under-calls: a read from the MUTANT allele covering the dupC
unit shows the 8C only with probability `d` (per-read detection rate — measured ~0.16 in one sample,
~0.50 in another: it DEPENDS on the patient/context/basecaller). A healthy read shows a false 8C with
probability `f` (~0.02, homopolymer over-call noise). Question: how many mutant-allele reads are needed
to CALL the dupC at a given power? → one-sided binomial test vs the noise `f`.

Used to size an **adaptive sampling** run: `min_reads(d)` → `coverage_for(...)` gives the target
on-target coverage. ⚠ `d` must be CALIBRATED on a well-covered carrier (otherwise 3/19 → a very wide
CI): `estimate_d(k, n)` gives d + a Wilson CI, to propagate into the sizing.
Pure stdlib. A better statistic (shifted C-tract length distribution, not just the 8C bin) would raise
the effective `d` → fewer reads; see `docs/…` (upcoming upgrade).
"""
from __future__ import annotations
from math import lgamma, log, log1p, exp


def _tail(n: int, k: int, p: float) -> float:
    """P(X >= k) for X ~ Binom(n, p).

    Sum IN LOG-SPACE (lgamma + logsumexp): the naive form `comb(n, i) * p**i` overflows
    (`OverflowError: int too large to convert to float`) as soon as `comb(n, i)` exceeds ~1e308,
    i.e. for an `n` on the order of 10^3 (deep coverage: e.g. LR-PCR data).
    We walk the terms from `k`, stopping once PAST the mode and the term is negligible
    (< max − 45 in log ≈ e^−45), which bounds the cost even for a very large `n`.
    """
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    logp, logq = log(p), log1p(-p)
    lg_n1 = lgamma(n + 1)
    mode = (n + 1) * p                       # peak of the binomial pmf
    terms = []
    cur_max = float("-inf")
    for i in range(k, n + 1):
        lt = lg_n1 - lgamma(i + 1) - lgamma(n - i + 1) + i * logp + (n - i) * logq
        terms.append(lt)
        if lt > cur_max:
            cur_max = lt
        if i > mode and lt < cur_max - 45.0:  # upper tail exhausted -> stop
            break
    return min(1.0, exp(cur_max) * sum(exp(t - cur_max) for t in terms))


def crit_k(n: int, f: float = 0.02, alpha: float = 0.05) -> int:
    """Smallest number of observed 8C declaring the dupC significant vs the noise `f` (alpha threshold)."""
    for k in range(1, n + 1):
        if _tail(n, k, f) < alpha:
            return k
    return n + 1                       # never significant at this N


def power(n: int, d: float, f: float = 0.02, alpha: float = 0.05) -> float:
    """Power: probability of calling the dupC with `n` mutant reads (true rate `d`) vs noise `f`."""
    kc = crit_k(n, f, alpha)
    return _tail(n, kc, d) if kc <= n else 0.0


def min_reads(d: float, f: float = 0.02, alpha: float = 0.05, target: float = 0.95,
              nmax: int = 2000) -> int | None:
    """Mutant-allele reads (covering the dupC unit) required to reach the power `target`."""
    for n in range(3, nmax):
        if power(n, d, f, alpha) >= target:
            return n
    return None


def coverage_for(n_mut: int, spanning: float = 0.8, tag: float = 0.6) -> float:
    """useful-mutant-reads → TOTAL on-target coverage required.

    Factor 2 = diploid (the mutant allele ≈ half); `spanning` = fraction of reads covering the dupC
    unit end-to-end (depends on read length vs VNTR length);
    `tag` = fraction phased by whatshap (HP:i).
    """
    return n_mut * 2 / (spanning * tag)


def estimate_d(k: int, n: int, z: float = 1.96) -> dict:
    """Estimate `d` (= k 8C over n spanning mutant reads) + a Wilson CI (calibration on a carrier)."""
    if n == 0:
        return {"d": None, "ci": (0.0, 1.0), "k": k, "n": n}
    p = k / n
    z2 = z * z
    center = (p + z2 / (2 * n)) / (1 + z2 / n)
    half = z * ((p * (1 - p) / n + z2 / (4 * n * n)) ** 0.5) / (1 + z2 / n)
    return {"d": round(p, 3), "ci": (round(max(0.0, center - half), 3),
                                     round(min(1.0, center + half), 3)), "k": k, "n": n}


def plan(d: float, *, f: float = 0.02, alpha: float = 0.05, targets=(0.95, 0.99),
         spanning: float = 0.8, tag: float = 0.6) -> dict:
    """Full sizing for a rate `d`: mutant reads + on-target coverage per target."""
    out = {"d": d, "f": f, "alpha": alpha, "spanning": spanning, "tag": tag, "targets": {}}
    for t in targets:
        n = min_reads(d, f, alpha, t)
        out["targets"][t] = {"n_mut": n,
                             "on_target_cov": round(coverage_for(n, spanning, tag)) if n else None}
    return out


if __name__ == "__main__":
    import argparse
    import json
    ap = argparse.ArgumentParser(prog="muc1_analyzer.dupc_power",
                                 description="Adaptive-sampling sizing for the ONT dupC call")
    ap.add_argument("--d", type=float, help="8C/read detection rate (otherwise sweep a grid)")
    ap.add_argument("--calibrate", nargs=2, type=int, metavar=("K", "N"),
                    help="estimate d + CI from K observed 8C over N spanning mutant reads (carrier)")
    ap.add_argument("-f", "--noise", type=float, default=0.02)
    ap.add_argument("--spanning", type=float, default=0.8)
    ap.add_argument("--tag", type=float, default=0.6)
    args = ap.parse_args()
    if args.calibrate:
        est = estimate_d(*args.calibrate)
        lo, hi = est["ci"]
        est["plan_d"] = plan(est["d"], f=args.noise, spanning=args.spanning, tag=args.tag)
        est["plan_ci_low"] = plan(max(lo, 0.03), f=args.noise, spanning=args.spanning, tag=args.tag)
        print(json.dumps(est, indent=2, ensure_ascii=False))
    elif args.d is not None:
        print(json.dumps(plan(args.d, f=args.noise, spanning=args.spanning, tag=args.tag),
                         indent=2, ensure_ascii=False))
    else:
        print(f"{'d':>6} {'N_mut(95%)':>11} {'cov(95%)':>9} {'N_mut(99%)':>11} {'cov(99%)':>9}")
        for d in (0.10, 0.16, 0.20, 0.30, 0.40, 0.50):
            p = plan(d, f=args.noise, spanning=args.spanning, tag=args.tag)
            t95, t99 = p["targets"][0.95], p["targets"][0.99]
            print(f"{d:>6.2f} {str(t95['n_mut']):>11} {str(t95['on_target_cov']):>9} "
                  f"{str(t99['n_mut']):>11} {str(t99['on_target_cov']):>9}")
