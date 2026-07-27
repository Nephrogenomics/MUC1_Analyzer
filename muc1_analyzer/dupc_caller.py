"""MULTI-RESOLUTION dupC caller — automates the de-dilution (for some samples) + PoN null + posterior.

We showed (for some sample) that a dupC can be DROWNED in a context shared at `ctx_units=3` (3/298 ≈ null)
and CONCENTRATE at `ctx_units=4` (3/9, p=1.4e-4). This caller scans several context resolutions,
tests each candidate context of the mutant allele vs the per-context PoN null (`dupc_pon.f_for_context`),
corrects the number of tests (Bonferroni), keeps the most significant, and returns the calibrated posterior.

`d_assumed` = assumed 8C/read detection rate IF carrier (prior, ~0.16-0.5 measured) for the posterior.
⚠ The null `f` must come from a PoN **independent** of the tested set (otherwise circular — cf. `dupc_detection.md`).
"""
from __future__ import annotations
from math import comb

from .dupc_power import _tail                     # P(X>=k | Binom(n,f)) = one-sided binomial test
from .dupc_prior import posterior as _posterior
from .detectors.dupc_pon import f_for_context

# Decoding ARTIFACT contexts (`match_motifs` mis-decode), not biology: `58_59delCC` fragments
# mint rare contexts that concentrate spurious 8C (cf. spec. 0.66 → 32 FP, `docs/dupc_detection.md`).
# `?` = undecoded units. We refuse to call on them.
def _is_artifact_ctx(ctx: str, mut_len: int = 8) -> bool:
    """Decoding ARTIFACT context to mask. `?` = undecoded unit → always unreliable.
    `58_59delCC` in the context = a neighbor decoded as delCC: for the dupC (8C) this is a sign
    of local decoding instability → mask; for the **delCC (5C) it IS THE SIGNAL** (a real delCC
    decodes this way, cf. one sample: without this guard, `tier=none`; with it, p_bonf 2e-15) → do NOT mask."""
    if "?" in ctx:
        return True
    return mut_len != 5 and "58_59delCC" in ctx


def _fisher_right(a: int, b: int, c: int, d: int) -> float:
    """One-sided p (Fisher exact) that the `mut` row [a 8C, b 7C] is ENRICHED in 8C vs `healthy`
    [c 8C, d 7C]. Paired intra-patient null (same chemistry/coverage/context). No dependency."""
    r1, r2, c1, n = a + b, c + d, a + c, a + b + c + d
    if r1 == 0 or r2 == 0 or c1 == 0 or c1 == n:
        return 1.0                                  # degenerate test (one allele with no read/no 8C)
    denom = comb(n, c1)
    lo, hi = max(0, c1 - r2), min(r1, c1)           # admissible a'
    return sum(comb(r1, x) * comb(r2, c1 - x) for x in range(a, hi + 1)) / denom


def _evaluate(scans_by_ctx: dict, *, hp_mut=None, pon=None, prior: float = 0.3,
              d_assumed: float = 0.3, alpha: float = 0.05, default_f: float = 0.02,
              min_n8c: int = 4, paired_alpha: float = 0.05, require_paired: bool = True,
              drop_artifact_ctx: bool = True, min_report: int = 2,
              borderline_min: float = 0.9, mut_len: int = 8) -> dict:
    """Pure scoring (testable): dict {ctx_units: scan_dupc_context output} → dupC call + tier.

    A candidate is CALLED (`tier="called"`) if it passes three gates (calibrated on spec. 0.66→1.0,
    cf. `docs/dupc_detection.md`): (a) absolute support `n8C >= min_n8c` (a real dupC gives many:
    e.g. 10; `2/2`=noise); (b) PAIRED enrichment mut-HP vs healthy-HP same context (Fisher,
    `paired_alpha`); (c) non-artifact context (`58_59delCC`/`?`). Then test vs PoN + Bonferroni.

    A candidate below the support threshold but with a **Bayesian posterior** ≥ `borderline_min` (hotspot
    prior × likelihood) comes out as `tier="borderline"` (e.g. one sample: 3 8C, high posterior but
    PRIOR-dominated) — so as not to hide a real-but-below-threshold signal behind a silent `False`.
    ⚠ a borderline's posterior is driven by the prior (few reads): read `n8C`/`tot` with it.
    `hp_mut` None → test both HP (discovery / negative control).
    """
    cand = []                                       # all non-artifact candidates, support >= min_report
    for k, res in scans_by_ctx.items():
        for c in res.get("dupC_candidates", []):
            h = c["hp_8C"]
            if hp_mut in ("1", "2") and h != hp_mut:
                continue
            if drop_artifact_ctx and _is_artifact_ctx(c["context"], mut_len):
                continue                            # gate (c): decoding artifact (target-aware)
            other = "2" if h == "1" else "1"
            n8, n7 = c["n8C"][h], c["n7C"][h]
            n8o, n7o = c["n8C"][other], c["n7C"][other]
            tot = n8 + n7
            if tot == 0 or n8 < min_report:          # too few to even report
                continue
            f = f_for_context(pon, c["context"], default=default_f, target=mut_len) if pon else default_f
            cand.append({"ctx_units": k, "context": c["context"], "hp": h, "n8C": n8, "tot": tot,
                         "n8C_healthy": n8o, "tot_healthy": n8o + n7o,
                         "frac_8C": round(n8 / tot, 3), "f": f, "p": _tail(tot, n8, f),
                         "p_paired": _fisher_right(n8, n7, n8o, n7o),
                         "posterior": _posterior(prior, n8, tot, d=d_assumed, f=f)})
    # CALL gates: absolute support -> Bonferroni over these tests only -> paired
    callable_ = [t for t in cand if t["n8C"] >= min_n8c]
    n_tests = len(callable_)
    for t in callable_:
        t["p_bonf"] = min(1.0, t["p"] * n_tests)
        t["pass_paired"] = (not require_paired) or (t["p_paired"] <= paired_alpha)
    callable_.sort(key=lambda t: (t["p"], -t["n8C"]))
    best = callable_[0] if callable_ else None
    called = bool(best and best["p_bonf"] <= alpha and best["pass_paired"])
    # borderline: best NON-called candidate by posterior (suggestive signal, prior-dominated)
    pool = [t for t in cand if not (called and t is best)]
    borderline = max(pool, key=lambda t: (t["posterior"], t["n8C"]), default=None)
    if called:
        tier = "called"
    elif borderline and borderline["posterior"] >= borderline_min:
        tier = "borderline"
    else:
        tier = "none"
    cand.sort(key=lambda t: (-t["posterior"], t["p"]))
    return {"hp_mut": hp_mut, "mut_len": mut_len, "n_tests": n_tests, "tier": tier,
            "called": called, "best": best,
            "borderline": borderline if tier == "borderline" else None,
            "candidates": cand[:8]}


def _region_has_hp(bam: str, region: str = None, genome_ref: str = None, n: int = 500) -> bool:
    """True if ≥1 read in the region carries an HP tag (PHASED bam). Used to decide whether the
    context caller (HP-dependent) is usable, or whether to route to the SEQUENCE-anchored positional
    path (unphased / off-hg38). Any error (missing contig, invalid region) → False → route seq-anchor (safe)."""
    import pysam
    from .config import GRCh38
    if region is None:
        loc = GRCh38["LOCUS"]
        region = f"{loc.chrom}:{loc.start}-{loc.end}"
    chrom, coords = region.split(":")
    start, end = (int(x) for x in coords.replace(",", "").split("-"))
    mode = "rc" if str(bam).endswith(".cram") else "rb"
    kw = {"reference_filename": genome_ref} if (mode == "rc" and genome_ref) else {}
    try:
        with pysam.AlignmentFile(bam, mode, **kw) as af:
            for i, r in enumerate(af.fetch(chrom, start, end)):
                if r.has_tag("HP"):
                    return True
                if i >= n:
                    break
    except Exception:
        return False
    return False


def call_dupc(bam: str, *, hp_mut=None, pon=None, ctx_range=(3, 4, 5), genome_ref: str = None,
              region: str = None, min_support: int = 2, prior: float = 0.3,
              d_assumed: float = 0.3, alpha: float = 0.05, default_f: float = 0.02,
              min_n8c: int = 4, paired_alpha: float = 0.05, require_paired: bool = True,
              drop_artifact_ctx: bool = True, borderline_min: float = 0.9,
              mut_len: int = 8, wt_len: int = 7, seq_anchor: bool = False,
              probe: str = "ctract") -> dict:
    """Call a C-tract variant by scanning `ctx_range` + PoN null + Bonferroni + posterior.

    `mut_len` = mutant C-tract length (8 = dupC, 5 = delCC), `wt_len` = WT (7). The probe reads the
    LITERAL C-tract (not `match_motifs`). Specificity gates + borderline tier: see `_evaluate`.
    `probe="unitlen"`: switches to the UNIT-LENGTH probe (del8_27 / body deletions), sequence-anchored
    (unphased/T2T/urine OK) — a separate path, does not touch the C-tract.
    """
    if probe == "unitlen":
        from .detectors.vntr_dupc import scan_del_positional_seq
        pos = scan_del_positional_seq(bam, region=region, genome_ref=genome_ref)
        return {"tier": "positional" if pos.get("called") else "none", "called": pos.get("called"),
                "best": pos.get("best"), "probe": "unitlen",
                "positional": {"called": pos.get("called"), "best": pos.get("best"),
                               "n_span": pos.get("n_span"), "anchor": pos.get("anchor"),
                               "frame": pos.get("frame"), "frames": pos.get("frames")}}
    # `--seq-anchor` = we explicitly want the sequence-anchored positional path → SKIP the context scan
    # (the per-read `match_motifs`, costly: timeout on login for a large CRAM). Otherwise normal scan.
    if seq_anchor:
        out = {"tier": "none", "called": False, "ctx_range": list(ctx_range), "context_skipped": True}
    else:
        from .detectors.vntr_dupc import scan_dupc_context
        scans = {k: scan_dupc_context(bam, region=region, genome_ref=genome_ref, hp_mut=hp_mut,
                                      ctx_units=k, min_support=min_support, mut_len=mut_len,
                                      wt_len=wt_len) for k in ctx_range}
        out = _evaluate(scans, hp_mut=hp_mut, pon=pon, prior=prior, d_assumed=d_assumed,
                        alpha=alpha, default_f=default_f, min_n8c=min_n8c, paired_alpha=paired_alpha,
                        require_paired=require_paired, drop_artifact_ctx=drop_artifact_ctx,
                        min_report=min_support, borderline_min=borderline_min, mut_len=mut_len)
        out["ctx_range"] = list(ctx_range)

    # HP-AGNOSTIC POSITION-ANCHORED FALLBACK (dupC only): the per-HP caller fails when whatshap phasing
    # collapses (length-homozygous: a starved HP) -> no candidate. `scan_dupc_positional` ignores the HP,
    # anchors on the flanks, tests the 8C by POSITION (binomial + Bonferroni). Validated on the cohort:
    # recovers homozygous/het carriers, 0 FP (negatives rejected). We only enable it if the
    # context called NOTHING, with a distinct `positional` tier (transparent, less context).
    # `mut_len in (8, 5)`: the sequence-anchored positional works for any C-tract (8=dupC, 5=delCC) — the
    # delCC re-opens the anchor to non-dupC (⚠ stronger ONT −2 noise → watch specificity on negatives).
    if mut_len in (8, 5) and out.get("tier") != "called":
        try:
            # AUTO-ROUTING (resolves "hg38-native + HP-dependent"): SEQUENCE-flank anchoring (T2T-native,
            # no HP nor hg38 anchors) if requested OR if the bam is not phased (no HP tag) → an unphased /
            # off-hg38 bam (e.g. URINE on 'muc1win') yields a real result instead of None. Otherwise
            # (phased hg38, mute context like a homozygous-length sample) → historical hg38 genomic anchoring.
            use_seq = seq_anchor or not _region_has_hp(bam, region, genome_ref)
            if use_seq:
                from .detectors.vntr_dupc import scan_dupc_positional_seq
                pos = scan_dupc_positional_seq(bam, region=region, genome_ref=genome_ref, mut_len=mut_len)
            else:
                from .detectors.vntr_dupc import scan_dupc_positional
                pos = scan_dupc_positional(bam, genome_ref=genome_ref, mut_len=mut_len)
            out["positional"] = {"called": pos.get("called"), "best": pos.get("best"),
                                 "n_span": pos.get("n_span"), "anchor": pos.get("anchor", "genomic"),
                                 "frame": pos.get("frame"), "frames": pos.get("frames")}
            if pos.get("called"):
                out["tier"] = "positional"
                out["called"] = True
                out["best"] = out.get("best") or {"context": "positional", "hp": None,
                                                  "posterior": None, **(pos.get("best") or {})}
        except Exception as e:                      # a fallback failure does not lose the context result
            out["positional"] = {"error": str(e)}
    return out


if __name__ == "__main__":
    import argparse
    import json
    ap = argparse.ArgumentParser(prog="muc1_analyzer.dupc_caller",
                                 description="Multi-resolution dupC caller (auto de-dilution + PoN)")
    ap.add_argument("-b", "--bam", required=True)
    ap.add_argument("--genome-ref", default=None)
    ap.add_argument("--region", default=None,
                    help="region to scan (default = GRCh38 chr1 dupC locus). For a BAM aligned on "
                         "the T2T 'muc1win' contig: --region muc1win:2328605-2330538")
    ap.add_argument("--hp-mut", choices=["1", "2"], default=None,
                    help="expected mutant HP; omitted = test both HP (discovery / negative control)")
    ap.add_argument("--pon", default=None, help="PoN table (pon_dupc.json) for the per-context null")
    ap.add_argument("--ctx", default="3,4,5", help="context resolutions to scan")
    ap.add_argument("--prior", type=float, default=0.3, help="mutation (hotspot) prior for the posterior")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--min-n8c", type=int, default=4,
                    help="minimum absolute 8C on the mutant HP to CALL (specificity gate; 2/2 = noise)")
    ap.add_argument("--paired-alpha", type=float, default=0.05,
                    help="threshold of the paired mut-HP vs healthy-HP test (Fisher)")
    ap.add_argument("--no-paired", action="store_true",
                    help="disable the paired intra-patient null gate (not recommended)")
    ap.add_argument("--mut-len", type=int, default=8,
                    help="mutant C-tract length: 8 = dupC (default), 5 = delCC")
    ap.add_argument("--wt-len", type=int, default=7, help="WT C-tract length (default 7)")
    ap.add_argument("--min-support", type=int, default=2,
                    help="minimum mut C-tracts to keep a context candidate (scan)")
    ap.add_argument("--no-artifact-mask", action="store_true",
                    help="do NOT mask the `58_59delCC`/`?` contexts (⚠ for delCC, the signal IS there)")
    ap.add_argument("--seq-anchor", action="store_true",
                    help="SEQUENCE-anchored positional fallback (alignment-free, no HP nor hg38 anchors) "
                         "— for an UNPHASED bam or one aligned off-hg38 (e.g. urine on 'muc1win')")
    ap.add_argument("--probe", choices=["ctract", "unitlen"], default="ctract",
                    help="ctract = C-tract (dupC/delCC, default); unitlen = UNIT LENGTH (del8_27 / "
                         "motif-body deletions), sequence-anchored")
    args = ap.parse_args()
    pon = json.load(open(args.pon)) if args.pon else None
    res = call_dupc(args.bam, hp_mut=args.hp_mut, pon=pon,
                    ctx_range=tuple(int(x) for x in args.ctx.split(",")),
                    genome_ref=args.genome_ref, region=args.region, prior=args.prior, alpha=args.alpha,
                    min_n8c=args.min_n8c, paired_alpha=args.paired_alpha, min_support=args.min_support,
                    require_paired=not args.no_paired, mut_len=args.mut_len, wt_len=args.wt_len,
                    drop_artifact_ctx=not args.no_artifact_mask, seq_anchor=args.seq_anchor,
                    probe=args.probe)
    print(json.dumps(res, indent=2, ensure_ascii=False))
