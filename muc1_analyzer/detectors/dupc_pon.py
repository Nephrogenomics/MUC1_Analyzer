"""PoN of the ONT homopolymer error at the MUC1 dupC — distribution of C-tracts BY CONTEXT, over a
cohort of non-carriers (non-MUC1 CKD / healthy controls). Alignment-free (reuses `vntr_dupc`).

Delivers, per context (local topology = types of the neighboring units):
  - the false `8C` rate (= `f`, the null of `dupc_power`), robust because aggregated over N patients;
  - a DIFFICULTY MAP by topology (fraction of C-tracts != 7C = ONT noise of the context);
  - the basis for a pre-test probability / a population LoD.

Two modes:
  1. PER-PATIENT profile (run on each CRAM, on the cluster) — emits a JSON `{context: {ctract_len: n}}`;
  2. AGGREGATION of patient profiles -> PoN table by context (f, inter-patient mean/sd/max,
     difficulty). No patient identifier retained (aggregate).
Warning: CKD cohort may contain undiagnosed MUC1 carriers -> for an ERROR model this is negligible
(a few carriers/100 do not shift the distribution), but cross-check against healthy controls.
"""
from __future__ import annotations
import collections
import json
import sys

import pysam

from ..config import GRCh38
from .vntr_dupc import _CTX, gene_oriented


def _base(name: str) -> str:
    return name.split("-")[0] if name else "?"


def context_ctract_profile(bam: str, *, region: str = None, genome_ref: str = None,
                           ctx_units: int = 3) -> dict:
    """{context: Counter{ctract_len: n}} POOLED over all primary reads (unphased)."""
    from .vntr_scaffold import _import_analyzer
    M = _import_analyzer()
    loc = GRCh38["LOCUS"]
    region = region or f"{loc.chrom}:{loc.start}-{loc.end}"
    chrom, coords = region.split(":")
    start, end = (int(x) for x in coords.replace(",", "").split("-"))
    mode = "rc" if bam.endswith(".cram") else "rb"
    kw = {"reference_filename": genome_ref} if (mode == "rc" and genome_ref) else {}

    prof = collections.defaultdict(collections.Counter)
    n_reads = 0
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        for rd in af.fetch(chrom, start, end):
            if rd.is_secondary or rd.is_supplementary or rd.is_unmapped or not rd.query_sequence:
                continue
            n_reads += 1
            g = gene_oriented(rd.query_sequence)
            for m in _CTX.finditer(g):
                clen = len(m.group(1))
                win = g[max(0, m.start() - 60 * ctx_units): m.start() + 11]   # preceding units
                matched, _ = M.match_motifs(win)
                ctx = "-".join(_base(x["name"]) for x in matched[-ctx_units:]) if matched else "?"
                prof[ctx][clen] += 1
    return {"n_reads": n_reads, "n_contexts": len(prof),
            "profile": {k: dict(v) for k, v in prof.items()}}


def aggregate(profiles: list) -> dict:
    """Pool of patient profiles -> PoN by context (f, inter-patient stats, difficulty)."""
    pooled = collections.defaultdict(collections.Counter)
    per_ctx_f8 = collections.defaultdict(list)          # context -> [frac_8C per patient]
    for p in profiles:
        for ctx, hist in p.get("profile", {}).items():
            h = collections.Counter({int(k): v for k, v in hist.items()})
            pooled[ctx].update(h)
            tot = sum(h.values())
            if tot:
                per_ctx_f8[ctx].append(h.get(8, 0) / tot)

    out = {}
    for ctx, h in pooled.items():
        tot = sum(h.values())
        fr = per_ctx_f8[ctx]
        mean = sum(fr) / len(fr) if fr else 0.0
        sd = (sum((x - mean) ** 2 for x in fr) / len(fr)) ** 0.5 if len(fr) > 1 else 0.0
        out[ctx] = {
            "n_patients": len(fr), "n_ctracts": tot,
            "pooled_7C": h.get(7, 0), "pooled_8C": h.get(8, 0),
            "counts": {str(k): v for k, v in sorted(h.items())},   # FULL C-tract histogram (5C/6C/8C...)
            "f_8C_pooled": round(h.get(8, 0) / tot, 4) if tot else None,
            "f_8C_mean": round(mean, 4), "f_8C_sd": round(sd, 4),
            "f_8C_max": round(max(fr), 4) if fr else None,
            "difficulty": round(1 - h.get(7, 0) / tot, 4) if tot else None,   # frac != 7C = context noise
        }
    # frequent contexts first (best null); the most 'difficult' can be surfaced by external sorting
    return {"n_profiles": len(profiles),
            "contexts": dict(sorted(out.items(), key=lambda kv: -kv[1]["n_ctracts"]))}


def f_for_context(pon: dict, context: str, default: float = 0.02, min_ctracts: int = 30,
                  target: int = 8) -> float:
    """`f` (null of a false C-tract of length `target`) for a context, from the aggregated PoN table.
    `default` if absent or too little data (< `min_ctracts`). Pass to `dupc_power`/`dupc_prior.posterior`/
    `dupc_caller` for a calibrated null. `target=8` (dupC) -> ~0.012 ; `target=5` (delCC) -> f of 5C per
    context. Requires a PoN re-aggregated with `counts` (full histogram); falls back to `f_8C_pooled` for target 8."""
    c = (pon.get("contexts") or {}).get(context)
    if not c or c.get("n_ctracts", 0) < min_ctracts:
        return default
    counts = c.get("counts")
    if counts is not None:
        tot = c.get("n_ctracts") or sum(counts.values())
        return round(counts.get(str(target), 0) / tot, 4) if tot else default
    if target == 8 and c.get("f_8C_pooled") is not None:   # compat old PoN (without `counts`)
        return c["f_8C_pooled"]
    return default


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(prog="muc1_analyzer.detectors.dupc_pon",
                                 description="dupC homopolymer-error PoN by context (cohort)")
    ap.add_argument("-b", "--bam", help="profile ONE patient (emits the profile JSON)")
    ap.add_argument("--genome-ref", default=None)
    ap.add_argument("--aggregate", nargs="+", metavar="PROFILE.json",
                    help="aggregate patient profiles (JSON) -> PoN table")
    args = ap.parse_args()
    if args.aggregate:
        profs, n_skipped = [], 0
        for f in args.aggregate:
            try:
                with open(f) as fh:
                    d = json.load(fh)
            except (json.JSONDecodeError, OSError):
                d = None
            if d and d.get("profile"):
                profs.append(d)
            else:
                n_skipped += 1                      # empty (0 read / module failure) or unreadable
        res = aggregate(profs)
        res["n_skipped"] = n_skipped
        print(json.dumps(res, indent=2, ensure_ascii=False))
        if n_skipped:
            print(f"# {n_skipped} profiles skipped (empty/unreadable)", file=sys.stderr)
    elif args.bam:
        print(json.dumps(context_ctract_profile(args.bam, genome_ref=args.genome_ref),
                         indent=2, ensure_ascii=False))
    else:
        ap.error("provide --bam (patient profile) or --aggregate (cohort)")
