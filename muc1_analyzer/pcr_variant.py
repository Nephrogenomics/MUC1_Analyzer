#!/usr/bin/env python3
"""pcr_variant — GENERIC per-allele frameshift caller for HIGH-DEPTH PCR/amplicon MUC1 data.

Generalises `pcr_dupc` beyond the C-tract. The whole amplicon scaffolding — split alleles by length, scan
PER allele (de-diluted from the PCR amplification skew), call by EFFECT SIZE with the internal null +
coverage/focality gate + per-strand support — is **variant-agnostic**; only the per-unit FEATURE was
dupC-specific. Here the feature is `vntr_segment.classify_segment` (reads the LITERAL inter-anchor segment,
no decode artifact) → it types every array unit as WT / dupC / delCC / delC / dupCCCC / dupG / **insG** from
the KNOWN_REPEATS dictionary. So the same machinery calls ANY dictionary variant, not just dupC.

Scope (honest, per docs/clef_generique_nonDupC.md + decisions):
  * insG / dupG — a G INSERTED in/around the C-run is a DISTINCT sequence event (not an ONT homopolymer
    slip) → expected to call as well as, or better than, dupC at amplicon depth. Ground truth: the dupG and insG carriers.
  * dupC / delCC — also typed here, but delCC's ceiling is PHYSICAL (the ONT 7C→5C −2 under-call is common
    and read-level-indistinguishable, caps ~0.85); the dedicated C-tract `pcr_dupc` stays the reference for
    dupC. Ground truth for delCC: the delCC carrier.
  * del8_27 — `--variant del8_27` switches the feature to UNIT LENGTH (`_scan_unitlen`). ⚠ UNVALIDATED: on the
    the del8_27 ground truth it does NOT fire — del8_27 produces no short C-tract gap (its unit's C-tract is
    missed → the gap SPANS the unit ~120 bp, not ~40) → del8_27 stays length-caller + visual-review. Kept for
    other body-deletions that DO shorten the gap, but do not rely on it for del8_27.

Reuses the validated decision core (`pcr_dupc._decide`) and length split (`pcr_dupc.split_by_length`) so the
specificity behaviour (internal null, focality gate, strand, four-way status) is IDENTICAL to the dupC mode.

  python3 -m muc1_analyzer.pcr_variant -b sample.chr1.bam --region chr1:155188000-155192000 --variant insG
"""
from __future__ import annotations
import argparse
import collections
import json
from typing import Optional

from .pcr_dupc import split_by_length, _decide

VARIANTS = ("insG", "dupG", "dupC", "delCC", "delC", "dupCCCC", "del8_27")
UNITLEN_VARIANTS = ("del8_27",)   # body deletions: the feature is UNIT LENGTH (≤ short_max), not the C-tract segment


# ── pure cores (pysam-free, unit-tested) ──────────────────────────────────────

def _best_variant_index(idxmap: dict, target: str, min_tot: int, max_other_frac: float = 0.5):
    """Best array index for `target` variant type: {index,n_var,tot,frac} over indices covered ≥ min_tot.
    BOUNDARY GUARD — skip indices dominated by `other` (≥ max_other_frac): the array-START primer region
    (idx≈3, and sometimes idx5/6) is universally `other`-saturated (0.5-0.8 across ALL samples) because the
    inter-anchor segment there is not a clean 60 bp VNTR unit → segment typing is junk. A clean unit index is
    mostly WT with a MINORITY variant (a real dupG at idx4 has other=0.06; a delCC at idx41 other=0.21)."""
    best = None
    for i, c in idxmap.items():
        tot = sum(c.values())
        if tot >= min_tot and c.get("other", 0) / tot < max_other_frac:
            fr = c.get(target, 0) / tot
            if best is None or fr > best["frac"]:
                best = {"index": i, "n_var": c.get(target, 0), "tot": tot, "frac": round(fr, 3)}
    return best


def _idx_variant(reads, target: str, min_tot: int, max_other_frac: float = 0.5):
    """From reads = [(n_units, strand, [seg_type per unit]),…] build {index: Counter{type}} (both strands)
    → best index for `target`, plus that index's per-strand target fraction. Pure."""
    both = collections.defaultdict(collections.Counter)
    per_strand = {"fwd": collections.defaultdict(collections.Counter),
                  "rev": collections.defaultdict(collections.Counter)}
    for _n, s, types in reads:
        for i, t in enumerate(types):
            both[i][t] += 1
            if s in per_strand:
                per_strand[s][i][t] += 1
    best = _best_variant_index(both, target, min_tot, max_other_frac)
    strand_frac = None
    if best:
        idx = best["index"]
        strand_frac = {}
        for s in ("fwd", "rev"):
            c = per_strand[s].get(idx, collections.Counter())
            tot = sum(c.values())
            strand_frac[s] = round(c.get(target, 0) / tot, 3) if tot else None
    return best, strand_frac


def call_variant_from_reads(reads, *, variant: str = "insG", floor: float = 0.10, min_abs: float = 0.07,
                            min_ratio: float = 1.4, min_tot: int = 20, min_sep: int = 5,
                            min_strand_frac: float = 0.04, min_cand: int = 5,
                            min_cov_frac: float = 0.03, max_other_frac: float = 0.5) -> dict:
    """Full per-allele call for ANY dictionary `variant` from segment-typed reads (no pysam) — testable.
    reads = [(n_units, strand, [seg_type per unit]), …]. Same internal-null / focality / strand / status
    machinery as `pcr_dupc.call_from_reads`; only the per-index feature is the variant type, not the 8C count."""
    split = split_by_length([n for n, _, _ in reads], min_sep=min_sep)
    if split["het"]:
        thr = split["threshold"]
        bins = {"short": [r for r in reads if r[0] <= thr], "long": [r for r in reads if r[0] > thr]}
    else:
        bins = {"single": reads}
    per = {}
    for name, rs in bins.items():
        eff_min_tot = max(min_tot, int(min_cov_frac * len(rs)))   # coverage/focality gate (see pcr_dupc)
        best, sf = _idx_variant(rs, variant, eff_min_tot, max_other_frac)
        ok_cov = bool(best and best.get("tot", 0) >= eff_min_tot)
        frac = best["frac"] if (best and ok_cov) else None
        cov = "ok" if ok_cov else ("low" if len(rs) >= min_cand else "none")
        strand_ok = True
        if frac is not None:
            vals = [v for v in (sf or {}).values() if v is not None]
            strand_ok = len(vals) >= 2 and min(vals) >= min_strand_frac
        per[name] = {"n_reads": len(rs), "best": best, "frac": frac, "min_tot_eff": eff_min_tot,
                     "strand_frac": sf, "strand_ok": strand_ok, "cov": cov}
    out = {"n_reads": len(reads), "length_split": split, "variant": variant, "alleles": per}
    out.update(_decide(per, floor=floor, min_abs=min_abs, min_ratio=min_ratio))
    return out


# ── pysam reader ──────────────────────────────────────────────────────────────

def _scan_segments(bam: str, chrom: str, start: int, end: int, ref: Optional[str] = None):
    """Per read spanning [start,end): (n_units, strand, [segment TYPE per unit]) via classify_segment."""
    import pysam
    from .detectors.vntr_dupc import gene_oriented
    from .detectors.vntr_segment import _SEG, classify_segment
    mode = "rc" if str(bam).endswith(".cram") else "rb"
    kw = {"reference_filename": ref} if (mode == "rc" and ref) else {}
    out = []
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        for r in af.fetch(chrom, start, end):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or r.query_sequence is None:
                continue
            seq = gene_oriented(r.query_sequence.upper())
            types = [classify_segment(m.group(1))["type"] for m in _SEG.finditer(seq)]
            if not types:
                continue
            out.append((len(types), "rev" if r.is_reverse else "fwd", types))
    return out


def _scan_unitlen(bam: str, chrom: str, start: int, end: int, ref: Optional[str] = None,
                  short_max: int = 50):
    """Per read: (n_units, strand, ['del8_27' if unit ≤ short_max else 'WT', …]). Body deletions (del8_27)
    leave the terminal C-tract anchor intact but shorten the unit ~60→~40 bp → the feature is the
    UNIT LENGTH, not the C-tract composition. gene5-anchored for a consistent per-unit index across reads."""
    import pysam
    from .detectors.vntr_dupc import unitlens_from_gene5
    mode = "rc" if str(bam).endswith(".cram") else "rb"
    kw = {"reference_filename": ref} if (mode == "rc" and ref) else {}
    out = []
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        for r in af.fetch(chrom, start, end):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or r.query_sequence is None:
                continue
            gaps = unitlens_from_gene5(r.query_sequence.upper())
            if not gaps:
                continue
            classes = ["del8_27" if g <= short_max else "WT" for g in gaps]
            out.append((len(gaps), "rev" if r.is_reverse else "fwd", classes))
    return out


def pcr_variant(bam: str, *, chrom: str, start: int, end: int, variant: str = "insG",
                ref: Optional[str] = None, floor: float = 0.10, min_tot: int = 20, min_sep: int = 5,
                min_cov_frac: float = 0.03, max_other_frac: float = 0.5) -> dict:
    """Per-allele generic-variant call on a deep amplicon BAM. Feature = segment composition for C-tract
    variants; UNIT LENGTH for body deletions (`del8_27`). Same per-allele decision core either way."""
    reads = (_scan_unitlen(bam, chrom, start, end, ref) if variant in UNITLEN_VARIANTS
             else _scan_segments(bam, chrom, start, end, ref))
    res = call_variant_from_reads(reads, variant=variant, floor=floor, min_tot=min_tot, min_sep=min_sep,
                                  min_cov_frac=min_cov_frac, max_other_frac=max_other_frac)
    return {"region": f"{chrom}:{start}-{end}", **res}


def _parse_region(s: str):
    chrom, rng = s.split(":")
    a, b = rng.replace(",", "").split("-")
    return chrom, int(a), int(b)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="muc1_analyzer.pcr_variant",
                                 description="Generic per-allele frameshift caller for high-depth PCR/amplicon MUC1")
    ap.add_argument("-b", "--bam", required=True)
    ap.add_argument("--region", required=True, help="e.g. chr1:155188000-155192000 or muc1win:2328605-2330538")
    ap.add_argument("--variant", default="insG", choices=VARIANTS, help="dictionary variant type to call")
    ap.add_argument("--ref", default=None, help="reference FASTA (required for a CRAM)")
    ap.add_argument("--floor", type=float, default=0.10)
    ap.add_argument("--min-tot", type=int, default=20)
    ap.add_argument("--min-sep", type=int, default=5)
    ap.add_argument("--min-cov-frac", type=float, default=0.03, help="coverage/focality gate (0 = off)")
    ap.add_argument("--max-other-frac", type=float, default=0.5,
                    help="boundary guard: skip array indices where 'other' segments dominate (1 = off)")
    a = ap.parse_args(argv)
    chrom, start, end = _parse_region(a.region)
    out = pcr_variant(a.bam, chrom=chrom, start=start, end=end, variant=a.variant, ref=a.ref,
                      floor=a.floor, min_tot=a.min_tot, min_sep=a.min_sep, min_cov_frac=a.min_cov_frac,
                      max_other_frac=a.max_other_frac)
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
