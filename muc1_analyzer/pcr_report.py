#!/usr/bin/env python3
"""pcr_report — clinical-coordinate REPORTING layer for confirmed MUC1 amplicon frameshifts.

The fast per-allele callers (`pcr_dupc` / `pcr_variant`) DETECT a frameshift (spec-validated) but report
INTERNAL coordinates — a C-tract segment-index and a C-tract-count length bin — which are (a) not clinical and
(b) an unreliable allele proxy (the C-tract count does not track array length → a 52dupG sample landed on the
wrong length bin). This layer runs on a CONFIRMED positive and produces the clinical report by FULL per-read
motif decomposition (`MUC1_Analyzer.match_motifs` → `detectors.vntr.frameshift_variants`), which yields:
  - the **motif name** already in the community `pristanna/muc1repeats` convention (`X-52dupG`, `X-58_59insG`,
    `X-59dupC`, `X-8_27del…` — all present in `KNOWN_REPEATS`), rendered as Kmoch `X:52dupG`;
  - the **array-repeat position** (`repeat_index`, 1-based over ALL units, not just C-tract segments);
  - the **carrier allele length** = the TRUE decomposition length of the reads carrying the variant (fixes the
    segment-count mis-assignment: that sample's carriers cluster at ~70 units, not the '80' the segment split implied).

    → reported as `X:52dupG @ repeat 11 on the ~70-unit allele`, matching the clinical ground truth.

This complements (does not replace) the fast detectors: DETECT with pcr_dupc/pcr_variant (spec-validated),
then REPORT here. del8_27 (`X-8_27del…`) is in the dict, but the greedy matcher tends to absorb a body
deletion as WT → it stays a length/visual-review case; this layer surfaces it only when a read decodes it.

⚠ pysam + MUC1_Analyzer (compute; match_motifs is per-read). The aggregation core is pure + unit-tested.

  python3 -m muc1_analyzer.pcr_report -b sample.chr1.bam --region chr1:155188000-155192000
"""
from __future__ import annotations
import argparse
import collections
import json
import statistics as st
from typing import Optional


def _kmoch(name: str) -> str:
    """`X-52dupG` → `X:52dupG` (pristanna/muc1repeats community format: `<UnitLetter>:<pos><type><bases>`)."""
    return name.replace("-", ":", 1) if name and "-" in name else name


# ── pure aggregation core (pysam-free, unit-tested) ───────────────────────────

# Fraction of carrier reads that must agree on the SAME index for the repeat number to be reported as
# exact. Below it the reads are truncated and the mode is a starting offset, not a position.
MIN_REPEAT_CONCENTRATION = 0.5


def aggregate_frameshifts(records, *, min_reads: int = 5, min_strand: int = 2,
                          focal_max_units: int = 3) -> dict:
    """records = [(total_units:int, frameshifts:list[dict], strand:str), …] where each frameshift dict is a
    `detectors.vntr.frameshift_variants` item {repeat_index, motif, delta_bp, exact, known_pathogenic}.
    Aggregate per motif name → a clinical variant report. Pure.

    Per-read focal guard: a read decoding into MANY frameshift units (> focal_max_units) is decomposition
    noise (array-wide artifact), dropped. A real variant is FOCAL → many reads agree on one motif+repeat."""
    clean = [(tu, fv, s) for tu, fv, s in records if 0 < len(fv) <= focal_max_units]
    all_len = [tu for tu, _, _ in records if tu > 0]
    by_motif = collections.defaultdict(list)          # motif -> [(repeat_index, total_units, strand, exact)]
    for tu, fv, s in clean:
        for v in fv:
            by_motif[v["motif"]].append((v["repeat_index"], tu, s, bool(v.get("exact", False))))
    variants = []
    for motif, hits in by_motif.items():
        if len(hits) < min_reads:
            continue
        # The per-read index is the unit's position in THAT READ's decomposition, not in the allele, so a
        # read truncated at the 5' end shifts every index down. On full-length amplicons all reads start at
        # the same place and the mode IS the repeat number; on truncated reads the mode is just the most
        # common starting offset. Measured on a real sample: 147 carrier reads, median read 25 units of a
        # 66-unit allele, index spread over 45 distinct values, mode holding 8 % of the reads — the reported
        # "repeat 60" was noise. Ship the concentration so the repeat is never read as exact when it is not.
        idx_counts = collections.Counter(h[0] for h in hits)
        repeat_index, idx_mode_n = idx_counts.most_common(1)[0]
        repeat_conc = idx_mode_n / len(hits)
        strands = sorted({h[2] for h in hits})
        carrier_units = [h[1] for h in hits]
        variants.append({
            "motif": motif, "kmoch": _kmoch(motif), "repeat_index": repeat_index,
            "repeat_concentration": round(repeat_conc, 3),
            "repeat_index_distinct": len(idx_counts),
            "repeat_confident": bool(repeat_conc >= MIN_REPEAT_CONCENTRATION),
            "n_reads": len(hits), "n_exact": sum(1 for h in hits if h[3]),
            "strands": strands, "strand_ok": len(strands) >= min_strand,
            "allele_units": int(st.median(carrier_units)),
            "allele_units_range": [min(carrier_units), max(carrier_units)],
        })
    variants.sort(key=lambda v: -v["n_reads"])
    # rough two allele-length modes from ALL reads (context; the length caller is authoritative)
    modes = [m for m, _ in collections.Counter(round(x / 2) * 2 for x in all_len).most_common(2)] if all_len else []
    return {"n_reads": len(records), "n_clean": len(clean),
            "allele_length_modes": sorted(modes, reverse=True), "variants": variants}


def format_report(res: dict) -> str:
    """One clinical line per detected variant (+ a note when none pass)."""
    if not res.get("variants"):
        return "no focal frameshift decoded (detection is by pcr_dupc/pcr_variant; may be length/review-only)"
    lines = []
    for v in res["variants"]:
        strand = "both strands" if v["strand_ok"] else "one strand"
        rep = (f"repeat {v['repeat_index']}" if v.get("repeat_confident", True) else
               f"repeat ~{v['repeat_index']} (POSITION UNRELIABLE: only "
               f"{v.get('repeat_concentration', 0):.0%} of carrier reads agree, index spread over "
               f"{v.get('repeat_index_distinct', 0)} values — reads truncated)")
        lines.append(f"{v['kmoch']} @ {rep} on the ~{v['allele_units']}-unit allele "
                     f"({v['n_reads']} reads, {v['n_exact']} exact, {strand})")
    return " ; ".join(lines)


# ── pysam + MUC1_Analyzer reader (fetch = IO, decode = CPU → parallel) ─────────

def _decode_seq(arg):
    """Worker (top-level → picklable for multiprocessing): raw read seq → (total_units, frameshifts, strand).
    match_motifs is pure-Python CPU-bound → parallelise across PROCESSES (the GIL makes threads useless)."""
    seq, strand = arg
    from .detectors.vntr_scaffold import _import_analyzer
    from .detectors.vntr import frameshift_variants
    from .detectors.vntr_dupc import gene_oriented
    motifs, _ = _import_analyzer().match_motifs(gene_oriented(seq))
    return (len(motifs), frameshift_variants(motifs), strand)


def _decompose_reads(bam: str, chrom: str, start: int, end: int, ref: Optional[str] = None,
                     max_reads: Optional[int] = None, jobs: int = 1):
    """Fetch read sequences (pysam, IO), then decompose them (match_motifs, CPU) — parallel over `jobs`."""
    import pysam
    mode = "rc" if str(bam).endswith(".cram") else "rb"
    kw = {"reference_filename": ref} if (mode == "rc" and ref) else {}
    seqs = []
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        for r in af.fetch(chrom, start, end):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or r.query_sequence is None:
                continue
            seqs.append((r.query_sequence.upper(), "rev" if r.is_reverse else "fwd"))
            if max_reads and len(seqs) >= max_reads:
                break
    if jobs and jobs > 1 and len(seqs) > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            return list(ex.map(_decode_seq, seqs, chunksize=4))
    return [_decode_seq(s) for s in seqs]


def pcr_report(bam: str, *, chrom: str, start: int, end: int, ref: Optional[str] = None,
               min_reads: int = 5, max_reads: Optional[int] = None, jobs: int = 1) -> dict:
    """Clinical-coordinate frameshift report for a deep amplicon (full per-read decomposition + aggregation)."""
    records = _decompose_reads(bam, chrom, start, end, ref, max_reads=max_reads, jobs=jobs)
    res = aggregate_frameshifts(records, min_reads=min_reads)
    res["report"] = format_report(res)
    return {"region": f"{chrom}:{start}-{end}", **res}


def _parse_region(s: str):
    chrom, rng = s.split(":")
    a, b = rng.replace(",", "").split("-")
    return chrom, int(a), int(b)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="muc1_analyzer.pcr_report",
                                 description="Clinical-coordinate frameshift report for high-depth PCR/amplicon MUC1")
    ap.add_argument("-b", "--bam", required=True)
    ap.add_argument("--region", required=True, help="e.g. chr1:155188000-155192000")
    ap.add_argument("--ref", default=None, help="reference FASTA (required for a CRAM)")
    ap.add_argument("--min-reads", type=int, default=5, help="min carrier reads to report a variant")
    ap.add_argument("--max-reads", type=int, default=None, help="cap reads decomposed (match_motifs is per-read)")
    ap.add_argument("--jobs", type=int, default=1, help="parallel processes for match_motifs (CPU-bound)")
    a = ap.parse_args(argv)
    chrom, start, end = _parse_region(a.region)
    out = pcr_report(a.bam, chrom=chrom, start=start, end=end, ref=a.ref,
                     min_reads=a.min_reads, max_reads=a.max_reads, jobs=a.jobs)
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
