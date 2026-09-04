"""VNTR detector : length per haplotype, dupC presence, mutant vs healthy allele, ratio.

Relies on MUC1_Analyzer.py (already validated). Two intended modes :
  - parse of an existing MUC1_Analyzer JSON (fast path, implemented),
  - orchestration Bam_cleaner -> realignment -> MUC1_Analyzer (Phase 1, TODO).
"""
from __future__ import annotations
import json
import re
from typing import Optional

import numpy as np
import pysam


UNIT_BP = 60  # length of a normal VNTR unit (20 codons, in frame)

# Focal guard (thresholds in config.py, no heavy dependency) : a real ADTKD frameshift hits
# ONE unit. If a large fraction of the array carries the SAME shift -> decomposition artifact
# (ref motif mis-aligned at low coverage), not a mutation -> flag `array_wide_frameshift_artifact`.
from ..config import ARTIFACT_MIN_UNITS, FOCAL_MAX_FRACTION  # noqa: E402

# ── Two classes of frameshift motif — the distinction `docs/muc1repeats_contribution.md` already
# ── states in prose, and which the code did not know (2026-07-29).
#
# The motif dictionary registers a unit for TWO different reasons, and conflating them is what put 3 of
# 8 negatives over the carrier flag in the release gate:
#   1. ATTESTED — a variant clinically attested as ADTKD-MUC1 (ours + the literature).
#   2. DECODING FORM — "registered to avoid mis-decomposition" (the catalogue's own words). `60dupA` is
#      the type case: it decodes 100 % exactly in 9-37 reads across 4 unrelated subjects at early repeats
#      (2026-07-26), and a NEGATIVE control was reported a "60dupA severe carrier" before the
#      AS path gated it. `run` never had that gate, so the same artifact resurfaced on N4/N5 at index 10.
# A unit whose length shifts the frame is a frameshift either way; only class 1 may set the CARRIER FLAG.
# Class 2 stays visible as an alarm — suppressing it would hide a real novel variant.
_ATTESTED_RE = re.compile(r"(?:59dupC(?!C)|58_59delCC|52dupG|58_59insG|8_27del)")

# Kept as the historical name/behaviour for the dupC family alone (used for the caller's own label).
_KNOWN_PATHOGENIC_RE = re.compile(r"dupC(?!C)")


def strand_balance(n_fwd, n_rev):
    """Fraction of the consensus's reads on its MINORITY strand: 0.5 = balanced, 0.0 = single-stranded.
    `None` when strand support was not recorded (an older analyzer JSON). Pure (testable).

    A single-stranded pile-up is the signature of the decode artifacts: a negative control's false
    `60dupA` rested
    on 5 primary reads, 3 forward / 2 reverse, once the secondary alignments of a multi-contig smear were
    dropped. Reported here, not yet gating — a threshold must be priced on an independent set first."""
    if n_fwd is None or n_rev is None:
        return None
    tot = n_fwd + n_rev
    return round(min(n_fwd, n_rev) / tot, 3) if tot else None


def is_attested(name: str) -> bool:
    """Is this motif a CLINICALLY ATTESTED ADTKD-MUC1 variant (vs a registered decoding form)? Pure.

    Attested in our ONT cohort: 59dupC (canonical), 58_59delCC, 52dupG, 8_27del,
    58_59insG. `dupC(?!C)` keeps `56_59dupCCCC` out — a longer C-run is not the canonical variant.
    ⚠ Deliberately NOT a "looks like a frameshift" test: that is the question the carrier flag was asking,
    and it is the wrong one."""
    return bool(_ATTESTED_RE.search(str(name or "")))

# ── Auto-resolution of genome + contig for the MUC1 fetch window ──────────────────────────────────
# The chr1 LENGTH distinguishes the build (the NAME is not enough : hg38 and T2T both use
# "chr1"). MUC1 is at a DIFFERENT locus on each build -> we choose the contig name AND the window.
#   GRCh38 chr1 = 248,956,422 · GRCh37 = 249,250,621 · T2T-CHM13 = 248,387,328
_CHR1_LEN = {248_956_422: "hg38", 249_250_621: "hg38", 248_387_328: "t2t"}
_MUC1_WINDOW = {"hg38": (155_185_000, 155_192_000),      # = config.GRCh38["LOCUS"]
                "t2t":  (154_324_904, 154_335_104)}       # NC_060925.1 MUC1 window (− strand)
_CHR1_NAMES = ("chr1", "1", "CM000663.2", "NC_000001.11", "NC_060925.1", "CP068277.2")


def _region_from_sqs(sqs) -> Optional[str]:
    """PURE (testable) : from the @SQ pairs (name, length), choose the chr1 contig (by name, fallback
    by length) and the MUC1 window matching the build (hg38/T2T via the chr1 length).
    Returns 'contig:start-end', or None if no chr1 is identifiable."""
    lens = dict(sqs)
    contig = next((c for c in _CHR1_NAMES if c in lens), None)
    if contig is None:                                    # atypical name -> fallback on chr1 length
        contig = next((n for n, ln in sqs if ln in _CHR1_LEN), None)
    if contig is None:
        return None
    build = _CHR1_LEN.get(lens.get(contig), "hg38")       # default hg38 coords if length unknown
    s, e = _MUC1_WINDOW[build]
    return f"{contig}:{s}-{e}"


def _resolve_region(bam: str, genome_ref: str = None) -> Optional[str]:
    """Read the BAM/CRAM header and return the auto-detected MUC1 fetch window (build + contig),
    or None if unreadable (the caller falls back to the hg38 default). Robust to 'chr1' vs '1' naming and
    to an input aligned on T2T (different locus) -> avoids silently extracting 0 reads."""
    try:
        mode = "rc" if str(bam).endswith(".cram") else "rb"
        kw = {"reference_filename": genome_ref} if (mode == "rc" and genome_ref) else {}
        with pysam.AlignmentFile(bam, mode, **kw) as af:
            sqs = list(zip(af.references, af.lengths))
    except Exception:
        return None
    return _region_from_sqs(sqs)


def _contig_repeat_count(contig: str) -> Optional[int]:
    """Extract N from 'MUC1_VNTR_Nrepeats'."""
    m = re.search(r"_(\d+)repeats?", contig or "")
    return int(m.group(1)) if m else None


def frameshift_variants(motifs: list) -> list:
    """PATHOGENIC (frameshift) units of a haplotype + the NUMBER of the affected repeat.

    Biological criterion : a unit whose length differs from 60 bp by an amount NOT a multiple
    of 3 shifts the reading frame -> toxic MUC1-fs protein. Covers dupC, delCC, insG, dupA,
    dupCCCC... IN-FRAME indels (delta multiple of 3 : VNTR copy-number variation, 33_34ins18...)
    are benign and excluded. `known_pathogenic` = the variant is a reported ADTKD-MUC1 (dupC family).

    `motifs` = MUC1_Analyzer list [{name,pos,length,mismatches,...}]. `repeat_index` =
    position (1-based) of the motif in the array = number of the repeat carrying the anomaly.
    """
    out = []
    for idx, m in enumerate(motifs, 1):
        name = m.get("name") or ""
        delta = int(m.get("length", UNIT_BP)) - UNIT_BP
        if delta != 0 and delta % 3 != 0:               # frame shift
            out.append({"repeat_index": idx, "motif": name, "delta_bp": delta,
                        "exact": int(m.get("mismatches", 1)) == 0,
                        "attested": is_attested(name),
                        "known_pathogenic": bool(_KNOWN_PATHOGENIC_RE.search(name))})
    return out


def from_analyzer_json(path: str) -> dict:
    """Parse a MUC1_Analyzer JSON -> VNTR parameters + localized pathogenic variant.

    ADTKD-MUC1 = dominant -> heterozygous : ONE single allele carries the frameshift anomaly.
    We classify each allele (frameshift = pathogenic), force heterozygosity (mutant =
    allele with the most exact variants) and report the **precise variant + the number
    of the affected repeat** (generalizes beyond dupC).
    """
    with open(path) as fh:
        data = json.load(fh)
    haps = data.get("haplotypes", [])
    parsed = []
    for h in haps:
        variants = frameshift_variants(h.get("motifs", []))
        n_rep = _contig_repeat_count(h.get("contig", ""))
        n_fs = len(variants)
        # array-wide = the SAME shift recurring on a large fraction of the array -> decomposition
        # artifact, not a mutation. We call ONLY the focal frameshift.
        array_wide = bool(n_fs >= ARTIFACT_MIN_UNITS and n_rep and n_fs > FOCAL_MAX_FRACTION * n_rep)
        focal = [] if array_wide else variants
        exact = [v for v in focal if v["exact"]]
        # PROMOTABLE = what the consensus alone may turn into a clinical verdict: an ATTESTED variant,
        # decoded EXACTLY. `exact` matters because the consensus is documented to mint approximate
        # `delCC` tokens, and delCC is attested — so the class test alone would let those through.
        promotable = [v for v in focal if v["attested"] and v["exact"]]
        # `coverage` = {n_reads, n_fwd, n_rev} of the reads that actually BUILT this consensus (primary
        # only). `read_count` is the contig's raw tally and over-states it in a multi-contig smear —
        # measured on one negative: 12 apparent, 5 primary, 3 fwd / 2 rev. A promotion rule must read `coverage`.
        cov = h.get("coverage") or {}
        parsed.append({"contig": h.get("contig"), "n_repeats": n_rep,
                       "reads": h.get("read_count"), "n_approx": h.get("n_approx", 0),
                       "coverage": cov or None,
                       "n_primary": cov.get("n_reads"), "n_fwd": cov.get("n_fwd"),
                       "n_rev": cov.get("n_rev"),
                       "strand_balance": strand_balance(cov.get("n_fwd"), cov.get("n_rev")),
                       "has_mut": bool(focal), "variants": focal,
                       "has_attested": bool(promotable), "n_attested": len(promotable),
                       "n_frameshift_exact": len(exact), "n_frameshift_raw": n_fs,
                       "array_wide_artifact": array_wide, "nomenclature": h.get("nomenclature", "")})

    def _n_known(p):
        return sum(1 for v in p["variants"] if v["known_pathogenic"])

    # A frameshift that is only a registered DECODING FORM (60dupA…) stays visible in `variants` and in
    # the notes — it is an ALARM — but it no longer decides `carrier=`. Measured on the release gate:
    # N4 and N5, two unrelated negatives, both carried `60dupA` at the same repeat index 10.
    carriers = [p for p in parsed if p["has_attested"]]
    alarms = [p for p in parsed if p["has_mut"] and not p["has_attested"]]
    mutation_present = bool(carriers) if parsed else None
    # mutant = most strongly carrying allele : KNOWN variants first, then exact, then count
    mut = max(carriers, key=lambda p: (_n_known(p), p["n_frameshift_exact"], len(p["variants"]))) if carriers else None
    healthy = next((p for p in parsed if p is not mut), None)

    mut_variant = None
    if mut and mut["variants"]:
        # representative : attested first (it is what promoted this haplotype), then known, then exact
        v = sorted(mut["variants"],
                   key=lambda x: (not (x["attested"] and x["exact"]),
                                  not x["known_pathogenic"], not x["exact"]))[0]
        mut_variant = {"motif": v["motif"], "repeat_index": v["repeat_index"],
                       "delta_bp": v["delta_bp"], "exact": v["exact"],
                       "known_pathogenic": v["known_pathogenic"]}

    # ── the DEDICATED caller's verdict, which the consensus alone used to hide ──────────────────
    # `mutation_present` above comes from the CONSENSUS motifs only. The dedicated dupC caller
    # (statistically gated, PoN/in-sample null, Bonferroni) writes its verdict into the same JSON and was
    # never read here — so a CONFIRMED carrier whose consensus missed the frameshift came out
    # `carrier=False`, in the summary line AND in the PDF, with the JSON saying the opposite two levels
    # down. Measured 2026-07-29 on VNTRtools testdata3_index: dupC CONFIRMED (p_bonf 1.4e-04, 27/38 reads,
    # carrier_contig MUC1_VNTR_82repeats), reported as a non-carrier. A tool that detects and does not
    # report is worse than one that does not detect.
    # The consensus is documented depth-fragile (it also mints spurious `delCC` tokens), so it is the
    # weaker of the two witnesses — it must not be able to veto the stronger one.
    _dres = ((data.get("dupc") or {}).get("result") or {})
    dupc_called = bool(_dres.get("called"))
    mutation_source = ("consensus" if mutation_present else None)
    if dupc_called:
        mutation_source = "both" if mutation_present else "dupc_caller"
        mutation_present = True
        if mut_variant is None:                      # consensus saw nothing → take the caller's position
            _v = _dres.get("variant") or {}
            _label = _v.get("label")
            if _v.get("repeat") is not None or _label:
                mut_variant = {"motif": None, "repeat_index": _v.get("repeat"), "delta_bp": None,
                               "exact": None,
                               "known_pathogenic": bool(_KNOWN_PATHOGENIC_RE.search(str(_label or ""))),
                               "label": _label, "source": "dupc_caller"}

    return {
        "mutation_present": mutation_present,
        "mutation_source": mutation_source,          # consensus | dupc_caller | both | None — auditable
        "dupc_caller_called": dupc_called,
        "vntr_len_mut": mut["n_repeats"] if mut else None,
        "vntr_len_healthy": healthy["n_repeats"] if healthy else None,
        "mutation_variant": mut_variant,     # pathogenic motif + repeat number + known
        # An unattested frameshift is REPORTED, never silenced: it may be a novel pathogenic variant, and
        # this is the only place a reviewer would see it once it stops driving `carrier=`.
        "frameshift_alarms": [{"contig": p["contig"], "reads": p["reads"],
                               "variants": [v["motif"] for v in p["variants"]]} for p in alarms],
        "notes": {"vntr_haplotypes": parsed,
                  "array_wide_frameshift_artifact": [p["contig"] for p in parsed if p["array_wide_artifact"]],
                  "dupc_caller": {"called": dupc_called,
                                  "interpretation": _dres.get("interpretation"),
                                  "carrier_contig": _dres.get("carrier_contig")}},
    }


def _run(cmd: list, **kw):
    """subprocess.run with check + readable error message."""
    import subprocess
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise RuntimeError(f"failed: {' '.join(cmd[:3])}... (rc={r.returncode})\n{r.stderr[-800:]}")
    return r


# ── VNTR length per read via flank anchoring (port of script 11) ───────────────
# Anchors in UNIQUE sequence flanking the VNTR (GRCh38 "chr"). Robust : the array
# itself is unalignable, but the flanks are. O(reads), no re-alignment.
_VNTR_ANCHORS = dict(chrom="chr1", l_anchor=155_188_600, r_anchor=155_191_300,
                     flank_const=700, unit=60, del_size=364)


def _qpos_at_ref(read, target: int):
    """Read coordinate aligned to reference position `target` (CIGAR walk)."""
    rpos = read.reference_start
    qpos = 0
    for op, ln in (read.cigartuples or []):
        if op in (0, 7, 8):                 # M/=/X
            if rpos <= target < rpos + ln:
                return qpos + (target - rpos)
            rpos += ln
            qpos += ln
        elif op in (2, 3):                  # D/N
            if rpos <= target < rpos + ln:
                return qpos
            rpos += ln
        elif op in (1, 4):                  # I/S
            qpos += ln
        # H(5)/P(6) : nothing
    return None


def measure_vntr_lengths(bam: str, ref: str = None, *, chrom: str = None,
                         l_anchor: int = None, r_anchor: int = None,
                         flank_const: int = None, unit: int = None,
                         margin: int = 60, min_mq: int = 20, min_reads: int = 3,
                         del_genotype_hp: str = None, del_size: int = None) -> dict:
    """Measure the VNTR length per haplotype (read span between 2 unique anchors).

    Returns {hp1:{span,copies,n,iqr}, hp2:{...}, unphased:{...}, alleles:[...]}.
    copies = (span - flank_const) / unit. Optional DEL correction if the HP carrying
    the DEL is provided (`del_genotype_hp` in '1'/'2') : we add del_size to that allele.
    """
    import statistics as st
    a = _VNTR_ANCHORS
    chrom = chrom or a["chrom"]
    l_anchor = l_anchor if l_anchor is not None else a["l_anchor"]
    r_anchor = r_anchor if r_anchor is not None else a["r_anchor"]
    flank_const = flank_const if flank_const is not None else a["flank_const"]
    unit = unit if unit is not None else a["unit"]
    del_size = del_size if del_size is not None else a["del_size"]

    mode = "rc" if str(bam).endswith(".cram") else "rb"
    kw = {"reference_filename": ref} if (mode == "rc" and ref) else {}
    by = {"1": [], "2": [], "U": []}
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        for r in af.fetch(chrom, l_anchor - 200, r_anchor + 200):
            if r.is_unmapped or r.is_secondary or r.is_supplementary \
                    or r.is_duplicate or r.mapping_quality < min_mq:
                continue
            if r.reference_start > l_anchor - margin:      # must extend past the left anchor
                continue
            if (r.reference_end or 0) < r_anchor + margin:  # and the right anchor
                continue
            qL = _qpos_at_ref(r, l_anchor)
            qR = _qpos_at_ref(r, r_anchor)
            if qL is None or qR is None:
                continue
            span = qR - qL
            if span <= 0 or span > 12000:
                continue
            hp = "U"
            if r.has_tag("HP"):
                t = str(r.get_tag("HP"))
                hp = t if t in ("1", "2") else "U"
            by[hp].append(span)

    def summ(v, hp):
        if len(v) < min_reads:
            return None
        span = st.median(v)
        if del_genotype_hp == hp:       # DEL correction : add the footprint on the carrier allele
            span += del_size
        iqr = int(np.subtract(*np.percentile(v, [75, 25]))) if len(v) >= 2 else 0
        return {"span": round(span, 1), "copies": round((span - flank_const) / unit, 1),
                "n": len(v), "iqr": iqr}

    hp1, hp2, u = summ(by["1"], "1"), summ(by["2"], "2"), summ(by["U"], "U")
    alleles = [x["copies"] for x in (hp1, hp2) if x]
    return {"hp1": hp1, "hp2": hp2, "unphased": u,
            "allele_copies": sorted(alleles) if alleles else None,
            "n_reads": {"hp1": len(by["1"]), "hp2": len(by["2"]), "unphased": len(by["U"])}}


def detect(bam: str, vntr_ref: str, sample: str = "", *,
           genome_ref: str = None, region: str = None, workdir: str = None,
           preset: str = "map-ont", threads: int = 4, min_mapq: int = 0,
           keep: bool = False) -> dict:
    """VNTR length PER HAPLOTYPE via the best contig-length (Laurent's logic).

    For each haplotype (HP tag), we extract its reads from the locus, map them on the
    multi-contig VNTR reference, and the `MUC1_VNTR_Nrepeats` contig capturing the most
    reads gives the length N of that allele.

    Key optimizations :
      - PHASE PRESERVED : we split by HP tag BEFORE mapping (the fastq lost the
        phase, hence the 1-read dispersion of the 1st attempt).
      - PAF, NOT SAM : to *count* which contig wins, we don't need the base-level
        alignment (`-a`) — the DP explodes on the tandem VNTR. `minimap2 -x`
        (PAF, approximate mapping only) gives the best contig per read, without DP -> fast.
      - `--secondary=no` : one primary contig per read.

    `genome_ref` required for a CRAM. Returns the lengths per haplotype + diagnostics.
    """
    import os
    import sys
    import shutil
    import tempfile
    import subprocess
    import time
    from collections import Counter

    from ..config import GRCh38

    for tool in ("samtools", "minimap2"):
        if shutil.which(tool) is None:
            raise RuntimeError(f"required tool missing from PATH : {tool}")

    loc = GRCh38["LOCUS"]
    region = region or f"{loc.chrom}:{loc.start}-{loc.end}"
    tmp = workdir or tempfile.mkdtemp(prefix=f"muc1_vntr_{sample or 'sample'}_")
    os.makedirs(tmp, exist_ok=True)

    _t0 = [time.time()]

    def _lap(label):
        now = time.time()
        print(f"[vntr] {label}: {now - _t0[0]:.1f}s", file=sys.stderr, flush=True)
        _t0[0] = now

    # cached minimap2 index (.mmi) — built once
    mmi = vntr_ref + ".mmi"
    if not os.path.exists(mmi):
        try:
            _run(["minimap2", "-x", preset, "-d", mmi, vntr_ref])
        except Exception:
            mmi = vntr_ref
    ref_mm2 = mmi if os.path.exists(mmi) else vntr_ref

    def _map_hp(hp):
        """Extract the reads of the given HP (or all if hp=None), map them to PAF,
        return (Counter{length:reads}, n_reads, n_mapped, top)."""
        tag = os.path.join(tmp, f"hp{hp or 'all'}")
        hp_bam, hp_fq, paf = tag + ".bam", tag + ".fq", tag + ".paf"
        view = ["samtools", "view", "-b"]
        if genome_ref:
            view += ["-T", genome_ref]
        if hp in ("1", "2"):
            view += ["-d", f"HP:{hp}"]
        view += ["-o", hp_bam, bam, region]
        _run(view)
        _run(["samtools", "fastq", "-0", hp_fq, "-n", hp_bam])
        try:
            with open(hp_fq) as fh:
                n_reads = sum(1 for _ in fh) // 4
        except OSError:
            n_reads = 0
        if n_reads == 0:
            return Counter(), 0, 0, []
        with open(paf, "wb") as pafh, open(tag + ".log", "wb") as errfh:
            subprocess.run(["minimap2", "-x", preset, "--secondary=no",
                            "-t", str(threads), ref_mm2, hp_fq], stdout=pafh, stderr=errfh)
        counts = Counter()
        n_mapped = 0
        with open(paf) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 12:
                    continue
                if int(f[11]) < min_mapq:          # col 12 = mapping quality
                    continue
                n = _contig_repeat_count(f[5])     # col 6 = target contig
                if n is not None:
                    counts[n] += 1
                    n_mapped += 1
        return counts, n_reads, n_mapped, counts.most_common(3)

    def _call_allele(counts):
        """Allele length = majority contig, aggregated over a +-2 window
        (ONT indel noise disperses +-a few units)."""
        if not counts:
            return None, 0
        best_n, best_win = 0, None
        for center in counts:
            w = sum(counts[c] for c in range(center - 2, center + 3))
            if w > best_n:
                best_n, best_win = w, center
        # refined center = exact mode within the winning window
        mode = max(range(best_win - 2, best_win + 3), key=lambda c: counts.get(c, 0))
        return mode, best_n

    per_hp = {}
    for hp in ("1", "2"):
        counts, n_reads, n_mapped, top = _map_hp(hp)
        length, win_n = _call_allele(counts)
        per_hp[hp] = {"length": length, "reads_in_window": win_n, "n_reads": n_reads,
                      "n_mapped": n_mapped, "top3_lengths": top}
        _lap(f"HP{hp}: {n_reads} reads, {n_mapped} mapped -> {length} rep ({win_n} reads +-2)")

    lens = [per_hp[hp]["length"] for hp in ("1", "2") if per_hp[hp]["length"]]
    return {"vntr_len_hp1": per_hp["1"]["length"], "vntr_len_hp2": per_hp["2"]["length"],
            "allele_lengths": sorted(lens) if lens else None,
            "notes": {"per_hp": per_hp, "region": region,
                      "workdir": tmp if (keep or workdir) else None}}


def detect_analyzer(bam: str, vntr_ref: str, sample: str = "", *,
                    genome_ref: str = None, region: str = None, workdir: str = None,
                    analyzer: str = None, preset: str = "map-ont", threads: int = 4,
                    min_mq: int = 0, keep: bool = True) -> dict:
    """ORIGINAL path (slow but complete) : multi-contig SAM re-alignment ->
    MUC1_Analyzer -> length per haplotype + **dupC**.

    Key : MUC1_Analyzer is launched with **--min-mq {min_mq}** (def. 0), because against 150
    near-identical contigs minimap2 sets MAPQ 0 ; the default 20 would filter out all the reads
    (that is what broke the rewrite). Slow (chaining 150 contigs) — reserve for the
    few MUC1 patients ; for length alone and fast, use detect().
    """
    import os
    import sys
    import shutil
    import tempfile
    import subprocess
    import time

    for tool in ("samtools", "minimap2"):
        if shutil.which(tool) is None:
            raise RuntimeError(f"required tool missing from PATH : {tool}")
    from ..config import GRCh38
    if region is None:                                  # auto-detect build (hg38/T2T) + contig chr1/1
        region = _resolve_region(bam, genome_ref)
    if region is None:                                  # unreadable header -> hg38 default fallback
        loc = GRCh38["LOCUS"]
        region = f"{loc.chrom}:{loc.start}-{loc.end}"
    print(f"[vntrA] fetch region = {region}", file=sys.stderr, flush=True)
    analyzer = analyzer or os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "caller.py")
    tmp = workdir or tempfile.mkdtemp(prefix=f"muc1_vntrA_{sample or 'sample'}_")
    os.makedirs(tmp, exist_ok=True)
    _t0 = [time.time()]

    def _lap(label):
        now = time.time(); print(f"[vntrA] {label}: {now-_t0[0]:.1f}s", file=sys.stderr, flush=True); _t0[0] = now

    fq = os.path.join(tmp, "reg.fq")
    region_bam = os.path.join(tmp, "reg.bam")
    sam = os.path.join(tmp, "realigned.sam")
    realigned = os.path.join(tmp, "realigned.bam")
    out_json = os.path.join(tmp, f"{sample or 'sample'}.analyzer.json")

    view = ["samtools", "view", "-b", "-o", region_bam]
    if genome_ref:
        view += ["-T", genome_ref]
    view += [bam, region]
    _run(view)
    _run(["samtools", "fastq", "-0", fq, "-n", region_bam])
    _lap("extraction")

    mmi = vntr_ref + ".mmi"
    ref_mm2 = mmi if os.path.exists(mmi) else vntr_ref
    with open(sam, "wb") as sh, open(os.path.join(tmp, "mm2.log"), "wb") as eh:
        rc = subprocess.run(["minimap2", "-ax", preset, "-t", str(threads), ref_mm2, fq],
                            stdout=sh, stderr=eh)
    if rc.returncode != 0:
        raise RuntimeError("minimap2 failed (see mm2.log)")
    _run(["samtools", "sort", "-@", str(threads), "-o", realigned, sam])
    _run(["samtools", "index", realigned])
    _lap("minimap2+sort")

    _run([sys.executable, analyzer, "-b", realigned, "-r", vntr_ref, "-s", sample or "sample",
          "--min-mq", str(min_mq), "--threads", str(threads), "--json", out_json])
    _lap("MUC1_Analyzer (--min-mq %d)" % min_mq)

    frag = from_analyzer_json(out_json)
    frag.setdefault("notes", {})["vntr_paths"] = {"workdir": tmp, "realigned_bam": realigned,
                                                  "analyzer_json": out_json, "region": region}
    return frag


# ── Anomaly <-> DEL <-> rs4072037 phasing ──────────────────────────────────────
def phase_config(hp_mut, hp_del=None, hp_alt=None) -> dict:
    """Phase of the VNTR anomaly vs DEL and SNP from the carrier haplotypes.

    cis_mut = the variant (DEL or SNP ALT) is on the SAME haplotype as the anomaly
    (aggravating for the DEL : amplifies the mutant allele) ; trans = on the healthy allele.
    None if one of the haplotypes is not resolved.
    """
    def rel(hp_carrier):
        if hp_mut not in ("1", "2") or hp_carrier not in ("1", "2"):
            return None
        return "cis_mut" if hp_carrier == hp_mut else "trans"
    return {"del_phase": rel(hp_del), "snp_phase": rel(hp_alt)}


def resolve_hp_mut(per_hp: dict) -> tuple:
    """Specificity filter by PHASING : which HP carries the mutation ? -> (hp_mut, conflict).

    A dominant **heterozygous** ADTKD mutation is on ONE SINGLE haplotype : after splitting by
    genomic HP tag, all the mutant reads fall on one HP, the other is clean. If BOTH HP
    are carriers (`has_mut`) -> biologically incoherent (non-viable homozygote / mis-phasing /
    bilateral homopolymer error at low cov.) -> **conflict** : suspect, unphasable (hp_mut=None).
    Filter orthogonal to the focal guard (intra-allele) : this one is inter-allele.
    """
    carriers = [hp for hp, d in per_hp.items() if d.get("has_mut")]
    if len(carriers) >= 2:
        return None, True                # both HP mutant -> conflict, we do not phase
    return (carriers[0] if carriers else None), False


def _extract_realign_analyze(bam, vntr_ref, region, tmp, tag, *, genome_ref, analyzer,
                             preset, threads, min_mq, hp=None, top_n=None):
    """Extract (opt. per HP) -> re-align -> MUC1_Analyzer -> parse. (None,0) if 0 read."""
    import os
    import sys
    import subprocess
    hp_bam = os.path.join(tmp, tag + ".bam")
    fq = os.path.join(tmp, tag + ".fq")
    sam = os.path.join(tmp, tag + ".sam")
    bamr = os.path.join(tmp, tag + ".realigned.bam")
    out_json = os.path.join(tmp, tag + ".analyzer.json")
    view = ["samtools", "view", "-b"]
    if genome_ref:
        view += ["-T", genome_ref]
    if hp in ("1", "2"):
        view += ["-d", f"HP:{hp}"]
    view += ["-o", hp_bam, bam, region]
    _run(view)
    _run(["samtools", "fastq", "-0", fq, "-n", hp_bam])
    try:
        with open(fq) as f:
            n = sum(1 for _ in f) // 4
    except OSError:
        n = 0
    if n == 0:
        return None, 0
    mmi = vntr_ref + ".mmi"
    ref = mmi if os.path.exists(mmi) else vntr_ref
    with open(sam, "wb") as sh, open(os.path.join(tmp, tag + ".mm2.log"), "wb") as eh:
        rc = subprocess.run(["minimap2", "-ax", preset, "-t", str(threads), ref, fq],
                            stdout=sh, stderr=eh)
    if rc.returncode != 0:
        raise RuntimeError(f"minimap2 failed ({tag})")
    _run(["samtools", "sort", "-@", str(threads), "-o", bamr, sam])
    _run(["samtools", "index", bamr])
    cmd = [sys.executable, analyzer, "-b", bamr, "-r", vntr_ref, "-s", tag,
           "--min-mq", str(min_mq), "--threads", str(threads), "--json", out_json]
    if top_n:
        cmd += ["--top-n", str(top_n)]
    _run(cmd)
    return from_analyzer_json(out_json), n


def detect_analyzer_phased(bam, vntr_ref, sample="", *, genome_ref=None, region=None,
                           workdir=None, analyzer=None, preset="map-ont", threads=4,
                           min_mq=0, hp_del=None, hp_alt=None, keep=True) -> dict:
    """PHASED VNTR anomaly : run MUC1_Analyzer PER HAPLOTYPE (extraction by HP tag)
    to tie the anomaly to genome HP1/HP2, then phase it with the DEL (`hp_del`) and
    rs4072037 (`hp_alt`). ~2x the cost of detect_analyzer (one re-alignment per HP).
    """
    import os
    import shutil
    import tempfile
    for tool in ("samtools", "minimap2"):
        if shutil.which(tool) is None:
            raise RuntimeError(f"required tool missing from PATH : {tool}")
    from ..config import GRCh38
    if region is None:                                  # auto-detect build (hg38/T2T) + contig chr1/1
        region = _resolve_region(bam, genome_ref)
    if region is None:
        loc = GRCh38["LOCUS"]
        region = f"{loc.chrom}:{loc.start}-{loc.end}"
    analyzer = analyzer or os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "caller.py")
    tmp = workdir or tempfile.mkdtemp(prefix=f"muc1_vntrP_{sample or 'sample'}_")
    os.makedirs(tmp, exist_ok=True)

    per_hp = {}
    for hp in ("1", "2"):
        frag, n = _extract_realign_analyze(bam, vntr_ref, region, tmp, f"{sample or 's'}_hp{hp}",
                                           genome_ref=genome_ref, analyzer=analyzer, preset=preset,
                                           threads=threads, min_mq=min_mq, hp=hp, top_n=1)
        if frag is None:
            per_hp[hp] = {"n_reads": 0, "length": None, "has_mut": None, "variants": [],
                          "mutation_variant": None}
            continue
        haps = frag.get("notes", {}).get("vntr_haplotypes", [])
        top = haps[0] if haps else {}
        per_hp[hp] = {"n_reads": n, "length": top.get("n_repeats"),
                      "has_mut": top.get("has_mut"), "variants": top.get("variants", []),
                      "mutation_variant": frag.get("mutation_variant"),
                      "array_wide": bool(frag.get("notes", {}).get("array_wide_frameshift_artifact"))}

    carriers = [hp for hp, d in per_hp.items() if d.get("has_mut")]
    hp_mut, mut_hp_conflict = resolve_hp_mut(per_hp)     # inter-allele specificity filter (1 HP only)
    phases = phase_config(hp_mut, hp_del, hp_alt)
    # representative variant : if conflict (hp_mut=None), we still show the called one, for review.
    rep_hp = hp_mut or (carriers[0] if carriers else None)
    return {
        "mutation_present": bool(carriers),
        "hp_mut": hp_mut,
        "mut_hp_conflict": mut_hp_conflict,   # both HP mutant -> suspect, unphased (see resolve_hp_mut)
        "mutation_variant": per_hp[rep_hp]["mutation_variant"] if rep_hp else None,
        "vntr_len_hp1": per_hp["1"]["length"], "vntr_len_hp2": per_hp["2"]["length"],
        "del_phase": phases["del_phase"], "snp_phase": phases["snp_phase"],
        "array_wide_artifact": [hp for hp, d in per_hp.items() if d.get("array_wide")],
        "notes": {"per_hp": per_hp, "region": region, "hp_del": hp_del, "hp_alt": hp_alt,
                  "workdir": tmp if (keep or workdir) else None},
    }
