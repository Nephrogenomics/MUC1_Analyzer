"""3a — Per-allele scaffold/consensus reconstructor (rescue by DENSE consensus).

Diagnosis (measurement 2026-07-04) : ~120 primary reads in the MUC1 region, but the call only uses
~4-8/allele -> the reads **scatter** over the 150 near-identical contigs (MAPQ 0), so the consensus
of the winning contig is built on ~4 reads -> "-2 on every unit" artifact (observed in one sample) that masks the
true dupC.

Fix (no POA, no new dependency) : **bin** the reads by allele (whatshap HP tag), **re-align them on
ONE SINGLE contig** of the right length (instead of the 150) -> (a) all the reads of the allele concentrate
on that contig (dense consensus ~60 reads), (b) the ambiguity disappears so **MAPQ becomes informative
again**. Then consensus (`build_consensus_from_readset`) -> decomposition (`match_motifs`) ->
`frameshift_variants` + focal guard. Should kill the artifact + rescue the affected sample.

Reuses `MUC1_Analyzer.py` (already validated). The allele-length choice comes upstream
(`detect_analyzer_phased` : winning contig per HP); here we only **densify** the consensus.
"""
from __future__ import annotations
import os
import sys
import tempfile

import pysam

from ..config import GRCh38, ARTIFACT_MIN_UNITS, FOCAL_MAX_FRACTION
from .vntr import _run, frameshift_variants


def _import_analyzer():
    """Import the VNTR caller (`muc1_analyzer.caller`) as a module (consensus/motifs functions)."""
    from .. import caller
    return caller


def _focal_call(variants: list, n_rep: int) -> dict:
    """Apply the focal guard (identical to from_analyzer_json) to a list of variants of an allele."""
    n_fs = len(variants)
    array_wide = bool(n_fs >= ARTIFACT_MIN_UNITS and n_rep and n_fs > FOCAL_MAX_FRACTION * n_rep)
    focal = [] if array_wide else variants
    mv = None
    if focal:
        v = sorted(focal, key=lambda x: (not x["known_pathogenic"], not x["exact"]))[0]
        mv = {k: v[k] for k in ("motif", "repeat_index", "delta_bp", "exact", "known_pathogenic")}
    return {"has_mut": bool(focal), "variants": focal, "n_frameshift_raw": n_fs,
            "array_wide_artifact": array_wide, "mutation_variant": mv}


def _aligned_span(rd) -> int:
    """Length of the REFERENCE segment covered by the alignment (reference_end - reference_start).
    In re-aligned VNTR space, this is the true coverage of the contig by the read (no GRCh38 soft-clip)."""
    if rd.reference_start is None or rd.reference_end is None:
        return 0
    return rd.reference_end - rd.reference_start


def _softclip_len(rd) -> int:
    c = rd.cigartuples
    n = 0
    if c:
        if c[0][0] == 4:
            n += c[0][1]
        if c[-1][0] == 4:
            n += c[-1][1]
    return n


def _read_passes(rd, min_span: int, max_softclip) -> bool:
    """Span/soft-clip filter IN VNTR SPACE (Bam_cleaner on OUTPUT, pre-consensus) : keeps the reads
    that cover >= `min_span` bp of the contig and soft-clip <= `max_softclip` (None = no bound)."""
    return _aligned_span(rd) >= min_span and (max_softclip is None or _softclip_len(rd) <= max_softclip)


def _poa_consensus(seqs: list) -> str:
    """POA consensus (partial-order alignment, `pyabpoa`). Aligns the reads TO EACH OTHER (not
    against a fixed reference) -> better realignment of the tandem units at low depth than the
    pileup-to-reference. Lazy import (pyabpoa optional)."""
    try:
        import pyabpoa as pa
    except ImportError as e:
        raise RuntimeError("pyabpoa missing : `pip install --user pyabpoa`") from e
    seqs = [s for s in seqs if s]
    if not seqs:
        return ""
    res = pa.msa_aligner().msa(seqs, out_cons=True, out_msa=False)
    return res.cons_seq[0] if getattr(res, "cons_seq", None) else ""


def scaffold_hp(bam: str, vntr_ref: str, length: int, *, region: str = None, hp: str = None,
                genome_ref: str = None, workdir: str = None, preset: str = "map-ont",
                threads: int = 4, min_mq: int = 0, max_mismatch: int = 3, method: str = "pileup",
                keep: bool = True, min_span_frac: float = 0.0, max_softclip=None) -> dict:
    """DENSE consensus of an allele : reads (optionally HP) re-aligned on the SINGLE contig of length
    `length`, consensus over ALL these reads, decomposition + focal guard.

    `method` : "pileup" (`build_consensus_from_readset`, default) or "poa" (`pyabpoa`, aligns the reads
    to each other -> better realignment at low depth; the sequences are taken REFERENCE-ORIENTED
    from the aligned BAM, otherwise POA would mix the strands).

    `min_span_frac`/`max_softclip` : span/soft-clip filter **in re-aligned VNTR space** (Bam_cleaner on
    OUTPUT, pre-consensus) — keeps the reads covering >= `min_span_frac x contig_length` bp and
    soft-clipping <= `max_softclip`. Default off (0.0/None) -> no change. Purpose : remove the partial
    reads that shift the consensus (array-wide artifact). Warning: measure before making it default.

    Returns {length, contig, method, n_reads, n_reads_seen, consensus_len, has_mut, variants,
              mutation_variant, array_wide_artifact, n_frameshift_raw}.
    """
    import shutil
    for tool in ("samtools", "minimap2"):
        if shutil.which(tool) is None:
            raise RuntimeError(f"required tool missing from PATH : {tool}")
    M = _import_analyzer()
    loc = GRCh38["LOCUS"]
    region = region or f"{loc.chrom}:{loc.start}-{loc.end}"
    tmp = workdir or tempfile.mkdtemp(prefix="muc1_scaffold_")
    os.makedirs(tmp, exist_ok=True)
    tag = f"hp{hp}" if hp in ("1", "2") else "all"
    contig = f"MUC1_VNTR_{length}repeats"

    # 1. extract the reads (optionally per HP) -> fastq
    hp_bam = os.path.join(tmp, f"{tag}.bam")
    fq = os.path.join(tmp, f"{tag}.fq")
    view = ["samtools", "view", "-b"]
    if genome_ref:
        view += ["-T", genome_ref]
    if hp in ("1", "2"):
        view += ["-d", f"HP:{hp}"]
    view += ["-o", hp_bam, bam, region]
    _run(view)
    _run(["samtools", "fastq", "-0", fq, "-n", hp_bam])

    # 2. ref = ONE SINGLE contig (the right length) -> MAPQ becomes informative again
    #    samtools faidx writes the extracted contig to stdout (captured by _run, text=True).
    single = os.path.join(tmp, f"{contig}.fa")
    r = _run(["samtools", "faidx", vntr_ref, contig])
    with open(single, "w") as fh:
        fh.write(r.stdout)
    _run(["samtools", "faidx", single])

    # 3. align the allele's reads on this single contig -> sorted bam
    sam = os.path.join(tmp, f"{tag}.sam")
    bamr = os.path.join(tmp, f"{tag}.realigned.bam")
    with open(sam, "wb") as sh, open(os.path.join(tmp, f"{tag}.mm2.log"), "wb") as eh:
        import subprocess
        rc = subprocess.run(["minimap2", "-ax", preset, "-t", str(threads), single, fq],
                            stdout=sh, stderr=eh)
    if rc.returncode != 0:
        raise RuntimeError(f"minimap2 failed (scaffold {tag})")
    _run(["samtools", "sort", "-@", str(threads), "-o", bamr, sam])
    _run(["samtools", "index", bamr])

    # 4. DENSE consensus (all the allele's reads) -> decomposition -> focal guard
    if method == "poa":
        seqs = []                                   # REFERENCE-ORIENTED sequences (pysam)
        with pysam.AlignmentFile(bamr, "rb") as af:
            for rd in af.fetch(until_eof=True):
                if rd.is_unmapped or rd.is_secondary or rd.is_supplementary:
                    continue
                if rd.query_sequence:
                    seqs.append(rd.query_sequence)
        n_reads = n_seen = len(seqs)                 # POA : no span filter (branch removed anyway)
        consensus = _poa_consensus(seqs)
    else:
        read_names = set()
        n_seen = 0
        with pysam.AlignmentFile(bamr, "rb") as af:
            contig_len = af.get_reference_length(contig)
            min_span = int(min_span_frac * contig_len)      # span filter IN VNTR SPACE (pre-consensus)
            for rd in af.fetch(until_eof=True):
                if rd.is_unmapped:
                    continue
                n_seen += 1
                if (min_span or max_softclip is not None) and not _read_passes(rd, min_span, max_softclip):
                    continue
                read_names.add(rd.query_name)
        n_reads = len(read_names)
        ref_fa = pysam.FastaFile(single)
        consensus = M.build_consensus_from_readset(bamr, ref_fa, contig, read_names, min_mq=min_mq)
        ref_fa.close()
    matched, _ = M.match_motifs(consensus, max_mismatch=max_mismatch)
    call = _focal_call(frameshift_variants(matched), length)

    out = {"length": length, "contig": contig, "method": method, "n_reads": n_reads,
           "n_reads_seen": n_seen, "min_span_frac": min_span_frac,
           "consensus_len": len(consensus), **call,
           "workdir": tmp if (keep or workdir) else None}
    if not keep and not workdir:
        shutil.rmtree(tmp, ignore_errors=True)
    return out


# ── 3a+ : diplotype — assignment of reads by HP **and by LENGTH** ───────────────
# Motivation (one sample) : whatshap only tags a fraction of the reads; the mutant allele (hp2) dropped
# to 15 reads -> artifact. Here we ALSO recover the UN-tagged reads by routing them by length
# (competitive alignment on the 2 allele contigs) -> dense consensus of BOTH alleles.

def _region_fastq_and_hp(bam, region, genome_ref, tmp):
    """Extract ALL reads of the region -> fastq + map {query_name: HP '1'/'2'/None}."""
    region_bam = os.path.join(tmp, "region.bam")
    fq = os.path.join(tmp, "region.fq")
    view = ["samtools", "view", "-b"]
    if genome_ref:
        view += ["-T", genome_ref]
    view += ["-o", region_bam, bam, region]
    _run(view)
    hp_map = {}
    with pysam.AlignmentFile(region_bam, "rb") as af:
        for rd in af.fetch(until_eof=True):
            if rd.is_secondary or rd.is_supplementary:
                continue
            hp_map[rd.query_name] = str(rd.get_tag("HP")) if rd.has_tag("HP") else None
    _run(["samtools", "fastq", "-0", fq, "-n", region_bam])
    return fq, hp_map


def _align_single(fq, vntr_ref, length, tmp, threads, preset):
    """Align the fastq on THE single contig of length `length` -> (sorted+indexed bam, single.fa, contig)."""
    import subprocess
    contig = f"MUC1_VNTR_{length}repeats"
    single = os.path.join(tmp, f"{contig}.fa")
    r = _run(["samtools", "faidx", vntr_ref, contig])
    with open(single, "w") as fh:
        fh.write(r.stdout)
    _run(["samtools", "faidx", single])
    sam = os.path.join(tmp, f"{length}.sam")
    bamr = os.path.join(tmp, f"{length}.bam")
    with open(sam, "wb") as sh, open(os.path.join(tmp, f"{length}.mm2.log"), "wb") as eh:
        rc = subprocess.run(["minimap2", "-ax", preset, "-t", str(threads), single, fq],
                            stdout=sh, stderr=eh)
    if rc.returncode != 0:
        raise RuntimeError(f"minimap2 failed (L={length})")
    _run(["samtools", "sort", "-@", str(threads), "-o", bamr, sam])
    _run(["samtools", "index", bamr])
    return bamr, single, contig


def _scores(bam):
    """{query_name: AS alignment score} for the mapped primary reads."""
    out = {}
    with pysam.AlignmentFile(bam, "rb") as af:
        for rd in af.fetch(until_eof=True):
            if rd.is_unmapped or rd.is_secondary or rd.is_supplementary:
                continue
            out[rd.query_name] = rd.get_tag("AS") if rd.has_tag("AS") else 0
    return out


def scaffold_diplotype(bam, vntr_ref, len_a: int, len_b: int, *, region=None, genome_ref=None,
                       workdir=None, preset="map-ont", threads=4, min_mq=0, max_mismatch=3,
                       keep=True) -> dict:
    """Dense consensus of BOTH alleles, recovering the un-tagged reads by LENGTH.

    Warning: **FAILED APPROACH for long VNTRs** (tested on one sample 2026-07-04) : competitive alignment
    does NOT distinguish the length of a tandem (a repetitive read aligns just as well on any contig
    long enough -> ~equal AS score) -> all the reads stick on the `len_a` contig, the mutant allele is
    drowned. The true length signal = flank-to-flank span, rare for long VNTRs (soft-clip). For
    that sample, the un-tagged reads are unassignable -> prefer `scaffold_hp` (per HP) + 3b to densify
    the mutant allele. Kept for reference / short-VNTR cases where the reads span.

    1. extract all region reads (+ HP map); 2. align on the 2 contigs (len_a, len_b);
    3. route each read to an allele : by HP if tagged (majority length of that HP), else by
    best alignment score (competitive); 4. dense consensus per allele (all its reads).
    """
    import shutil
    for tool in ("samtools", "minimap2"):
        if shutil.which(tool) is None:
            raise RuntimeError(f"required tool missing from PATH : {tool}")
    M = _import_analyzer()
    loc = GRCh38["LOCUS"]
    region = region or f"{loc.chrom}:{loc.start}-{loc.end}"
    tmp = workdir or tempfile.mkdtemp(prefix="muc1_diplo_")
    os.makedirs(tmp, exist_ok=True)

    fq, hp_map = _region_fastq_and_hp(bam, region, genome_ref, tmp)
    bam_a, fa_a, contig_a = _align_single(fq, vntr_ref, len_a, tmp, threads, preset)
    bam_b, fa_b, contig_b = _align_single(fq, vntr_ref, len_b, tmp, threads, preset)
    sc_a, sc_b = _scores(bam_a), _scores(bam_b)

    # majority length per HP (among the tagged reads) -> establishes the HP<->length correspondence
    hp_len = {}
    for hp in ("1", "2"):
        tagged = [q for q, h in hp_map.items() if h == hp]
        va = sum(1 for q in tagged if sc_a.get(q, -1) >= sc_b.get(q, -1) and q in sc_a)
        vb = sum(1 for q in tagged if q in sc_b and sc_b.get(q, -1) > sc_a.get(q, -1))
        hp_len[hp] = len_a if va >= vb else len_b

    # route each read to an allele : HP if tagged, else best score
    reads = {len_a: set(), len_b: set()}
    for q in set(sc_a) | set(sc_b):
        h = hp_map.get(q)
        if h in ("1", "2"):
            reads[hp_len[h]].add(q)
        else:
            reads[len_a if sc_a.get(q, -1) >= sc_b.get(q, -1) else len_b].add(q)

    def _call(length, bamr, fa, contig, rset):
        ref_fa = pysam.FastaFile(fa)
        cons = M.build_consensus_from_readset(bamr, ref_fa, contig, rset, min_mq=min_mq)
        ref_fa.close()
        matched, _ = M.match_motifs(cons, max_mismatch=max_mismatch)
        call = _focal_call(frameshift_variants(matched), length)
        n_tag = sum(1 for q in rset if hp_map.get(q) in ("1", "2"))
        return {"length": length, "contig": contig, "n_reads": len(rset),
                "n_tagged": n_tag, "n_untagged": len(rset) - n_tag,
                "consensus_len": len(cons), **call}

    allele_a = _call(len_a, bam_a, fa_a, contig_a, reads[len_a])
    allele_b = _call(len_b, bam_b, fa_b, contig_b, reads[len_b])
    carriers = [a for a in (allele_a, allele_b) if a["has_mut"]]
    conflict = len(carriers) >= 2
    mut = None if conflict else (carriers[0] if carriers else None)

    out = {"mutation_present": bool(carriers), "mut_hp_conflict": conflict,
           "mutation_variant": mut["mutation_variant"] if mut else None,
           "mut_length": mut["length"] if mut else None,
           "hp_len": hp_len, "allele_a": {k: v for k, v in allele_a.items() if k != "variants"},
           "allele_b": {k: v for k, v in allele_b.items() if k != "variants"},
           "workdir": tmp if (keep or workdir) else None}
    if not keep and not workdir:
        shutil.rmtree(tmp, ignore_errors=True)
    return out


if __name__ == "__main__":   # test 1 patient (e.g. a sample) : scaffold_hp OR diplotype
    import argparse
    import json
    ap = argparse.ArgumentParser(prog="muc1_analyzer.detectors.vntr_scaffold",
                                 description="Dense allele consensus (allele rescue)")
    ap.add_argument("-b", "--bam", required=True)
    ap.add_argument("-r", "--vntr-ref", required=True)
    ap.add_argument("-L", "--length", type=int, required=True, help="allele length (n repeats)")
    ap.add_argument("--len-b", type=int, default=None,
                    help="2nd length -> DIPLOTYPE mode (recovers the un-tagged reads by length)")
    ap.add_argument("--hp", choices=["1", "2"], default=None)
    ap.add_argument("--method", choices=["pileup", "poa"], default="pileup",
                    help="poa = consensus via pyabpoa (better realignment at low depth)")
    ap.add_argument("--genome-ref", default=None)
    ap.add_argument("--min-mq", type=int, default=0)
    ap.add_argument("--min-span-frac", type=float, default=0.0,
                    help="span filter IN VNTR SPACE (pre-consensus) : keep the reads covering >= frac x contig length")
    ap.add_argument("--max-softclip", type=int, default=None,
                    help="max soft-clip (VNTR space); None = no bound")
    args = ap.parse_args()
    if args.len_b:
        res = scaffold_diplotype(args.bam, args.vntr_ref, args.length, args.len_b,
                                 genome_ref=args.genome_ref, min_mq=args.min_mq)
    else:
        res = scaffold_hp(args.bam, args.vntr_ref, args.length, hp=args.hp,
                          genome_ref=args.genome_ref, min_mq=args.min_mq, method=args.method,
                          min_span_frac=args.min_span_frac, max_softclip=args.max_softclip)
        res.pop("variants", None)
    print(json.dumps(res, indent=2, ensure_ascii=False))
