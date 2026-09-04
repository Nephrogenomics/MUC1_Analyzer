#!/usr/bin/env python3
"""vntr_scaffold (3a) — a DENSE per-allele consensus for the fragmented long/mutant allele.

Why
---
The caller ranks ~150 near-identical length-contigs and consensuses only the top-2; a real allele whose
spanning reads scatter across many contigs (typically the LONG allele) is never consensused densely, so
its frameshift is invisible to the plain call — only a supplied `--mut-variant` recovers it, and the
two-axis COVERAGE-GATE then (correctly) flags the plain call as unreliable. But the reads EXIST: the
alignment-free length arbiter bins them fine (on a real allele-specific cohort, a depleted ~90-copy
allele still carried ~18 spanning reads), they are just fragmented across length-contigs.

What 3a does
------------
1. BIN the region reads by their alignment-free flank-to-flank copy number (`vntr_raw_length.read_length`),
   keeping those within `tol` of the target allele.
2. RE-ALIGN that pooled read set to a SINGLE representative contig (the N-repeat idealised contig closest
   to the allele length) with minimap2 (`prepare.align`) → one dense pile-up instead of ~150 thin ones.
3. CONSENSUS + DECODE with the unchanged maison layer (`caller.build_consensus_for_contig` + `match_motifs`)
   → the frameshift on the depleted allele, now on a dense bi-strand consensus.

The frameshift decoder is unchanged; 3a only feeds it a denser consensus. Alignment is via **mappy** (the
minimap2 Python binding — no CLI, so no login-node watchdog and no external tool). `scaffold_both_alleles`
bins+decodes BOTH alleles in one pass (merging Ilias' span_bin_call design), and `caller.call` AUTO-TRIGGERS
it when it detects a multi-contig SMEAR (flat contig ranking) — so the smeared per-allele nomenclature +
frameshift are recovered end-to-end, no separate manual step. Still runnable standalone (`--both`).

STATUS (validated on a real allele-specific cohort): the densification MECHANICALLY works (a depleted
~90-copy allele → a dense ~93-unit bi-strand consensus; reproduces a known short-allele 59dupC at the
right repeat; a non-carrier → none), but it adds NO onset on that cohort. 3a can only help when the
mutant allele is fragmented OUT of the caller's top-2 contigs AND its frameshift is ONT-decodable
(≥~10 bp) — and those never co-occur there: the fragmented (depleted long) alleles all carry 1-bp indels
that are BASECALLING-limited even at dense coverage, while a decodable dupC sits on a SHORT allele that is
already a top-2 contig (the plain call already decodes it). Kept as a conditional densification tool for a
future sample with a large frameshift on a fragmented allele; NOT wired into the onset path. See
docs/VNTR_read_recovery_design.md.
"""
from __future__ import annotations

import argparse
import os
import sys

# GAP / read_length live in the top-level arbiter (no pysam import cost beyond what the caller already pays).
import vntr_raw_length as V


def copies_of(seq: str, maxmm: int, offset: int):
    """Alignment-free copy number of a read sequence (flank-to-flank), or None if it does not span."""
    g = V.read_length(seq, maxmm)
    if g is None:
        return None
    return round((g - V.GAP) / 60) + offset


def select_for_allele(reads, target, *, tol=4, maxmm=9, offset=0, copies_fn=copies_of):
    """Pick reads whose copy number is within `tol` of `target`.

    `reads` = iterable of (name, seq, qual). Returns [(name, seq, qual, copies)] — the pooled allele bin.
    `copies_fn` is injectable so the binning window can be unit-tested without real flank anchors."""
    out = []
    for name, seq, qual in reads:
        c = copies_fn(seq, maxmm, offset)
        if c is not None and abs(c - target) <= tol:
            out.append((name, seq, qual, c))
    return out


def _iter_primary(bam_path):
    """Yield (name, seq, qual_string) for each primary read of a BAM, in the read's NATIVE (basecalled)
    orientation. Using get_forward_sequence() (not query_sequence, which the original alignment already
    re-oriented to its contig's + strand) lets minimap RE-DERIVE the strand on the scaffold contig, so the
    pooled consensus recovers true bi-strand support instead of collapsing to one strand."""
    import pysam
    seen = set()
    with pysam.AlignmentFile(bam_path) as af:
        for r in af.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.query_name is None:
                continue
            seq = r.get_forward_sequence()
            if seq is None:
                continue
            if r.query_name in seen:
                continue
            seen.add(r.query_name)
            fq = r.get_forward_qualities()
            qual = pysam.qualities_to_qualitystring(fq) if fq is not None else "I" * len(seq)
            yield r.query_name, seq.upper(), qual


def _contig_for_copies(vntr_ref, target_copies):
    """(name, length) of the reference contig whose repeat count is closest to `target_copies`.

    The VNTR ref is a family `MUC1_VNTR_<N>repeats` (each = flank + N×60 bp). We pick by the arbiter's
    COPY NUMBER (parsed from the contig NAME), NOT by read length: the prep reads carry ~10 kb of extra
    genomic flank beyond the AL/AH anchors, so their full length massively over-states N — while N×60
    alone under-states it (ignores the ~7500 bp reference flank). The read's genomic overhang simply
    soft-clips onto the correctly-sized contig."""
    import re
    import pysam
    fa = pysam.FastaFile(vntr_ref)
    try:
        pairs = list(zip(fa.references, fa.lengths))
    finally:
        fa.close()
    parsed = [(name, ln, int(m.group(1))) for name, ln in pairs
              if (m := re.search(r"(\d+)\s*repeats", name))]
    if parsed:
        name, ln, _n = min(parsed, key=lambda p: abs(p[2] - target_copies))
        return name, ln
    # fallback: no parseable repeat count → nearest by (flank-agnostic) length proxy N×60
    return min(pairs, key=lambda p: abs(p[1] - target_copies * 60))


def _write_single_contig_ref(vntr_ref, contig, out_fa):
    """Extract one contig from `vntr_ref` into its own FASTA (+ .fai) so minimap forces all reads onto it."""
    import pysam
    fa = pysam.FastaFile(vntr_ref)
    try:
        seq = fa.fetch(contig)
    finally:
        fa.close()
    with open(out_fa, "w") as fh:
        fh.write(f">{contig}\n")
        for i in range(0, len(seq), 70):
            fh.write(seq[i:i + 70] + "\n")
    pysam.faidx(out_fa)
    return out_fa


def _align_bin_mappy(binned, mini_ref, contig, contig_len, out_bam, *, preset="map-ont", threads=4):
    """Align a length-bin to its SINGLE representative contig with **mappy** (the minimap2 Python binding) —
    no minimap2 CLI, so it never trips the login-node watchdog and needs no external tool. Writes a
    sorted+indexed BAM (one contig) that `build_consensus_for_contig` piles up directly. Each read keeps its
    NATIVE orientation input → mappy re-derives the strand on the scaffold → true bi-strand consensus."""
    import mappy
    import pysam
    aligner = mappy.Aligner(mini_ref, preset=preset, n_threads=threads)
    if not aligner:
        raise RuntimeError(f"mappy could not index {mini_ref}")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": contig, "LN": contig_len}]}
    tmp = out_bam + ".unsorted.bam"
    n = 0
    with pysam.AlignmentFile(tmp, "wb", header=header) as out:
        for name, seq, qual, _c in binned:
            hit = None
            for h in aligner.map(seq):           # mappy yields the primary first
                if h.is_primary:
                    hit = h
                    break
                hit = hit or h
            if hit is None:
                continue
            a = pysam.AlignedSegment(out.header)
            a.query_name = name
            if hit.strand == -1:
                s, q = mappy.revcomp(seq), qual[::-1]
                lclip, rclip = len(seq) - hit.q_en, hit.q_st           # clips swap in revcomp space
            else:
                s, q = seq, qual
                lclip, rclip = hit.q_st, len(seq) - hit.q_en
            a.query_sequence = s
            a.query_qualities = pysam.qualitystring_to_array(q)
            a.flag = 16 if hit.strand == -1 else 0
            a.reference_id = 0
            a.reference_start = hit.r_st
            a.mapping_quality = min(hit.mapq, 60)
            cig = ([(4, lclip)] if lclip > 0 else []) \
                + [(op, ln) for ln, op in hit.cigar] \
                + ([(4, rclip)] if rclip > 0 else [])
            a.cigartuples = cig
            out.write(a)
            n += 1
    pysam.sort("-o", out_bam, tmp)
    pysam.index(out_bam)
    os.remove(tmp)
    return n


def _strand_counts(bam_path, contig, min_mq):
    from .caller import _strand_counts_for_contig
    return _strand_counts_for_contig(bam_path, contig, min_mq)


def scaffold_allele(prep_bam, target, vntr_ref, *, workdir, tol=4, maxmm=9, offset=0,
                    min_depth=3, min_bq=20, min_mq=0, max_mismatch=3,
                    ins_threshold=0.5, long_ins_threshold=0.4, del_threshold=0.6,
                    threads=4, sample="scaffold", pacbio=False, binned=None, verbose=False):
    """Pool the reads of one allele → re-align to a representative contig (via mappy) → dense consensus →
    decode. Returns a JSON-friendly dict. `n_binned` = reads pooled by length; `n_aligned`/coverage = the
    dense pile-up the frameshift was decoded on; `indel` = the frameshift summary (or '' if none).
    `binned` may be pre-supplied (both-allele path bins once and reuses); else it is computed here."""
    import pysam
    from . import caller as C

    os.makedirs(workdir, exist_ok=True)
    if binned is None:
        binned = select_for_allele(_iter_primary(prep_bam), target, tol=tol, maxmm=maxmm, offset=offset)
    if not binned:
        return {"target": target, "n_binned": 0, "available": False,
                "note": "no spanning reads within tol of the target allele"}

    # representative contig: chosen by the arbiter COPY NUMBER (parsed from the contig name), not read
    # length — the reads over-hang the anchors by ~10 kb of genomic flank (that overhang soft-clips).
    import statistics
    median_len = int(statistics.median(len(seq) for _n, seq, _q, _c in binned))
    contig, contig_len = _contig_for_copies(vntr_ref, target)

    mini_ref = _write_single_contig_ref(vntr_ref, contig, os.path.join(workdir, f"{sample}_{contig}.fa"))
    out_bam = os.path.join(workdir, f"{sample}_allele{target}.bam")
    preset = "map-hifi" if pacbio else "map-ont"
    n_aligned = _align_bin_mappy(binned, mini_ref, contig, contig_len, out_bam, preset=preset, threads=threads)
    ref_fa = pysam.FastaFile(mini_ref)
    try:
        consensus = C.build_consensus_for_contig(
            bam_path=out_bam, ref_fa=ref_fa, contig=contig,
            min_depth=min_depth, min_bq=min_bq, min_mq=min_mq,
            ins_threshold=ins_threshold, long_ins_threshold=long_ins_threshold,
            del_threshold=del_threshold, verbose=verbose)
    finally:
        ref_fa.close()
    matched, unmatched = C.match_motifs(consensus, max_mismatch=max_mismatch, verbose=verbose)
    cov = _strand_counts(out_bam, contig, min_mq)

    return {
        "target": target, "contig": contig, "contig_len": contig_len, "median_read_len": median_len,
        "n_binned": len(binned), "n_aligned": n_aligned, "available": True, "coverage": cov,
        "nomenclature": C.haplotype_string(matched), "n_motifs": len(matched),
        "indel": C._indel_summary(matched), "consensus_len": len(consensus),
    }


def _phasing_concordance(bam, contig, snps, A, B, min_grp_support=5, maj_frac=0.8):
    """Count SNP positions that TRULY discriminate the two phased groups: A's majority base != B's, each
    backed by >= maj_frac of its group. Shared systematic errors look like SNPs but are the SAME in both
    groups -> they do not discriminate. So this count is ~0 for a genuine homozygote (a false split on
    noise) and equals the real differences for a true same/near-length heterozygote."""
    import collections
    import pysam
    A, B, snpset, discrim = set(A), set(B), set(snps), 0
    with pysam.AlignmentFile(bam) as af:
        for col in af.pileup(contig, truncate=True, min_base_quality=0):
            if col.reference_pos not in snpset:
                continue
            ca, cb = collections.Counter(), collections.Counter()
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip or pr.query_position is None:
                    continue
                nm = pr.alignment.query_name
                base = pr.alignment.query_sequence[pr.query_position]
                if nm in A:
                    ca[base] += 1
                elif nm in B:
                    cb[base] += 1
            na, nb = sum(ca.values()), sum(cb.values())
            if na < min_grp_support or nb < min_grp_support:
                continue
            (ba, va), (bb, vb) = ca.most_common(1)[0], cb.most_common(1)[0]
            if ba != bb and va >= maj_frac * na and vb >= maj_frac * nb:
                discrim += 1
    return discrim


def _snp_phase_alleles(all_binned, vntr_ref, dominant_len, workdir, *, pacbio, min_mq=0,
                       min_snps=3, min_reads=5, threads=4):
    """FOREIGN/PacBio only — try to split a LENGTH-COLLAPSED read cluster into two alleles by SEQUENCE.

    When `call_alleles` returns a single length (two near-length alleles merged by `sep`, e.g. 79|84 or
    46|47), the length cannot separate them and pooling both onto one consensus DILUTES a carrier's +1
    below the 50 % rule -> a false negative. Instead: realign every spanning read onto the ONE dominant
    ruler contig, find phasing SNPs, phase, and gate the split on `_phasing_concordance` (>= `min_snps`
    truly-discriminating positions, each group >= `min_reads`) so ONT/HiFi noise never fabricates a split.
    Returns [(lenA, groupA_binned), (lenB, groupB_binned)] with each len = the MODAL copy number of its
    group, or None (caller falls back to length binning). Fail-closed: any error -> None."""
    import os
    import collections
    import pysam
    from . import caller as C
    try:
        if len(all_binned) < 2 * min_reads:
            return None
        contig, contig_len = _contig_for_copies(vntr_ref, dominant_len)
        mini_ref = _write_single_contig_ref(vntr_ref, contig, os.path.join(workdir, f"phase_{contig}.fa"))
        out_bam = os.path.join(workdir, "phase_all.bam")
        preset = "map-hifi" if pacbio else "map-ont"
        n_aln = _align_bin_mappy(all_binned, mini_ref, contig, contig_len, out_bam,
                                 preset=preset, threads=threads)
        if n_aln < 2 * min_reads:
            return None
        ref_len = pysam.FastaFile(mini_ref).get_reference_length(contig)
        # het positions: minor allele 30-70 % (two equimolar alleles), min_mq inherited from the caller
        snps = C.find_phasing_snps(out_bam, contig, ref_len, min_mq=min_mq, min_depth=max(5, 2 * min_reads),
                                   min_af=0.30, max_af=0.70)
        if not snps or len(snps) < min_snps:
            return None
        rA, rB = C.phase_reads(out_bam, contig, ref_len, snps, min_mq=min_mq)
        names = {r[0] for r in all_binned}
        A = [n for n in rA if n in names]
        B = [n for n in rB if n in names]
        if len(A) < min_reads or len(B) < min_reads:
            return None
        if _phasing_concordance(out_bam, contig, snps, A, B) < min_snps:
            return None                                       # split not backed by real differences -> homozygote
        Aset, Bset = set(A), set(B)
        grpA = [r for r in all_binned if r[0] in Aset]
        grpB = [r for r in all_binned if r[0] in Bset]
        _mode = lambda g: collections.Counter(r[3] for r in g).most_common(1)[0][0]
        pair = sorted([(_mode(grpA), grpA), (_mode(grpB), grpB)], key=lambda x: x[0])
        return pair
    except Exception:
        return None


def scaffold_both_alleles(prep_bam, vntr_ref, *, workdir, tol=4, maxmm=9, offset=0, min_mq=0,
                          threads=4, sample="scaffold", pacbio=False, alleles=None, verbose=False,
                          copies_fn=copies_of, snp_phase=False, min_snps=3):
    """Length-drive the WHOLE call on a SMEARED sample (Ilias' span_bin idea, merged): measure each read's
    copy number ONCE (alignment-free), call the two alleles, bin each, and scaffold+decode BOTH → a dense
    per-allele consensus + frameshift where the flat multi-contig ranking gives only generic 'X'. `alleles`
    may be supplied (the caller already has the arbiter result) to skip re-measuring. Returns
    {available, alleles, per_allele:[scaffold_allele dict, …]}.

    `copies_fn(seq, maxmm, offset)` is injectable: it defaults to the FLANK-to-flank measure (`copies_of`,
    valid on our native LR-PCR), but a FOREIGN amplicon lacking our AL/AH anchors (e.g. VNTRtools on PacBio
    HiFi) can pass the amplicon-agnostic CASSETTE measure so the binning still works. Only the initial
    binning uses it — per-allele reads are then pre-supplied to `scaffold_allele`.

    `snp_phase` (FOREIGN/PacBio path only): when the length caller collapses to a SINGLE length, try to
    split two near-length alleles by SNP first (``_snp_phase_alleles``, gated on `min_snps` discriminating
    positions). Fail-closed → any failure keeps the current length-binning. Never enabled for native ONT."""
    reads = list(_iter_primary(prep_bam))
    data = [(copies_fn(seq, maxmm, offset), name, seq, qual) for name, seq, qual in reads]
    data = [d for d in data if d[0] is not None]
    if not data:
        return {"available": False, "note": "no reads span the VNTR anchors (flank or cassette)", "per_allele": []}
    if alleles is None:
        alleles = V.call_alleles(sorted(c for c, *_ in data)).get("alleles", [])
    # SNP-phase a length-collapsed cluster (near-length het that `sep` merged) — foreign/PacBio only.
    if snp_phase and len(alleles) == 1:
        all_binned = [(n, s, q, c) for c, n, s, q in data]
        phased = _snp_phase_alleles(all_binned, vntr_ref, alleles[0], workdir, pacbio=pacbio,
                                    min_mq=min_mq, min_snps=min_snps, threads=threads)
        if phased:
            out = []
            for tgt, grp in phased:
                out.append(scaffold_allele(prep_bam, tgt, vntr_ref, workdir=workdir, tol=tol, maxmm=maxmm,
                                           offset=offset, min_mq=min_mq, threads=threads, sample=sample,
                                           pacbio=pacbio, binned=grp, verbose=verbose))
            return {"available": bool(out), "alleles": [t for t, _ in phased],
                    "n_spanning": len(data), "per_allele": out, "phased_by_snp": True}
    out = []
    for target in alleles:
        bin_reads = [(n, s, q, c) for c, n, s, q in data if abs(c - target) <= tol]
        res = scaffold_allele(prep_bam, target, vntr_ref, workdir=workdir, tol=tol, maxmm=maxmm,
                              offset=offset, min_mq=min_mq, threads=threads, sample=sample,
                              pacbio=pacbio, binned=bin_reads, verbose=verbose)
        out.append(res)
    return {"available": bool(out), "alleles": alleles, "n_spanning": len(data), "per_allele": out}


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="muc1_analyzer.vntr_scaffold",
        description="3a rescue: dense per-allele consensus of a fragmented/depleted VNTR allele (mappy).")
    ap.add_argument("-b", "--prep-bam", required=True, help="region reads (from `prepare`)")
    ap.add_argument("-r", "--vntr-ref", required=True, help="the 150-contig VNTR reference FASTA")
    ap.add_argument("--copies", type=int, default=None,
                    help="target allele copies to rescue (default: the LONG allele from the arbiter)")
    ap.add_argument("--tol", type=int, default=4, help="copy-number window around the target (default 4)")
    ap.add_argument("--offset", type=int, default=0, help="copies added per read (0 AS/WGS/T2T, 4 PCR)")
    ap.add_argument("--maxmm", type=int, default=9, help="max mismatches for the flank anchors")
    ap.add_argument("--min-mq", type=int, default=0)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--both", action="store_true",
                    help="bin+scaffold BOTH alleles (length-driven per-allele call — the smear fix), not just the long")
    ap.add_argument("--pacbio", action="store_true", help="PacBio HiFi reads (mappy map-hifi; default map-ont)")
    ap.add_argument("-s", "--sample", default="scaffold")
    ap.add_argument("--workdir", default="scaffold_out")
    ap.add_argument("--json", default=None, help="write the result dict as JSON")
    a = ap.parse_args(argv)

    if a.both:
        res = scaffold_both_alleles(a.prep_bam, a.vntr_ref, workdir=a.workdir, tol=a.tol, maxmm=a.maxmm,
                                    offset=a.offset, min_mq=a.min_mq, threads=a.threads, sample=a.sample,
                                    pacbio=a.pacbio)
        print(f"[scaffold] alleles={res.get('alleles')} (n_spanning={res.get('n_spanning')})", file=sys.stderr)
        for r in res.get("per_allele", []):
            print(f"  allele {r.get('target')}: binned={r.get('n_binned')} aligned={r.get('n_aligned')} "
                  f"| {r.get('nomenclature', '')}  →  FRAMESHIFT: {r.get('indel') or '(none)'}")
        if a.json:
            import json
            with open(a.json, "w") as jf:
                json.dump(res, jf, indent=2)
            print(f"[scaffold] JSON → {a.json}", file=sys.stderr)
        return 0

    target = a.copies
    if target is None:
        cops = [c for c, _ in V.copies_from_bam(a.prep_bam, chrom=None, offset=a.offset, maxmm=a.maxmm)]
        call = V.call_alleles(cops)
        if not call["alleles"]:
            print("[scaffold] the arbiter found no allele peak — pass --copies explicitly.", file=sys.stderr)
            return 1
        target = max(call["alleles"])          # the LONG allele is the one the plain call fragments/misses
        print(f"[scaffold] arbiter alleles={call['alleles']} counts={call['counts']} → rescuing {target}",
              file=sys.stderr)

    res = scaffold_allele(a.prep_bam, target, a.vntr_ref, workdir=a.workdir, tol=a.tol, maxmm=a.maxmm,
                          offset=a.offset, min_mq=a.min_mq, threads=a.threads, sample=a.sample,
                          pacbio=a.pacbio)
    print(f"[scaffold] allele {res.get('target')} : binned={res.get('n_binned')} "
          f"aligned={res.get('n_aligned')} coverage={res.get('coverage')}", file=sys.stderr)
    print(f"[scaffold] nomenclature : {res.get('nomenclature', '')}")
    print(f"[scaffold] FRAMESHIFT   : {res.get('indel') or '(none decoded)'}")
    if a.json:
        import json
        with open(a.json, "w") as jf:
            json.dump(res, jf, indent=2)
        print(f"[scaffold] JSON → {a.json}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
