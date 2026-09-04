#!/usr/bin/env python3
"""consensus_arbiter.py — decide one-vs-two alleles for CLOSE lengths (Δ≤1) by COMPARING CONSENSUS
sequences at the BASE level.

WHY (the length arbiter alone is not enough). Two allele lengths that differ by ≥2 copies are two alleles
outright — a single allele's flank-length measurement does not scatter by 2. The ambiguous case is Δ≤1:
either two genuine alleles that happen to sit 0–1 copy apart (a 43|44 pair), or ONE allele whose measured
length scatters over adjacent bins (a length-homozygote, 44|44). Length cannot separate these; the
SEQUENCE can. Two distinct MUC1 alleles carry several base differences along the array, whereas one allele
split by measurement noise gives ~identical consensuses.

Method (the user's, at the BASE level rather than the motif level — a motif call can flip on a single
artefactual SNP, whereas a ≥3-base-difference threshold is robust to that):
  1. group the reads by measured length;
  2. align each group to ITS OWN length-contig (xavier) and build a majority consensus (the caller's rule);
  3. align the two consensuses against EACH OTHER (mappy) — the (L2−L1) extra repeats show up as a clean
     insertion GAP, and true allele differences show up as MISMATCHES outside the gap;
  4. count mismatches outside indels; ≥ min_diff ⇒ two distinct alleles, else one (homozygous / noise split).

The reference VNTR contigs are an ARTIFICIAL 1→150-times tandem of a single reference repeat — not the
patient's true sequence — so we never compare a group to the reference; we compare the two GROUP consensuses
to each other. The contig is only a scaffold to lay reads down and read a per-position majority.
"""
import collections
import os

import pysam
import mappy

import vntr_raw_length as V

XAVIER_KW = dict(preset="lr:hq", scoring=[2, 1, 10, 50, 10, 50],
                 bw=200, bw_long=7000, best_n=1, min_chain_score=100)


def _extract_contig(ref, name, out_fa):
    with pysam.FastaFile(ref) as fa:
        seq = fa.fetch(name)
    with open(out_fa, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + "\n")
    pysam.faidx(out_fa)
    return seq.upper()


def _align_group_to_contig(reads, names, contig_fa, out_bam):
    """Align each named read (reads[name]=seq or (seq,qual)) to the single contig in contig_fa (xavier),
    write a sorted+indexed BAM. Returns the contig name."""
    al = mappy.Aligner(contig_fa, **XAVIER_KW)
    if not al:
        raise RuntimeError(f"mappy failed to index {contig_fa}")
    ctg = al.seq_names[0]
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": ctg, "LN": len(al.seq(ctg))}]}
    H = pysam.AlignmentHeader.from_dict(header)
    tmp = out_bam + ".uns.bam"
    with pysam.AlignmentFile(tmp, "wb", header=header) as out:
        for name in names:
            r = reads[name]
            seq = r[0] if isinstance(r, tuple) else r
            best = None
            for h in al.map(seq):
                if h.is_primary:
                    best = h
                    break
                if best is None or h.mlen > best.mlen:
                    best = h
            if best is None:
                continue
            seg = pysam.AlignedSegment(H)
            seg.query_name = name
            if best.strand == 1:
                qs, lclip, rclip, flag = seq, best.q_st, len(seq) - best.q_en, 0
            else:
                qs, lclip, rclip, flag = V.rc(seq), len(seq) - best.q_en, best.q_st, 16
            cig = ([(4, lclip)] if lclip else []) \
                + [(op, l) for (l, op) in best.cigar] \
                + ([(4, rclip)] if rclip else [])
            seg.flag, seg.reference_id, seg.reference_start = flag, 0, best.r_st
            seg.mapping_quality, seg.query_sequence, seg.cigartuples = 60, qs, cig
            out.write(seg)
    pysam.sort("-o", out_bam, tmp)
    pysam.index(out_bam)
    os.remove(tmp)
    return ctg


def _consensus(bam, ctg, ref_seq, *, min_depth=3, min_bq=20, min_mq=20,
               ins_threshold=0.5, long_ins_threshold=0.25, del_threshold=0.6):
    """Majority consensus over the reads in `bam` aligned to `ctg` (the caller's build_consensus rule,
    verbatim: per-position majority base; consensual deletion drops the position; short/long insertions
    added above their fraction thresholds). Low-coverage positions fall back to the reference base."""
    INS_SHORT_MAX = 3
    ref_len = len(ref_seq)
    bases_at = collections.defaultdict(list)
    ins_after = collections.defaultdict(list)
    with pysam.AlignmentFile(bam, "rb") as af:
        for pcol in af.pileup(ctg, 0, ref_len, min_base_quality=min_bq,
                              min_mapping_quality=min_mq, stepper="all", truncate=True):
            pos = pcol.reference_pos
            for pr in pcol.pileups:
                if pr.is_refskip:
                    continue
                if pr.is_del:
                    bases_at[pos].append("-")
                else:
                    b = pr.alignment.query_sequence[pr.query_position].upper()
                    bases_at[pos].append(b)
                    if pr.indel > 0:
                        q0 = pr.query_position + 1
                        ins_after[pos].append(pr.alignment.query_sequence[q0:q0 + pr.indel].upper())
    parts = []
    for pos in range(ref_len):
        bases = bases_at[pos]
        depth = len(bases)
        if depth < min_depth:
            parts.append(ref_seq[pos])
        else:
            ctr = collections.Counter(bases)
            top_base, top_count = ctr.most_common(1)[0]
            if top_base == "-" and top_count / depth >= del_threshold:
                pass
            else:
                non_del = [b for b in bases if b != "-"]
                if non_del:
                    parts.append(collections.Counter(non_del).most_common(1)[0][0])
                elif top_base != "-":
                    parts.append(top_base)
        ins_list = ins_after.get(pos, [])
        if ins_list:
            best_ins, best_count = collections.Counter(ins_list).most_common(1)[0]
            ins_len = len(best_ins)
            if ins_len <= INS_SHORT_MAX:
                if ins_len == 1:
                    for delta in (-1, 1):
                        best_count += sum(1 for s in ins_after.get(pos + delta, []) if s == best_ins)
                    ref_depth = max(len(bases_at.get(max(0, pos - 1), [])), depth,
                                    len(bases_at.get(min(ref_len - 1, pos + 1), [])))
                    frac = best_count / max(ref_depth, 1)
                else:
                    frac = best_count / max(depth, 1)
                threshold = ins_threshold
            else:
                frac = best_count / max(depth, 1)
                threshold = long_ins_threshold
            if frac >= threshold:
                parts.append(best_ins)
    return "".join(parts)


def _count_mismatches(cons_a, cons_b):
    """Align cons_b to cons_a (mappy), count SUBSTITUTIONS outside indels. The (len difference) extra
    repeats appear as a clean insertion/deletion in the CIGAR and are NOT counted; only base-for-base
    substitutions in aligned blocks are. Returns (n_mismatch, n_aligned, detail)."""
    tmpfa = "_consA.fa"
    with open(tmpfa, "w") as f:
        f.write(">A\n")
        for i in range(0, len(cons_a), 60):
            f.write(cons_a[i:i + 60] + "\n")
    try:
        al = mappy.Aligner(tmpfa, preset="asm5")   # asm5: close sequences, honest mismatch/indel split
        if not al:
            al = mappy.Aligner(seq=cons_a, preset="asm10")
        hit = None
        for h in al.map(cons_b):
            if h.is_primary:
                hit = h
                break
            hit = hit or h
        if hit is None:
            return None, 0, {"note": "no alignment between consensuses"}
        # NM = edit distance; subtract indel bases to get substitutions
        indel_bases = sum(l for (l, op) in hit.cigar if op in (1, 2))
        nm = hit.NM if hasattr(hit, "NM") else None
        if nm is None:
            # fall back: count from cs/cigar mismatches
            mism = getattr(hit, "mlen", 0)
            return None, hit.blen, {"note": "no NM tag"}
        mism = nm - indel_bases
        return max(0, mism), hit.blen, {"NM": nm, "indel_bases": indel_bases,
                                        "aligned_len": hit.blen, "mlen": hit.mlen}
    finally:
        for e in ("", ".fai"):
            try:
                os.remove(tmpfa + e)
            except OSError:
                pass


def same_or_two_alleles(reads, names_a, names_b, contig_a, contig_b, ref, *,
                        min_diff=3, min_depth=3, outdir=None, return_detail=False):
    """Compare the two GROUPS' consensuses at the base level. Returns True if they are TWO distinct
    alleles (≥ min_diff substitutions between the consensuses), False if ONE allele (homozygous / noise
    split). `contig_a`/`contig_b` are each group's own length-contig scaffold.

      reads   : {name: seq} or {name: (seq, qual)}
      names_a : reads of allele-length A;  names_b : reads of allele-length B
      contig_a/contig_b : e.g. 'MUC1_VNTR_43repeats' / 'MUC1_VNTR_44repeats'
    """
    work = outdir or "."
    os.makedirs(work, exist_ok=True)
    fa_a = os.path.join(work, "_consarb_A.fa")
    fa_b = os.path.join(work, "_consarb_B.fa")
    bam_a = os.path.join(work, "_consarb_A.bam")
    bam_b = os.path.join(work, "_consarb_B.bam")
    detail = {"n_a": len(names_a), "n_b": len(names_b), "n_mismatch": None, "min_diff": min_diff}
    try:
        ref_a = _extract_contig(ref, contig_a, fa_a)
        ref_b = _extract_contig(ref, contig_b, fa_b)
        ctg_a = _align_group_to_contig(reads, list(names_a), fa_a, bam_a)
        ctg_b = _align_group_to_contig(reads, list(names_b), fa_b, bam_b)
        cons_a = _consensus(bam_a, ctg_a, ref_a, min_depth=min_depth)
        cons_b = _consensus(bam_b, ctg_b, ref_b, min_depth=min_depth)
        detail["len_cons_a"], detail["len_cons_b"] = len(cons_a), len(cons_b)
        mism, aln, mdet = _count_mismatches(cons_a, cons_b)
        detail.update(mdet)
        detail["n_mismatch"] = mism
        verdict = (mism is not None) and (mism >= min_diff)
    finally:
        for base in (fa_a, fa_b):
            for e in ("", ".fai"):
                try:
                    os.remove(base + e)
                except OSError:
                    pass
        for base in (bam_a, bam_b):
            for e in ("", ".bai"):
                try:
                    os.remove(base + e)
                except OSError:
                    pass
    return (verdict, detail) if return_detail else verdict
