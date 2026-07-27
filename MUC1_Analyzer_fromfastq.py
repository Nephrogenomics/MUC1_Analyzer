#!/usr/bin/env python3
"""MUC1_Analyzer_fromfastq.py — alignment-first VNTR haplotype caller for MUC1 (length → per-allele call).

STAGE 1 (length): align every read to ONE long "ruler" contig (the faked ref's MUC1_VNTR_150repeats,
which already carries the unique flanks AL/AH AND a 150-unit tandem ≥ any allele → a shorter allele is
DELETION-encoded, so the CIGAR count is exact and never capped by soft-clipped excess copies). A single
contig = no 150-way smear. Per read, transfer the two ref tandem edges (AL end, AH start) to read coords
THROUGH THE CIGAR and count read bases between them → exact copies for a read carrying BOTH its own flanks
(SPAN), a LOWER BOUND for one flank (LEFT/RIGHT). Call the 2 alleles from the FEW-but-ACCURATE SPAN reads
(low-count clustering; no min_reads=6 floor).

STAGE 2 (per-allele biology): bin reads to their allele (SPAN set the length; LEFT/RIGHT densify), re-align
each bin to its single MUC1_VNTR_<N>repeats contig, and run `muc1_analyzer call` per allele → per-allele
consensus + motifs + frameshift + PDF. The two per-allele calls are then MERGED and scored with the
project's own two-axis model (`detectors.vntr.from_analyzer_json` → `clinical_call.score_from_fields`),
giving the combined MUC1_Score. `--stage1-only` stops after the length call.

Run from the repo root (so `import vntr_raw_length` and `python -m muc1_analyzer` resolve). Deps: pysam, mappy.
"""
import argparse
import json
import os
import re
import subprocess
import sys
import statistics
import threading
from array import array
from collections import Counter, defaultdict

import pysam
import mappy

try:
    import vntr_raw_length as V   # AL, AH, GAP, fuzzy_find, rc — a repo-root module (as used by span_bin_call.py)
except ModuleNotFoundError:
    try:
        from muc1_analyzer import vntr_raw_length as V   # fallback: if it lives inside the package
    except ModuleNotFoundError:
        sys.exit("[MUC1_Analyzer_fromfastq] cannot import 'vntr_raw_length'. It must be importable — a "
                 "`vntr_raw_length.py` at the repository root (next to this script, as span_bin_call.py "
                 "also uses it) or inside the muc1_analyzer package. Run from the repo root, or restore it: "
                 "`git checkout <branch-that-has-it> -- vntr_raw_length.py`.")

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
_UNIT = 60
_RULER_CONTIG = "MUC1_VNTR_150repeats"
_CONTIG_TEMPLATE = "MUC1_VNTR_{n}repeats"


# ── ruler reference ─────────────────────────────────────────────────────────────
def build_ruler(faked_ref, contig, outdir, maxmm=9):
    """Extract the ruler contig, ORIENTED so AL precedes AH (the faked ref stores it revcomp'd), indexed.
    Orientation matters: the alignment coordinates and the AL/AH edge coordinates must share one frame."""
    if not os.path.exists(faked_ref + ".fai"):
        pysam.faidx(faked_ref)
    with pysam.FastaFile(faked_ref) as fa:
        if contig not in fa.references:
            raise SystemExit(f"[cigar-len] ruler contig {contig} absent from {faked_ref}")
        seq = fa.fetch(contig).upper()
    if not (V.fuzzy_find(V.AL, seq, maxmm) and V.fuzzy_find(V.AH, seq, maxmm)):
        rseq = V.rc(seq)
        if V.fuzzy_find(V.AL, rseq, maxmm) and V.fuzzy_find(V.AH, rseq, maxmm):
            seq = rseq
    ruler_fa = os.path.join(outdir, "ruler." + contig + ".fa")
    with open(ruler_fa, "w") as f:
        f.write(f">{contig}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + "\n")
    pysam.faidx(ruler_fa)
    return ruler_fa


def locate_tandem_edges(ref_fa, contig, maxmm):
    """AL/AH in the ruler → (gAL_end, gAH_start) [half-open tandem span in ref coords]."""
    with pysam.FastaFile(ref_fa) as fa:
        seq = fa.fetch(contig).upper()
    al, ah = V.fuzzy_find(V.AL, seq, maxmm), V.fuzzy_find(V.AH, seq, maxmm)
    if not (al and ah):
        rseq = V.rc(seq)
        al2, ah2 = V.fuzzy_find(V.AL, rseq, maxmm), V.fuzzy_find(V.AH, rseq, maxmm)
        if al2 and ah2:
            al, ah = al2, ah2
    if not (al and ah):
        raise SystemExit("[cigar-len] AL/AH anchors not found in the ruler contig.")
    gAL_end, gAH_start = al[0] + len(V.AL), ah[0]
    if gAH_start <= gAL_end:
        raise SystemExit("[cigar-len] AH is not above AL in the ruler — unexpected orientation.")
    ruler_copies = round((gAH_start - gAL_end - V.GAP) / _UNIT)
    return gAL_end, gAH_start, {"gAL_end": gAL_end, "gAH_start": gAH_start,
                                "anchor_span_bp": gAH_start - gAL_end, "ruler_copies": ruler_copies}


# ── read loading (qualities kept COMPACTLY as array('B'), 1 byte/base, not list(q)) ──
def iter_reads(path, ref_cram=None):
    """Yield (name, seq, quals) ONE read at a time — STREAMING, never holding the whole file in memory
    (an AS FASTQ can be tens of millions of reads / tens of GB). fastq(.gz)/fasta OR uBAM/BAM/CRAM
    (primary only). quals is an array('B') of Phred values or None. De-duplication of kept reads happens
    downstream by name, so no big `seen` set is needed here."""
    if str(path).endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz", ".fasta", ".fa", ".fa.gz")):
        with pysam.FastxFile(path) as fh:
            for e in fh:
                q = array("B", (ord(c) - 33 for c in e.quality)) if e.quality else None
                yield e.name, e.sequence.upper(), q
    else:
        mode = "rc" if str(path).endswith(".cram") else "rb"
        kw = {"reference_filename": ref_cram} if (mode == "rc" and ref_cram) else {}
        with pysam.AlignmentFile(path, mode, check_sq=False, **kw) as af:
            for r in af.fetch(until_eof=True):
                if r.is_secondary or r.is_supplementary or r.query_sequence is None:
                    continue
                q = r.query_qualities
                yield r.query_name, r.query_sequence.upper(), array("B", q) if q is not None else None


# ── per-read length from CIGAR + spanning gate from the read's own flanks ────────
def _q_in_window(r_st, cigar_lenop, ws, we):
    ref, q_in = r_st, 0
    for ln, op in cigar_lenop:
        if op in (0, 7, 8):
            a, b = ref, ref + ln
            q_in += max(0, min(b, we) - max(a, ws))
            ref = b
        elif op in (2, 3):
            ref += ln
        elif op == 1:
            if ws <= ref < we:
                q_in += ln
    return q_in, ref


def _has_anchor(seq, probe, maxmm):
    return bool(V.fuzzy_find(probe, seq, maxmm) or V.fuzzy_find(probe, V.rc(seq), maxmm))


def _seg(hdr, name, seq, quals, hit):
    s = pysam.AlignedSegment(hdr)
    s.query_name = name
    if hit.strand == 1:
        qseq, qq, lclip, rclip, flag = seq, quals, hit.q_st, len(seq) - hit.q_en, 0
    else:
        qseq = V.rc(seq)
        qq = quals[::-1] if quals is not None else None
        lclip, rclip, flag = len(seq) - hit.q_en, hit.q_st, 16
    cig = ([(4, lclip)] if lclip else []) + [(op, ln) for (ln, op) in hit.cigar] \
        + ([(4, rclip)] if rclip else [])
    s.flag, s.reference_id, s.reference_start = flag, 0, hit.r_st
    s.mapping_quality = min(hit.mapq, 60)
    s.query_sequence = qseq
    if qq is not None:
        s.query_qualities = array("B", qq)
    s.cigartuples = cig
    return s


def measure_reads(read_iter, ruler_fa, preset, ws, we, offset, maxmm, out_bam=None, threads=1, keep_thr=1):
    """Stream reads, align each to the ruler, and classify. To survive a tens-of-millions-of-reads AS
    FASTQ, only the reads that actually COVER the VNTR (copies >= keep_thr) are kept in memory (their
    sequence, record, and BAM segment); everything else (the 0-copy / off-target flood) is counted for
    the histogram and dropped. Returns (records, kept_reads, hist, n_total, n_aligned)."""
    aln = mappy.Aligner(ruler_fa, preset=preset, n_threads=max(1, threads))
    if not aln:
        raise SystemExit(f"[cigar-len] mappy failed to index {ruler_fa}")
    ctg = aln.seq_names[0]
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": ctg, "LN": len(aln.seq(ctg))}]}
    hdr = pysam.AlignmentHeader.from_dict(header)
    writer = tmp = None
    if out_bam:
        tmp = out_bam + ".uns.bam"
        writer = pysam.AlignmentFile(tmp, "wb", header=header)
    records, kept_reads = {}, {}
    hist = defaultdict(Counter)
    n_total = n_aligned = 0

    _tl = threading.local()

    def _work(item):
        """Align + classify ONE read (worker thread; mappy releases the GIL, own ThreadBuffer per thread).
        Returns a light tuple; the read's sequence rides along ONLY if it covers the VNTR (kept)."""
        name, seq, quals = item
        buf = getattr(_tl, "buf", None)
        if buf is None:
            buf = _tl.buf = mappy.ThreadBuffer()
        best = None
        for h in aln.map(seq, buf=buf):
            if h.is_primary:
                best = h
                break
            if best is None or h.mlen > best.mlen:
                best = h
        if best is None:
            return None
        q_in, _ = _q_in_window(best.r_st, best.cigar, ws, we)
        has_al, has_ah = _has_anchor(seq, V.AL, maxmm), _has_anchor(seq, V.AH, maxmm)
        if has_al and has_ah:
            cls, exact, copies = "SPAN", True, round((q_in - V.GAP) / _UNIT) + offset
        elif has_al:
            cls, exact, copies = "LEFT", False, max(0, round((q_in - V.GAP) / _UNIT))
        elif has_ah:
            cls, exact, copies = "RIGHT", False, max(0, round((q_in - V.GAP) / _UNIT))
        else:
            cls, exact, copies = "INT", False, max(0, round(q_in / _UNIT))
        rec = {"copies": copies, "class": cls, "exact": exact, "q_in": q_in, "mapq": min(best.mapq, 60)}
        if copies >= keep_thr:                                # covers the VNTR → keep seq/segment
            seg = _seg(hdr, name, seq, quals, best) if writer is not None else None
            return name, rec, seg, seq, quals
        return name, rec, None, None, None                    # off-target → record for histogram only

    def _consume(result):
        nonlocal n_aligned
        if result is None:
            return
        name, rec, seg, seq, quals = result
        n_aligned += 1
        hist[rec["class"]][rec["copies"]] += 1
        if seq is not None:
            records[name] = rec
            kept_reads[name] = (seq, quals)
            if writer is not None and seg is not None:
                writer.write(seg)

    if threads and threads > 1:
        from concurrent.futures import ThreadPoolExecutor
        from itertools import islice
        with ThreadPoolExecutor(max_workers=threads) as ex:
            while True:                                       # batch to bound memory (one batch at a time)
                chunk = list(islice(read_iter, 50_000))
                if not chunk:
                    break
                n_total += len(chunk)
                for result in ex.map(_work, chunk, chunksize=max(1, len(chunk) // (threads * 4))):
                    _consume(result)
    else:
        for item in read_iter:
            n_total += 1
            _consume(_work(item))

    if writer is not None:
        writer.close()
        pysam.sort("-o", out_bam, tmp)
        pysam.index(out_bam)
        os.remove(tmp)
    return records, kept_reads, hist, n_total, n_aligned


def call_alleles_lowcount(values, sep, merge_tol, min_reads_per_allele):
    if not values:
        return [], []
    vals = sorted(values)
    clusters = [[vals[0]]]
    for v in vals[1:]:
        (clusters[-1].append(v) if v - clusters[-1][-1] <= merge_tol else clusters.append([v]))
    cl = [{"center": int(round(statistics.median(c))), "support": len(c)} for c in clusters]
    cl = [c for c in cl if c["support"] >= min_reads_per_allele]
    if not cl:
        return [], []
    cl.sort(key=lambda c: (-c["support"], c["center"]))
    chosen = [cl[0]]
    for c in cl[1:]:
        if all(abs(c["center"] - k["center"]) >= sep for k in chosen):
            chosen.append(c)
        if len(chosen) == 2:
            break
    chosen.sort(key=lambda c: c["center"])
    return [c["center"] for c in chosen], [c["support"] for c in chosen]


def collapse_scatter_allele(alleles, support, sep, ratio_near, ratio_far):
    """Fold a parasitic second length-allele into the dominant one when it is likely scatter, not a true
    allele. Reads mis-measure onto NEIGHBOURING lengths (stutter around the true peak), so a second peak
    close to the first is far more suspect than a distant one → a distance-dependent relative threshold:
    a neighbouring second peak (Δ < sep) must reach `ratio_near` (default 0.5) of the dominant's support,
    a distant one only `ratio_far` (default 0.25). Below threshold → treated as scatter → length-homozygous
    (the folded reads stay in the pool for phasing). The absolute floor against tiny distant noise is
    `--min-reads-per-allele`, applied upstream in call_alleles_lowcount (never the ratio), so a true but
    shallow distant allele — e.g. 5 reads — is preserved. Returns (alleles, support, fold_info|None)."""
    if len(alleles) != 2:
        return alleles, support, None
    (L1, s1), (L2, s2) = sorted(zip(alleles, support), key=lambda x: -x[1])   # dominant first
    delta = abs(L1 - L2)
    thr = ratio_near if delta < sep else ratio_far
    if s1 > 0 and (s2 / s1) < thr:
        info = {"folded_allele": L2, "folded_support": s2, "kept_allele": L1, "kept_support": s1,
                "ratio": round(s2 / s1, 3), "threshold": thr, "neighbouring": bool(delta < sep), "delta": delta}
        return [L1], [s1], info
    return alleles, support, None


def print_histogram(hist):
    by = hist
    print("\n[cigar-len] per-read copy histogram  (SPAN = both flanks = LENGTH; LEFT/RIGHT = one flank "
          "= lower bound / DEPTH; INT = internal)", file=sys.stderr)
    print(f"{'copies':>7} | {'SPAN':>5} {'LEFT':>5} {'RIGHT':>6} {'INT':>5}", file=sys.stderr)
    for c in sorted({c for cnt in by.values() for c in cnt}):
        print(f"{c:>7} | {by['SPAN'][c]:>5} {by['LEFT'][c]:>5} {by['RIGHT'][c]:>6} {by['INT'][c]:>5}"
              f"  {'#'*min(by['SPAN'][c],50)}", file=sys.stderr)
    return sum(by["SPAN"].values())


# ── STAGE 2 : per-allele single-contig re-align + muc1_analyzer call + combined score ──
def available_contigs(ref, template=_CONTIG_TEMPLATE):
    if not os.path.exists(ref + ".fai"):
        pysam.faidx(ref)
    pat = re.compile(re.escape(template).replace(r"\{n\}", r"(\d+)"))
    out = {}
    with open(ref + ".fai") as fh:
        for line in fh:
            m = pat.fullmatch(line.split("\t")[0])
            if m:
                out[int(m.group(1))] = line.split("\t")[0]
    return out


def pick_contig(copies, contigs, contig_offset):
    n = min(contigs, key=lambda k: abs(k - (copies + contig_offset)))
    return n, contigs[n]


def _indel_burden(reads, names, contig_fa, preset):
    """Total inserted+deleted bp when the bin's reads are aligned to `contig_fa`. Minimal when the
    scaffold length matches the reads' true VNTR length: too short → reads' excess units pile up as
    insertions; too long → big deletions. A far more reliable length than the ruler AL→AH count, whose
    anchors sit a few units inside the tandem (the ruler self-measures ~146 for the 150-mer)."""
    aln = mappy.Aligner(contig_fa, preset=preset)
    tot = 0
    for name in names:
        seq = reads[name][0]
        best = None
        for h in aln.map(seq):
            if h.is_primary:
                best = h
                break
            if best is None or h.mlen > best.mlen:
                best = h
        if best:
            tot += sum(length for length, op in best.cigar if op in (1, 2))
    return tot


def parse_clair3_vcf(vcf_path, allele_contigs, vntr_start=4574):
    """Ingest a Clair3 VCF (called against the faked VNTR reference) and extract, per allele contig, the
    first PASS frameshift indel in the tandem — mirroring the reference pipeline's awk. Clair3 is a
    deep-learning ONT caller that models homopolymer error, so it is BOTH sensitive and specific for a
    dupC where a hand-rolled 8C-fraction test cannot be (it correctly calls the 59dupC that our
    statistical caller, held to specificity, misses at low depth). repeat index = (POS-start)//60,
    unit position = (POS-start)%60; a G>GC at unit position 20 is the canonical MUC1 dupC.
    Returns {scaffold_contig: {repeat, unit_pos, ref, alt, label, is_dupc} | None}."""
    import gzip
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    contig_n = {c: L for L, c in allele_contigs.items()}       # contig → scaffold length
    hits = {c: None for c in contig_n}
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            chrom, pos, _id, ref, alt, _qual, filt = f[:7]
            if chrom not in contig_n or filt != "PASS":
                continue
            pos = int(pos)
            nb = contig_n[chrom]
            if not (vntr_start <= pos < vntr_start + 60 * nb):
                continue
            if len(ref) == len(alt):                            # indel only (frameshift)
                continue
            if hits[chrom] is not None:                         # first PASS indel per allele (like awk `exit`)
                continue
            rep = (pos - vntr_start) // 60
            upos = (pos - vntr_start) % 60
            is_dupc = (upos == 20 and ref == "G" and alt == "GC")
            label = "59dupC" if is_dupc else f"{upos}:{ref}>{alt}"
            hits[chrom] = {"repeat": rep, "unit_pos": upos, "ref": ref, "alt": alt,
                           "label": label, "is_dupc": is_dupc, "pos": pos}
    return hits


def assign_reads_to_alleles(reads, allele_contigs, known_allele_of, ref, preset, outdir):
    """Assign every read to an allele, RELIABLY BY CONSTRUCTION. A read whose allele is DETERMINABLE
    (spanning → exact length; or a one-flank read whose lower bound already exceeds the short allele →
    must be the long one) is placed by that KNOWN allele and never second-guessed by alignment score.
    Only GENUINELY AMBIGUOUS reads — partial internal fragments covering just the shared tandem middle,
    whose allele no method can resolve — fall back to the minimum-burden (best-fit) contig. This recovers
    depth at the variant locus without ever mis-placing an allele-informative read. Returns
    ({L: names}, {L: {'confident': n, 'ambiguous': n}})."""
    import mappy
    aligners = {}
    for L, cname in allele_contigs.items():
        cfa = os.path.join(outdir, f"_assign_{L}.fa")
        extract_contig(ref, cname, cfa)
        aligners[L] = mappy.Aligner(cfa, preset=preset)
    out = {L: [] for L in allele_contigs}
    counts = {L: {"confident": 0, "ambiguous": 0} for L in allele_contigs}
    for name, (seq, _) in reads.items():
        kL = known_allele_of.get(name)
        if kL in allele_contigs:                               # allele known → trust it, no burden
            out[kL].append(name)
            counts[kL]["confident"] += 1
            continue
        best_L, best_burden = None, None                       # ambiguous → minimum-burden best fit
        for L, aln in aligners.items():
            hit = None
            for h in aln.map(seq):
                if h.is_primary:
                    hit = h
                    break
                if hit is None or h.mlen > hit.mlen:
                    hit = h
            if hit is None:
                continue
            b = sum(l for l, op in hit.cigar if op in (1, 2)) + hit.q_st + (len(seq) - hit.q_en)
            if best_burden is None or b < best_burden:
                best_burden, best_L = b, L
        if best_L is not None:
            out[best_L].append(name)
            counts[best_L]["ambiguous"] += 1
    for L in allele_contigs:
        for ext in ("", ".fai"):
            try:
                os.remove(os.path.join(outdir, f"_assign_{L}.fa" + ext))
            except OSError:
                pass
    return out, counts


def refine_scaffold(reads, names, contigs, center, ref, preset, outdir, window=6):
    """Pick the scaffold contig that MINIMISES the bin's indel burden, scanning ±window around the
    stage-1 estimate. Returns (n, contig_name, burden_table). This is what makes a 3′ frameshift (e.g.
    59dupC) recoverable: on a scaffold even 2–4 units too short the carrier units are soft-clipped or
    buried in large insertions, and the variant never reaches the consensus."""
    cand = sorted(n for n in contigs if center - window <= n <= center + window)
    if not cand:
        n, name = pick_contig(center, contigs, 0)
        return n, name, []
    burdens = []
    for n in cand:
        cfa = os.path.join(outdir, f"_sweep_{n}.fa")
        extract_contig(ref, contigs[n], cfa)
        burdens.append((_indel_burden(reads, names, cfa, preset), n))
        for ext in ("", ".fai"):
            try:
                os.remove(cfa + ext)
            except OSError:
                pass
    burdens.sort()
    best_n = burdens[0][1]
    return best_n, contigs[best_n], burdens


def extract_contig(ref, name, out_fa):
    with pysam.FastaFile(ref) as fa:
        seq = fa.fetch(name)
    with open(out_fa, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + "\n")
    pysam.faidx(out_fa)
    return out_fa


def realign_bin(reads, names, contig_fa, preset, out_bam):
    aln = mappy.Aligner(contig_fa, preset=preset)
    if not aln:
        raise RuntimeError(f"mappy failed to index {contig_fa}")
    ctg = aln.seq_names[0]
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": ctg, "LN": len(aln.seq(ctg))}]}
    hdr = pysam.AlignmentHeader.from_dict(header)
    tmp = out_bam + ".uns.bam"
    n_ok = n_no = 0
    with pysam.AlignmentFile(tmp, "wb", header=header) as out:
        for name in names:
            seq, quals = reads[name]
            best = None
            for h in aln.map(seq):
                if h.is_primary:
                    best = h
                    break
                if best is None or h.mlen > best.mlen:
                    best = h
            if best is None:
                n_no += 1
                continue
            out.write(_seg(hdr, name, seq, quals, best))
            n_ok += 1
    pysam.sort("-o", out_bam, tmp)
    pysam.index(out_bam)
    os.remove(tmp)
    return n_ok, n_no


def _phasing_concordance(bam, contig, snps, A, B, min_grp_support=5, maj_frac=0.8):
    """Count SNP positions that TRULY discriminate the two phased groups: A's majority base ≠ B's, each
    backed by ≥ maj_frac of its group. ONT systematic errors (homopolymers, etc.) look like SNPs but are
    SHARED by both groups → they don't discriminate. So this count is ~0 for a genuine homozygote (a false
    split on noise) and equals the real differences for a true same-length heterozygote."""
    import collections
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
                (ca if nm in A else cb if nm in B else collections.Counter())[base] += 1
            na, nb = sum(ca.values()), sum(cb.values())
            if na < min_grp_support or nb < min_grp_support:
                continue
            (ba, va), (bb, vb) = ca.most_common(1)[0], cb.most_common(1)[0]
            if ba != bb and va >= maj_frac * na and vb >= maj_frac * nb:
                discrim += 1
    return discrim


def phase_hom_reads(cname, names, reads, ref, preset, outdir, min_snps=3, min_reads=3,
                    min_depth=5, snp_af=0.30):
    """Length-homozygous case: the two alleles have the SAME length and pile on one scaffold, so the
    length can't separate them — split them by SEQUENCE instead, reusing the reference caller's
    find_phasing_snps + phase_reads (which operate on a single contig). A concordance gate then rejects
    a false split driven by ONT noise (true homozygote). Returns (namesA, namesB, n_discriminating_snps)
    when the two groups genuinely differ, else None."""
    try:
        from muc1_analyzer.caller import find_phasing_snps, phase_reads
    except Exception:
        return None
    ref_len = pysam.FastaFile(ref).get_reference_length(cname)
    cfa = os.path.join(outdir, "_phase.fa")
    cbam = os.path.join(outdir, "_phase.bam")
    extract_contig(ref, cname, cfa)
    realign_bin(reads, names, cfa, preset, cbam)
    result = None
    try:
        snps = find_phasing_snps(cbam, cname, ref_len, min_depth=max(min_depth, 5),
                                 min_af=snp_af, max_af=1.0 - snp_af)
        if snps:
            rA, rB = phase_reads(cbam, cname, ref_len, snps)
            A = [n for n in rA if n in names]
            B = [n for n in rB if n in names]
            if len(A) >= min_reads and len(B) >= min_reads:
                discrim = _phasing_concordance(cbam, cname, snps, A, B)
                if discrim >= min_snps:                        # the two groups really differ → real het
                    result = (A, B, discrim)
    except Exception as e:
        print(f"[cigar-len] phasing skipped ({e})", file=sys.stderr)
    for pth, exts in ((cfa, ("", ".fai")), (cbam, ("", ".bai"))):
        for e in exts:
            try:
                os.remove(pth + e)
            except OSError:
                pass
    return result


def _consensus_snp_diffs(hapA, hapB):
    """Number of SNP (nucleotide) differences between the two haplotype consensuses. The two motif
    sequences are aligned by motif order (difflib, so a length difference of one repeat is an INDEL, not
    counted), then aligned motif pairs are compared base by base. This is the clinician's criterion — the
    number of SNPs between the alleles — read off the reconstructed consensus (robust: per-read errors are
    averaged out). ≤ threshold → the two 'alleles' are effectively identical → homozygote."""
    import difflib
    ma, mb = hapA.get("motifs") or [], hapB.get("motifs") or []
    na = [m.get("name", "?") for m in ma]
    nb = [m.get("name", "?") for m in mb]
    snps = 0
    for op, i1, i2, j1, j2 in difflib.SequenceMatcher(None, na, nb, autojunk=False).get_opcodes():
        if op in ("equal", "replace"):
            for k in range(min(i2 - i1, j2 - j1)):
                sa = ma[i1 + k].get("sequence", "") or ""
                sb = mb[j1 + k].get("sequence", "") or ""
                L = min(len(sa), len(sb))
                snps += sum(1 for x in range(L) if sa[x] != sb[x]) + abs(len(sa) - len(sb))
        # 'insert'/'delete' = a whole extra/missing repeat (length difference) → NOT a SNP, not counted
    return snps


def call_unit(u, reads, ref, outdir, preset, min_depth, ins_threshold, keep_pdf):
    """Realign one haplotype unit's reads to its scaffold and run `muc1_analyzer call`. Returns everything
    the caller needs to record the haplotype (hap dict, summary bin, consensus/novel paths, intermediates)."""
    cname, names, meta, tag, scaffold_L = u["cname"], u["names"], u["meta"], u["tag"], u["scaffold_L"]
    cfa = os.path.join(outdir, f"{tag}.contig.fa")
    cbam = os.path.join(outdir, f"{tag}.bam")
    extract_contig(ref, cname, cfa)
    n_ok, n_no = realign_bin(reads, names, cfa, preset, cbam)
    cjson, cpdf, rc, cons_fa, novel_tsv = run_call(cbam, cfa, tag, outdir, min_depth=min_depth,
                                                   ins_threshold=ins_threshold, want_pdf=keep_pdf)
    hap = top_haplotype(cjson)
    if hap is not None:
        hap.setdefault("contig", cname)
        hap["tag"] = tag
        hap["hap_label"] = u["hap"]
    bin_dict = {"allele_copies": meta["stage1_allele"], "scaffold_copies": scaffold_L,
                "stage1_allele": meta["stage1_allele"], "bin_median_copies": meta["bin_median"],
                "scaffold_sweep_burden": meta["burden"], "scaffold_contig": cname,
                "tag": tag, "hap_label": u["hap"], "reads_binned": len(names), "reads_aligned": n_ok,
                "bam": cbam, "analyzer_json": cjson, "call_returncode": rc}
    return {"hap": hap, "bin": bin_dict, "cons_fa": cons_fa, "novel_tsv": novel_tsv, "n_ok": n_ok, "n_no": n_no,
            "rc": rc, "cname": cname, "scaffold_L": scaffold_L, "hap_label": u["hap"], "n_names": len(names),
            "ints": [cfa, cfa + ".fai", cbam, cbam + ".bai", cjson, cpdf, cons_fa, novel_tsv]}


def run_call(bam, contig_fa, tag, outdir, min_depth=3, ins_threshold=0.5,
             fasta_consensus=True, novel_motifs=True, want_pdf=False):
    cjson = os.path.join(outdir, f"{tag}.analyzer.json")
    cpdf = os.path.join(outdir, f"{tag}.report.pdf") if want_pdf else None
    cmd = [sys.executable, "-m", "muc1_analyzer", "call", "-b", bam, "-r", contig_fa,
           "-s", tag, "--min-mq", "0", "--min-depth", str(min_depth),
           "--ins-threshold", str(ins_threshold), "--json", cjson]
    if cpdf:
        cmd += ["--pdf", cpdf]
    cfa = os.path.join(outdir, f"{tag}.consensus.fa") if fasta_consensus else None
    ctsv = os.path.join(outdir, f"{tag}.novel_motifs.tsv") if novel_motifs else None
    if cfa:
        cmd += ["--fasta-consensus", cfa]
    if ctsv:
        cmd += ["--novel-motifs", ctsv]
    rc = subprocess.run(cmd, cwd=REPO_ROOT)
    return cjson, cpdf, rc.returncode, cfa, ctsv


def merge_allele_bams(bins, out_bam):
    """Combine the per-allele BAMs into ONE coordinate-sorted, indexed BAM holding BOTH haplotypes — each
    read still on its own allele's scaffold contig (header carries both contigs). Handy for IGV: load it
    against the faked reference and switch between the two length-contigs. Returns the path or None."""
    entries = [b for b in sorted(bins, key=lambda b: b["scaffold_copies"]) if b.get("bam") and os.path.exists(b["bam"])]
    if not entries:
        return None
    sq, seen = [], set()
    for b in entries:
        with pysam.AlignmentFile(b["bam"]) as af:
            for s in af.header.to_dict().get("SQ", []):
                if s["SN"] not in seen:
                    seen.add(s["SN"])
                    sq.append({"SN": s["SN"], "LN": s["LN"]})
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}
    H = pysam.AlignmentHeader.from_dict(header)
    tmp = out_bam + ".unsorted.bam"
    with pysam.AlignmentFile(tmp, "wb", header=header) as out:
        for b in entries:
            with pysam.AlignmentFile(b["bam"]) as af:
                for r in af.fetch(until_eof=True):
                    out.write(pysam.AlignedSegment.from_dict(r.to_dict(), H))   # remap ref to merged header
    pysam.sort("-o", out_bam, tmp)
    pysam.index(out_bam)
    os.remove(tmp)
    return out_bam


def merge_consensus_fasta(fastas, out_path):
    """Combine the per-allele consensus FASTAs into one file (both haplotypes, in length order)."""
    kept = []
    with open(out_path, "w") as out:
        for fa in fastas:
            if fa and os.path.exists(fa):
                out.write(open(fa).read().rstrip("\n") + "\n")
                kept.append(fa)
    return out_path if kept else None


def merge_novel_motifs(tsvs, out_path, sample=""):
    """Merge the per-allele novel-motif TSVs, re-aggregating by sequence across both haplotypes (sum
    occurrences, union haplotypes/positions, keep the closest known motif) — the same cross-haplotype
    semantics as the reference caller's `collect_novel_60bp_motifs`, but assembled from the per-allele runs."""
    header = ["sample", "occurrences", "haplotypes", "motif_positions",
              "closest_motif", "closest_dist", "diff_positions", "sequence"]
    agg = {}
    for tsv in tsvs:
        if not (tsv and os.path.exists(tsv)):
            continue
        with open(tsv) as fh:
            rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
        if not rows or rows[0][0] == "sample":
            rows = rows[1:]
        for r in rows:
            if len(r) < 8:
                continue
            smp, occ, haps, pos, cm, cd, diff, seq = r[:8]
            e = agg.setdefault(seq, {"occ": 0, "haps": [], "pos": [], "cm": cm, "cd": int(cd), "diff": diff})
            e["occ"] += int(occ)
            e["haps"] += [h for h in haps.split(";") if h and h not in e["haps"]]
            e["pos"] += [p for p in pos.split(";") if p]
            if int(cd) < e["cd"]:
                e["cd"], e["cm"], e["diff"] = int(cd), cm, diff
    cands = sorted(agg.items(), key=lambda kv: (-kv[1]["occ"], kv[1]["cd"]))
    with open(out_path, "w") as out:
        out.write("\t".join(header) + "\n")
        for seq, e in cands:
            out.write("\t".join([sample, str(e["occ"]), ";".join(e["haps"]), ";".join(e["pos"]),
                                 e["cm"], str(e["cd"]), e["diff"], seq]) + "\n")
    return out_path, len(cands)


def top_haplotype(cjson):
    """The most-supported haplotype from a per-allele call JSON (guards against the caller splitting a
    single-length bin into noise haplotypes: we keep the densest one)."""
    if not os.path.exists(cjson):
        return None
    haps = json.load(open(cjson)).get("haplotypes", [])
    return max(haps, key=lambda h: h.get("read_count", 0)) if haps else None


def _variant_short(motif):
    """'X-59dupC' → '59dupC' (drop the leading motif-letter family prefix) for the clinical string."""
    return re.sub(r"^[A-Za-z0-9]+-", "", motif or "") or (motif or "?")


def clinical_string(frag, haps=None):
    """'36 repeat | 70 repeat (59dupC ~ repeat 20)' — short|long, parenthetical on the MUTANT allele.
    Both allele lengths come from the haplotypes when given (so they always show, carrier or not)."""
    from muc1_analyzer.detectors.vntr import _contig_repeat_count
    if haps:
        lens = sorted(n for n in (_contig_repeat_count(h.get("contig", "")) for h in haps) if n)
    else:
        lens = sorted(x for x in (frag.get("vntr_len_healthy"), frag.get("vntr_len_mut")) if x is not None)
    if not lens:
        return "?"
    mv, mut_len = frag.get("mutation_variant"), frag.get("vntr_len_mut")
    ann = f" ({_variant_short(mv.get('motif'))} ~ repeat {mv.get('repeat_index')})" if mv else ""
    parts, annotated = [], False
    for L in lens:
        if mv and L == mut_len and not annotated:        # only ONE allele carries it (same-length safe)
            parts.append(f"{L} repeat{ann}")
            annotated = True
        else:
            parts.append(f"{L} repeat")
    return " | ".join(parts)


def _cohort_ref():
    """Cohort calibration of the onset positioning — mean/std of frac_to_mut (= position/mut_len) from
    the n=27 survey stored in muc1_analyzer.config.SEVERITY_CAL. Falls back to those values if the import
    fails so the plot still renders."""
    try:
        from muc1_analyzer.config import SEVERITY_CAL as S
        return float(S["mean_frac"]), float(S["std_frac"]), 27
    except Exception:
        return 0.4625, 0.2061, 27


def _cohort_plot(onset_index, mean_oi, std_oi, z):
    """A blue Gaussian of the calibration cohort on the onset_index axis (0→1) with a red dot at the
    tested sample. onset_index = fraction of the mutant VNTR downstream of the frameshift (translated as
    the MUC1fs neoprotein); mean_oi = 1 − mean_frac."""
    from reportlab.graphics.shapes import Drawing, Line, PolyLine, Circle, String
    from reportlab.lib import colors
    import math
    W, Hd, lx, rx, by, ty = 470, 150, 45, 440, 48, 126
    d = Drawing(W, Hd)
    X = lambda oi: lx + max(0.0, min(1.0, oi)) * (rx - lx)
    G = lambda oi: math.exp(-0.5 * ((oi - mean_oi) / std_oi) ** 2)
    Y = lambda g: by + g * (ty - by)
    pts = []
    for i in range(81):
        oi = i / 80.0
        pts += [X(oi), Y(G(oi))]
    d.add(PolyLine(pts, strokeColor=colors.HexColor("#185FA5"), strokeWidth=1.8))
    d.add(Line(lx, by, rx, by, strokeColor=colors.HexColor("#888"), strokeWidth=1))
    for oi, lab in [(0.0, "0"), (mean_oi, "cohort mean"), (1.0, "1")]:
        x = X(oi)
        d.add(Line(x, by, x, by - 4, strokeColor=colors.HexColor("#888")))
        d.add(String(x, by - 14, lab, fontName="Helvetica", fontSize=7,
                     fillColor=colors.HexColor("#666"), textAnchor="middle"))
    # direction of effect: onset_index↑ (right) = earlier frameshift → more neoprotein → earlier onset
    d.add(String(lx, by - 26, "\u25c4 later onset", fontName="Helvetica-Oblique", fontSize=7,
                 fillColor=colors.HexColor("#888"), textAnchor="start"))
    d.add(String(rx, by - 26, "earlier onset \u25ba", fontName="Helvetica-Oblique", fontSize=7,
                 fillColor=colors.HexColor("#888"), textAnchor="end"))
    xm = X(mean_oi)
    d.add(Line(xm, by, xm, Y(G(mean_oi)), strokeColor=colors.HexColor("#9BB6D6"),
              strokeWidth=0.8, strokeDashArray=[2, 2]))
    xs, ys = X(onset_index), Y(G(onset_index))
    d.add(Line(xs, by, xs, ys, strokeColor=colors.red, strokeWidth=0.8, strokeDashArray=[2, 2]))
    d.add(Circle(xs, ys, 4, fillColor=colors.red, strokeColor=colors.red))
    d.add(String(xs, ys + 8, f"sample (z = {z:+.2f})", fontName="Helvetica-Bold", fontSize=7.5,
                 fillColor=colors.red, textAnchor="middle"))
    d.add(String((lx + rx) / 2, Hd - 9, "Onset position in the calibration cohort (n=27)",
                 fontName="Helvetica", fontSize=8, fillColor=colors.HexColor("#444"), textAnchor="middle"))
    d.add(String((lx + rx) / 2, 6, "onset_index = fraction of the mutant VNTR translated into neoprotein",
                 fontName="Helvetica-Oblique", fontSize=6.5, fillColor=colors.HexColor("#999"), textAnchor="middle"))
    return d


def write_merged_pdf(path, sample, frag, score, haps, rel_dupc=None, rs_detected=None, rs_alert=None):
    """One PDF for both haplotypes — DELEGATED to `caller.write_pdf_report`, the single renderer.

    This used to be a second, parallel PDF writer. Two writers meant a clinician could be handed two
    different-looking reports for the same sample, and a field fixed on one side stayed wrong on the
    other. Its distinctive blocks (the blue clinical-call box, the pathogenic motif in red, the
    power-aware dupC wording, the cohort onset plot) now live in `muc1_analyzer.report_blocks` and are
    rendered by the shared writer, which additionally carries the arbiter length, the T2T-native
    rs4072037 and the LAYER CONGRUENCE block this one never had.

    Falls back to the legacy body below only if the shared renderer is unavailable."""
    try:
        from muc1_analyzer.caller import write_pdf_report
        from muc1_analyzer.detectors.vntr import _contig_repeat_count
    except Exception:
        write_pdf_report = None
    if write_pdf_report is not None:
        var = (rel_dupc or {}).get("variant") or {}
        carrier_contig = (rel_dupc or {}).get("carrier_contig")
        results = []
        for i, h in enumerate(sorted(haps, key=lambda x: _contig_repeat_count(x.get("contig", "")) or 0), 1):
            n = _contig_repeat_count(h.get("contig", "")) or 0
            motifs = h.get("motifs") or []
            results.append({"rank": i, "contig": h.get("contig", "?"),
                            "read_count": h.get("read_count", 0), "contig_length": n * 60,
                            "n_motifs": len(motifs) or n,
                            "n_exact": sum(1 for m in motifs if not m.get("mismatches")),
                            "n_approx": sum(1 for m in motifs if m.get("mismatches")),
                            "motifs": motifs})
        rs_block = None
        if rs_detected and rs_detected.get("base"):
            rs_block = {"available": True, "genotype": rs_detected["base"],
                        "note": f"read off the mutant allele (C:{rs_detected.get('C')} / "
                                f"T:{rs_detected.get('T')}, depth {rs_detected.get('depth')})"}
        try:
            write_pdf_report(results, path, sample=sample, dupc={"result": rel_dupc} if rel_dupc else None,
                             rs4072037=rs_block,
                             clinical_call=clinical_string(frag, haps),
                             onset_index=score.get("onset_index"),
                             variant_repeat=var.get("repeat"), variant_label=var.get("label"),
                             carrier_contig=carrier_contig,
                             score=score, rs_alert=rs_alert)
            if rs_alert:
                print(f"[cigar-len] {rs_alert}", file=sys.stderr)
            return path
        except Exception as e:
            print(f"[cigar-len] shared renderer failed ({e}) — falling back to the legacy writer",
                  file=sys.stderr)
    # ── legacy fallback ──────────────────────────────────────────────────────────
    try:
        from reportlab.lib.pagesizes import A4
        from reportlab.lib import colors
        from reportlab.lib.units import cm
        from reportlab.lib.styles import ParagraphStyle
        from reportlab.lib.enums import TA_LEFT
        from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, HRFlowable
    except ImportError:
        print("[cigar-len] reportlab unavailable — skipping merged PDF.", file=sys.stderr)
        return None
    H1 = ParagraphStyle("h1", fontName="Helvetica-Bold", fontSize=15, spaceAfter=8)
    BOX = ParagraphStyle("box", fontName="Helvetica-Bold", fontSize=13, textColor=colors.HexColor("#0C447C"),
                         leading=17)
    BODY = ParagraphStyle("body", fontName="Helvetica", fontSize=10.5, leading=15)
    META = ParagraphStyle("meta", fontName="Helvetica", fontSize=9.5, leading=13, textColor=colors.HexColor("#444"))
    GENO = ParagraphStyle("geno", fontName="Helvetica-Bold", fontSize=13, leading=17, spaceBefore=4, alignment=TA_LEFT)

    carrier = bool(frag.get("mutation_present"))
    call_str = clinical_string(frag, haps)
    story = [Paragraph(f"MUC1 — merged result · {sample}", H1)]

    # ── blue reliable-length box (clinical call string) ───────────────────────────
    box = Paragraph(f"VNTR — clinical call&nbsp;: {call_str}", BOX)
    tbl = Table([[box]], colWidths=[17 * cm])
    tbl.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#E6F1FB")),
                             ("BOX", (0, 0), (-1, -1), 1.2, colors.HexColor("#185FA5")),
                             ("LEFTPADDING", (0, 0), (-1, -1), 10), ("RIGHTPADDING", (0, 0), (-1, -1), 10),
                             ("TOPPADDING", (0, 0), (-1, -1), 9), ("BOTTOMPADDING", (0, 0), (-1, -1), 9)]))
    story += [tbl, Spacer(1, 0.35 * cm)]

    # ── pathogenic variant presence + MUC1_Score ─────────────────────────────────
    if carrier and frag.get("mutation_variant"):
        mv = frag["mutation_variant"]
        which = "long" if score.get("mut_is_long") else "short"
        pres = (f"<b>Pathogenic variant (drives the MUC1_Score)&nbsp;: PRESENT</b> — "
                f"{_variant_short(mv.get('motif'))} on the {which} allele "
                f"({frag.get('vntr_len_mut')} rep.), repeat {mv.get('repeat_index')}"
                f"{'' if mv.get('known_pathogenic') else ' (not catalogued)'}")
    else:
        pres = "<b>Pathogenic variant (drives the MUC1_Score)&nbsp;: ABSENT</b> — no focal frameshift detected"
    story += [Paragraph(pres, BODY), Spacer(1, 0.2 * cm)]

    # trustworthy dupC verdict (specificity-calibrated caller), separate from the depth-fragile token
    # trustworthy dupC verdict — delegated to the genome-anchored caller (separate from the fragile token)
    if rel_dupc is not None:
        if not rel_dupc.get("assessed"):
            line = "dupC (reliable verdict)&nbsp;: <b>NOT ASSESSED</b> — " + str(rel_dupc.get("reason", ""))
        else:
            interp = rel_dupc.get("interpretation")
            pw = rel_dupc.get("power_by_allele", {})
            cov = ", ".join(f"{c}: {t['n_cover']} covering reads"
                            f"{' (powered)' if t['adequately_powered'] else ' (underpowered)'}"
                            for c, t in pw.items())
            if interp == "CONFIRMED":
                meth = {"clair3": "Clair3", "runlen_shift": "homopolymer run-length shift, calibrated specificity"}.get(
                    rel_dupc.get("method"), "context+paired, calibrated specificity")
                var = rel_dupc.get("variant", {})
                extra = f", repeat {var.get('repeat')}" if var.get("repeat") is not None else ""
                line = (f"<b>dupC&nbsp;: CONFIRMED</b> on {rel_dupc.get('carrier_contig')}{extra} ({meth})")
            elif interp == "NEGATIVE":
                line = (f"<b>dupC&nbsp;: NEGATIVE (reliable)</b> — not detected, with adequate power "
                        f"[{cov}]")
            elif interp == "NEGATIVE_BORDERLINE":
                line = (f"<b>dupC&nbsp;: NEGATIVE (borderline power)</b> — not detected&nbsp;; coverage at the edge "
                        f"of detectability [{cov}]&nbsp;: a low-d carrier could be missed")
            else:
                need = max((t["reads_for_power0.90"].get("d=0.3", 12) for t in pw.values()), default=12)
                line = (f"<b>dupC&nbsp;: INDETERMINATE</b> — not detected but coverage underpowered, "
                        f"a dupC cannot be excluded [{cov}]&nbsp;; ~{need} reads/allele required (d=0.3)")
        story.append(Paragraph(line, META))
        story.append(Spacer(1, 0.2 * cm))

    story += [HRFlowable(width="100%", color=colors.HexColor("#CCC")),
              Spacer(1, 0.25 * cm), Paragraph("Haplotypes", H1)]

    # ── the two haplotypes (length + reads + nomenclature; variant motif in red as X(59dupC)) ──────
    from muc1_analyzer.detectors.vntr import _contig_repeat_count
    carrier = (rel_dupc or {}).get("carrier_tag") or (rel_dupc or {}).get("carrier_contig")
    var = (rel_dupc or {}).get("variant") or {}
    var_rep, var_lbl = var.get("repeat"), var.get("label")
    for h in sorted(haps, key=lambda x: _contig_repeat_count(x.get("contig", "")) or 0):
        n = _contig_repeat_count(h.get("contig", ""))
        hlab = f"  ·  hap {h.get('hap_label')}" if h.get("hap_label") else ""
        story.append(Paragraph(f"<b>{n} repeats</b>{hlab} — {h.get('contig')} "
                               f"({h.get('read_count', '?')} reads)", BODY))
        # rebuild the nomenclature from the motif list so we can recolour the variant motif exactly;
        # the variant repeat is a 1-based motif index, so motif #var_rep gets "<name>(label)" in red.
        motifs = h.get("motifs") or []
        is_carrier = (h.get("tag") == carrier) or (h.get("hap_label") is None and h.get("contig") == carrier)
        if motifs:
            parts = []
            for idx, m in enumerate(motifs):
                name = m.get("name", "?")
                if is_carrier and var_rep is not None and idx + 1 == var_rep:
                    # the consensus may ALREADY name the motif with the frameshift (e.g. "X-59dupC");
                    # only append "(label)" when it doesn't, so the variant is never written twice.
                    disp = name if (var_lbl and var_lbl in name) else f"{name}({var_lbl})"
                    parts.append(f'<font color="#C00000"><b>{disp}</b></font>')
                else:
                    parts.append(name)
            nom = "-".join(parts)
        else:
            nom = h.get("nomenclature", "")
        if nom:
            story.append(Paragraph(nom, META))
        story.append(Spacer(1, 0.2 * cm))

    # ── SECTION 2 — MUC1_Score (values + cohort position + genotype phenotype flag) ──────────────
    story += [Spacer(1, 0.15 * cm), HRFlowable(width="100%", color=colors.HexColor("#CCC")),
              Spacer(1, 0.25 * cm), Paragraph("MUC1_Score", H1)]
    sev = score.get("severity_score")
    if rs_detected and rs_detected.get("base"):
        src = f"read off the mutant allele — C:{rs_detected['C']} / T:{rs_detected['T']}, depth {rs_detected['depth']}"
    elif rs_detected and rs_detected.get("label") == "low_depth":
        src = f"insufficient depth (C:{rs_detected['C']} / T:{rs_detected['T']})"
    else:
        src = "provided"
    sev_txt = (f"{sev:+g} (rs4072037-{score.get('splice_base')} → "
               f"{'severe' if sev and sev > 0 else 'protective'} ; {src})") if sev is not None \
        else "n/a — rs4072037 undetermined (no usable read and no --rs4072037)"
    rows = [["MUC1_Score — ONSET axis", f"{score.get('onset_score')}  (onset_index {score.get('onset_index')})"],
            ["MUC1_Score — SEVERITY axis", sev_txt],
            ["VNTR ratio (mut/healthy)", f"{score.get('ratio')}"]]
    st = Table(rows, colWidths=[6 * cm, 11 * cm])
    st.setStyle(TableStyle([("FONT", (0, 0), (-1, -1), "Helvetica", 9.5),
                            ("FONT", (0, 0), (0, -1), "Helvetica-Bold", 9.5),
                            ("LINEBELOW", (0, 0), (-1, -2), 0.4, colors.HexColor("#DDD")),
                            ("TOPPADDING", (0, 0), (-1, -1), 4), ("BOTTOMPADDING", (0, 0), (-1, -1), 4)]))
    story += [st, Spacer(1, 0.3 * cm)]

    # cohort calibration curve + red sample dot (only meaningful when there is an onset position)
    oi = score.get("onset_index")
    if oi is not None:
        mean_frac, std_frac, n = _cohort_ref()
        mean_oi, z = 1.0 - mean_frac, (oi - (1.0 - mean_frac)) / std_frac
        try:
            story += [_cohort_plot(oi, mean_oi, std_frac, z), Spacer(1, 0.25 * cm)]
        except Exception as e:
            story.append(Paragraph(f"<i>cohort curve unavailable ({e})</i>", META))

    # genotype → phenotype flag, driven by the severity axis
    if sev is not None and sev > 0:
        story.append(Paragraph('<font color="#C00000"><b>GENOTYPE ASSOCIATED WITH A SEVERE PHENOTYPE</b></font>', GENO))
    elif sev is not None and sev < 0:
        story.append(Paragraph('<font color="#1B7A2F"><b>GENOTYPE ASSOCIATED WITH A PROTECTIVE PHENOTYPE</b></font>', GENO))
    else:
        story.append(Paragraph('<font color="#888"><b>SEVERITY UNDETERMINED — rs4072037 NOT DETERMINED</b></font>', GENO))

    if rs_alert:
        story.append(Spacer(1, 0.15 * cm))
        story.append(Paragraph(f'<font color="#B8860B"><b>⚠ {rs_alert}</b></font>', META))

    SimpleDocTemplate(path, pagesize=A4, topMargin=1.4 * cm, bottomMargin=1.4 * cm,
                      leftMargin=1.6 * cm, rightMargin=1.6 * cm).build(story)
    return path


def _motif_index_at(analyzer_json, refpos):
    """1-based index of the motif in the consensus nomenclature that spans `refpos` — i.e. WHICH repeat
    carries the variant, counted from the first motif exactly like the total VNTR length. This is the
    clinically meaningful repeat number (consistent with the analyzer and with the allele length), not a
    genomic-coordinate count from an arbitrary origin."""
    try:
        haps = json.load(open(analyzer_json)).get("haplotypes", [])
        motifs = haps[0].get("motifs", []) if haps else []
        for idx, m in enumerate(motifs):
            p, L = m.get("pos"), m.get("length", 60)
            if p is not None and p <= refpos < p + L:
                return idx + 1
    except Exception:
        return None
    return None


def runlen_shift_dupc(bins, ref, vntr_start=4574, repeat_offset=4, mut_len=8, wt_len=7,
                      alpha=0.001, min_k=4, min_frac=0.25, min_reads=6):
    """General, self-calibrated homopolymer-frameshift caller — the right statistic for a dupC.

    Why not count reads showing exactly `mut_len` C's: ONT miscounts a homopolymer, so a TRUE 8C run
    reads as 6/7/8/9. Counting `==8` keeps only ~1/3 of a carrier's reads (the 3/9 puzzle) and discards
    the 9's (over-called 8's). Instead we score the fraction of reads whose C-tract is `>= mut_len`, and
    we test it against an IN-SAMPLE EMPIRICAL NULL: how this very sample's ONT reads its true `wt_len`
    runs (measured over every other C-run) — so the test adapts to the patient/basecaller and needs no
    PoN. One-sided binomial per C-run position, Bonferroni over positions. Detects ANY C-run expansion
    (not just the canonical 59dupC), unlike a Clair3 signature match. Validated on a real carrier:
    calls the dupC at p_bonf ~5e-6 (6/12 reads ≥8 vs a 0.02 null) where the `==8` test gave p_bonf 0.3."""
    from muc1_analyzer.dupc_power import _tail
    import re as _re
    import collections as _c
    fa = pysam.FastaFile(ref)
    out = {"assessed": True, "method": "runlen_shift", "called": False, "per_allele": {}}
    best_overall = None
    for b in bins:
        contig = b["scaffold_contig"]
        key = b.get("tag", contig)                             # unique per haplotype (2 phased haps share a contig)
        hap = b.get("hap_label")
        refseq = fa.fetch(contig).upper()
        runpos = [m.start() for m in _re.finditer("C{%d}A" % wt_len, refseq)]
        per = _c.defaultdict(list)
        with pysam.AlignmentFile(b["bam"]) as af:
            for r in af.fetch(contig):
                if r.is_unmapped or not r.query_sequence:
                    continue
                r2q = {rp: qp for qp, rp in r.get_aligned_pairs() if rp is not None and qp is not None}
                s = r.query_sequence
                for idx, rs in enumerate(runpos):
                    q = r2q.get(rs)
                    if q is None:
                        continue
                    i = q
                    while i < len(s) and s[i] == "C":
                        i += 1
                    j = q
                    while j > 0 and s[j - 1] == "C":
                        j -= 1
                    rl = i - j
                    if rl >= wt_len - 3:                       # drop alignment artifacts (run not captured)
                        per[idx].append(rl)
        allrl = [x for v in per.values() for x in v]
        if not allrl:
            out["per_allele"][key] = {"tested": 0, "note": "no C-runs measured", "contig": contig, "hap": hap}
            continue
        f_null = max(sum(x >= mut_len for x in allrl) / len(allrl), 1e-4)
        tested = [i for i in per if len(per[i]) >= min_reads]
        scored = []
        for i in tested:
            v = per[i]
            k = sum(x >= mut_len for x in v)
            p = _tail(len(v), k, f_null)
            # repeat number = the 1-based MOTIF INDEX carrying the variant (consistent with the length),
            # from the per-allele call's motif list; fall back to the genomic count only if unavailable.
            midx = _motif_index_at(b.get("analyzer_json"), runpos[i])
            repeat = midx if midx is not None else (runpos[i] - vntr_start) // 60 + repeat_offset
            scored.append({"run_index": i, "repeat": repeat,
                           "k_ge_mut": k, "n": len(v), "frac": round(k / len(v), 3),
                           "p": p, "p_bonf": min(1.0, p * max(1, len(tested)))})
        scored.sort(key=lambda d: d["p_bonf"])
        cand = scored[0] if scored else None
        called = bool(cand and cand["p_bonf"] <= alpha and cand["k_ge_mut"] >= min_k
                      and cand["frac"] >= min_frac)
        out["per_allele"][key] = {"f_null": round(f_null, 3), "n_tested": len(tested),
                                  "candidate": cand, "called": called, "contig": contig, "hap": hap}
        if called and (best_overall is None or cand["p_bonf"] < best_overall[1]["candidate"]["p_bonf"]):
            best_overall = (key, out["per_allele"][key], contig, hap)
    if best_overall:
        key, info, contig, hap = best_overall
        out["called"] = True
        out["carrier_tag"] = key
        out["carrier_contig"] = contig
        out["carrier_hap"] = hap
        out["variant"] = {"repeat": info["candidate"]["repeat"], "label": "59dupC",
                          "k_ge_mut": info["candidate"]["k_ge_mut"], "n": info["candidate"]["n"],
                          "p_bonf": info["candidate"]["p_bonf"]}
    out["interpretation"] = "CONFIRMED" if out["called"] else "NEGATIVE"
    return out


def reliable_dupc_context(reads, bins, ruler_fa, ruler_contig, preset, outdir):
    """Specificity-controlled dupC via the CONTEXT+PAIRED caller, self-contained (no external genome).

    The insight that makes this work without T2T: the context caller (`scan_dupc_context`) needs only a
    common reference to fetch reads over + HP tags to contrast the two alleles — and we already have both.
    We align every binned read to the RULER (the 150-mer common reference, so partial reads that merely
    cover the dupC context contribute — no full-span requirement), carry each read's LENGTH bin over as its
    HP tag (bin→HP1/HP2 = the phasing the paired test needs), and run `call_dupc` in paired mode. This beats
    the ~240× multiple-testing (context grouping) AND uses partial reads AND controls specificity in-sample
    (mutant vs its own healthy allele). Returns the caller verdict. Needs the 2-allele (HET) case for the
    pairing; HOM or single-bin → abstains (no internal healthy control)."""
    binlist = sorted(bins, key=lambda b: b["scaffold_copies"])
    if len(binlist) < 2:
        return {"assessed": False, "reason": "single allele (HOM) — no internal healthy control for the paired dupC test"}
    try:
        import mappy
        from muc1_analyzer.dupc_caller import call_dupc
    except Exception as e:
        return {"assessed": False, "reason": f"dupC context caller unavailable: {e}"}
    hp = {}
    for i, b in enumerate(binlist[:2]):                      # shortest→HP1, longest→HP2 (labels only)
        with pysam.AlignmentFile(b["bam"]) as af:
            for r in af.fetch(until_eof=True):
                hp.setdefault(r.query_name, str(i + 1))
    aln = mappy.Aligner(ruler_fa, preset=preset)
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"},
           "SQ": [{"SN": ruler_contig, "LN": len(aln.seq(ruler_contig))}]}
    H = pysam.AlignmentHeader.from_dict(hdr)
    tmp = os.path.join(outdir, "_dupc_ctx.uns.bam")
    with pysam.AlignmentFile(tmp, "wb", header=hdr) as out:
        for name, tag in hp.items():
            if name not in reads:
                continue
            seq, q = reads[name]
            best = None
            for h in aln.map(seq):
                if h.is_primary:
                    best = h
                    break
                if best is None or h.mlen > best.mlen:
                    best = h
            if not best:
                continue
            seg = _seg(H, name, seq, q, best)
            seg.set_tag("HP", int(tag))
            out.write(seg)
    bam = os.path.join(outdir, "_dupc_ctx.bam")
    pysam.sort("-o", bam, tmp)
    pysam.index(bam)
    os.remove(tmp)
    region = f"{ruler_contig}:1-{len(aln.seq(ruler_contig))}"
    res = call_dupc(bam, region=region, hp_mut=None, require_paired=True, mut_len=8, min_support=2)
    # which HP (allele) carries it, if called
    best = res.get("best") or {}
    hp_carrier = best.get("hp") or (best.get("hp_8C") if isinstance(best, dict) else None)
    carrier_contig = None
    if hp_carrier in ("1", "2"):
        carrier_contig = binlist[int(hp_carrier) - 1]["scaffold_contig"]
    return {"assessed": True, "method": "context+paired (ruler-anchored, length-HP)",
            "called": bool(res.get("called")), "tier": res.get("tier"),
            "carrier_contig": carrier_contig, "best": best, "region_reads": len(hp)}


def _tandem_core_depth(bam, contig):
    """Median depth over the central 40–60% of the scaffold — where a 59dupC-type frameshift sits — as the
    n that drives detection power. Central (not 25–75%): the tandem ends are sparse in AS and a mid-tandem
    dupC is what we power for."""
    with pysam.AlignmentFile(bam) as af:
        L = af.get_reference_length(contig)
        lo, hi = int(0.40 * L), int(0.60 * L)
        try:
            cov = af.count_coverage(contig, lo, hi, quality_threshold=0)
        except Exception:
            return 0
    depths = [cov[0][i] + cov[1][i] + cov[2][i] + cov[3][i] for i in range(hi - lo)]
    depths.sort()
    return depths[len(depths) // 2] if depths else 0


def dupc_power_tier(n):
    """Turn a mutant-allele depth `n` into a detection-power verdict, so a 'not called' is read correctly.
    Power depends on the per-read 8C rate `d` (patient/basecaller-specific, ~0.16 pessimistic … 0.50
    optimistic); we report the range and tier on the central estimate d=0.30: ADEQUATE (power ≥ 0.90,
    n ≳ 12), BORDERLINE (0.75–0.90), POOR (< 0.75)."""
    try:
        from muc1_analyzer.dupc_power import power, crit_k, min_reads
    except Exception:
        return None
    p30 = power(n, 0.30)
    tier = "adequate" if p30 >= 0.90 else ("borderline" if p30 >= 0.75 else "poor")
    return {"n_cover": n, "crit_k": crit_k(n), "tier": tier,
            "power": {f"d={d}": round(power(n, d), 2) for d in (0.16, 0.30, 0.50)},
            "adequately_powered": p30 >= 0.90,
            "reads_for_power0.90": {f"d={d}": min_reads(d, target=0.90) for d in (0.16, 0.30, 0.50)}}


def dupc_interpretation(rel_dupc, bins):
    """Attach per-allele detection power to the dupC verdict and derive an INTERPRETATION:
    CONFIRMED / NEGATIVE (reliable) / NEGATIVE (borderline power) / INDETERMINATE (underpowered)."""
    if not rel_dupc.get("assessed"):
        return rel_dupc
    per_allele, tiers = {}, []
    for b in bins:
        n = _tandem_core_depth(b["bam"], b["scaffold_contig"])
        t = dupc_power_tier(n)
        if t is not None:
            per_allele[b.get("tag", b["scaffold_contig"])] = t
            tiers.append(t["tier"])
    rel_dupc["power_by_allele"] = per_allele
    if rel_dupc.get("called"):
        rel_dupc["interpretation"] = "CONFIRMED"
    elif not tiers or "poor" in tiers:
        rel_dupc["interpretation"] = "INDETERMINATE"        # some allele too shallow → cannot exclude
    elif "borderline" in tiers:
        rel_dupc["interpretation"] = "NEGATIVE_BORDERLINE"  # negative, but power at the edge
    else:
        rel_dupc["interpretation"] = "NEGATIVE"             # negative, adequately powered
    return rel_dupc


def dupc_verdict(genomic_bam=None, genome_ref=None, pon=None):
    """Specificity-controlled dupC verdict, DELEGATED to muc1_analyzer's genome-anchored caller.

    Why delegate (not call it off the faked-contig bins): a homopolymer dupC is one +1 C among ~240
    near-identical C-runs, so a position-blind scan pays a ~240× multiple-testing penalty that Bonferroni
    turns into p≈1 at AS depth. The ONLY validated way to beat that is `dupc_dispatch`'s motif-CONTEXT
    grouping (+ PoN null / paired), which needs reads anchored on the genome (or HP-phased) to build the
    contexts — the length-first faked-contig bins provide neither, and the caller then silently degrades to
    a span-requiring positional scan that abstains. So: if a chr1-aligned BAM is given, run the real caller;
    otherwise report NOT ASSESSED rather than a misleading silent negative. Length stays a faked-contig job;
    the dupC is a genome-anchored job."""
    if not genomic_bam:
        return {"assessed": False,
                "reason": ("dupC not assessed on the faked-contig bins (a specificity-controlled call needs "
                           "motif-context grouping over ~240 C-runs, hence genome/HP-anchored reads). Pass "
                           "--genomic-bam (chr1-aligned) [+ --dupc-pon] to run muc1_analyzer.dupc_dispatch.")}
    try:
        from muc1_analyzer.dupc_dispatch import dispatch_dupc
        from muc1_analyzer.detectors.vntr import _VNTR_ANCHORS as A
    except Exception as e:
        return {"assessed": False, "reason": f"dupc_dispatch unavailable: {e}"}
    pon_d = json.load(open(pon)) if (pon and isinstance(pon, str)) else pon
    out = dispatch_dupc(genomic_bam, chrom=A["chrom"], start=A["l_anchor"], end=A["r_anchor"],
                        ref=genome_ref, pon=pon_d, fast=(pon is None))
    out["assessed"] = True
    return out


def detect_rs4072037(bins, rel_dupc, ref, sample_level=False):
    """Read rs4072037 NATIVELY off the reads (offset 4265 on the faked ref): count C vs T, majority = the
    base. By default this is done on the MUTANT allele's reads only (the cis base the severity axis needs);
    with sample_level=True it pools all bins (sample genotype, used on a negative only if opted in).
    Returns {base, C, T, depth, frac_T, label, contig} or None (no carrier / no anchor / low depth)."""
    try:
        from muc1_analyzer.detectors.splice_snp import genotype_rs4072037_vntr_ref
    except Exception:
        return None
    bam = contig = None
    if sample_level:
        # any bin's bam covers the 5' flank anchor; use the first available
        for b in bins:
            if b.get("bam") and os.path.exists(b["bam"]):
                bam, contig = b["bam"], b.get("scaffold_contig"); break
    else:
        if not (rel_dupc and rel_dupc.get("called")):
            return None
        ctag, cctg = rel_dupc.get("carrier_tag"), rel_dupc.get("carrier_contig")
        for b in bins:
            if (ctag and b.get("tag") == ctag) or (not ctag and b.get("scaffold_contig") == cctg):
                bam, contig = b.get("bam"), b.get("scaffold_contig"); break
    if not bam or not os.path.exists(bam):
        return None
    try:
        g = genotype_rs4072037_vntr_ref(bam, ref, min_mq=0)
    except Exception:
        return None
    if not g:
        return None
    C, T = g.get("ref_count", 0), g.get("alt_count", 0)
    depth = C + T
    if depth < 5 or g.get("gt_label") == "low_depth":
        return {"base": None, "C": C, "T": T, "depth": depth, "frac_T": g.get("alt_frac"),
                "label": "low_depth", "contig": contig}
    return {"base": ("T" if T > C else "C"), "C": C, "T": T, "depth": depth,
            "frac_T": g.get("alt_frac"), "label": g.get("gt_label"), "contig": contig}


def resolve_rs4072037(provided, detected):
    """Reconcile the CLI override with the native detection. The provided value WINS (clinician's call),
    but a disagreement raises an alert for the report. Returns (final_base, source, alert|None)."""
    det = (detected or {}).get("base")
    if provided is not None:
        if det and det != provided:
            alert = (f"rs4072037 provided ({provided}) ≠ detected on the mutant allele "
                     f"({det} ; C:{detected['C']} / T:{detected['T']}, depth {detected['depth']}) — "
                     f"provided value kept")
            return provided, "provided (override)", alert
        return provided, "provided", None
    if det:
        return det, f"detected (C:{detected['C']} / T:{detected['T']})", None
    return None, "undetermined", None


def assemble_score(haps, sample, outdir, rs4072037=None):
    """Merge the per-allele top haplotypes → the project's own two-axis MUC1_Score. Returns
    (merged_json_path, fragments, score) using detectors.vntr + clinical_call (lazy import)."""
    from muc1_analyzer.detectors import vntr as VD
    from muc1_analyzer.clinical_call import score_from_fields
    merged = os.path.join(outdir, f"{sample}.merged_haplotypes.json")
    json.dump({"sample": sample, "haplotypes": haps}, open(merged, "w"), indent=2)
    frag = VD.from_analyzer_json(merged)
    pos = frag["mutation_variant"]["repeat_index"] if frag.get("mutation_variant") else None
    score = score_from_fields(frag.get("vntr_len_mut"), frag.get("vntr_len_healthy"), pos,
                              rs4072037_mut=rs4072037)
    return merged, frag, score


# ── main ─────────────────────────────────────────────────────────────────────────
def main(argv=None):
    ap = argparse.ArgumentParser(description="Alignment-first VNTR caller (length → per-allele call + score).")
    ap.add_argument("-b", "--bam", required=True, help="reads: fastq(.gz), uBAM, or any BAM/CRAM (seq only)")
    ap.add_argument("-r", "--ref", required=True, help="multi-contig faked VNTR reference (MUC1_VNTR_Nrepeats)")
    ap.add_argument("-s", "--sample", default="sample")
    ap.add_argument("-o", "--outdir", default="cigar_out")
    ap.add_argument("--ruler-contig", default=_RULER_CONTIG)
    ap.add_argument("--ref-cram", default=None, help="reference FASTA if --bam is a CRAM")
    ap.add_argument("--pacbio", action="store_true", help="PacBio HiFi (mappy map-hifi; default map-ont)")
    ap.add_argument("--maxmm", type=int, default=9, help="max mismatches for AL/AH anchor match")
    ap.add_argument("--offset", type=int, default=0, help="copies added to each SPAN read (calibration)")
    ap.add_argument("--sep", type=int, default=8, help="min copies between the two alleles")
    ap.add_argument("--merge-tol", type=int, default=3, help="agglomerate SPAN copies within this many units")
    ap.add_argument("--min-reads-per-allele", type=int, default=2, help="min SPAN reads to accept an allele")
    ap.add_argument("--contig-offset", type=int, default=0, help="shift between measured copies and scaffold contig")
    ap.add_argument("--rs4072037", choices=["C", "T"], default=None,
                    help="OVERRIDE the rs4072037 severity base of the mutant allele (by default it is read "
                         "natively off the mutant allele's reads). If given, it wins over the detection but a "
                         "disagreement is flagged on the report")
    ap.add_argument("--no-rs4072037-detect", action="store_true",
                    help="disable native rs4072037 detection; rely only on --rs4072037")
    ap.add_argument("--rs4072037-on-negative", action="store_true",
                    help="also read rs4072037 (sample-level) when NO frameshift is found (informational; "
                         "severity stays null without a mutant allele)")
    ap.add_argument("--call-min-depth", type=int, default=3,
                    help="min depth for a consensus base in `call` (lower to 2 for a shallow LONG allele, "
                         "so a germline frameshift seen by all its reads is not floored back to reference)")
    ap.add_argument("--call-ins-threshold", type=float, default=0.5,
                    help="min read fraction to keep a ≤3 bp insertion in the CONSENSUS nomenclature (kept at the "
                         "safe 0.5). NB: a homopolymer dupC is NOT reliably callable from this per-position "
                         "fraction (its value tracks ONT homopolymer error, not VAF) — the trustworthy verdict "
                         "comes from the specificity-calibrated statistical caller, reported separately below)")
    ap.add_argument("--no-scaffold-refine", action="store_true",
                    help="disable the indel-burden scaffold sweep; scaffold straight on the stage-1 allele "
                         "(the sweep corrects the ruler's few-unit under-count, e.g. 66→70, which is what "
                         "lets a 3′ frameshift reach the consensus)")
    ap.add_argument("--scaffold-window", type=int, default=6,
                    help="± units scanned around the stage-1 allele for the indel-burden scaffold sweep")
    ap.add_argument("--genomic-bam", default=None,
                    help="chr1-aligned BAM/CRAM for the specificity-controlled dupC call (delegated to "
                         "muc1_analyzer.dupc_dispatch at the MUC1 locus). Without it the dupC is reported NOT "
                         "ASSESSED — it cannot be called from the faked-contig bins at AS depth")
    ap.add_argument("--genome-ref-dupc", default=None, help="genome FASTA if --genomic-bam is a CRAM")
    ap.add_argument("--dupc-pon", default=None, help="panel-of-normals JSON for the dupC caller (sharpens the null)")
    ap.add_argument("--no-recover-reads", dest="recover_reads", action="store_false",
                    help="disable read recovery; use only the flank-anchored spanning reads per allele "
                         "(the default recovers ALL reads by best-fit allele — no spanning gate — which "
                         "lifts per-allele depth at the variant locus, e.g. 8→12 on a shallow long allele)")
    ap.set_defaults(recover_reads=True)
    ap.add_argument("--min-vntr-copies", type=int, default=1,
                    help="a read must cover at least this many VNTR copies to be assigned to an allele "
                         "(default 1 = exclude the 0-copy / off-target reads that only touch the locus "
                         "flanks; raise to ~5–10 to also drop very short tandem fragments)")
    ap.add_argument("--clair3-vcf", default=None,
                    help="Clair3 VCF (called vs the faked VNTR reference) → authoritative dupC/frameshift "
                         "verdict per allele (deep-learning ONT caller: sensitive AND specific, unlike the "
                         "consensus token or the statistical caller). Mirrors the reference pipeline.")
    ap.add_argument("--vntr-start", type=int, default=4574,
                    help="tandem start offset on the faked contigs for Clair3 coordinates (repeat = "
                         "(POS-start)//60); default 4574 matches MUC1_faked…Kirby.fa")
    ap.add_argument("--repeat-offset", type=int, default=4,
                    help="added to the reported dupC repeat index so it matches your clinical convention "
                         "(the analyzer counts repeats from the FIRST motif ~pos 4302, i.e. +4 vs the awk's "
                         "4574 origin; default 4 → e.g. raw 15 → repeat 19). Set 0 for the raw (POS-4574)//60")
    ap.add_argument("-t", "--threads", type=int, default=1,
                    help="threads for the ruler alignment (the bottleneck on large AS fastqs; mappy releases "
                         "the GIL so this scales well — e.g. -t 20). Default 1")
    ap.add_argument("--hom-ratio-near", type=float, default=0.5,
                    help="a NEIGHBOURING second length-allele (Δ < --sep) is folded as scatter when its SPAN "
                         "support is below this fraction of the dominant allele's (default 0.5)")
    ap.add_argument("--hom-ratio-far", type=float, default=0.25,
                    help="a DISTANT second length-allele (Δ ≥ --sep) is folded as scatter below this fraction "
                         "of the dominant's support (default 0.25; the absolute floor is --min-reads-per-allele)")
    ap.add_argument("--min-snp-diff", type=int, default=3,
                    help="a phasing split is kept as a true heterozygote only if the two reconstructed "
                         "haplotype consensuses differ by MORE than this many SNPs; at or below it (≤3 by "
                         "default) the sample is treated as homozygous and collapsed to a single haplotype")
    ap.add_argument("--phase-close-within", type=int, default=3,
                    help="when TWO alleles are called with a length gap ≤ this many repeats, separate them by "
                         "SNP phasing instead of length-binning (measurement noise can't split close lengths); "
                         "each haplotype is then re-measured for its own length. Lower --min-phasing-snps if the "
                         "two close alleles share few substitutional SNPs (default 3)")
    ap.add_argument("--no-phase-hom", action="store_true",
                    help="disable SNP phasing of a length-homozygous sample (two same-length alleles). By "
                         "default, when one length-allele is called, the reads are split by phasing SNPs into "
                         "two same-length haplotypes so the dupC can be attributed to one of them")
    ap.add_argument("--min-phasing-snps", type=int, default=2,
                    help="min DISCRIMINATING pileup SNPs (≥5 reads/group, ≥80%% within-group agreement) just "
                         "to ATTEMPT a phasing split; the real heterozygote/homozygote decision is then made on "
                         "the reconstructed consensus via --min-snp-diff (default 2)")
    ap.add_argument("--snp-min-af", type=float, default=0.30,
                    help="min minor-allele frequency for a candidate phasing SNP (default 0.30; max-af = 1 - this)")
    ap.add_argument("--keep-intermediates", action="store_true",
                    help="keep the per-allele intermediate files (BAMs, per-allele JSON/consensus/novel and "
                         "the per-allele report PDFs). By default only the final both-allele outputs are kept "
                         "(the single-allele PDFs are misleading — a carrier allele's own report omits the dupC)")
    ap.add_argument("--stage1-only", action="store_true", help="stop after the length call (no per-allele call/score)")
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    preset = "map-hifi" if a.pacbio else "map-ont"

    # Resolve inputs to ABSOLUTE paths up front. The script and its subprocess (`muc1_analyzer call`, run
    # with cwd=REPO_ROOT) must resolve a relative -r / -b the same way; an absolute path removes any
    # working-directory ambiguity and gives a clear early error instead of a confusing mid-run one.
    a.ref = os.path.abspath(os.path.expanduser(a.ref))
    if not os.path.isfile(a.ref):
        sys.exit(f"[cigar-len] reference not found: {a.ref}\n"
                 f"  → check the -r path is correct and readable from where you launch the script.")
    if os.path.exists(a.bam):
        a.bam = os.path.abspath(os.path.expanduser(a.bam))
    if not os.path.exists(a.ref + ".fai"):                    # index once, next to the reference
        try:
            pysam.faidx(a.ref)
        except Exception as e:
            sys.exit(f"[cigar-len] cannot build the .fai index for {a.ref} ({e}).\n"
                     f"  → the reference's directory is likely read-only. Pre-index it once with "
                     f"`samtools faidx {a.ref}`, or copy the reference to a writable location.")

    ruler_fa = build_ruler(a.ref, a.ruler_contig, a.outdir, a.maxmm)
    gAL, gAH, diag = locate_tandem_edges(ruler_fa, a.ruler_contig, a.maxmm)
    print(f"[cigar-len] ruler {a.ruler_contig}: tandem [{gAL},{gAH})  span {diag['anchor_span_bp']} bp "
          f"→ {diag['ruler_copies']} copies (ceiling)", file=sys.stderr)

    print(f"[cigar-len] streaming reads from {a.bam}; aligning to the ruler ({preset}, -t {a.threads}) …",
          file=sys.stderr)
    aln_bam = os.path.join(a.outdir, f"{a.sample}.ruler.bam")
    keep_thr = max(1, a.min_vntr_copies)                      # keep sequences only for VNTR-covering reads
    records, reads, hist, n_total, n_aligned = measure_reads(
        iter_reads(a.bam, a.ref_cram), ruler_fa, preset, gAL, gAH, a.offset, a.maxmm,
        out_bam=aln_bam, threads=a.threads, keep_thr=keep_thr)
    if not n_total:
        sys.exit("[cigar-len] no reads in the input.")
    print(f"[cigar-len] {n_total} reads processed; {n_aligned} aligned to the ruler; "
          f"{len(reads)} cover the VNTR (kept)", file=sys.stderr)
    if not records:
        sys.exit("[cigar-len] no read covers the VNTR (no flank/tandem signal).")
    n_span = print_histogram(hist)
    span_vals = [r["copies"] for r in records.values() if r["class"] == "SPAN"]
    print(f"\n[cigar-len] {n_span} SPAN reads / {n_aligned} aligned", file=sys.stderr)

    alleles, support = call_alleles_lowcount(span_vals, a.sep, a.merge_tol, a.min_reads_per_allele)
    alleles, support, folded = collapse_scatter_allele(alleles, support, a.sep, a.hom_ratio_near, a.hom_ratio_far)
    if folded:
        print(f"[cigar-len] second allele {folded['folded_allele']} ({folded['folded_support']} SPAN reads) "
              f"folded into {folded['kept_allele']} as scatter (ratio {folded['ratio']} < {folded['threshold']}, "
              f"Δ={folded['delta']} {'neighbouring' if folded['neighbouring'] else 'distant'}) → length-homozygous",
              file=sys.stderr)
        summary["scatter_fold"] = folded
    kind = "HET" if len(alleles) == 2 else ("HOM" if len(alleles) == 1 else "NONE")
    print(f"[cigar-len] alleles ({kind}): {alleles} copies  (SPAN support {support})", file=sys.stderr)

    summary = {"sample": a.sample, "ruler_contig": a.ruler_contig, "tandem": diag,
               "n_reads": len(reads), "n_aligned": len(records), "n_span": n_span,
               "alleles": alleles, "allele_support": support, "kind": kind,
               "span_histogram": dict(sorted(Counter(span_vals).items())),
               "aligned_bam": aln_bam, "bins": []}

    if not a.stage1_only and alleles:
        contigs = available_contigs(a.ref)
        if not contigs:
            print("[cigar-len] no MUC1_VNTR_Nrepeats contigs → skipping per-allele call.", file=sys.stderr)
        else:
            # bin: SPAN by nearest allele; LEFT/RIGHT (lower bound) → nearest allele long enough to hold it
            bins = defaultdict(list)
            for name, rec in records.items():
                c, cls = rec["copies"], rec["class"]
                if cls == "SPAN":
                    tgt = min(alleles, key=lambda al: abs(al - c))
                elif cls in ("LEFT", "RIGHT"):
                    cand = [al for al in alleles if al >= c - a.merge_tol] or alleles
                    tgt = min(cand, key=lambda al: abs(al - c))
                else:
                    continue
                bins[tgt].append(name)

            # PASS 1 — determine each allele's scaffold length/contig from the SPANNING bins (accurate length)
            allele_contigs, bin_meta = {}, {}
            for L in sorted(bins):
                names = bins[L]
                meas = [records[n]["copies"] for n in names if records[n]["class"] == "SPAN"]
                bin_median = int(round(statistics.median(meas))) if meas else None
                if a.no_scaffold_refine:
                    scaffold_L, cname = pick_contig(L, contigs, a.contig_offset)
                    burden = None
                else:
                    scaffold_L, cname, burden_tbl = refine_scaffold(
                        reads, names, contigs, L + a.contig_offset, a.ref, preset, a.outdir, a.scaffold_window)
                    burden = burden_tbl[:5]
                    if burden_tbl:
                        print(f"[cigar-len] allele {L}: scaffold sweep → {scaffold_L} "
                              f"(burden {burden_tbl[0][0]}; stage-1 said {L})", file=sys.stderr)
                allele_contigs[scaffold_L] = cname
                bin_meta[scaffold_L] = {"stage1_allele": L, "bin_median": bin_median, "burden": burden}

            # RECOVER READS — reassign reads to alleles, RELIABLY: a read whose allele is determinable
            # (spanning, or a one-flank lower bound already past the short allele) keeps its KNOWN allele;
            # only genuinely ambiguous internal partials fall back to best-fit. Mirrors the reference
            # pipeline's all-reads architecture but without ever mis-placing an allele-informative read.
            # Length still comes from the spanning reads (pass 1).
            if a.recover_reads and len(allele_contigs) >= 1:
                stage1_to_scaffold = {m["stage1_allele"]: L for L, m in bin_meta.items()}
                short_allele = min(alleles)
                # Only VNTR-covering reads are eligible. The 0-copy / off-target flood was already dropped
                # while streaming (kept_thr), so `reads` here is exactly the VNTR-covering set.
                vntr_reads = {name: reads[name] for name, rec in records.items()
                              if rec["copies"] >= a.min_vntr_copies and name in reads}
                print(f"[cigar-len] recovery pool: {len(vntr_reads)} reads cover the VNTR "
                      f"(≥{a.min_vntr_copies} copies; off-target reads already dropped upstream)", file=sys.stderr)
                known_allele_of = {}
                for name, rec in records.items():
                    if name not in vntr_reads:
                        continue
                    c, cls = rec["copies"], rec["class"]
                    if cls == "SPAN":
                        known_allele_of[name] = stage1_to_scaffold.get(min(alleles, key=lambda al: abs(al - c)))
                    elif cls in ("LEFT", "RIGHT") and (c - a.merge_tol) > short_allele:
                        known_allele_of[name] = stage1_to_scaffold.get(max(alleles))
                    # else: ambiguous → resolved by burden inside the assigner
                recovered, rel_counts = assign_reads_to_alleles(
                    vntr_reads, allele_contigs, known_allele_of, a.ref, preset, a.outdir)
                summary["read_assignment"] = {allele_contigs[L]: rel_counts[L] for L in allele_contigs}
                print("[cigar-len] read recovery (reliable: known allele by length, burden only for "
                      "ambiguous internals):", file=sys.stderr)
                for L in sorted(allele_contigs):
                    rc_ = rel_counts[L]
                    print(f"[cigar-len]   {allele_contigs[L]}: {len(recovered[L])} reads "
                          f"({rc_['confident']} allele-sure + {rc_['ambiguous']} ambiguous)", file=sys.stderr)
            else:
                recovered = {L: bins[bin_meta[L]["stage1_allele"]] for L in allele_contigs}

            # Build the HAPLOTYPE UNITS to call. Normally one per length-allele. If the sample is
            # length-homozygous (one allele), try SNP phasing on the scaffold to split the pool into two
            # SAME-LENGTH haplotypes — otherwise the two alleles would be merged into one consensus and the
            # dupC could not be attributed to a haplotype. Each unit has a UNIQUE tag (two phased haplotypes
            # share the contig, so the contig alone is not an identity).
            haplo_units = []
            phase_pool = None                                  # set when a phasing split is provisional (may collapse)
            if len(allele_contigs) == 1 and not a.no_phase_hom:
                L = next(iter(allele_contigs)); cname = allele_contigs[L]; names = recovered[L]
                phased = phase_hom_reads(cname, names, reads, a.ref, preset, a.outdir,
                                         a.min_phasing_snps, a.min_reads_per_allele, a.call_min_depth, a.snp_min_af)
                if phased:
                    A, B, ndisc = phased
                    print(f"[cigar-len] one length-cluster ({L} rep) phased by sequence: {ndisc} discriminating "
                          f"SNP(s) → 2 haplotypes (hap A {len(A)} reads, hap B {len(B)} reads); re-measuring each",
                          file=sys.stderr)
                    summary["phasing"] = {"length_homozygous_cluster": L, "n_discriminating_snps": ndisc,
                                          "reads_hapA": len(A), "reads_hapB": len(B)}
                    for lbl, grp in (("A", A), ("B", B)):        # re-scaffold each phased haplotype on its OWN
                        meas = [records[n]["copies"] for n in grp if n in records and records[n]["class"] == "SPAN"]
                        stage1 = int(round(statistics.median(meas))) if meas else L
                        if a.no_scaffold_refine:
                            sL, cn = pick_contig(stage1, contigs, a.contig_offset); bd = None
                        else:
                            sL, cn, bt = refine_scaffold(reads, grp, contigs, stage1 + a.contig_offset,
                                                         a.ref, preset, a.outdir, a.scaffold_window)
                            bd = bt[:5] if bt else None
                        haplo_units.append({"tag": f"{a.sample}_{sL}rep_hap{lbl}", "cname": cn, "scaffold_L": sL,
                                            "names": grp, "meta": {"stage1_allele": stage1, "bin_median": stage1,
                                                                   "burden": bd}, "hap": lbl})
                    phase_pool = {"names": names, "cname": cname, "scaffold_L": L, "meta": bin_meta[L]}
                else:
                    print(f"[cigar-len] length-homozygous ({L} rep): no informative phasing SNP → single "
                          f"haplotype (true homozygote or unphaseable)", file=sys.stderr)
                    summary["phasing"] = {"length_homozygous": True, "n_phasing_snps": 0}
            elif (len(allele_contigs) == 2 and not a.no_phase_hom
                  and abs(max(allele_contigs) - min(allele_contigs)) <= a.phase_close_within):
                # Two REAL but very close-length alleles (Δ ≤ phase_close_within): length-binning can't
                # separate them (measurement noise ±1-2 mixes the reads), so phase them by SEQUENCE instead.
                # Pool both bins on the LARGER scaffold (short-allele reads then carry a clean ~60 bp deletion),
                # phase by SNP, then RE-MEASURE each haplotype's own length so 43 | 44 can still show.
                Ls = sorted(allele_contigs)
                big = allele_contigs[Ls[-1]]
                pool = recovered[Ls[0]] + recovered[Ls[1]]
                phased = phase_hom_reads(big, pool, reads, a.ref, preset, a.outdir,
                                         a.min_phasing_snps, a.min_reads_per_allele, a.call_min_depth, a.snp_min_af)
                if phased:
                    A, B, ndisc = phased
                    print(f"[cigar-len] close-length alleles {Ls} (Δ={Ls[-1]-Ls[0]}): {ndisc} discriminating "
                          f"SNP(s) → phased by sequence (hap A {len(A)} reads, hap B {len(B)} reads)", file=sys.stderr)
                    summary["phasing"] = {"close_length_alleles": Ls, "n_discriminating_snps": ndisc,
                                          "reads_hapA": len(A), "reads_hapB": len(B)}
                    for lbl, grp in (("A", A), ("B", B)):        # re-scaffold each haplotype on its OWN length
                        meas = [records[n]["copies"] for n in grp if n in records and records[n]["class"] == "SPAN"]
                        stage1 = int(round(statistics.median(meas))) if meas else Ls[-1]
                        if a.no_scaffold_refine:
                            sL, cn = pick_contig(stage1, contigs, a.contig_offset); bd = None
                        else:
                            sL, cn, bt = refine_scaffold(reads, grp, contigs, stage1 + a.contig_offset,
                                                         a.ref, preset, a.outdir, a.scaffold_window)
                            bd = bt[:5] if bt else None
                        haplo_units.append({"tag": f"{a.sample}_{sL}rep_hap{lbl}", "cname": cn, "scaffold_L": sL,
                                            "names": grp, "meta": {"stage1_allele": stage1, "bin_median": stage1,
                                                                   "burden": bd}, "hap": lbl})
                    dom = max(allele_contigs, key=lambda L: len(recovered[L]))
                    phase_pool = {"names": pool, "cname": allele_contigs[dom], "scaffold_L": dom,
                                  "meta": bin_meta[dom]}
                else:
                    print(f"[cigar-len] close-length alleles {Ls} (Δ={Ls[-1]-Ls[0]}): no discriminating SNP → "
                          f"cannot phase; falling back to length bins (they may be poorly separated)", file=sys.stderr)
                    summary["phasing"] = {"close_length_alleles": Ls, "n_discriminating_snps": 0, "phased": False}
            if not haplo_units:
                for L in sorted(allele_contigs):
                    haplo_units.append({"tag": f"{a.sample}_{L}rep", "cname": allele_contigs[L], "scaffold_L": L,
                                        "names": recovered[L], "meta": bin_meta[L], "hap": None})

            # PASS 2 — per haplotype: realign the read set, call consensus/motifs/frameshift
            per_allele_haps = []
            per_allele_fastas, per_allele_novels = [], []
            intermediates = [ruler_fa, ruler_fa + ".fai", aln_bam, aln_bam + ".bai"]
            results = []
            for u in haplo_units:
                res = call_unit(u, reads, a.ref, a.outdir, preset, a.call_min_depth, a.call_ins_threshold,
                                a.keep_intermediates)
                results.append(res)
                intermediates += res["ints"]
                lbl = f" [hap {res['hap_label']}]" if res["hap_label"] else ""
                print(f"[cigar-len] haplotype {res['scaffold_L']}rep{lbl}: {res['n_names']} reads → {res['n_ok']} "
                      f"on {res['cname']} ({res['n_no']} unmapped) → call rc={res['rc']}", file=sys.stderr)

            # BIOLOGICAL GATE — decide heterozygote vs homozygote for a phasing split. The two haplotypes
            # are genuinely distinct if they differ by LENGTH (different re-measured scaffold) OR by more
            # than --min-snp-diff consensus SNPs. Only when the length does NOT separate them (same
            # scaffold) AND they differ by ≤ --min-snp-diff SNPs are they effectively identical (a true
            # homozygote, or an ONT-noise split) → collapse to a single pooled haplotype. So a real
            # 43 | 44 carrier (distinct lengths, few substitutional SNPs) is KEPT, while a same-length
            # near-identical split is folded. This is the clinician's rule, applied only where length can't.
            haps_ok = [r for r in results if r["hap"] is not None]
            if phase_pool is not None and len(haps_ok) == 2:
                same_length = haps_ok[0]["scaffold_L"] == haps_ok[1]["scaffold_L"]
                nsnp = _consensus_snp_diffs(haps_ok[0]["hap"], haps_ok[1]["hap"])
                summary.setdefault("phasing", {})["consensus_snp_diffs"] = nsnp
                summary["phasing"]["same_length"] = bool(same_length)
                if same_length and nsnp <= a.min_snp_diff:
                    print(f"[cigar-len] phased haplotypes: same length and only {nsnp} SNP(s) "
                          f"(≤ {a.min_snp_diff}) → homozygous; collapsing to a single haplotype", file=sys.stderr)
                    summary["phasing"]["collapsed_to_homozygous"] = True
                    single = {"tag": f"{a.sample}_{phase_pool['scaffold_L']}rep", "cname": phase_pool["cname"],
                              "scaffold_L": phase_pool["scaffold_L"], "names": phase_pool["names"],
                              "meta": phase_pool["meta"], "hap": None}
                    res = call_unit(single, reads, a.ref, a.outdir, preset, a.call_min_depth,
                                    a.call_ins_threshold, a.keep_intermediates)
                    intermediates += res["ints"]
                    print(f"[cigar-len] single haplotype {res['scaffold_L']}rep: {res['n_names']} reads → "
                          f"{res['n_ok']} on {res['cname']} → call rc={res['rc']}", file=sys.stderr)
                    results = [res]
                else:
                    why = ("distinct lengths" if not same_length else f"{nsnp} SNP(s) > {a.min_snp_diff}")
                    print(f"[cigar-len] phased haplotypes kept ({why}) → confirmed heterozygote", file=sys.stderr)

            for res in results:                                # commit the final haplotype set
                if res["hap"] is not None:
                    per_allele_haps.append(res["hap"])
                per_allele_fastas.append(res["cons_fa"])
                per_allele_novels.append(res["novel_tsv"])
                summary["bins"].append(res["bin"])

            # combined both-allele BAM (final output, kept) — one BAM holding both haplotypes for IGV
            vntr_bam = os.path.join(a.outdir, f"{a.sample}.vntr.bam")
            if merge_allele_bams(summary["bins"], vntr_bam):
                summary["vntr_bam"] = vntr_bam

            # combined consensus FASTA (both haplotypes) + merged novel-motif TSV (candidates not yet in
            # KNOWN_REPEATS), aggregated across the two alleles — reusing the caller's own consensus/motif logic.
            cons_out = os.path.join(a.outdir, f"{a.sample}.consensus.fa")
            if merge_consensus_fasta(per_allele_fastas, cons_out):
                summary["consensus_fasta"] = cons_out
            novel_out = os.path.join(a.outdir, f"{a.sample}.novel_motifs.tsv")
            novel_path, n_novel = merge_novel_motifs(per_allele_novels, novel_out, a.sample)
            summary["novel_motifs_tsv"] = novel_path
            summary["n_novel_motifs"] = n_novel
            print(f"[cigar-len] consensus FASTA → {os.path.basename(cons_out)} ; "
                  f"novel motifs → {os.path.basename(novel_out)} ({n_novel} candidat(s))", file=sys.stderr)

            # trustworthy dupC verdict. PRIORITY: --clair3-vcf (deep-learning ONT caller — sensitive AND
            # specific, the reference-pipeline authority; correctly calls a low-depth dupC our statistical
            # test would miss). Else --genomic-bam → dupc_dispatch (PoN null). Else the self-contained
            # context+paired caller. The consensus X-59dupC token in the nomenclature stays labelled fragile.
            clair3_frag = None
            contig_L = {c: L for L, c in allele_contigs.items()}
            if a.clair3_vcf:
                c3 = parse_clair3_vcf(a.clair3_vcf, allele_contigs, a.vntr_start)
                summary["clair3_frameshifts"] = {c: v for c, v in c3.items()}
                carrier = None
                for cname, v in c3.items():
                    if v is not None and (carrier is None or v["is_dupc"]):
                        carrier = (cname, v)
                lens = sorted(allele_contigs)
                if carrier:
                    cname, v = carrier
                    mut_L = contig_L[cname]
                    healthy = [L for L in lens if L != mut_L] or [mut_L]
                    clair3_frag = {"mutation_present": True, "vntr_len_mut": mut_L,
                                   "vntr_len_healthy": healthy[0],
                                   "mutation_variant": {"motif": v["label"], "repeat_index": v["repeat"],
                                                        "known_pathogenic": v["is_dupc"]}}
                    rel_dupc = {"assessed": True, "called": True, "interpretation": "CONFIRMED",
                                "method": "clair3", "carrier_contig": cname, "variant": v}
                    print(f"[cigar-len] dupC (Clair3): CONFIRMED — {v['label']} on {cname}, repeat {v['repeat']}",
                          file=sys.stderr)
                else:
                    clair3_frag = {"mutation_present": False,
                                   "vntr_len_mut": max(lens), "vntr_len_healthy": min(lens)}
                    rel_dupc = {"assessed": True, "called": False, "interpretation": "NEGATIVE",
                                "method": "clair3"}
                    print("[cigar-len] dupC (Clair3): no PASS frameshift in the tandem → NEGATIVE", file=sys.stderr)
            elif a.genomic_bam:
                rel_dupc = dupc_interpretation(dupc_verdict(a.genomic_bam, a.genome_ref_dupc, a.dupc_pon),
                                               summary["bins"])
            else:
                # DEFAULT: run-length-shift caller (general homopolymer-expansion detector, self-calibrated).
                rel_dupc = runlen_shift_dupc(summary["bins"], a.ref, a.vntr_start, a.repeat_offset)
                if rel_dupc.get("called") and rel_dupc.get("variant"):
                    cname = rel_dupc["carrier_contig"]
                    mut_L = contig_L[cname]
                    lens = sorted(allele_contigs)
                    healthy = [L for L in lens if L != mut_L] or [mut_L]
                    v = rel_dupc["variant"]
                    clair3_frag = {"mutation_present": True, "vntr_len_mut": mut_L,
                                   "vntr_len_healthy": healthy[0],
                                   "mutation_variant": {"motif": v["label"], "repeat_index": v["repeat"],
                                                        "known_pathogenic": True}}
                    print(f"[cigar-len] dupC (run-length shift): CONFIRMED — repeat {v['repeat']} on "
                          f"{cname} ({v['k_ge_mut']}/{v['n']} reads ≥{8}C, p_bonf {v['p_bonf']:.1e})",
                          file=sys.stderr)
                else:                                           # negative → annotate with detection power
                    rel_dupc = dupc_interpretation(rel_dupc, summary["bins"])
            summary["reliable_dupc"] = rel_dupc
            if not rel_dupc.get("assessed"):
                print(f"[cigar-len] dupC: NOT ASSESSED — {rel_dupc.get('reason')}", file=sys.stderr)
            elif rel_dupc.get("power_by_allele"):
                interp = rel_dupc.get("interpretation")
                powers = "; ".join(f"{c.split('_')[-1]}:n{t['n_cover']}"
                                   f"{'✓' if t['adequately_powered'] else '✗'}"
                                   for c, t in rel_dupc.get("power_by_allele", {}).items())
                print(f"[cigar-len] dupC verdict: {interp} (called={rel_dupc.get('called')}; "
                      f"power {powers})", file=sys.stderr)
                if interp == "INDETERMINATE":
                    need = max((t["reads_for_power0.90"]["d=0.3"]
                                for t in rel_dupc.get("power_by_allele", {}).values()), default=12)
                    print(f"[cigar-len]   → underpowered: ~{need} mutant-allele reads needed for power 0.90 "
                          f"(d=0.3); current coverage cannot exclude a dupC.", file=sys.stderr)

            # rs4072037 severity base — read NATIVELY off the mutant allele's reads (majority C/T at
            # offset 4265), not just taken from the CLI. Only on a carrier (severity is cis to the mutant
            # allele); on a negative, read it only if explicitly asked. The CLI value, if given, overrides
            # the detection but a disagreement is flagged on the report.
            rs_detected = None
            if not a.no_rs4072037_detect:
                if rel_dupc and rel_dupc.get("called"):
                    rs_detected = detect_rs4072037(summary["bins"], rel_dupc, a.ref)
                elif a.rs4072037_on_negative:
                    rs_detected = detect_rs4072037(summary["bins"], rel_dupc, a.ref, sample_level=True)
            rs_final, rs_source, rs_alert = resolve_rs4072037(a.rs4072037, rs_detected)
            if rs_detected:
                summary["rs4072037_detection"] = rs_detected
            summary["rs4072037_used"] = {"base": rs_final, "source": rs_source, "alert": rs_alert}
            if rs_detected and rs_detected.get("base"):
                print(f"[cigar-len] rs4072037 (mutant allele): {rs_detected['base']} "
                      f"(C:{rs_detected['C']} / T:{rs_detected['T']}, depth {rs_detected['depth']}) → "
                      f"{'severe' if rs_final == 'C' else 'protective' if rs_final == 'T' else 'n/a'} [{rs_source}]",
                      file=sys.stderr)
            if rs_alert:
                print(f"[cigar-len]   ⚠ {rs_alert}", file=sys.stderr)

            if per_allele_haps or clair3_frag:
                try:
                    from muc1_analyzer.clinical_call import score_from_fields
                    if per_allele_haps:
                        merged, frag_c, score_c = assemble_score(per_allele_haps, a.sample, a.outdir, rs_final)
                        summary["merged_haplotypes_json"] = merged
                    else:
                        frag_c, score_c = {}, {}
                    frag = clair3_frag if clair3_frag is not None else frag_c
                    if clair3_frag is not None:                 # score from the Clair3-authoritative variant
                        pos = frag["mutation_variant"]["repeat_index"] if frag.get("mutation_variant") else None
                        score = score_from_fields(frag.get("vntr_len_mut"), frag.get("vntr_len_healthy"),
                                                  pos, rs4072037_mut=rs_final)
                    else:
                        score = score_c
                    summary["vntr_fragments"] = frag
                    oi = score.get("onset_index")
                    if oi is not None:
                        mean_frac, std_frac, n_cal = _cohort_ref()
                        score["cohort_onset"] = {"z": round((oi - (1.0 - mean_frac)) / std_frac, 3),
                                                 "sample_onset_index": oi,
                                                 "cohort_mean_onset_index": round(1.0 - mean_frac, 4),
                                                 "cohort_std": std_frac, "n": n_cal}
                    summary["muc1_two_axis_score"] = score
                    summary["clinical_call"] = clinical_string(frag, per_allele_haps)
                    mpdf = os.path.join(a.outdir, f"{a.sample}.MUC1_score.pdf")
                    if write_merged_pdf(mpdf, a.sample, frag, score, per_allele_haps, rel_dupc,
                                        rs_detected=rs_detected, rs_alert=rs_alert):
                        summary["merged_pdf"] = mpdf
                    print(f"[cigar-len] MUC1_Score: {clinical_string(frag, per_allele_haps)}  carrier={frag.get('mutation_present')} "
                          f"onset={score.get('onset_score')} severity={score.get('severity_score')}", file=sys.stderr)
                except Exception as e:
                    summary["score_error"] = str(e)
                    print(f"[cigar-len] score assembly failed: {e}", file=sys.stderr)

    # Keep ONLY the final, both-allele outputs; drop the per-allele intermediates (their single-allele PDFs
    # are misleading — e.g. a carrier allele's own report never shows the dupC). --keep-intermediates keeps them.
    if not a.keep_intermediates:
        removed = 0
        for p in dict.fromkeys(x for x in locals().get("intermediates", []) if x):
            try:
                os.remove(p)
                removed += 1
            except OSError:
                pass
        summary.pop("aligned_bam", None)                     # strip references to now-deleted files
        for b in summary.get("bins", []):
            b.pop("bam", None)
            b.pop("analyzer_json", None)
        if removed:
            print(f"[cigar-len] cleaned {removed} per-allele intermediate file(s); kept the final outputs "
                  f"({a.sample}.MUC1_score.pdf, .cigar.summary.json, .vntr.bam, .consensus.fa, "
                  f".novel_motifs.tsv, .merged_haplotypes.json)", file=sys.stderr)

    sjson = os.path.join(a.outdir, f"{a.sample}.cigar.summary.json")
    json.dump(summary, open(sjson, "w"), indent=2, default=str)
    print(f"[cigar-len] wrote {sjson}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
