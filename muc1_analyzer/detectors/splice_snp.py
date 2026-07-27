"""Detector for the MUC1 exon 2 splice SNP: rs4072037 (chr1:155,192,276, A/G).

Functional splice-regulatory variant (MUC1 isoform) + renal/urate GWAS signal.
Genotyped by targeted pileup, phased by HP tag. Its proximity (~3 kb from the DEL, ~170 bp
from the 5'UTR CpG window) requires disentangling its LD/phase with the DEL -> see
muc1_analyzer/del_snp_matrix.py.

Target build: GRCh38 "chr". Handles BAM + CRAM (reference_filename).
"""
from __future__ import annotations
from typing import Optional

import pysam

from ..config import GRCh38, SNP_RS4072037, MIN_MQ, MIN_BQ

_SNP = GRCh38["SPLICE_SNP_EXON2"]


def _open(bam: str, ref: Optional[str]):
    mode = "rc" if bam.endswith(".cram") else "rb"
    kw = {"reference_filename": ref} if (mode == "rc" and ref) else {}
    return pysam.AlignmentFile(bam, mode, **kw)


def _query_index_at(read, target: int):
    """Read base index aligned to reference position `target` (0-based).

    CIGAR walk (like the DEL detector) — reliable on BAM and CRAM, unlike
    pileup(). Returns None if `target` falls inside a deletion/skip of the read.
    """
    qpos = 0
    rpos = read.reference_start
    for op, ln in (read.cigartuples or []):
        if op in (0, 7, 8):            # M/=/X : consume read + ref
            if rpos <= target < rpos + ln:
                return qpos + (target - rpos)
            qpos += ln
            rpos += ln
        elif op == 1:                  # I : consumes the read
            qpos += ln
        elif op in (2, 3):             # D/N : consume the ref
            if rpos <= target < rpos + ln:
                return None            # position deleted in this read
            rpos += ln
        elif op == 4:                  # S : consumes the read
            qpos += ln
        # H (5) / P (6) : nothing
    return None


def genotype_snp(bam: str, ref: Optional[str] = None, *, chrom: str = None,
                 pos1: int = None, ref_allele: str = None, alt_allele: str = None,
                 min_mq: int = MIN_MQ, min_bq: int = MIN_BQ,
                 het_lo: float = 0.20, het_hi: float = 0.80, min_depth: int = 5) -> dict:
    """Genotype rs4072037 by pileup; count ref/alt globally and per HP.

    Returns {snp_genotype: 0|1|2|None, gt_label, alt_frac, depth,
              ref_count, alt_count, hp_alt: '1'|'2'|None,
              per_hp: {...}}.
    `hp_alt` = the haplotype carrying the ALT allele majoritarily (for the DEL x SNP phase).
    """
    chrom = chrom or _SNP.chrom
    pos1 = pos1 if pos1 is not None else _SNP.start
    pos0 = pos1 - 1

    # ── Count of the 4 bases (forward reference orientation), global + per HP ──
    # NB : MUC1 is on the − strand, so the forward alleles = complement of A/G (= T/C).
    # We stay strand-agnostic : REF comes from the FASTA, ALT = the other segregating allele.
    base_cnt = {b: 0 for b in "ACGT"}
    hp_base = {"1": {b: 0 for b in "ACGT"}, "2": {b: 0 for b in "ACGT"}}

    with _open(bam, ref) as af:
        for rd in af.fetch(chrom, pos0, pos0 + 1):
            if rd.is_unmapped or rd.is_secondary or rd.is_supplementary \
                    or rd.is_duplicate or rd.mapping_quality < min_mq:
                continue
            qi = _query_index_at(rd, pos0)
            if qi is None or rd.query_sequence is None:
                continue
            if rd.query_qualities is not None and rd.query_qualities[qi] < min_bq:
                continue
            base = rd.query_sequence[qi].upper()
            if base not in base_cnt:
                continue
            base_cnt[base] += 1
            if rd.has_tag("HP"):
                hp = str(rd.get_tag("HP"))
                if hp in ("1", "2"):
                    hp_base[hp][base] += 1

    # ── REF : explicit override, else base from the FASTA, else majority allele ──
    if ref_allele:
        ref_allele = ref_allele.upper()
    elif ref:
        try:
            fa = pysam.FastaFile(ref)
            ref_allele = fa.fetch(chrom, pos0, pos0 + 1).upper()
            fa.close()
        except Exception:
            ref_allele = None
    if not ref_allele or ref_allele not in base_cnt:
        ref_allele = max(base_cnt, key=base_cnt.get)          # fallback : most frequent

    # ── ALT : override, else most frequent segregating allele != REF ──
    if alt_allele:
        alt_allele = alt_allele.upper()
    else:
        others = sorted((b for b in "ACGT" if b != ref_allele),
                        key=lambda b: base_cnt[b], reverse=True)
        alt_allele = others[0] if others else "N"

    cnt = {ref_allele: base_cnt.get(ref_allele, 0), alt_allele: base_cnt.get(alt_allele, 0),
           "other": sum(v for b, v in base_cnt.items() if b not in (ref_allele, alt_allele))}
    hp_cnt = {hp: {ref_allele: hp_base[hp].get(ref_allele, 0),
                   alt_allele: hp_base[hp].get(alt_allele, 0)} for hp in ("1", "2")}

    depth = cnt[ref_allele] + cnt[alt_allele]
    alt_frac = (cnt[alt_allele] / depth) if depth else None
    if depth < min_depth:
        gt_label, gt = "NA(lowcov)", None
    elif alt_frac <= (1 - het_hi):
        gt_label, gt = f"{ref_allele}/{ref_allele}", 0
    elif alt_frac >= het_hi:
        gt_label, gt = f"{alt_allele}/{alt_allele}", 2
    elif het_lo <= alt_frac <= het_hi:
        gt_label, gt = f"{ref_allele}/{alt_allele}", 1
    else:
        gt_label, gt = "ambiguous", None

    # HP carrying the ALT (alt fraction per HP)
    hp_alt = None
    fa = {}
    for hp in ("1", "2"):
        d = hp_cnt[hp][ref_allele] + hp_cnt[hp][alt_allele]
        fa[hp] = (hp_cnt[hp][alt_allele] / d) if d else None
    if fa["1"] is not None and fa["2"] is not None:
        if fa["1"] > 0.5 and fa["2"] < 0.5:
            hp_alt = "1"
        elif fa["2"] > 0.5 and fa["1"] < 0.5:
            hp_alt = "2"

    return {"snp_genotype": gt, "gt_label": gt_label,
            "ref_allele": ref_allele, "alt_allele": alt_allele,
            "alt_frac": None if alt_frac is None else round(alt_frac, 3),
            "depth": depth, "ref_count": cnt[ref_allele], "alt_count": cnt[alt_allele],
            "other_count": cnt["other"], "base_counts": base_cnt, "hp_alt": hp_alt,
            "per_hp": {"1": {**hp_cnt["1"], "alt_frac": fa["1"]},
                       "2": {**hp_cnt["2"], "alt_frac": fa["2"]}},
            "rsid": SNP_RS4072037["rsid"]}


# ══════════════════════════════════════════════════════════════════════════════
# rs4072037 — T2T-NATIVE genotyping, directly on the VNTR-reference alignment
# ══════════════════════════════════════════════════════════════════════════════
# The multi-contig VNTR reference differs between contigs ONLY in the number of repeat
# units: the 150 contigs share a 4636 bp identical 5' prefix (and a 3462 bp identical
# suffix), the array varying in between. rs4072037 sits in that constant 5' flank, so it
# is at the SAME offset on every contig — genotypeable on the already-realigned BAM.
# That removes the GRCh38 dependency: `prepare → call → score` never needs a second
# alignment space just to read this SNP.
#
# Verified on MUC1_fakedVNTR1to150revcomplKirby.fa: the 41 bp anchor below occurs in
# 150/150 contigs at a single offset, placing the SNP at offset 4265 (0-based), i.e.
# inside the constant prefix.
#
# STRAND: the VNTR reference is the reverse complement of the GRCh38 forward strand (it is
# in GENE orientation). On it the alleles are G/A; on GRCh38-forward they are C/T. Everything
# is reported in the GRCh38-forward alphabet (C/T) so it stays interchangeable with
# `genotype_snp` and with the `--rs4072037-mut C|T` severity input.

#: 41 bp anchor centred on rs4072037, in VNTR-reference (gene-strand) orientation.
RS4072037_VNTR_ANCHOR = "AAACCCGCAACAGTTGTTACGGGTTCTGGTCATGCAAGCTC"
RS4072037_ANCHOR_SNP_INDEX = 20          # the SNP is the middle base of the anchor
RS4072037_VNTR_OFFSET = 4265             # 0-based, verified constant across the 150 contigs
#: gene-strand (VNTR reference) -> GRCh38-forward alphabet
_GENE_TO_FWD = {"G": "C", "A": "T", "C": "G", "T": "A"}


def find_rs4072037_offset(vntr_ref_fasta: str) -> int | None:
    """Locate rs4072037 in a VNTR reference by ANCHOR SEARCH, and check the offset is constant.

    Anchoring rather than trusting the hardcoded constant means a VNTR reference rebuilt with
    a different flank length still works, and a reference where the offset is NOT constant is
    REFUSED instead of silently genotyping the wrong base. Returns the 0-based offset, or None.
    """
    offsets, name, seq = set(), None, []

    def _flush():
        if name is None:
            return
        s = "".join(seq)
        i = s.find(RS4072037_VNTR_ANCHOR)
        offsets.add(None if i < 0 else i + RS4072037_ANCHOR_SNP_INDEX)

    with open(vntr_ref_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                _flush()
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
        _flush()

    if not offsets or None in offsets or len(offsets) != 1:
        return None                       # absent from some contig, or not at a single offset
    return offsets.pop()


def genotype_rs4072037_vntr_ref(bam: str, vntr_ref: str, *, min_mq: int = 0,
                                min_bq: int = MIN_BQ, min_depth: int = 5,
                                het_lo: float = 0.20, het_hi: float = 0.80) -> dict | None:
    """Genotype rs4072037 on a BAM aligned to the multi-contig VNTR reference (no GRCh38).

    Returns the same shape as `genotype_snp` (alleles in the GRCh38-forward alphabet, and
    `snp_genotype` = T dosage 0/1/2 as the severity axis expects), or None if the reference
    has no usable constant-offset anchor.

    `min_mq` defaults to 0 on purpose: 150 near-identical contigs make almost every read
    MAPQ 0, so the usual MAPQ filter would discard the whole locus (cf. CLAUDE.md).
    """
    offset = find_rs4072037_offset(vntr_ref)
    if offset is None:
        return None

    counts = {"C": 0, "T": 0, "other": 0}
    hp_counts = {"1": {"C": 0, "T": 0}, "2": {"C": 0, "T": 0}}
    with pysam.AlignmentFile(bam, "rb") as af:
        for read in af.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mq:
                continue
            qi = _query_index_at(read, offset)
            if qi is None:
                continue
            quals = read.query_qualities
            if quals is not None and quals[qi] < min_bq:
                continue
            fwd = _GENE_TO_FWD.get(read.query_sequence[qi].upper(), "other")
            key = fwd if fwd in ("C", "T") else "other"
            counts[key] += 1
            if key != "other":
                hp = read.get_tag("HP") if read.has_tag("HP") else None
                if str(hp) in hp_counts:
                    hp_counts[str(hp)][key] += 1

    depth = counts["C"] + counts["T"]
    if depth < min_depth:
        return {"snp_genotype": None, "gt_label": "low_depth", "depth": depth,
                "ref_allele": "C", "alt_allele": "T", "alt_frac": None,
                "ref_count": counts["C"], "alt_count": counts["T"], "other_count": counts["other"],
                "offset": offset, "source": "vntr_ref_native", "rsid": SNP_RS4072037["rsid"]}

    t_frac = counts["T"] / depth
    if t_frac >= het_hi:
        gt, label = 2, "T/T"
    elif t_frac <= (1 - het_hi):
        gt, label = 0, "C/C"
    elif het_lo <= t_frac <= het_hi:
        gt, label = 1, "C/T"
    else:
        gt, label = None, "ambiguous"

    return {"snp_genotype": gt, "dosage": gt, "genotype": label, "gt_label": label,
            "ref_allele": "C", "alt_allele": "T", "alt_frac": round(t_frac, 3),
            "depth": depth, "ref_count": counts["C"], "alt_count": counts["T"],
            "other_count": counts["other"], "per_hp": hp_counts, "offset": offset,
            "source": "vntr_ref_native", "rsid": SNP_RS4072037["rsid"]}
