#!/usr/bin/env python3
"""
VNTR Haplotype Caller
=====================
Logic:
  - The FASTA reference is dedicated to the VNTR: each contig represents
    one possible VNTR length.
  - The BAM contains the region's reads aligned against these contigs.
  - The two contigs receiving the most reads correspond to the patient's
    two haplotypes.
  - For each, a consensus is reconstructed, then the sequence is translated
    into motif nomenclature (greedy exact then approximate matching).

Dependencies: pysam  (pip install pysam)

Usage:
    python vntr_haplotype_caller.py -b reads.bam -r vntr_ref.fa [options]

Options:
    -b / --bam           BAM file (must be .bai indexed)
    -r / --ref           Custom FASTA reference (one contig = one VNTR length)
    -o / --output        Text output file (default: stdout)
    -s / --sample        Sample name
    --top-n              Number of haplotypes to report (default: 2)
    --min-depth          Min depth to call a consensus base (default: 3)
    --min-bq             Min base quality (default: 20)
    --min-mq             Min mapping quality (default: 20)
    --max-mismatch       Max mismatches for approximate matching (default: 3)
    --ins-threshold      Min fraction of reads carrying an insertion (default: 0.5)
    --long-ins-threshold Min fraction to include a long insertion > 3 bp (default: 0.4)
    --del-threshold      Min fraction of reads carrying a deletion (default: 0.6)
    --dump-consensus     Write the consensus sequences as FASTA to stderr
    --fasta-consensus    Path of the output FASTA file for the consensus sequences
    --pdf                Path of the output PDF file (visual haplotype report)
    --json               JSON file for the detailed export
    --homo-threshold     Min fraction of reads on the dominant contig to trigger
                         intra-contig phasing (default: 0.65)
    --snp-min-af         Min minor-allele frequency for a position to be
                         considered an informative phasing SNP (default: 0.15)
    -v / --verbose       Detailed matching output
"""

import pysam
import argparse
import bisect
import sys
import json
import re
from collections import Counter, defaultdict

# ══════════════════════════════════════════════════════════════════════════════
# Known motif dictionary
# Key = sequence (60 bp for standard motifs, variable length for indels)
# Value = motif name
# ══════════════════════════════════════════════════════════════════════════════
KNOWN_REPEATS = {
    'GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCGCCCGCA': 'A',
    'GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA': 'B',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCAA': 'C',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCCGCCCCGGGCTCCACCGCCCCCCCA': 'D',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCCGCCCCGGGCTCCACCGCGCCCGCA': 'E',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCACA': 'F',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCGCCCGCA': 'G',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCGCCCCA': 'H',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCGCCCCCA': 'I',
    'GCCCACGGTGTCACCTCGGCACCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCGCCCGCA': 'J',
    'GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCTGGGCTCCACCGCCCCCCCA': 'K',
    'GCCCACGGTGTCACCTCAGCCCCGGACACCAGGCCGGCCCCGGGCTCCGCCGCCGCCCCA': 'L',
    'GCCCACGGTGTCACCTCGGCACCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA': 'M',
    'GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCACGGGCTCCACCGCCCCCCCA': 'N',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCAGGCTCCACCGCCCCCCCA': 'O',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCAGCCCCGGGCTCCACCGCCCCCCCA': 'P',
    'GCCCACAGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA': 'Q',
    'GCCCACGGTGTCACCTCGGCCACGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA': 'R',
    'GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCCCCCGCA': 'S',
    'GCCCACGGTGTCACCTCAGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCGCCCCA': 'T',
    'GCCCACGGTGTCACCTCGGCCCCGGAGACCAGGCCGGCCCCGGGCTCCACCGCGCCCGCA': 'U',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCACCCCCA': 'V',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCG': 'W',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA': 'X',
    'GCCCACGGTGTCACCTCGGCCCCGGAGACCAGGCCGGCCCCGGGCTCCACCGCGCCCCCA': 'Y',
    'GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCGCCCCCA': 'Z',
    'GCCCACGGTGTCACCTCGGCCCCGGAGACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA': 'aA',
    'CCCACGGTGTCACCATCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCCCGGCCCCG': 'aB',
    'GCCCACGATGTCACCTCAGCCCCGGACAACAAGCCAGCCCCGGGCTCCATCGCCCCCCCA': 'aC',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCCGCCCCAGGCTCCACCGCCCCCCCA': 'aD',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCGGGGCTCCACCGCCCCCCCA': 'aE',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCTGGGCTCCACCGCCCCCCCA': 'aF',
    'GCCCACGGTGTCACCTCGGCCCCGGACAGCAGGCCGGCCCCGGGCTCCACCGCGCCCGCA': 'aG',
    'GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCCTCCCCA': 'aH',
    'GCCCACGGTGTCACCTCGGCCCTGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA': 'aI',
    'GCCCACGGTGTCACCTTGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCAA': 'aJ',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCCGCCCCGGGCTCCACCGCGCCCCCA': 'aK',
    'GCCCACGGTGTCACCTCGGCCCCGGAGACCAGGCCGGCCCCGGGCTCCACCGCCCCCGCA': 'aL',
    'GCCCACGGTGTTACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCGCCCGCA': 'aM',
    'GCCCATGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA': 'aN',
    'GGCTCCACCGCACCCCCAGCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCG': 'aO',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGACCCGCCCCGGGCTCCACCGCGCCCGCA': 'aP',
    'GCCCACGGTGTCACCTCGGCCCCGGACAGCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA': 'aQ',
    'GCCCACGGTGTCACCTCGGCCCCGGATACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA': 'aR',
    'GGCTCCACCGCCCCCCAAGCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCG': 'aS',
    'GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACTGCGCCCGCA': 'aT',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCGCA': 'aU',
    'GCCCACGGTGTCACCTCAGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA': 'aV',
    'GCCCACGGTGTCACCTCGGCCCCGGACACAAGGCCGGCCCCGGGCTCCACCGCCCCCCCA': 'aW',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCAGGCTCCACCGCCCCCCAA': 'aX',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCAGGCTCCACCGCGCCCGCA': 'aY',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCACCCCCCCA': 'aZ',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCACGCCCGCA': 'bA',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCCCCCCCCCA': 'bB',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCGCCA': 'bC',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCGGGCCCCGGGCTCCACCCGGGCCCCG': 'bD',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCTGGCCCCGGGCTCCACCGCCCCCCCA': 'bE',
    'GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCCCCCCAA': 'bF',
    'GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCTGGGCTCCACCGCGCCCGCA': 'bG',
    'GCCCACGGTGTCACCTCGTCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA': 'bH',
    'GCCCACGGTGTCACCTTGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA': 'bI',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCCGCCCCGGGCTCCACCGCCCCCGCA': 'bJ',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGTCCCCCCA': 'bK',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCTGGCCCCGGGCTCCACCGCGCCCGCA': 'bL',
    # ── Motifs with indels (variable lengths) ──────────────────────────
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCGGGCCCCGGGCTCCACCCCGGCCCCGGGCTCCACCGCCCCCCCA': 'X-33_34insCGGGCCCCGGGCTCCACC',
    'AAGGAGACTTCGGCTACCCAGAGAAGTTCAGTGCCCAGCTCTACTGAGAAGAATGCTGTG': '1',
    'AGTATGACCAGCAGCGTACTCTCCAGCCACAGCCCCGGTTCAGGCTCCTCCACCACTCAG': '2',
    'GGACAGGATGTCACTCTGGCCCCGGCCACGGAACCAGCTTCAGGTTCAGCTGCCACCTGG': '3',
    'GGACAGGATGTCACCTCGGTCCCAGTCACCAGGCCAGCCCTGGGCTCCACCACCCCGCCA': '4',
    'GGACAGGATGTCACCTCGGTCCCAGTCACCAGGCCAGCCCTGGGCTCCACCACCCCACCA': '4+',
    'GGACAGGATGTCACCTCCGTCCCAGTCACCAGGCCAGCCCTGGGCTCCACCACCCCGCCA': '4++',
    'GGACAGGATGTCACCTCGGTCCCAGTCACCAGGCCAGCACTGGGCTCCACCACCCCGCCA': '4+++',
    'GCCCACGATGTCACCTCAGCCCCGGACAACAAGCCAGCCCCGGGCTCCACCGCCCCCCCA': '5',
    'GCCCACGATGTCACCTCAGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCAA': '5C',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCGGGCCCCGGGCTCCACCCCGGCCCCG': '6',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCCCGGCCCCG': '6+',
    'GGCTCCACCGCCCCCCCAGCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCG': '7',
    'GGCTCCACCGCCTCCCCAGCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCG': '7+',
    'GGCTCCACCGCCCCCCCAGCCCATGGTGTCACCTCGGCCCCGGACAACAGGCCCGCCTTG': '8',
    'GGCTCCACCGCCCCTCCAGTCCACAATGTCACCTCGGCCTCAGGCTCTGCATCAGGCTCA': '9',
    'GCCCACGGTGTCACCTCGGCCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCGCCCGCA': 'A-23dupC',
    'GCCCACGGTGTCACCTCGGCCCCCGGACACCAGGCCCGCCCCGGGCTCCACCGCGCCCGCA': 'E-23dupC',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCCA': 'X-59dupC',
    'GCCCACGATGTCACCTCAGCCCCGGACAACAAGCCAGCCCCGGGCTCCACCGCCCCCCCCA': '5-59dupC',
    'GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCCCCCCCCA': 'B-59dupC',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCCGCCCCGGGCTCCACCGCCCCCCCCA': 'D-59dupC',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCAA': 'X-60dupA',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCAAA': 'C-60dupA',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCGCCGGCCCCGGGCTCCACCGCCCCCCCA': 'X-33_48dupGCCGGCCCCGGGCTCC',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCCCCCA': 'X-56_59dupCCCC',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCGGGCTCCACCGCCCCCCAA': 'C-42_57dupGGGCTCCACCGCCCCC',
    'GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCCATCCCCA': 'aH-54_55insA',
    'GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCCCCCCGCA': 'B-58_59insG',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCCGCA': 'X-58_59insG',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCCGCCCCGGGCTCCACCGCCGCCCGCA': 'bJ-54_55insG',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCCCCCA': 'X-58_59delCC',
    'GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCCCCCA': 'B-58_59delCC',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCCGCCCCGGGCTCCACCGCCCCCA': 'D-58_59delCC',
    'GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCGCCCA': 'A-58_59delCC',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCCGCCCCGGGCTCCACCGCGCCCA': 'E-58_59delCC',
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGCGCCCA': 'G-58_59delCC',
    'GCCCACGGTGTCACCTCGGCCGGCCCCGGGCTCCACCGCCCCCCCA': 'X-18_31delGGCCCCGGACACCA',
    # ── Additions 2026-07-12: dupG(52) + del8_27 (CANONICAL base-X forms — ⚠ base TO BE VERIFIED vs reads
    #    WOO_Jer / INF_Lio). NB del8_27 is ABSORBED by the tandem (INF_Lio reads decoded into clean
    #    60 bp units, cf. MUC1_log 2026-07-12) → nomenclature entry, will NOT match its reads in practice.
    'GCCCACGGTGTCACCTCGGCCCCGGACACCAGGCCGGCCCCGGGCTCCACCGGCCCCCCCA': 'X-52dupG',
    'GCCCACGACCAGGCCGGCCCCGGGCTCCACCGCCCCCCCA': 'X-8_27delGTGTCACCTCGGCCCCGGAC',
}

# Pre-sort by decreasing length (priority to motifs with long indels)
_SORTED_MOTIFS = sorted(KNOWN_REPEATS.items(), key=lambda x: len(x[0]), reverse=True)
_EXACT_INDEX   = {seq: name for seq, name in KNOWN_REPEATS.items()}
_MOTIF_LENGTHS = sorted({len(k) for k in KNOWN_REPEATS}, reverse=True)

# Block prefilter (pigeonhole) for the Hamming step over 60 bp motifs.
# A motif at Hamming distance <= 3 necessarily shares >= 1 of the 4 EXACT 15 bp blocks (3
# mismatches can touch at most 3 blocks). We index the 60 bp motifs by
# (block_index, 15-mer); at scan time, we only Hamming-test candidates sharing a block
# — a few instead of 79. Exact for max_mismatch <= 3 (otherwise fallback to full scan).
# _MOTIF60 preserves the order of _SORTED_MOTIFS → decoding tie-break preserved.
_MOTIF60 = [(s, n) for s, n in _SORTED_MOTIFS if len(s) == 60]
_HAM_BLK = 15
_HAM_NBLK = 4
_MOTIF60_BLOCK_IDX = {}
for _i60, (_s60, _n60) in enumerate(_MOTIF60):
    for _b60 in range(_HAM_NBLK):
        _MOTIF60_BLOCK_IDX.setdefault(
            (_b60, _s60[_b60 * _HAM_BLK:(_b60 + 1) * _HAM_BLK]), []).append(_i60)

# For the Levenshtein step: deletions (<60bp) tested FIRST,
# then insertions (>60bp). At equal distance, the deletion thus wins
# over the corresponding insertion (which shares the same 60bp prefix).
_SORTED_INDEL_DEL = [(s, n) for s, n in KNOWN_REPEATS.items() if len(s) < 60]
_SORTED_INDEL_INS = sorted(
    [(s, n) for s, n in KNOWN_REPEATS.items() if len(s) > 60],
    key=lambda x: len(x[0]), reverse=True,
)
_SORTED_INDELS = _SORTED_INDEL_DEL + _SORTED_INDEL_INS

# Motif names involved in the terminal correction
_MOTIF_33_34INS = 'X-33_34insCGGGCCCCGGGCTCCACC'
_MOTIF_33_34INS_LEN = 78
_MOTIFS_6  = {'6', '6+'}
_MOTIFS_7  = {'7', '7+'}


def fix_terminal_6_7(
    matched: list,
    unmatched: list,
    consensus: str,
    max_mismatch: int,
) -> tuple:
    """
    Post-correction: if the LAST X-33_34insCGGGCCCCGGGCTCCACC in the
    nomenclature is a terminal false positive (= should be 6/6+ + 7/7+),
    replace it.

    Mechanism:
      Locate the last X-33_34ins in `matched`. Re-match the consensus
      from its position, forcing first 6/6+, then 7/7+, then letting
      the greedy continue freely. If this hypothesis yields fewer unmatched
      bases than the current hypothesis over the same region, apply it.

    Returns (matched_corrected, unmatched_corrected).
    """
    # Find the last X-33_34ins
    last_ins_idx = None
    for i in range(len(matched) - 1, -1, -1):
        if matched[i]['name'] == _MOTIF_33_34INS:
            last_ins_idx = i
            break

    if last_ins_idx is None:
        return matched, unmatched

    pos_ins = matched[last_ins_idx]['pos']
    n = len(consensus)

    # ── Find the best 6/6+ at pos_ins ───────────────────────────────
    best6, best6_dist = None, float('inf')
    for motif_seq, motif_name in KNOWN_REPEATS.items():
        if motif_name not in _MOTIFS_6:
            continue
        mlen = len(motif_seq)
        if pos_ins + mlen > n:
            continue
        window = consensus[pos_ins: pos_ins + mlen]
        dist = 0 if window == motif_seq else hamming(window, motif_seq)
        if dist < best6_dist:
            best6_dist = dist
            best6 = {'name': motif_name, 'pos': pos_ins, 'length': mlen,
                     'mismatches': dist, 'sequence': window}

    if best6 is None or best6_dist > max_mismatch:
        return matched, unmatched

    # ── Find the best 7/7+ right after ─────────────────────────────
    pos7 = pos_ins + best6['length']
    best7, best7_dist = None, float('inf')
    for motif_seq, motif_name in KNOWN_REPEATS.items():
        if motif_name not in _MOTIFS_7:
            continue
        mlen = len(motif_seq)
        if pos7 + mlen > n:
            continue
        window = consensus[pos7: pos7 + mlen]
        dist = 0 if window == motif_seq else hamming(window, motif_seq)
        if dist < best7_dist:
            best7_dist = dist
            best7 = {'name': motif_name, 'pos': pos7, 'length': mlen,
                     'mismatches': dist, 'sequence': window}

    if best7 is None or best7_dist > max_mismatch:
        return matched, unmatched

    # ── Continue the greedy from pos_after_6_7 ──────────────────────────
    pos_after = pos7 + best7['length']
    tail_matched, tail_unmatched = _greedy_from(consensus, pos_after, max_mismatch)

    # ── Compare unmatched bases: current vs 6+7 hypothesis ─────────
    # Current unmatched in the region [pos_ins, end]
    current_unmatched_len = sum(
        len(b['sequence']) for b in unmatched if b['pos'] >= pos_ins
    )
    # Unmatched in the 6+7 hypothesis
    new_unmatched_len = sum(len(b['sequence']) for b in tail_unmatched)

    if new_unmatched_len > current_unmatched_len:
        return matched, unmatched  # the 6+7 hypothesis is worse → keep X-33_34ins

    # ── Apply the correction ────────────────────────────────────────────
    corrected_matched   = matched[:last_ins_idx] + [best6, best7] + tail_matched
    corrected_unmatched = [b for b in unmatched if b['pos'] < pos_ins] + tail_unmatched
    return corrected_matched, corrected_unmatched


def _greedy_from(consensus: str, start: int, max_mismatch: int) -> tuple:
    """
    Run the greedy matching from `start` to the end of the consensus.
    Lightweight version (without the post-correction, to avoid recursion).
    """
    pos      = start
    n        = len(consensus)
    matched  = []
    leftover = []

    while pos < n:
        best      = None
        best_dist = float('inf')

        for length in _MOTIF_LENGTHS:
            if pos + length > n:
                continue
            window = consensus[pos: pos + length]
            if window in _EXACT_INDEX:
                best = {'name': _EXACT_INDEX[window], 'pos': pos, 'length': length,
                        'mismatches': 0, 'sequence': window}
                best_dist = 0
                break

        if best_dist > 0 and pos + 60 <= n:
            window = consensus[pos: pos + 60]
            if max_mismatch <= _HAM_NBLK - 1:
                cand = set()
                for b in range(_HAM_NBLK):
                    hit = _MOTIF60_BLOCK_IDX.get((b, window[b * _HAM_BLK:(b + 1) * _HAM_BLK]))
                    if hit:
                        cand.update(hit)
                order = sorted(cand)
            else:
                order = range(len(_MOTIF60))
            for i in order:
                motif_seq, motif_name = _MOTIF60[i]
                dist = hamming(window, motif_seq, cap=best_dist)
                if dist < best_dist:
                    best_dist = dist
                    best = {'name': motif_name, 'pos': pos, 'length': 60,
                            'mismatches': dist, 'sequence': window}
                if best_dist == 0:
                    break

        if best_dist > 0 and pos + min(_MOTIF_LENGTHS) <= n:
            for motif_seq, motif_name in _SORTED_MOTIFS:
                mlen = len(motif_seq)
                if mlen == 60:
                    continue
                if pos + mlen > n + max_mismatch:
                    continue
                end   = min(pos + mlen, n)
                chunk = consensus[pos: end]
                dist  = levenshtein(chunk, motif_seq, max_dist=min(best_dist, max_mismatch + 1) - 1)
                if dist < best_dist:
                    best_dist = dist
                    best = {'name': motif_name, 'pos': pos, 'length': mlen,
                            'mismatches': dist, 'sequence': chunk}

        if best and best_dist <= max_mismatch:
            matched.append(best)
            pos += best['length']
        else:
            leftover.append({'pos': pos, 'base': consensus[pos]})
            pos += 1

    unmatched_blocks = []
    if leftover:
        blk = {'pos': leftover[0]['pos'], 'sequence': leftover[0]['base']}
        for seg in leftover[1:]:
            if seg['pos'] == blk['pos'] + len(blk['sequence']):
                blk['sequence'] += seg['base']
            else:
                unmatched_blocks.append(blk)
                blk = {'pos': seg['pos'], 'sequence': seg['base']}
        unmatched_blocks.append(blk)

    return matched, unmatched_blocks


# ══════════════════════════════════════════════════════════════════════════════
# Utilitaires
# ══════════════════════════════════════════════════════════════════════════════

def hamming(s1: str, s2: str, cap=None) -> int:
    """Hamming distance (strings of equal length).

    If `cap` is provided, stops as soon as the distance REACHES `cap` and returns `cap`: the result
    is then only valid as a ">= cap" bound. Usage: `hamming(w, m, cap=best_dist)` — a motif is only
    kept if `dist < best_dist`, so once the distance reaches `best_dist` there is no point counting
    further. Without `cap`, behavior unchanged (exact distance). Gain on motifs that diverge early
    (the majority).
    """
    if cap is None:
        return sum(a != b for a, b in zip(s1, s2))
    d = 0
    for a, b in zip(s1, s2):
        if a != b:
            d += 1
            if d >= cap:
                return d
    return d


def levenshtein(s1: str, s2: str, max_dist: int = None) -> int:
    """Edit distance. With `max_dist`, **BANDED** DP (width 2·max_dist+1 around the
    diagonal): any out-of-band cell has distance > max_dist, so ignoring it does not change
    the result when the true distance is ≤ max_dist; otherwise we return max_dist+1. Cost
    O(m·max_dist) instead of O(m·n) — this is the hot spot of `match_motifs` on a divergent
    consensus (80k calls). Without `max_dist`, classic full DP.

    Interchangeable with the old version for `match_motifs` usage: the latter only reads
    the value via `dist < best_dist` and only accepts a match if `dist ≤ max_mismatch`, and the
    band is exact over all `dist ≤ max_dist` and "> threshold" otherwise → decisions unchanged.
    """
    m, n = len(s1), len(s2)
    if max_dist is None:                                   # full DP (exact distance, unbounded)
        prev = list(range(n + 1))
        for i in range(1, m + 1):
            curr = [i] + [0] * n
            for j in range(1, n + 1):
                curr[j] = prev[j - 1] if s1[i - 1] == s2[j - 1] \
                    else 1 + min(prev[j], curr[j - 1], prev[j - 1])
            prev = curr
        return prev[n]

    if abs(m - n) > max_dist:
        return max_dist + 1
    d = max_dist
    INF = d + 1                                            # out of band / > threshold
    prev = [j if j <= d else INF for j in range(n + 1)]    # row 0: j deletions, cap out of band
    for i in range(1, m + 1):
        curr = [INF] * (n + 1)
        if i <= d:
            curr[0] = i
        lo = i - d if i - d > 1 else 1
        hi = i + d if i + d < n else n
        row_min = curr[0]
        s1c = s1[i - 1]
        for j in range(lo, hi + 1):
            cost = 0 if s1c == s2[j - 1] else 1
            v = prev[j - 1] + cost
            dj = prev[j] + 1
            if dj < v:
                v = dj
            ij = curr[j - 1] + 1
            if ij < v:
                v = ij
            curr[j] = v
            if v < row_min:
                row_min = v
        if row_min > d:
            return d + 1
        prev = curr
    return prev[n] if prev[n] <= d else d + 1


# ══════════════════════════════════════════════════════════════════════════════
# Read counting per contig
# ══════════════════════════════════════════════════════════════════════════════

def count_reads_per_contig(bam_path: str, min_mq: int = 20) -> dict:
    """
    Returns a dict {contig_name: read_count} for reads
    correctly mapped (non-supplementary, non-duplicate, MQ ≥ min_mq).
    """
    counts = Counter()
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam.fetch(until_eof=True):
            if (read.is_unmapped
                    or read.is_secondary
                    or read.is_supplementary
                    or read.is_duplicate
                    or read.mapping_quality < min_mq):
                continue
            counts[read.reference_name] += 1
    return counts


# ══════════════════════════════════════════════════════════════════════════════
# Consensus construction for a given contig
# ══════════════════════════════════════════════════════════════════════════════

def _tally_bases_cigar(bam_path, contig, ref_len, min_bq, min_mq, read_names=None):
    """Per-ref-position base tallies for the consensus builders, via fetch() + CIGAR walk over ONLY the
    (optional) `read_names` subset. Replaces a whole-contig `pysam.pileup()` that processed EVERY read on
    the contig even when just a 200-read sample was wanted (the deep-HiFi `call` hang: pileup builds a
    column from all reads, then the caller discards the ones not in read_names). Also honours the project
    rule "genotype via fetch()+CIGAR, never pileup()".

    Returns (bases_at, ins_after): bases_at[pos] = ['A'|'C'|'G'|'T'|'-', ...], ins_after[pos] = [ins_seq,…].
    Matches the old `pileup(stepper='all', min_base_quality, min_mapping_quality, truncate=True)`: skip
    unmapped/secondary/qcfail/duplicate reads, require mapq ≥ min_mq, drop bases with bq < min_bq, record a
    deletion as '-', and attach an insertion to the last kept aligned base before it (iff that base passed).
    """
    bases_at = defaultdict(list)
    ins_after = defaultdict(list)
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for r in bam.fetch(contig):
            if (r.is_unmapped or r.is_secondary or r.is_qcfail or r.is_duplicate
                    or r.mapping_quality < min_mq):
                continue
            if read_names is not None and r.query_name not in read_names:
                continue
            seq = r.query_sequence
            cig = r.cigartuples
            if seq is None or cig is None:
                continue
            quals = r.query_qualities
            qpos = 0
            rpos = r.reference_start
            prev_rp = None                                   # ref pos of the last kept aligned base
            prev_kept = False
            for op, length in cig:
                if op == 0 or op == 7 or op == 8:            # M / = / X  (aligned)
                    for k in range(length):
                        rp = rpos + k
                        if 0 <= rp < ref_len and (quals is None or quals[qpos + k] >= min_bq):
                            bases_at[rp].append(seq[qpos + k].upper())
                    end = rpos + length - 1
                    prev_rp = end
                    prev_kept = (0 <= end < ref_len) and (quals is None or quals[qpos + length - 1] >= min_bq)
                    qpos += length
                    rpos += length
                elif op == 1:                                # I  → after the last kept aligned base
                    if prev_kept and prev_rp is not None and 0 <= prev_rp < ref_len:
                        ins_after[prev_rp].append(seq[qpos: qpos + length].upper())
                    qpos += length
                elif op == 2:                                # D  → '-' at each deleted ref position
                    for k in range(length):
                        rp = rpos + k
                        if 0 <= rp < ref_len:
                            bases_at[rp].append('-')
                    rpos += length
                    prev_rp, prev_kept = None, False
                elif op == 3:                                # N  (ref skip)
                    rpos += length
                    prev_rp, prev_kept = None, False
                elif op == 4:                                # S  (soft clip)
                    qpos += length
                # H (5) / P (6): consume nothing
    return bases_at, ins_after


def build_consensus_for_contig(
    bam_path: str,
    ref_fa: pysam.FastaFile,
    contig: str,
    min_depth: int = 3,
    min_bq: int = 20,
    min_mq: int = 20,
    ins_threshold: float = 0.5,
    long_ins_threshold: float = 0.25,
    del_threshold: float = 0.6,
    verbose: bool = False,
) -> str:
    """
    Builds the consensus sequence for a contig (= a candidate VNTR haplotype).
    Low-coverage positions are filled in from the reference sequence.

    Two insertion thresholds:
      ins_threshold      (default 0.50) – insertions ≤ 3 bp.
      long_ins_threshold (default 0.25) – insertions > 3 bp (e.g. X-33_34ins 18 bp).
        Long insertions appear only in reads covering the carrier motif
        → fraction often well below 0.40 in a long VNTR.
      For 1 bp insertions in a homopolymer context (e.g. X-60dupA), the
        signal is aggregated over a ±1 bp window before comparison to the threshold.
    """
    INS_SHORT_MAX = 3   # max length to consider an insertion "short"
    ref_seq = ref_fa.fetch(contig).upper()
    ref_len = len(ref_seq)

    # pos → [base,...] / pos → [inserted_seq,...], via fetch()+CIGAR (was a whole-contig pileup).
    bases_at, ins_after = _tally_bases_cigar(bam_path, contig, ref_len, min_bq, min_mq)

    consensus_parts = []
    low_cov = 0

    for pos in range(ref_len):
        bases = bases_at[pos]
        depth = len(bases)

        if depth < min_depth:
            consensus_parts.append(ref_seq[pos])
            low_cov += 1
        else:
            ctr = Counter(bases)
            top_base, top_count = ctr.most_common(1)[0]

            if top_base == '-' and top_count / depth >= del_threshold:
                pass  # consensus deletion → position removed
            else:
                non_del = [b for b in bases if b != '-']
                if non_del:
                    consensus_parts.append(Counter(non_del).most_common(1)[0][0])
                elif top_base != '-':
                    consensus_parts.append(top_base)

        # ── Consensus insertion after this position ───────────────────
        ins_list = ins_after.get(pos, [])
        if ins_list:
            ctr_ins = Counter(ins_list)
            best_ins, best_count = ctr_ins.most_common(1)[0]
            ins_len = len(best_ins)

            if ins_len <= INS_SHORT_MAX:
                # Short insertions: aggregate over ±1 bp (homopolymers)
                if ins_len == 1:
                    for delta in (-1, 1):
                        adj = ins_after.get(pos + delta, [])
                        best_count += sum(1 for s in adj if s == best_ins)
                    ref_depth = max(
                        len(bases_at.get(max(0, pos - 1), [])),
                        depth,
                        len(bases_at.get(min(ref_len - 1, pos + 1), [])),
                    )
                    frac = best_count / max(ref_depth, 1)
                else:
                    frac = best_count / max(depth, 1)
                threshold = ins_threshold
            else:
                # Long insertions (e.g. X-33_34ins): lower threshold
                frac = best_count / max(depth, 1)
                threshold = long_ins_threshold

            if frac >= threshold:
                consensus_parts.append(best_ins)
                if verbose and ins_len > INS_SHORT_MAX:
                    print(
                        f"  [ins] pos {pos:>5}  +{ins_len}pb"
                        f"  frac={frac:.2f}  '{best_ins[:25]}'",
                        file=sys.stderr,
                    )

    consensus = ''.join(consensus_parts)

    if verbose:
        pct = 100 * low_cov / ref_len if ref_len else 0
        print(
            f"  [consensus] {contig}: {len(consensus)} bp "
            f"({low_cov} low-coverage pos, {pct:.1f}%)",
            file=sys.stderr,
        )

    return consensus


# ══════════════════════════════════════════════════════════════════════════════
# Greedy matching consensus → motifs
# ══════════════════════════════════════════════════════════════════════════════

def match_motifs(
    consensus: str,
    max_mismatch: int = 3,
    verbose: bool = False,
) -> tuple:
    """
    Greedy traversal: at each position, look for the best known motif.

    Priorities:
      1. Exact match (longest length first → captures indels)
      2. Approximate match by Hamming on the 60 bp motifs
      3. Approximate match by Levenshtein on motifs of length ≠ 60 bp

    Returns (matched, unmatched_blocks).
    """
    pos      = 0
    n        = len(consensus)
    matched  = []
    leftover = []   # characters with no match

    while pos < n:
        best      = None
        best_dist = float('inf')

        # ── 1. Match exact ──────────────────────────────────────────────
        for length in _MOTIF_LENGTHS:
            if pos + length > n:
                continue
            window = consensus[pos: pos + length]
            if window in _EXACT_INDEX:
                best = {
                    'name': _EXACT_INDEX[window],
                    'pos': pos,
                    'length': length,
                    'mismatches': 0,
                    'sequence': window,
                }
                best_dist = 0
                break

        # ── 2. Approximate Hamming match (60 bp motifs) — block prefilter ──
        if best_dist > 0 and pos + 60 <= n:
            window = consensus[pos: pos + 60]
            if max_mismatch <= _HAM_NBLK - 1:
                # pigeonhole: candidates = motifs sharing >= 1 exact block (covers all <=3 mm)
                cand = set()
                for b in range(_HAM_NBLK):
                    hit = _MOTIF60_BLOCK_IDX.get((b, window[b * _HAM_BLK:(b + 1) * _HAM_BLK]))
                    if hit:
                        cand.update(hit)
                order = sorted(cand)                       # _SORTED_MOTIFS order → tie-break preserved
            else:
                order = range(len(_MOTIF60))               # exact fallback: all 60 bp motifs
            for i in order:
                motif_seq, motif_name = _MOTIF60[i]
                dist = hamming(window, motif_seq, cap=best_dist)
                if dist < best_dist:
                    best_dist = dist
                    best = {
                        'name': motif_name,
                        'pos': pos,
                        'length': 60,
                        'mismatches': dist,
                        'sequence': window,
                    }
                if best_dist == 0:
                    break

        # ── 3. Approximate Levenshtein match (motifs of length ≠ 60) ────────
        # Always attempted as soon as an approximate 60 bp match exists (best_dist > 0)
        # so that motifs with indels compete with the Hamming match.
        # Order: deletions (<60 bp) FIRST, then insertions (>60 bp, desc).
        # At equal distance, the deletion thus wins over the corresponding
        # insertion, which shares the same 60 bp prefix.
        if best_dist > 0 and pos + min(_MOTIF_LENGTHS) <= n:
            for motif_seq, motif_name in _SORTED_INDELS:
                mlen = len(motif_seq)
                if pos + mlen > n + max_mismatch:
                    continue
                end   = min(pos + mlen, n)
                chunk = consensus[pos: end]
                dist  = levenshtein(chunk, motif_seq, max_dist=min(best_dist, max_mismatch + 1) - 1)
                if dist < best_dist:
                    best_dist = dist
                    best = {
                        'name': motif_name,
                        'pos': pos,
                        'length': mlen,
                        'mismatches': dist,
                        'sequence': chunk,
                    }

        # (X-33_34ins vs 6+7 arbitration handled in post-correction, cf. fix_terminal_6_7)

        # ── Decision ────────────────────────────────────────────────────
        if best and best_dist <= max_mismatch:
            matched.append(best)
            if verbose:
                tag = 'exact' if best_dist == 0 else f"{best_dist} mm"
                print(
                    f"    pos {pos:>5}  len {best['length']:>3}  "
                    f"[{tag:>5}]  {best['name']}",
                    file=sys.stderr,
                )
            pos += best['length']
        else:
            leftover.append({'pos': pos, 'base': consensus[pos]})
            if verbose:
                print(f"    pos {pos:>5}  UNMATCHED: {consensus[pos]}", file=sys.stderr)
            pos += 1

    # Consolidate runs of unmatched bases
    unmatched_blocks = []
    if leftover:
        blk = {'pos': leftover[0]['pos'], 'sequence': leftover[0]['base']}
        for seg in leftover[1:]:
            if seg['pos'] == blk['pos'] + len(blk['sequence']):
                blk['sequence'] += seg['base']
            else:
                unmatched_blocks.append(blk)
                blk = {'pos': seg['pos'], 'sequence': seg['base']}
        unmatched_blocks.append(blk)

    # ── Post-correction: terminal X-33_34ins → 6/6+ + 7/7+ ─────────────
    matched, unmatched_blocks = fix_terminal_6_7(
        matched, unmatched_blocks, consensus, max_mismatch
    )

    return matched, unmatched_blocks


# ══════════════════════════════════════════════════════════════════════════════
# Deconvolution of same-length haplotypes (size-homozygous case)
# ══════════════════════════════════════════════════════════════════════════════

def is_same_length_case(
    read_counts: dict,
    top_contigs: list,
    homo_threshold: float = 0.65,
    min_reads_dominant: int = 5,
) -> bool:
    """
    Detects whether the two haplotypes probably have the same VNTR length.

    Criterion: the ratio reads(top-2) / reads(top-1) is strictly below 0.25,
    i.e. the second best-represented contig receives less than a quarter
    of the reads of the first. This ratio is robust even when reads are scattered
    across many minority contigs (which would lower the global fraction of
    top-1 and mislead a criterion based on the total).

    A minimum of min_reads_dominant reads on the top-1 contig is required to
    avoid calls on too poorly covered data.

    The homo_threshold parameter is kept for CLI compatibility but is no
    longer used in the computation (the ratio is more robust).
    """
    if len(top_contigs) < 2:
        return False
    c1, c2 = top_contigs[0], top_contigs[1]
    n1 = read_counts.get(c1, 0)
    n2 = read_counts.get(c2, 0)
    if n1 < min_reads_dominant:
        return False
    return (n2 / n1) < 0.25


def find_phasing_snps(
    bam_path: str,
    contig: str,
    ref_len: int,
    min_bq: int = 20,
    min_mq: int = 20,
    min_depth: int = 5,
    min_af: float = 0.15,
    max_af: float = 0.85,
    verbose: bool = False,
) -> list:
    """
    Identifies bi-allelic positions useful for phasing.

    Returns a list of positions (int) where the minor allele frequency
    is between min_af and max_af with depth ≥ min_depth.
    These positions serve as SNP markers to separate the two haplotypes.
    """
    # fetch()+CIGAR (was a whole-contig pileup — same deep-data hang class as the consensus builders).
    # SNP phasing looks only at MATCHED bases: drop the deletion markers ('-') the tally records, exactly
    # as the old pileup skipped `is_del`. Insertions are irrelevant here and ignored.
    bases_all, _ = _tally_bases_cigar(bam_path, contig, ref_len, min_bq, min_mq)

    snp_positions = []
    for p, blist in bases_all.items():
        bases = [b for b in blist if b != '-']
        depth = len(bases)
        if depth < min_depth:
            continue
        ctr   = Counter(bases)
        top2  = ctr.most_common(2)
        if len(top2) < 2:
            continue
        af_minor = top2[1][1] / depth
        if min_af <= af_minor <= max_af:
            snp_positions.append(p)
            if verbose:
                print(
                    f"  [SNP] pos {p:>6}  depth={depth:>4}  "
                    f"{top2[0][0]}={top2[0][1]}  {top2[1][0]}={top2[1][1]}  "
                    f"AF_minor={af_minor:.2f}",
                    file=sys.stderr,
                )
    return sorted(snp_positions)


def _min_variance_split_threshold(scores: list) -> float:
    """PURE, O(n log n): the threshold that splits `scores` (per-read minor-allele fractions) into two
    groups minimizing the total within-group variance — the 1-D 2-means gap used to phase reads into two
    haplotypes. Candidate thresholds are `s + 0.001` for each score in ascending order; prefix sums of
    s and s² give each split's variance in O(1) and `bisect_left` locates the split point. The result is
    IDENTICAL to the naive O(n²) scan it replaces (same thresholds, same `var = Σs² − (Σs)²/n`, same
    strict `<` tie-break) — but that scan was quadratic and HUNG on a deep HiFi smear. Returns 0.5 when
    no non-trivial split exists (fewer than 2 scores, or every threshold leaves one side empty)."""
    s = sorted(scores)
    n = len(s)
    if n < 2:
        return 0.5
    P1 = [0.0] * (n + 1)
    P2 = [0.0] * (n + 1)
    for i, v in enumerate(s):
        P1[i + 1] = P1[i] + v
        P2[i + 1] = P2[i] + v * v
    best_gap, best_var = 0.5, float('inf')
    for v in s[:-1]:
        t = v + 0.001
        k = bisect.bisect_left(s, t)                        # |grpA| = #{x < t}
        if k == 0 or k == n:
            continue
        var_a = P2[k] - P1[k] * P1[k] / k
        var_b = (P2[n] - P2[k]) - (P1[n] - P1[k]) ** 2 / (n - k)
        if var_a + var_b < best_var:
            best_var, best_gap = var_a + var_b, t
    return best_gap


def phase_reads(
    bam_path: str,
    contig: str,
    ref_len: int,
    snp_positions: list,
    min_bq: int = 20,
    min_mq: int = 20,
    verbose: bool = False,
) -> tuple[list, list]:
    """
    Separates reads into two groups (haplotype A and haplotype B) by SNP phasing.

    Algorithm:
      1. For each read, extract the bases at the covered SNP positions.
      2. Compute a binary vector: major_allele → 0, minor_allele → 1.
      3. Simple hierarchical clustering (initialized on the most extreme read,
         then assigning remaining reads to the nearest centroid).
      4. Return (read_names_A, read_names_B).

    Returns two lists of read names (query_name).
    """
    # ── 1. Determine the major/minor alleles at each SNP position ──
    allele_major = {}   # pos → major base
    allele_minor = {}   # pos → minor base
    bases_at_snp  = defaultdict(list)

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for pcol in bam.pileup(
            contig, 0, ref_len,
            min_base_quality=min_bq,
            min_mapping_quality=min_mq,
            stepper='all', truncate=True,
        ):
            p = pcol.reference_pos
            if p not in snp_positions:
                continue
            for pr in pcol.pileups:
                if pr.is_refskip or pr.is_del:
                    continue
                b = pr.alignment.query_sequence[pr.query_position].upper()
                bases_at_snp[p].append(b)

    for p in snp_positions:
        bases = bases_at_snp[p]
        if not bases:
            continue
        top2 = Counter(bases).most_common(2)
        allele_major[p] = top2[0][0]
        allele_minor[p] = top2[1][0] if len(top2) > 1 else top2[0][0]

    informative = [p for p in snp_positions if p in allele_major and allele_major[p] != allele_minor[p]]
    if not informative:
        return [], []

    # ── 2. Extract each read's vector ──────────────────────────
    read_vectors = {}  # query_name → list of (pos, value) where value=0 major, 1=minor

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam.fetch(contig):
            if (read.is_unmapped or read.is_secondary or read.is_supplementary
                    or read.is_duplicate or read.mapping_quality < min_mq):
                continue
            seq = read.query_sequence
            if seq is None:
                continue
            # ref_pos → query index, built ONCE per read. The old code recomputed
            # get_reference_positions() and a linear scan for EACH SNP → O(SNPs·read_len) per
            # read, which hangs on deep HiFi (14k long reads × many noisy SNPs). Keep the FIRST
            # query index for a ref pos to match the old `next(...)` semantics (equivalent anyway:
            # full_length ref positions are monotone/unique). Result is byte-identical.
            try:
                qpos = read.get_reference_positions(full_length=True)
            except Exception:
                continue
            ref2q = {}
            for i, rp in enumerate(qpos):
                if rp is not None and rp not in ref2q:
                    ref2q[rp] = i
            vec = []
            for p in informative:
                idx = ref2q.get(p)
                if idx is None:
                    continue
                b = seq[idx].upper()
                if b == allele_major[p]:
                    vec.append((p, 0))
                elif b == allele_minor[p]:
                    vec.append((p, 1))
            if vec:
                read_vectors[read.query_name] = vec

    if not read_vectors:
        return [], []

    # ── 3. Simple clustering ─────────────────────────────────────────────
    # Represent each read by its fraction of minor alleles over the covered SNPs
    def snp_score(vec):
        return sum(v for _, v in vec) / len(vec) if vec else 0.5

    scored = {name: snp_score(vec) for name, vec in read_vectors.items()}

    # Initialize the two centroids on the most extreme reads
    sorted_reads = sorted(scored.items(), key=lambda x: x[1])
    if len(sorted_reads) < 2:
        return [r for r, _ in sorted_reads], []

    scores = [s for _, s in sorted_reads]                  # sorted ascending (sorted_reads is)

    # Natural gap that minimizes the within-group variance. Extracted to a PURE O(n log n) helper —
    # the old inline scan was O(n²) and hung on a deep HiFi smear (thousands of phased reads).
    best_gap = _min_variance_split_threshold(scores)

    reads_A = [n for n, s in scored.items() if s <  best_gap]
    reads_B = [n for n, s in scored.items() if s >= best_gap]

    # Reads without an SNP vector → assign to the larger group
    all_read_names = set()
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam.fetch(contig):
            if (read.is_unmapped or read.is_secondary or read.is_supplementary
                    or read.is_duplicate or read.mapping_quality < min_mq):
                continue
            all_read_names.add(read.query_name)

    unphased = all_read_names - set(read_vectors.keys())
    if unphased:
        if len(reads_A) >= len(reads_B):
            reads_A.extend(unphased)
        else:
            reads_B.extend(unphased)

    if verbose:
        print(
            f"  [phasing] {len(informative)} informative SNPs  "
            f"→ {len(reads_A)} reads hap-A / {len(reads_B)} reads hap-B  "
            f"({len(unphased)} unphased reads distributed)",
            file=sys.stderr,
        )

    return reads_A, reads_B


def build_consensus_from_readset(
    bam_path: str,
    ref_fa: pysam.FastaFile,
    contig: str,
    read_names: set,
    min_depth: int = 3,
    min_bq: int = 20,
    min_mq: int = 20,
    ins_threshold: float = 0.5,
    long_ins_threshold: float = 0.25,
    del_threshold: float = 0.6,
    verbose: bool = False,
) -> str:
    """
    Identical to build_consensus_for_contig but restricted to a subset
    of reads (identified by their query_name).

    Used to build the separate consensuses after intra-contig phasing.
    Positions without sufficient coverage in the subset are
    filled in from the reference sequence.
    """
    INS_SHORT_MAX = 3
    ref_seq = ref_fa.fetch(contig).upper()
    ref_len = len(ref_seq)

    # Restricted to the `read_names` subset — fetch()+CIGAR touches only those reads (was a whole-contig
    # pileup that processed EVERY read then filtered, the deep-HiFi hang).
    bases_at, ins_after = _tally_bases_cigar(bam_path, contig, ref_len, min_bq, min_mq,
                                             read_names=read_names)

    consensus_parts = []
    low_cov = 0

    for pos in range(ref_len):
        bases = bases_at[pos]
        depth = len(bases)

        if depth < min_depth:
            consensus_parts.append(ref_seq[pos])
            low_cov += 1
        else:
            ctr = Counter(bases)
            top_base, top_count = ctr.most_common(1)[0]
            if top_base == '-' and top_count / depth >= del_threshold:
                pass
            else:
                non_del = [b for b in bases if b != '-']
                if non_del:
                    consensus_parts.append(Counter(non_del).most_common(1)[0][0])
                elif top_base != '-':
                    consensus_parts.append(top_base)

        ins_list = ins_after.get(pos, [])
        if ins_list:
            ctr_ins   = Counter(ins_list)
            best_ins, best_count = ctr_ins.most_common(1)[0]
            ins_len   = len(best_ins)
            if ins_len <= INS_SHORT_MAX:
                if ins_len == 1:
                    for delta in (-1, 1):
                        adj = ins_after.get(pos + delta, [])
                        best_count += sum(1 for s in adj if s == best_ins)
                    ref_depth = max(
                        len(bases_at.get(max(0, pos - 1), [])),
                        depth,
                        len(bases_at.get(min(ref_len - 1, pos + 1), [])),
                    )
                    frac = best_count / max(ref_depth, 1)
                else:
                    frac = best_count / max(depth, 1)
                threshold = ins_threshold
            else:
                frac      = best_count / max(depth, 1)
                threshold = long_ins_threshold
            if frac >= threshold:
                consensus_parts.append(best_ins)

    consensus = ''.join(consensus_parts)
    if verbose:
        pct = 100 * low_cov / ref_len if ref_len else 0
        print(
            f"  [consensus subset] {contig}: {len(consensus)} bp "
            f"({low_cov} low-coverage pos, {pct:.1f}%)",
            file=sys.stderr,
        )
    return consensus


# ══════════════════════════════════════════════════════════════════════════════
# Detection and export of candidate novel motifs (60 bp, approximate matches)
# ══════════════════════════════════════════════════════════════════════════════

def collect_novel_60bp_motifs(all_results: list) -> list:
    """
    Collects the 60 bp sequences assigned to a known motif with at least
    1 mismatch (approximate matches). These are candidates for new motifs
    not yet described in KNOWN_REPEATS.

    For each distinct sequence:
      - closest_motif  : known motif to which it was assigned
      - closest_dist   : number of mismatches (Hamming or Levenshtein distance)
      - diff_positions : list of positions (0-indexed) and bases differing from the
                         closest motif, format "pos:consensus>motif"
      - occurrences    : total number of appearances across all haplotypes/samples
      - haplotypes     : haplotypes in which this sequence was observed
      - motif_positions : positions in the VNTR in "rank/total" format
                          (e.g. "44/64"), one entry per occurrence

    The result is sorted by decreasing occurrences then increasing distance.
    """
    seen: dict = {}  # sequence → aggregated dict

    for h in all_results:
        hap_label  = f"hap{h['rank']}_{h['contig']}"
        total_motifs = len(h['motifs'])
        for motif_idx, m in enumerate(h['motifs']):
            # Only the 60 bp approximate matches
            if m['mismatches'] == 0 or m['length'] != 60:
                continue
            seq      = m['sequence']
            position = f"{motif_idx + 1}/{total_motifs}"

            if seq not in seen:
                # Compute the exact difference positions against the assigned motif
                ref_seq = next(
                    (s for s, n in KNOWN_REPEATS.items() if n == m['name']), None
                )
                diff_pos = []
                if ref_seq and len(ref_seq) == 60:
                    diff_pos = [
                        f"{i}:{seq[i]}>{ref_seq[i]}"
                        for i in range(60)
                        if seq[i] != ref_seq[i]
                    ]
                seen[seq] = {
                    'sequence':        seq,
                    'closest_motif':   m['name'],
                    'closest_dist':    m['mismatches'],
                    'diff_positions':  diff_pos,
                    'occurrences':     0,
                    'haplotypes':      [],
                    'motif_positions': [],
                }
            rec = seen[seq]
            rec['occurrences'] += 1
            rec['motif_positions'].append(position)
            # Keep the match with the smallest distance
            if m['mismatches'] < rec['closest_dist']:
                rec['closest_dist']  = m['mismatches']
                rec['closest_motif'] = m['name']
                ref_seq = next(
                    (s for s, n in KNOWN_REPEATS.items() if n == m['name']), None
                )
                if ref_seq and len(ref_seq) == 60:
                    rec['diff_positions'] = [
                        f"{i}:{seq[i]}>{ref_seq[i]}"
                        for i in range(60)
                        if seq[i] != ref_seq[i]
                    ]
            if hap_label not in rec['haplotypes']:
                rec['haplotypes'].append(hap_label)

    candidates = list(seen.values())
    candidates.sort(key=lambda r: (-r['occurrences'], r['closest_dist']))
    return candidates


def write_novel_motifs(candidates: list, path: str, sample: str = '') -> None:
    """
    Writes the TSV file of candidate novel motifs (60 bp, approximate matches).

    Columns:
      sample           Sample name
      occurrences      Number of times the exact sequence was observed
      haplotypes       Haplotypes containing this sequence (separated by ';')
      motif_positions  Position(s) in the VNTR in "rank/total" format
                       (e.g. "44/64"), separated by ';' if multiple occurrences
      closest_motif    Known motif to which the sequence was assigned
      closest_dist     Number of mismatches with that motif
      diff_positions   Positions and bases differing from the motif (e.g. "12:A>G;45:C>T")
      sequence         Candidate nucleotide sequence (60 bp)
    """
    header = '\t'.join([
        'sample', 'occurrences', 'haplotypes', 'motif_positions',
        'closest_motif', 'closest_dist', 'diff_positions', 'sequence',
    ])
    lines = [header]
    for c in candidates:
        lines.append('\t'.join([
            sample or '',
            str(c['occurrences']),
            ';'.join(c['haplotypes']),
            ';'.join(c['motif_positions']),
            c['closest_motif'],
            str(c['closest_dist']),
            ';'.join(c['diff_positions']),
            c['sequence'],
        ]))
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


# ══════════════════════════════════════════════════════════════════════════════
# Export FASTA consensus
# ══════════════════════════════════════════════════════════════════════════════

def write_fasta_consensus(results: list, path: str, sample: str = '') -> None:
    """
    Writes the consensus sequences of each haplotype to a FASTA file.
    The header uses exactly the contig name (e.g. MUC1_VNTR_77repeats),
    prefixed with the sample name if provided.
    """
    with open(path, 'w') as fh:
        for h in results:
            header = h['contig']
            if sample:
                header = f"{sample}_{header}"
            fh.write(f">{header}\n")
            seq = h['consensus']
            for i in range(0, len(seq), 80):
                fh.write(seq[i: i + 80] + '\n')


# ══════════════════════════════════════════════════════════════════════════════
# PDF visual report export
# ══════════════════════════════════════════════════════════════════════════════

# Names of motifs carrying an indel (length ≠ 60 bp) — used for
# differential coloring in the PDF.
_INDEL_MOTIF_NAMES = frozenset(
    name for seq, name in KNOWN_REPEATS.items() if len(seq) != 60
)
# ADTKD-MUC1 is caused by a FRAMESHIFT: only a unit whose length differs from 60 by a NON-multiple of 3
# shifts the reading frame and creates the pathogenic MUC1-fs neoprotein. An in-frame indel (Δ divisible
# by 3) changes the protein locally but keeps the frame → it is NOT a MUC1-fs variant and must never drive
# the onset axis. (Today exactly one dictionary unit is in-frame: X-33_34insCGGGCCCCGGGCTCCACC, +18 bp —
# which is also the recurring partial-motif-duplication consensus artefact seen on unrelated samples.)
_FRAMESHIFT_MOTIF_NAMES = frozenset(
    name for seq, name in KNOWN_REPEATS.items() if (len(seq) - 60) % 3 != 0
)
_INFRAME_INDEL_MOTIF_NAMES = _INDEL_MOTIF_NAMES - _FRAMESHIFT_MOTIF_NAMES


def _indel_summary(motifs: list) -> str:
    """
    Builds the parenthesized part of the molecular result line — FRAMESHIFTS ONLY (the pathogenic class,
    which drives the onset axis). E.g.: '[X-59dupC ~ repeat 12]' or empty string if none.
    Multiple frameshifts are separated by ' ; '. In-frame indels are reported by `_inframe_summary`.
    """
    parts = []
    for i, m in enumerate(motifs, 1):
        if m['name'] in _FRAMESHIFT_MOTIF_NAMES:
            parts.append(f"{m['name']} ~ repeat {i}")
    return ' ; '.join(parts)


def _inframe_summary(motifs: list) -> str:
    """In-frame indel units (Δ length divisible by 3) — reported for transparency but NOT as MUC1-fs."""
    parts = []
    for i, m in enumerate(motifs, 1):
        if m['name'] in _INFRAME_INDEL_MOTIF_NAMES:
            parts.append(f"{m['name']} ~ repeat {i}")
    return ' ; '.join(parts)


def write_pdf_report(results: list, path: str, sample: str = '', two_axis: dict = None,
                     length_warning: str = None, arbiter: dict = None, rs4072037: dict = None,
                     dupc: dict = None, congruence: dict = None, clinical_call: str = None,
                     onset_index: float = None, variant_repeat: int = None,
                     variant_label: str = None, carrier_contig: str = None,
                     score: dict = None, rs_alert: str = None) -> None:
    """
    Generates a PDF report with, for each haplotype:
      - The contig name
      - The colored motif nomenclature (indels in orange)
      - At the bottom: the final molecular result in the form
          X motifs ([indel] ~ repeat m) | Y motifs ([indel] ~ repeat n)
    """
    try:
        from reportlab.lib.pagesizes import A4
        from reportlab.lib import colors
        from reportlab.lib.units import cm
        from reportlab.lib.styles import ParagraphStyle
        from reportlab.platypus import (
            SimpleDocTemplate, Paragraph, Spacer, HRFlowable, KeepTogether
        )
        from reportlab.lib.enums import TA_LEFT, TA_CENTER
    except ImportError:
        print("[ERROR] reportlab not available. Install it with: pip install reportlab",
              file=sys.stderr)
        return

    # ── Styles ────────────────────────────────────────────────────────────
    COLOR_NORMAL = '#1a1a2e'
    COLOR_INDEL  = '#c0392b'   # brick red for motifs with an indel
    COLOR_HEADER = '#2c3e6b'
    COLOR_RULE   = colors.HexColor('#2c3e6b')

    style_title = ParagraphStyle(
        'title',
        fontName='Helvetica-Bold',
        fontSize=14,
        textColor=colors.HexColor(COLOR_HEADER),
        spaceAfter=6,
        alignment=TA_CENTER,
    )
    style_contig = ParagraphStyle(
        'contig',
        fontName='Helvetica-Bold',
        fontSize=11,
        textColor=colors.HexColor(COLOR_HEADER),
        spaceBefore=14,
        spaceAfter=4,
    )
    style_nomen = ParagraphStyle(
        'nomen',
        fontName='Helvetica',
        fontSize=7.5,
        leading=11,
        textColor=colors.HexColor(COLOR_NORMAL),
        spaceAfter=8,
        wordWrap='CJK',   # coupe partout, pas seulement aux espaces
    )
    style_result = ParagraphStyle(
        'result',
        fontName='Helvetica-Bold',
        fontSize=10,
        textColor=colors.HexColor(COLOR_HEADER),
        spaceBefore=18,
        spaceAfter=6,
        alignment=TA_CENTER,
        borderPad=6,
        backColor=colors.HexColor('#eaf0fb'),
        borderColor=colors.HexColor(COLOR_HEADER),
        borderWidth=1,
        borderRadius=4,
    )
    style_meta = ParagraphStyle(
        'meta',
        fontName='Helvetica',
        fontSize=8,
        textColor=colors.HexColor('#555555'),
        spaceAfter=2,
    )
    style_section = ParagraphStyle(
        'section',
        fontName='Helvetica-Bold',
        fontSize=13,
        textColor=colors.HexColor(COLOR_HEADER),
        spaceBefore=6,
        spaceAfter=6,
    )

    # ── Document construction ──────────────────────────────────────────
    doc = SimpleDocTemplate(
        path,
        pagesize=A4,
        leftMargin=2 * cm,
        rightMargin=2 * cm,
        topMargin=2 * cm,
        bottomMargin=2 * cm,
    )
    story = []

    # Title
    import re as _re
    _disp_sample = _re.sub(r'[-_]+vntrcigarlength$', '', sample or '').rstrip('-_ ')
    title_text = f"VNTR Report — {_disp_sample}" if _disp_sample else "VNTR Report"
    story.append(Paragraph(title_text, style_title))
    story.append(HRFlowable(width='100%', thickness=1.5,
                             color=COLOR_RULE, spaceAfter=10))

    # ── Clinical call, boxed at the top — the result the reader came for ─────
    if clinical_call:
        try:
            from .report_blocks import clinical_call_box
            story.append(clinical_call_box(clinical_call))
            story.append(Spacer(1, 0.35 * cm))
        except Exception as e:
            print(f"[WARN] clinical box skipped: {e}", file=sys.stderr)

    # ── Reliable VNTR length (alignment-free arbiter) — the number to read ───
    if arbiter and arbiter.get("available"):
        _al, _ct = arbiter.get("alleles", []), arbiter.get("counts", [])
        _lc = " — LONG allele low-confidence / PCR-depleted" if arbiter.get("long_low_conf") else ""
        if len(_al) == 2:
            _body = f"{_al[0]} / {_al[1]} copies (n={_ct[0]}/{_ct[1]}){_lc}"
        elif len(_al) == 1:
            _body = f"~{_al[0]} copies (homozygous)"
        else:
            _body = "no clear peak — too shallow"
        ok_style = ParagraphStyle('arblen', fontName='Helvetica-Bold', fontSize=10,
                                  textColor=colors.HexColor('#1b5e20'), spaceAfter=4,
                                  backColor=colors.HexColor('#e8f5e9'), borderPad=6,
                                  borderColor=colors.HexColor('#1b5e20'), borderWidth=1)
        story.append(Paragraph(f"VNTR LENGTH — {arbiter['method']} (alignment-free): {_body}", ok_style))
        # technology + (if a PCR amplicon) the product band size
        _plat = arbiter.get("platform", "unknown")
        if arbiter.get("is_amplicon"):
            _prod = arbiter.get("product_bp", [])
            _bands = " / ".join(f"~{round(p / 1000, 1)} kb" for p in _prod) if _prod else "?"
            _ctx = f"Technology: {_plat}  ·  PCR amplicon detected — product {_bands} (flanks + VNTR)"
        else:
            _ctx = f"Technology: {_plat}  ·  adaptive sampling / WGS (not a PCR amplicon)"
        story.append(Paragraph(_ctx, ParagraphStyle('arbctx', fontName='Helvetica', fontSize=8,
                     textColor=colors.HexColor('#555555'), spaceAfter=4)))
    # ── rs4072037 and the dupC verdict ──────────────────────────────────────
    # NOT nested under the arbiter block, where they used to live: a sample whose length could not be
    # called would then lose its VARIANT verdict from the report entirely — the one line a clinician
    # cannot do without. The length failing says nothing about whether the dupC was found.
    _rs = _rs4072037_line(rs4072037)
    if _rs:
        story.append(Paragraph(_rs, ParagraphStyle('arbrs', fontName='Helvetica', fontSize=9,
                     textColor=colors.HexColor('#1b5e20'), spaceAfter=4)))
    # power-aware wording when the verdict carries an interpretation; else our compact line. A bare
    # "not called" cannot tell an adequately-powered negative from a sample too thin to exclude.
    _dc = ""
    try:
        from .report_blocks import dupc_verdict_html
        _dc = dupc_verdict_html((dupc or {}).get("result", dupc))
    except Exception:
        _dc = ""
    _dc = _dc or _dupc_line(dupc)
    if _dc:
        _res = (dupc or {}).get("result", dupc) or {}
        _pos = str(_res.get("status", "")).lower() == "positive" or bool(_res.get("called")) \
            or _res.get("interpretation") == "CONFIRMED"
        story.append(Paragraph(_dc, ParagraphStyle('arbdc', fontName='Helvetica-Bold', fontSize=9,
                     textColor=colors.HexColor('#8a1f11' if _pos else '#1b5e20'), spaceAfter=10)))

    # ── Layer congruence ────────────────────────────────────────────────────
    # The cross-layer check used to reach the terminal and the JSON only. A DIVERGENCE the clinician never
    # sees is a divergence that does not exist, so the PDF carries it too — same lines, one renderer.
    _cong_lines = []
    try:
        from .congruence import render_lines as _cong_render
        _cong_lines = _cong_render(congruence)
    except Exception:                                    # never let the report die on its own footer
        _cong_lines = []
    if _cong_lines:
        _COL = {"head": '#2c3e6b', "agree": '#1b5e20', "disagree": '#8a1f11', "info": '#4a4a4a'}
        _block = []
        for _i, (_txt, _kind) in enumerate(_cong_lines):
            _block.append(Paragraph(
                _txt,
                ParagraphStyle(f'cong{_i}', fontName='Helvetica-Bold' if _i == 0 else 'Helvetica',
                               fontSize=9 if _i == 0 else 8, leftIndent=0 if _i == 0 else 10,
                               textColor=colors.HexColor(_COL.get(_kind, '#4a4a4a')), spaceAfter=2)))
        story.append(KeepTogether(_block))
        story.append(Spacer(1, 0.25 * cm))

    # ── Length-unreliable banner (multi-contig smear / thin support) ─────────
    if length_warning:
        warn_style = ParagraphStyle('warn', fontName='Helvetica-Bold', fontSize=9,
                                    textColor=colors.HexColor('#8a1f11'), spaceAfter=10,
                                    backColor=colors.HexColor('#fdecea'), borderPad=6,
                                    borderColor=colors.HexColor('#8a1f11'), borderWidth=1)
        story.append(Paragraph("⚠ LENGTH UNRELIABLE — " + length_warning, warn_style))

    # ── One block per haplotype ─────────────────────────────────────────────
    story.append(Paragraph("Haplotypes", style_section))
    for h in results:
        block = []

        # Haplotype header
        block.append(Paragraph(
            f"Haplotype {h['rank']} — <font color='{COLOR_HEADER}'>{h['contig']}</font>",
            style_contig,
        ))
        block.append(Paragraph(
            f"{h['read_count']} reads  |  {h['contig_length']} bp reference  |  "
            f"{h['n_motifs']} motifs identified "
            f"(exact: {h['n_exact']} / approximate: {h['n_approx']})",
            style_meta,
        ))

        # Colored nomenclature: normal motifs in black, indels in red, and — on the CARRIER haplotype —
        # the pathogenic unit itself picked out, so the eye lands on the repeat that causes the disease
        # instead of on every indel-looking motif equally.
        _is_carrier = bool(carrier_contig) and h.get('contig') == carrier_contig
        tokens = []
        for i, m in enumerate(h['motifs']):
            sep = '-' if i > 0 else ''
            if _is_carrier and variant_repeat is not None and i + 1 == variant_repeat:
                _disp = m['name'] if (variant_label and variant_label in m['name']) \
                    else f"{m['name']}({variant_label})"
                tokens.append(f"{sep}<font color='#C00000'><b>{_disp}</b></font>")
            elif m['name'] in _INDEL_MOTIF_NAMES:
                tokens.append(
                    f"{sep}<font color='{COLOR_INDEL}'><b>{m['name']}</b></font>"
                )
            else:
                tokens.append(f"{sep}{m['name']}")
        nomen_html = ''.join(tokens)
        block.append(Paragraph(nomen_html, style_nomen))

        story.append(KeepTogether(block))
        story.append(HRFlowable(width='100%', thickness=0.5,
                                  color=colors.HexColor('#cccccc'), spaceAfter=4))

    # ── Final molecular-result line — only when there is NO clinical-call box above (else the same
    #     genotype would appear twice, in two blue boxes; the top box is the one the reader came for) ──
    if not clinical_call:
        parts_mol = []
        for h in results:
            n = h['n_motifs']
            indel_str = _indel_summary(h['motifs'])
            parts_mol.append(f"{n} motifs ({indel_str})" if indel_str else f"{n} motifs")
        story.append(Spacer(1, 0.3 * cm))
        story.append(Paragraph(" | ".join(parts_mol), style_result))

    # ── MUC1 Score — the onset / severity detail table (shown when the caller forwards its score) ────
    if score:
        from reportlab.platypus import Table, TableStyle
        story.append(Spacer(1, 0.15 * cm))
        story.append(HRFlowable(width='100%', thickness=0.5, color=colors.HexColor('#cccccc'), spaceAfter=6))
        story.append(Paragraph("MUC1 Score", style_section))
        _sev = score.get("severity_score")
        if rs4072037 and rs4072037.get("genotype"):
            _src = rs4072037.get("note") or f"rs4072037-{rs4072037.get('genotype')}"
        else:
            _src = "provided"
        if _sev is not None:
            _sev_txt = (f"{_sev:+g} (rs4072037-{score.get('splice_base')} → "
                        f"{'severe' if _sev > 0 else 'protective'}; {_src})")
        else:
            _sev_txt = "n/a — rs4072037 undetermined (no usable read and no --rs4072037)"
        _rows = [["MUC1_Score — ONSET axis",
                  f"{score.get('onset_score')}  (onset_index {score.get('onset_index')})"],
                 ["MUC1_Score — SEVERITY axis", _sev_txt],
                 ["VNTR ratio (mut/healthy)", f"{score.get('ratio')}"]]
        _st = Table(_rows, colWidths=[6 * cm, 11 * cm])
        _st.setStyle(TableStyle([("FONT", (0, 0), (-1, -1), "Helvetica", 9.5),
                                 ("FONT", (0, 0), (0, -1), "Helvetica-Bold", 9.5),
                                 ("LINEBELOW", (0, 0), (-1, -2), 0.4, colors.HexColor("#DDDDDD")),
                                 ("TOPPADDING", (0, 0), (-1, -1), 4), ("BOTTOMPADDING", (0, 0), (-1, -1), 4)]))
        story.append(_st)
        story.append(Spacer(1, 0.3 * cm))

    # ── Two-axis clinical summary (optional; the genomic pipeline's own prognosis block) ─────────────
    if two_axis is not None:
        from . import two_axis_report as TA
        story.extend(TA.render_pdf_blocks(two_axis, {"result": style_result, "meta": style_meta}))

    # Where this sample's onset sits in the calibration cohort. Kept OUTSIDE the two-axis block: the
    # figure is about the onset axis, which predicts better than the absolute length, and it should not
    # disappear because some other summary field was missing.
    _oi = onset_index if onset_index is not None else (two_axis or {}).get("onset_index")
    if _oi is not None:
        try:
            from .report_blocks import cohort_plot
            _plot = cohort_plot(_oi)
            if _plot is not None:
                story.append(Spacer(1, 0.2 * cm))
                story.append(_plot)
        except Exception as e:
            print(f"[WARN] cohort plot skipped: {e}", file=sys.stderr)

    # ── Genotype → phenotype flag, under the cohort curve (driven by the severity axis) ──────────────
    if score:
        _sev = score.get("severity_score")
        _geno = ParagraphStyle('geno', fontName='Helvetica-Bold', fontSize=13, spaceBefore=6, leading=17)
        if _sev is not None and _sev > 0:
            story.append(Paragraph('<font color="#C00000"><b>GENOTYPE ASSOCIATED WITH A SEVERE PHENOTYPE</b></font>', _geno))
        elif _sev is not None and _sev < 0:
            story.append(Paragraph('<font color="#1B7A2F"><b>GENOTYPE ASSOCIATED WITH A PROTECTIVE PHENOTYPE</b></font>', _geno))
        else:
            story.append(Paragraph('<font color="#888888"><b>SEVERITY UNDETERMINED — rs4072037 NOT DETERMINED</b></font>', _geno))
        if rs_alert:
            story.append(Spacer(1, 0.1 * cm))
            story.append(Paragraph(f'<font color="#B8860B"><b>&#9888; {rs_alert}</b></font>', style_meta))

    doc.build(story)


# ══════════════════════════════════════════════════════════════════════════════
# Text output formatting
# ══════════════════════════════════════════════════════════════════════════════

def haplotype_string(matched: list) -> str:
    """Return the compact nomenclature: A-C-Z-X-X-X-G-A-..."""
    return '-'.join(m['name'] for m in matched)


def format_report(results: list, sample: str = '', length_warning: str = None) -> str:
    """
    Format the full report for all haplotypes.
    `results` is a list of per-haplotype dicts.
    """
    lines = []
    sep   = '═' * 80

    if sample:
        lines.append(sep)
        lines.append(f"  Sample: {sample}")
    lines.append(sep)
    if length_warning:
        lines.append("  ⚠ LENGTH UNRELIABLE — " + length_warning)
        lines.append(sep)

    for i, h in enumerate(results, 1):
        lines.append(f"\n{'─'*80}")
        lines.append(
            f"  Haplotype {i} — contig: {h['contig']}  "
            f"({h['read_count']} reads, {h['contig_length']} bp reference)"
        )
        lines.append(f"{'─'*80}")
        lines.append(f"\n  Nomenclature:\n  {h['nomenclature']}\n")
        lines.append(
            f"  Motifs identified: {h['n_motifs']}  "
            f"(exact: {h['n_exact']}  |  approximate: {h['n_approx']})"
        )
        if h['unmatched']:
            lines.append(f"  ⚠ Unmatched bases ({len(h['unmatched'])} block(s)):")
            for b in h['unmatched']:
                lines.append(f"      pos {b['pos']:>5} : {b['sequence']}")

        lines.append(f"\n  Details:")
        lines.append(
            f"  {'Pos':>6}  {'Len.':>5}  {'Motif':<38}  {'Match':>6}  Sequence"
        )
        lines.append(f"  {'─'*6}  {'─'*5}  {'─'*38}  {'─'*6}  {'─'*20}")
        for m in h['motifs']:
            tag = 'exact' if m['mismatches'] == 0 else f"{m['mismatches']} mm"
            lines.append(
                f"  {m['pos']:>6}  {m['length']:>5}  {m['name']:<38}  "
                f"{tag:>6}  {m['sequence']}"
            )

    lines.append(f"\n{sep}")
    return '\n'.join(lines)


# ══════════════════════════════════════════════════════════════════════════════
# Entry point
# ══════════════════════════════════════════════════════════════════════════════

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description='VNTR Haplotype Caller — multi-contig BAM → motif nomenclature',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('-b', '--bam',     required=True, help='BAM file (.bai index required)')
    p.add_argument('-r', '--ref',     required=True, help='FASTA reference (one contig = one VNTR length)')
    p.add_argument('-o', '--output',  default='-',   help='Text output file (- = stdout)')
    p.add_argument('-s', '--sample',  default='',    help="Sample name")
    p.add_argument('--top-n',    type=int,   default=2,   help="Number of haplotypes to report")
    p.add_argument('--min-depth', type=int,  default=3,   help='Min depth to call a consensus base')
    p.add_argument('--min-bq',   type=int,   default=20,  help='Min base quality')
    p.add_argument('--min-mq',   type=int,   default=20,  help='Min mapping quality')
    p.add_argument('--max-mismatch', type=int, default=3, help='Max mismatches (approximate matching)')
    p.add_argument('--threads', type=int, default=1,
                   help='Number of processes for per-contig decomposition (default 1 = sequential)')
    p.add_argument('--max-depth', type=int, default=None,
                   help='Subsample to N reads/contig (deterministic + faster consensus on deep PCR)')
    p.add_argument('--pcr', action='store_true',
                   help='PCR amplicon mode: adapted defaults (--max-depth 200 if unset) + slippage note. '
                        'NB: the analyzer consumes an already-aligned BAM; the FROZEN PCR minimap params '
                        '(set B: -z 600,200 -r 2000,20000) live at ALIGNMENT — config.py::PCR_MINIMAP_FLAGS, '
                        'applied at the alignment step.')
    p.add_argument('--ins-threshold', type=float, default=0.5,
                   help='Min fraction to include an insertion ≤ 3 bp (homopolymer)')
    p.add_argument('--long-ins-threshold', type=float, default=0.4,
                   help='Min fraction to include a long insertion > 3 bp (e.g. X-33_34ins)')
    p.add_argument('--del-threshold', type=float, default=0.6, help='Min fraction to apply a deletion')
    p.add_argument('--dump-consensus', action='store_true', help='Write consensus FASTA to stderr')
    p.add_argument('--fasta-consensus', default=None,
                   help='Output FASTA file for the consensus sequences')
    p.add_argument('--pdf', default=None,
                   help='Output PDF file — visual report with colored nomenclature')
    p.add_argument('--json', default=None, help='JSON file for detailed export')
    p.add_argument('--homo-threshold', type=float, default=0.65,
                   help='Min fraction of reads on the dominant contig to trigger intra-contig phasing')
    p.add_argument('--snp-min-af', type=float, default=0.15,
                   help='Min minor-allele frequency for a phasing SNP')
    # ── Two-axis clinical summary (optional): arbiter length + rs4072037 per haplotype ──
    # The VNTR-ref caller cannot see rs4072037 (outside the tandem) and under-calls the long allele on
    # LR-PCR. Point --genomic-bam at the chr1 BAM/CRAM to add the arbiter length + rs4072037/hap section.
    p.add_argument('--genomic-bam', default=None,
                   help='chr1 genomic/PCR BAM or CRAM → adds a two-axis summary (arbiter length + '
                        'rs4072037 phased per haplotype) to the report/PDF/JSON')
    p.add_argument('--genome-ref', default=None, help='genome FASTA (required if --genomic-bam is a CRAM)')
    p.add_argument('--dupc-json', default=None,
                   help='reliable dupC verdict (dupc_dispatch JSON) → rendered in the report/PDF/JSON as THE '
                        'trustworthy dupC (distinct from the depth-fragile consensus X-59dupC token). `run` wires this.')
    p.add_argument('--report-json', default=None,
                   help='pcr_report JSON (clinical coordinates of a CONFIRMED positive) → enables the '
                        'LAYER CONGRUENCE block, which cross-checks variant/repeat/carrier-allele between '
                        'the independent layers instead of printing them side by side. `run` wires this.')
    p.add_argument('--snp-chrom', default='chr1', help='contig name of rs4072037 in the genomic BAM')
    p.add_argument('--snp-pos', type=int, default=None,
                   help='rs4072037 1-based position (default = GRCh38 155,192,276)')
    p.add_argument('--copy-offset', type=int, default=None,
                   help='flank→physical copy offset (default 4 ALWAYS: the AL/AH anchors sit ~4 units inside the array)')
    # Supply a VALIDATED mutation (lab nomenclature / careful per-allele caller run) for the ONSET axis —
    # authoritative when the top-N decomposition missed it (a PCR-depleted long mutant). The report only
    # displays it (attaches to the nearest allele + promotes it); it never re-runs a caller → no false call.
    p.add_argument('--mut-variant', default=None, help='validated mutation name, e.g. 59dupC / 58_59insG / del8_27')
    p.add_argument('--mut-repeat', type=int, default=None, help='repeat position of the supplied mutation')
    p.add_argument('--mut-allele-len', type=int, default=None,
                   help='copies of the mutant allele (attaches the supplied mutation to it)')
    # Two-axis COVERAGE-GATE (only with --genomic-bam): a plain-call frameshift decoded on a consensus
    # with < min-reads reads, or fewer than min-strand reads on either strand, is flagged unreliable.
    p.add_argument('--gate-min-reads', type=int, default=5,
                   help='two-axis onset coverage-gate: min reads backing the mutant consensus (default 5)')
    p.add_argument('--gate-min-strand', type=int, default=2,
                   help='two-axis onset coverage-gate: min reads PER STRAND (bi-strand, default 2)')
    p.add_argument('--novel-motifs', default=None,
                   help='TSV file listing the 60 bp sequences assigned with mismatches '
                        '— candidates for new motifs not yet described in KNOWN_REPEATS')
    p.add_argument('-v', '--verbose', action='store_true', help='Verbose mode')
    return p.parse_args(argv)


def _read_names_for_contig(bam_path: str, contig: str, min_mq: int) -> list:
    """Read names (unique, stable order) mapped on `contig` above `min_mq`."""
    seen = {}
    with pysam.AlignmentFile(bam_path) as bam:
        for r in bam.fetch(contig):
            if r.is_unmapped or r.mapping_quality < min_mq or r.query_name is None:
                continue
            seen[r.query_name] = None
    return list(seen)


def _strand_counts_for_contig(bam_path: str, contig: str, min_mq: int) -> dict:
    """Bi-strand read support for a whole contig — {'n_reads','n_fwd','n_rev'} (unique PRIMARY reads).

    Feeds the two-axis COVERAGE-GATE: a frameshift decoded on a consensus with too few reads or only
    one strand is a decode artefact. Uses the SAME filter as count_reads_per_contig / the consensus
    builders (drops secondary/supplementary/duplicate) so the count matches the reads that actually
    built the consensus — in a multi-contig smear, fetch() returns a read's SECONDARY alignments on
    every length-contig, which would otherwise inflate the support and pollute the strand tally."""
    seen = {}
    with pysam.AlignmentFile(bam_path) as bam:
        for r in bam.fetch(contig):
            if (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate
                    or r.mapping_quality < min_mq or r.query_name is None):
                continue
            seen[r.query_name] = r.is_reverse
    n_rev = sum(1 for v in seen.values() if v)
    return {"n_reads": len(seen), "n_fwd": len(seen) - n_rev, "n_rev": n_rev}


def _strand_counts_for_readset(bam_path: str, contig: str, names, min_mq: int) -> dict:
    """Bi-strand read support restricted to `names` (a phased read group) on `contig` — PRIMARY only."""
    names = set(names)
    seen = {}
    with pysam.AlignmentFile(bam_path) as bam:
        for r in bam.fetch(contig):
            if (r.query_name not in names or r.is_unmapped or r.is_secondary
                    or r.is_supplementary or r.is_duplicate or r.mapping_quality < min_mq):
                continue
            seen[r.query_name] = r.is_reverse
    n_rev = sum(1 for v in seen.values() if v)
    return {"n_reads": len(seen), "n_fwd": len(seen) - n_rev, "n_rev": n_rev}


def _process_contig(a: dict) -> dict:
    """Worker (picklable): consensus + matching for ONE contig. Opens its own pysam.
    If `max_depth` is set, subsamples the reads (deterministic) via
    build_consensus_from_readset — fast/consistent on deep PCR. The result has the
    same shape as the sequential branch (all_results)."""
    import random
    ref_fa = pysam.FastaFile(a['ref'])
    try:
        md = a.get('max_depth')
        if md:
            names = _read_names_for_contig(a['bam'], a['contig'], a['min_mq'])
            names = (set(random.Random(1234).sample(names, md)) if len(names) > md else set(names))
            consensus = build_consensus_from_readset(
                bam_path=a['bam'], ref_fa=ref_fa, contig=a['contig'], read_names=names,
                min_depth=a['min_depth'], min_bq=a['min_bq'], min_mq=a['min_mq'],
                ins_threshold=a['ins_threshold'], long_ins_threshold=a['long_ins_threshold'],
                del_threshold=a['del_threshold'], verbose=False)
        else:
            consensus = build_consensus_for_contig(
                bam_path=a['bam'], ref_fa=ref_fa, contig=a['contig'],
                min_depth=a['min_depth'], min_bq=a['min_bq'], min_mq=a['min_mq'],
                ins_threshold=a['ins_threshold'], long_ins_threshold=a['long_ins_threshold'],
                del_threshold=a['del_threshold'], verbose=False)
        matched, unmatched = match_motifs(consensus, max_mismatch=a['max_mismatch'], verbose=False)
        return {
            'rank': a['rank'], 'contig': a['contig'],
            'contig_length': ref_fa.get_reference_length(a['contig']),
            'read_count': a['read_count'], 'consensus': consensus,
            'nomenclature': haplotype_string(matched), 'n_motifs': len(matched),
            'n_exact': sum(1 for m in matched if m['mismatches'] == 0),
            'n_approx': sum(1 for m in matched if m['mismatches'] > 0),
            'motifs': matched, 'unmatched': unmatched, 'phased': False, 'same_length': False,
        }
    finally:
        ref_fa.close()


def _rs4072037_line(rs) -> str | None:
    """One-line rs4072037 genotype read VNTR-NATIVELY off the multi-contig bam (offset 4265, no GRCh38) —
    so it resolves from the caller's own BAM without a separate chr1/--genomic-bam. None if unavailable."""
    if not rs:
        return None
    gt = rs.get("genotype") or rs.get("gt_label") or "?"
    if rs.get("snp_genotype") is None:
        return f"rs4072037: undetermined ({gt}, depth {rs.get('depth', 0)})"
    src = rs.get("source") or "VNTR-native"
    return f"rs4072037: {gt}  (T-fraction {rs.get('alt_frac')}, depth {rs.get('depth')}) — {src}"


try:                                  # script-safe (this module is also run directly)
    from .config import STRAND_SKEW_WARN
except ImportError:
    from config import STRAND_SKEW_WARN


def _dupc_line(dupc) -> str | None:
    """One-line RELIABLE dupC/frameshift verdict from `dupc_dispatch` (the depth-routed caller: per-allele
    `pcr_dupc` on amplicons, statistical/positional on WGS/AS). None if not run. This is the TRUSTWORTHY dupC —
    distinct from the consensus-nomenclature `X-59dupC` token, which is depth-fragile on a smear."""
    if not dupc:
        return None
    res = dupc.get("result", dupc) or {}
    regime, caller = dupc.get("regime", "?"), dupc.get("caller", "?")
    status = res.get("status")                       # pcr_dupc / pcr_variant → positive/negative/undetermined
    if status is not None:
        allele = res.get("mut_allele")
        crit = res.get("criterion", "")
        return (f"59dupC: {str(status).upper()} — {caller} [{regime}]"
                f"{f', mutant allele={allele}' if allele else ''}{f' · {crit}' if crit else ''}")
    best = res.get("best") or {}
    frac, idx = best.get("frac"), best.get("index")
    if res.get("called"):
        return f"59dupC: CALLED — {caller} [{regime}] (tier {res.get('tier','')}, index {idx}, frac {frac})"
    tail = f" (best frac {frac} @ index {idx})" if frac is not None else ""
    return f"59dupC: not called — {caller} [{regime}]{tail}"


def _scaffold_block(sc) -> str:
    """Per-allele nomenclature + frameshift from the length-driven scaffold (mappy) — the smear-robust
    consensus that replaces the generic 'X' the flat multi-contig ranking produces. Empty if not run."""
    if not sc or not sc.get("per_allele"):
        return ""
    lines = ["=" * 68,
             "PER-ALLELE CONSENSUS — length-driven scaffold (smear-robust, replaces the generic 'X')"]
    for r in sc["per_allele"]:
        if not r.get("available"):
            lines.append(f"  allele ~{r.get('target')} cp : {r.get('note', 'no spanning reads')}")
            continue
        lines.append(f"  allele ~{r.get('target')} cp  (binned {r.get('n_binned')}, aligned "
                     f"{r.get('n_aligned')}):  frameshift = {r.get('indel') or 'none decoded'}")
        if r.get("nomenclature"):
            lines.append(f"      {r['nomenclature']}")
    lines.append("=" * 68)
    return "\n".join(lines) + "\n\n"


def _arbiter_length_block(arb, rs=None, dupc=None) -> str:
    """A prominent, FIRST-in-report block for the alignment-free arbiter length + the VNTR-native rs4072037
    genotype + the depth-routed RELIABLE dupC verdict — the numbers the user should read (the contig ranking
    below is unreliable on a smear). Empty string when the length is unavailable."""
    if not arb or not arb.get("available"):
        return ""
    al, ct = arb.get("alleles", []), arb.get("counts", [])
    lc = "  [LONG allele low-confidence / PCR-depleted]" if arb.get("long_low_conf") else ""
    if len(al) == 2:
        body = f"HET: {al[0]} / {al[1]} copies  (n={ct[0]}/{ct[1]}){lc}"
    elif len(al) == 1:
        body = f"HOM: ~{al[0]} copies  (n={ct[0] if ct else '?'})"
    else:
        body = "no clear peak — too shallow to length"
    plat = arb.get("platform", "unknown")
    if arb.get("is_amplicon"):
        prod = arb.get("product_bp", [])
        bands = " / ".join(f"~{round(p / 1000, 1)}kb" for p in prod) if prod else "?"
        ctx = f"  Technology: {plat}  ·  PCR amplicon — product {bands} (flanks+VNTR)"
    else:
        ctx = f"  Technology: {plat}  ·  adaptive sampling / WGS (not PCR)"
    # STRAND BALANCE: surface a skew at SAMPLE level. A variant carried by the weak strand is under-powered
    # even at huge depth, and our variant gates demand BOTH strands — so a skewed library silently costs
    # sensitivity. Report it here rather than letting it surface (or not) inside the variant caller.
    sb = ""
    _mf = arb.get("strand_minor_frac")
    if _mf is not None and arb.get("n_fwd") is not None:
        _tag = ("  ⚠ STRAND-SKEWED — a variant on the weak strand is under-powered; bi-strand gates may "
                "fail despite depth" if _mf < STRAND_SKEW_WARN else "")
        sb = f"  Strand balance: {arb['n_fwd']} fwd / {arb['n_rev']} rev (minor {_mf:.0%}){_tag}\n"
    rs_line = _rs4072037_line(rs)
    rs_block = f"  {rs_line}\n" if rs_line else ""
    dupc_line = _dupc_line(dupc)
    dupc_block = f"  {dupc_line}\n" if dupc_line else ""
    bar = "=" * 68
    return (f"{bar}\nVNTR LENGTH — alignment-free arbiter ({arb['method']}) — THE RELIABLE LENGTH\n"
            f"  {body}\n{ctx}\n{sb}{rs_block}{dupc_block}{bar}\n\n")


def main(argv=None):
    args = parse_args(argv)
    # lazy + script-safe: this module is also executed directly (no package context) by the harness
    try:
        from .config import FLANK_TO_PHYSICAL_OFFSET
    except ImportError:
        from config import FLANK_TO_PHYSICAL_OFFSET
    if args.pcr and args.max_depth is None:
        args.max_depth = 200                     # deep PCR: default cap (fast + stable consensus)
    if args.pcr:
        print("[PCR] amplicon mode: --max-depth=%d ; ⚠ PCR slippage can bias VNTR length "
              "(±1 unit) — the chr1 alignment must use the FROZEN PCR params "
              "(set B: -z 600,200 -r 2000,20000 ; cf. config.PCR_MINIMAP_FLAGS)." % args.max_depth,
              file=sys.stderr)

    # ── 1. Count reads per contig ──────────────────────────────
    print("[STEP 1] Counting reads per reference contig…", file=sys.stderr)
    read_counts = count_reads_per_contig(args.bam, min_mq=args.min_mq)

    if not read_counts:
        print("[ERROR] No mapped reads found in the BAM.", file=sys.stderr)
        sys.exit(1)

    # Descending sort by read count
    ranked = sorted(read_counts.items(), key=lambda x: x[1], reverse=True)

    print(f"  Contigs ranked by read count:", file=sys.stderr)
    for rank_i, (contig, cnt) in enumerate(ranked[:max(args.top_n * 2, 10)]):
        marker = ' ◄ candidate haplotype' if rank_i < args.top_n else ''
        print(f"    {contig:<40} {cnt:>8} reads{marker}", file=sys.stderr)

    top_contigs = [c for c, _ in ranked[: args.top_n]]

    # ── LENGTH GUARDRAIL (cf. memory muc1-vntr-length-pitfall) ────────
    # Contig counting UNDER-DETECTS LONG alleles (few full-span reads) and
    # can give a FALSE-SHORT. Warn if top contig weakly supported OR flat distribution (smear).
    _top_cnt = ranked[0][1]
    _tail = ranked[min(9, len(ranked) - 1)][1]
    _flat = len(ranked) >= 5 and _top_cnt < 1.5 * _tail
    length_warning = None
    if _top_cnt < 20 or _flat:
        length_warning = (
            ("top contig supported by < 20 reads. " if _top_cnt < 20 else "")
            + ("FLAT read distribution across contigs = multi-contig smear (reads do not map uniquely to "
               "one length-contig, so the top-N ranking and the per-contig consensus are UNRELIABLE — the "
               "nomenclature will collapse to generic 'X' motifs). " if _flat else "")
            + "LONG alleles are under-detected here — DO NOT conclude the length from this report. Use the "
              "alignment-free arbiter `vntr_raw_length.py` for the length instead.")
        print("[⚠ LENGTH UNRELIABLE] " + length_warning, file=sys.stderr)

    # ── 2. Load the FASTA reference ───────────────────────────────────
    ref_fa = pysam.FastaFile(args.ref)

    # ── 3. For each candidate haplotype: consensus + matching ────────
    all_results = []
    # per-haplotype bi-strand read support (rank → {'n_reads','n_fwd','n_rev'}); only filled when the
    # genomic two-axis add-on runs, and used there to COVERAGE-GATE a plain-call frameshift onset.
    hap_coverage = {}

    # Detect the length-homozygous case
    same_length = is_same_length_case(read_counts, top_contigs, args.homo_threshold)
    if same_length:
        print(
            f"\n[INFO] Probable length-homozygous case detected: "
            f"contig '{top_contigs[0]}' concentrates ≥ {args.homo_threshold*100:.0f}% of the reads. "
            f"Intra-contig phasing enabled.",
            file=sys.stderr,
        )

    if same_length and args.top_n >= 2:
        # ── Homozygous case: intra-contig phasing on top-1 ──────────────
        dominant_contig = top_contigs[0]
        ref_len_dom = ref_fa.get_reference_length(dominant_contig)

        print(f"\n[STEP 2] Intra-contig phasing on {dominant_contig}…", file=sys.stderr)

        # Identify the phasing SNPs
        snp_pos = find_phasing_snps(
            bam_path  = args.bam,
            contig    = dominant_contig,
            ref_len   = ref_len_dom,
            min_bq    = args.min_bq,
            min_mq    = args.min_mq,
            min_depth = max(args.min_depth, 5),
            min_af    = args.snp_min_af,
            verbose   = args.verbose,
        )
        print(f"  {len(snp_pos)} phasing SNP(s) identified.", file=sys.stderr)

        if len(snp_pos) == 0:
            # No SNP → probably true homozygous (same sequences)
            print(
                "  No informative SNP — haplotypes probably identical in sequence. "
                "Reconstructing a single consensus.",
                file=sys.stderr,
            )
            # Build a single consensus and duplicate it
            consensus = build_consensus_for_contig(
                bam_path=args.bam, ref_fa=ref_fa, contig=dominant_contig,
                min_depth=args.min_depth, min_bq=args.min_bq, min_mq=args.min_mq,
                ins_threshold=args.ins_threshold, long_ins_threshold=args.long_ins_threshold,
                del_threshold=args.del_threshold, verbose=args.verbose,
            )
            hap_consensuses = [(dominant_contig, consensus, True),
                               (dominant_contig, consensus, True)]  # flag True = duplicated
            hap_readnames = [None, None]                            # None = whole-contig support
        else:
            # Phase the reads into two groups
            reads_A, reads_B = phase_reads(
                bam_path=args.bam, contig=dominant_contig, ref_len=ref_len_dom,
                snp_positions=snp_pos, min_bq=args.min_bq, min_mq=args.min_mq,
                verbose=args.verbose,
            )
            if not reads_A or not reads_B:
                print("  Phasing failed (empty groups). Falling back to global consensus.",
                      file=sys.stderr)
                consensus = build_consensus_for_contig(
                    bam_path=args.bam, ref_fa=ref_fa, contig=dominant_contig,
                    min_depth=args.min_depth, min_bq=args.min_bq, min_mq=args.min_mq,
                    ins_threshold=args.ins_threshold, long_ins_threshold=args.long_ins_threshold,
                    del_threshold=args.del_threshold, verbose=args.verbose,
                )
                hap_consensuses = [(dominant_contig, consensus, False),
                                   (dominant_contig, consensus, False)]
                hap_readnames = [None, None]                        # phasing failed → whole-contig support
            else:
                print(f"  Haplotype A: {len(reads_A)} reads — Haplotype B: {len(reads_B)} reads",
                      file=sys.stderr)
                cons_A = build_consensus_from_readset(
                    bam_path=args.bam, ref_fa=ref_fa, contig=dominant_contig,
                    read_names=set(reads_A), min_depth=args.min_depth,
                    min_bq=args.min_bq, min_mq=args.min_mq,
                    ins_threshold=args.ins_threshold, long_ins_threshold=args.long_ins_threshold,
                    del_threshold=args.del_threshold, verbose=args.verbose,
                )
                cons_B = build_consensus_from_readset(
                    bam_path=args.bam, ref_fa=ref_fa, contig=dominant_contig,
                    read_names=set(reads_B), min_depth=args.min_depth,
                    min_bq=args.min_bq, min_mq=args.min_mq,
                    ins_threshold=args.ins_threshold, long_ins_threshold=args.long_ins_threshold,
                    del_threshold=args.del_threshold, verbose=args.verbose,
                )
                hap_consensuses = [(dominant_contig, cons_A, False),
                                   (dominant_contig, cons_B, False)]
                hap_readnames = [set(reads_A), set(reads_B)]        # per-haplotype phased read groups

        # Build the results for the 2 phased haplotypes
        for rank, (contig, consensus, is_dup) in enumerate(hap_consensuses, 1):
            if args.genomic_bam:
                rn = hap_readnames[rank - 1]
                hap_coverage[rank] = (
                    _strand_counts_for_readset(args.bam, contig, rn, args.min_mq) if rn is not None
                    else _strand_counts_for_contig(args.bam, contig, args.min_mq))
            if args.dump_consensus:
                tag = args.sample or 'consensus'
                label = f'hap{rank}_phased{"_dup" if is_dup else ""}_{contig}'
                print(f">{tag}_{label}", file=sys.stderr)
                for i in range(0, len(consensus), 80):
                    print(consensus[i: i + 80], file=sys.stderr)

            matched, unmatched = match_motifs(
                consensus, max_mismatch=args.max_mismatch, verbose=args.verbose,
            )
            all_results.append({
                'rank':          rank,
                'contig':        contig,
                'contig_length': ref_fa.get_reference_length(contig),
                'read_count':    read_counts.get(contig, 0),
                'consensus':     consensus,
                'nomenclature':  haplotype_string(matched),
                'n_motifs':      len(matched),
                'n_exact':       sum(1 for m in matched if m['mismatches'] == 0),
                'n_approx':      sum(1 for m in matched if m['mismatches'] > 0),
                'motifs':        matched,
                'unmatched':     unmatched,
                'phased':        not is_dup,
                'same_length':   True,
            })
            dup_tag = ' [identical — probable true homozygous]' if is_dup else ''
            print(f"  → Haplotype {rank}{dup_tag}: {haplotype_string(matched)}", file=sys.stderr)

    else:
        # ── Standard case: haplotypes of different lengths ───────────
        # PER-CONTIG decomposition (independent) → parallelizable (--threads) + optional
        # subsampling (--max-depth / --pcr). The worker `_process_contig` opens its own pysam.
        tasks = [{
            'bam': args.bam, 'ref': args.ref, 'contig': contig, 'rank': rank,
            'read_count': read_counts[contig], 'max_depth': args.max_depth,
            'min_depth': args.min_depth, 'min_bq': args.min_bq, 'min_mq': args.min_mq,
            'ins_threshold': args.ins_threshold, 'long_ins_threshold': args.long_ins_threshold,
            'del_threshold': args.del_threshold, 'max_mismatch': args.max_mismatch,
        } for rank, contig in enumerate(top_contigs, 1)]

        nproc = min(args.threads, len(tasks))
        if nproc > 1:
            print(f"\n[STEP 2] Decomposing {len(tasks)} contigs across {nproc} processes…", file=sys.stderr)
            import multiprocessing as mp
            with mp.Pool(nproc) as pool:
                results = pool.map(_process_contig, tasks)
        else:
            results = [_process_contig(t) for t in tasks]

        for res in sorted(results, key=lambda r: r['rank']):
            all_results.append(res)
            if args.genomic_bam:
                hap_coverage[res['rank']] = _strand_counts_for_contig(
                    args.bam, res['contig'], args.min_mq)
            if args.dump_consensus:
                tag = args.sample or 'consensus'
                print(f">{tag}_hap{res['rank']}_{res['contig']}", file=sys.stderr)
                for i in range(0, len(res['consensus']), 80):
                    print(res['consensus'][i: i + 80], file=sys.stderr)
            print(f"  → Haplotype {res['rank']} ({res['contig']}, {res['read_count']} reads): "
                  f"{res['nomenclature']}", file=sys.stderr)

    ref_fa.close()

    # ── 3bis. Two-axis clinical summary (optional; needs the genomic chr1 BAM) ──
    two_axis = None
    if args.genomic_bam:
        from . import two_axis_report as TA
        offset = args.copy_offset if args.copy_offset is not None else FLANK_TO_PHYSICAL_OFFSET
        haps = [{'rank': h['rank'], 'length': h['n_motifs'], 'indel': _indel_summary(h['motifs']),
                 'coverage': hap_coverage.get(h['rank'])}
                for h in all_results]
        mutation = ({'variant': args.mut_variant, 'repeat': args.mut_repeat,
                     'mutant_length': args.mut_allele_len} if args.mut_variant else None)
        try:
            two_axis = TA.build_summary(args.genomic_bam, haps, chrom=args.snp_chrom,
                                        snp_pos=args.snp_pos, copy_offset=offset, ref=args.genome_ref,
                                        mutation=mutation, gate_min_reads=args.gate_min_reads,
                                        gate_min_strand=args.gate_min_strand,
                                        smear=length_warning is not None)
        except Exception as e:                       # never let the add-on sink the core call
            print(f"[WARN] two-axis summary failed: {e}", file=sys.stderr)
            two_axis = {"available": False, "note": f"error: {e}"}

    # ── 3ter. AUTHORITATIVE arbiter length (alignment-free) — END-TO-END, no separate command needed.
    # The contig ranking above is UNRELIABLE on a multi-contig smear; THIS is the length to report. Auto-
    # selects flank (native LR-PCR/AS/WGS) vs cassette 1→9 (foreign/short PCR). Never sinks the core call.
    arbiter_len = None
    try:
        import vntr_raw_length as V
        _off = args.copy_offset if args.copy_offset is not None else FLANK_TO_PHYSICAL_OFFSET
        arbiter_len = V.auto_length(args.bam, chrom=None, offset=_off)
    except Exception as e:
        print(f"[WARN] arbiter length failed: {e}", file=sys.stderr)

    # AUTO-SCAFFOLD on a multi-contig SMEAR: the top-N ranking + per-contig consensus above collapse to
    # generic 'X', but the alignment-free arbiter binned the alleles fine. Length-drive a DENSE per-allele
    # consensus (Ilias' span_bin idea, merged into vntr_scaffold via mappy) → recover the motifs + frameshift
    # the smear hides. Only fires on a smear (else the plain top-N call is already dense/authoritative).
    # Trigger whenever the multi-contig length is UNRELIABLE — FLAT smear OR a THIN top contig (< 20 reads):
    # both mean the per-contig consensus is starved (7/4 reads on the wrong contigs here), which the
    # length-driven scaffold densifies by pooling ALL of an allele's spanning reads onto one contig.
    scaffold = None
    if length_warning and arbiter_len and arbiter_len.get("alleles"):
        try:
            import tempfile
            from . import vntr_scaffold as SC
            scaffold = SC.scaffold_both_alleles(
                args.bam, args.ref, workdir=tempfile.mkdtemp(prefix="muc1_scaffold_"), offset=_off,
                alleles=arbiter_len["alleles"], pacbio=(arbiter_len.get("platform") == "PacBio HiFi"),
                sample=args.sample or "sample", threads=args.threads)
            print("[scaffold] smear → length-driven per-allele consensus:", file=sys.stderr)
            for r in (scaffold.get("per_allele") or []):
                print(f"[scaffold]   allele {r.get('target')}: {r.get('nomenclature', '')} → "
                      f"{r.get('indel') or '(none)'}", file=sys.stderr)
        except Exception as e:
            print(f"[WARN] auto-scaffold failed: {e}", file=sys.stderr)

    # rs4072037 read VNTR-NATIVELY off the caller's own multi-contig BAM (offset 4265, no GRCh38 / no
    # --genomic-bam needed) — so the SNP resolves end-to-end for AS/PCR alike, not "n=0" when the genomic
    # BAM is the VNTR one. (The two-axis `--genomic-bam` path stays available for phased/severity work.)
    rs4072037 = None
    try:
        from .detectors import splice_snp as SS
        rs4072037 = SS.genotype_rs4072037_vntr_ref(args.bam, args.ref, min_mq=args.min_mq)
    except Exception as e:
        print(f"[WARN] rs4072037 (VNTR-native) failed: {e}", file=sys.stderr)
    # Fallback (LM 2026-07-25): the VNTR-native genotype is starved when the multi-contig alignment SMEARS
    # (reads spread ~1/contig → no depth at the fixed 5' offset). If a genomic hg38 BAM/CRAM is on hand
    # (--genomic-bam, which `run` already forwards for an aligned input), read rs4072037 directly at the
    # GRCh38 locus chr1:155,192,276 — those reads DO cover it. Only kicks in when our own method returns
    # undetermined, so it never overrides a good VNTR-native call.
    if (rs4072037 is None or rs4072037.get("snp_genotype") is None) and args.genomic_bam:
        try:
            from .detectors import splice_snp as SS
            g = SS.genotype_snp(args.genomic_bam, args.genome_ref, chrom=args.snp_chrom,
                                pos1=args.snp_pos, min_mq=args.min_mq)
            if g.get("snp_genotype") is not None:
                g["source"] = "genomic-fallback"
                rs4072037 = g
                print(f"[rs4072037] VNTR-native undetermined (multi-contig smear) → genomic fallback "
                      f"@ {args.snp_chrom}:{args.snp_pos or 155192276} = {g.get('gt_label')}", file=sys.stderr)
        except Exception as e:
            print(f"[WARN] rs4072037 genomic fallback failed: {e}", file=sys.stderr)

    # reliable dupC verdict (depth-routed caller), wired in by `run` as a JSON — distinct from the
    # consensus-nomenclature X-59dupC token (which is depth-fragile on a smear).
    dupc_verdict = None
    if args.dupc_json:
        try:
            dupc_verdict = json.load(open(args.dupc_json))
        except Exception as e:
            print(f"[WARN] --dupc-json unreadable: {e}", file=sys.stderr)

    report_verdict = None
    if args.report_json:
        try:
            report_verdict = json.load(open(args.report_json))
        except Exception as e:
            print(f"[WARN] --report-json unreadable: {e}", file=sys.stderr)

    # CONGRUENCE across the independent layers — a disagreement about the pathogenic variant must be SAID,
    # not left for the reader to notice (on a real case three layers looked contradictory and only a manual
    # read-fraction computation settled it).
    cong = None
    try:
        from . import congruence as CG
        _hv, _hr = None, None
        for _h in all_results:                      # the caller's own nomenclature call, if it decoded one
            _ind = _indel_summary(_h.get('motifs', []))
            if _ind:
                _hv = _ind.split(' ~ ')[0]
                _m = re.search(r'repeat\s*(\d+)', _ind)
                _hr = int(_m.group(1)) if _m else None
                break
        cong = CG.check(dupc=dupc_verdict, report=report_verdict, arbiter=arbiter_len,
                        hap_variant=_hv, hap_repeat=_hr)
    except Exception as e:
        print(f"[WARN] congruence check failed: {e}", file=sys.stderr)

    # ── 4. Text report ─────────────────────────────────────────────────
    report = format_report(all_results, sample=args.sample, length_warning=length_warning)
    report = (_arbiter_length_block(arbiter_len, rs4072037, dupc_verdict)
              + (CG.render(cong) if cong else "") + _scaffold_block(scaffold) + report)
    if two_axis is not None:
        from . import two_axis_report as TA
        report += "\n" + TA.render_text(two_axis)

    if args.output == '-':
        print(report)
    else:
        with open(args.output, 'w') as fh:
            fh.write(report + '\n')
        print(f"\n[INFO] Report written to: {args.output}", file=sys.stderr)

    # ── 5. Consensus FASTA export ────────────────────────────────────────
    if args.fasta_consensus:
        write_fasta_consensus(all_results, args.fasta_consensus, sample=args.sample)
        print(f"[INFO] Consensus FASTA written to: {args.fasta_consensus}", file=sys.stderr)

    # ── 6. PDF export ────────────────────────────────────────────────────
    if args.pdf:
        write_pdf_report(all_results, args.pdf, sample=args.sample, two_axis=two_axis,
                         length_warning=length_warning, arbiter=arbiter_len, rs4072037=rs4072037,
                         dupc=dupc_verdict, congruence=cong)
        print(f"[INFO] PDF report written to: {args.pdf}", file=sys.stderr)

    # ── 7. JSON export ───────────────────────────────────────────────────
    if args.json:
        payload = {'sample': args.sample, 'haplotypes': list(all_results)}
        if arbiter_len is not None:
            payload['arbiter_length'] = arbiter_len
        if rs4072037 is not None:
            payload['rs4072037'] = rs4072037
        if dupc_verdict is not None:
            payload['dupc'] = dupc_verdict
        if scaffold is not None:
            payload['scaffold'] = scaffold
        if cong is not None:
            payload['congruence'] = cong
        if report_verdict is not None:
            payload['pcr_report'] = report_verdict
        if two_axis is not None:
            payload['two_axis'] = two_axis
        with open(args.json, 'w') as jf:
            json.dump(payload, jf, indent=2)
        print(f"[INFO] JSON export: {args.json}", file=sys.stderr)

    # ── 8. Export candidate novel motifs (60 bp, approximate matches) ───
    if args.novel_motifs:
        candidates = collect_novel_60bp_motifs(all_results)
        write_novel_motifs(candidates, args.novel_motifs, sample=args.sample)
        print(
            f"[INFO] Candidate novel motifs ({len(candidates)} unique sequences) "
            f"→ {args.novel_motifs}",
            file=sys.stderr,
        )

    print("\n[DONE]", file=sys.stderr)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

