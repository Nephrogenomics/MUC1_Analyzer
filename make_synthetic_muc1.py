#!/usr/bin/env python3
"""Génère un jeu de test SYNTHÉTIQUE pour MUC1_Analyzer (aucune donnée patient).

Crée :
  - vntr_ref.fa (+ .fai) : contigs `MUC1_VNTR_Nrepeats`, chacun = N x motif A (60 pb)
  - synthetic.bam (+ .bai) : reads mappés, hétérozygote en longueur (20 vs 25 répétitions)

Le but est un smoke-test : prouver que le pipeline tourne bout-en-bout localement.
"""
import pysam, os, sys


def main(out=None):
    """Write the synthetic reference and BAM into `out`. Called from __main__ — importing this
    module must NOT write anything: it ships as a user tool, and a release gate that merely
    imports it once committed a generated `vntr_ref.fa` into the public repository."""
    OUT = out if out is not None else (sys.argv[1] if len(sys.argv) > 1 else ".")
    os.makedirs(OUT, exist_ok=True)

    # Motif A canonique (60 pb) tiré de KNOWN_REPEATS dans muc1_analyzer/caller.py
    MOT_A = "GCCCACGGTGTCACCTCGGCCCCGGAGAGCAGGCCGGCCCCGGGCTCCACCGCGCCCGCA"
    assert len(MOT_A) == 60

    # Contigs : quelques longueurs de VNTR (dont les 2 vrais haplotypes 20 & 25)
    LENGTHS = [18, 20, 22, 25, 30]
    contigs = {f"MUC1_VNTR_{n}repeats": MOT_A * n for n in LENGTHS}

    # --- 1. Référence FASTA ---
    ref_path = os.path.join(OUT, "vntr_ref.fa")
    with open(ref_path, "w") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i+70] + "\n")
    pysam.faidx(ref_path)
    print(f"[OK] {ref_path}  ({len(contigs)} contigs)")

    # --- 2. BAM synthétique ---
    # En-tête : une ligne SQ par contig
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
             "SQ": [{"SN": name, "LN": len(seq)} for name, seq in contigs.items()]}

    # Nb de reads par contig : haplotypes 25 (12 reads) et 20 (10 reads) dominent,
    # quelques reads "bruit" sur les leurres pour un ranking réaliste.
    read_plan = {"MUC1_VNTR_25repeats": 12, "MUC1_VNTR_20repeats": 10,
                 "MUC1_VNTR_18repeats": 2, "MUC1_VNTR_22repeats": 2,
                 "MUC1_VNTR_30repeats": 1}

    name2tid = {name: i for i, name in enumerate(contigs)}
    bam_unsorted = os.path.join(OUT, "synthetic.unsorted.bam")
    with pysam.AlignmentFile(bam_unsorted, "wb", header=header) as bam:
        for contig, nreads in read_plan.items():
            seq = contigs[contig]
            L = len(seq)
            for k in range(nreads):
                a = pysam.AlignedSegment()
                a.query_name = f"read_{contig}_{k}"
                a.query_sequence = seq
                a.flag = 0                       # mappé, brin +, primaire
                a.reference_id = name2tid[contig]
                a.reference_start = 0
                a.mapping_quality = 60
                a.cigartuples = [(0, L)]         # L match (M)
                a.query_qualities = pysam.qualitystring_to_array("I" * L)  # Q40
                bam.write(a)

    sorted_bam = os.path.join(OUT, "synthetic.bam")
    pysam.sort("-o", sorted_bam, bam_unsorted)
    pysam.index(sorted_bam)
    os.remove(bam_unsorted)
    print(f"[OK] {sorted_bam}  ({sum(read_plan.values())} reads sur {len(read_plan)} contigs)")
    print(f"[OK] plan de reads : {read_plan}")


if __name__ == "__main__":
    main()
