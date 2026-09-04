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
    # This is the FIRST command in the public quickstart, so a mistyped invocation is a new user's first
    # contact with the tool. Without this, `--help` was taken as an output directory and the run died in
    # a samtools usage traceback (audit, 2026-08-07).
    argv = sys.argv[1:]
    if out is None and argv and argv[0] in ("-h", "--help"):
        print(f"usage: python3 {os.path.basename(__file__)} [OUTDIR]\n\n"
              "Writes a SYNTHETIC MUC1 test set (no patient data) into OUTDIR (default: .):\n"
              "  vntr_ref.fa (+ .fai)   contigs MUC1_VNTR_Nrepeats, each N x the 60 bp motif A\n"
              "  synthetic.bam (+ .bai) mapped reads, length-heterozygous (20 vs 25 repeats)\n\n"
              "Then:\n"
              "  python3 -m muc1_analyzer run -i OUTDIR/synthetic.bam -r OUTDIR/vntr_ref.fa "
              "-s DEMO -o OUTDIR/out")
        return
    OUT = out if out is not None else (argv[0] if argv else ".")
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
