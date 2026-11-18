# MUC1_Analyzer — tutorial

From a runnable synthetic quickstart (no patient data) to the real-data workflow and the
`MUC1_Score` (frameshift tail). Pair with the [README](../README.md).

External tools on `PATH`: **samtools**, **minimap2** (see `environment.yml`).

---

## 1. Quickstart on synthetic data (no patient material)

```bash
# (a) synthetic VNTR reference + reads
python make_synthetic_muc1.py /tmp/demo               # -> /tmp/demo/vntr_ref.fa , synthetic.bam

# (b) call the two VNTR haplotypes.  --min-mq 0 is REQUIRED: the length-contigs are near-identical,
#     so minimap2 assigns MAPQ 0. The tool selects the top-N contigs by mapping density.
python -m muc1_analyzer call -b /tmp/demo/synthetic.bam -r /tmp/demo/vntr_ref.fa -s DEMO --min-mq 0 \
    --json /tmp/demo/DEMO.json --pdf /tmp/demo/DEMO.pdf     # text + JSON + PDF report

# (c) MUC1_Score from the analyzer output
python -m muc1_analyzer score --analyzer-json /tmp/demo/DEMO.json --snp-genotype 1
```

The synthetic reads are pure motif A (no frameshift), so no mutation is called and no score is
produced (there is no frameshift tail) — that is correct, it just proves the pipeline runs end to end.

---

## 2. Real data — detect, then score

**Inputs**: FASTQ or unaligned BAM (uBAM), ONT or PacBio HiFi, from LR-PCR, adaptive sampling, or WGS.
An aligned BAM/CRAM is refused (v1.0.1): extract the MUC1 reads first (below).
The MUC1 coding VNTR collapses on GRCh38, so length/structure are resolved by re-aligning the reads to
the multi-contig VNTR reference (one contig = one length). Chromosome naming (`chr1` / `1` /
`CM000663.2` / `NC_000001.11`) and genome (hg38 / T2T) are auto-detected.

```bash
# from an aligned BAM/CRAM: extract the MUC1 window (GRCh38) as FASTQ, keeping the methylation tags
samtools view -b -F 0x904 -T GRCh38.fa input.cram chr1:155180000-155200000 \
    | samtools fastq -T MM,ML - > SAMPLE.MUC1.fq

# everything in one command: align + call + score
python -m muc1_analyzer run -i SAMPLE.MUC1.fq --input-type wgs-ont \
    -r MUC1_fakedVNTR1to150revcomplKirby.fa -s SAMPLE -o out/
#   --input-type as-ont for adaptive sampling, pcr-ont for an LR-PCR amplicon, *-pacbio for HiFi.
```

That is the recommended path. The steps are also available individually, which is what you want when
you already have an alignment, or when you want to cross-check a length:

```bash
# 1) align only (methylation tags MM/ML are preserved)
python -m muc1_analyzer prepare -i SAMPLE.MUC1.fq \
    -r MUC1_fakedVNTR1to150revcomplKirby.fa -o muc1.vntr.bam -s SAMPLE

# 2) call haplotypes + motifs + frameshift -> JSON
python -m muc1_analyzer call -b muc1.vntr.bam -r MUC1_fakedVNTR1to150revcomplKirby.fa -s SAMPLE \
    --min-mq 0 --json SAMPLE.json --pdf SAMPLE.pdf

# 3) VNTR length cross-check (alignment-free, robust to collapse)
python vntr_raw_length.py --bam muc1.vntr.bam --name SAMPLE

# 4) MUC1_Score (frameshift tail; rs4072037 is informational: --snp-genotype 0/1/2 or --rs4072037-mut C/T)
python -m muc1_analyzer score --analyzer-json SAMPLE.json --snp-genotype 1 -o SAMPLE.score.json
```

**On input types.** Only FASTQ and uBAM are accepted, and they are taken whole, so they should already
be MUC1-enriched (an LR-PCR amplicon, a targeted uBAM, or reads extracted at the locus as above).
Feeding a whole-run FASTQ aligns *everything* to the VNTR reference. Do not read the alleles from the
number of reads per contig of the VNTR BAM: a read that does not span the VNTR lands on an arbitrary
length-contig. The lengths come from the alignment-free arbiter.

Behaviour worth knowing (automatic):

- **LR-PCR amplicons (`--pcr`, auto-detected).** A VNTR-spanning amplicon read must chain across the
  whole tandem; the default minimap2 bandwidth fragments it and you lose length discrimination.
  `prepare` applies frozen chaining flags (`-z 600,200 -r 2000,20000`) — **auto-enabled** when the
  locus is amplicon-scale (≥ 2000 reads that actually *align* to the VNTR reference, so a busy
  adaptive-sampling window is not mistaken for an amplicon). Force with `--pcr`, disable auto-detection
  with `--no-pcr-autodetect`.

**On rs4072037.** `run` genotypes it directly on the VNTR-reference alignment — the 150 contigs share
an identical 5′ flank, so the SNP sits at a fixed offset. No GRCh38 alignment is needed for it. Note
this is the *genotype*; the score reports, for information, the base on the **mutant** allele (see §6).

**Manual / visual scoring** — for variants finalized by eye (e.g. `del8_27`, invisible to the automated
caller) or a coverage-floored heterozygote:

```bash
python -m muc1_analyzer score \
    --clinical-call "64 repeats | 80 repeats (del8_27 @ repeat 37)" \
    --rs4072037-mut C -s SAMPLE
#   -> fs_tail (80 − 37 = 43 repeats), percentile, category; rs4072037 informational; provenance tagged
```

The score is the frameshift tail placed on the calibration cohort's Gaussian (n = 33): percentile and
one of four categories (very mild < 10th, mild < 50th, severe ≤ 90th, very severe > 90th). rs4072037 is
reported for information and carries no weight.

---

## 3. Detecting the frameshift at low depth (dupC caller + PoN)

The statistical dupC caller compares each C-tract context to a **Panel-of-Normals** (ONT homopolymer
error per context) with a binomial test and a Bayesian posterior. A **calibrated PoN is shipped**
(`pon_results/pon_dupc.json`, 71 public 1000 Genomes LCLs, ONT R10.4.1 WGS):

```bash
python -m muc1_analyzer.dupc_caller -b muc1.hg38.bam --genome-ref GRCh38.fa \
    --pon pon_results/pon_dupc.json --mut-len 8              # 8 = dupC ; 5 = delCC ; --seq-anchor for T2T/urine
```

Check what you were shipped, and its provenance (`meta` block):

```bash
python verify_pon.py pon_results/pon_dupc.json
```

**When to build your own instead.** `f` is sensitive to chemistry and basecaller version (both set
homopolymer accuracy). The shipped panel mixes dorado versions, so it is conservative rather than
matched to a specific run. If your chemistry/basecaller differs materially — and always for the
anti-circularity requirement (never test a cohort against a PoN built from itself) — build one on
YOUR negative controls:

```bash
# per-sample profile (one JSON per control), then aggregate, then verify before sharing
python -m muc1_analyzer.detectors.dupc_pon --bam control_i.bam --genome-ref GRCh38.fa --out prof_i.json
python -m muc1_analyzer.detectors.dupc_pon --aggregate prof_*.json --out my_pon.json
python verify_pon.py my_pon.json
```

---

## 4. Visual review (IGV)

```bash
# per-haplotype BAM + reference + GFF3 track for eyeballing a variant in IGV
python -m muc1_analyzer.review_bundle --cram SAMPLE.cram --ref GRCh38.fa --vntr-ref MUC1_fakedVNTR1to150revcomplKirby.fa --out bundle/
```

Load the per-haplotype BAM with `MUC1_fakedVNTR1to150revcomplKirby.fa` + `.gff3` in IGV. The
`muc1_analyzer call --pdf` report is a standalone visual summary of the two haplotypes and their motifs.

---

## 4b. Reading VNTR length correctly at low coverage (the `LENGTH UNRELIABLE` warning)

The caller separates the two alleles by **read count per length-contig**: the top-N contigs by
supporting reads become the haplotypes. At **low coverage this split is not trustworthy** — reads
scatter across near-identical-length contigs, and a 1-copy difference between the top two is within
the noise. The caller says so explicitly:

```
[⚠ LENGTH UNRELIABLE] top contig supported by < 20 reads. LONG alleles are under-detected here —
DO NOT conclude the length. Cross-check with sv8533_length on the WHOLE TANDEM …
```

**Worked example — a quasi-homozygous long sample at ~37× WGS.** The contig ranking was:

```
MUC1_VNTR_80repeats  11 reads  ◄ candidate haplotype
MUC1_VNTR_79repeats   5 reads  ◄ candidate haplotype
MUC1_VNTR_81repeats   5 · 83:5 · 82:4 · 85:3 · 116:3 …
```

The caller reported haplotypes **80 and 79** — but this is **not** a real 79-vs-80 heterozygosity:
with 11 + 5 reads, the tallest contigs (79/80/81/82/83…) all carry a handful of reads, so the
1-copy difference is noise. The sample is **length-homozygous ~80**. Reading it as het would be a
mistake — and it is exactly why it is not a phasing/ASE candidate (no length handle to separate the
alleles).

**You don't have to do anything — `call` already reports the reliable length for you.** Since the
arbiter is now computed inside `call`, its length (and the rs4072037 genotype) sits at the TOP of the
report/PDF/JSON; the per-contig ranking below is the detail, not the number to quote. The banner reads:

```
VNTR LENGTH — flank AL→AH (alignment-free): 45 / 82 copies
Technology: PacBio HiFi  ·  PCR amplicon detected — product ~7.8 / ~10.7 kb (flanks + VNTR)
rs4072037: C/T  (T-fraction 0.49, depth 210) — VNTR-native
```

**The arbiter is amplicon-agnostic.** It auto-selects the method: the unique **flank** anchors (AL+AH)
for our native LR-PCR / adaptive-sampling / WGS, and the invariant **cassette** motifs `1→9` for a
*foreign or shorter* MUC1 PCR whose reads do not span our flanks — so length works on **any** protocol.
A PCR-depleted long allele is surfaced as `low-confidence`. To inspect it directly (or on a bare BAM):

```bash
python vntr_raw_length.py --bam muc1.vntr.bam --name SAMPLE                 # auto flank/cassette
python vntr_raw_length.py --bam muc1.vntr.bam --amplicon --name SAMPLE      # amplicon signature (kind, flanks, product kb, platform)
```

A genuine length-heterozygote shows a **clear bimodal peak** in the raw histogram (well-separated modes
with real read support), not two adjacent contigs a single copy apart.

## 5. Regenerating the VNTR reference from T2T-CHM13

The shipped `MUC1_fakedVNTR1to150revcomplKirby.fa` (+ `.gff3`) is generated from **T2T-CHM13v2.0**
(NCBI `GCA_009914755.4`). To rebuild it (e.g. a different unit or length range):

```bash
python -m muc1_analyzer.build_vntr_ref --genome chm13v2.0.fa --region <MUC1 array coords> --revcomp \
    --nmin 1 --nmax 150
```

---

## 6. Phasing rs4072037 with the mutant allele (informational)

The rs4072037 base reported with the score is the base **of the mutant allele** — not the genotype. It
does not change the score. A homozygote is unambiguous, a heterozygote needs the phase; without it the
output falls back to the cis rule and records that it did.

```bash
# single-molecule phasing: group full-length reads by their rs4072037 base, read each
# allele's VNTR length off the same molecules, then take the base of the mutant length
python phase_fs_snp.py --bam SAMPLE.chr1.bam --mut-len 44 --name SAMPLE
```

Why one read is enough: rs4072037 lies ~370 bp 5′ of the tandem, so a full-length read spans the SNP,
the whole array, and the frameshift. Requirements: reads that span SNP → VNTR → frameshift, and enough
depth **per allele** (LR-PCR qualifies; adaptive sampling often does, but a long allele sampled thin
will not resolve — the tool reports `resolved: False` rather than guessing).

Precedence for the reported base, highest first:

| Source | How | `splice_source` |
|---|---|---|
| Observed on the mutant allele | `--rs4072037-mut C\|T` (visual review, dRNA-seq) | `observed` |
| Phased from single molecules | `phase_fs_snp.py` → pass the base via `--rs4072037-mut` | `observed` |
| Inferred from genotype + cis rule | `--snp-genotype 0\|1\|2` (T is cis with the short VNTR, ~94%) | `inferred(dosage/cis)` |

A length-homozygous heterozygote has no cis handle: the base is left **unset** rather than guessed.
Check `splice_source` in the output before treating it as observed.

## Command reference

| Task | Command |
|---|---|
| Call (length + rs4072037 + amplicon, end-to-end) | `python -m muc1_analyzer call -b B.bam -r vntr_ref.fa --min-mq 0 --json out.json` |
| VNTR length only (auto flank/cassette) | `python vntr_raw_length.py --bam B.bam --name S` |
| Amplicon signature (kind, flanks, product kb, platform) | `python vntr_raw_length.py --bam B.bam --amplicon --name S` |
| dupC statistical call | `python -m muc1_analyzer.dupc_caller -b B.bam --genome-ref g.fa --pon PON.json --mut-len 8` |
| MUC1_Score (auto) | `python -m muc1_analyzer score --analyzer-json S.json --snp-genotype 1` |
| MUC1_Score (visual) | `python -m muc1_analyzer score --clinical-call "…(variant @ repeat N)" --rs4072037-mut C` |
| Phase rs4072037 (single-molecule) | `python phase_fs_snp.py --bam S.chr1.bam --mut-len N --name S` |
| Full pipeline | `python -m muc1_analyzer run -i input -r vntr_ref.fa -s S -o out/` |
| Align only | `python -m muc1_analyzer prepare -i S.MUC1.fq -r vntr_ref.fa -o S.vntr.bam` |
| Align, LR-PCR amplicon (usually auto) | `python -m muc1_analyzer prepare -i amplicon.fq.gz -r vntr_ref.fa -o S.vntr.bam --pcr` |
| Extract MUC1 reads from an aligned CRAM | `samtools view -b -F 0x904 -T g.fa S.cram chr1:155180000-155200000 \| samtools fastq -T MM,ML - > S.MUC1.fq` |
| IGV bundle | `python -m muc1_analyzer.review_bundle --cram S.cram --ref g.fa --vntr-ref vntr_ref.fa --out bundle/` |
