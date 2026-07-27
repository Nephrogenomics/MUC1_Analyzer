# MUC1_Analyzer

Long-read resolution of the **MUC1** coding VNTR for ADTKD-*MUC1*, from Oxford Nanopore or PacBio
reads (LR-PCR, adaptive sampling, or WGS). It calls the two VNTR haplotypes, decodes the repeat-motif
structure, detects the pathogenic frameshift (canonical **59dupC** and non-canonical variants), and
computes a **two-axis `MUC1_Score`** for genotype–phenotype work.

> Inputs: uBAM / BAM / CRAM / FASTQ · ONT or PacBio HiFi · PCR, adaptive-sampling, or WGS depth.
> The coding VNTR is GC-rich and collapses on GRCh38, so length/structure are resolved on a
> **synthetic multi-contig reference** (one contig = one VNTR length) built from **T2T-CHM13**.

## The two-axis `MUC1_Score`

The score is deliberately **two orthogonal axes**, reported separately (there is no single composite —
each axis maps to a different clinical dimension):

- **Axis A — onset** (*when*): `onset_index = 1 − position/L_mut`, the fraction of the mutant allele
  translated downstream of the frameshift as MUC1fs neoprotein (an earlier frameshift → longer
  neoprotein tail → earlier onset).
- **Axis B — severity** (*how severe*): the splice haplotype **rs4072037** of the mutant allele —
  C/MUC1-TR retains the VNTR (severe), T/MUC1-Y splices it out (protected).

plus the descriptor **VNTR ratio** = mutant / normal allele length (which allele is short, which is long,
which carries the mutation). Weights are direction-of-effect and **not calibrated** (small cohort).

## Install

```bash
conda env create -f environment.yml && conda activate muc1        # pins python + pysam + samtools + minimap2
# or: pip install -r requirements.txt   (and provide samtools + minimap2 on PATH)
python -m muc1_analyzer --help                                     # sanity: list the subcommands
```

External tools on `PATH`: **samtools** and **minimap2**. Optional: `reportlab` for the PDF report.

## Quickstart

One command takes any input to a score — it extracts the MUC1 reads, aligns them to the VNTR
reference, calls the haplotypes and scores the result:

```bash
python -m muc1_analyzer run -i reads.fastq.gz -r MUC1_fakedVNTR1to150revcomplKirby.fa -s SAMPLE -o out/
```

`-i` accepts **FASTQ, uBAM, or an aligned BAM/CRAM** (GRCh38 or T2T) — the input type and reference
are detected automatically, and an input already aligned to the VNTR reference is reused as-is. The
**sequencing platform is auto-detected from the read names** (PacBio `…/ccs` → `minimap2 map-hifi`,
ONT UUID → `map-ont`), so you rarely need `--pacbio`/`--ont`; pass one to override. `--ref genome.fa`
is needed only to read a CRAM.

`prepare` adapts to the input **automatically** — you rarely need these flags, but they exist:

| Situation | What `prepare` does (flag) |
|---|---|
| **LR-PCR amplicon** | Applies frozen minimap2 chaining flags (`-z 600,200 -r 2000,20000`) so a VNTR-spanning amplicon read chains across the tandem instead of fragmenting. **Auto-enabled** when the locus is amplicon-scale (≥ 2000 reads that actually *align* to the VNTR reference); force with `--pcr`, disable auto-detection with `--no-pcr-autodetect`. |
| **hg38-aligned BAM/CRAM** | **`--recover-unmapped` is ON by default**: GRCh38 collapses the coding VNTR, so VNTR-heavy reads are flagged *unmapped* and a plain window fetch misses them; they are pulled back and re-aligned (validated up to ~35×, e.g. 264 → 9324 reads). Disable with `--no-recover-unmapped`. |
| **T2T-CHM13 input** (full-genome or CRAM) | Restricted to the MUC1 VNTR window like an hg38 input (no whole-run over-selection). T2T represents the VNTR faithfully, so recover-unmapped is low-yield there (harmless). |

Outputs in `out/`: `SAMPLE.vntr.bam` (reads on the VNTR reference), `SAMPLE.analyzer.json`
(haplotypes, motif nomenclature, frameshift), `SAMPLE.score.json` (the two-axis score).

### What `call` reports — end to end, no extra command

`call` (and therefore `run`) puts the numbers you actually need **first** in the report, on the PDF and
in the JSON — computed from the VNTR-reference BAM alone (no `--genomic-bam`, no separate length tool):

- **VNTR length** — the alignment-free **arbiter** length (auto-selects the flank anchors for our native
  LR-PCR / AS / WGS, or the amplicon-agnostic **cassette** motifs `1→9` for a *foreign/shorter* MUC1 PCR
  whose reads don't span our flanks). This is the length to trust — the per-contig ranking below it is
  unreliable on a multi-contig smear. A PCR-depleted long allele is flagged *low-confidence*.
- **rs4072037** — genotyped VNTR-natively off the same BAM (C/C · C/T · T/T, with T-fraction and depth).
- **Technology** (ONT / PacBio) and, when the reads look like a PCR amplicon (fixed primer flanks), the
  **product band size in kb** (flanks + VNTR) — your gel band, and a dropout tell if it comes out short.

`vntr_raw_length.py --amplicon --bam SAMPLE.vntr.bam` prints that amplicon signature on its own; the tool
is amplicon-agnostic, so it works on **any** MUC1 PCR protocol, not only ours.

### Try it without any data

```bash
python make_synthetic_muc1.py /tmp/demo                          # synthetic reference + reads, no PHI
python -m muc1_analyzer call -b /tmp/demo/synthetic.bam -r /tmp/demo/vntr_ref.fa -s DEMO \
    --min-mq 0 --json /tmp/demo/DEMO.json
python -m muc1_analyzer score --analyzer-json /tmp/demo/DEMO.json --snp-genotype 1
```

### The subcommands

| Subcommand | What it does |
|---|---|
| `run` | the whole pipeline: `prepare` → `call` → `score` |
| `prepare` | any input → BAM on the VNTR reference (extract + align; methylation tags preserved) |
| `call` | VNTR haplotype caller: the reliable **arbiter length** + **rs4072037** (VNTR-native) + technology/amplicon up front, then motif nomenclature + frameshift |
| `score` | the two-axis `MUC1_Score` |
| `bundle` | per-haplotype BAM + reference + GFF3 for IGV visual review |

`python -m muc1_analyzer <subcommand> --help` for the options of each.

## Scoring a carrier

Two paths to the two-axis score:

```bash
# auto — from a MUC1_Analyzer JSON (frameshift + lengths detected)
python -m muc1_analyzer score --analyzer-json PATIENT.json --snp-genotype 1

# manual/visual — from a finalized clinical call (variants confirmed by eye, e.g. del8_27)
python -m muc1_analyzer score \
    --clinical-call "64 repeats | 80 repeats (del8_27 @ repeat 37)" \
    --rs4072037-mut C
#   -> onset_index, severity (splice), ratio; provenance is tagged (auto vs clinical/visual vs dRNA-seq)
```

`--clinical-call` parses the lab nomenclature: the parenthetical `(variant @ repeat N)` marks the mutant
allele (`@` canonical, `~` accepted). rs4072037 comes from the mutant-allele base (`--rs4072037-mut`,
observed visually or by dRNA-seq — highest precedence) or from a dosage (`--snp-genotype`) + the cis rule.

## References (genomes)

Not shipped (too large) — reference externally, or regenerate the VNTR reference:

- **T2T-CHM13v2.0** (VNTR reference source): NCBI `GCA_009914755.4`.
- **GRCh38** — *not required*. rs4072037 (`chr1:155,192,276`, C/T) is genotyped directly on the VNTR
  reference: the 150 contigs share an identical 5′ flank, so the SNP sits at a fixed offset and is read
  off the same alignment. A GRCh38 build is only needed to read a CRAM input (`--ref`).
- Rebuild the shipped VNTR reference: `python -m muc1_analyzer.build_vntr_ref --genome chm13v2.0.fa --region <MUC1 array>`.
  The pre-built reference (`MUC1_fakedVNTR1to150revcomplKirby.fa` + `.gff3` for IGV) is included.

## Detecting the frameshift (dupC) at low depth — Panel-of-Normals

The statistical dupC caller models the ONT homopolymer error per sequence context against a
**Panel-of-Normals** (PoN) — a `{sequence context → C-tract length histogram}` table. This is an
aggregate error model (the analogue of a technical error profile or gnomAD allele frequencies): it
contains no reads, no genotypes, and no sample identifiers.

**A calibrated PoN is shipped** as `pon_results/pon_dupc.json`, built on **71 public 1000 Genomes
LCLs** sequenced **ONT R10.4.1 (LSK114), WGS**. Its measured null over the common `7-7-7` context is
**f ≈ 0.0099** for dupC (8C) and **f ≈ 0.0038** for delCC (5C); pooled over all 314 contexts,
≈ 0.017 (8C) and ≈ 0.026 (5C). You can re-check exactly what you were given with the bundled
`verify_pon.py`, and read the provenance in the file's `meta` block.

> **Chemistry / basecaller caveat.** `f` is robust to coverage and protocol but **sensitive to
> chemistry and basecaller version** — both set homopolymer accuracy, which is precisely what this
> model measures. The shipped panel was basecalled with **mixed dorado versions**, so its `f` is a
> *mixture* over basecalling regimes: conservative (it includes older, noisier basecallers) rather
> than matched to any single one. If your chemistry or basecaller differs materially, **build your
> own PoN on your negative controls** (this is also the anti-circularity requirement — never test a
> cohort against a PoN built from it).

A small **synthetic demo PoN** (`pon_results/pon_dupc_demo.json`, no data) is also included as a
format example. To build your own: per-sample profiles → `muc1_analyzer.detectors.dupc_pon`
aggregation (see the tutorial), then `python verify_pon.py your_pon.json` before sharing it.

## Phasing rs4072037 with the mutant allele (preferred, not required)

The severity axis needs rs4072037 **on the mutant allele**, not just the genotype. Three ways, in
decreasing order of directness:

1. **Observed directly** — `--rs4072037-mut C|T` when the base has been read off the mutant allele
   (visual review, dRNA-seq). Highest precedence.
2. **Phased from single molecules** — `phase_fs_snp.py` groups full-length reads by their rs4072037
   base and reads each allele's VNTR length off the same molecules. rs4072037 sits ~370 bp 5′ of the
   tandem, so one read can span the SNP *and* the whole array; the mutant allele is then the length
   carrying the frameshift, and its base is **observed** rather than inferred.

   ```bash
   python phase_fs_snp.py --bam SAMPLE.chr1.bam --mut-len 44 --name SAMPLE
   ```

   This needs reads that span SNP → VNTR → frameshift (LR-PCR does; adaptive sampling often does,
   given enough per-allele depth).
3. **Inferred** — a genotype alone (`--snp-genotype 0|1|2`) plus the cis rule (rs4072037-T is cis
   with the short VNTR in ~94% of alleles). Homozygotes are unambiguous; a heterozygote whose two
   alleles are the same length has no cis handle and is left unresolved rather than guessed.

The score records which of the three was used (`splice_source`), so an inferred severity is never
silently reported as an observed one.

## Components

| File | Role |
|---|---|
| `muc1_analyzer/caller.py` (`… call`) | VNTR haplotype caller: per-contig read counts (= length), 2-haplotype consensus, motif nomenclature, indel/frameshift detection. Text / JSON / FASTA / PDF output. |
| `vntr_raw_length.py` | alignment-free VNTR length arbiter (flank-anchored on the read sequence; robust to tandem collapse). |
| `Bam_cleaner.py` | upstream span/soft-clip read filter (handles CRAM). |
| `muc1_analyzer/clinical_call`, `haplotype_score` | the **two-axis** score (onset ← position, severity ← rs4072037). |
| `muc1_analyzer/dupc_caller` + `dupc_dispatch` + `dupc_pon`/`dupc_prior`/`dupc_power` + `pcr_dupc`/`pcr_variant` | statistical frameshift caller (PoN null → binomial test → Bayesian posterior), routed by locus depth (PCR vs AS/WGS). |
| `muc1_analyzer/detectors/{vntr, splice_snp, vntr_dupc, vntr_scaffold, vntr_segment}` | VNTR length/motifs; rs4072037 genotyping; targeted dupC/segment probes. |
| `muc1_analyzer/build_vntr_ref` | (re)generate the multi-contig VNTR reference + GFF3 from a genome. |
| `muc1_analyzer/prepare` (`… prepare`) | universal input → VNTR-reference BAM: auto-detects FASTQ/uBAM/hg38/T2T, preserves MM/ML methylation tags, bypasses an input already aligned. |
| `verify_pon.py` | check a PoN is aggregate-only (no per-sample data) before sharing it. |
| `phase_fs_snp.py` | single-molecule phasing of rs4072037 with the VNTR allele (severity axis, observed). |
| `muc1_analyzer/review_bundle` | per-haplotype BAM + reference + GFF3 for **IGV visual review** of a variant. |

## Method note (MAPQ 0)

Length is resolved by re-aligning each read to the multi-contig reference (one contig per VNTR length,
generated from T2T-CHM13). Because the contigs are near-identical, minimap2 assigns **MAPQ 0**, so all
steps use `--min-mq 0` plus a top-*N* selection by mapping density. This is expected, not a defect.

## Code / Data availability · License · Citation

Patient sequencing/clinical data are potentially identifying and are **not** included; they are available
from the corresponding author on reasonable request, subject to ethics approval. Aggregate results and the
synthetic VNTR reference are in this repository.

Licensed under the **GNU Affero General Public License v3.0** (`LICENSE`). If you use this software,
please cite the MUC1_Analyzer manuscript (details on acceptance).
