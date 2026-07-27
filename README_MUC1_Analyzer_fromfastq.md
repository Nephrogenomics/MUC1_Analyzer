# MUC1_Analyzer_fromfastq.py

*Alignment-first MUC1 VNTR caller that runs **directly from a FASTQ** — length of both alleles (including two same-length or one-repeat-apart alleles, separated by SNP phasing), homopolymer **59dupC / frameshift** detection, **native rs4072037 severity**, and the two-axis **MUC1_Score** with cohort positioning, in a single self-contained pass.*

This tool is a **companion to `MUC1_Analyzer`** (the `muc1_analyzer` package). It reuses that package's validated internals (consensus reconstruction, motif nomenclature, frameshift detection, novel-motif discovery, SNP phasing, rs4072037 genotyping, two-axis clinical score) and adds everything needed to go **from raw adaptive-sampling long reads to a both-allele clinical result without a pre-built, locus-aligned BAM**.

> Runs from the repository root (it imports `muc1_analyzer` and `vntr_raw_length`). It lives on the branch **`muc1-analyzer-fromfastq`** (based on `Dev-Muc1-score`, which carries these dependencies).

---

## 1. Rationale — why a FASTQ entry point

`MUC1_Analyzer`'s `call` subcommand expects a BAM already aligned to the *faked* multi-length reference (`MUC1_VNTR_Nrepeats` contigs). That is ideal in a lab pipeline but makes three things hard when you only have reads:

1. **Length calling from adaptive-sampling (AS) reads.** On a near-identical tandem, aligning reads to a *multi-contig* reference smears them across neighbouring lengths (MAPQ-0 / split alignments), so a depth-per-contig peak is fragile and long alleles are under-detected.
2. **The homopolymer 59dupC.** ONT miscounts the C-tract, so a consensus/insertion-fraction call at a single position is not simultaneously sensitive and specific at AS depth.
3. **Reproducibility from a single input.** A full both-allele result required chaining several external tools.

`MUC1_Analyzer_fromfastq.py` takes a **FASTQ (or uBAM/BAM/CRAM — sequences only)** and does the whole thing internally: it measures length **alignment-first against a "ruler"**, fixes the physical length with an **indel-burden scaffold sweep**, recovers reads reliably, calls each allele's consensus/motifs/frameshift, detects the dupC, reads **rs4072037 natively** off the mutant allele, computes the **MUC1_Score**, and emits **merged, both-allele outputs**.

---

## 2. What it does (pipeline)

```
FASTQ / uBAM / BAM / CRAM  (sequences only)
        |
        v
[1] STREAMING length - read one at a time, align each to the RULER
        (MUC1_VNTR_150repeats); classify SPAN (both flanks -> exact length),
        LEFT/RIGHT (one flank -> lower bound), INT (internal). Only reads that
        cover the VNTR are kept in memory -> 30 M+ reads without OOM.
        |
        v
[2] Allele calling - cluster SPAN copies -> HET / HOM
        + scatter-fold: a parasitic second length-allele (measurement
          scatter around the true peak) is folded into the dominant one
          (distance-dependent relative threshold)
        |
        v
[3] Scaffold sweep - pick, per allele, the faked contig that MINIMISES the
        bin's indel burden = the true physical length (corrects the ruler's
        ~4-unit AL->AH under-count, e.g. measured 66 -> true 70)
        |
        v
[4] Read recovery - assign EVERY VNTR-covering read to an allele, reliably
        (known allele by length; only ambiguous internal partials by best-fit)
        |
        v
[5] Haplotype units
        - distinct-length alleles -> one unit each
        - same-length OR <=3-repeats-apart alleles -> PHASE by SNP into two
          haplotypes, then RE-MEASURE each haplotype's own length (recovers
          e.g. 43 | 44 even when the length step clustered them together)
        - het/hom gate: two haplotypes are kept if they differ in LENGTH or
          by > --min-snp-diff SNPs; only two SAME-length haplotypes differing
          by <= that are collapsed to one (a genuine homozygote / noise split)
        |
        v
[6] Per-haplotype call - consensus, motif nomenclature, frameshift,
        novel-motif discovery (muc1_analyzer call)
        |
        v
[7] 59dupC - statistical run-length-shift caller by default; optional
        --clair3-vcf or --genomic-bam authorities; power-aware verdict
        |
        v
[8] rs4072037 severity - read NATIVELY off the MUTANT allele's reads
        (majority C/T at the constant flank offset; no GRCh38/T2T needed)
        |
        v
[9] MUC1_Score (two axes) + cohort z-score -> merged both-allele outputs
```

---

## 3. Installation / requirements

* **Python 3.9+**
* `pysam`, `mappy` (minimap2 binding), `reportlab`
* The **`muc1_analyzer` package** and **`vntr_raw_length.py`** must be importable -> **run from the repository root** (branch `Dev-Muc1-score` / `muc1-analyzer-fromfastq`, which carry them).
* The **faked reference** `MUC1_fakedVNTR1to150revcomplKirby.fa` (150 contigs `MUC1_VNTR_1repeats ... _150repeats`), passed with `-r` (absolute path recommended). Not versioned in the repo - supply it yourself.

```bash
pip install pysam mappy reportlab
```

---

## 4. Usage

```bash
python MUC1_Analyzer_fromfastq.py \
    -b  reads.fastq.gz \
    -r  /abs/path/MUC1_fakedVNTR1to150revcomplKirby.fa \
    -s  SAMPLE_ID \
    -o  VNTR_results/SAMPLE_ID \
    --call-min-depth 2 \
    -t 20
```

* `-b` accepts a **FASTQ(.gz)**, uBAM, or any BAM/CRAM (only read *sequences* are used).
* **rs4072037 is now detected automatically** off the mutant allele's reads - do **not** pass `--rs4072037` unless you want to override it.
* `-r` is resolved to an absolute path and indexed once at start-up (clear early error if missing / read-only).
* The output directory is created automatically. **Quote** paths with spaces or `[]`, e.g. `-o "VNTR_results/[AS]SAMPLE"`.

### Options (defaults in brackets)

**Core**

| Option | Purpose |
|---|---|
| `-b, --bam` *(req.)* | reads: fastq(.gz), uBAM, or any BAM/CRAM (sequences only) |
| `-r, --ref` *(req.)* | faked multi-length reference (`MUC1_VNTR_Nrepeats`) |
| `-s, --sample` `[sample]` / `-o, --outdir` `[.]` | sample name / output directory |
| `-t, --threads` `[1]` | threads for the ruler alignment (bottleneck on large AS FASTQs; scales well, e.g. `-t 20`) |
| `--pacbio` | PacBio HiFi (`map-hifi`; default `map-ont`) |
| `--ref-cram` | reference FASTA if `--bam` is a CRAM |

**rs4072037 severity (new)**

| Option | Purpose |
|---|---|
| `--rs4072037 {C,T}` | **override** the mutant allele's rs4072037 base (default: detected). If given it wins, but a disagreement with the detected base is flagged on the report |
| `--no-rs4072037-detect` | disable native detection; rely only on `--rs4072037` |
| `--rs4072037-on-negative` | also read rs4072037 (sample-level, informational) when no frameshift is found |

**Length & alleles**

| Option | Purpose |
|---|---|
| `--sep` `[8]` / `--merge-tol` `[3]` / `--min-reads-per-allele` `[2]` | allele clustering |
| `--hom-ratio-near` `[0.5]` | a NEIGHBOURING second length-allele (delta < `--sep`) is folded as scatter below this fraction of the dominant's support |
| `--hom-ratio-far` `[0.25]` | a DISTANT second length-allele (delta >= `--sep`) folded as scatter below this fraction (absolute floor = `--min-reads-per-allele`) |
| `--no-scaffold-refine` / `--scaffold-window` `[6]` | the indel-burden scaffold sweep (true physical length) |
| `--maxmm` `[9]` / `--offset` `[0]` / `--contig-offset` `[0]` / `--ruler-contig` | length calibration knobs |

**Same-length / close-length phasing (new)**

| Option | Purpose |
|---|---|
| `--phase-close-within` `[3]` | when two alleles are called <= this many repeats apart, separate them by SNP phasing (not length-binning), then re-measure each haplotype's own length |
| `--no-phase-hom` | disable SNP phasing of a length-homozygous sample |
| `--min-snp-diff` `[3]` | two SAME-length phased haplotypes are collapsed to one (homozygous) when their consensuses differ by <= this many SNPs; haplotypes of DIFFERENT length are always kept, whatever the SNP count |
| `--min-phasing-snps` `[2]` / `--snp-min-af` `[0.30]` | pileup pre-gate to *attempt* a split (the het/hom decision is `--min-snp-diff`) |

**Read recovery & dupC**

| Option | Purpose |
|---|---|
| `--min-vntr-copies` `[1]` | min VNTR copies for a read to be assigned (excludes 0-copy off-target; raise to ~5-10 to drop tiny fragments) |
| `--no-recover-reads` | use only flank-anchored spanning reads per allele |
| `--call-min-depth` `[3]` | min depth for a consensus base (use `2` for a shallow long allele) |
| `--call-ins-threshold` `[0.5]` | min read fraction to keep a <=3 bp insertion in the consensus nomenclature |
| `--clair3-vcf` | Clair3 VCF (vs the faked ref) -> authoritative dupC verdict per allele |
| `--genomic-bam` / `--genome-ref-dupc` / `--dupc-pon` | delegate the dupC call to `muc1_analyzer.dupc_dispatch` at the MUC1 locus |
| `--vntr-start` `[4574]` / `--repeat-offset` `[4]` | coordinate conventions for the reported repeat index |

**Misc**

| Option | Purpose |
|---|---|
| `--keep-intermediates` | keep the per-allele files (BAMs, per-allele JSON/consensus/PDFs) for debugging |
| `--stage1-only` | stop after the length call |

---

## 5. Output files

By default **only the final, both-allele files are kept** (per-allele intermediates - including the misleading single-allele PDFs - are removed; `--keep-intermediates` retains them).

| File | Content |
|---|---|
| `SAMPLE.MUC1_score.pdf` | **The clinical report** - clinical call, variant motif **in red** (`X(59dupC)` / `X-59dupC`, never doubled), both haplotypes, and a **MUC1_Score** section: score table + cohort calibration curve with the sample's **z-score** (onset tardif <- -> precoce) + a **GENOTYPE ASSOCIE A UN PHENOTYPE DE SEVERITE / PROTECTEUR** flag, the rs4072037 source (detected counts or override), and a warning if an override disagrees with the detected base |
| `SAMPLE.vntr.bam` (+ `.bai`) | one coordinate-sorted BAM with **both** haplotypes - load in IGV against the faked reference |
| `SAMPLE.consensus.fa` | consensus of **both** alleles |
| `SAMPLE.novel_motifs.tsv` | 60 bp motifs not yet in `KNOWN_REPEATS`, aggregated |
| `SAMPLE.merged_haplotypes.json` | merged analyzer JSON (both haplotypes), input to the score |
| `SAMPLE.cigar.summary.json` | full machine-readable summary (alleles, scatter-fold, phasing + `consensus_snp_diffs`, dupC verdict, `rs4072037_detection` / `rs4072037_used`, score) |

---

## 6. Toy example (a control carrier)

Input: a MUC1 carrier - short allele 36 repeats, long allele 70 repeats carrying the 59dupC.

```
36 M reads processed; ~1 M aligned to the ruler; ~50 cover the VNTR (kept)
alleles (HET): [32, 66] copies
allele 32 -> scaffold 36 ; allele 66 -> scaffold 70
dupC (run-length shift): CONFIRME - repetition 20 on MUC1_VNTR_70repeats
rs4072037 (mutant allele): T (C:0 / T:10, depth 10) -> protecteur [detecte]
MUC1_Score: 36 repeat | 70 repeat (59dupC ~ repeat 20)  onset=1.429 severity=-1.0
```

> Note: rs4072037 is read from the mutant allele's reads (here **T -> protective**). Passing `--rs4072037 C` would override to severe **and** print a disagreement alert on the PDF.

`cigar.summary.json` (excerpt):

```json
{
  "alleles": [32, 66],
  "clinical_call": "36 repeat | 70 repeat (59dupC ~ repeat 20)",
  "reliable_dupc": {"called": true, "interpretation": "CONFIRMED",
                    "carrier_contig": "MUC1_VNTR_70repeats",
                    "variant": {"repeat": 20, "label": "59dupC", "p_bonf": 4.66e-06}},
  "rs4072037_detection": {"base": "T", "C": 0, "T": 10, "depth": 10},
  "rs4072037_used": {"base": "T", "source": "detecte (C:0 / T:10)"},
  "muc1_two_axis_score": {"mut_len": 70, "healthy_len": 36,
                          "onset_score": 1.429, "severity_score": -1.0}
}
```

---

## 7. Method notes (design rationale)

* **Ruler-based length.** Length comes only from reads carrying **both** unique flanks - exact and immune to the near-identical-tandem smear. One-flank reads are a lower bound (depth), never a length.
* **Scaffold by indel burden.** The `AL->AH` anchors sit ~4 units inside the tandem, so the ruler under-reads; each allele is scaffolded on the faked contig its reads fit with the lowest indel burden - the true physical length.
* **Scatter-fold.** A second length-peak close to the dominant one is likely measurement scatter, so it must reach a higher relative support (`--hom-ratio-near`, 0.5) than a distant one (`--hom-ratio-far`, 0.25) to be kept as a real allele; the absolute floor against tiny distant noise is `--min-reads-per-allele`.
* **Same-length / close-length phasing.** When the length can't separate two alleles (identical, or one repeat apart clustered together), they are split by **SNP phasing** on the scaffold (the reference caller's `find_phasing_snps` + `phase_reads`), and each phased haplotype is **re-measured for its own length** - so `43 | 44` shows even when the length step merged them. The 59dupC is then attributed to the carrier haplotype.
* **Het vs homozygote — length first, then SNPs.** Two phased haplotypes are distinct if they differ in **re-measured length** OR by **more than `--min-snp-diff` (3) SNPs** between their consensuses. Only when the length does NOT separate them (same scaffold) AND they differ by <= that are they treated as one allele (a true homozygote, or a noise split) and collapsed. So a real close-length carrier (e.g. 43 | 44, few substitutional SNPs) is kept by its length difference, while a same-length near-identical split is folded. Read off the consensus, the SNP count is robust to per-read ONT error.
* **Run-length-shift dupC.** ONT can't count a homopolymer, so the caller scores the fraction of reads whose C-tract is **>= mut_len** against an **in-sample empirical null**, with a binomial + Bonferroni test - general, self-calibrated, needing no PoN or Clair3.
* **rs4072037 read natively.** The SNP sits in the constant 5' flank present in every faked contig, so its base is read straight off the mutant allele's reads (majority C/T at a constant offset, gene-strand mapped back to C/T) - **no GRCh38/T2T alignment needed**. On a negative there is no cis allele, so it is read only with `--rs4072037-on-negative`; below ~5 reads it is reported indeterminate rather than defaulted to severe.
* **Variant repeat number = motif index.** The reported repeat is the 1-based motif index carrying the variant, consistent with the total length count.
* **MUC1_Score.** Two axes - onset (`onset_index` = fraction of the mutant VNTR downstream of the frameshift = fraction translated as the MUC1fs neoprotein) and severity (rs4072037: C -> +1, T -> -1) - plus a cohort z-score positioning the sample against the n=27 calibration survey.

---

## 8. Notes & caveats

* The **rs4072037 severity axis weights and the cohort calibration are provisional** (n=27, LOO rho ~ 0.52) - a direction-of-effect signal, not a validated predictor.
* **Phasing sensitivity**: two SAME-length alleles differing by <= `--min-snp-diff` substitutional SNPs are reported as a single (homozygous) haplotype by design (their dupC, if any, is still detected in the pool but not attributed). Alleles of different length are always separated. Tune `--min-snp-diff` / `--min-phasing-snps` on known cases.
* **PCR long-read amplicons** are out of scope: stutter/chimeras break the burden sweep (use the dedicated tool for those).
* The homopolymer caller's **specificity** should still be confirmed on a well-covered true negative.
* Everything above is validated on real and synthetic controls; per-sample thresholds may need tuning on your cohort.

---

## 9. Where it fits in the lab pipeline

Use `MUC1_Analyzer_fromfastq.py` as the **FASTQ entry point** in place of the
`minimap2 (faked ref) -> mosdepth -> Clair3 -> vntr_haplotype_caller` chain when you
start from adaptive-sampling reads and want a single both-allele result. If you
already produce a Clair3 VCF against the faked reference, pass it with
`--clair3-vcf` to use it as the authoritative dupC verdict; the rest of the
report is produced identically.
