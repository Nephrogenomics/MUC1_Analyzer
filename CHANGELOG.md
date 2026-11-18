# Changelog

## v1.0.1

Fixes from the first external user reports (adaptive sampling and WGS run from GRCh38 CRAMs).
Detection, length arbitration and the MUC1_Score are unchanged.

### Fixed
- **Relative output folder crashed `from-fastq`** (`FileNotFoundError` on the per-allele BAMs). `-o` and
  the optional input paths are now resolved to absolute paths before anything runs.
- **A failed per-allele call no longer produces a verdict.** The run stops with status `FAILED` in
  `.cigar.summary.json`, keeps its intermediates for diagnosis and exits non-zero, instead of writing a
  PDF and a verdict from the other allele only.
- **Aligned BAM/CRAM inputs are no longer accepted.** `run`, `run_pcr`, `prepare` and
  `MUC1_Analyzer_fromfastq.py` take FASTQ or an unaligned BAM (uBAM) only, and refuse an aligned file
  with the command that extracts the MUC1 reads from it. Re-extraction from a GRCh38 CRAM had pulled
  every unplaced read of the file: on an adaptive-sampling CRAM, ~81 000 rejected 300–500 bp fragments
  buried 98 real locus reads and made the sample FLAT; on a WGS CRAM it produced a 569 MB read file.
- **"PDF report written" was printed when reportlab was missing.** Success is now reported only when the
  file exists; otherwise the run says that no PDF was written and why.
- **`vntr_raw_length.py --help` crashed** on an unescaped `%` in a help string.

### Changed
- The former alias of the `muc1` alignment profile is removed: `--align-profile` accepts `setB`, `muc1`
  and `muc1_scoring` only (`muc1` is the default of `run_pcr`).
- Under `⚠ LENGTH UNRELIABLE`, the PDF withholds the per-haplotype motif nomenclature unless the variant
  caller is positive (a consensus on non-spanning reads displayed a spurious `58_59delCC` on a negative
  report). The text report labels its per-motif table "expert review only".
- The read concentration on the two busiest length-contigs is reported (stderr, PDF, `analyzer.json`).
- Every output carries a `provenance` block (tool, version, git commit); the PDF prints it in its footer.
- README and tutorial: FASTQ/uBAM input, extraction command for aligned files, good-practice section.
- README, tutorial and the from-fastq README now describe the **single-axis `MUC1_Score`** (frameshift
  tail, four categories, rs4072037 informational). They still described the former two-axis model.
