# Changelog

## 3.0.0 — 2026-09-08

### Added
- Native MIT/Hsu pair scoring and guide-level specificity.
- Optional Doench Rule Set 2 on-target scoring provider.
- Optional CFD off-target scoring provider.
- Multi-contig FASTA-aware local-reference screening.
- Contig-aware off-target reporting.
- GuideScan2 whole-genome backend integration.
- Whole-genome/local-reference mode separation in the Streamlit UI.
- Regression tests for MIT scoring and FASTA contig preservation.

### Changed
- Renamed the old activity value as a transparent legacy heuristic rather than implying it is Doench Rule Set 2.
- Replaced the old custom specificity formula as the primary local specificity metric with MIT/Hsu specificity.
- Expanded export fields for scientific provenance.
- Updated documentation for optional genome-indexed workflows and current limitations.

### Fixed
- Multi-FASTA references are no longer concatenated across record boundaries.
- Additional exact target copies remain critical off-targets while one intended exact match is excluded.

### Known limitations
- Gene lookup is still transcript/cDNA-first and requires genomic exon/build validation.
- CRISPRi mode is not yet true TSS-aware design.
- Local-reference mode does not model bulges, variants, or chromatin.
- Whole-genome GuideScan2 requires a local/server installation and prebuilt genome index.
