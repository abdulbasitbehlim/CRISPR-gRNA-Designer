# Changelog

## 3.1.0 — 2026-09-08

### Added
- True Ensembl canonical-transcript **TSS-aware CRISPRi** design.
- Strand-aware genomic sequence retrieval around the TSS.
- CRISPRi candidate filtering from **−50 to +300 bp** relative to TSS.
- Preferred CRISPRi placement band at **+50 to +100 bp**.
- Genomic assembly, chromosome, transcript, TSS coordinate, genomic guide coordinates and TSS distance in CRISPRi output.
- Explicit TSS coordinate input for pasted-sequence CRISPRi mode.
- Regression tests for TSS bands, plus/minus strand genomic coordinate conversion and TSS validation.

### Changed
- Restored the polished v2-style Streamlit dashboard while retaining the v3 scientific backend.
- CRISPRi no longer uses a transcript/cDNA first-400-bp positional proxy.
- CRISPRi ranking now prioritizes TSS placement, then specificity/on-target metrics.
- Updated README and in-app methods/limitations to distinguish solved and remaining limitations.

### Fixed
- Removed the misleading limitation that CRISPRi is not TSS-aware.
- Prevented pasted-sequence CRISPRi from silently assuming a TSS.

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

### Fixed
- Multi-FASTA references are no longer concatenated across record boundaries.
- Additional exact target copies remain critical off-targets while one intended exact match is excluded.
