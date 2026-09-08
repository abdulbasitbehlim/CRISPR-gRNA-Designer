# Changelog

## 3.2.0 — 2026-09-08

### Added
- Third target-input mode: **Accession ID**.
- Online nucleotide/stable-ID resolution through **NCBI Nucleotide** and **Ensembl**.
- Auto database routing for accession lookup.
- Accession provenance: resolved accession, database, record type and sequence length.
- Professional dark/light UI palette with explicit readable styling for Streamlit inputs, textareas, selects, popovers, tabs and uploaders.
- Dedicated accession lookup module: `accession_lookup.py`.
- Accession normalization and database-routing regression tests.

### Changed
- Repository documentation and architecture now reflect the current four-file application structure: `app.py`, `grna_designer.py`, `crispri.py`, and `accession_lookup.py`.
- Dark/light UI styling was refined for consistent contrast and visibility.
- Validation is presented as readable status cards rather than raw JSON.
- Specificity states distinguish **not screened** from a failed/low-specificity result.
- Accession-only CRISPRi is explicitly blocked from pretending to be TSS-aware; users are directed to Gene lookup or pasted genomic sequence with an explicit TSS.

### Preserved
- Native MIT/Hsu scoring.
- Optional Doench Rule Set 2 and CFD providers.
- Multi-contig local-reference screening.
- GuideScan2 whole-genome adapter.
- TSS-aware CRISPRi for gene lookup.

## 3.1.1 — 2026-09-08

### Changed
- Polished guide-detail dashboard.
- Removed raw validation JSON from the user interface.
- Improved dark/light theme coverage of Streamlit controls.
- Added readable guide-quality cards and styled cloning-oligo presentation.

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
- Restored the polished dashboard while retaining the v3 scientific backend.
- CRISPRi no longer uses a transcript/cDNA first-400-bp positional proxy.
- CRISPRi ranking now prioritizes TSS placement, then specificity/on-target metrics.

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
- Renamed the old activity value as a transparent heuristic rather than implying it is Doench Rule Set 2.
- Replaced the old custom specificity formula as the primary local specificity metric with MIT/Hsu specificity.

### Fixed
- Multi-FASTA references are no longer concatenated across record boundaries during local-reference screening.
- Additional exact target copies remain critical off-targets while one intended exact match is excluded.
