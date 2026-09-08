# CRISPR Studio v3.2.0

Release date: 2026-09-08

## Summary

CRISPR Studio v3.2.0 consolidates the current application architecture and adds online accession-based sequence retrieval while preserving the v3 scientific backend and TSS-aware CRISPRi workflow.

## Highlights

- Professional, high-contrast dark and light dashboard themes.
- Three input modes: Gene lookup, Accession ID and Paste sequence.
- NCBI Nucleotide and Ensembl accession/stable-ID retrieval with provenance.
- True Ensembl canonical-transcript TSS-aware CRISPRi.
- Native MIT/Hsu specificity.
- Optional Doench Rule Set 2 and CFD providers.
- Multi-contig local-reference off-target screening.
- GuideScan2 whole-genome adapter.
- CSV, FASTA and JSON exports.
- CI coverage on Python 3.10, 3.11 and 3.12.

## Architecture

- `app.py` — Streamlit UI and orchestration
- `grna_designer.py` — guide discovery/scoring/specificity
- `crispri.py` — TSS-aware CRISPRi
- `accession_lookup.py` — NCBI/Ensembl accession retrieval

## Important scientific notes

- Knockout gene lookup remains representative transcript/cDNA-first; final guides should be mapped to the intended genomic coding exon and assembly.
- Accession-only CRISPRi is not inferred because arbitrary nucleotide records may not define a valid TSS.
- Whole-genome specificity requires GuideScan2 plus a matching external index.
- Doench Rule Set 2 and CFD remain optional provider-backed metrics.
- Computational shortlisting does not replace experimental validation.

## Suggested GitHub release title

`CRISPR Studio v3.2.0`

## Suggested tag

`v3.2.0`
