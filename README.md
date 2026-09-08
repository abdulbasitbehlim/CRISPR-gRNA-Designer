# CRISPR Studio — gRNA Designer v3.1

An open-source Python/Streamlit workbench for designing and ranking **SpCas9 (20 nt + NGG)** guide RNAs for knockout and CRISPRi experiments.

## What v3.1 adds

- Restored polished dark/light Streamlit dashboard with guide cards, plots, detailed inspection, exports and methods/limitations views.
- **True TSS-aware CRISPRi gene design** using the canonical Ensembl transcript rather than a cDNA 5′-position proxy.
- Strand-aware genomic TSS-window retrieval from Ensembl.
- CRISPRi filtering to **−50 to +300 bp relative to the annotated TSS**, with **+50 to +100 bp** treated as the preferred placement band for dCas9-KRAB-style repression.
- Pasted-sequence CRISPRi requires an explicit TSS coordinate instead of guessing one.
- Genomic assembly, chromosome, transcript, TSS coordinate, guide genomic coordinates and TSS distance are exported for gene-based CRISPRi designs.
- Native MIT/Hsu off-target scoring and guide-level specificity.
- Optional Doench Rule Set 2 on-target scoring and Doench CFD off-target scoring through a compatible GuideMaker installation.
- Multi-contig FASTA-safe local-reference screening.
- Optional whole-genome off-target analysis through a local/prebuilt GuideScan2 index.

## Scientific score separation

CRISPR Studio intentionally keeps different biological concepts separate:

- **Heuristic score** — transparent sequence-quality ranking retained for continuity.
- **Doench Rule Set 2** — on-target cutting/activity prediction when the optional provider is installed.
- **MIT/Hsu specificity** — off-target specificity calculated natively.
- **CFD** — off-target cleavage likelihood/specificity when the optional provider is installed.

The app does not average these into a misleading single score.

## CRISPRi workflow

For **Gene lookup → CRISPRi repression**, CRISPR Studio:

1. Resolves the gene in Ensembl.
2. Selects the canonical transcript when available.
3. Determines the transcript strand and true genomic TSS.
4. Retrieves genomic DNA around that TSS in transcriptional orientation.
5. Finds SpCas9 NGG sites on both strands.
6. Keeps candidates whose spacer midpoint falls from **−50 to +300 bp** relative to the TSS.
7. Prioritizes the **+50 to +100 bp** band, then ranks using specificity and on-target metrics.
8. Reports assembly, chromosome, genomic coordinates, strand, transcript and TSS distance.

For **Paste sequence → CRISPRi repression**, the user must provide the TSS position. The pasted sequence must be genomic DNA in 5′→3′ transcriptional orientation.

## Installation

```bash
git clone https://github.com/abdulbasitbehlim/CRISPR-gRNA-Designer.git
cd CRISPR-gRNA-Designer
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell
# .venv\Scripts\Activate.ps1
pip install -r requirements.txt
streamlit run app.py
```

## Whole-genome mode with GuideScan2

Whole-genome analysis requires the external `guidescan` executable and a matching prebuilt genome index.

```bash
conda install -c bioconda guidescan
guidescan index genome.fa
```

Then choose **Specificity analysis → Whole genome (GuideScan2)** and provide the index path available to the server running Streamlit.

Large genome FASTA/index files should not be committed to this repository.

## Local-reference mode

Local-reference analysis is useful for bacterial genomes, viral genomes, plasmids, amplicons, contigs, paralog panels and other custom references. Multi-FASTA records are preserved independently so artificial contig-boundary targets are not created.

## Remaining scientific boundaries

- Knockout gene lookup still begins from a representative transcript/cDNA; shortlisted guides should be mapped to the intended genomic coding exon and assembly before experimental use.
- Variant-aware filtering, chromatin-state integration and DNA/RNA bulge modeling require additional external genomic resources/workflows.
- GuideScan2 whole-genome search requires a server/local environment where its executable and genome index are installed; ordinary Streamlit Community Cloud does not automatically provide these large indexes.
- Computational scores prioritize candidates but do not replace experimental validation.

## Tests

The repository includes core and Streamlit smoke tests. GitHub Actions runs the test suite on Python 3.10, 3.11 and 3.12.

## Primary references

- Jinek et al. 2012 — programmable Cas9 cleavage.
- Hsu et al. 2013 — SpCas9 specificity / MIT scoring.
- Gilbert et al. 2014 — CRISPRi/CRISPRa genome-scale regulation.
- Horlbeck et al. 2016 — compact CRISPRi/a libraries and TSS-position rules.
- Doench et al. 2016 — Rule Set 2 and CFD.
- GuideScan2 — genome-wide CRISPR guide specificity analysis.

## License

MIT License.
