# CRISPR Studio — gRNA Designer v3.2.0

CRISPR Studio is an open-source Python/Streamlit workbench for designing and ranking **SpCas9 (20 nt + NGG)** guide RNAs for knockout and TSS-aware CRISPRi workflows.

## Live app

**Launch CRISPR Studio:** https://crispr-grna-designer-v6mhgxd4o3eqbhgur3anvh.streamlit.app/

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://crispr-grna-designer-v6mhgxd4o3eqbhgur3anvh.streamlit.app/)

> The hosted app runs from the repository's `main` branch. If Streamlit is waking from sleep, the first load may take a short moment.

## v3.2.0 highlights

- Professional dark and light Streamlit dashboard with readable controls, metric cards, plots, validation cards and exports.
- Three target-input modes:
  - **Gene lookup** through NCBI or Ensembl.
  - **Accession ID** lookup through NCBI Nucleotide or Ensembl stable IDs.
  - **Paste sequence** for DNA/RNA/FASTA input.
- **True TSS-aware CRISPRi** for gene lookup using the canonical Ensembl transcript.
- Native **MIT/Hsu** off-target scoring and guide-level specificity.
- Optional **Doench Rule Set 2** on-target scoring and **CFD** off-target scoring through a compatible GuideMaker installation.
- Multi-contig FASTA-safe local-reference screening.
- Optional whole-genome analysis through a local/prebuilt **GuideScan2** index.

## Architecture

```text
app.py
  ├─ grna_designer.py      # SpCas9 discovery, heuristic, MIT, local screening, GuideScan2 adapter
  ├─ crispri.py            # Ensembl TSS retrieval and TSS-aware CRISPRi ranking
  └─ accession_lookup.py   # NCBI / Ensembl accession resolution
```

The Streamlit UI coordinates the workflow but does not reimplement the scientific core.

## Scientific score separation

CRISPR Studio keeps distinct biological concepts separate:

- **Heuristic score** — transparent sequence-quality ranking retained for continuity.
- **Doench Rule Set 2** — on-target activity prediction when the optional provider is installed.
- **MIT/Hsu specificity** — off-target specificity calculated natively.
- **CFD** — per-hit off-target cleavage-risk scoring when the optional provider is installed.

These values are not averaged into one misleading score.

## Target input modes

### Gene lookup

Enter a gene symbol and organism and choose NCBI or Ensembl. Knockout discovery uses a representative transcript/cDNA. For **CRISPRi**, the app switches to Ensembl genomic annotation, resolves the canonical transcript TSS and retrieves strand-aware genomic DNA around that TSS.

### Accession ID

Choose **Accession ID** and enter a record such as:

- `NM_000546.6`
- `NC_000017.11`
- `ENST00000269305`
- `ENSG00000141510`

Auto mode routes `ENS...` identifiers to Ensembl and other nucleotide accessions to NCBI Nucleotide. The app reports the resolved accession, database, record type and sequence length.

Accession-only CRISPRi is intentionally not inferred because a nucleotide accession does not necessarily define a reliable regulatory TSS. Use **Gene lookup** for true TSS-aware CRISPRi, or paste a genomic sequence with an explicit TSS.

### Paste sequence

Paste plain DNA, RNA or FASTA. For CRISPRi, provide the 1-based TSS position and ensure the sequence is genomic DNA in 5′→3′ transcriptional orientation.

## CRISPRi workflow

For **Gene lookup → CRISPRi repression**, CRISPR Studio:

1. resolves the gene in Ensembl;
2. selects the canonical transcript when available;
3. determines transcript strand and genomic TSS;
4. retrieves genomic DNA around the TSS in transcriptional orientation;
5. finds NGG-compatible SpCas9 sites on both strands;
6. keeps guides whose spacer midpoint falls from **−50 to +300 bp** relative to TSS;
7. prioritizes the **+50 to +100 bp** placement band;
8. reports assembly, chromosome, genomic coordinates, transcript and TSS distance.

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

For development:

```bash
pip install -r requirements-dev.txt
pytest -v
```

GitHub Actions tests Python 3.10, 3.11 and 3.12.

## Specificity modes

### Local reference

Useful for plasmids, amplicons, bacterial/viral genomes, contigs and paralog panels. Multi-FASTA contigs remain separate, preventing artificial junction sites. The app scans PAM-compatible sites, reports near matches and calculates MIT/Hsu specificity.

### Whole genome with GuideScan2

Whole-genome analysis requires an external `guidescan` executable plus a matching prebuilt genome index available to the machine running Streamlit.

```bash
conda install -c bioconda guidescan
```

Large genome FASTA/index files should not be committed to this repository.

## Optional Doench Rule Set 2 / CFD provider

The core attempts to use a compatible GuideMaker installation for Doench Rule Set 2 and CFD. When it is not installed, those fields are reported as unavailable rather than being replaced with a mislabeled heuristic.

## Exports

The dashboard exports:

- CSV ranked guide table;
- FASTA synthesis-ready spacer records;
- JSON structured analysis metadata and guide results.

CRISPRi exports also include TSS/genomic provenance when available.

## Current scientific boundaries

- Knockout **gene lookup** still begins from representative transcript/cDNA. Shortlisted guides must be mapped to the intended genomic assembly and coding exon before experimental use.
- Accession lookup retrieves the record requested; it does not automatically infer exon, CDS or promoter biology from arbitrary nucleotide records.
- Variant-aware filtering, chromatin-state integration and DNA/RNA bulge modeling require additional external genome-aware resources.
- GuideScan2 whole-genome search requires an installed executable and prebuilt index.
- Computational scores prioritize candidates but do not replace genome-aware review or experimental validation.

## Primary references

- Jinek et al. 2012 — programmable Cas9 cleavage.
- Hsu et al. 2013 — SpCas9 specificity / MIT scoring.
- Gilbert et al. 2014 — CRISPRi/CRISPRa genome-scale regulation.
- Horlbeck et al. 2016 — CRISPRi/a libraries and TSS-position rules.
- Doench et al. 2016 — Rule Set 2 and CFD.
- GuideScan2 — genome-wide CRISPR guide specificity analysis.

## License

MIT License.
