# CRISPR Studio — gRNA Designer v3.3.0

CRISPR Studio is a Python and Streamlit application for helping researchers shortlist **SpCas9 guide RNAs**.

The program looks for **20-nucleotide guide sequences next to an NGG PAM** and supports both:

- CRISPR knockout guide design;
- human/mouse TSS-aware CRISPRi guide design.

The tool is meant to support research and learning. It does **not** guarantee that a guide will work experimentally, and its local specificity screen is not a substitute for complete genome-wide validation.

**Live app:** https://crispr-grna-designer-v6mhgxd4o3eqbhgur3anvh.streamlit.app/

---

## What this tool does

In simple terms, the program follows this workflow:

1. receive a gene, accession, or DNA sequence;
2. find possible 20 nt guides next to an NGG PAM;
3. check basic sequence quality;
4. rank the candidates;
5. optionally compare them with a supplied reference sequence for local near-matches;
6. export the selected results for review.

The software also keeps sequence provenance and validation information so that the result can be traced back to the exact biological input used.

---

## Main workflows

### 1. Knockout design

For NCBI gene lookup, the program retrieves genomic sequence and keeps candidate cut sites inside annotated coding sequence.

It does not create an artificial sequence by joining exons together.

### 2. CRISPRi design

For human and mouse CRISPRi, the program uses a genomic window around a transcription start site and applies a dCas9-KRAB-oriented placement heuristic.

This is a design aid, not a direct measurement of transcriptional repression.

### 3. Manual sequence design

A user can also paste one contiguous genomic DNA sequence and search it directly.

In this case, the user is responsible for knowing where that sequence came from and whether it represents the intended biological locus.

---

## Important behaviour

The software uses several safeguards so that a simple-looking result is not mistaken for stronger evidence than it really provides.

- The intended genomic site is excluded from off-target counting only when the user gives an explicit record, coordinate, and strand that can be verified.
- Duplicate or empty FASTA identifiers are rejected.
- Local specificity counts use all qualifying hits, even when the screen only displays part of the hit table.
- A result with incomplete screening remains marked for review.
- A high-risk non-intended local hit cannot be hidden by a good aggregate specificity score.
- Large inputs are rejected instead of being silently truncated.

---

## Input options

| Input type | What the program does |
|---|---|
| Gene lookup + NCBI + knockout | Uses genomic CDS-aware guide discovery. |
| Gene lookup + CRISPRi | Uses a genomic TSS window for human/mouse CRISPRi. |
| NCBI accession | Accepts genomic DNA accessions. RNA/cDNA accessions are rejected. |
| Ensembl accession | Retrieves contiguous genomic sequence for supported IDs. |
| Paste sequence | Searches one user-supplied contiguous genomic DNA region. |

Accession-only CRISPRi is rejected because a nucleotide accession by itself does not define the intended transcription start site.

---

## Local specificity screen

The optional local screen checks the supplied **linear reference sequence** on both strands.

It currently supports:

- NGG PAMs;
- 0–4 substitution mismatches;
- MIT/Hsu-style specificity calculations;
- bundled CFD scoring;
- explicit intended-site exclusion.

It does **not** model:

- bulges;
- alternate PAMs;
- sample variants;
- chromatin state;
- circular sequence junctions;
- a complete genome unless the user supplies that genome as the reference.

The local screen is therefore best understood as a **reference-limited prioritization step**.

---

## Understanding PASS and REVIEW

The program does not use PASS to mean “experimentally proven”.

A local result can receive **LOCAL CHECKS MET** only when the configured local checks are sufficiently complete and no Critical or High non-intended hit is found.

A result remains **REVIEW** when, for example:

- the mismatch search is too shallow;
- the intended locus has not been declared correctly;
- the reference is ambiguous;
- a Critical or High local hit exists.

---

## Scores

The software uses several scores for different purposes.

### Sequence heuristic

This helps rank knockout candidates using sequence-based properties.

It is not a trained editing-success probability.

### MIT specificity

This summarizes near-match risk in the supplied reference.

It should not be interpreted as a probability of safety.

### CFD

CFD is also used for local near-match assessment.

The bundled weight data and their provenance are documented in the repository.

### Optional Rule Set 2

GuideMaker Rule Set 2 can be reported separately when the required sequence context is available.

It is not used to predict CRISPRi repression.

---

## Installation

Clone the repository:

\`\`\`bash
git clone https://github.com/abdulbasitbehlim/CRISPR-gRNA-Designer.git
cd CRISPR-gRNA-Designer
\`\`\`

Create a virtual environment:

\`\`\`bash
python -m venv .venv
\`\`\`

Activate it.

### Windows PowerShell

\`\`\`powershell
.venv\Scripts\Activate.ps1
\`\`\`

### Linux/macOS

\`\`\`bash
source .venv/bin/activate
\`\`\`

Install the required packages:

\`\`\`bash
python -m pip install -r requirements.txt
\`\`\`

Run the app:

\`\`\`bash
python -m streamlit run app.py
\`\`\`

---

## Running the tests

Install the development requirements:

\`\`\`bash
python -m pip install -r requirements-dev.txt
\`\`\`

Run the test suite:

\`\`\`bash
python -m pytest -q
\`\`\`

The CI workflow also checks supported Python versions automatically.

---

## Main files

The project is divided into smaller modules so that each part has a clear job.

- \`app.py\` — Streamlit user interface.
- \`grna_designer.py\` — main guide discovery and ranking logic.
- \`crispri.py\` — CRISPRi-specific logic.
- \`knockout.py\` — knockout-related genomic handling.
- \`accession_lookup.py\` — biological accession retrieval.
- \`sequence_io.py\` — FASTA and sequence input handling.
- \`models.py\` — shared data structures.
- \`reporting.py\` — result/export helpers.
- \`network.py\` — network-related helper functions.
- \`tests/\` — automated tests.

For a deeper technical explanation, see:

- [SCIENTIFIC_AUDIT.md](SCIENTIFIC_AUDIT.md)
- [USER_GUIDE.md](USER_GUIDE.md)
- [ARCHITECTURE.md](ARCHITECTURE.md)
- [CHANGELOG.md](CHANGELOG.md)

---

## Exports

The application can export reviewable results in formats such as:

- CSV;
- spacer FASTA;
- JSON with settings and provenance.

The JSON export records details such as sequence/reference fingerprints, selected settings, screening scope, hit counts, and genomic/TSS provenance when available.

---

## Whole-genome adapter

The repository includes an optional GuideScan2 adapter.

GuideScan2 requires:

- the external \`guidescan\` executable;
- a matching prebuilt genome index.

Its output is kept separate because GuideScan2 uses its own scoring system and should not be confused with the local MIT/CFD screen.

---

## Scientific limitations

This program helps **prioritize** candidates.

It does not guarantee:

- loss of function;
- frameshift formation;
- activity in a specific cell type;
- complete isoform coverage;
- chromatin accessibility;
- genome-wide safety;
- experimental editing efficiency.

Any final experimental guide should still be checked using appropriate genome-aware tools and biological validation.

---

## Primary references

- Hsu et al. (2013). *DNA targeting specificity of RNA-guided Cas9 nucleases.* DOI: 10.1038/nbt.2647
- Doench et al. (2016). *Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9.* DOI: 10.1038/nbt.3437
- Gilbert et al. (2014). *Genome-scale CRISPR-mediated control of gene repression and activation.* DOI: 10.1016/j.cell.2014.09.029
- Horlbeck et al. (2016). *Compact and highly active next-generation libraries for CRISPR-mediated gene repression and activation.* DOI: 10.7554/eLife.19760
- Poudel et al. (2022). *GuideMaker: Software to design CRISPR-Cas guide RNA pools in non-model genomes.* DOI: 10.1093/gigascience/giac007

GuideScan2 implementation: https://github.com/pritykinlab/guidescan-cli

---

## Citation and license

Software citation metadata is available in [CITATION.cff](CITATION.cff).

Zenodo-ready metadata is available in [.zenodo.json](.zenodo.json). A DOI should only be added after an actual archived release has received one.

The original software is available under the **MIT License**.

CFD weight data keep their separately documented CC0 dedication.
