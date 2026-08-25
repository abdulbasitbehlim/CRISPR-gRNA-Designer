# CRISPR Studio â€” gRNA Designer

An interactive Python and Streamlit workbench for discovering, ranking, reviewing, and exporting **SpCas9 (NGG)** guide RNA candidates.

**Live app:** [Open CRISPR Studio](https://crispr-grna-designer-v6mhgxd4o3eqbhgur3anvh.streamlit.app/)

> Research and educational use only. This application creates a candidate shortlist; it does not replace genome-wide specificity analysis, genomic/exon mapping, or experimental validation.

## What is new in version 2.0

- Polished responsive dashboard with a built-in **dark/light mode toggle**
- Gene lookup from **Ensembl** or **NCBI**, plus DNA/RNA/FASTA paste mode
- Both-strand SpCas9 discovery with a transparent **0-100 activity ranking**
- Interactive guide landscape, GC distribution, strand balance, and score gauge
- Guide-by-guide score explanation and validation checklist
- Optional **PAM-aware local-reference similarity screen** on both DNA strands
- Specificity score, mismatch positions, seed-region mismatches, and risk tiers
- Downloadable **CSV, Excel, FASTA, and JSON** outputs
- Clear scientific limitations and experimental reminders inside the app
- Offline unit tests for sequence cleaning, PAM scanning, scoring, validation, and off-target logic

## Analysis workflow

```text
Gene symbol + organism OR pasted sequence
                    â”‚
                    â–¼
        NCBI / Ensembl / manual input
                    â”‚
                    â–¼
        Scan both strands for 20 nt + NGG
                    â”‚
                    â–¼
      Activity ranking + validation checks
                    â”‚
           optional local reference
                    â–¼
      PAM-aware near-match specificity screen
                    â”‚
                    â–¼
    Interactive report + CSV/XLSX/FASTA/JSON
```

## Run locally

Python 3.10 or newer is recommended.

```bash
git clone https://github.com/abdulbasitbehlim/CRISPR-gRNA-Designer.git
cd CRISPR-gRNA-Designer

python -m venv .venv
```

Activate the environment:

```bash
# Windows PowerShell
.venv\Scripts\Activate.ps1

# macOS / Linux
source .venv/bin/activate
```

Install and start the app:

```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
streamlit run app.py
```

Open the local URL printed by Streamlit, normally `http://localhost:8501`.

## Use the dashboard

### Gene lookup

1. Choose **Gene lookup**.
2. Enter an official gene symbol such as `TP53`, `BRCA1`, or `ACTB`.
3. Enter the organism, for example `Homo sapiens`.
4. Choose Ensembl or NCBI and select Knockout or Knockdown/CRISPRi.
5. Adjust the score threshold and result count in the sidebar.
6. Select **Design and analyze guides**.

### Paste a sequence

Choose **Paste sequence** and enter plain DNA, RNA, or FASTA. RNA `U` bases are converted to `T`; standard ambiguous IUPAC bases are represented as `N`. The hosted app accepts target sequences up to 50,000 bp.

### Local-reference off-target screen

Enable **Add local-reference off-target screen**, then upload or paste a FASTA/text reference up to 250,000 bp. The app:

- finds NGG-compatible candidate sites on both strands;
- compares each reference spacer with the query spacer;
- excludes one exact site as the presumed intended target;
- reports additional exact or near matches;
- shows mismatch positions and mismatches in positions 13-20 (the PAM-proximal seed region);
- assigns an explainable risk tier and local-reference specificity score.

This screen is useful for a plasmid, amplicon, contig, paralog panel, or small reference region. It is intentionally **not presented as a whole-genome off-target search**. Confirm final guides with CRISPOR, CHOPCHOP, Cas-OFFinder, GuideScan, or an equivalent genome-indexed workflow.

## Cross-check results with independent online tools

CRISPR Studio produces a transparent candidate shortlist. Before ordering any guide, cross-check the shortlisted candidates with at least one genome-indexed platform using the same organism, genome assembly, nuclease (**SpCas9**), PAM (**NGG**), and target locus.

> **Important:** This is independent computational cross-checking, not experimental validation. Different tools use different reference data, algorithms, and score scales, so their numerical scores are not directly interchangeable.

| Online tool | Useful cross-check |
|---|---|
| [CRISPOR](https://crispor.org/) | Compare guide placement and on-target ranking, then inspect genome-context and predicted off-target results. |
| [CHOPCHOP](https://chopchop.cbu.uib.no/) | Redesign from a gene, genomic coordinates, or sequence and compare candidate placement and predicted off-targets. |
| [CRISPick](https://portals.broadinstitute.org/gppx/crispick/public) | Obtain an independent candidate ranking for CRISPR knockout, activation, or interference workflows. |
| [GuideScan2](https://guidescan.com/py/) | Check genome-aware specificity for supported assemblies or search directly for exported gRNA sequences. |
| [CRISPRdirect](https://crispr.dbcls.jp/) | Perform an additional sequence-based guide-selection check focused on reducing unintended targets. |
| [Cas-OFFinder](https://www.rgenome.net/cas-offinder/) | Search a selected reference genome for potential off-target sites with a configurable mismatch limit. |

### Suggested cross-check procedure

1. Export the ranked guide table as CSV or FASTA.
2. Confirm the intended organism, genome assembly, transcript isoform, coding exon, and genomic locus. If CRISPR Studio retrieved cDNA, map each guide to genomic DNA first because an exon-junction candidate may not exist as a continuous genomic target.
3. Submit the exact 20 nt spacer in 5-prime to 3-prime orientation and select SpCas9 with an NGG PAM. When a service requests a target sequence, use the relevant genomic sequence rather than relying only on cDNA.
4. Confirm that the PAM, strand, and genomic coordinate identify the same candidate before comparing results.
5. Review predicted off-targets, especially additional exact matches, low-mismatch sites, coding-region hits, and repetitive or multi-mapping candidates.
6. Prefer multiple independent guides that remain acceptable after genome-aware review, then validate editing efficiency and specificity experimentally.

Agreement across independent tools can strengthen confidence in a computational shortlist, but disagreement should be investigated. No online score proves that a guide will be efficient or safe in a biological experiment.

## Python API

```python
from grna_designer import design_from_gene

accession, description, sequence, guides = design_from_gene(
    gene_name="TP53",
    organism="Homo sapiens",
    source="ensembl",
    application="knockout",
    min_score=40,
    max_guides=10,
)

for guide in guides:
    print(guide.to_dict())
```

Use your own sequence:

```python
from grna_designer import design_guides

guides = design_guides(
    sequence=my_target_sequence,
    application="knockout",
    min_score=35,
    max_guides=20,
    genome_context=optional_reference_sequence,
    max_mismatches=3,
)
```

Inspect one local-reference screen directly:

```python
from grna_designer import analyze_offtargets

report = analyze_offtargets(
    spacer="GCTAGCTAGCTAGCTAGCTA",
    whole_genome=my_reference_sequence,
    max_mismatches=3,
)

print(report.specificity_score)
for hit in report.hits:
    print(hit.to_dict())
```

## How the activity ranking works

The activity value is a transparent heuristic inspired by published SpCas9 design observations. It combines:

- GC content, preferring approximately 40-70%;
- selected position-specific nucleotide preferences near the PAM;
- PAM context;
- poly-G/poly-T and homopolymer penalties;
- a simple self-complementarity penalty;
- an application-position adjustment.

Knockout mode mildly favors candidates in the early 30% of the input sequence. Knockdown/CRISPRi mode favors the first 400 bp as a simple 5-prime proxy. Every non-zero contribution can be inspected in the dashboard.

The result is **not the trained Doench Rule Set 2/Azimuth model**, even though some feature choices are literature-inspired. The value should be used to rank candidates inside this app, not compared directly with scores produced by other tools.

## Scientific limitations

1. Database lookup retrieves a representative transcript/cDNA. A candidate may cross an exon-exon junction or lack genomic, isoform, and regulatory context. Map every spacer back to the intended genome assembly and coding exon.
2. CRISPRi normally requires an experimentally relevant transcription start site window. Transcript position is only an approximation.
3. The lightweight reference screen does not model DNA/RNA bulges, chromatin accessibility, genetic variants, non-NGG PAMs, or genome-scale repeats.
4. A high activity or specificity value does not guarantee editing performance or safety.
5. Use multiple independent guides, non-targeting controls, and orthogonal validation such as amplicon sequencing, ICE/TIDE, RT-qPCR, or Western blot as appropriate.

## Project structure

```text
CRISPR-gRNA-Designer/
â”œâ”€â”€ app.py                     # Streamlit dashboard and export workflow
â”œâ”€â”€ grna_designer.py           # Sequence fetching, design, score, validation, screening
â”œâ”€â”€ test_grna_designer.py      # Offline unit tests
â”œâ”€â”€ requirements.txt           # Streamlit Cloud/runtime dependencies
â”œâ”€â”€ requirements-dev.txt       # Runtime + development dependencies
â”œâ”€â”€ .streamlit/
â”‚   â””â”€â”€ config.toml            # Theme, security, and upload settings
â”œâ”€â”€ .github/workflows/
â”‚   â””â”€â”€ tests.yml              # Continuous integration
â”œâ”€â”€ ARCHITECTURE.md
â”œâ”€â”€ USER_GUIDE.md
â”œâ”€â”€ CONTRIBUTING.md
â””â”€â”€ LICENSE
```

## Tests

```bash
pip install -r requirements-dev.txt
pytest -q
```

All sequence-processing tests run without NCBI or Ensembl network access.

## Deploy on Streamlit Community Cloud

1. Push the updated files to the `main` branch of this GitHub repository.
2. In Streamlit Community Cloud, select **Create app**.
3. Choose this repository and branch.
4. Set the main file path to `app.py`.
5. Deploy or reboot the existing app.

The existing public URL stays the same when the connected repository and app entry point are unchanged.

NCBI asks API clients to identify themselves. For production use, set an `NCBI_EMAIL` environment variable to a monitored contact email. No API key is required for normal low-volume use.

## Primary references

- Jinek M. et al. *Science* (2012), [doi:10.1126/science.1225829](https://doi.org/10.1126/science.1225829)
- Hsu P.D. et al. *Nature Biotechnology* (2013), [doi:10.1038/nbt.2647](https://doi.org/10.1038/nbt.2647)
- Doench J.G. et al. *Nature Biotechnology* (2016), [doi:10.1038/nbt.3437](https://doi.org/10.1038/nbt.3437)
- Moreno-Mateos M.A. et al. *Nature Methods* (2015), [doi:10.1038/nmeth.3543](https://doi.org/10.1038/nmeth.3543)

## Contributing

Contributors are welcome. Please read [CONTRIBUTING.md](CONTRIBUTING.md) for the development setup, test commands, issue-reporting guidance, and pull-request workflow.

Useful contributions include bug fixes, additional tests, documentation improvements, carefully validated scoring methods, support for additional nucleases/PAMs, and optional genome-aware integrations. Scientific changes should clearly state their assumptions, cite the underlying method or data source, include appropriate tests, and preserve the project's limitation notices.

## License

Released under the [MIT License](LICENSE).
