# CRISPR Studio — gRNA Designer v3

An open-source Python/Streamlit workbench for designing and ranking **SpCas9 (20 nt + NGG)** guide RNAs.

Version 3 separates three concepts that should not be conflated:

- **on-target activity** — legacy transparent heuristic, plus optional **Doench Rule Set 2**
- **off-target pair risk** — **MIT/Hsu** and optional **CFD** scores
- **off-target search** — either a user-supplied **local reference** or an indexed **whole genome with GuideScan2**

> Research-use software. Always validate genomic coordinates, genome assembly, exon/TSS context, off-targets, and experimental controls before ordering guides.

## What is new in v3.0.0

- Native MIT/Hsu off-target scoring using the published positional mismatch weights.
- Guide-level MIT specificity aggregation.
- Optional Doench 2016 Rule Set 2 scoring through a compatible GuideMaker provider when installed.
- Optional Doench 2016 CFD off-target scoring through GuideMaker when installed.
- Multi-FASTA local references preserve contig/chromosome boundaries instead of concatenating them.
- Local-reference hits now report contig, position, strand, mismatch positions, MIT risk, optional CFD risk, and risk tier.
- Whole-genome mode uses a local/prebuilt **GuideScan2** index rather than pretending that a short sequence scan is genome-wide.
- Clear separation of Local Reference and Whole Genome analysis in the UI.
- Expanded CSV/FASTA/JSON export fields.
- New regression tests for MIT scoring and multi-contig references.

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

The default install is intentionally lightweight. MIT scoring and local-reference analysis work without additional CRISPR packages.

## Optional advanced scoring

The v3 core can use the USDA **GuideMaker** implementations of Doench Rule Set 2 and CFD when that package and its model/data files are available in the Python environment. If the provider is unavailable, the application reports `N/A` and **does not relabel the legacy heuristic as Doench**.

The relevant scientific methods are:

- Hsu et al. 2013 — MIT/Hsu off-target scoring
- Doench et al. 2016 — Rule Set 2 on-target activity and CFD off-target activity

## Whole-genome mode with GuideScan2

Whole-genome analysis requires the external `guidescan` executable and a prebuilt genome index.

Typical installation with Bioconda:

```bash
conda install -c bioconda guidescan
```

Build an index from a genome FASTA:

```bash
guidescan index genome.fa
```

Then launch CRISPR Studio and choose:

**Specificity analysis → Whole genome (GuideScan2)**

Enter the path/prefix of the GuideScan2 index. The v3 backend runs `guidescan enumerate` for the current candidate guides. Prebuilt indexes can also be used when available.

Large human/mouse genome indexes should live outside this Git repository. Do not commit genome FASTA files or indexes to GitHub.

## Local-reference mode

Use Local Reference mode for:

- bacterial genomes
- viral genomes
- plasmids
- synthetic constructs
- amplicons
- paralog panels
- selected genomic regions

FASTA records are preserved separately:

```text
>chr1
...
>chr2
...
```

v3 never creates an artificial `chr1 → chr2` junction.

For each retained near-match, the report can contain:

- contig
- 1-based position
- strand
- candidate spacer
- PAM
- mismatch count
- mismatch positions
- seed mismatches
- MIT off-target score
- CFD off-target score, when provider is installed
- risk tier

## MIT score in v3

For a guide/off-target pair, v3 uses the Hsu/MIT positional mismatch-weight formulation. The pair score is between 0 and 1, where larger values indicate a more concerning off-target pair.

A guide-level specificity score is then computed as:

```text
MIT specificity = 100 / (1 + sum(pair MIT scores))
```

Higher guide-level specificity is better.

An additional perfect genomic/reference copy therefore strongly lowers specificity.

## Doench Rule Set 2

Rule Set 2 is an **on-target efficiency** model, not an off-target score. It requires sequence context around the protospacer. v3 constructs a 30-nt context when enough flanking sequence is available and passes it to the optional Rule Set 2 provider.

When the provider is absent or the required context is unavailable, the value remains `N/A`.

## CFD

CFD is a **guide/off-target cleavage-risk** score from Doench et al. 2016. v3 can calculate it through the optional GuideMaker provider. The local report also derives a CFD-based specificity summary when CFD scores are available.

## Current workflow

```text
Gene / pasted DNA
        |
        v
20 nt + NGG candidate discovery on both strands
        |
        +--> legacy transparent activity heuristic
        |
        +--> optional Doench Rule Set 2
        |
        v
candidate guides
        |
        +--> Local Reference
        |       |
        |       +--> MIT/Hsu
        |       +--> optional CFD
        |
        +--> Whole Genome (GuideScan2 index)
                |
                +--> genome-wide indexed enumeration
```

## Important limitations that remain

Version 3 fixes several major scoring/reference issues, but it does **not** claim to solve all CRISPR-design biology.

1. **Gene lookup currently starts from representative transcript/cDNA.** A candidate can in principle cross an exon-exon junction. Shortlisted guides must be mapped back to the intended genome build and exon before experimental use.
2. **CRISPRi is not yet truly TSS-aware.** The current knockdown mode still uses a simple 5-prime positional preference. A future version should use genomic TSS annotations and CRISPRi-specific activity models.
3. **Variant-aware filtering is not yet included.** dbSNP/VCF population variation can change both on-target and off-target sites.
4. **Bulges are not modeled by the lightweight local scanner.** GuideScan2 can support richer genome-wide search options outside the lightweight hosted workflow.
5. **Chromatin, cell type, epigenetic state, essential domains, frameshift probability, and repair outcome models are not currently incorporated.**
6. **A high computational score does not guarantee biological activity or safety.** Experimental validation is required.

## Testing

Run the pure-Python core tests with:

```bash
pytest -q test_grna_designer.py
```

The v3 core regression suite covers candidate generation, MIT pair scoring, guide-level MIT specificity, exact-match exclusion, duplicate perfect targets, and multi-FASTA contig preservation.

Run the Streamlit smoke tests with:

```bash
pytest -q test_app.py
```

## Files

- `app.py` — Streamlit interface
- `grna_designer.py` — design/scoring/off-target core
- `test_grna_designer.py` — offline scientific/core tests
- `test_app.py` — Streamlit smoke tests
- `ARCHITECTURE.md` — project architecture notes
- `USER_GUIDE.md` — user documentation
- `requirements.txt` — default Python dependencies

## Scientific references

- Jinek M, et al. (2012). *A programmable dual-RNA-guided DNA endonuclease in adaptive bacterial immunity.* Science. https://doi.org/10.1126/science.1225829
- Hsu PD, et al. (2013). *DNA targeting specificity of RNA-guided Cas9 nucleases.* Nature Biotechnology 31:827–832. https://doi.org/10.1038/nbt.2647
- Doench JG, et al. (2016). *Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9.* Nature Biotechnology 34:184–191. https://doi.org/10.1038/nbt.3437
- Schmidt H, et al. GuideScan2 / GuideScan software for genome-indexed CRISPR off-target enumeration. https://github.com/pritykinlab/guidescan-cli
- USDA-ARS GuideMaker, including a maintained ONNX implementation of the Doench 2016 model. https://github.com/USDA-ARS-GBRU/GuideMaker

## License

See [LICENSE](LICENSE).

## Version

Current release line: **3.0.0**
