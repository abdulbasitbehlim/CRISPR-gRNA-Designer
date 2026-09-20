# CRISPR Studio — gRNA Designer v3.3.0

A Python/Streamlit workbench for shortlisting **SpCas9 20 nt + NGG** guides for
knockout and human/mouse TSS-aware CRISPRi.

[Launch CRISPR Studio](https://crispr-grna-designer-v6mhgxd4o3eqbhgur3anvh.streamlit.app/)

**Keywords:** `bioinformatics` · `crispr` · `grna` · `spcas9` · `crispri` · `streamlit` · `computational-biology` · `genome-editing`

## What changed

- **Explicit intended locus:** no first-exact-match removal. Exclusion requires a
  declared contig, position and strand, verified against the reference.
- **Complete local scoring:** every qualifying hit contributes to MIT/CFD and
  counts. Display limits cap rows only.
- **Genomic knockout discovery:** NCBI gene lookup retrieves contiguous genomic
  DNA and keeps cuts inside an annotated CDS. It never joins exons.
- **Strict FASTA:** duplicate/empty identifiers or records are errors. Target
  input accepts one record; reference input preserves multiple contigs.
- **Honest validation:** unscreened, ambiguous-reference and undeclared-locus
  results cannot receive the local-checks-met status.
- **Scope-aware validation:** `LOCAL CHECKS MET` requires a search through at
  least three mismatches and no non-intended Critical or High hit. Aggregate
  MIT specificity alone cannot override a high-risk individual site.
- Bundled CFD weights, corrected optional Rule Set 2 interface, reproducible JSON
  exports, bounded local screening, accession-version checks and additional tests.

See [the scientific audit](SCIENTIFIC_AUDIT.md), [user guide](USER_GUIDE.md),
[architecture](ARCHITECTURE.md) and [changelog](CHANGELOG.md).

## Installation

```bash
git clone https://github.com/abdulbasitbehlim/CRISPR-gRNA-Designer.git
cd CRISPR-gRNA-Designer
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

Tests: `python -m pip install -r requirements-dev.txt`, then `python -m pytest -q`.
The CI matrix covers Python 3.10, 3.11 and 3.12.

## Choosing a workflow

| Input | Behavior and boundary |
|---|---|
| Gene lookup + NCBI + knockout | Genomic CDS-aware design. Choose a RefSeq protein/transcript isoform or use the stated representative selection. |
| Gene lookup + CRISPRi | Ensembl genomic TSS window; canonical or explicitly selected transcript. Human/mouse dCas9-KRAB placement heuristic. |
| Ensembl gene-symbol knockout | Disabled until an annotated genomic implementation is available. Select NCBI. |
| NCBI accession | Genomic DNA only. RNA/cDNA accessions are rejected. Large chromosome records exceed the 2 Mb input cap. |
| Ensembl accession | Gene, transcript or exon IDs retrieve contiguous **genomic** sequence, including introns where applicable. No CDS filter in this mode. |
| Paste sequence | Supply one contiguous genomic DNA region. Its genomic provenance and coding annotation are the user's responsibility. |

Accession-only CRISPRi is rejected because a nucleotide identifier alone does not
establish the intended TSS. For custom CRISPRi, supply a 1-based TSS in genomic DNA
oriented in the direction of transcription.

## Local specificity

The search covers the supplied **linear** reference, both strands, NGG PAMs and
0–4 substitution mismatches. No alternate PAMs, bulges, variants or circular
junctions are searched. Ambiguous sites are skipped and flagged.

For an intended locus, enter the FASTA record ID, the 1-based left edge of the
**whole input region** within that reference and its orientation. The full region
must match; individual guide coordinates are then mapped explicitly. Leave the
record ID blank to retain all exact matches. A wrong declaration is an error.

Twenty exact copies produce MIT specificity **4.76** with no intended locus,
or **5.00** after one verified intended locus is excluded. Both values remain
unchanged when showing 0, 1, 10 or all hit rows. These are reference-limited
prioritization scores, not probabilities or a genome-wide safety certificate.

Every local report records the searched mismatch radius, complete Critical,
High, Moderate and Low hit counts, and maximum per-site MIT and CFD risk. These
summaries use every qualifying hit, including hits omitted from the displayed
detail table. Searches limited to zero, one or two mismatches remain `REVIEW`.
Any non-intended Critical or High hit also requires `REVIEW`, even when the
aggregate MIT specificity is 50 or higher. No universal CFD pass cutoff is used.

Limits: 2 Mb target, 5 Mb local reference, 500,000 indexed NGG sites, 50,000
candidates and 100 million candidate/site comparisons. Oversized work is rejected
with guidance instead of silently truncated screening.

## Scores and exports

- The sequence heuristic ranks knockout candidates; local MIT specificity breaks
  ties. It is not a trained activity model or editing-success probability.
- Optional GuideMaker Rule Set 2 uses an oriented 30-mer and is reported
  separately. Missing context/provider gives N/A. CRISPRi does not use this
  nuclease-activity model to predict repression.
- MIT and bundled CFD aggregate **all** qualifying local hits. CFD data
  provenance and CC0 dedication are in [third_party](third_party/README.md).
- CSV and spacer FASTA support review. JSON includes settings, sequence/reference
  SHA-256 fingerprints, genomic/TSS provenance, full hit counts and capped details.
  Spacer FASTA does not include a scaffold or vector-specific cloning validation.

## Whole-genome adapter

GuideScan2 requires an external executable (`guidescan`) and matching prebuilt
index on the app host. The adapter exposes its raw per-hit CSV output separately;
it does not rerank candidates or convert it into a local-checks-met status.
GuideScan2's native specificity scale and formula differ from the local MIT/CFD
aggregation. No-match output is not proof of a verified intended target. This
release has offline adapter contract tests, not a live genome-index benchmark.

## Scientific boundaries

CDS placement does not guarantee loss of function, frameshift, coverage of every
isoform or experimental activity. Representative selection prefers an NP_ protein,
then longest CDS; it is not a MANE/canonical assertion. CRISPRi does not measure
active cellular TSS or chromatin. No variant-aware or experimental off-target
validation is supplied. `LOCAL CHECKS MET` requires radius >=3 and no Critical or
High local hit, but it still describes only the supplied reference and model scope.
See the audit for test evidence and remaining limitations.

## Citation

Software citation metadata is provided in [CITATION.cff](CITATION.cff). Zenodo-ready release metadata is provided in [.zenodo.json](.zenodo.json); add a DOI only after Zenodo actually archives a release and mints one.

## Primary references

- Hsu et al. (2013). DNA targeting specificity of RNA-guided Cas9 nucleases.
  [DOI 10.1038/nbt.2647](https://doi.org/10.1038/nbt.2647).
- Doench et al. (2016). Optimized sgRNA design to maximize activity and minimize
  off-target effects of CRISPR-Cas9. [DOI 10.1038/nbt.3437](https://doi.org/10.1038/nbt.3437).
- Gilbert et al. (2014). Genome-scale CRISPR-mediated control of gene repression
  and activation. [DOI 10.1016/j.cell.2014.09.029](https://doi.org/10.1016/j.cell.2014.09.029).
- Horlbeck et al. (2016). Compact and highly active next-generation libraries for
  CRISPR-mediated gene repression and activation.
  [DOI 10.7554/eLife.19760](https://doi.org/10.7554/eLife.19760).
- Poudel et al. (2022). GuideMaker: Software to design CRISPR-Cas guide RNA pools
  in non-model genomes. [DOI 10.1093/gigascience/giac007](https://doi.org/10.1093/gigascience/giac007).
- [GuideScan2 official implementation](https://github.com/pritykinlab/guidescan-cli).

MIT License; CFD weight data retain their separately documented CC0 dedication.
