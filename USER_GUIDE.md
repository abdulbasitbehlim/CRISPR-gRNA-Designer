# CRISPR Studio v3.2.0 — user guide

## Before you begin

CRISPR Studio helps create and compare a shortlist of **SpCas9 20 nt + NGG** guide RNAs. It supports knockout and TSS-aware CRISPRi workflows, optional specificity analysis and online sequence retrieval.

It does not prove that a guide is safe or experimentally effective. Confirm genomic context, off-target profile and experimental controls before ordering guides.

## 1. Open the app

Use the live Streamlit app or run locally:

```bash
streamlit run app.py
```

Use the **Dark mode** switch in the sidebar. Both themes are designed for readable inputs, dropdowns, cards and plots.

## 2. Choose a target input

### Gene lookup

1. Select **Gene lookup**.
2. Enter a gene symbol such as `TP53`.
3. Enter an organism such as `Homo sapiens`.
4. Choose Ensembl or NCBI for knockout candidate discovery.

For **CRISPRi repression**, gene lookup automatically uses Ensembl genomic TSS annotation rather than the transcript 5′ end.

### Accession ID

1. Select **Accession ID**.
2. Enter an accession/stable ID.
3. Choose **Auto**, **NCBI Nucleotide**, or **Ensembl**.

Examples:

- `NM_000546.6`
- `NC_000017.11`
- `ENST00000269305`
- `ENSG00000141510`

Auto mode routes `ENS...` IDs to Ensembl and other nucleotide accessions to NCBI. The app reports the resolved accession, source database, record type and sequence length.

Accession-only CRISPRi is not inferred because arbitrary sequence records may not define a biologically reliable TSS. Use **Gene lookup** for TSS-aware CRISPRi or **Paste sequence** with an explicit TSS.

### Paste sequence

Paste plain DNA, RNA or FASTA. The app removes FASTA headers/whitespace, converts `U` to `T`, and converts supported ambiguous IUPAC bases to `N`.

For pasted-sequence CRISPRi, provide the **1-based TSS position** and use genomic DNA in 5′→3′ transcriptional orientation.

## 3. Choose the design intent

### Knockout

Scans both strands for 20 nt spacers adjacent to NGG PAMs and applies a mild early-target preference in the transparent heuristic. Confirm final guides in the intended genomic coding exon and assembly.

### CRISPRi repression

For gene lookup, the app resolves the canonical Ensembl transcript and true genomic TSS. It keeps guides from **−50 to +300 bp** relative to TSS and prioritizes **+50 to +100 bp**.

## 4. Adjust design settings

Sidebar controls include:

- **Maximum guides**
- **Minimum heuristic score**
- **Off-target mismatches**

A practical starting point is 20 guides, heuristic threshold 35 and three mismatches.

## 5. Specificity analysis

### None

No off-target specificity score is calculated. The dashboard clearly reports **Not screened** rather than treating this as a failed result.

### Local reference

Upload or paste a FASTA/reference sequence. Multi-FASTA contigs are preserved independently. The app scans NGG-compatible sites on both strands and calculates native MIT/Hsu specificity; CFD is shown when its optional provider is available.

The hosted local-reference limit is 5,000,000 bp total.

### Whole genome (GuideScan2)

Requires an installed `guidescan` executable and matching prebuilt genome index on the machine running the app. This is not automatically available on a normal hosted Streamlit instance.

## 6. Read the report

### Ranked guides

The table includes spacer, PAM, strand, coordinates, GC, heuristic score, optional Doench Rule Set 2, MIT specificity, CFD specificity and notes. CRISPRi rows additionally include TSS/genomic provenance.

### Design landscape

For knockout, the plot shows guide position versus heuristic score. For CRISPRi, the x-axis is TSS distance and the preferred +50 to +100 bp region is highlighted.

### Guide details

The guide detail dashboard shows:

- synthesis-oriented spacer + PAM;
- heuristic, Doench RS2, GC, MIT and CFD cards;
- TSS/genomic metadata for CRISPRi;
- readable validation cards instead of raw JSON;
- heuristic contribution table;
- example BbsI cloning oligos.

### Off-target screen

When specificity analysis was run, inspect additional exact and near matches. Additional perfect copies are critical. Local-reference results apply only to the sequence supplied.

### Export

- **CSV** — ranked table
- **FASTA** — synthesis-ready spacers
- **JSON** — structured metadata and guide results

## 7. Understand the scores

- **Heuristic**: transparent sequence-ranking score; not a calibrated editing probability.
- **Doench Rule Set 2**: on-target activity model when the optional compatible provider is installed.
- **MIT/Hsu**: native off-target specificity score.
- **CFD**: optional off-target cleavage-risk score.

Do not average these metrics mentally into one number; they answer different questions.

## 8. Recommended shortlisting workflow

1. Confirm the spacer exists in the intended genome assembly.
2. Confirm exon/regulatory context and intended transcript.
3. Reject exon-junction-only knockout candidates.
4. Review GC, poly-T and other validation flags.
5. Run genome-wide specificity analysis for final candidates.
6. Check relevant population/strain variants when appropriate.
7. Select multiple independent guides where practical.
8. Validate experimentally with suitable positive, negative and non-targeting controls.

## 9. Common problems

### Gene not found

Check the official symbol and organism name, switch database, use an accession ID or paste a verified sequence.

### Accession not found

Verify the accession/version and choose the correct database explicitly. Public APIs may temporarily throttle or fail.

### No guides pass

Lower the heuristic threshold, use a longer target or verify that the sequence contains NGG PAMs and few ambiguous bases.

### Doench / CFD shows N/A

Those models are optional providers. The app does not substitute the internal heuristic and call it Doench.

### MIT shows Not screened

Choose Local reference or Whole genome mode. A guide that was not screened is not the same as a guide that failed specificity.

## 10. Research-use reminder

CRISPR Studio is a research and educational candidate-shortlisting tool. It does not replace genome-aware design review, biological judgment or experimental validation.
