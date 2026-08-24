# CRISPR Studio user guide

## Before you begin

CRISPR Studio helps you create and compare a **shortlist** of SpCas9 guide RNAs. It supports 20 nt spacers beside an NGG PAM and scans both strands.

It does not prove that a guide is safe or experimentally effective. Always confirm the genomic sequence, exon or regulatory context, whole-genome off-target profile, and experimental controls.

## 1. Open the app

Use the [live Streamlit app](https://crispr-grna-designer-v6mhgxd4o3eqbhgur3anvh.streamlit.app/) or start it locally with:

```bash
streamlit run app.py
```

Use the **Dark mode** switch at the top of the sidebar to choose the appearance you prefer.

## 2. Choose a target input

### Gene lookup

1. Select **Gene lookup**.
2. Enter an official gene symbol, for example `TP53`.
3. Enter the scientific organism name, for example `Homo sapiens`.
4. Choose **Ensembl** or **NCBI**.

The database supplies a representative transcript. If one service is temporarily unavailable or returns no record, try the other service or paste a sequence.

### Paste sequence

1. Select **Paste sequence**.
2. Paste plain DNA, RNA, or FASTA.
3. Make sure the cleaned sequence is at least 50 bp and no more than 50,000 bp.

FASTA headers, spaces, and line numbers are handled automatically. `U` is converted to `T`, and ambiguous IUPAC bases become `N`.

## 3. Choose the design intent

- **Knockout** adds a mild preference for guides in the early 30% of the input sequence. For a real knockout, confirm that the selected site is in a shared coding exon and that indels are likely to disrupt the reading frame.
- **Knockdown / CRISPRi** favors the first 400 bp as a simple 5-prime proxy. Real CRISPRi design should use a verified transcription start site and an appropriate promoter window.

## 4. Adjust the filters

The sidebar controls:

- **Maximum guides:** how many ranked candidates to retain;
- **Minimum activity score:** discard candidates below this heuristic value;
- **Reference-screen mismatches:** maximum substitutions accepted during an optional local-reference screen.

Start with 20 guides, a minimum score of 35, and three mismatches. Tighten the filters after inspecting whether the target contains enough PAM sites.

## 5. Optional local-reference screen

Enable **Add local-reference off-target screen** to compare candidates with a plasmid, amplicon, contig, paralog collection, or small genomic region.

Upload a `.fa`, `.fasta`, `.fna`, or `.txt` file, or paste the sequence. The cleaned reference can contain up to 250,000 bp in the hosted dashboard.

The screen:

- examines NGG-compatible sites on both strands;
- compares their 20 nt spacers with each guide;
- excludes one exact match as the intended target;
- reports additional exact and near matches;
- identifies PAM-proximal seed mismatches;
- provides an internal local-reference specificity score.

This is not a replacement for CRISPOR, CHOPCHOP, Cas-OFFinder, GuideScan, or another genome-indexed workflow.

## 6. Read the report

### Ranked guides

The main table includes:

| Column | Meaning |
|---|---|
| Rank | Ordering after all filters |
| Spacer | 20 nt sequence normally used in guide synthesis |
| PAM | Effective NGG PAM adjacent to the target; do not include it in the spacer oligo unless your protocol says otherwise |
| Strand | Input-sequence strand containing the guide target |
| Start / End | One-based displayed target coordinates |
| GC% | Spacer GC content |
| Score | Explainable 0-100 activity ranking heuristic |
| Specificity | Local-reference score, shown only when a reference was supplied |
| Off-target hits | Retained additional PAM-compatible near matches in that reference |
| Notes | Preferences and issues that deserve review |

### Design landscape

Use the position-versus-score plot to see whether strong candidates cluster in one part of the target. Point size represents GC content and color represents strand. The other plots summarize GC and strand balance.

### Guide details

Select any guide to see:

- spacer and PAM orientation;
- activity, GC, strand, and coordinates;
- validation checklist;
- local sequence context;
- example BbsI cloning oligos;
- score gauge and every non-zero score contribution.

The BbsI overhangs are only an example. Verify the exact vector protocol before ordering oligos.

### Off-target screen

When a reference is present, this tab shows guide-level specificity and lets you inspect every retained hit. An additional exact match is **Critical**. Sites with few mismatches and a conserved PAM-proximal seed are ranked as higher risk.

A blank result only means no near matches were found in the sequence you supplied under the selected mismatch rule.

### Export

- **CSV:** ranked guide table;
- **Excel:** ranked guides, validation, metadata, and reference hits;
- **FASTA:** 20 nt spacers with score/PAM/strand headers;
- **JSON:** structured report for scripts or pipelines.

## 7. How to shortlist candidates

For routine research planning:

1. confirm the spacer exists exactly in the intended genome assembly;
2. reject exon-junction candidates for genomic cutting;
3. prioritize shared coding exons for knockout, unless an isoform-specific design is intentional;
4. prefer acceptable GC, no poly-T motif, and a strong activity rank;
5. run a real whole-genome specificity analysis;
6. check population or strain variants at the spacer and PAM;
7. choose at least two independent guides where practical;
8. include non-targeting and positive controls;
9. validate editing and phenotype with an orthogonal assay.

## 8. Understand the score

The dashboard score is transparent and rule-based. It uses GC range, selected bases near the PAM, PAM context, homopolymers, a simple self-complementarity check, and the design-intent position preference.

It is not the trained Doench Rule Set 2/Azimuth score. Use it to compare candidates generated in the same run; do not treat it as an editing percentage.

## 9. Common problems

### Gene not found

- verify the official gene symbol;
- use the scientific organism name;
- try the other database;
- paste a verified sequence.

### No guides pass

- lower the minimum score;
- use a longer sequence;
- verify that the input contains NGG PAMs;
- check whether many bases were converted to `N`.

### Database timeout

Public services can throttle or briefly fail. Wait and retry, switch database, or use pasted sequence mode.

### Too many reference hits

Reduce the mismatch limit, use a more relevant reference, or move to a genome-indexed off-target tool. Repeated sequences and paralogs can legitimately produce many sites.

## 10. Research-use reminder

CRISPR Studio is an educational design aid. Do not use its results as clinical guidance or as the only evidence for guide safety, specificity, or efficacy.

