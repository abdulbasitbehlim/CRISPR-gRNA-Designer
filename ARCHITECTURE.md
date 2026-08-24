# Architecture and scientific design notes

## 1. Scope

CRISPR Studio is a lightweight SpCas9 guide-candidate workbench. Its job is to make four operations transparent and easy to audit:

1. obtain or accept a target sequence;
2. find 20 nt spacers beside NGG PAMs on both strands;
3. rank candidates with an explainable activity heuristic;
4. optionally screen those candidates against a small user-supplied reference.

It is deliberately not a genome browser, clinical decision system, trained Rule Set 2 implementation, or genome-indexed off-target aligner.

## 2. Components

```text
┌────────────────────────────────────────────────────────────┐
│ app.py                                                     │
│ Streamlit state, form validation, dashboard, plots, export │
└────────────────────────────┬───────────────────────────────┘
                             │ pure Python calls
┌────────────────────────────▼───────────────────────────────┐
│ grna_designer.py                                           │
│ clean → fetch → scan → score → validate → reference screen │
└───────────────┬───────────────────────────────┬────────────┘
                │                               │
       ┌────────▼────────┐             ┌────────▼────────┐
       │ NCBI Entrez     │             │ Ensembl REST   │
       │ representative │             │ canonical cDNA │
       │ transcript      │             │ where present  │
       └─────────────────┘             └─────────────────┘
```

The Streamlit layer never reimplements guide logic. This keeps the core importable, testable without a browser, and suitable for a later CLI or API.

## 3. Core data models

### `GuideRNA`

Stores the spacer, effective NGG PAM, strand, zero-based leftmost input-sequence coordinate, end coordinate, GC content, activity score, notes, application, and optional reference-screen results.

`GuideRNA.to_dict()` converts internal coordinates to a one-based display start and creates a stable tabular export schema.

### `OffTargetHit`

Stores a single PAM-compatible reference near-match:

- leftmost zero-based reference coordinate;
- strand;
- candidate spacer and effective PAM;
- total mismatches;
- one-based mismatch positions in guide 5-prime to 3-prime orientation;
- mismatch count in spacer positions 13-20;
- risk tier.

### `OffTargetReport`

Contains the retained hits, local-reference specificity score, PAM-site count, and whether one exact match was excluded as the presumed intended target.

## 4. Sequence normalization

`clean_dna_sequence()` accepts plain DNA, RNA, or FASTA text.

- FASTA header lines are removed.
- Whitespace and numeric line labels are removed.
- `U` is converted to `T`.
- standard ambiguous IUPAC symbols are converted to `N`.
- unexpected punctuation raises `ValueError` instead of being silently deleted.

Guide discovery skips any 20 nt spacer containing `N`, preserving predictable coordinates and avoiding false certainty at ambiguous bases.

## 5. Sequence retrieval

### NCBI

The NCBI path searches Gene with `gene symbol + organism`, follows nucleotide links, prefers RefSeq RNA records, and selects a record described as mRNA/transcript or containing a CDS feature. Requests retry up to three times with a short backoff.

### Ensembl

The Ensembl path maps common organism names to species slugs, expands a gene-symbol lookup, chooses a canonical transcript where marked, and downloads cDNA.

### Important consequence

Both lookup paths return a representative transcript/cDNA, not an annotated genomic exon model. A 20 nt candidate can cross an exon-exon junction and therefore fail to exist in genomic DNA. Final candidates must be mapped to the intended genome assembly, transcript isoform, and coding exon before ordering.

## 6. SpCas9 discovery

Forward candidates use:

```text
5′ — [20 nt spacer] [NGG PAM] — 3′
```

Reverse candidates are recognized as `CCN` in the input sequence. The downstream 20 nt are reverse-complemented, and the `CCN` is reported as its effective `NGG` PAM so every exported spacer is in the synthesis-ready guide orientation.

The scanner uses overlapping regular-expression lookaheads, so nearby PAM sites are not skipped.

## 7. Activity ranking

`score_breakdown()` starts at 50 and records every adjustment. The final value is clamped to 0-100.

| Feature | Adjustment |
|---|---:|
| GC 40-70% | +15 |
| GC 30-40% or 70-80% | +5 |
| GC outside those ranges | -15 |
| G at spacer position 20 | +10 |
| C at spacer position 20 | -8 |
| T at spacer position 20 | -12 |
| A/G at position 19 | +5 |
| C at position 16 | +3 |
| C at position 18 | +3 |
| CGG PAM | +5 |
| TGG PAM | -5 |
| GGG or TTT in spacer | -8 |
| GGGG or TTTT in spacer | additional -20 |
| high simple reverse-complement match | -10 |
| knockout in early 30% | +8 |
| knockdown/CRISPRi in first 400 bp | +15 |
| knockdown/CRISPRi after first 400 bp | -5 |

This is a literature-inspired linear ranking heuristic. It is not the trained Doench Rule Set 2/Azimuth model, and its values should not be interpreted as calibrated editing probabilities.

## 8. Validation checklist

`validate_guide()` reports:

- spacer length is 20 nt;
- PAM ends in GG;
- GC lies in preferred and acceptable windows;
- no homopolymer of five or more bases;
- activity score is at least 40;
- no `TTTT` polymerase III termination motif;
- whether reference screening ran;
- whether local-reference specificity is at least 50.

The overall pass requires the structural checks, acceptable GC, no long homopolymer/poly-T motif, score threshold, and an acceptable specificity result when available.

## 9. Local-reference similarity screening

`analyze_offtargets()` performs a bounded, PAM-aware screen:

1. normalize the reference;
2. enumerate forward NGG and reverse CCN sites;
3. orient every candidate as a guide spacer plus effective NGG PAM;
4. calculate spacer Hamming distance without indels;
5. retain sites within the configured mismatch limit;
6. exclude the first exact site as the presumed intended target;
7. rank additional hits by mismatch count, seed mismatches, and coordinate.

Positions 13-20 are treated as the PAM-proximal seed region. Risk tiers are intentionally simple:

- **Critical:** an additional exact site;
- **High:** one mismatch, or two mismatches with no more than one seed mismatch;
- **Moderate:** up to three mismatches with no more than one seed mismatch;
- **Low:** retained sites with more seed disruption.

The specificity calculation is:

```text
specificity = 100 / (1 + weighted_hit_risk / 20)
```

Base weights for 0, 1, 2, and 3 mismatches are 100, 25, 8, and 2. A fully conserved seed multiplies risk by 1.5; one seed mismatch uses 1.0; two or more seed mismatches use 0.5.

This score is designed for relative sorting inside a supplied reference. It is not CFD, MIT, cutting-frequency, or genome-wide specificity.

### Complexity and limits

Reference PAM sites are enumerated for each screened guide. For `G` guides and `P` PAM sites, runtime is approximately `O(G × P × 20)`. The dashboard limits custom targets to 50,000 bp and reference input to 250,000 bp so a free hosted instance remains responsive. Each report retains at most 250 hits.

## 10. Dashboard state and caching

- Public-database responses are cached for one hour.
- Pure guide-design results are cached by all input parameters.
- Completed analysis is placed in `st.session_state` so changing result tabs does not discard the report.
- Dark/light mode is CSS-variable based and persists in the Streamlit session.
- Plotly charts use a matching dark or light template.

No target sequence is intentionally written to server disk. Export files are generated in memory.

## 11. Export design

- CSV: ranked guide table.
- Excel: ranked guides, validation, metadata, and reference hits when present.
- FASTA: synthesis-ready 20 nt spacer records with score/PAM/strand metadata.
- JSON: metadata plus the ranked table for downstream automation.

The example BbsI cloning oligos are explicitly labeled vector-specific because overhangs differ among plasmid systems.

## 12. Error handling and deployment

- Network requests use timeouts.
- NCBI requests retry with backoff.
- Unknown source/application values raise clear errors.
- Unsupported nuclease/PAM options fail explicitly rather than appearing to work.
- Upload size is limited by `.streamlit/config.toml`, while cleaned biological sequence lengths are checked in `app.py`.
- CORS and XSRF protection remain enabled for deployment.

Production deployments should set `NCBI_EMAIL` to a monitored contact address.

## 13. Test strategy

The offline suite covers:

- FASTA/RNA/IUPAC normalization and invalid input;
- GC and homopolymer helpers;
- scoring range and application clamp;
- forward/reverse PAM discovery and ranking;
- guide validation;
- intended-target exclusion;
- exact duplicate and mismatch detection;
- attachment of screening results to guide objects.

Streamlit's application test runner is used during release checks to load the dashboard and execute a complete pasted-sequence workflow without calling an external database.

## 14. Recommended future work

1. Fetch exon/CDS annotations and prevent exon-junction candidates.
2. Integrate a validated Rule Set 2 or DeepSpCas9 implementation with model/version metadata.
3. Add genome-build-aware Bowtie2 or Cas-OFFinder execution for real genome-wide specificity.
4. Add CRISPRi TSS annotation and strand-aware activity windows.
5. Add alternate nucleases, PAMs, base-editor windows, and prime-editing workflows as separate validated models.
6. Add batch design with provenance-aware Excel output.

