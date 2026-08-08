# Architecture & Design Notes

This document describes the internal design of the CRISPR gRNA Designer:
why it was built this way, the scientific rules it implements, and where
its boundaries are. It's aimed at contributors and anyone who wants to
extend or audit the tool.

---

## 1. Motivation

Designing CRISPR guide RNAs usually means: look up the gene sequence in
NCBI or Ensembl, scan it by hand (or with a separate tool) for PAM sites,
and cross-reference scoring rules from the literature. Excellent web
tools already do this well (see below), but there wasn't a small,
self-contained, Python-native tool that goes straight from a **gene
symbol** to a **ranked guide list** with an explicit Knockout/Knockdown
mode, runnable locally with no account or server setup.

That gap is what this project targets — a lightweight dashboard suitable
for teaching, quick prototyping, or embedding inside a larger pipeline.

---

## 2. Where this fits among existing tools

| Tool | Type | Strengths | Notes |
|------|------|-----------|-------|
| [CRISPOR](https://crispor.tefor.net) | Web | Off-target search, many scoring models | Gold standard for single guides |
| [CHOPCHOP](https://chopchop.cbu.uib.no) | Web | KO / KD / activation modes, many organisms | Very user-friendly |
| crisprVerse (Bioconductor) | R packages | Most comprehensive modern ecosystem | Requires R |
| GuideScan2 | CLI + web | Genome-wide, custom genomes | C++ heavy |
| CRISPRdirect | Web + API | Fast | Good API |
| Benchling | Commercial web | Sequence import from NCBI/Ensembl | Closed source |
| GuideMaker | Python + web | Non-model organisms | |
| CRISPRware | Python | Context-aware libraries | |

This project doesn't try to replace any of the above, particularly for
off-target analysis. It fills a narrower niche: a single-file-installable
Streamlit dashboard with a pure-Python core, transparent scoring, and an
explicit Knockout-vs-Knockdown bias.

---

## 3. Scientific basis

**Core mechanism**
- Jinek et al., *Science* 2012 — Cas9 is guided by a 20-nt spacer plus a
  PAM (NGG for SpCas9).

**On-target activity heuristics implemented here**
- Doench et al., *Nature Biotechnology* 2014 & 2016 (Rule Set 1 / Rule
  Set 2) — position-specific nucleotide preferences, GC content 40–70 %,
  G preferred at position 20, avoid T at position 20, avoid poly-G/T runs.
- Moreno-Mateos et al., *Nature Methods* 2015 (CRISPRscan) — similar
  positional features validated in zebrafish / in-vitro T7 systems.
- Hsu et al., *Nature Biotechnology* 2013 — mismatch tolerance and the
  importance of the PAM-proximal seed region.

**Specificity / off-target**
- Mismatches in the seed region (the ~10 nt closest to the PAM) are more
  disruptive to binding than mismatches further away (Hsu et al. 2013).
- This tool deliberately does **not** implement a genome-wide off-target
  search — that requires a full genome index (Bowtie/BWA) which would
  make the tool heavy and installation-unfriendly. Instead it flags this
  clearly and points users to CRISPOR / Cas-OFFinder for that step.

**Application-specific targeting**
- **Knockout (NHEJ):** guides are biased toward the early part of the
  fetched sequence, as a proxy for early exons — a frameshift there is
  more likely to trigger nonsense-mediated decay and fully disrupt the
  protein.
- **Knockdown (CRISPRi-style):** guides are biased toward the first
  ~400 bp, as a proxy for proximity to the transcription start site,
  where dCas9/KRAB-dCas9 binding is most effective at blocking
  transcription.

**Practical cloning notes**
- U6-driven guides should avoid long T-stretches (RNA Pol III
  terminator), which is why guides ending in "TTTT" or similar are
  penalized.
- The app includes an example BbsI/U6-style oligo (`CACC…` /
  `AAAC…`) for the top-ranked guide, matching the common pX330 /
  pSpCas9(BB) cloning scheme.

---

## 4. Architecture

```
crispr-grna-designer/
├── app.py                  # Streamlit dashboard (UI layer)
├── grna_designer.py        # Core library: fetch, design, score, validate
├── test_grna_designer.py   # Unit tests
├── requirements.txt        # Runtime dependencies
├── requirements-dev.txt    # + pytest, for running the test suite
├── USER_GUIDE.md           # Plain-language walkthrough
├── ARCHITECTURE.md         # This file
└── .github/workflows/
    └── tests.yml           # CI: runs pytest on every push/PR
```

### Key functions in `grna_designer.py`

| Function | Role |
|----------|------|
| `fetch_gene_ncbi(gene, organism)` | Entrez Gene → Nuccore → preferred RefSeq transcript |
| `fetch_gene_ensembl(gene, species)` | Ensembl REST `/lookup/symbol` + `/sequence/id` |
| `fetch_sequence(...)` | Unified dispatcher between the two sources |
| `design_guides(sequence, application, ...)` | Scans both strands for NGG/CCN, scores, ranks |
| `_doench_like_score(spacer, pam)` | Transparent heuristic score, 0–100 |
| `validate_guide(guide)` | Boolean checklist for downstream QC |
| `design_from_gene(...)` | One-call, end-to-end convenience API |

### Scoring features currently used
- GC content window (40–70 % optimal)
- Position 20 = G bonus / T penalty
- Position 19 A/G preference
- PAM identity (CGG slightly preferred)
- Homopolymer penalties (e.g. GGGG, TTTT)
- Simple self-complementarity penalty
- Application positional bias (early region for KO, 5′-proximal for KD)

---

## 5. Validation strategy

1. **Automated checks** (`validate_guide`): exactly 20 nt, PAM ends in
   GG, GC between 30–80 %, no ≥5-nt homopolymer, score ≥ 40.
2. **Biological sanity check:** for KO, confirm the guide falls inside a
   coding exon (genome browser / Ensembl gene view); for KD, confirm
   proximity to the TSS.
3. **Cross-check against a gold-standard tool:** paste the same sequence
   into CRISPOR or CHOPCHOP and compare the ranking of the top 3–5
   guides, including off-target counts (which this tool does not compute).
4. **Experimental confirmation:** T7E1 / ICE / TIDE / NGS after
   transfection; Western blot or a functional assay for KO; RT-qPCR for KD.
5. **Suggested regression test for CI:** run a sequence containing a
   published, well-characterized guide (e.g. an EMX1 or VEGFA guide from
   the literature) and confirm it scores near the top of the ranking.

---

## 6. Known limitations

- No genome-wide off-target search (by design, to keep the tool
  dependency-light — use CRISPOR / Cas-OFFinder for that).
- Transcript selection is heuristic (prefers RefSeq mRNA / Ensembl
  canonical transcript); alternative isoforms may be more appropriate
  for some genes.
- The score is a transparent, linear heuristic — not the full Doench
  Rule Set 2 gradient-boosted model, nor a deep-learning model such as
  DeepHF or CRISPRon.
- NCBI Entrez can rate-limit repeated requests; the Ensembl REST API is
  generally more robust for higher-throughput use.

---

## 7. Possible extensions

1. Integrate a proper ML on-target model (e.g. a Rule Set 2 / Azimuth
   port) alongside the current transparent heuristic.
2. Optional Bowtie2 / Cas-OFFinder off-target module when a genome FASTA
   is supplied locally.
3. Support additional systems: Cas12a (TTTV PAM), Cas13 (RNA targeting),
   base editors, and prime-editing pegRNAs.
4. Batch mode: upload a list of gene symbols and export a full guide
   library to Excel.
5. A CI regression test against a small set of published, validated
   guides (see Section 5).
6. Optional Docker image for one-command deployment.

---

## 8. License & citation

Released under the [MIT License](LICENSE). If you use or extend this
tool for published research, please cite the original scoring papers
(Doench 2016, Hsu 2013, Moreno-Mateos 2015) and the public databases used
(NCBI, Ensembl).
