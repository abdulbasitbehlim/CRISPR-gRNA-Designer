# CRISPR Studio v3.3.0 architecture

| Module | Responsibility |
|---|---|
| app.py | Inputs, workflow routing, explicit intended region, reports and exports |
| sequence_io.py | Strict single/multi-FASTA parsing and DNA normalization |
| knockout.py | NCBI genomic retrieval, CDS/isoform annotations and cut filtering |
| grna_designer.py | Both-strand NGG discovery, coordinates, heuristic, MIT, local index and GuideScan2 adapter |
| crispri.py | Ensembl TSS context, strand mapping, position filtering and ranking |
| accession_lookup.py | Version-checked nucleotide retrieval; RNA/protein safeguards |
| models.py | Bundled CFD and optional GuideMaker Rule Set 2 provider |
| network.py | Bounded transient retries and process-local NCBI pacing |
| reporting.py | Strict JSON, settings and sequence/reference fingerprints |

## Coordinate contract

Internal intervals are zero-based, half-open. A GuideRNA start/end interval
covers the 23-base spacer-plus-PAM locus at its left edge on the input sequence,
including on the minus strand. Public tables use 1-based inclusive positions.
The plus-strand spacer spans start to end-3; minus spans start+3 to end. Predicted
cut boundaries are start+17 and start+6, respectively. A cut boundary identifies
the number of bases before the break; it is not another nucleotide coordinate.

TargetLocus(contig, start, strand) uses a zero-based 23-base left edge. Legacy
three-element tuples remain accepted. A declared region is checked against the
entire input sequence in the requested orientation before guide loci are mapped.
An exact sequence match somewhere else is never proof of the intended locus.

## Knockout pipeline

NCBI Gene search must resolve one exact symbol and one genomic mapping. Fetch a
forward-orientation chromosome slice and complete local GenBank CDS annotations.
Choose an explicitly requested isoform or prefer NP_ proteins, then longest CDS
and accession for deterministic selection. Scan contiguous genomic DNA and retain
cuts strictly inside a selected CDS segment. Introns are present during scanning;
CDS/exon segments are never concatenated. Annotation reports compatible coding
isoforms and strand-aware CDS progress without claiming a functional knockout.

Ensembl gene-symbol knockout is rejected. Ensembl accession mode requests
`type=genomic`, verifies interval length and reports that no CDS filter was applied.
NCBI accession mode requires genomic biomolecule metadata; RNA/cDNA is rejected.
Legacy transcript retrieval helpers exist for callers needing sequences, but no
application or design_from_gene knockout path uses them.

## Local screening contract

ReferenceIndex holds separate normalized contigs and reusable NGG site vectors.
Mismatch enumeration is chunked. Counts and MIT/CFD risk sums are accumulated over
all qualifying sites; only detail storage is bounded. Exact-match count includes
a declared intended site, whereas total_hits excludes it after verification.
No-match, undeclared and verified-locus states are distinct. A missing declared
site raises an error. Intended mapping and full scoring precede top-N selection.

GuideRNA retains total count, details, status, intended locus, reference
ambiguity, searched mismatch radius, complete risk-category counts and maximum
per-site MIT/CFD risk. These summaries are accumulated before display truncation.
Validation requires sequence checks, a verified intended locus, exactly one
reference exact match, no reference ambiguity, local MIT >=50, a search radius
of at least three mismatches and zero non-intended Critical/High hits. Moderate
and Low hits remain governed by the aggregate evidence and user review. This
heuristic status is named LOCAL CHECKS MET, never experimental validation.

## Models and ranking

CFD uses the product of published mismatch and PAM weights; the shipped JSON
comes from the documented GuideMaker revision. The local enumerator remains NGG
only. Optional RS2 expects a 30-mer (4 upstream bases, 20-base spacer, PAM, 3
downstream bases), oriented to the guide. Failures/missing context return N/A with
provider warnings where appropriate. Regression output is displayed at 100x,
without calling it probability or clipping it. It is not mixed with the heuristic.

Knockout ranking is heuristic then local specificity. CRISPRi first filters to the
TSS window, then ranks position band, local specificity, heuristic and proximity
to +75 bp. RS2 is not calculated for CRISPRi. GuideScan2 output stays separate.

## Reproducibility and tests

Exports preserve UTC time, version, parameter settings, SHA-256 fingerprints,
reference contig identities, genomic/TSS provenance, full hit counts and detail
truncation. Nonfinite JSON values become null. Tests include the four reported
bugs, API routing, strand conversion, scoring weights, accession errors, saved
ACTB/GAPDH contexts and Streamlit workflows. CI runs the same offline suite.

External services can change. Network retries handle transient failures but do
not validate annotation correctness. Live Ensembl and GuideScan2 genome-index
behavior are not established by mock tests. See SCIENTIFIC_AUDIT.md.
