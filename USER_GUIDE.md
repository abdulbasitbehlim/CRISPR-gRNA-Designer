# CRISPR Studio v3.3.0 user guide

## What the software helps you do

CRISPR Studio makes a shortlist of possible SpCas9 guides. A spacer is the 20-base
sequence directing Cas9; an NGG PAM is a three-base motif beside the genomic
site. The software does not demonstrate that a guide edits cells successfully.

## Genomic knockout design

1. Choose Gene lookup, enter an exact gene symbol and scientific organism name.
2. Select NCBI and Knockout. For a quick human example, use ACTB and Homo sapiens.
3. Optionally enter a RefSeq protein or compatible transcript accession for the
   intended isoform. Leaving it blank uses the selection rule shown in the result.
4. Choose the minimum heuristic and number of returned guides.
5. Add a local reference if available, or leave specificity as not screened.
6. Run the analysis. Check the genomic accession, selected protein, CDS position,
   spacer sequence and predicted cut coordinates.
7. Review specificity evidence, biological context and experimental controls.

The software searches actual genomic DNA with introns present. A candidate is
kept only when its predicted cut lies inside a coding segment of the selected
isoform. A guide's complete 23-base site can legitimately extend beyond that
coding segment when it exists contiguously in genomic DNA.

Ensembl gene-symbol knockout is disabled. Use NCBI for this annotated workflow.

## Other input modes

NCBI accession mode accepts records verified as genomic DNA. NM_/XM_/NR_/XR_
transcripts and other RNA/cDNA records are rejected. Ensembl gene/transcript/exon
IDs retrieve the genomic interval, never the spliced transcript. Protein IDs are
rejected. The 2 Mb limit means whole human chromosome accessions are too large.
Use gene lookup or a smaller genomic region.

Paste sequence accepts one genomic DNA record. Multiple target records are an
error, not joined sequences. The app cannot prove the provenance of pasted DNA;
do not paste a spliced cDNA for genomic editing. Coding annotation is not applied
in accession/paste modes, so a candidate can lie in an intron or noncoding region.

## Declaring your intended target in a local reference

Upload/paste the correct reference, with a unique ID for every FASTA record.
If the target input appears inside reference contig chr1 beginning at base 101,
enter chr1, 101 and +. For a reverse-complement input region, enter its leftmost
reference coordinate and -. The whole input region must match at that location.
The program then maps each candidate to its own intended locus before ranking.

If the reference consists of exactly your target under >reference, enter
reference, 1 and +. This is useful for testing or a construct, but it only screens
that supplied reference. It does not examine other genome sequences.

Leave the intended-record field blank when you do not know the location. Every
exact copy is then retained, including the possible intended copy. A high score
with no exact match requires investigation; it can indicate the wrong reference.

## Reading scores and coordinates

| Field | Meaning |
|---|---|
| Heuristic | Transparent sequence preference; not a validated success percentage |
| Doench RS2 | Optional nuclease-activity prediction; N/A if provider/context unavailable |
| MIT specificity | Reference-limited aggregation of MIT pair scores |
| CFD specificity | Reference-limited aggregation of published CFD pair scores |
| Off-target hits | All qualifying hits after an explicitly verified intended exclusion |
| Exact reference matches | All exact NGG matches, including a verified intended match |
| Start and End | 1-based complete spacer-plus-PAM interval on input DNA |
| Spacer start and end | 1-based spacer interval without PAM |
| Cut after base | Predicted between-base cleavage position |

Details are capped at 250 rows per guide by default; all hits still contribute to
counts and scores. LOCAL CHECKS MET means the stated local checks passed. REVIEW
also covers unscreened, unverified-locus and ambiguous-reference cases.

## CRISPRi repression

Gene lookup uses Ensembl for human/mouse dCas9-KRAB design, regardless of the
sequence-database radio selection. Confirm the canonical transcript or provide
an explicit Ensembl transcript ID; the tool will not silently choose an arbitrary
transcript if no canonical transcript is supplied. The annotation TSS may differ
from the active TSS in your cell type.

The app keeps spacer midpoints -50 to +300 bases relative to TSS and prioritizes
+50 to +100. These are placement heuristics, not universal biological rules.
Custom input must be genomic DNA in transcriptional orientation with a 1-based
TSS. The app does not apply these rules to gene-based plant/bacterial designs.
It does not use the RS2 nuclease model to predict repression.

## Exports and external analysis

Download CSV for the guide table, FASTA for spacer sequences and JSON for settings,
reference fingerprints, provenance and local evidence. Keep the original reference
and exact genome assembly alongside the JSON; a hash identifies a sequence but
cannot reconstruct it. FASTA is not a complete expression construct.

GuideScan2 requires an installed executable and a compatible genome index on the
server. Its raw per-hit table uses native scoring and is not folded into local
validation or ranking. This release tests the adapter interface with simulated
output; users must validate their actual external installation/index.

## Troubleshooting

- Duplicate FASTA identifier: give every contig a unique name; do not delete a
  relevant contig merely to silence the error.
- Declared target mismatch: check assembly, contig naming, region coordinate and
  orientation. Do not exclude a different exact match to force a better score.
- No complete CDS: confirm gene symbol, organism and selected coding isoform.
- RNA/cDNA accession rejected: use NCBI Gene lookup or genomic input.
- Rule Set 2 N/A: optional provider missing, invalid 30-mer context or provider
  failure. The heuristic remains clearly separate.
- Limits exceeded: reduce the analysis region or use an indexed genome tool.
- Network failure: retry later or use verified saved genomic input. Provider
  outages must not be interpreted as a biological absence of candidate guides.
