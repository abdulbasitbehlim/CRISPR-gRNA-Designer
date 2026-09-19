# CRISPR Studio v3.3.0 scientific and technical audit

Review date: 19 September 2026. Original baseline: commit
`8235cd88e9cf8f6b1e943441166bf5a51105fda2`. The original 24 tests passed during the
initial review; the expanded final suite contains 82 passing tests. Passing tests
establish the behavior covered below, not freedom from every possible software
error or experimental validation of the guides.

## Four demonstrated failures and their resolution

| Failure | Consequence | Implemented resolution | Regression evidence |
|---|---|---|---|
| First exact match removed automatically | An unknown exact site could vanish and yield a misleading score of 100 | Explicit TargetLocus with contig, zero-based full-site left edge and strand; exact site must be verified | One undeclared exact site remains one hit, MIT/CFD 50; wrong locus raises an error |
| Hit-display cap changed specificity | Identical biology produced different scores | Accumulate counts and risk sums across all qualifying hits; cap details only | Twenty copies: 4.76 with none excluded, 5.00 with one declared exclusion, invariant at caps 0/1/10/100; 301-copy stress case also invariant |
| Spliced transcript scanned for knockout | Artificial exon-junction spacer/PAM could be ranked | NCBI genomic retrieval and CDS-cut filter; no exon concatenation; unsafe gene-symbol Ensembl route disabled; RNA accession rejection; Ensembl accession uses genomic interval | A20 followed by exon-2 AGG is found in the deliberately spliced control but absent when the intervening intron is retained |
| Duplicate FASTA identifiers overwritten | A contig containing real matches could be lost | Reject duplicate identifiers, empty IDs/records, malformed preamble and multiple target records | Parser and Streamlit regression tests; previous results are cleared when new input fails |

## Additional changes

1. Both spacer and PAM must be unambiguous A/C/G/T. Ambiguous reference windows
   are skipped and reported; a reference containing ambiguity cannot meet the
   app's local quality status.
2. Intended regions may occur inside larger contigs and in reverse orientation.
   Full-region sequence identity is checked before mapping each guide locus.
   Exclusion and scoring run before top-N selection, including CRISPRi ranking.
3. Unscreened or unverified results cannot receive LOCAL CHECKS MET. That label
   requires the sequence checks, a verified locus, one exact reference match,
   no reference ambiguity, local MIT >=50, a search through at least three
   mismatches and no non-intended Critical or High hit. The aggregate threshold
   cannot override high-risk per-site evidence.
4. Reports preserve complete Critical, High, Moderate and Low counts and maximum
   per-site MIT/CFD risk before hit-detail truncation. No universal CFD cutoff was
   invented; the empirical per-site evidence is exported for review.
5. CFD weights are bundled from a pinned GuideMaker source revision, with data
   license/provenance. Exact-match and known mismatch/PAM weights are tested.
6. The optional RS2 interface uses the provider's sequence-array call, checks a
   correctly oriented 30-mer, reports availability/failure and keeps its scale
   separate. Missing edge context is N/A. No repression predictions use RS2.
7. NCBI records must match requested versions. Ensembl versions and genomic
   interval lengths are checked; protein translation IDs are rejected before
   normalization can misinterpret a protein as DNA.
8. CRISPRi filters the TSS window before output truncation; reports half-base
   spacer-midpoint distances correctly; does not silently substitute an arbitrary
   transcript when a canonical transcript is absent. Gene-based use is restricted
   to the human/mouse dCas9-KRAB placement scope.
9. A reusable, chunked local NGG index and explicit workload limits replace
   unbounded repeated searches. Limits cause errors, never partial biological
   scores presented as complete.
10. Public-data retrieval has transient-error retries and process-local NCBI
   pacing. This is not a distributed quota controller for many server replicas.
11. JSON records sequence/reference SHA-256, settings, coordinates, intended
    loci, provenance, full hit totals and capped evidence. NaN becomes null.
12. GuideScan2 input now supplies a numeric coordinate placeholder, as required
    by the upstream CSV reader. Missing output, schema errors and missing guide
    IDs raise errors. Multiple hit rows per guide remain valid. Native output
    remains separate because its score convention differs from local scoring.

## Genomic knockout interpretation

The selected NCBI chromosome slice is forward-oriented and contiguous. Complete
CDS feature segments guide filtering; the predicted cut must fall strictly inside
one segment of the selected isoform. Scanning a real genomic sequence also allows
real exon/intron-border sites when the cut is coding. It does not require the
entire spacer plus PAM to fit within an exon, and does not invent junction sites.

Automatic isoform selection prefers an NP_ RefSeq protein then longest CDS and
accession as a tie-breaker. This is a deterministic representative, not a claim
of MANE/canonical or cell-specific expression. Users can provide a compatible
protein/transcript accession. Tables report coding isoforms at the predicted cut;
this is annotation overlap, not experimental evidence of editing every isoform.

## Verification performed

- `python -m pytest -q`: **82 passed** in the final local Python 3.12 run. No
  optional genome-index or activity-provider integrations were needed for this run.
- Seven Streamlit tests cover loading, pasted design, short-input errors,
  duplicate references/stale result removal, intended-region screening,
  genomic gene routing and blocked Ensembl transcript knockout.
- Live NCBI ACTB retrieval after the fixes: NC_000007.14:5527144-5530605,
  3,462 bases; selected NP_001092.1; 20 returned candidates, all coding cuts.
- Saved same-session GAPDH context: NC_000012.12:6534513-6538375,
  3,863 bases; selected NP_001276674.1; 20 returned candidates.
- Both genomic snapshots are included in `tests/fixtures`; deterministic tests
  check each full target against genomic DNA and each cut against selected CDS.
- Minus-strand mapping is covered by ACTB and synthetic both-orientation tests.
- The Ensembl path is covered by mocked response and coordinate tests. Live
  Ensembl requests timed out during review, so live service behavior is not certified.
- GuideScan2 argument/CSV contracts are tested with simulated output. No matching
  genome index was run. RS2 provider was absent; wrapper/context behavior is
  checked but end-to-end model predictions were not validated in this environment.
- Source compilation and `git diff --check` completed without errors.

The scope-aware validation tests reproduce both follow-up findings. A verified
radius-zero screen with MIT 100 remains REVIEW because the search scope is
incomplete. A one-mismatch position-1 site produces MIT specificity 50 and CFD
specificity about 52.63, is counted as High, and also remains REVIEW. A clean
verified radius-three control still receives LOCAL CHECKS MET.

## Literature and source interpretation

Hsu et al. (2013), DOI [10.1038/nbt.2647](https://doi.org/10.1038/nbt.2647),
characterized mismatch-position-dependent SpCas9 specificity in human cells.
This supports reporting a mismatch-aware score while preserving its model scope.
It does not support assuming that the first exact reference match is intended.

Doench et al. (2016), DOI [10.1038/nbt.3437](https://doi.org/10.1038/nbt.3437),
described activity and off-target models. The implementation separates RS2 from
CFD and from the app's hand-coded sequence heuristic. The bundled weights and
independent scoring implementation are documented in `third_party/README.md`.

Gilbert et al. (2014), DOI [10.1016/j.cell.2014.09.029](https://doi.org/10.1016/j.cell.2014.09.029),
and Horlbeck et al. (2016), DOI [10.7554/eLife.19760](https://doi.org/10.7554/eLife.19760),
provide the mammalian CRISPRi/a context. This review treats TSS bands as a placement
heuristic, not a universal rule for plants, bacteria, every repressor or cell type.

Poudel et al. (2022), DOI [10.1093/gigascience/giac007](https://doi.org/10.1093/gigascience/giac007),
and the [GuideMaker implementation](https://github.com/USDA-ARS-GBRU/GuideMaker)
were used to check the optional RS2 interface and CFD data provenance.
[GuideScan2 source](https://github.com/pritykinlab/guidescan-cli) was checked for
its numeric input-coordinate contract and native per-hit output conventions.

These sources explain model origins and boundaries. Regression examples supply
the evidence that the four software defects have been fixed; papers do not
independently validate this application or its selected guides.

## Remaining limitations

- Local NGG-only substitution screening excludes alternate PAMs, DNA/RNA bulges,
  variants, rearrangements, chromatin and circular-origin targets. It covers only
  the supplied sequences and selected mismatch radius.
- A high local score is not a genome-wide specificity guarantee. Local checks can
  be met on an incomplete reference; users must supply a suitable reference.
- Pasted DNA cannot be authenticated as genomic from sequence text alone.
- Accession/paste modes do not infer coding annotation, frameshift probability,
  essential domains or loss of function. The gene workflow filters cuts to CDS
  but still does not predict these outcomes.
- Ranking remains heuristic. Experimental activity, specificity and phenotypes
  require independent biological validation.
- Ensembl gene-symbol knockout is deliberately unavailable pending an annotated
  genomic implementation. Ensembl genomic accession mode lacks CDS filtering.
- API availability and annotation releases can change. The included snapshots
  establish reproducible examples, not universal species coverage.
- External GuideScan2 and optional RS2 installations/indexes need environment-
  specific validation; offline adapter tests are not a live integration benchmark.
