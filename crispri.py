# ============================================================================
# CRISPRI
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Contains CRISPR interference (CRISPRi) helper logic used when designing repression-oriented guides.
#
# HOW TO READ THIS FILE:
# 1. Start with the imports and constants.
# 2. Read each function separately; every function performs one part of the workflow.
# 3. Follow the function calls from the application/workflow rather than trying to
#    understand the entire file at once.
# 4. Scientific equations, thresholds, validation rules and public function names
#    are intentionally kept unchanged while the code is being humanized.
#
# MAIN TOP-LEVEL PARTS IN THIS FILE:
# - function: species_slug
# - class: TSSContext
# - class: CRISPRiGuide
# - function: _canonical_transcript
# - function: fetch_ensembl_tss_context
# - function: tss_band
# - function: _spacer_interval
# - function: _genomic_interval
# - function: _genomic_guide_strand
# - function: design_crispri_guides
# - function: design_crispri_from_sequence
# ============================================================================

"""TSS-aware CRISPRi design utilities for CRISPR Studio.

Gene-based CRISPRi uses the canonical Ensembl transcript TSS and retrieves a
strand-aware genomic sequence window. Candidate guides are filtered to the
well-established dCas9-KRAB CRISPRi window from -50 to +300 bp relative to
that TSS, with +50 to +100 bp marked as the preferred placement band.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

import requests
from network import get

from grna_designer import GuideRNA, clean_dna_sequence, design_guides, screen_guides

ENSEMBL_REST = "https://rest.ensembl.org"
CRISPRI_MIN = -50
CRISPRI_MAX = 300
CRISPRI_OPTIMAL_MIN = 50
CRISPRI_OPTIMAL_MAX = 100

SPECIES_ALIASES = {
    "homo sapiens": "homo_sapiens",
    "human": "homo_sapiens",
    "mus musculus": "mus_musculus",
    "mouse": "mus_musculus",
    "rattus norvegicus": "rattus_norvegicus",
    "rat": "rattus_norvegicus",
    "danio rerio": "danio_rerio",
    "zebrafish": "danio_rerio",
    "drosophila melanogaster": "drosophila_melanogaster",
    "arabidopsis thaliana": "arabidopsis_thaliana",
    "saccharomyces cerevisiae": "saccharomyces_cerevisiae",
}



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: species_slug
# ----------------------------------------------------------------------------
def species_slug(organism: str) -> str:
    key = organism.strip().lower()
    return SPECIES_ALIASES.get(key, key.replace(" ", "_"))


@dataclass(frozen=True)
class TSSContext:
    gene: str
    organism: str
    species: str
    gene_id: str
    transcript_id: str
    assembly: str
    chromosome: str
    transcript_strand: int
    tss_coordinate: int
    region_start: int
    region_end: int
    tss_offset: int
    sequence: str
    description: str

    @property
    def strand_label(self) -> str:
        return "+" if self.transcript_strand == 1 else "-"


@dataclass
class CRISPRiGuide:
    guide: GuideRNA
    tss_distance: float
    tss_band: str
    tss_priority: int
    genomic_start: Optional[int]
    genomic_end: Optional[int]
    genomic_strand: Optional[str]
    chromosome: Optional[str]
    assembly: Optional[str]
    transcript_id: Optional[str]
    tss_coordinate: Optional[int]

    def to_dict(self) -> Dict[str, object]:
        row = self.guide.to_dict()
        row.update(
            {
                "TSS distance (bp)": self.tss_distance,
                "TSS band": self.tss_band,
                "Chromosome": self.chromosome,
                "Genomic start": self.genomic_start,
                "Genomic end": self.genomic_end,
                "Genomic strand": self.genomic_strand,
                "Assembly": self.assembly,
                "Transcript": self.transcript_id,
                "TSS coordinate": self.tss_coordinate,
            }
        )
        return row



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _canonical_transcript
# ----------------------------------------------------------------------------
def _canonical_transcript(gene_record: Dict[str, object]) -> Dict[str, object]:
    transcripts = gene_record.get("Transcript", []) or []
    if not transcripts:
        raise ValueError("Ensembl returned no transcripts for this gene.")
    canonical = str(gene_record.get("canonical_transcript", "")).split(".")[0]
    for transcript in transcripts:
        tid = str(transcript.get("id", "")).split(".")[0]
        if transcript.get("is_canonical") or (canonical and tid == canonical):
            return transcript
    raise ValueError("No canonical transcript is annotated; provide an explicit transcript ID.")



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: fetch_ensembl_tss_context
# ----------------------------------------------------------------------------
def fetch_ensembl_tss_context(
    gene_name: str,
    organism: str,
    fetch_upstream: int = 100,
    fetch_downstream: int = 350,
    transcript_id: str = "",
) -> TSSContext:
    """Fetch the canonical-transcript genomic TSS window in transcriptional orientation."""
    species = species_slug(organism)
    if species not in {'homo_sapiens', 'mus_musculus'}:
        raise ValueError('The built-in CRISPRi placement heuristic supports human/mouse dCas9-KRAB, not other organisms or repressors.')
    if not 0 <= fetch_upstream <= 100_000 or not 0 <= fetch_downstream <= 100_000:
        raise ValueError('TSS flanks must be between 0 and 100,000 bp.')
    headers = {"Accept": "application/json", "User-Agent": "CRISPR_gRNA_Designer"}
    lookup = get(
        f"{ENSEMBL_REST}/lookup/symbol/{species}/{gene_name}?expand=1",
        headers=headers,
        timeout=30,
    )
    if not lookup.ok:
        detail = lookup.text[:180].replace("\n", " ")
        raise ValueError(f"Ensembl TSS lookup failed ({lookup.status_code}): {detail}")

    gene = lookup.json()
    transcript = next((t for t in gene.get('Transcript', []) if t.get('id') == transcript_id), None) if transcript_id else _canonical_transcript(gene)
    if transcript is None:
        raise ValueError('Requested transcript is not annotated for this gene.')
    strand = int(transcript.get("strand") or gene.get("strand") or 0)
    if strand not in (-1, 1):
        raise ValueError("Ensembl did not return a valid transcript strand.")

    chromosome = str(transcript.get("seq_region_name") or gene.get("seq_region_name") or "")
    if not chromosome:
        raise ValueError("Ensembl did not return a chromosome/contig for the transcript.")

    tx_start = int(transcript["start"])
    tx_end = int(transcript["end"])
    tss = tx_start if strand == 1 else tx_end

    if strand == 1:
        region_start = max(1, tss - fetch_upstream)
        region_end = tss + fetch_downstream
        tss_offset = tss - region_start
    else:
        region_start = max(1, tss - fetch_downstream)
        region_end = tss + fetch_upstream
        tss_offset = region_end - tss

    region = f"{chromosome}:{region_start}..{region_end}:{strand}"
    sequence_response = get(
        f"{ENSEMBL_REST}/sequence/region/{species}/{region}",
        headers={"Accept": "text/plain", "User-Agent": "CRISPR_gRNA_Designer"},
        timeout=30,
    )
    if not sequence_response.ok:
        detail = sequence_response.text[:180].replace("\n", " ")
        raise ValueError(
            f"Could not fetch genomic TSS window ({sequence_response.status_code}): {detail}"
        )

    sequence = clean_dna_sequence(sequence_response.text)
    expected = region_end - region_start + 1
    if len(sequence) != expected:
        raise ValueError(
            f"Unexpected Ensembl sequence length: got {len(sequence)}, expected {expected}."
        )

    transcript_id = str(transcript.get("id", ""))
    assembly = str(gene.get("assembly_name") or transcript.get("assembly_name") or "Unknown")
    description = (
        f"Ensembl annotated transcript TSS | {transcript_id} | {assembly} {chromosome}:{tss} "
        f"({ '+' if strand == 1 else '-' } strand)"
    )
    return TSSContext(
        gene=gene_name,
        organism=organism,
        species=species,
        gene_id=str(gene.get("id", "")),
        transcript_id=transcript_id,
        assembly=assembly,
        chromosome=chromosome,
        transcript_strand=strand,
        tss_coordinate=tss,
        region_start=region_start,
        region_end=region_end,
        tss_offset=tss_offset,
        sequence=sequence,
        description=description,
    )



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: tss_band
# ----------------------------------------------------------------------------
def tss_band(distance: int) -> tuple[str, int]:
    if CRISPRI_OPTIMAL_MIN <= distance <= CRISPRI_OPTIMAL_MAX:
        return "Preferred (+50 to +100)", 3
    if 0 <= distance < CRISPRI_OPTIMAL_MIN or CRISPRI_OPTIMAL_MAX < distance <= 200:
        return "High-priority", 2
    return "CRISPRi window", 1



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _spacer_interval
# ----------------------------------------------------------------------------
def _spacer_interval(guide: GuideRNA) -> tuple[int, int]:
    if guide.strand == "+":
        return guide.start, guide.end - 3
    return guide.start + 3, guide.end



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _genomic_interval
# ----------------------------------------------------------------------------
def _genomic_interval(context: TSSContext, start: int, end: int) -> tuple[int, int]:
    """Map a 0-based half-open interval in transcription-oriented sequence to 1-based genome coords."""
    if context.transcript_strand == 1:
        g_start = context.region_start + start
        g_end = context.region_start + end - 1
    else:
        g_start = context.region_end - end + 1
        g_end = context.region_end - start
    return min(g_start, g_end), max(g_start, g_end)



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _genomic_guide_strand
# ----------------------------------------------------------------------------
def _genomic_guide_strand(context: TSSContext, guide_strand: str) -> str:
    if context.transcript_strand == 1:
        return guide_strand
    return "+" if guide_strand == "-" else "-"



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: design_crispri_guides
# ----------------------------------------------------------------------------
def design_crispri_guides(
    context: TSSContext,
    max_guides: int = 20,
    min_score: float = 30.0,
    genome_context=None,
    max_mismatches: int = 3,
    intended_region=None,
) -> List[CRISPRiGuide]:
    """Design and rank SpCas9 CRISPRi guides by true TSS-relative placement."""
    if type(max_guides) is not int or max_guides < 1:
        raise ValueError('Maximum guides must be a positive integer.')
    if context.species not in {'homo_sapiens', 'mus_musculus', 'custom'}:
        raise ValueError('Built-in placement rules are restricted to human/mouse dCas9-KRAB.')
    candidates = design_guides(
        context.sequence,
        application="crispri",
        min_score=min_score,
        max_guides=50_000,
        prefer_5prime=False,
        genome_context=None,
        max_mismatches=max_mismatches,
    )

    annotated: List[CRISPRiGuide] = []
    for guide in candidates:
        spacer_start, spacer_end = _spacer_interval(guide)
        midpoint = (spacer_start + spacer_end - 1) / 2
        distance = midpoint - context.tss_offset
        if not (CRISPRI_MIN <= distance <= CRISPRI_MAX):
            continue

        band, priority = tss_band(distance)
        genomic_start, genomic_end = _genomic_interval(context, spacer_start, spacer_end)
        guide.application = "crispri"
        guide.notes.append(f"TSS-aware: {distance:+g} bp")
        if priority == 3:
            guide.notes.append("Preferred CRISPRi placement band")

        annotated.append(
            CRISPRiGuide(
                guide=guide,
                tss_distance=distance,
                tss_band=band,
                tss_priority=priority,
                genomic_start=genomic_start,
                genomic_end=genomic_end,
                genomic_strand=_genomic_guide_strand(context, guide.strand),
                chromosome=context.chromosome,
                assembly=context.assembly,
                transcript_id=context.transcript_id,
                tss_coordinate=context.tss_coordinate,
            )
        )

    if genome_context is not None:
        screen_guides([item.guide for item in annotated], genome_context, max_mismatches, intended_region=intended_region, target_sequence=context.sequence)
    annotated.sort(
        key=lambda item: (
            item.tss_priority,
            item.guide.specificity_score if item.guide.specificity_score is not None else -1,
            item.guide.score,
            -abs(item.tss_distance - 75),
        ),
        reverse=True,
    )
    return annotated[:max_guides]



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: design_crispri_from_sequence
# ----------------------------------------------------------------------------
def design_crispri_from_sequence(
    sequence: str,
    tss_position_1based: int,
    max_guides: int = 20,
    min_score: float = 30.0,
    genome_context=None,
    max_mismatches: int = 3,
    intended_region=None,
) -> List[CRISPRiGuide]:
    """TSS-aware CRISPRi for a user-provided transcription-oriented genomic sequence."""
    sequence = clean_dna_sequence(sequence)
    if type(tss_position_1based) is not int or not 1 <= tss_position_1based <= len(sequence):
        raise ValueError("TSS position must fall inside the pasted sequence.")
    context = TSSContext(
        gene="Custom target",
        organism="User supplied",
        species="custom",
        gene_id="",
        transcript_id="CUSTOM_TSS",
        assembly="Custom",
        chromosome="custom",
        transcript_strand=1,
        tss_coordinate=tss_position_1based,
        region_start=1,
        region_end=len(sequence),
        tss_offset=tss_position_1based - 1,
        sequence=sequence,
        description="User-provided genomic sequence with declared TSS; sequence assumed 5'→3' in transcriptional orientation.",
    )
    guides = design_crispri_guides(
        context,
        max_guides=max_guides,
        min_score=min_score,
        genome_context=genome_context,
        max_mismatches=max_mismatches,
        intended_region=intended_region,
    )
    for item in guides:
        item.chromosome = None
        item.assembly = None
        item.genomic_start = None
        item.genomic_end = None
        item.genomic_strand = None
        item.tss_coordinate = tss_position_1based
    return guides
