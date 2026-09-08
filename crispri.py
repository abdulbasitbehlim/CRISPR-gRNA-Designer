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

from grna_designer import GuideRNA, clean_dna_sequence, design_guides

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
    tss_distance: int
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


def _canonical_transcript(gene_record: Dict[str, object]) -> Dict[str, object]:
    transcripts = gene_record.get("Transcript", []) or []
    if not transcripts:
        raise ValueError("Ensembl returned no transcripts for this gene.")
    canonical = str(gene_record.get("canonical_transcript", "")).split(".")[0]
    for transcript in transcripts:
        tid = str(transcript.get("id", "")).split(".")[0]
        if transcript.get("is_canonical") or (canonical and tid == canonical):
            return transcript
    protein_coding = [t for t in transcripts if t.get("biotype") == "protein_coding"]
    return protein_coding[0] if protein_coding else transcripts[0]


def fetch_ensembl_tss_context(
    gene_name: str,
    organism: str,
    fetch_upstream: int = 100,
    fetch_downstream: int = 350,
) -> TSSContext:
    """Fetch the canonical-transcript genomic TSS window in transcriptional orientation."""
    species = species_slug(organism)
    headers = {"Accept": "application/json", "User-Agent": "CRISPR_gRNA_Designer"}
    lookup = requests.get(
        f"{ENSEMBL_REST}/lookup/symbol/{species}/{gene_name}?expand=1",
        headers=headers,
        timeout=30,
    )
    if not lookup.ok:
        detail = lookup.text[:180].replace("\n", " ")
        raise ValueError(f"Ensembl TSS lookup failed ({lookup.status_code}): {detail}")

    gene = lookup.json()
    transcript = _canonical_transcript(gene)
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
    sequence_response = requests.get(
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
        f"Ensembl canonical TSS | {transcript_id} | {assembly} {chromosome}:{tss} "
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


def tss_band(distance: int) -> tuple[str, int]:
    if CRISPRI_OPTIMAL_MIN <= distance <= CRISPRI_OPTIMAL_MAX:
        return "Preferred (+50 to +100)", 3
    if 0 <= distance < CRISPRI_OPTIMAL_MIN or CRISPRI_OPTIMAL_MAX < distance <= 200:
        return "High-priority", 2
    return "CRISPRi window", 1


def _spacer_interval(guide: GuideRNA) -> tuple[int, int]:
    if guide.strand == "+":
        return guide.start, guide.end - 3
    return guide.start + 3, guide.end


def _genomic_interval(context: TSSContext, start: int, end: int) -> tuple[int, int]:
    """Map a 0-based half-open interval in transcription-oriented sequence to 1-based genome coords."""
    if context.transcript_strand == 1:
        g_start = context.region_start + start
        g_end = context.region_start + end - 1
    else:
        g_start = context.region_end - end + 1
        g_end = context.region_end - start
    return min(g_start, g_end), max(g_start, g_end)


def _genomic_guide_strand(context: TSSContext, guide_strand: str) -> str:
    if context.transcript_strand == 1:
        return guide_strand
    return "+" if guide_strand == "-" else "-"


def design_crispri_guides(
    context: TSSContext,
    max_guides: int = 20,
    min_score: float = 30.0,
    genome_context=None,
    max_mismatches: int = 3,
) -> List[CRISPRiGuide]:
    """Design and rank SpCas9 CRISPRi guides by true TSS-relative placement."""
    candidates = design_guides(
        context.sequence,
        application="knockout",
        min_score=min_score,
        max_guides=200,
        prefer_5prime=False,
        genome_context=genome_context,
        max_mismatches=max_mismatches,
    )

    annotated: List[CRISPRiGuide] = []
    for guide in candidates:
        spacer_start, spacer_end = _spacer_interval(guide)
        midpoint = int(round((spacer_start + spacer_end - 1) / 2))
        distance = midpoint - context.tss_offset
        if not (CRISPRI_MIN <= distance <= CRISPRI_MAX):
            continue

        band, priority = tss_band(distance)
        genomic_start, genomic_end = _genomic_interval(context, spacer_start, spacer_end)
        guide.application = "crispri"
        guide.notes.append(f"TSS-aware: {distance:+d} bp")
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

    annotated.sort(
        key=lambda item: (
            item.tss_priority,
            item.guide.specificity_score if item.guide.specificity_score is not None else -1,
            item.guide.doench_score if item.guide.doench_score is not None else item.guide.score,
            -abs(item.tss_distance - 75),
        ),
        reverse=True,
    )
    return annotated[:max_guides]


def design_crispri_from_sequence(
    sequence: str,
    tss_position_1based: int,
    max_guides: int = 20,
    min_score: float = 30.0,
    genome_context=None,
    max_mismatches: int = 3,
) -> List[CRISPRiGuide]:
    """TSS-aware CRISPRi for a user-provided transcription-oriented genomic sequence."""
    sequence = clean_dna_sequence(sequence)
    if not 1 <= tss_position_1based <= len(sequence):
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
    )
    for item in guides:
        item.chromosome = None
        item.assembly = None
        item.genomic_start = None
        item.genomic_end = None
        item.genomic_strand = None
        item.tss_coordinate = tss_position_1based
    return guides
