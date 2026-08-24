#!/usr/bin/env python3
"""Core sequence retrieval, SpCas9 guide design, and local off-target screening.

The scoring model in this project is intentionally transparent. It combines a
small set of published sequence preferences into a 0-100 ranking heuristic; it
is not the trained Doench Rule Set 2 model and it is not a replacement for a
genome-indexed off-target tool.
"""

from __future__ import annotations

import os
import re
import time
from dataclasses import dataclass, field
from io import StringIO
from typing import Any, Dict, List, Optional, Tuple

import requests
from Bio import SeqIO
from Bio.Seq import Seq


# NCBI asks callers to identify themselves. Deployers can set NCBI_EMAIL in
# Streamlit Cloud secrets/environment variables without changing this file.
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "crispr.grna.designer@example.com")
NCBI_TOOL = "CRISPR_gRNA_Designer"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

PAM_PATTERN = re.compile(r"(?=([ATCG]GG))", re.IGNORECASE)
PAM_PATTERN_REV = re.compile(r"(?=(CC[ATCG]))", re.IGNORECASE)
DNA_ALPHABET = frozenset("ATGCN")
IUPAC_AMBIGUOUS = re.compile(r"[RYSWKMBDHVX]")


@dataclass(frozen=True)
class OffTargetHit:
    """One PAM-compatible near-match in a user-supplied reference sequence."""

    start: int
    strand: str
    sequence: str
    pam: str
    mismatches: int
    mismatch_positions: Tuple[int, ...]
    seed_mismatches: int
    risk: str

    def to_dict(self) -> Dict[str, Any]:
        return {
            "Position": self.start + 1,
            "Strand": self.strand,
            "Candidate spacer": self.sequence,
            "PAM": self.pam,
            "Mismatches": self.mismatches,
            "Mismatch positions": ", ".join(map(str, self.mismatch_positions)) or "Exact",
            "Seed mismatches": self.seed_mismatches,
            "Risk": self.risk,
        }


@dataclass(frozen=True)
class OffTargetReport:
    """Summary of a local-reference similarity screen for one spacer."""

    hits: Tuple[OffTargetHit, ...]
    specificity_score: float
    pam_sites_scanned: int
    on_target_excluded: bool


@dataclass
class GuideRNA:
    """A ranked 20 nt SpCas9 spacer and its PAM."""

    sequence: str
    pam: str
    strand: str
    start: int
    end: int
    gc_content: float
    score: float
    notes: List[str] = field(default_factory=list)
    application: str = "knockout"
    specificity_score: Optional[float] = None
    off_target_count: Optional[int] = None
    off_target_details: List[Dict[str, Any]] = field(default_factory=list)

    @property
    def full_target(self) -> str:
        return self.sequence + self.pam

    def to_dict(self) -> Dict[str, Any]:
        return {
            "Spacer (20 nt)": self.sequence,
            "PAM": self.pam,
            "Strand": self.strand,
            "Start": self.start + 1,
            "End": self.end,
            "GC%": round(self.gc_content, 1),
            "Score": round(self.score, 1),
            "Specificity": (
                round(self.specificity_score, 1)
                if self.specificity_score is not None
                else None
            ),
            "Off-target hits": self.off_target_count,
            "Application": self.application,
            "Notes": "; ".join(self.notes) if self.notes else "No flags",
        }


# ---------------------------------------------------------------------------
# Sequence handling and retrieval
# ---------------------------------------------------------------------------
def clean_dna_sequence(raw_sequence: str) -> str:
    """Normalize plain DNA/RNA or FASTA text while preserving coordinates.

    Whitespace and line numbers are removed, U is converted to T, and standard
    ambiguous IUPAC bases are converted to N. Unexpected symbols raise a clear
    error instead of being silently deleted.
    """

    if not raw_sequence:
        return ""

    sequence_lines = [
        line.strip()
        for line in raw_sequence.splitlines()
        if line.strip() and not line.lstrip().startswith(">")
    ]
    sequence = "".join(sequence_lines).upper().replace("U", "T")
    sequence = re.sub(r"[\s\d]", "", sequence)
    sequence = IUPAC_AMBIGUOUS.sub("N", sequence)
    invalid = sorted(set(sequence) - DNA_ALPHABET)
    if invalid:
        raise ValueError(
            "Sequence contains unsupported character(s): " + ", ".join(invalid)
        )
    return sequence


def fetch_gene_ncbi(
    gene_name: str,
    organism: str,
    max_retries: int = 3,
) -> Tuple[str, str, str]:
    """Fetch a representative RefSeq transcript from NCBI Gene."""

    query = f"{gene_name}[Gene Name] AND {organism}[Organism] AND alive[prop]"
    for attempt in range(max_retries):
        try:
            common_params = {"tool": NCBI_TOOL, "email": NCBI_EMAIL}
            search_response = requests.get(
                f"{NCBI_EUTILS}/esearch.fcgi",
                params={
                    **common_params,
                    "db": "gene",
                    "term": query,
                    "retmax": 5,
                    "retmode": "json",
                },
                timeout=30,
            )
            search_response.raise_for_status()
            gene_ids = search_response.json().get("esearchresult", {}).get("idlist", [])
            if not gene_ids:
                broader_query = f"{gene_name} AND {organism}[Organism]"
                search_response = requests.get(
                    f"{NCBI_EUTILS}/esearch.fcgi",
                    params={
                        **common_params,
                        "db": "gene",
                        "term": broader_query,
                        "retmax": 5,
                        "retmode": "json",
                    },
                    timeout=30,
                )
                search_response.raise_for_status()
                gene_ids = search_response.json().get("esearchresult", {}).get(
                    "idlist", []
                )
            if not gene_ids:
                raise ValueError(
                    f"No gene found for '{gene_name}' in '{organism}' on NCBI."
                )

            gene_id = gene_ids[0]

            def linked_nucleotide_ids(link_name: Optional[str]) -> List[str]:
                params: Dict[str, Any] = {
                    **common_params,
                    "dbfrom": "gene",
                    "db": "nuccore",
                    "id": gene_id,
                    "retmode": "json",
                }
                if link_name:
                    params["linkname"] = link_name
                link_response = requests.get(
                    f"{NCBI_EUTILS}/elink.fcgi",
                    params=params,
                    timeout=30,
                )
                link_response.raise_for_status()
                linksets = link_response.json().get("linksets", [])
                if not linksets:
                    return []
                link_databases = linksets[0].get("linksetdbs", [])
                if not link_databases:
                    return []
                return list(link_databases[0].get("links", []))

            nucleotide_ids = linked_nucleotide_ids("gene_nuccore_refseqrna")

            if not nucleotide_ids:
                nucleotide_ids = linked_nucleotide_ids(None)

            if not nucleotide_ids:
                raise ValueError("Could not resolve a nucleotide record for this gene.")

            fetch_response = requests.get(
                f"{NCBI_EUTILS}/efetch.fcgi",
                params={
                    **common_params,
                    "db": "nuccore",
                    "id": ",".join(nucleotide_ids[:30]),
                    "rettype": "gb",
                    "retmode": "text",
                },
                timeout=30,
            )
            fetch_response.raise_for_status()
            records = list(SeqIO.parse(StringIO(fetch_response.text), "genbank"))
            if not records:
                raise ValueError("NCBI returned no readable nucleotide records.")

            def transcript_priority(candidate: Any) -> Tuple[int, int, str]:
                accession = candidate.id.upper()
                description = candidate.description.upper()
                record_class = (
                    0
                    if accession.startswith("NM_")
                    else 1
                    if accession.startswith("XM_")
                    else 2
                )
                variant_priority = (
                    0
                    if "TRANSCRIPT VARIANT 1" in description
                    else 1
                    if "TRANSCRIPT VARIANT" not in description
                    else 2
                )
                return record_class, variant_priority, accession

            coding_records = [
                candidate
                for candidate in records
                if "MRNA" in candidate.description.upper()
                or any(feature.type == "CDS" for feature in candidate.features)
            ]
            selected = min(coding_records or records, key=transcript_priority)
            sequence = clean_dna_sequence(str(selected.seq))
            return selected.id, selected.description, sequence
        except Exception:
            if attempt == max_retries - 1:
                raise
            time.sleep(1 + attempt)

    raise RuntimeError("NCBI lookup failed after retries.")


def fetch_gene_ensembl(
    gene_name: str,
    species: str = "homo_sapiens",
) -> Tuple[str, str, str]:
    """Fetch the canonical Ensembl transcript cDNA for a gene symbol."""

    server = "https://rest.ensembl.org"
    species = species.lower().strip().replace(" ", "_")
    headers = {"Accept": "application/json", "User-Agent": NCBI_TOOL}
    lookup_url = f"{server}/lookup/symbol/{species}/{gene_name}?expand=1"
    response = requests.get(lookup_url, headers=headers, timeout=30)
    if not response.ok:
        detail = response.text[:180].replace("\n", " ")
        raise ValueError(f"Ensembl lookup failed ({response.status_code}): {detail}")

    data = response.json()
    transcripts = data.get("Transcript", [])
    if not transcripts:
        raise ValueError("No transcripts were returned for this gene by Ensembl.")

    canonical_id = str(data.get("canonical_transcript", "")).split(".")[0]
    transcript = next(
        (
            item
            for item in transcripts
            if item.get("is_canonical") or item.get("id") == canonical_id
        ),
        transcripts[0],
    )
    transcript_id = transcript["id"]

    sequence_response = requests.get(
        f"{server}/sequence/id/{transcript_id}?type=cdna",
        headers={"Accept": "text/plain", "User-Agent": NCBI_TOOL},
        timeout=30,
    )
    if not sequence_response.ok:
        raise ValueError(f"Could not fetch cDNA sequence for {transcript_id}.")

    sequence = clean_dna_sequence(sequence_response.text)
    description = (
        f"Ensembl {transcript_id} | {data.get('display_name', gene_name)} | "
        f"{data.get('description', '')}"
    ).strip(" |")
    return transcript_id, description, sequence


def fetch_sequence(
    gene_name: str,
    organism: str,
    source: str = "ncbi",
) -> Tuple[str, str, str]:
    """Fetch a representative transcript from NCBI or Ensembl."""

    source = source.lower().strip()
    if source not in {"ncbi", "ensembl"}:
        raise ValueError("Source must be either 'ncbi' or 'ensembl'.")
    if source == "ncbi":
        return fetch_gene_ncbi(gene_name, organism)

    species_map = {
        "homo sapiens": "homo_sapiens",
        "human": "homo_sapiens",
        "mus musculus": "mus_musculus",
        "mouse": "mus_musculus",
        "rattus norvegicus": "rattus_norvegicus",
        "rat": "rattus_norvegicus",
        "danio rerio": "danio_rerio",
        "zebrafish": "danio_rerio",
        "saccharomyces cerevisiae": "saccharomyces_cerevisiae",
        "arabidopsis thaliana": "arabidopsis_thaliana",
        "drosophila melanogaster": "drosophila_melanogaster",
    }
    species = species_map.get(organism.lower(), organism.lower().replace(" ", "_"))
    return fetch_gene_ensembl(gene_name, species)


# ---------------------------------------------------------------------------
# Transparent guide scoring
# ---------------------------------------------------------------------------
def _gc_content(sequence: str) -> float:
    sequence = sequence.upper()
    return (
        100.0 * (sequence.count("G") + sequence.count("C")) / len(sequence)
        if sequence
        else 0.0
    )


def _has_homopolymer(sequence: str, max_run: int = 4) -> bool:
    pattern = r"(A{%d,}|T{%d,}|G{%d,}|C{%d,})" % (
        max_run,
        max_run,
        max_run,
        max_run,
    )
    return bool(re.search(pattern, sequence.upper()))


def score_breakdown(
    spacer: str,
    pam: str,
    application: Optional[str] = None,
    start: Optional[int] = None,
    sequence_length: Optional[int] = None,
) -> Dict[str, float]:
    """Return every contribution to the project's ranking heuristic."""

    spacer = spacer.upper()
    pam = pam.upper()
    if len(spacer) != 20 or len(pam) != 3:
        return {"Final score": 0.0}

    components: Dict[str, float] = {"Baseline": 50.0}
    gc = _gc_content(spacer)
    if 40 <= gc <= 70:
        components["GC content"] = 15.0
    elif 30 <= gc < 40 or 70 < gc <= 80:
        components["GC content"] = 5.0
    else:
        components["GC content"] = -15.0

    components["Position 20"] = 10.0 if spacer[19] == "G" else 0.0
    if spacer[19] == "C":
        components["Position 20"] = -8.0
    elif spacer[19] == "T":
        components["Position 20"] = -12.0
    components["Position 19"] = 5.0 if spacer[18] in "AG" else 0.0
    components["Position 16"] = 3.0 if spacer[15] == "C" else 0.0
    components["Position 18"] = 3.0 if spacer[17] == "C" else 0.0
    components["PAM context"] = 5.0 if pam == "CGG" else 0.0
    if pam.startswith("T"):
        components["PAM context"] = -5.0

    homopolymer_penalty = 0.0
    if "GGGG" in spacer or "TTTT" in spacer:
        homopolymer_penalty -= 20.0
    if "GGG" in spacer or "TTT" in spacer:
        homopolymer_penalty -= 8.0
    components["Homopolymers"] = homopolymer_penalty

    reverse_complement = str(Seq(spacer).reverse_complement())
    complement_matches = sum(
        left == right for left, right in zip(spacer, reverse_complement)
    )
    components["Self-complementarity"] = -10.0 if complement_matches > 12 else 0.0

    position_adjustment = 0.0
    if application and start is not None and sequence_length:
        if application.lower() == "knockdown":
            position_adjustment = 15.0 if start < 400 else -5.0
        elif application.lower() == "knockout" and start < sequence_length * 0.3:
            position_adjustment = 8.0
    components["Application position"] = position_adjustment

    raw_total = sum(components.values())
    components["Final score"] = round(max(0.0, min(100.0, raw_total)), 1)
    return components


def _doench_like_score(spacer: str, pam: str) -> float:
    """Backward-compatible entry point for the sequence-only score."""

    return score_breakdown(spacer, pam)["Final score"]


def _guide_notes(
    spacer: str,
    gc_content: float,
    application: str,
    start: int,
    sequence_length: int,
) -> List[str]:
    notes: List[str] = []
    if not 40 <= gc_content <= 70:
        notes.append("GC outside preferred 40-70%")
    if gc_content < 20 or gc_content > 80:
        notes.append("Extreme GC")
    if _has_homopolymer(spacer, 4):
        notes.append("Homopolymer >=4")
    if spacer.endswith("T"):
        notes.append("Ends with T (U6 risk)")
    if application == "knockdown" and start < 400:
        notes.append("5-prime-proximal preference")
    elif application == "knockout" and start < sequence_length * 0.3:
        notes.append("Early-region preference")
    return notes


# ---------------------------------------------------------------------------
# Local-reference off-target screen
# ---------------------------------------------------------------------------
def _risk_label(mismatches: int, seed_mismatches: int) -> str:
    if mismatches == 0:
        return "Critical"
    if mismatches == 1 or (mismatches == 2 and seed_mismatches <= 1):
        return "High"
    if mismatches <= 3 and seed_mismatches <= 1:
        return "Moderate"
    return "Low"


def _risk_weight(mismatches: int, seed_mismatches: int) -> float:
    base_weight = {0: 100.0, 1: 25.0, 2: 8.0, 3: 2.0}.get(mismatches, 0.5)
    seed_multiplier = 1.5 if seed_mismatches == 0 else (1.0 if seed_mismatches == 1 else 0.5)
    return base_weight * seed_multiplier


def analyze_offtargets(
    spacer: str,
    whole_genome: Optional[str],
    max_mismatches: int = 3,
    max_hits: int = 250,
) -> OffTargetReport:
    """Screen a spacer against PAM-compatible sites in a supplied reference.

    Both DNA strands are scanned. One exact match is treated as the intended
    target and excluded; additional exact matches remain critical hits. The
    final score is a transparent ranking aid, not a validated CFD score.
    """

    spacer = clean_dna_sequence(spacer)
    if len(spacer) != 20:
        raise ValueError("Off-target screening requires a 20 nt spacer.")
    if not 0 <= max_mismatches <= 6:
        raise ValueError("max_mismatches must be between 0 and 6.")
    if not whole_genome:
        return OffTargetReport((), 0.0, 0, False)

    reference = clean_dna_sequence(whole_genome)
    candidates: List[Tuple[int, str, str, str]] = []

    for match in PAM_PATTERN.finditer(reference):
        pam_start = match.start()
        spacer_start = pam_start - 20
        if spacer_start < 0:
            continue
        candidate = reference[spacer_start:pam_start]
        pam = reference[pam_start : pam_start + 3]
        if "N" not in candidate and len(pam) == 3:
            candidates.append((spacer_start, "+", candidate, pam))

    for match in PAM_PATTERN_REV.finditer(reference):
        pam_start = match.start()
        spacer_end = pam_start + 23
        if spacer_end > len(reference):
            continue
        candidate = str(Seq(reference[pam_start + 3 : spacer_end]).reverse_complement())
        pam = str(Seq(reference[pam_start : pam_start + 3]).reverse_complement())
        if "N" not in candidate:
            candidates.append((pam_start, "-", candidate, pam))

    exact_match_excluded = False
    hits: List[OffTargetHit] = []
    for start, strand, candidate, pam in candidates:
        mismatch_positions = tuple(
            index + 1
            for index, (query_base, candidate_base) in enumerate(zip(spacer, candidate))
            if query_base != candidate_base
        )
        mismatch_count = len(mismatch_positions)
        if mismatch_count > max_mismatches:
            continue
        if mismatch_count == 0 and not exact_match_excluded:
            exact_match_excluded = True
            continue

        seed_mismatches = sum(position >= 13 for position in mismatch_positions)
        hits.append(
            OffTargetHit(
                start=start,
                strand=strand,
                sequence=candidate,
                pam=pam,
                mismatches=mismatch_count,
                mismatch_positions=mismatch_positions,
                seed_mismatches=seed_mismatches,
                risk=_risk_label(mismatch_count, seed_mismatches),
            )
        )

    hits.sort(key=lambda hit: (hit.mismatches, hit.seed_mismatches, hit.start))
    retained_hits = hits[:max_hits]
    risk_sum = sum(
        _risk_weight(hit.mismatches, hit.seed_mismatches) for hit in retained_hits
    )
    specificity = 100.0 if not retained_hits else 100.0 / (1.0 + risk_sum / 20.0)
    return OffTargetReport(
        hits=tuple(retained_hits),
        specificity_score=round(specificity, 1),
        pam_sites_scanned=len(candidates),
        on_target_excluded=exact_match_excluded,
    )


def calculate_offtarget_score(
    spacer: str,
    pam: str,
    whole_genome: Optional[str] = None,
    max_mismatches: int = 3,
) -> Tuple[int, float, List[str]]:
    """Compatibility wrapper returning count, score, and readable hit lines."""

    del pam  # PAM compatibility is evaluated at every reference candidate site.
    if not whole_genome:
        return 0, 0.0, ["Genome not provided"]
    report = analyze_offtargets(spacer, whole_genome, max_mismatches=max_mismatches)
    descriptions = [
        (
            f"Position {hit.start + 1} | {hit.strand} strand | "
            f"{hit.sequence}{hit.pam} | {hit.mismatches} mismatch(es) | {hit.risk} risk"
        )
        for hit in report.hits
    ]
    return len(report.hits), report.specificity_score, descriptions


# ---------------------------------------------------------------------------
# Guide discovery and validation
# ---------------------------------------------------------------------------
def design_guides(
    sequence: str,
    application: str = "knockout",
    pam: str = "NGG",
    spacer_len: int = 20,
    min_score: float = 30.0,
    max_guides: int = 50,
    prefer_5prime: bool = True,
    genome_context: Optional[str] = None,
    max_mismatches: int = 3,
) -> List[GuideRNA]:
    """Scan both strands for SpCas9 NGG sites and return ranked guides."""

    if pam.upper() != "NGG" or spacer_len != 20:
        raise ValueError("This release supports SpCas9 with a 20 nt spacer and NGG PAM.")
    application = application.lower().strip()
    if application not in {"knockout", "knockdown"}:
        raise ValueError("Application must be 'knockout' or 'knockdown'.")

    sequence = clean_dna_sequence(sequence)
    guides: List[GuideRNA] = []

    def add_guide(
        spacer: str,
        pam_sequence: str,
        strand: str,
        start: int,
        end: int,
    ) -> None:
        gc_content = _gc_content(spacer)
        scoring_application = application if prefer_5prime else None
        score = score_breakdown(
            spacer,
            pam_sequence,
            application=scoring_application,
            start=start,
            sequence_length=len(sequence),
        )["Final score"]
        if score < min_score:
            return

        notes = _guide_notes(spacer, gc_content, application, start, len(sequence))
        guide = GuideRNA(
            sequence=spacer,
            pam=pam_sequence,
            strand=strand,
            start=start,
            end=end,
            gc_content=gc_content,
            score=score,
            notes=notes,
            application=application,
        )
        if genome_context:
            report = analyze_offtargets(
                spacer,
                genome_context,
                max_mismatches=max_mismatches,
            )
            guide.specificity_score = report.specificity_score
            guide.off_target_count = len(report.hits)
            guide.off_target_details = [hit.to_dict() for hit in report.hits]
            if report.hits:
                notes.append(f"Reference screen: {len(report.hits)} near-match(es)")
        guides.append(guide)

    for match in PAM_PATTERN.finditer(sequence):
        pam_start = match.start()
        spacer_start = pam_start - spacer_len
        if spacer_start < 0:
            continue
        spacer = sequence[spacer_start:pam_start]
        pam_sequence = sequence[pam_start : pam_start + 3]
        if "N" not in spacer and len(pam_sequence) == 3:
            add_guide(
                spacer,
                pam_sequence,
                "+",
                spacer_start,
                pam_start + 3,
            )

    for match in PAM_PATTERN_REV.finditer(sequence):
        pam_start = match.start()
        spacer_end = pam_start + 3 + spacer_len
        if spacer_end > len(sequence):
            continue
        genomic_spacer = sequence[pam_start + 3 : spacer_end]
        spacer = str(Seq(genomic_spacer).reverse_complement())
        pam_sequence = str(
            Seq(sequence[pam_start : pam_start + 3]).reverse_complement()
        )
        if "N" not in spacer:
            add_guide(spacer, pam_sequence, "-", pam_start, spacer_end)

    guides.sort(
        key=lambda guide: (
            guide.score,
            guide.specificity_score if guide.specificity_score is not None else -1,
            -guide.start,
        ),
        reverse=True,
    )
    return guides[:max_guides]


def validate_guide(
    guide: GuideRNA,
    genome_context: Optional[str] = None,
) -> Dict[str, Any]:
    """Return transparent quality checks for a guide candidate."""

    specificity = guide.specificity_score
    if genome_context and specificity is None:
        specificity = analyze_offtargets(guide.sequence, genome_context).specificity_score

    checks: Dict[str, Any] = {
        "length_ok": len(guide.sequence) == 20,
        "pam_ok": len(guide.pam) == 3 and guide.pam.upper().endswith("GG"),
        "gc_in_preferred_range": 40 <= guide.gc_content <= 70,
        "gc_in_acceptable_range": 30 <= guide.gc_content <= 80,
        "no_extreme_homopolymer": not _has_homopolymer(guide.sequence, 5),
        "score_above_threshold": guide.score >= 40,
        "no_poly_t": "TTTT" not in guide.sequence,
        "specificity_screened": specificity is not None,
        "specificity_ok": specificity is None or specificity >= 50,
    }
    checks["overall_pass"] = all(
        (
            checks["length_ok"],
            checks["pam_ok"],
            checks["gc_in_acceptable_range"],
            checks["no_extreme_homopolymer"],
            checks["score_above_threshold"],
            checks["no_poly_t"],
            checks["specificity_ok"],
        )
    )
    return checks


def design_from_gene(
    gene_name: str,
    organism: str,
    application: str = "knockout",
    source: str = "ncbi",
    max_guides: int = 20,
    min_score: float = 30.0,
    genome_context: Optional[str] = None,
    max_mismatches: int = 3,
) -> Tuple[str, str, str, List[GuideRNA]]:
    """Fetch a representative transcript and design ranked guide candidates."""

    accession, description, sequence = fetch_sequence(
        gene_name,
        organism,
        source=source,
    )
    guides = design_guides(
        sequence,
        application=application,
        min_score=min_score,
        max_guides=max_guides,
        genome_context=genome_context,
        max_mismatches=max_mismatches,
    )
    return accession, description, sequence, guides

