# ============================================================================
# ACCESSION LOOKUP
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Handles accession-based sequence lookup and retrieval used by the CRISPR application.
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
# - class: AccessionRecord
# - function: normalize_accession
# - function: detect_database
# - function: _ncbi_get
# - function: fetch_ncbi_accession
# - function: _ensembl_lookup
# - function: _ensembl_sequence
# - function: fetch_ensembl_accession
# - function: fetch_accession
# ============================================================================

"""Online accession lookup for CRISPR Studio.

Supports nucleotide records from NCBI Nucleotide/RefSeq and stable IDs from
Ensembl. The module returns DNA plus provenance so the Streamlit UI can show
exactly where a target sequence came from.
"""
from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
import re
from typing import Optional

import requests
from network import get
from Bio import SeqIO

from grna_designer import NCBI_EMAIL, NCBI_EUTILS, NCBI_TOOL, clean_dna_sequence

ENSEMBL_REST = "https://rest.ensembl.org"
MAX_ONLINE_RECORD_BP = 2_000_000


@dataclass(frozen=True)
class AccessionRecord:
    query: str
    accession: str
    database: str
    description: str
    sequence: str
    record_url: str
    object_type: str = "nucleotide"

    @property
    def length(self) -> int:
        return len(self.sequence)



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: normalize_accession
# ----------------------------------------------------------------------------
def normalize_accession(value: str) -> str:
    accession = value.strip()
    if not accession:
        raise ValueError("Enter an accession or stable ID.")
    if not re.fullmatch(r"[A-Za-z][A-Za-z_]*[0-9]+(?:\.[0-9]+)?", accession):
        raise ValueError("Enter one valid nucleotide accession or Ensembl stable ID without spaces or query syntax.")
    return accession



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: detect_database
# ----------------------------------------------------------------------------
def detect_database(accession: str) -> str:
    """Return 'ensembl' for ENS stable IDs; use NCBI otherwise."""
    accession = normalize_accession(accession).upper()
    return "ensembl" if accession.startswith("ENS") else "ncbi"



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _ncbi_get
# ----------------------------------------------------------------------------
def _ncbi_get(path: str, params: dict) -> requests.Response:
    response = get(
        f"{NCBI_EUTILS}/{path}",
        params={"tool": NCBI_TOOL, "email": NCBI_EMAIL, **params},
        timeout=35,
    )
    response.raise_for_status()
    return response



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: fetch_ncbi_accession
# ----------------------------------------------------------------------------
def fetch_ncbi_accession(accession: str, max_bp: int = MAX_ONLINE_RECORD_BP) -> AccessionRecord:
    accession = normalize_accession(accession)

    search = _ncbi_get(
        "esearch.fcgi",
        {"db": "nuccore", "term": f"{accession}[Accession]", "retmax": 5, "retmode": "json"},
    )
    ids = search.json().get("esearchresult", {}).get("idlist", [])
    if len(ids) > 1:
        raise ValueError("Accession query matched multiple records; specify the exact version.")
    if not ids:
        raise ValueError(f"NCBI Nucleotide could not find '{accession}'.")

    summary = _ncbi_get(
        "esummary.fcgi",
        {"db": "nuccore", "id": ids[0], "retmode": "json"},
    ).json()
    row = summary.get("result", {}).get(str(ids[0]), {})
    biomol = str(row.get('biomol', '')).lower()
    if accession.upper().startswith(('NM_', 'XM_', 'NR_', 'XR_')) or 'rna' in biomol or biomol == 'cdna':
        raise ValueError('RNA/cDNA accessions cannot define contiguous genomic knockout targets. Use NCBI Gene lookup or a genomic DNA accession.')
    if biomol != 'genomic':
        raise ValueError('NCBI record is not verified as genomic DNA. Use Gene lookup or supply a verified genomic region.')
    length = int(row.get("slen") or 0)
    if length and length > max_bp:
        raise ValueError(
            f"The NCBI record is {length:,} bp. The hosted accession limit is {max_bp:,} bp; "
            "use Gene lookup or paste a smaller genomic target region."
        )

    fasta = _ncbi_get(
        "efetch.fcgi",
        {"db": "nuccore", "id": ids[0], "rettype": "fasta", "retmode": "text"},
    ).text
    records = list(SeqIO.parse(StringIO(fasta), "fasta"))
    if len(records) != 1:
        raise ValueError(f"NCBI must return exactly one nucleotide record for '{accession}'.")
    record = records[0]
    sequence = clean_dna_sequence(str(record.seq))
    if not sequence:
        raise ValueError("NCBI returned an empty nucleotide sequence.")
    if len(sequence) > max_bp:
        raise ValueError(f"Fetched NCBI sequence exceeds the {max_bp:,} bp hosted limit.")

    resolved = record.id
    expected = accession.upper()
    if (resolved.upper() != expected if '.' in expected else resolved.split('.')[0].upper() != expected):
        raise ValueError(f"Requested accession {accession} does not match returned record {resolved}; no version fallback was used.")
    title = str(row.get("title") or record.description or resolved)
    return AccessionRecord(
        query=accession,
        accession=resolved,
        database="NCBI Nucleotide",
        description=title,
        sequence=sequence,
        record_url=f"https://www.ncbi.nlm.nih.gov/nuccore/{resolved}",
        object_type="genomic DNA (coding annotation not applied)",
    )



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _ensembl_lookup
# ----------------------------------------------------------------------------
def _ensembl_lookup(accession: str) -> dict:
    response = get(
        f"{ENSEMBL_REST}/lookup/id/{accession}?expand=1",
        headers={"Accept": "application/json", "User-Agent": NCBI_TOOL},
        timeout=35,
    )
    if not response.ok:
        raise ValueError(f"Ensembl could not resolve '{accession}' ({response.status_code}).")
    return response.json()



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _ensembl_sequence
# ----------------------------------------------------------------------------
def _ensembl_sequence(stable_id: str, sequence_type: Optional[str] = None) -> str:
    params = f"?type={sequence_type}" if sequence_type else ""
    response = get(
        f"{ENSEMBL_REST}/sequence/id/{stable_id}{params}",
        headers={"Accept": "text/plain", "User-Agent": NCBI_TOOL},
        timeout=35,
    )
    if not response.ok:
        raise ValueError(f"Ensembl could not fetch sequence for '{stable_id}' ({response.status_code}).")
    return clean_dna_sequence(response.text)



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: fetch_ensembl_accession
# ----------------------------------------------------------------------------
def fetch_ensembl_accession(accession: str, max_bp: int = MAX_ONLINE_RECORD_BP) -> AccessionRecord:
    accession = normalize_accession(accession)
    base, _, requested_version = accession.partition('.')
    data = _ensembl_lookup(base)
    if requested_version and str(data.get('version', '')) != requested_version:
        raise ValueError('The requested Ensembl version does not match the current record.')
    object_type = str(data.get("object_type") or "").lower()
    resolved = str(data.get("id") or accession)
    description = str(data.get("description") or data.get("display_name") or resolved)

    if object_type not in {'gene', 'transcript', 'exon'}:
        raise ValueError("Only nucleotide gene, transcript or exon IDs are accepted; protein translation IDs cannot be used as DNA.")
    expected_length = int(data['end']) - int(data['start']) + 1
    if expected_length < 1 or expected_length > max_bp:
        raise ValueError(f'Genomic interval exceeds the {max_bp:,} bp limit or has invalid coordinates.')
    sequence = _ensembl_sequence(resolved, 'genomic')
    if len(sequence) != expected_length:
        raise ValueError('Ensembl genomic sequence length does not match the annotated interval.')
    description = f"{description} | genomic {data.get('assembly_name', 'unknown assembly')} {data.get('seq_region_name', '')}:{data['start']}-{data['end']} strand {data.get('strand', '')}"
    object_type = f'{object_type} genomic DNA (coding annotation not applied)'

    if not sequence:
        raise ValueError("Ensembl returned an empty DNA sequence.")
    if len(sequence) > max_bp:
        raise ValueError(
            f"The Ensembl sequence is {len(sequence):,} bp. The hosted accession limit is {max_bp:,} bp."
        )

    return AccessionRecord(
        query=accession,
        accession=resolved,
        database="Ensembl",
        description=description,
        sequence=sequence,
        record_url=f"https://www.ensembl.org/id/{accession}",
        object_type=object_type,
    )



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: fetch_accession
# ----------------------------------------------------------------------------
def fetch_accession(accession: str, database: str = "Auto") -> AccessionRecord:
    """Fetch an online nucleotide target using Auto, NCBI, or Ensembl routing."""
    accession = normalize_accession(accession)
    choice = database.strip().lower()
    if choice == "auto":
        choice = detect_database(accession)
    if choice.startswith("ncbi"):
        return fetch_ncbi_accession(accession)
    if choice.startswith("ensembl"):
        return fetch_ensembl_accession(accession)
    raise ValueError("Database must be Auto, NCBI Nucleotide, or Ensembl.")
