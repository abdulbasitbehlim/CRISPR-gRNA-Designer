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


def normalize_accession(value: str) -> str:
    accession = value.strip()
    if not accession:
        raise ValueError("Enter an accession or stable ID.")
    if any(ch.isspace() for ch in accession):
        raise ValueError("Accession IDs cannot contain spaces.")
    return accession


def detect_database(accession: str) -> str:
    """Return 'ensembl' for ENS stable IDs; use NCBI otherwise."""
    accession = normalize_accession(accession).upper()
    return "ensembl" if accession.startswith("ENS") else "ncbi"


def _ncbi_get(path: str, params: dict) -> requests.Response:
    response = requests.get(
        f"{NCBI_EUTILS}/{path}",
        params={"tool": NCBI_TOOL, "email": NCBI_EMAIL, **params},
        timeout=35,
    )
    response.raise_for_status()
    return response


def fetch_ncbi_accession(accession: str, max_bp: int = MAX_ONLINE_RECORD_BP) -> AccessionRecord:
    accession = normalize_accession(accession)

    search = _ncbi_get(
        "esearch.fcgi",
        {"db": "nuccore", "term": f"{accession}[Accession]", "retmax": 5, "retmode": "json"},
    )
    ids = search.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        # Version-less RefSeq/GenBank accessions are commonly entered by users.
        base = accession.split(".")[0]
        search = _ncbi_get(
            "esearch.fcgi",
            {"db": "nuccore", "term": base, "retmax": 5, "retmode": "json"},
        )
        ids = search.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        raise ValueError(f"NCBI Nucleotide could not find '{accession}'.")

    summary = _ncbi_get(
        "esummary.fcgi",
        {"db": "nuccore", "id": ids[0], "retmode": "json"},
    ).json()
    row = summary.get("result", {}).get(str(ids[0]), {})
    length = int(row.get("slen") or 0)
    if length and length > max_bp:
        raise ValueError(
            f"The NCBI record is {length:,} bp. The hosted accession limit is {max_bp:,} bp; "
            "use a gene/transcript accession or paste a smaller target region."
        )

    fasta = _ncbi_get(
        "efetch.fcgi",
        {"db": "nuccore", "id": ids[0], "rettype": "fasta", "retmode": "text"},
    ).text
    records = list(SeqIO.parse(StringIO(fasta), "fasta"))
    if not records:
        raise ValueError(f"NCBI returned no nucleotide sequence for '{accession}'.")
    record = records[0]
    sequence = clean_dna_sequence(str(record.seq))
    if not sequence:
        raise ValueError("NCBI returned an empty nucleotide sequence.")
    if len(sequence) > max_bp:
        raise ValueError(f"Fetched NCBI sequence exceeds the {max_bp:,} bp hosted limit.")

    resolved = record.id
    title = str(row.get("title") or record.description or resolved)
    return AccessionRecord(
        query=accession,
        accession=resolved,
        database="NCBI Nucleotide",
        description=title,
        sequence=sequence,
        record_url=f"https://www.ncbi.nlm.nih.gov/nuccore/{resolved}",
        object_type="nucleotide",
    )


def _ensembl_lookup(accession: str) -> dict:
    response = requests.get(
        f"{ENSEMBL_REST}/lookup/id/{accession}?expand=1",
        headers={"Accept": "application/json", "User-Agent": NCBI_TOOL},
        timeout=35,
    )
    if not response.ok:
        raise ValueError(f"Ensembl could not resolve '{accession}' ({response.status_code}).")
    return response.json()


def _ensembl_sequence(stable_id: str, sequence_type: Optional[str] = None) -> str:
    params = f"?type={sequence_type}" if sequence_type else ""
    response = requests.get(
        f"{ENSEMBL_REST}/sequence/id/{stable_id}{params}",
        headers={"Accept": "text/plain", "User-Agent": NCBI_TOOL},
        timeout=35,
    )
    if not response.ok:
        raise ValueError(f"Ensembl could not fetch sequence for '{stable_id}' ({response.status_code}).")
    return clean_dna_sequence(response.text)


def fetch_ensembl_accession(accession: str, max_bp: int = MAX_ONLINE_RECORD_BP) -> AccessionRecord:
    accession = normalize_accession(accession)
    data = _ensembl_lookup(accession)
    object_type = str(data.get("object_type") or "").lower()
    resolved = str(data.get("id") or accession)
    description = str(data.get("description") or data.get("display_name") or resolved)

    if object_type == "gene":
        transcripts = data.get("Transcript", []) or []
        canonical = str(data.get("canonical_transcript") or "").split(".")[0]
        transcript = next(
            (
                tx for tx in transcripts
                if tx.get("is_canonical") or str(tx.get("id", "")).split(".")[0] == canonical
            ),
            None,
        )
        if transcript is None and transcripts:
            transcript = transcripts[0]
        if transcript is None:
            raise ValueError(f"Ensembl gene '{accession}' has no transcript sequence to design from.")
        resolved = str(transcript.get("id"))
        sequence = _ensembl_sequence(resolved, "cdna")
        description = f"{description} | canonical/representative transcript {resolved}"
        object_type = "gene → transcript cDNA"
    elif object_type == "transcript":
        sequence = _ensembl_sequence(resolved, "cdna")
        object_type = "transcript cDNA"
    elif object_type in {"exon", "translation"}:
        sequence = _ensembl_sequence(resolved)
    else:
        sequence = _ensembl_sequence(resolved)
        object_type = object_type or "stable ID"

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
