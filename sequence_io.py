"""Strict sequence input handling; FASTA records are never silently joined."""
import re


def _normalize(raw: str) -> str:
    sequence = re.sub(r"[\s\d]", "", raw).upper().replace("U", "T")
    sequence = re.sub(r"[RYSWKMBDHVX]", "N", sequence)
    invalid = sorted(set(sequence) - set("ACGTN"))
    if invalid:
        raise ValueError("Sequence contains unsupported character(s): " + ", ".join(invalid))
    return sequence


def parse_reference(raw: str) -> dict[str, str]:
    """Read plain DNA or strict multi-FASTA, preserving record identity."""
    if not raw or not raw.strip():
        return {}
    raw = raw.lstrip("\ufeff")
    if not any(line.lstrip().startswith(">") for line in raw.splitlines()):
        sequence = _normalize(raw)
        return {"reference": sequence} if sequence else {}
    records, name, parts = {}, None, []

    def save():
        if name is not None:
            sequence = _normalize("".join(parts))
            if not sequence:
                raise ValueError(f"FASTA record '{name}' is empty.")
            records[name] = sequence

    for line in raw.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            save()
            name = line[1:].split()[0] if line[1:].split() else ""
            if not name:
                raise ValueError("Every FASTA record needs an identifier.")
            if name in records:
                raise ValueError(f"Duplicate FASTA identifier '{name}'; use unique names.")
            parts = []
        elif name is None:
            raise ValueError("Sequence occurs before the first FASTA header.")
        else:
            parts.append(line)
    save()
    return records


def clean_dna_sequence(raw: str) -> str:
    """Normalize a single target. Ambiguity stays N and is not a design base."""
    records = parse_reference(raw)
    if len(records) > 1:
        raise ValueError("Target input must contain one FASTA record; analyze targets separately.")
    return next(iter(records.values()), "")
