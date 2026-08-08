#!/usr/bin/env python3
"""
CRISPR Guide RNA Designer Core Module
=====================================
Fetches gene sequences from NCBI / Ensembl and designs SpCas9 gRNAs
for Knockout (NHEJ) or Knockdown (CRISPRi-style targeting).

Design principles based on:
- Jinek et al., Science 2012 (Cas9 + PAM NGG)
- Doench et al., Nature Biotechnology 2016 (Rule Set 2 features)
- Hsu et al., Nature Biotechnology 2013 (off-target considerations)
- Moreno-Mateos et al., Nature Methods 2015 (CRISPRscan)
- CRISPOR / CHOPCHOP community best practices
"""

from __future__ import annotations

import re
import time
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict, Any

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import requests

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
Entrez.email = "crispr.designer@example.com"  # NCBI requires a valid email
Entrez.tool = "CRISPR_gRNA_Designer"

PAM_PATTERN = re.compile(r"(?=([ATCG]GG))", re.IGNORECASE)  # NGG (lookahead)
PAM_PATTERN_REV = re.compile(r"(?=(CC[ATCG]))", re.IGNORECASE)  # CCN for reverse


@dataclass
class GuideRNA:
    """Represents a single designed guide RNA."""
    sequence: str                 # 20 nt spacer (DNA form)
    pam: str                      # 3 nt PAM
    strand: str                   # '+' or '-'
    start: int                    # 0-based start on input sequence
    end: int
    gc_content: float
    score: float                  # composite on-target score 0-100
    notes: List[str] = field(default_factory=list)
    application: str = "knockout"  # knockout | knockdown

    @property
    def full_target(self) -> str:
        return self.sequence + self.pam

    def to_dict(self) -> Dict[str, Any]:
        return {
            "Spacer (20 nt)": self.sequence,
            "PAM": self.pam,
            "Strand": self.strand,
            "Start": self.start + 1,  # 1-based for user
            "End": self.end,
            "GC%": round(self.gc_content, 1),
            "Score": round(self.score, 1),
            "Application": self.application,
            "Notes": "; ".join(self.notes) if self.notes else "OK",
        }


# ---------------------------------------------------------------------------
# Sequence fetching
# ---------------------------------------------------------------------------
def fetch_gene_ncbi(gene_name: str, organism: str, max_retries: int = 3) -> Tuple[str, str, str]:
    """
    Fetch gene sequence from NCBI using gene name + organism.
    Returns (accession, description, sequence).
    """
    query = f"{gene_name}[Gene Name] AND {organism}[Organism] AND alive[prop]"
    for attempt in range(max_retries):
        try:
            # Search Gene database
            handle = Entrez.esearch(db="gene", term=query, retmax=5)
            record = Entrez.read(handle)
            handle.close()
            if not record["IdList"]:
                # fallback broader search
                query2 = f"{gene_name} AND {organism}[Organism]"
                handle = Entrez.esearch(db="gene", term=query2, retmax=5)
                record = Entrez.read(handle)
                handle.close()
            if not record["IdList"]:
                raise ValueError(f"No gene found for '{gene_name}' in '{organism}' on NCBI.")

            gene_id = record["IdList"][0]
            # Get nucleotide links
            handle = Entrez.elink(dbfrom="gene", db="nuccore", id=gene_id, linkname="gene_nuccore_refseqrna")
            links = Entrez.read(handle)
            handle.close()

            nuc_ids = []
            if links and links[0]["LinkSetDb"]:
                nuc_ids = [link["Id"] for link in links[0]["LinkSetDb"][0]["Link"]]

            if not nuc_ids:
                # try genomic
                handle = Entrez.elink(dbfrom="gene", db="nuccore", id=gene_id)
                links = Entrez.read(handle)
                handle.close()
                if links and links[0]["LinkSetDb"]:
                    nuc_ids = [link["Id"] for link in links[0]["LinkSetDb"][0]["Link"]][:3]

            if not nuc_ids:
                raise ValueError("Could not resolve nucleotide accession for this gene.")

            # Prefer RefSeq NM_ or XR_ transcripts
            best_id = nuc_ids[0]
            handle = Entrez.efetch(db="nuccore", id=",".join(nuc_ids[:5]), rettype="gb", retmode="text")
            records = list(SeqIO.parse(handle, "genbank"))
            handle.close()

            # Prefer mRNA / complete CDS
            selected = None
            for rec in records:
                desc = rec.description.upper()
                if "MRNA" in desc or "TRANSCRIPT" in desc or any(f.type == "CDS" for f in rec.features):
                    selected = rec
                    break
            if selected is None:
                selected = records[0]

            seq = str(selected.seq).upper().replace("U", "T")
            return selected.id, selected.description, seq

        except Exception as e:
            if attempt == max_retries - 1:
                raise
            time.sleep(1 + attempt)
    raise RuntimeError("Failed after retries")


def fetch_gene_ensembl(gene_name: str, species: str = "homo_sapiens") -> Tuple[str, str, str]:
    """
    Fetch canonical transcript sequence from Ensembl REST API.
    species examples: homo_sapiens, mus_musculus, danio_rerio, saccharomyces_cerevisiae
    """
    server = "https://rest.ensembl.org"
    # lookup gene
    ext = f"/lookup/symbol/{species}/{gene_name}?expand=1"
    r = requests.get(server + ext, headers={"Content-Type": "application/json"}, timeout=30)
    if r.status_code == 400:
        # try with different species formatting
        species_map = {
            "human": "homo_sapiens",
            "mouse": "mus_musculus",
            "rat": "rattus_norvegicus",
            "zebrafish": "danio_rerio",
            "yeast": "saccharomyces_cerevisiae",
            "fly": "drosophila_melanogaster",
            "arabidopsis": "arabidopsis_thaliana",
        }
        sp = species_map.get(species.lower(), species.lower().replace(" ", "_"))
        ext = f"/lookup/symbol/{sp}/{gene_name}?expand=1"
        r = requests.get(server + ext, headers={"Content-Type": "application/json"}, timeout=30)

    if not r.ok:
        raise ValueError(f"Ensembl lookup failed ({r.status_code}): {r.text[:200]}")

    data = r.json()
    gene_id = data["id"]
    # get canonical transcript
    transcripts = data.get("Transcript", [])
    if not transcripts:
        raise ValueError("No transcripts found for this gene on Ensembl.")

    # prefer canonical
    canon = None
    for t in transcripts:
        if t.get("is_canonical"):
            canon = t
            break
    if canon is None:
        canon = transcripts[0]

    tx_id = canon["id"]
    # fetch sequence
    ext2 = f"/sequence/id/{tx_id}?type=cdna"
    r2 = requests.get(server + ext2, headers={"Content-Type": "text/plain"}, timeout=30)
    if not r2.ok:
        # genomic sequence fallback
        ext2 = f"/sequence/id/{gene_id}"
        r2 = requests.get(server + ext2, headers={"Content-Type": "text/plain"}, timeout=30)
    if not r2.ok:
        raise ValueError(f"Could not fetch sequence from Ensembl for {tx_id}")

    seq = r2.text.upper().replace("U", "T")
    desc = f"Ensembl {tx_id} | {data.get('display_name', gene_name)} | {data.get('description', '')}"
    return tx_id, desc, seq


def fetch_sequence(gene_name: str, organism: str, source: str = "ncbi") -> Tuple[str, str, str]:
    """Unified fetcher."""
    source = source.lower()
    if source == "ensembl":
        # map common names
        species_map = {
            "homo sapiens": "homo_sapiens",
            "human": "homo_sapiens",
            "mus musculus": "mus_musculus",
            "mouse": "mus_musculus",
            "rattus norvegicus": "rattus_norvegicus",
            "danio rerio": "danio_rerio",
            "zebrafish": "danio_rerio",
            "saccharomyces cerevisiae": "saccharomyces_cerevisiae",
            "arabidopsis thaliana": "arabidopsis_thaliana",
        }
        sp = species_map.get(organism.lower(), organism.lower().replace(" ", "_"))
        return fetch_gene_ensembl(gene_name, sp)
    else:
        return fetch_gene_ncbi(gene_name, organism)


# ---------------------------------------------------------------------------
# gRNA design algorithms (research-backed)
# ---------------------------------------------------------------------------
def _gc_content(seq: str) -> float:
    seq = seq.upper()
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def _has_homopolymer(seq: str, max_run: int = 4) -> bool:
    return bool(re.search(r"(A{%d,}|T{%d,}|G{%d,}|C{%d,})" % (max_run, max_run, max_run, max_run), seq.upper()))


def _doench_like_score(spacer: str, pam: str) -> float:
    """
    Simplified Doench Rule Set 2 inspired score (0-100).
    Real RS2 is a machine-learning model; here we implement the most
    important sequence features reported in Doench et al. 2016 and
    subsequent literature for educational / practical ranking.
    """
    s = spacer.upper()
    if len(s) != 20:
        return 0.0
    score = 50.0  # baseline

    # GC content optimal 40-70 %
    gc = _gc_content(s)
    if 40 <= gc <= 70:
        score += 15
    elif 30 <= gc < 40 or 70 < gc <= 80:
        score += 5
    else:
        score -= 15

    # Prefer G at position 20 (most important for activity)
    if s[19] == "G":
        score += 10
    elif s[19] == "C":
        score -= 8

    # Prefer A/G at position 19
    if s[18] in "AG":
        score += 5

    # Avoid T at position 20 (U6 terminator related)
    if s[19] == "T":
        score -= 12

    # Prefer C at position 16-18 region (approximate)
    if s[15] == "C":
        score += 3
    if s[17] == "C":
        score += 3

    # PAM: CGG preferred over others
    if pam.upper() == "CGG":
        score += 5
    elif pam.upper()[0] == "T":
        score -= 5

    # penalize poly-G / poly-T
    if "GGGG" in s or "TTTT" in s:
        score -= 20
    if "GGG" in s or "TTT" in s:
        score -= 8

    # self-complementarity rough penalty (simple)
    rev = str(Seq(s).reverse_complement())
    matches = sum(1 for a, b in zip(s, rev) if a == b)
    if matches > 12:
        score -= 10

    return max(0.0, min(100.0, score))


def design_guides(
    sequence: str,
    application: str = "knockout",
    pam: str = "NGG",
    spacer_len: int = 20,
    min_score: float = 30.0,
    max_guides: int = 50,
    prefer_5prime: bool = True,
    genome_context: Optional[str] = None,
) -> List[GuideRNA]:
    """
    Scan sequence for NGG PAMs and design ranked gRNAs.

    application:
      - "knockout"  : prefer early region of sequence (proxy for early exons)
      - "knockdown" : prefer first ~200-500 bp (proxy for near TSS / promoter proximal)
    """
    sequence = sequence.upper().replace("U", "T")
    # clean non-ATGC
    sequence = re.sub(r"[^ATGC]", "N", sequence)
    guides: List[GuideRNA] = []

    # Forward strand NGG
    for m in PAM_PATTERN.finditer(sequence):
        pam_start = m.start()
        spacer_start = pam_start - spacer_len
        if spacer_start < 0:
            continue
        spacer = sequence[spacer_start:pam_start]
        if "N" in spacer:
            continue
        pam_seq = sequence[pam_start:pam_start + 3]
        if len(pam_seq) < 3:
            continue

        notes = []
        gc = _gc_content(spacer)
        if gc < 20 or gc > 80:
            notes.append("Extreme GC")
        if _has_homopolymer(spacer, 4):
            notes.append("Homopolymer ≥4")
        if spacer.endswith("T"):
            notes.append("Ends with T (U6 risk)")

        score = _doench_like_score(spacer, pam_seq)

        # application bias
        if application.lower() == "knockdown":
            # boost guides in first 400 bp
            if spacer_start < 400:
                score += 15
                notes.append("Near 5' (KD preferred)")
            else:
                score -= 5
        else:  # knockout – mild preference for earlier positions
            if prefer_5prime and spacer_start < len(sequence) * 0.3:
                score += 8
                notes.append("Early region (KO preferred)")

        # === OFF-TARGET ANALYSIS ADDED (NEW) ===
        if genome_context:
            off_count, off_score, off_list = calculate_offtarget_score(
                spacer, pam_seq, genome_context, max_mismatches=3
            )
            notes.append(f"OT: {off_count} off-target(s) | Score: {off_score}")
            notes.extend(off_list[:3])  # only show top 3 off-targets
        else:
            off_count, off_score = 0, 100.0

        if score >= min_score:
            guides.append(
                GuideRNA(
                    sequence=spacer,
                    pam=pam_seq,
                    strand="+",
                    start=spacer_start,
                    end=pam_start + 3,
                    gc_content=gc,
                    score=score,
                    notes=notes,
                    application=application.lower(),
                )
            )

    # Reverse strand (CCN)
    for m in PAM_PATTERN_REV.finditer(sequence):
        pam_start = m.start()
        # spacer is downstream of CCN on reverse
        spacer_end = pam_start + 3 + spacer_len
        if spacer_end > len(sequence):
            continue
        spacer = sequence[pam_start + 3 : spacer_end]
        if "N" in spacer:
            continue
        # reverse complement to get the actual spacer that will be in the gRNA
        spacer_rc = str(Seq(spacer).reverse_complement())
        pam_seq = sequence[pam_start:pam_start + 3]
        # PAM on reverse is CCN; the effective PAM recognized is NGG on the opposite strand
        effective_pam = str(Seq(pam_seq).reverse_complement())

        notes = []
        gc = _gc_content(spacer_rc)
        if gc < 20 or gc > 80:
            notes.append("Extreme GC")
        if _has_homopolymer(spacer_rc, 4):
            notes.append("Homopolymer ≥4")
        if spacer_rc.endswith("T"):
            notes.append("Ends with T (U6 risk)")

        score = _doench_like_score(spacer_rc, effective_pam)

        if application.lower() == "knockdown":
            if pam_start < 400:
                score += 15
                notes.append("Near 5' (KD preferred)")
            else:
                score -= 5
        else:
            if prefer_5prime and pam_start < len(sequence) * 0.3:
                score += 8
                notes.append("Early region (KO preferred)")

        # === OFF-TARGET ANALYSIS ADDED (NEW) ===
        if genome_context:
            off_count, off_score, off_list = calculate_offtarget_score(
                spacer_rc, effective_pam, genome_context, max_mismatches=3
            )
            notes.append(f"OT: {off_count} off-target(s) | Score: {off_score}")
            notes.extend(off_list[:3])  # only show top 3 off-targets
        else:
            off_count, off_score = 0, 100.0

        if score >= min_score:
            guides.append(
                GuideRNA(
                    sequence=spacer_rc,
                    pam=effective_pam,
                    strand="-",
                    start=pam_start,
                    end=spacer_end,
                    gc_content=gc,
                    score=score,
                    notes=notes,
                    application=application.lower(),
                )
            )

    # rank
    guides.sort(key=lambda g: g.score, reverse=True)
    return guides[:max_guides]


def validate_guide(guide: GuideRNA, genome_context: Optional[str] = None) -> Dict[str, Any]:
    """
    Basic validation checks that a user / downstream AI can use to judge whether a guide is reasonable.
    """
    checks = {
        "length_ok": len(guide.sequence) == 20,
        "pam_ok": guide.pam.upper().endswith("GG"),
        "gc_in_range": 30 <= guide.gc_content <= 80,
        "no_extreme_homopolymer": not _has_homopolymer(guide.sequence, 5),
        "score_above_threshold": guide.score >= 40,
        "no_polyT_end": not guide.sequence.endswith("TT"),
    }
    checks["overall_pass"] = all(
        [
            checks["length_ok"],
            checks["pam_ok"],
            checks["gc_in_range"],
            checks["no_extreme_homopolymer"],
            checks["score_above_threshold"],
        ]
    )
    return checks


# ---------------------------------------------------------------------------
# Convenience high-level API
# ---------------------------------------------------------------------------
def design_from_gene(
    gene_name: str,
    organism: str,
    application: str = "knockout",
    source: str = "ncbi",
    max_guides: int = 20,
) -> Tuple[str, str, str, List[GuideRNA]]:
    """
    End-to-end: fetch → design.
    Returns (accession, description, sequence, list_of_guides)
    """
    acc, desc, seq = fetch_sequence(gene_name, organism, source=source)
    guides = design_guides(seq, application=application, max_guides=max_guides)
    return acc, desc, seq, guides


# ---------------------------------------------------------------------------
# Off-target analysis (added in your update)
# ---------------------------------------------------------------------------
def calculate_offtarget_score(
    spacer: str,
    pam: str,
    whole_genome: Optional[str] = None,
    max_mismatches: int = 3,
) -> Tuple[int, float, List[str]]:
    """
    Calculates off-target statistics by scanning `whole_genome` for windows
    that match the full target (spacer + PAM) within `max_mismatches`
    substitutions (simple Hamming-distance search — no indels).

    This intentionally stays dependency-free for portability. For production
    / genome-scale use, swap this out for a proper aligner-backed search:
    - Cas-OFFinder (faster, better — the standard tool for this)
    - CRISPOR / CHOPCHOP API
    - Bowtie2 seed-and-extend

    Returns (num_offtargets, specificity_score_0_100, list_of_hit_descriptions).
    """
    if not whole_genome:
        return 0, 0.0, ["Genome not provided"]

    full_target = (spacer + pam.upper()).upper()
    genome = whole_genome.upper()
    target_len = len(full_target)

    total_offtargets = 0
    off_target_list: List[str] = []

    i = 0
    last_pos = len(genome) - target_len
    while i <= last_pos:
        window = genome[i : i + target_len]
        mismatches = sum(1 for a, b in zip(window, full_target) if a != b)
        if mismatches <= max_mismatches:
            total_offtargets += 1
            plural = "es" if mismatches != 1 else ""
            off_target_list.append(f"Position {i} : {window} ({mismatches} mismatch{plural})")
            # Jump past this hit so one repeated locus isn't counted once
            # per overlapping shift of the sliding window.
            i += target_len
        else:
            i += 1

    # Specificity score: more / closer off-targets pull the score down.
    if total_offtargets == 0:
        off_target_score = 100.0
    else:
        off_target_score = max(0.0, 100.0 - (total_offtargets * (max_mismatches / 3.0)))

    return total_offtargets, round(off_target_score, 1), off_target_list