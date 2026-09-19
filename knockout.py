"""Genomic knockout discovery using NCBI annotated coding regions.

All discovery uses a contiguous chromosome slice, never concatenated CDS/cDNA.
Coordinates inside this module are zero-based half-open relative to that slice.
"""
from dataclasses import dataclass
from io import StringIO
from urllib.parse import quote

from Bio import SeqIO
from Bio.SeqFeature import ExactPosition

from grna_designer import NCBI_EUTILS, NCBI_EMAIL, NCBI_TOOL, clean_dna_sequence, design_guides, screen_guides
from network import get


@dataclass(frozen=True)
class CodingIsoform:
    protein_id: str
    transcript_ids: tuple[str, ...]
    strand: int
    segments: tuple[tuple[int, int], ...]

    @property
    def coding_length(self):
        return sum(end - start for start, end in self.segments)


@dataclass(frozen=True)
class KnockoutContext:
    gene: str
    organism: str
    gene_id: str
    accession: str
    chromosome: str
    region_start: int  # 1-based genomic coordinate of first input base
    sequence: str  # always forward genomic orientation
    isoforms: tuple[CodingIsoform, ...]
    selected_protein: str
    selection_rule: str
    description: str


def _get(path, **params):
    response = get(f'{NCBI_EUTILS}/{path}', params={'tool': NCBI_TOOL, 'email': NCBI_EMAIL, **params}, timeout=35)
    response.raise_for_status()
    return response


def context_from_genbank(text, gene, organism, gene_id, accession, chromosome, region_start, protein_id=''):
    record = SeqIO.read(StringIO(text), 'genbank')
    if record.id != accession:
        raise ValueError('Fetched genomic accession differs from the requested reference.')
    sequence = clean_dna_sequence(str(record.seq))
    features = [f for f in record.features if f.qualifiers.get('gene', [''])[0].upper() == gene.upper()]
    transcripts = [f for f in features if f.type == 'mRNA' and f.location is not None]
    isoforms = []
    for feature in features:
        if feature.type != 'CDS' or not feature.qualifiers.get('protein_id') or feature.location is None:
            continue
        parts = feature.location.parts
        if any(not isinstance(p.start, ExactPosition) or not isinstance(p.end, ExactPosition) or p.ref for p in parts):
            continue  # Do not call a partial/remote feature a complete coding annotation.
        strand = feature.location.strand
        if strand not in (-1, 1):
            continue
        segments = tuple(sorted((int(p.start), int(p.end)) for p in parts))
        if any(a < 0 or b > len(sequence) or a >= b for a, b in segments):
            raise ValueError('Coding annotation is outside the fetched genomic region.')
        tids = []
        for tx in transcripts:
            if tx.location.strand == strand and all(any(int(p.start) <= a and b <= int(p.end) for p in tx.location.parts) for a, b in segments):
                tids.extend(tx.qualifiers.get('transcript_id', []))
        isoforms.append(CodingIsoform(feature.qualifiers['protein_id'][0], tuple(sorted(set(tids))), strand, segments))
    if not isoforms:
        raise ValueError('No complete protein-coding annotation for this gene in the fetched genomic record. Noncoding genes need a different design strategy.')
    # Deterministic representative; explicitly not claimed to be MANE or canonical.
    isoforms.sort(key=lambda x: (not x.protein_id.startswith('NP_'), -x.coding_length, x.protein_id))
    selected = next((x for x in isoforms if x.protein_id == protein_id or protein_id in x.transcript_ids), None) if protein_id else isoforms[0]
    if selected is None:
        choices = ', '.join(x.protein_id for x in isoforms)
        raise ValueError(f'Isoform accession is not in this gene annotation. Available proteins: {choices}')
    return KnockoutContext(
        gene, organism, str(gene_id), accession, chromosome, region_start, sequence,
        tuple(isoforms), selected.protein_id,
        'User-selected isoform' if protein_id else 'Reviewed RefSeq protein preferred, then longest CDS; not a canonical-transcript claim',
        f'NCBI genomic coding-region design | {gene} | {accession}:{region_start}-{region_start + len(sequence) - 1}',
    )


def fetch_ncbi_knockout_context(gene, organism, protein_id='', max_bp=2_000_000):
    if not gene.strip() or not organism.strip() or any(c in gene + organism for c in '[]"'):
        raise ValueError('Enter a gene symbol and organism without query syntax.')
    result = _get('esearch.fcgi', db='gene', term=f'{gene}[Gene Name] AND {organism}[Organism] AND alive[prop]', retmax=20, retmode='json').json()['esearchresult']
    ids = result.get('idlist', [])
    if not ids:
        raise ValueError(f'No NCBI gene found for {gene} in {organism}.')
    if int(result.get('count', len(ids))) > len(ids):
        raise ValueError('Gene query is too broad; use the exact gene symbol and organism.')
    summary = _get('esummary.fcgi', db='gene', id=','.join(ids), retmode='json').json().get('result', {})
    matches = [summary[i] for i in ids if summary.get(i, {}).get('name', '').upper() == gene.upper()]
    if len(matches) != 1:
        raise ValueError('Gene name is ambiguous or is an alias; use the exact NCBI gene symbol.')
    row = matches[0]
    loci = row.get('genomicinfo', [])
    if len(loci) != 1:
        raise ValueError('Gene has zero or multiple genomic mappings; select a reference region explicitly.')
    locus = loci[0]
    # Four bases on either side permit full 30-mer context when possible.
    start = max(1, min(int(locus['chrstart']), int(locus['chrstop'])) + 1 - 4)
    end = max(int(locus['chrstart']), int(locus['chrstop'])) + 1 + 4
    if end - start + 1 > max_bp:
        raise ValueError(f'Genomic gene region exceeds {max_bp:,} bp; use a smaller annotated region.')
    accession = locus['chraccver']
    text = _get('efetch.fcgi', db='nuccore', id=accession, seq_start=start, seq_stop=end, strand=1, rettype='gbwithparts', retmode='text').text
    context = context_from_genbank(text, row['name'], row.get('organism', {}).get('scientificname', organism), row['uid'], accession, str(locus['chrloc']), start, protein_id)
    if len(context.sequence) != end - start + 1:
        raise ValueError('Fetched genomic region length differs from requested coordinates.')
    return context


def design_knockout_guides(context, max_guides=20, min_score=30, genome_context=None, max_mismatches=3, intended_region=None):
    if type(max_guides) is not int or max_guides < 1:
        raise ValueError('Maximum guides must be a positive integer.')
    selected = next(x for x in context.isoforms if x.protein_id == context.selected_protein)
    candidates = design_guides(context.sequence, min_score=min_score, max_guides=50_000, prefer_5prime=False)
    guides = []
    ordered = sorted(selected.segments, reverse=selected.strand == -1)
    for guide in candidates:
        cut = guide.cut_boundary
        segment = next(((a, b) for a, b in selected.segments if a < cut < b), None)
        if segment is None:
            continue
        before = 0
        for a, b in ordered:
            if (a, b) == segment:
                before += cut - a if selected.strand == 1 else b - cut
                break
            before += b - a
        supported = [x.protein_id for x in context.isoforms if any(a < cut < b for a, b in x.segments)]
        guide.annotation = {
            'Reference accession': context.accession, 'Chromosome': context.chromosome,
            'Genomic spacer start': context.region_start + guide.spacer_start,
            'Genomic spacer end': context.region_start + guide.spacer_end - 1,
            'Genomic cut after base': context.region_start + cut - 1,
            'Genomic strand': guide.strand, 'Selected protein': selected.protein_id,
            'Compatible transcript annotations': '; '.join(selected.transcript_ids),
            'CDS position percent': round(100 * before / selected.coding_length, 2),
            'Coding isoforms at cut': len(supported), 'Annotated coding isoforms': len(context.isoforms),
        }
        guide.notes.append('Predicted cut inside annotated CDS; editing outcome unmodeled')
        guides.append(guide)
    if genome_context is not None:
        screen_guides(guides, genome_context, max_mismatches, intended_region=intended_region, target_sequence=context.sequence)
    guides.sort(key=lambda g: (-g.score, -(g.specificity_score if g.specificity_score is not None else -1), g.start, g.strand))
    return guides[:max_guides]
