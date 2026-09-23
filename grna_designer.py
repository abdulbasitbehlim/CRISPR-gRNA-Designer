#!/usr/bin/env python3

# ============================================================================
# GRNA DESIGNER
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Contains the main scientific logic for finding, checking and ranking gRNA candidates.
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
# - class: TargetLocus
# - function: intended_loci
# - class: OffTargetHit
# - class: OffTargetReport
# - class: GuideRNA
# - function: fetch_gene_ensembl
# - function: fetch_gene_ncbi
# - function: fetch_sequence
# - function: _gc_content
# - function: _has_homopolymer
# - function: score_breakdown
# - function: _doench_like_score
# - function: mit_offtarget_score
# - function: mit_specificity
# - function: _optional_cfd
# - function: _risk_label
# - function: iter_pam_sites
# - class: ReferenceIndex
# - function: analyze_offtargets
# - function: screen_guides
# - function: calculate_offtarget_score
# - function: _context30
# - function: doench_rs2_score
# - function: design_guides
# - function: guidescan2_available
# - plus 3 additional helper functions/classes
# ============================================================================

"""CRISPR Studio v3 core: SpCas9 design, MIT specificity, local multi-FASTA screening,
optional Doench Rule Set 2, bundled CFD weights, and GuideScan2 whole-genome integration.
"""
from __future__ import annotations
import csv, os, re, shutil, subprocess, tempfile, time, math
import numpy as np
from sequence_io import clean_dna_sequence, parse_reference
from models import cfd_score, rs2_score
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, NamedTuple
import requests
from network import get
from Bio import SeqIO
from Bio.Seq import Seq

NCBI_EMAIL=os.getenv('NCBI_EMAIL','crispr.grna.designer@example.com')
NCBI_TOOL='CRISPR_gRNA_Designer'
NCBI_EUTILS='https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
PAM_PATTERN=re.compile(r'(?=([ATCG]GG))',re.I)
PAM_PATTERN_REV=re.compile(r'(?=(CC[ATCG]))',re.I)
DNA_ALPHABET=frozenset('ATGCN')
IUPAC_AMBIGUOUS=re.compile(r'[RYSWKMBDHVX]')
MIT_WEIGHTS=(0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583)


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: TargetLocus
# ----------------------------------------------------------------------------
class TargetLocus(NamedTuple):
    """Explicit site: contig, zero-based left edge of spacer+PAM, strand."""
    contig: str
    start: int
    strand: str



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: intended_loci
# ----------------------------------------------------------------------------
def intended_loci(guides, sequence, index, region):
    """Map a user-declared genomic region after verifying its full sequence."""
    contig, offset, orientation = region
    if type(offset) is not int or offset < 0 or orientation not in {'+', '-'}:
        raise ValueError('Intended region needs a nonnegative integer offset and + or - strand.')
    reference = index.contigs.get(contig, '')
    expected = sequence if orientation == '+' else str(Seq(sequence).reverse_complement())
    if not expected or reference[offset:offset + len(sequence)] != expected:
        raise ValueError('Declared target region does not exactly match the full input sequence at this contig, position and orientation.')
    return {(g.start, g.strand): TargetLocus(contig,
            offset + (g.start if orientation == '+' else len(sequence) - g.end),
            g.strand if orientation == '+' else ('-' if g.strand == '+' else '+'))
            for g in guides}

@dataclass(frozen=True)
class OffTargetHit:
    start:int; strand:str; sequence:str; pam:str; mismatches:int
    mismatch_positions:Tuple[int,...]; seed_mismatches:int; risk:str
    contig:str='reference'; mit_score:float=0.0; cfd_score:Optional[float]=None
    def to_dict(self)->Dict[str,Any]:
        return {'Contig':self.contig,'Position':self.start+1,'Strand':self.strand,
        'Candidate spacer':self.sequence,'PAM':self.pam,'Mismatches':self.mismatches,
        'Mismatch positions':', '.join(map(str,self.mismatch_positions)) or 'Exact',
        'Seed mismatches':self.seed_mismatches,'MIT off-target':round(self.mit_score*100,2),
        'CFD off-target':None if self.cfd_score is None else round(self.cfd_score*100,2),'Risk':self.risk}

@dataclass(frozen=True)
class OffTargetReport:
    hits:Tuple[OffTargetHit,...]; specificity_score:float; pam_sites_scanned:int; on_target_excluded:bool
    cfd_specificity:Optional[float]=None; reference_contigs:int=1
    total_hits:int=0; exact_matches:int=0; truncated:bool=False
    intended_target_status:str='not_declared'
    ambiguous_bases:int=0
    screened_mismatch_radius:int=0
    risk_counts:Dict[str,int]=field(default_factory=lambda: {'Critical':0, 'High':0, 'Moderate':0, 'Low':0})
    max_mit_risk:float=0.0
    max_cfd_risk:float=0.0
    @property
    def mit_specificity(self): return self.specificity_score

@dataclass
class GuideRNA:
    sequence:str; pam:str; strand:str; start:int; end:int; gc_content:float; score:float
    notes:List[str]=field(default_factory=list); application:str='knockout'
    specificity_score:Optional[float]=None; off_target_count:Optional[int]=None
    off_target_details:List[Dict[str,Any]]=field(default_factory=list)
    doench_score:Optional[float]=None; cfd_specificity:Optional[float]=None
    scoring_method:str='Heuristic'; genome_specificity:Optional[float]=None
    genome_backend:Optional[str]=None
    screen_status:str='not_screened'
    exact_match_count:Optional[int]=None
    intended_target:Optional[Tuple[str,int,str]]=None
    reference_ambiguous_bases:int=0
    screened_mismatch_radius:Optional[int]=None
    risk_counts:Dict[str,int]=field(default_factory=lambda: {'Critical':0, 'High':0, 'Moderate':0, 'Low':0})
    max_mit_risk:Optional[float]=None
    max_cfd_risk:Optional[float]=None
    annotation:Dict[str,Any]=field(default_factory=dict)
    @property
    def full_target(self): return self.sequence+self.pam
    @property
    def spacer_start(self): return self.start if self.strand=='+' else self.start+3
    @property
    def spacer_end(self): return self.end-3 if self.strand=='+' else self.end
    @property
    def cut_boundary(self): return self.start+17 if self.strand=='+' else self.start+6
    def to_dict(self):
        return {'Spacer (20 nt)':self.sequence,'PAM':self.pam,'Strand':self.strand,
        'Start':self.start+1,'End':self.end,'GC%':round(self.gc_content,1),'Heuristic':round(self.score,1),
        'Doench RS2':None if self.doench_score is None else round(self.doench_score,1),
        'MIT specificity':None if self.specificity_score is None else round(self.specificity_score,1),
        'CFD specificity':None if self.cfd_specificity is None else round(self.cfd_specificity,1),
        'Genome specificity':self.genome_specificity,'Off-target hits':self.off_target_count,
        'Application':self.application,'Notes':'; '.join(self.notes) if self.notes else 'No flags',
        'Spacer start':self.spacer_start+1,'Spacer end':self.spacer_end,
        'Cut after base':self.cut_boundary,'Screen status':self.screen_status,
        'Exact reference matches':self.exact_match_count,
        'Screened mismatch radius':self.screened_mismatch_radius,
        'Critical hits':None if self.screened_mismatch_radius is None else self.risk_counts.get('Critical', 0),
        'High-risk hits':None if self.screened_mismatch_radius is None else self.risk_counts.get('High', 0),
        'Moderate-risk hits':None if self.screened_mismatch_radius is None else self.risk_counts.get('Moderate', 0),
        'Low-risk hits':None if self.screened_mismatch_radius is None else self.risk_counts.get('Low', 0),
        'Maximum per-site MIT risk':None if self.max_mit_risk is None else round(self.max_mit_risk*100, 2),
        'Maximum per-site CFD risk':None if self.max_cfd_risk is None else round(self.max_cfd_risk*100, 2),
        'Intended locus':None if self.intended_target is None else {'contig':self.intended_target[0], 'start_0based':self.intended_target[1], 'strand':self.intended_target[2]},
        'Ambiguous reference bases':self.reference_ambiguous_bases, **self.annotation}


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: fetch_gene_ensembl
# ----------------------------------------------------------------------------
def fetch_gene_ensembl(gene_name:str,species='homo_sapiens')->Tuple[str,str,str]:
    species=species.lower().strip().replace(' ','_'); server='https://rest.ensembl.org'; headers={'Accept':'application/json','User-Agent':NCBI_TOOL}
    r=get(f'{server}/lookup/symbol/{species}/{gene_name}?expand=1',headers=headers,timeout=30)
    if not r.ok: raise ValueError(f'Ensembl lookup failed ({r.status_code}).')
    data=r.json(); txs=data.get('Transcript',[])
    if not txs: raise ValueError('No transcripts returned by Ensembl.')
    cid=str(data.get('canonical_transcript','')).split('.')[0]
    tx=next((t for t in txs if t.get('is_canonical') or t.get('id')==cid),txs[0]); tid=tx['id']
    sr=get(f'{server}/sequence/id/{tid}?type=cdna',headers={'Accept':'text/plain','User-Agent':NCBI_TOOL},timeout=30)
    if not sr.ok: raise ValueError(f'Could not fetch cDNA for {tid}.')
    return tid,f"Ensembl {tid} | {data.get('display_name',gene_name)} | {data.get('description','')}".strip(' |'),clean_dna_sequence(sr.text)


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: fetch_gene_ncbi
# ----------------------------------------------------------------------------
def fetch_gene_ncbi(gene_name:str,organism:str,max_retries=3)->Tuple[str,str,str]:
    q=f'{gene_name}[Gene Name] AND {organism}[Organism] AND alive[prop]'; common={'tool':NCBI_TOOL,'email':NCBI_EMAIL}
    for attempt in range(max_retries):
        try:
            s=get(f'{NCBI_EUTILS}/esearch.fcgi',params={**common,'db':'gene','term':q,'retmax':5,'retmode':'json'},timeout=30); s.raise_for_status()
            ids=s.json().get('esearchresult',{}).get('idlist',[])
            if not ids: raise ValueError(f"No gene found for '{gene_name}' in '{organism}'.")
            l=get(f'{NCBI_EUTILS}/elink.fcgi',params={**common,'dbfrom':'gene','db':'nuccore','id':ids[0],'linkname':'gene_nuccore_refseqrna','retmode':'json'},timeout=30); l.raise_for_status()
            dbs=l.json().get('linksets',[{}])[0].get('linksetdbs',[]); nids=dbs[0].get('links',[]) if dbs else []
            if not nids: raise ValueError('Could not resolve a RefSeq RNA record.')
            f=get(f'{NCBI_EUTILS}/efetch.fcgi',params={**common,'db':'nuccore','id':','.join(nids[:20]),'rettype':'gb','retmode':'text'},timeout=30); f.raise_for_status()
            recs=list(SeqIO.parse(StringIO(f.text),'genbank')); selected=min(recs,key=lambda r:(0 if r.id.upper().startswith('NM_') else 1,r.id))
            return selected.id,selected.description,clean_dna_sequence(str(selected.seq))
        except Exception:
            if attempt==max_retries-1: raise
            time.sleep(1+attempt)
    raise RuntimeError('NCBI lookup failed.')


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: fetch_sequence
# ----------------------------------------------------------------------------
def fetch_sequence(gene_name,organism,source='ncbi'):
    if source.lower()=='ncbi': return fetch_gene_ncbi(gene_name,organism)
    m={'homo sapiens':'homo_sapiens','human':'homo_sapiens','mus musculus':'mus_musculus','mouse':'mus_musculus','rattus norvegicus':'rattus_norvegicus','rat':'rattus_norvegicus','danio rerio':'danio_rerio','zebrafish':'danio_rerio'}
    return fetch_gene_ensembl(gene_name,m.get(organism.lower(),organism.lower().replace(' ','_')))


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _gc_content
# ----------------------------------------------------------------------------
def _gc_content(s): return 100*(s.upper().count('G')+s.upper().count('C'))/len(s) if s else 0.0

# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _has_homopolymer
# ----------------------------------------------------------------------------
def _has_homopolymer(s,max_run=4): return bool(re.search(rf'(A{{{max_run},}}|T{{{max_run},}}|G{{{max_run},}}|C{{{max_run},}})',s.upper()))


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: score_breakdown
# ----------------------------------------------------------------------------
def score_breakdown(spacer,pam,application=None,start=None,sequence_length=None):
    if len(spacer)!=20 or len(pam)!=3:return {'Final score':0.0}
    gc=_gc_content(spacer); c={'Baseline':50.0,'GC content':15.0 if 40<=gc<=70 else 5.0 if 30<=gc<=80 else -15.0,
    'Position 20':10.0 if spacer[19]=='G' else -8.0 if spacer[19]=='C' else -12.0 if spacer[19]=='T' else 0.0,
    'Position 19':5.0 if spacer[18] in 'AG' else 0.0,'Position 16':3.0 if spacer[15]=='C' else 0.0,'Position 18':3.0 if spacer[17]=='C' else 0.0,
    'PAM context':5.0 if pam.upper()=='CGG' else -5.0 if pam.upper().startswith('T') else 0.0,
    'Homopolymers':-20.0 if ('GGGG' in spacer or 'TTTT' in spacer) else -8.0 if ('GGG' in spacer or 'TTT' in spacer) else 0.0}
    adj=0.0
    if application and start is not None and sequence_length:
        if application=='knockdown': adj=15.0 if start<400 else -5.0
        elif application=='knockout' and start<sequence_length*.3: adj=8.0
    c['Application position']=adj; c['Final score']=round(max(0,min(100,sum(c.values()))),1); return c


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _doench_like_score
# ----------------------------------------------------------------------------
def _doench_like_score(spacer,pam): return score_breakdown(spacer,pam)['Final score']


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: mit_offtarget_score
# ----------------------------------------------------------------------------
def mit_offtarget_score(spacer,offtarget,pam='NGG'):
    spacer=clean_dna_sequence(spacer); offtarget=clean_dna_sequence(offtarget)
    if len(spacer)!=20 or len(offtarget)!=20 or (set(spacer+offtarget)-set('ACGT')): raise ValueError('MIT scoring requires two 20 nt sequences.')
    pos=[i+1 for i,(a,b) in enumerate(zip(spacer,offtarget)) if a!=b]
    if not pos:return 1.0
    m=len(pos); d=19 if m==1 else (max(pos)-min(pos))/(m-1)
    t1=1.0
    for p in pos:t1*=1-MIT_WEIGHTS[p-1]
    return max(0.0,min(1.0,t1*(1/(m*m))*(1/(((19-d)/19)*4+1))))


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: mit_specificity
# ----------------------------------------------------------------------------
def mit_specificity(scores): return round(100/(1+sum(scores)),2)


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _optional_cfd
# ----------------------------------------------------------------------------
def _optional_cfd(spacer, offtarget):
    return cfd_score(spacer, offtarget)


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _risk_label
# ----------------------------------------------------------------------------
def _risk_label(mm,seed):
    if mm==0:return 'Critical'
    if mm==1 or (mm==2 and seed<=1):return 'High'
    if mm<=3 and seed<=1:return 'Moderate'
    return 'Low'


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: iter_pam_sites
# ----------------------------------------------------------------------------
def iter_pam_sites(sequence, contig="reference"):
    """Yield real, unambiguous 23-base targets; start is left edge on input strand."""
    for pattern, strand in ((PAM_PATTERN, "+"), (PAM_PATTERN_REV, "-")):
        for match in pattern.finditer(sequence):
            ps = match.start()
            start, end = (ps - 20, ps + 3) if strand == "+" else (ps, ps + 23)
            if start < 0 or end > len(sequence):
                continue
            target = sequence[start:end]
            if set(target) - set("ACGT"):
                continue
            if strand == "-":
                target = str(Seq(target).reverse_complement())
            yield (contig, start, strand, target[:20], target[20:])



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: ReferenceIndex
# ----------------------------------------------------------------------------
class ReferenceIndex:
    """Reusable, bounded local NGG site index with vectorized mismatch search."""
    def __init__(self, reference):
        records = parse_reference(reference) if isinstance(reference, str) else reference
        self.contigs = {str(k): clean_dna_sequence(v) for k, v in records.items()}
        if not self.contigs or any(not k or not v for k, v in self.contigs.items()):
            raise ValueError("Reference must contain nonempty, named sequences.")
        self.total_bases = sum(map(len, self.contigs.values()))
        if self.total_bases > 5_000_000:
            raise ValueError("Local reference exceeds 5,000,000 bp; use an indexed genome tool.")
        self.ambiguous_bases = sum(v.count("N") for v in self.contigs.values())
        self.sites = []
        for k, seq in self.contigs.items():
            for site in iter_pam_sites(seq, k):
                if len(self.sites) >= 500_000:
                    raise ValueError("Reference exceeds 500,000 NGG sites; use an indexed genome tool.")
                self.sites.append(site)
        self.encoded = np.frombuffer("".join(x[3] for x in self.sites).encode("ascii"), dtype=np.uint8).reshape(-1, 20)



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: analyze_offtargets
# ----------------------------------------------------------------------------
def analyze_offtargets(spacer, whole_genome, max_mismatches=3, max_hits=250, intended_target=None):
    """Score ALL local hits; cap only displayed details, never counts or risk sums.

    No intended target is guessed. To exclude one verified exact site, provide
    (contig, zero-based target-left-edge, strand). All exact copies otherwise remain.
    Only NGG PAMs, substitution mismatches, and the supplied reference are searched.
    """
    spacer = clean_dna_sequence(spacer)
    if len(spacer) != 20 or set(spacer) - set("ACGT"):
        raise ValueError("Off-target screening requires an unambiguous 20 nt spacer.")
    if type(max_mismatches) is not int or not 0 <= max_mismatches <= 4:
        raise ValueError("Mismatch limit must be an integer from 0 to 4.")
    if type(max_hits) is not int or max_hits < 0:
        raise ValueError("Display hit limit must be a nonnegative integer.")
    index = whole_genome if isinstance(whole_genome, ReferenceIndex) else ReferenceIndex(whole_genome)
    if intended_target is not None:
        try:
            intended_target = TargetLocus(*intended_target)
        except (TypeError, ValueError):
            raise ValueError('Intended target must contain contig, start and strand.') from None
        if not isinstance(intended_target.contig, str) or type(intended_target.start) is not int or intended_target.start < 0 or intended_target.strand not in {'+', '-'}:
            raise ValueError('Invalid intended target contig, start or strand.')
    query = np.frombuffer(spacer.encode("ascii"), dtype=np.uint8)
    hits, total, exact, mit_sum, cfd_sum, cfd_complete = [], 0, 0, 0.0, 0.0, True
    risk_counts = {'Critical':0, 'High':0, 'Moderate':0, 'Low':0}
    max_mit_risk, max_cfd_risk = 0.0, 0.0
    excluded = False
    for offset in range(0, len(index.sites), 50_000):
        chunk = index.encoded[offset:offset + 50_000]
        counts = np.count_nonzero(chunk != query, axis=1)
        for local in np.flatnonzero(counts <= max_mismatches):
            cname, start, strand, cand, pam = index.sites[offset + int(local)]
            mm = int(counts[local])
            exact += mm == 0
            if intended_target == (cname, start, strand) and mm == 0:
                excluded = True
                continue
            pos = tuple(i + 1 for i, (a, b) in enumerate(zip(spacer, cand)) if a != b)
            seed = sum(p >= 13 for p in pos)
            mit = mit_offtarget_score(spacer, cand, pam)
            cfd = _optional_cfd(spacer, cand)
            risk = _risk_label(mm, seed)
            total += 1
            mit_sum += mit
            cfd_complete = cfd_complete and cfd is not None
            cfd_sum += cfd or 0.0
            risk_counts[risk] += 1
            max_mit_risk = max(max_mit_risk, mit)
            max_cfd_risk = max(max_cfd_risk, cfd or 0.0)
            if max_hits:
                hits.append(OffTargetHit(start, strand, cand, pam, mm, pos, seed, risk, cname, mit, cfd))
                if len(hits) > max_hits * 2:
                    hits.sort(key=lambda h: (h.mismatches, h.seed_mismatches, -h.mit_score, h.contig, h.start))
                    hits = hits[:max_hits]
    if intended_target is not None and not excluded:
        raise ValueError("Declared intended target is not an exact NGG site in this reference. Check contig, position and strand.")
    hits.sort(key=lambda h: (h.mismatches, h.seed_mismatches, -h.mit_score, h.contig, h.start))
    status = "verified_locus" if excluded else ("exact_matches_retained" if exact else "no_exact_match")
    return OffTargetReport(
        hits=tuple(hits[:max_hits]), specificity_score=round(100 / (1 + mit_sum), 2),
        pam_sites_scanned=len(index.sites), on_target_excluded=excluded,
        cfd_specificity=round(100 / (1 + cfd_sum), 2) if cfd_complete else None,
        reference_contigs=len(index.contigs), total_hits=total, exact_matches=exact,
        truncated=total > max_hits, intended_target_status=status,
        ambiguous_bases=index.ambiguous_bases, screened_mismatch_radius=max_mismatches,
        risk_counts=risk_counts, max_mit_risk=max_mit_risk, max_cfd_risk=max_cfd_risk,
    )



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: screen_guides
# ----------------------------------------------------------------------------
def screen_guides(guides, reference, max_mismatches=3, intended_targets=None, intended_region=None, target_sequence=None):
    """Build once and reuse across guides; refuse excessive work explicitly."""
    index = reference if isinstance(reference, ReferenceIndex) else ReferenceIndex(reference)
    if intended_region is not None:
        intended_targets = intended_loci(guides, target_sequence, index, intended_region)
    unique = len({(g.sequence, (intended_targets or {}).get((g.start, g.strand))) for g in guides})
    if unique * len(index.sites) > 100_000_000:
        raise ValueError("Local screen exceeds the 100 million comparison budget. Use a smaller target/reference or an indexed genome tool.")
    reports = {}
    for guide in guides:
        intended = (intended_targets or {}).get((guide.start, guide.strand))
        key = (guide.sequence, intended)
        if key not in reports:
            reports[key] = analyze_offtargets(guide.sequence, index, max_mismatches, intended_target=intended)
        report = reports[key]
        guide.intended_target = intended
        guide.reference_ambiguous_bases = report.ambiguous_bases
        guide.screened_mismatch_radius = report.screened_mismatch_radius
        guide.risk_counts = dict(report.risk_counts)
        guide.max_mit_risk = report.max_mit_risk
        guide.max_cfd_risk = report.max_cfd_risk
        guide.specificity_score = report.mit_specificity
        guide.cfd_specificity = report.cfd_specificity
        guide.off_target_count = report.total_hits
        guide.exact_match_count = report.exact_matches
        guide.screen_status = report.intended_target_status
        guide.off_target_details = [h.to_dict() for h in report.hits]
        if report.ambiguous_bases:
            guide.notes.append("Ambiguous reference bases excluded from search")
    return guides



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: calculate_offtarget_score
# ----------------------------------------------------------------------------
def calculate_offtarget_score(spacer,pam,whole_genome=None,max_mismatches=3):
    if not whole_genome:return 0,0.0,['Genome not provided']
    r=analyze_offtargets(spacer,whole_genome,max_mismatches); return r.total_hits,r.specificity_score,[f'{h.contig}:{h.start+1} | {h.strand} | {h.sequence}{h.pam} | {h.mismatches} mismatch(es) | MIT {h.mit_score*100:.1f}' for h in r.hits]


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _context30
# ----------------------------------------------------------------------------
def _context30(sequence,start,end,strand):
    if strand=='+':
        if start<4 or end+3>len(sequence):return None
        return sequence[start-4:end+3]
    if start<3 or end+4>len(sequence):return None
    return str(Seq(sequence[start-3:end+4]).reverse_complement())


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: doench_rs2_score
# ----------------------------------------------------------------------------
def doench_rs2_score(context30):
    return rs2_score(context30)



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: design_guides
# ----------------------------------------------------------------------------
def design_guides(sequence, application='knockout', pam='NGG', spacer_len=20,
                  min_score=30, max_guides=50, prefer_5prime=False,
                  genome_context=None, max_mismatches=3, intended_region=None):
    if pam.upper() != 'NGG' or spacer_len != 20:
        raise ValueError('Only SpCas9 20 nt + NGG design is supported.')
    if application not in {'knockout', 'crispri'}:
        raise ValueError('Use knockout or the explicit TSS-aware CRISPRi workflow.')
    if not math.isfinite(float(min_score)) or not 0 <= min_score <= 100:
        raise ValueError('Heuristic threshold must be between 0 and 100.')
    if not isinstance(max_guides, int) or max_guides < 1:
        raise ValueError('Maximum guides must be a positive integer.')
    sequence = clean_dna_sequence(sequence)
    if len(sequence) > 2_000_000:
        raise ValueError('Target exceeds the 2,000,000 bp design limit.')
    guides = []
    for _, start, strand, spacer, target_pam in iter_pam_sites(sequence):
        score = score_breakdown(spacer, target_pam, application if prefer_5prime else None, start, len(sequence))['Final score']
        if score < min_score:
            continue
        guide = GuideRNA(spacer, target_pam, strand, start, start + 23, _gc_content(spacer), score, application=application)
        # A nuclease activity model does not predict dCas9-KRAB repression.
        if application == 'knockout':
            guide.doench_score = doench_rs2_score(_context30(sequence, start, start + 23, strand))
        guide.scoring_method = 'Sequence heuristic; RS2 reported separately when available'
        if not 40 <= guide.gc_content <= 70:
            guide.notes.append('GC outside heuristic 40-70% range')
        if _has_homopolymer(spacer, 4):
            guide.notes.append('Homopolymer >=4')
        if 'TTTT' in spacer:
            guide.notes.append('Poly-T may interfere with U6 expression')
        guides.append(guide)
        if len(guides) > 50_000:
            raise ValueError("Target produces more than 50,000 candidates; select a smaller region.")
    if genome_context is not None:
        screen_guides(guides, genome_context, max_mismatches, intended_region=intended_region, target_sequence=sequence)
    # Never compare different model scales or put missing-model edge guides last.
    guides.sort(key=lambda g: (-g.score, -(g.specificity_score if g.specificity_score is not None else -1), g.start, g.strand))
    return guides[:max_guides]



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: guidescan2_available
# ----------------------------------------------------------------------------
def guidescan2_available(): return shutil.which('guidescan') is not None


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: run_guidescan2
# ----------------------------------------------------------------------------
def run_guidescan2(guides,index_path,max_mismatches=4,alt_pam='NAG'):
    """Run GuideScan2 succinct CSV mode. Requires a prebuilt local genome index."""
    if not guidescan2_available(): raise RuntimeError('GuideScan2 executable not found. Install with bioconda: guidescan.')
    if not index_path: raise ValueError('GuideScan2 index path is required.')
    with tempfile.TemporaryDirectory() as td:
        inp=Path(td)/'kmers.csv'; out=Path(td)/'result.csv'
        with inp.open('w',newline='') as fh:
            w=csv.writer(fh); w.writerow(['id','sequence','pam','chromosome','position','sense'])
            for i,g in enumerate(guides,1): w.writerow([f'g{i}',g.sequence,'NGG','__unmapped__',1,g.strand])
        cmd=['guidescan','enumerate',str(index_path),'-f',str(inp),'-o',str(out),'--format','csv','--mode','succinct','-m',str(max_mismatches)]
        if alt_pam:cmd += ['-a',alt_pam]
        p=subprocess.run(cmd,capture_output=True,text=True,timeout=600)
        if p.returncode!=0: raise RuntimeError('GuideScan2 failed: '+p.stderr[-600:])
        if not out.exists():
            raise RuntimeError('GuideScan2 returned no output file; no successful screen is recorded.')
        with out.open() as handle:
            reader = csv.DictReader(handle)
            if not {'id', 'sequence', 'specificity'}.issubset(reader.fieldnames or []):
                raise RuntimeError('GuideScan2 output has an unsupported CSV schema.')
            rows = list(reader)
        if {row['id'] for row in rows} != {f'g{i}' for i in range(1, len(guides) + 1)}:
            raise RuntimeError('GuideScan2 output is incomplete or contains unexpected guide IDs.')
        return rows


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: validate_guide
# ----------------------------------------------------------------------------
def validate_guide(guide, genome_context=None, min_score=40):
    spec = guide.specificity_score
    screened = spec is not None
    if genome_context is not None and not screened:
        screen_guides([guide], genome_context)
        spec, screened = guide.specificity_score, True
    full_local_scope = screened and guide.screened_mismatch_radius is not None and guide.screened_mismatch_radius >= 3
    no_critical_or_high_hits = screened and guide.risk_counts.get('Critical', 0) == 0 and guide.risk_counts.get('High', 0) == 0
    aggregate_specificity_ok = screened and spec >= 50
    locus_evidence_ok = screened and guide.screen_status == 'verified_locus' and guide.exact_match_count == 1 and guide.reference_ambiguous_bases == 0
    checks = {
        'length_ok': len(guide.sequence) == 20,
        'unambiguous': not (set(guide.sequence) - set('ACGT')),
        'pam_ok': bool(re.fullmatch('[ACGT]GG', guide.pam)),
        'gc_in_preferred_range': 40 <= guide.gc_content <= 70,
        'gc_in_acceptable_range': 30 <= guide.gc_content <= 80,
        'no_extreme_homopolymer': not _has_homopolymer(guide.sequence, 5),
        'score_above_threshold': guide.score >= min_score,
        'no_poly_t': 'TTTT' not in guide.sequence,
        'specificity_screened': screened,
        'screen_scope_complete': full_local_scope,
        'no_critical_or_high_hits': no_critical_or_high_hits,
        'aggregate_specificity_ok': aggregate_specificity_ok,
        'locus_evidence_ok': locus_evidence_ok,
        'specificity_ok': full_local_scope and no_critical_or_high_hits and aggregate_specificity_ok and locus_evidence_ok,
    }
    checks['sequence_checks_pass'] = all(checks[k] for k in (
        'length_ok', 'unambiguous', 'pam_ok', 'gc_in_acceptable_range',
        'no_extreme_homopolymer', 'score_above_threshold', 'no_poly_t'))
    checks['overall_pass'] = checks['sequence_checks_pass'] and checks['specificity_ok']
    return checks



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: design_from_gene
# ----------------------------------------------------------------------------
def design_from_gene(gene_name,organism,application='knockout',source='ncbi',max_guides=20,min_score=30,genome_context=None,max_mismatches=3):
    if application != 'knockout':
        raise ValueError('Use the explicit TSS-aware CRISPRi workflow.')
    if source.lower() != 'ncbi':
        raise ValueError('Annotated knockout gene design requires NCBI genomic CDS lookup. Ensembl transcript discovery is disabled.')
    from knockout import fetch_ncbi_knockout_context, design_knockout_guides
    context = fetch_ncbi_knockout_context(gene_name, organism)
    guides = design_knockout_guides(context, max_guides, min_score, genome_context, max_mismatches)
    return context.accession, context.description, context.sequence, guides
