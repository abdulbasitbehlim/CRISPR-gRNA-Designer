#!/usr/bin/env python3
"""CRISPR Studio v3 core: SpCas9 design, MIT specificity, local multi-FASTA screening,
optional Doench Rule Set 2/CFD providers, and GuideScan2 whole-genome integration.
"""
from __future__ import annotations
import csv, os, re, shutil, subprocess, tempfile, time
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import requests
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
    @property
    def full_target(self): return self.sequence+self.pam
    def to_dict(self):
        return {'Spacer (20 nt)':self.sequence,'PAM':self.pam,'Strand':self.strand,
        'Start':self.start+1,'End':self.end,'GC%':round(self.gc_content,1),'Heuristic':round(self.score,1),
        'Doench RS2':None if self.doench_score is None else round(self.doench_score,1),
        'MIT specificity':None if self.specificity_score is None else round(self.specificity_score,1),
        'CFD specificity':None if self.cfd_specificity is None else round(self.cfd_specificity,1),
        'Genome specificity':self.genome_specificity,'Off-target hits':self.off_target_count,
        'Application':self.application,'Notes':'; '.join(self.notes) if self.notes else 'No flags'}

def clean_dna_sequence(raw:str)->str:
    if not raw:return ''
    lines=[x.strip() for x in raw.splitlines() if x.strip() and not x.lstrip().startswith('>')]
    seq=''.join(lines).upper().replace('U','T'); seq=re.sub(r'[\s\d]','',seq); seq=IUPAC_AMBIGUOUS.sub('N',seq)
    bad=sorted(set(seq)-DNA_ALPHABET)
    if bad: raise ValueError('Sequence contains unsupported character(s): '+', '.join(bad))
    return seq

def parse_reference(raw:str)->Dict[str,str]:
    """Preserve FASTA contigs; plain sequence becomes one contig."""
    if not raw:return {}
    if not any(line.lstrip().startswith('>') for line in raw.splitlines()): return {'reference':clean_dna_sequence(raw)}
    records={}
    for rec in SeqIO.parse(StringIO(raw),'fasta'):
        name=rec.id or f'contig_{len(records)+1}'; records[name]=clean_dna_sequence(str(rec.seq))
    if not records: raise ValueError('No readable FASTA records found.')
    return records

def fetch_gene_ensembl(gene_name:str,species='homo_sapiens')->Tuple[str,str,str]:
    species=species.lower().strip().replace(' ','_'); server='https://rest.ensembl.org'; headers={'Accept':'application/json','User-Agent':NCBI_TOOL}
    r=requests.get(f'{server}/lookup/symbol/{species}/{gene_name}?expand=1',headers=headers,timeout=30)
    if not r.ok: raise ValueError(f'Ensembl lookup failed ({r.status_code}).')
    data=r.json(); txs=data.get('Transcript',[])
    if not txs: raise ValueError('No transcripts returned by Ensembl.')
    cid=str(data.get('canonical_transcript','')).split('.')[0]
    tx=next((t for t in txs if t.get('is_canonical') or t.get('id')==cid),txs[0]); tid=tx['id']
    sr=requests.get(f'{server}/sequence/id/{tid}?type=cdna',headers={'Accept':'text/plain','User-Agent':NCBI_TOOL},timeout=30)
    if not sr.ok: raise ValueError(f'Could not fetch cDNA for {tid}.')
    return tid,f"Ensembl {tid} | {data.get('display_name',gene_name)} | {data.get('description','')}".strip(' |'),clean_dna_sequence(sr.text)

def fetch_gene_ncbi(gene_name:str,organism:str,max_retries=3)->Tuple[str,str,str]:
    q=f'{gene_name}[Gene Name] AND {organism}[Organism] AND alive[prop]'; common={'tool':NCBI_TOOL,'email':NCBI_EMAIL}
    for attempt in range(max_retries):
        try:
            s=requests.get(f'{NCBI_EUTILS}/esearch.fcgi',params={**common,'db':'gene','term':q,'retmax':5,'retmode':'json'},timeout=30); s.raise_for_status()
            ids=s.json().get('esearchresult',{}).get('idlist',[])
            if not ids: raise ValueError(f"No gene found for '{gene_name}' in '{organism}'.")
            l=requests.get(f'{NCBI_EUTILS}/elink.fcgi',params={**common,'dbfrom':'gene','db':'nuccore','id':ids[0],'linkname':'gene_nuccore_refseqrna','retmode':'json'},timeout=30); l.raise_for_status()
            dbs=l.json().get('linksets',[{}])[0].get('linksetdbs',[]); nids=dbs[0].get('links',[]) if dbs else []
            if not nids: raise ValueError('Could not resolve a RefSeq RNA record.')
            f=requests.get(f'{NCBI_EUTILS}/efetch.fcgi',params={**common,'db':'nuccore','id':','.join(nids[:20]),'rettype':'gb','retmode':'text'},timeout=30); f.raise_for_status()
            recs=list(SeqIO.parse(StringIO(f.text),'genbank')); selected=min(recs,key=lambda r:(0 if r.id.upper().startswith('NM_') else 1,r.id))
            return selected.id,selected.description,clean_dna_sequence(str(selected.seq))
        except Exception:
            if attempt==max_retries-1: raise
            time.sleep(1+attempt)
    raise RuntimeError('NCBI lookup failed.')

def fetch_sequence(gene_name,organism,source='ncbi'):
    if source.lower()=='ncbi': return fetch_gene_ncbi(gene_name,organism)
    m={'homo sapiens':'homo_sapiens','human':'homo_sapiens','mus musculus':'mus_musculus','mouse':'mus_musculus','rattus norvegicus':'rattus_norvegicus','rat':'rattus_norvegicus','danio rerio':'danio_rerio','zebrafish':'danio_rerio'}
    return fetch_gene_ensembl(gene_name,m.get(organism.lower(),organism.lower().replace(' ','_')))

def _gc_content(s): return 100*(s.upper().count('G')+s.upper().count('C'))/len(s) if s else 0.0
def _has_homopolymer(s,max_run=4): return bool(re.search(rf'(A{{{max_run},}}|T{{{max_run},}}|G{{{max_run},}}|C{{{max_run},}})',s.upper()))

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

def _doench_like_score(spacer,pam): return score_breakdown(spacer,pam)['Final score']

def mit_offtarget_score(spacer,offtarget,pam='NGG'):
    spacer=clean_dna_sequence(spacer); offtarget=clean_dna_sequence(offtarget)
    if len(spacer)!=20 or len(offtarget)!=20: raise ValueError('MIT scoring requires two 20 nt sequences.')
    pos=[i+1 for i,(a,b) in enumerate(zip(spacer,offtarget)) if a!=b]
    if not pos:return 1.0
    m=len(pos); d=19 if m==1 else (max(pos)-min(pos))/(m-1)
    t1=1.0
    for p in pos:t1*=1-MIT_WEIGHTS[p-1]
    return max(0.0,min(1.0,t1*(1/(m*m))*(1/(((19-d)/19)*4+1))))

def mit_specificity(scores): return round(100/(1+sum(scores)),2)

def _optional_cfd(spacer,offtarget):
    try:
        from guidemaker.cfd_score_calculator import calc_cfd
        return float(calc_cfd(spacer,offtarget))
    except Exception:return None

def _risk_label(mm,seed):
    if mm==0:return 'Critical'
    if mm==1 or (mm==2 and seed<=1):return 'High'
    if mm<=3 and seed<=1:return 'Moderate'
    return 'Low'

def analyze_offtargets(spacer,whole_genome,max_mismatches=3,max_hits=250):
    spacer=clean_dna_sequence(spacer)
    if len(spacer)!=20: raise ValueError('Off-target screening requires a 20 nt spacer.')
    contigs=parse_reference(whole_genome) if isinstance(whole_genome,str) else whole_genome
    candidates=[]
    for cname,ref in contigs.items():
        for m in PAM_PATTERN.finditer(ref):
            ps=m.start(); ss=ps-20
            if ss>=0:
                cand=ref[ss:ps]; pam=ref[ps:ps+3]
                if 'N' not in cand:candidates.append((cname,ss,'+',cand,pam))
        for m in PAM_PATTERN_REV.finditer(ref):
            ps=m.start(); end=ps+23
            if end<=len(ref):
                cand=str(Seq(ref[ps+3:end]).reverse_complement()); pam=str(Seq(ref[ps:ps+3]).reverse_complement())
                if 'N' not in cand:candidates.append((cname,ps,'-',cand,pam))
    excluded=False; hits=[]
    for cname,start,strand,cand,pam in candidates:
        pos=tuple(i+1 for i,(a,b) in enumerate(zip(spacer,cand)) if a!=b); mm=len(pos)
        if mm>max_mismatches:continue
        if mm==0 and not excluded: excluded=True; continue
        seed=sum(p>=13 for p in pos); mit=mit_offtarget_score(spacer,cand,pam); cfd=_optional_cfd(spacer,cand)
        hits.append(OffTargetHit(start,strand,cand,pam,mm,pos,seed,_risk_label(mm,seed),cname,mit,cfd))
    hits.sort(key=lambda h:(h.mismatches,h.seed_mismatches,-h.mit_score,h.contig,h.start)); hits=hits[:max_hits]
    mit_spec=mit_specificity([h.mit_score for h in hits]); cfd_vals=[h.cfd_score for h in hits if h.cfd_score is not None]
    cfd_spec=round(100/(1+sum(cfd_vals)),2) if cfd_vals else None
    return OffTargetReport(tuple(hits),mit_spec,len(candidates),excluded,cfd_spec,len(contigs))

def calculate_offtarget_score(spacer,pam,whole_genome=None,max_mismatches=3):
    if not whole_genome:return 0,0.0,['Genome not provided']
    r=analyze_offtargets(spacer,whole_genome,max_mismatches); return len(r.hits),r.specificity_score,[f'{h.contig}:{h.start+1} | {h.strand} | {h.sequence}{h.pam} | {h.mismatches} mismatch(es) | MIT {h.mit_score*100:.1f}' for h in r.hits]

def _context30(sequence,start,end,strand):
    if strand=='+':
        if start<4 or end+3>len(sequence):return None
        return sequence[start-4:end+3]
    if start<3 or end+4>len(sequence):return None
    return str(Seq(sequence[start-3:end+4]).reverse_complement())

def doench_rs2_score(context30):
    if not context30 or len(context30)!=30:return None
    try:
        import numpy as np
        from guidemaker.doench_predict import predict
        val=float(predict(np.array([context30]),num_threads=1)[0]); return round(max(0,min(1,val))*100,2)
    except Exception:return None

def design_guides(sequence,application='knockout',pam='NGG',spacer_len=20,min_score=30,max_guides=50,prefer_5prime=True,genome_context=None,max_mismatches=3):
    if pam.upper()!='NGG' or spacer_len!=20:raise ValueError('v3 supports SpCas9 20 nt + NGG.')
    application=application.lower().strip(); sequence=clean_dna_sequence(sequence); guides=[]
    def add(sp,p,strand,start,end):
        sc=score_breakdown(sp,p,application if prefer_5prime else None,start,len(sequence))['Final score']
        if sc<min_score:return
        g=GuideRNA(sp,p,strand,start,end,_gc_content(sp),sc,application=application); g.doench_score=doench_rs2_score(_context30(sequence,start,end,strand))
        g.scoring_method='Doench RS2 + heuristic' if g.doench_score is not None else 'Heuristic (Doench provider unavailable)'
        if genome_context:
            r=analyze_offtargets(sp,genome_context,max_mismatches); g.specificity_score=r.mit_specificity; g.cfd_specificity=r.cfd_specificity; g.off_target_count=len(r.hits); g.off_target_details=[h.to_dict() for h in r.hits]
        if not 40<=g.gc_content<=70:g.notes.append('GC outside preferred 40-70%')
        if _has_homopolymer(sp,4):g.notes.append('Homopolymer >=4')
        guides.append(g)
    for m in PAM_PATTERN.finditer(sequence):
        ps=m.start(); ss=ps-20
        if ss>=0:add(sequence[ss:ps],sequence[ps:ps+3],'+',ss,ps+3)
    for m in PAM_PATTERN_REV.finditer(sequence):
        ps=m.start(); end=ps+23
        if end<=len(sequence):add(str(Seq(sequence[ps+3:end]).reverse_complement()),str(Seq(sequence[ps:ps+3]).reverse_complement()),'-',ps,end)
    guides.sort(key=lambda g:(g.doench_score if g.doench_score is not None else g.score,g.specificity_score or -1,-g.start),reverse=True); return guides[:max_guides]

def guidescan2_available(): return shutil.which('guidescan') is not None

def run_guidescan2(guides,index_path,max_mismatches=4,alt_pam='NAG'):
    """Run GuideScan2 succinct CSV mode. Requires a prebuilt local genome index."""
    if not guidescan2_available(): raise RuntimeError('GuideScan2 executable not found. Install with bioconda: guidescan.')
    if not index_path: raise ValueError('GuideScan2 index path is required.')
    with tempfile.TemporaryDirectory() as td:
        inp=Path(td)/'kmers.csv'; out=Path(td)/'result.csv'
        with inp.open('w',newline='') as fh:
            w=csv.writer(fh); w.writerow(['id','sequence','pam','chromosome','position','sense'])
            for i,g in enumerate(guides,1): w.writerow([f'g{i}',g.sequence,'NGG','','',g.strand])
        cmd=['guidescan','enumerate',str(index_path),'-f',str(inp),'-o',str(out),'--format','csv','--mode','succinct','-m',str(max_mismatches)]
        if alt_pam:cmd += ['-a',alt_pam]
        p=subprocess.run(cmd,capture_output=True,text=True,timeout=600)
        if p.returncode!=0: raise RuntimeError('GuideScan2 failed: '+p.stderr[-600:])
        rows=list(csv.DictReader(out.open())) if out.exists() else []
        return rows

def validate_guide(guide,genome_context=None):
    spec=guide.specificity_score
    if genome_context and spec is None:spec=analyze_offtargets(guide.sequence,genome_context).specificity_score
    c={'length_ok':len(guide.sequence)==20,'pam_ok':len(guide.pam)==3 and guide.pam.upper().endswith('GG'),'gc_in_preferred_range':40<=guide.gc_content<=70,'gc_in_acceptable_range':30<=guide.gc_content<=80,'no_extreme_homopolymer':not _has_homopolymer(guide.sequence,5),'score_above_threshold':guide.score>=40,'no_poly_t':'TTTT' not in guide.sequence,'specificity_screened':spec is not None,'specificity_ok':spec is None or spec>=50}
    c['overall_pass']=all(c[k] for k in ['length_ok','pam_ok','gc_in_acceptable_range','no_extreme_homopolymer','score_above_threshold','no_poly_t','specificity_ok']); return c

def design_from_gene(gene_name,organism,application='knockout',source='ncbi',max_guides=20,min_score=30,genome_context=None,max_mismatches=3):
    a,d,s=fetch_sequence(gene_name,organism,source); return a,d,s,design_guides(s,application=application,max_guides=max_guides,min_score=min_score,genome_context=genome_context,max_mismatches=max_mismatches)
