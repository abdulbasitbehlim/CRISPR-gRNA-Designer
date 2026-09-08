#!/usr/bin/env python3
from __future__ import annotations
import json, os, re
from statistics import median
from typing import List
import pandas as pd
import streamlit as st
from grna_designer import (
    GuideRNA, analyze_offtargets, clean_dna_sequence, design_guides, fetch_sequence,
    guidescan2_available, parse_reference, run_guidescan2, score_breakdown, validate_guide,
)
APP_VERSION='3.0.0'
MAX_CUSTOM_BP=100_000
MAX_LOCAL_REFERENCE_BP=5_000_000
st.set_page_config(page_title='CRISPR Studio v3 | gRNA Designer',page_icon='🧬',layout='wide')
dark_mode=st.sidebar.toggle('Dark mode',value=True,key='dark_mode')
st.markdown('''<style>.block-container{max-width:1450px;padding-top:1.5rem}.hero{padding:1.6rem 1.8rem;border:1px solid #6b728055;border-radius:18px;margin-bottom:1rem}.seq{font-family:monospace;overflow-wrap:anywhere}</style>''',unsafe_allow_html=True)
st.sidebar.markdown(f'### CRISPR Studio v{APP_VERSION}')
max_guides=st.sidebar.slider('Maximum guides',5,50,20,5)
min_score=st.sidebar.slider('Minimum heuristic score',0,90,35,5)
max_mismatches=st.sidebar.slider('Off-target mismatches',0,4,3)
st.sidebar.caption('MIT/Hsu specificity is calculated natively. CFD and Doench Rule Set 2 are enabled when the optional GuideMaker scoring provider is installed. Whole-genome mode uses a local GuideScan2 index.')

@st.cache_data(ttl=3600,show_spinner=False)
def cached_fetch(gene,organism,source): return fetch_sequence(gene,organism,source)

def read_upload(f):
    if f is None:return ''
    b=f.getvalue()
    try:return b.decode('utf-8-sig')
    except UnicodeDecodeError:return b.decode('latin-1')

def safe_name(s):return re.sub(r'[^A-Za-z0-9._-]+','_',s).strip('_') or 'crispr_guides'

st.markdown('''<section class="hero"><h1>CRISPR Studio <b>v3</b></h1><p>SpCas9 guide design with MIT/Hsu specificity, optional Doench Rule Set 2 and CFD scoring, multi-contig local-reference screening, and GuideScan2 whole-genome support.</p></section>''',unsafe_allow_html=True)
st.subheader('Design workspace')
input_mode=st.radio('Target input',['Gene lookup','Paste sequence'],horizontal=True,label_visibility='collapsed')
with st.form('design_form',border=True):
    a,b=st.columns([1.25,1])
    with a:
        if input_mode=='Gene lookup':
            gene_name=st.text_input('Gene symbol',value='TP53'); organism=st.text_input('Organism',value='Homo sapiens')
            source=st.radio('Sequence database',['Ensembl','NCBI'],horizontal=True).lower(); custom_sequence=''
        else:
            custom_sequence=st.text_area('Target DNA or FASTA',height=180,placeholder='>target\nATG...')
            gene_name='Custom target'; organism='Not specified'; source='manual'
    with b:
        app_label=st.radio('Design intent',['Knockout','Knockdown / CRISPRi']); application='knockout' if app_label=='Knockout' else 'knockdown'
        screen_mode=st.selectbox('Specificity analysis',['None','Local reference','Whole genome (GuideScan2)'])
        ref_upload=None; ref_text=''; genome_index=''
        if screen_mode=='Local reference':
            ref_upload=st.file_uploader('Reference FASTA / text',type=['fa','fasta','fna','txt'])
            ref_text=st.text_area('Or paste reference sequence',height=90)
        elif screen_mode.startswith('Whole genome'):
            genome_index=st.text_input('GuideScan2 index path',value=os.getenv('GUIDESCAN_INDEX',''),help='Path to a prebuilt GuideScan2 genome index on the server/local machine.')
            st.caption('Backend detected: '+('yes' if guidescan2_available() else 'no'))
    submitted=st.form_submit_button('Design and analyze guides',type='primary',use_container_width=True)

if submitted:
    try:
        if input_mode=='Paste sequence':
            sequence=clean_dna_sequence(custom_sequence)
            if len(sequence)<50:raise ValueError('The target sequence must contain at least 50 bp.')
            if len(sequence)>MAX_CUSTOM_BP:raise ValueError(f'Target exceeds {MAX_CUSTOM_BP:,} bp limit.')
            accession='CUSTOM'; description='User-provided target sequence'
        else:
            if not gene_name.strip() or not organism.strip():raise ValueError('Enter both gene symbol and organism.')
            accession,description,sequence=cached_fetch(gene_name.strip(),organism.strip(),source)
        raw_reference=''; reference=None
        if screen_mode=='Local reference':
            raw_reference=read_upload(ref_upload) or ref_text
            reference=parse_reference(raw_reference)
            total=sum(map(len,reference.values()))
            if not total:raise ValueError('Local reference screening is enabled but no reference was supplied.')
            if total>MAX_LOCAL_REFERENCE_BP:raise ValueError(f'Local reference exceeds {MAX_LOCAL_REFERENCE_BP:,} bp hosted limit.')
        guides=design_guides(sequence,application=application,max_guides=max_guides,min_score=float(min_score),genome_context=reference,max_mismatches=max_mismatches)
        genome_rows=[]
        if screen_mode.startswith('Whole genome') and guides:
            genome_rows=run_guidescan2(guides,genome_index,max_mismatches=max_mismatches)
        st.session_state['analysis']={'guides':guides,'sequence':sequence,'reference':reference,'reference_raw':raw_reference,'genome_rows':genome_rows,
            'gene':gene_name.strip() or 'Custom target','organism':organism.strip(),'source':source,'application':application,'accession':accession,'description':description,'screen_mode':screen_mode,'max_mismatches':max_mismatches}
        st.success(f'Analysis complete: {len(guides)} guide candidates passed the filter.')
    except Exception as exc: st.error(f'Could not complete the analysis: {exc}')

if 'analysis' in st.session_state:
    r=st.session_state['analysis']; guides:List[GuideRNA]=r['guides']; sequence=r['sequence']; reference=r['reference']
    st.divider(); st.subheader(f"Analysis report · {r['gene']}"); st.caption(f"{r['description'][:220]} | Accession: {r['accession']} | Input length: {len(sequence):,} bp")
    if r['source']!='manual': st.warning('Gene lookup currently designs from a representative transcript/cDNA. For experimental use, map shortlisted spacers back to the intended genome assembly/exon and reject exon-junction-only sequences.')
    if not guides: st.warning('No candidates passed the current filter.')
    else:
        df=pd.DataFrame([g.to_dict() for g in guides]); df.insert(0,'Rank',range(1,len(df)+1))
        m=st.columns(5); m[0].metric('Guides retained',len(guides)); m[1].metric('Best heuristic',f'{max(g.score for g in guides):.1f}'); m[2].metric('Median GC',f'{median(g.gc_content for g in guides):.1f}%')
        ds=[g.doench_score for g in guides if g.doench_score is not None]; m[3].metric('Best Doench RS2',f'{max(ds):.1f}' if ds else 'Provider not installed')
        ms=[g.specificity_score for g in guides if g.specificity_score is not None]; m[4].metric('Best MIT specificity',f'{max(ms):.1f}' if ms else ('GuideScan2 run' if r['genome_rows'] else 'Not run'))
        tabs=st.tabs(['Ranked guides','Design landscape','Guide details','Off-target screen','Export','Methods & limits'])
        with tabs[0]: st.dataframe(df,use_container_width=True,hide_index=True)
        with tabs[1]:
            st.markdown('#### Candidate positions'); st.scatter_chart(df[['Rank','Start','Heuristic','GC%','Strand']],x='Start',y='Heuristic',size='GC%',color='Strand')
        with tabs[2]:
            idx=st.selectbox('Inspect a guide',range(len(guides)),format_func=lambda i:f'#{i+1} · {guides[i].sequence}')
            g=guides[idx]; st.markdown(f'<div class="seq">5′—{g.sequence}<b>{g.pam}</b>—3′</div>',unsafe_allow_html=True)
            q=st.columns(5); q[0].metric('Heuristic',f'{g.score:.1f}'); q[1].metric('Doench RS2',f'{g.doench_score:.1f}' if g.doench_score is not None else 'N/A'); q[2].metric('GC',f'{g.gc_content:.1f}%'); q[3].metric('MIT',f'{g.specificity_score:.1f}' if g.specificity_score is not None else 'N/A'); q[4].metric('CFD spec.',f'{g.cfd_specificity:.1f}' if g.cfd_specificity is not None else 'N/A')
            st.json(validate_guide(g)); st.markdown('##### Legacy heuristic breakdown'); st.json(score_breakdown(g.sequence,g.pam,application=g.application,start=g.start,sequence_length=len(sequence)))
        with tabs[3]:
            if reference:
                oi=st.selectbox('Inspect local off-targets',range(len(guides)),key='otguide',format_func=lambda i:f'#{i+1} · {guides[i].sequence}')
                rep=analyze_offtargets(guides[oi].sequence,reference,max_mismatches=r['max_mismatches'])
                c=st.columns(5); c[0].metric('MIT specificity',f'{rep.mit_specificity:.1f}'); c[1].metric('CFD specificity',f'{rep.cfd_specificity:.1f}' if rep.cfd_specificity is not None else 'Provider N/A'); c[2].metric('Near-matches',len(rep.hits)); c[3].metric('PAM sites',rep.pam_sites_scanned); c[4].metric('Contigs',rep.reference_contigs)
                st.dataframe(pd.DataFrame([h.to_dict() for h in rep.hits]),use_container_width=True,hide_index=True) if rep.hits else st.success('No additional PAM-compatible near-matches within the selected mismatch radius.')
            elif r['genome_rows']:
                st.markdown('#### GuideScan2 whole-genome results'); st.dataframe(pd.DataFrame(r['genome_rows']),use_container_width=True,hide_index=True)
            else: st.info('No specificity analysis was run. Choose Local reference or Whole genome (GuideScan2).')
        with tabs[4]:
            base=safe_name(r['gene']); csvb=df.to_csv(index=False).encode(); fasta='\n'.join(f'>{base}_g{i}|pam={g.pam}|strand={g.strand}\n{g.sequence}' for i,g in enumerate(guides,1)).encode(); js=json.dumps({'metadata':{k:v for k,v in r.items() if k not in {'guides','sequence','reference','reference_raw','genome_rows'}},'guides':df.where(pd.notna(df),None).to_dict('records')},indent=2).encode()
            d=st.columns(3); d[0].download_button('Download CSV',csvb,file_name=f'{base}_v3_guides.csv'); d[1].download_button('Download FASTA',fasta,file_name=f'{base}_v3_spacers.fasta'); d[2].download_button('Download JSON',js,file_name=f'{base}_v3_analysis.json')
        with tabs[5]:
            st.markdown('''#### v3 scoring
- **MIT/Hsu** pair scores and guide-level specificity are calculated natively from the published positional mismatch weights.
- **Doench Rule Set 2** and **CFD** use the optional USDA GuideMaker implementation when installed; v3 never relabels the legacy heuristic as Doench.
- **Local reference** preserves FASTA contigs and never concatenates chromosome/contig boundaries.
- **Whole genome** uses a prebuilt **GuideScan2** index and is intentionally separated from local screening.

#### Remaining limitations
- Gene lookup still begins from representative cDNA; validate genomic exon/assembly coordinates before experimental use.
- CRISPRi TSS-aware ranking, variant-aware filtering, bulges, and chromatin annotations require external genomic annotation/index resources.
- A clean in-silico score does not guarantee activity or safety; experimental validation remains required.''')
