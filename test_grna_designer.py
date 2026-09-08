import pytest
from grna_designer import (
    GuideRNA, analyze_offtargets, clean_dna_sequence, design_guides, validate_guide,
    _gc_content, _has_homopolymer, _doench_like_score, calculate_offtarget_score,
    mit_offtarget_score, mit_specificity, parse_reference,
)
from crispri import (
    TSSContext, design_crispri_guides, design_crispri_from_sequence, tss_band,
    _genomic_interval, _genomic_guide_strand,
)
from accession_lookup import detect_database, normalize_accession, fetch_accession

SAMPLE_SEQUENCE=(
    'ATGGCTAGCTAGCTAGGCTAGCATCGATCGATCGGATCGATCGATCGATCGGCTAGCTAGCTAGCTAGG'
    'CTAGCATCGATCGATCGGATCGATCGATCGATCGGCTAGCTAGCTAGCTAGGCTAGCATCGATCGATCG'
    'GATCGATCGATCGATCGGCTAGCTAGCTAGCTAGGCTAGCATCGATCGATCGGATCGATCGATCGATCG'
)

def test_clean_fasta_and_rna():
    assert clean_dna_sequence('>x\nAUGC RYSW\n12')=='ATGCNNNN'

def test_gc_content():
    assert _gc_content('GCGCATAT')==50.0

def test_homopolymer():
    assert _has_homopolymer('ACGTGGGGACGT',4)

def test_legacy_score_range():
    assert 0<=_doench_like_score('A'*20,'AGG')<=100

def test_design_guides():
    guides=design_guides(SAMPLE_SEQUENCE,min_score=0,max_guides=5)
    assert guides
    assert all(isinstance(g,GuideRNA) for g in guides)
    assert all(len(g.sequence)==20 and g.pam.endswith('GG') for g in guides)

def test_validate():
    g=GuideRNA('ACGTACGTACGTACGTACGT','AGG','+',0,23,50,70)
    assert validate_guide(g)['overall_pass']

def test_mit_exact_pair_is_one():
    assert mit_offtarget_score('GCTAGCTAGCTAGCTAGCTA','GCTAGCTAGCTAGCTAGCTA')==1.0

def test_mit_mismatch_reduces_pair_score():
    s='GCTAGCTAGCTAGCTAGCTA'; o='GCTAGCTAGCTAACTAGCTA'
    assert 0<mit_offtarget_score(s,o)<1

def test_mit_aggregate():
    assert mit_specificity([])==100.0
    assert mit_specificity([1.0])==50.0

def test_multifasta_preserves_contigs():
    assert parse_reference('>chr1\nAAAA\n>chr2\nCCCC')=={'chr1':'AAAA','chr2':'CCCC'}

def test_exact_intended_target_excluded():
    s='GCTAGCTAGCTAGCTAGCTA'; r=analyze_offtargets(s,s+'AGG')
    assert r.on_target_excluded and not r.hits and r.specificity_score==100.0

def test_additional_exact_copy_is_critical():
    s='GCTAGCTAGCTAGCTAGCTA'; r=analyze_offtargets(s,(s+'AGG')*2)
    assert len(r.hits)==1 and r.hits[0].risk=='Critical' and r.specificity_score<=50.0

def test_offtarget_wrapper_empty():
    n,score,notes=calculate_offtarget_score('GCTAGCTAGCTAGCTAGCTA','AGG',None)
    assert n==0 and score==0 and 'Genome not provided' in notes

def test_local_report_has_mit_fields():
    s='GCTAGCTAGCTAGCTAGCTA'; near='GCTAGCTAGCTAACTAGCTAAGG'
    r=analyze_offtargets(s,near)
    assert r.hits and r.hits[0].mit_score>0


def _ctx(strand=1):
    return TSSContext(
        gene='TEST', organism='Human', species='homo_sapiens', gene_id='G1', transcript_id='T1',
        assembly='GRCh38', chromosome='1', transcript_strand=strand, tss_coordinate=1000,
        region_start=900 if strand==1 else 650, region_end=1350 if strand==1 else 1100,
        tss_offset=100, sequence='A'*451, description='test'
    )

def test_tss_bands():
    assert tss_band(75)[0].startswith('Preferred')
    assert tss_band(25)[0]=='High-priority'
    assert tss_band(250)[0]=='CRISPRi window'

def test_plus_strand_genomic_interval_mapping():
    assert _genomic_interval(_ctx(1), 100, 120)==(1000,1019)
    assert _genomic_guide_strand(_ctx(1), '+')=='+'

def test_minus_strand_genomic_interval_mapping():
    c=_ctx(-1)
    assert _genomic_interval(c,100,120)==(981,1000)
    assert _genomic_guide_strand(c,'+')=='-'

def test_custom_crispri_requires_valid_tss_position():
    with pytest.raises(ValueError):
        design_crispri_from_sequence(SAMPLE_SEQUENCE,0)
    with pytest.raises(ValueError):
        design_crispri_from_sequence(SAMPLE_SEQUENCE,len(SAMPLE_SEQUENCE)+1)

def test_crispri_filters_to_tss_window_and_annotates():
    seq='A'*90 + 'GCTAGCTAGCTAGCTAGCTAAGG' + 'A'*340
    c=TSSContext('TEST','Human','homo_sapiens','G1','T1','GRCh38','1',1,1000,900,1352,100,seq,'test')
    rows=design_crispri_guides(c,max_guides=20,min_score=0)
    assert rows
    assert all(-50 <= x.tss_distance <= 300 for x in rows)
    assert all(x.guide.application=='crispri' for x in rows)
    assert all(x.transcript_id=='T1' for x in rows)


def test_accession_normalization_and_database_detection():
    assert normalize_accession('  NM_000546.6  ') == 'NM_000546.6'
    assert detect_database('ENST00000269305') == 'ensembl'
    assert detect_database('NM_000546.6') == 'ncbi'
    with pytest.raises(ValueError):
        normalize_accession('bad id')


def test_accession_router_uses_expected_backend(monkeypatch):
    calls = []
    class Dummy:
        accession='X'; sequence='ATGC'; database='dummy'; description='d'; record_url=''; object_type='nucleotide'
    def fake_ncbi(accession):
        calls.append(('ncbi', accession)); return Dummy()
    def fake_ensembl(accession):
        calls.append(('ensembl', accession)); return Dummy()
    monkeypatch.setattr('accession_lookup.fetch_ncbi_accession', fake_ncbi)
    monkeypatch.setattr('accession_lookup.fetch_ensembl_accession', fake_ensembl)
    fetch_accession('NM_000546.6', 'Auto')
    fetch_accession('ENST00000269305', 'Auto')
    assert calls == [('ncbi','NM_000546.6'), ('ensembl','ENST00000269305')]
