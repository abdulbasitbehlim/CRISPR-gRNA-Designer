import pytest
from grna_designer import (
    GuideRNA, analyze_offtargets, clean_dna_sequence, design_guides, validate_guide,
    _gc_content, _has_homopolymer, _doench_like_score, calculate_offtarget_score,
    mit_offtarget_score, mit_specificity, parse_reference,
)

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
