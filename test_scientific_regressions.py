# ============================================================================
# SCIENTIFIC REGRESSION TESTS
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Protects important scientific behaviours with regression tests so future edits do not change validated results unexpectedly.
#
# HOW TO READ THIS FILE:
# 1. Each test describes one behaviour the scientific code must continue to satisfy.
# 2. Test inputs are prepared first, then the function is called, then assertions
#    check that the output still matches the expected scientific behaviour.
# 3. The assertions and scientific expectations are intentionally unchanged.
#
# MAIN TOP-LEVEL TEST FUNCTIONS:
# - function: test_malformed_or_multiple_target_records_rejected
# - function: test_ambiguous_guides_and_pams_never_emitted
# - function: test_reverse_spacer_pam_and_cut_coordinates
# - function: test_no_synthetic_contig_junction_target
# - function: test_duplicate_reference_ids_rejected
# - function: test_first_exact_match_is_not_silently_excluded
# - function: test_explicit_target_exclusion_and_retained_duplicate
# - function: test_declared_but_missing_locus_rejected
# - function: test_all_hits_scored_even_when_details_are_capped
# - function: test_ambiguous_reference_and_missing_exact_target_reported
# - function: test_invalid_mismatch_radius_rejected
# - function: test_unscreened_and_unknown_target_are_not_passed
# - function: test_zero_mismatch_screen_cannot_receive_full_local_pass
# - function: test_high_risk_single_mismatch_requires_review
# - function: test_complete_clean_local_screen_can_receive_full_local_pass
# - function: test_cfd_golden_exact_and_pam_weight
# - function: test_context30_orientation_and_edge_behavior
# - function: genomic_fixture
# - function: test_genomic_cds_filter_and_absolute_coordinates
# - function: test_crispri_restricts_repressor_scope
# - function: test_crispri_filters_before_global_candidate_truncation
# - function: test_missing_canonical_transcript_does_not_guess
# - function: test_translation_id_rejected_before_protein_is_normalized
# - function: test_accession_version_mismatch_never_falls_back
# - function: test_json_is_standard_and_keeps_evidence
# - function: test_user_twenty_copy_example
# - function: test_user_spliced_exon_junction_is_not_genomic_candidate
# - function: test_declared_region_inside_contig_maps_both_guide_strands
# - function: test_wrong_declared_region_rejected_before_result
# - function: test_ambiguous_reference_never_gets_local_checks_met
# - plus 8 additional tests/helpers
# ============================================================================

"""Regression tests for scientific failures reproduced against v3.2.0."""
from dataclasses import replace
import json
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from Bio import SeqIO
from io import StringIO

import accession_lookup
from grna_designer import (clean_dna_sequence, parse_reference, design_guides,
    analyze_offtargets, ReferenceIndex, GuideRNA, validate_guide, _context30,
    screen_guides)
from knockout import context_from_genbank, design_knockout_guides
from crispri import design_crispri_from_sequence, design_crispri_guides, TSSContext, _canonical_transcript
from models import cfd_score
from reporting import export_json, make_metadata

SPACER = 'GCTAGCTAGCTAGCTAGCTA'
TARGET = SPACER + 'AGG'


@pytest.mark.parametrize('text', ['>a\nACGT\n>b\nGGGG', '>a\nACGT\n>a\nGGGG', 'ACGT\n>a\nGGGG', '>a\n'])
def test_malformed_or_multiple_target_records_rejected(text):
    with pytest.raises(ValueError):
        clean_dna_sequence(text)


def test_ambiguous_guides_and_pams_never_emitted():
    assert not design_guides('N' * 20 + 'AGG', min_score=0)
    assert not design_guides(SPACER + 'NGG', min_score=0)
    assert not design_guides('CCN' + 'N' * 20, min_score=0)


def test_reverse_spacer_pam_and_cut_coordinates():
    target = str(Seq(TARGET).reverse_complement())
    guide = next(g for g in design_guides('AAAA' + target + 'TTTT', min_score=0) if g.strand == '-')
    assert guide.sequence == SPACER and guide.pam == 'AGG'
    assert (guide.start, guide.end) == (4, 27)
    assert (guide.spacer_start, guide.spacer_end, guide.cut_boundary) == (7, 27, 10)


def test_no_synthetic_contig_junction_target():
    report = analyze_offtargets('A' * 20, '>left\n' + 'A' * 20 + '\n>right\nAGG')
    assert report.total_hits == 0


def test_duplicate_reference_ids_rejected():
    with pytest.raises(ValueError, match='Duplicate'):
        parse_reference('>same\nAAAA\n>same\nCCCC')


def test_first_exact_match_is_not_silently_excluded():
    report = analyze_offtargets(SPACER, TARGET)
    assert not report.on_target_excluded and report.total_hits == 1
    assert report.exact_matches == 1 and report.specificity_score == 50


def test_explicit_target_exclusion_and_retained_duplicate():
    report = analyze_offtargets(SPACER, {'a': TARGET, 'b': TARGET}, intended_target=('b', 0, '+'))
    assert report.on_target_excluded and report.total_hits == 1
    assert report.hits[0].contig == 'a' and report.exact_matches == 2
    assert report.specificity_score == 50 and report.cfd_specificity == 50


def test_declared_but_missing_locus_rejected():
    with pytest.raises(ValueError, match='Declared'):
        analyze_offtargets(SPACER, TARGET, intended_target=('wrong', 0, '+'))


def test_all_hits_scored_even_when_details_are_capped():
    index = ReferenceIndex({f'copy{i}': TARGET for i in range(301)})
    small = analyze_offtargets(SPACER, index, max_hits=1, intended_target=('copy0', 0, '+'))
    full = analyze_offtargets(SPACER, index, max_hits=500, intended_target=('copy0', 0, '+'))
    assert small.total_hits == full.total_hits == 300
    assert small.specificity_score == full.specificity_score == round(100/301, 2)
    assert small.cfd_specificity == full.cfd_specificity
    assert small.risk_counts == full.risk_counts == {'Critical':300, 'High':0, 'Moderate':0, 'Low':0}
    assert small.max_mit_risk == full.max_mit_risk == 1.0
    assert small.max_cfd_risk == full.max_cfd_risk == 1.0
    assert small.truncated and len(small.hits) == 1 and not full.truncated


def test_ambiguous_reference_and_missing_exact_target_reported():
    report = analyze_offtargets(SPACER, TARGET.replace('G', 'N', 1))
    assert report.ambiguous_bases == 1 and report.intended_target_status == 'no_exact_match'


@pytest.mark.parametrize('mm', [-1, 5, 1.5])
def test_invalid_mismatch_radius_rejected(mm):
    with pytest.raises(ValueError):
        analyze_offtargets(SPACER, TARGET, max_mismatches=mm)


def test_unscreened_and_unknown_target_are_not_passed():
    g = GuideRNA(SPACER, 'AGG', '+', 0, 23, 50, 70)
    assert validate_guide(g)['sequence_checks_pass']
    assert not validate_guide(g)['overall_pass']
    assert not validate_guide(g, TARGET)['overall_pass']


def test_zero_mismatch_screen_cannot_receive_full_local_pass():
    guide = GuideRNA(SPACER, 'AGG', '+', 0, 23, 50, 68)
    screen_guides(
        [guide], {'target': TARGET}, max_mismatches=0,
        intended_targets={(0, '+'): ('target', 0, '+')},
    )
    checks = validate_guide(guide, min_score=35)
    assert guide.specificity_score == 100
    assert guide.screened_mismatch_radius == 0
    assert checks['locus_evidence_ok'] and checks['aggregate_specificity_ok']
    assert not checks['screen_scope_complete']
    assert not checks['specificity_ok'] and not checks['overall_pass']


def test_high_risk_single_mismatch_requires_review():
    high_risk = 'A' + SPACER[1:] + 'AGG'
    guide = GuideRNA(SPACER, 'AGG', '+', 0, 23, 50, 68)
    screen_guides(
        [guide], {'target': TARGET, 'off_target': high_risk}, max_mismatches=3,
        intended_targets={(0, '+'): ('target', 0, '+')},
    )
    checks = validate_guide(guide, min_score=35)
    assert guide.specificity_score == 50
    assert guide.cfd_specificity == pytest.approx(52.63)
    assert guide.risk_counts == {'Critical':0, 'High':1, 'Moderate':0, 'Low':0}
    assert guide.max_mit_risk == 1.0
    assert guide.max_cfd_risk == pytest.approx(0.9)
    assert checks['screen_scope_complete'] and checks['aggregate_specificity_ok']
    assert not checks['no_critical_or_high_hits']
    assert not checks['specificity_ok'] and not checks['overall_pass']


def test_complete_clean_local_screen_can_receive_full_local_pass():
    guide = GuideRNA(SPACER, 'AGG', '+', 0, 23, 50, 68)
    screen_guides(
        [guide], {'target': TARGET}, max_mismatches=3,
        intended_targets={(0, '+'): ('target', 0, '+')},
    )
    checks = validate_guide(guide, min_score=35)
    assert guide.screened_mismatch_radius == 3
    assert guide.risk_counts == {'Critical':0, 'High':0, 'Moderate':0, 'Low':0}
    assert checks['specificity_ok'] and checks['overall_pass']


def test_cfd_golden_exact_and_pam_weight():
    assert cfd_score(SPACER, SPACER) == 1
    assert cfd_score(SPACER, SPACER, 'AAG') == pytest.approx(0.259259259)


def test_context30_orientation_and_edge_behavior():
    sequence = 'ACGT' + TARGET + 'TGA'
    assert _context30(sequence, 4, 27, '+') == sequence
    rc = str(Seq(sequence).reverse_complement())
    assert _context30(rc, 3, 26, '-') == sequence
    assert _context30(TARGET, 0, 23, '+') is None


def genomic_fixture(strand=1):
    seq = 'A' * 10 + TARGET + 'A' * 17 + TARGET + 'A' * 17 + TARGET + 'A' * 15
    rec = SeqRecord(Seq(seq), id='NC_TEST.1', name='fixture', description='synthetic regression fixture')
    rec.annotations['molecule_type'] = 'DNA'
    parts = [SimpleLocation(0, 40, strand=strand), SimpleLocation(80, 120, strand=strand)]
    rec.features = [SeqFeature(CompoundLocation(parts), type='CDS', qualifiers={'gene': ['TEST'], 'protein_id': ['NP_TEST.1']})]
    handle = StringIO(); SeqIO.write(rec, handle, 'genbank')
    return context_from_genbank(handle.getvalue(), 'TEST', 'Synthetic', '1', 'NC_TEST.1', 'test', 1001)


@pytest.mark.parametrize('strand', [1, -1])
def test_genomic_cds_filter_and_absolute_coordinates(strand):
    context = genomic_fixture(strand)
    guides = design_knockout_guides(context, min_score=0, max_guides=100)
    starts = {g.start for g in guides if g.sequence == SPACER}
    assert starts == {10, 90}  # middle target is intronic
    for g in guides:
        assert g.annotation['Genomic cut after base'] == 1000 + g.cut_boundary
        assert any(a < g.cut_boundary < b for a, b in context.isoforms[0].segments)
    first = next(g for g in guides if g.start == 10)
    assert first.annotation['CDS position percent'] == (33.75 if strand == 1 else 66.25)


def test_crispri_restricts_repressor_scope():
    context = TSSContext('TEST', 'plant', 'arabidopsis_thaliana', '', '', '', '', 1, 101, 1, 500, 100, 'A' * 500, '')
    with pytest.raises(ValueError, match='restricted'):
        design_crispri_guides(context)


def test_crispri_filters_before_global_candidate_truncation():
    seq = 'G' * 8000 + 'A' * 100 + TARGET + 'A' * 350
    result = design_crispri_from_sequence(seq, 8101, min_score=0)
    assert result and any(row.guide.sequence == SPACER for row in result)
    assert all(row.guide.doench_score is None for row in result)
    assert all(row.tss_distance % 1 == 0.5 for row in result)


def test_missing_canonical_transcript_does_not_guess():
    with pytest.raises(ValueError, match='canonical'):
        _canonical_transcript({'Transcript': [{'id': 'T1'}]})


def test_translation_id_rejected_before_protein_is_normalized(monkeypatch):
    monkeypatch.setattr(accession_lookup, '_ensembl_lookup', lambda _: {'object_type': 'Translation', 'id': 'ENSP00001'})
    monkeypatch.setattr(accession_lookup, '_ensembl_sequence', lambda *_: pytest.fail('must not fetch protein as DNA'))
    with pytest.raises(ValueError, match='protein'):
        accession_lookup.fetch_ensembl_accession('ENSP00001')


def test_accession_version_mismatch_never_falls_back(monkeypatch):
    def fake(path, params):
        if path == 'esearch.fcgi': return SimpleNamespace(json=lambda: {'esearchresult': {'idlist': ['1']}})
        if path == 'esummary.fcgi': return SimpleNamespace(json=lambda: {'result': {'1': {'slen': 23, 'biomol': 'genomic'}}})
        return SimpleNamespace(text='>NC_000001.2\n' + TARGET)
    monkeypatch.setattr(accession_lookup, '_ncbi_get', fake)
    with pytest.raises(ValueError, match='does not match'):
        accession_lookup.fetch_ncbi_accession('NC_000001.1')


def test_json_is_standard_and_keeps_evidence():
    g = GuideRNA(SPACER, 'AGG', '+', 0, 23, 50, 70)
    g.doench_score = float('nan')
    meta = make_metadata(TARGET, {'ref': TARGET}, mismatch_limit=3)
    text = export_json(meta, [g]).decode()
    assert 'NaN' not in text
    data = json.loads(text)
    assert len(data['metadata']['target_sha256']) == 64
    assert 'local_hit_details' in data and data['guides'][0]['Doench RS2'] is None
    evidence = data['local_hit_details'][0]
    assert 'screened_mismatch_radius' in evidence and 'risk_counts' in evidence
    assert 'maximum_per_site_mit_risk' in evidence and 'maximum_per_site_cfd_risk' in evidence


@pytest.mark.parametrize('limit', [0, 1, 10, 100])
@pytest.mark.parametrize('declare', [False, True])
def test_user_twenty_copy_example(limit, declare):
    reference = {f'copy{i}': TARGET for i in range(20)}
    report = analyze_offtargets(SPACER, reference, max_hits=limit,
        intended_target=('copy0', 0, '+') if declare else None)
    assert report.total_hits == (19 if declare else 20)
    assert report.specificity_score == (5.0 if declare else 4.76)
    assert report.cfd_specificity == report.specificity_score
    assert len(report.hits) == min(limit, report.total_hits)


def test_user_spliced_exon_junction_is_not_genomic_candidate():
    from knockout import KnockoutContext, CodingIsoform
    exon1, intron, exon2 = 'A' * 20, 'T' * 40, 'AGG' + 'C' * 25
    assert any(g.sequence == 'A' * 20 for g in design_guides(exon1 + exon2, min_score=0))
    context = KnockoutContext('TEST', 'synthetic', '1', 'TEST.1', '1', 1,
        exon1 + intron + exon2, (CodingIsoform('P1', (), 1, ((0, 20), (60, 88))),), 'P1', 'fixture', 'fixture')
    guides = design_knockout_guides(context, min_score=0)
    assert not any(g.sequence == 'A' * 20 and g.pam == 'AGG' for g in guides)
    for g in guides:
        target = context.sequence[g.start:g.end]
        assert g.full_target == (target if g.strand == '+' else str(Seq(target).reverse_complement()))


@pytest.mark.parametrize('orientation', ['+', '-'])
def test_declared_region_inside_contig_maps_both_guide_strands(orientation):
    from grna_designer import screen_guides, TargetLocus
    sequence = 'AAAA' + TARGET + 'AAAA' + str(Seq(TARGET).reverse_complement()) + 'AAAA'
    reference_sequence = sequence if orientation == '+' else str(Seq(sequence).reverse_complement())
    reference = {'chr1': 'T' * 13 + reference_sequence + 'A' * 10}
    guides = design_guides(sequence, min_score=0, max_guides=100)
    screen_guides(guides, reference, intended_region=('chr1', 13, orientation), target_sequence=sequence)
    assert {g.strand for g in guides} == {'+', '-'}
    for g in guides:
        assert g.screen_status == 'verified_locus'
        assert isinstance(g.intended_target, TargetLocus)
        locus = g.intended_target
        site = reference['chr1'][locus.start:locus.start + 23]
        assert (site if locus.strand == '+' else str(Seq(site).reverse_complement())) == g.full_target


def test_wrong_declared_region_rejected_before_result():
    with pytest.raises(ValueError, match='does not exactly match'):
        design_guides(TARGET, min_score=0, genome_context={'chr1': TARGET}, intended_region=('chr1', 1, '+'))


def test_ambiguous_reference_never_gets_local_checks_met():
    from grna_designer import screen_guides
    g = GuideRNA(SPACER, 'AGG', '+', 0, 23, 50, 70)
    screen_guides([g], {'chr1': TARGET + 'N'}, intended_targets={(0, '+'): ('chr1', 0, '+')})
    assert g.screen_status == 'verified_locus'
    assert not validate_guide(g)['overall_pass']


def test_ncbi_rna_rejected_before_sequence_download(monkeypatch):
    def fake(path, params):
        if path == 'esearch.fcgi': return SimpleNamespace(json=lambda: {'esearchresult': {'idlist': ['1']}})
        if path == 'esummary.fcgi': return SimpleNamespace(json=lambda: {'result': {'1': {'slen': 23, 'biomol': 'mRNA'}}})
        pytest.fail('RNA must not be downloaded for knockout design')
    monkeypatch.setattr(accession_lookup, '_ncbi_get', fake)
    with pytest.raises(ValueError, match='RNA/cDNA'):
        accession_lookup.fetch_ncbi_accession('NM_000546.6')


@pytest.mark.parametrize('kind', ['gene', 'transcript', 'exon'])
def test_ensembl_accession_fetches_genomic_never_cdna(monkeypatch, kind):
    monkeypatch.setattr(accession_lookup, '_ensembl_lookup', lambda _: {'id': 'ENST00001', 'object_type': kind, 'start': 1, 'end': 23, 'version': 1})
    def fetch(stable_id, sequence_type):
        assert sequence_type == 'genomic'
        return TARGET
    monkeypatch.setattr(accession_lookup, '_ensembl_sequence', fetch)
    record = accession_lookup.fetch_ensembl_accession('ENST00001.1')
    assert 'genomic DNA' in record.object_type and record.sequence == TARGET


def test_gene_design_api_cannot_use_old_transcript_path(monkeypatch):
    import grna_designer
    import knockout
    context = genomic_fixture()
    monkeypatch.setattr(knockout, 'fetch_ncbi_knockout_context', lambda *a: context)
    monkeypatch.setattr(grna_designer, 'fetch_sequence', lambda *a: pytest.fail('must never use spliced transcript'))
    accession, _, sequence, guides = grna_designer.design_from_gene('TEST', 'synthetic', min_score=0)
    assert accession == context.accession and sequence == context.sequence and guides
    with pytest.raises(ValueError, match='NCBI genomic'):
        grna_designer.design_from_gene('TEST', 'synthetic', source='ensembl')


def test_cfd_nonperfect_published_weight():
    # G->A at position 1 uses the published rG:dT,1 mismatch weight.
    assert cfd_score('G' + 'A' * 19, 'A' * 20) == pytest.approx(0.9)


def test_guidescan_contract_numeric_position_and_multiple_hit_rows(monkeypatch):
    import csv
    import grna_designer
    from pathlib import Path
    monkeypatch.setattr(grna_designer, 'guidescan2_available', lambda: True)
    def run(command, **kwargs):
        rows = list(csv.DictReader(Path(command[command.index('-f') + 1]).open()))
        assert int(rows[0]['position']) == 1
        assert rows[0]['chromosome'] == '__unmapped__'
        output = Path(command[command.index('-o') + 1])
        output.write_text('id,sequence,specificity\ng1,' + SPACER + ',0.5\ng1,' + SPACER + ',0.5\n')
        return SimpleNamespace(returncode=0, stderr='')
    monkeypatch.setattr(grna_designer.subprocess, 'run', run)
    assert len(grna_designer.run_guidescan2([GuideRNA(SPACER, 'AGG', '+', 0, 23, 50, 70)], 'test-index')) == 2


def test_guidescan_missing_output_is_error(monkeypatch):
    import grna_designer
    monkeypatch.setattr(grna_designer, 'guidescan2_available', lambda: True)
    monkeypatch.setattr(grna_designer.subprocess, 'run', lambda *a, **k: SimpleNamespace(returncode=0, stderr=''))
    with pytest.raises(RuntimeError, match='no output'):
        grna_designer.run_guidescan2([GuideRNA(SPACER, 'AGG', '+', 0, 23, 50, 70)], 'test-index')


@pytest.mark.parametrize('gene', ['ACTB', 'GAPDH'])
def test_saved_ncbi_genomic_examples(gene):
    from pathlib import Path
    from knockout import KnockoutContext, CodingIsoform
    data = json.loads((Path(__file__).parent / 'tests' / 'fixtures' / f'{gene}_NCBI_context.json').read_text())
    data['isoforms'] = tuple(CodingIsoform(x['protein_id'], tuple(x['transcript_ids']), x['strand'], tuple(map(tuple, x['segments']))) for x in data['isoforms'])
    context = KnockoutContext(**data)
    guides = design_knockout_guides(context, max_guides=20, min_score=30)
    assert len(guides) == 20
    selected = next(x for x in context.isoforms if x.protein_id == context.selected_protein)
    for guide in guides:
        target = context.sequence[guide.start:guide.end]
        assert guide.full_target == (target if guide.strand == '+' else str(Seq(target).reverse_complement()))
        assert any(a < guide.cut_boundary < b for a, b in selected.segments)


def test_rs2_provider_contract_and_missing_context(monkeypatch):
    import models
    calls = []
    def predict(array, num_threads):
        calls.append((array.tolist(), num_threads))
        return [1.2]  # regression output is not clipped to a probability
    monkeypatch.setattr(models, 'rs2_provider', lambda: predict)
    context = 'ACGT' + TARGET + 'TGA'
    assert models.rs2_score(context) == 120.0
    assert calls == [([context], 1)]
    assert models.rs2_score(None) is None
    monkeypatch.setattr(models, 'rs2_provider', lambda: None)
    assert models.rs2_score(context) is None
