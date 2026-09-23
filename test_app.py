# ============================================================================
# APP
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Tests the main application behaviour so interface and workflow changes do not silently break the tool.
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
# - function: test_dashboard_loads_without_exception
# - function: test_pasted_sequence_completes_full_report
# - function: test_short_sequence_shows_clear_error
# - function: test_duplicate_reference_is_clear_error_and_clears_old_results
# - function: test_local_declared_region_reaches_verified_status_in_app
# - function: test_gene_lookup_uses_genomic_context_in_app
# - function: test_ensembl_gene_knockout_cannot_rank_cdna
# ============================================================================

"""Offline Streamlit smoke tests for the main dashboard workflow."""

from pathlib import Path

from streamlit.testing.v1 import AppTest


APP_PATH = Path(__file__).with_name("app.py")
SAMPLE_SEQUENCE = (
    "ATGGCTAGCTAGCTAGGCTAGCATCGATCGATCGGATCGATCGATCGATCGGCTAGCTAGCTAGCTAGG"
    * 4
)



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: test_dashboard_loads_without_exception
# ----------------------------------------------------------------------------
def test_dashboard_loads_without_exception():
    app = AppTest.from_file(str(APP_PATH), default_timeout=30).run()

    assert not app.exception
    assert any(button.label == "Design and analyze guides" for button in app.button)
    assert any(toggle.label == "Dark mode" for toggle in app.toggle)



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: test_pasted_sequence_completes_full_report
# ----------------------------------------------------------------------------
def test_pasted_sequence_completes_full_report():
    app = AppTest.from_file(str(APP_PATH), default_timeout=30).run()
    app.radio[0].set_value("Paste sequence").run()
    app.text_area[0].input(SAMPLE_SEQUENCE)
    app.button[0].click().run()

    assert not app.exception
    assert len(app.metric) >= 5
    assert len(app.tabs) == 6
    assert any("Analysis complete" in message.value for message in app.success)



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: test_short_sequence_shows_clear_error
# ----------------------------------------------------------------------------
def test_short_sequence_shows_clear_error():
    app = AppTest.from_file(str(APP_PATH), default_timeout=30).run()
    app.radio[0].set_value("Paste sequence").run()
    app.text_area[0].input("ATGC")
    app.button[0].click().run()

    assert not app.exception
    assert any("at least 50 bp" in message.value for message in app.error)




# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: test_duplicate_reference_is_clear_error_and_clears_old_results
# ----------------------------------------------------------------------------
def test_duplicate_reference_is_clear_error_and_clears_old_results():
    app = AppTest.from_file(str(APP_PATH), default_timeout=30).run()
    app.radio[0].set_value('Paste sequence').run()
    app.text_area[0].input(SAMPLE_SEQUENCE)
    app.button[0].click().run()
    assert 'analysis' in app.session_state
    next(x for x in app.selectbox if x.label == 'Specificity analysis').set_value('Local reference').run()
    next(x for x in app.text_area if x.label == 'Or paste reference sequence').input('>chr1\n' + SAMPLE_SEQUENCE + '\n>chr1\nAAAA')
    app.button[0].click().run()
    assert not app.exception
    assert any('Duplicate FASTA' in x.value for x in app.error)
    assert 'analysis' not in app.session_state



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: test_local_declared_region_reaches_verified_status_in_app
# ----------------------------------------------------------------------------
def test_local_declared_region_reaches_verified_status_in_app():
    app = AppTest.from_file(str(APP_PATH), default_timeout=30).run()
    app.radio[0].set_value('Paste sequence').run()
    app.text_area[0].input(SAMPLE_SEQUENCE)
    next(x for x in app.selectbox if x.label == 'Specificity analysis').set_value('Local reference').run()
    next(x for x in app.text_area if x.label == 'Or paste reference sequence').input('>chr1\n' + 'A'*10 + SAMPLE_SEQUENCE)
    next(x for x in app.text_input if x.label == 'Intended target record ID optional').input('chr1')
    next(x for x in app.number_input if x.label == 'Target region start in reference (1-based left edge)').set_value(11)
    app.button[0].click().run()
    assert not app.exception and not app.error
    result = app.session_state['analysis']
    assert result['guides']
    assert all(g.screen_status == 'verified_locus' for g in result['guides'])
    assert result['metadata']['intended_region'] == ('chr1', 10, '+')



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: test_gene_lookup_uses_genomic_context_in_app
# ----------------------------------------------------------------------------
def test_gene_lookup_uses_genomic_context_in_app(monkeypatch):
    import knockout
    from test_scientific_regressions import genomic_fixture
    monkeypatch.setattr(knockout, 'fetch_ncbi_knockout_context', lambda *a, **k: genomic_fixture())
    app = AppTest.from_file(str(APP_PATH), default_timeout=30).run()
    app.button[0].click().run()
    assert not app.exception and not app.error
    result = app.session_state['analysis']
    assert result['knockout_context'] is not None
    assert all('Selected protein' in g.annotation for g in result['guides'])



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: test_ensembl_gene_knockout_cannot_rank_cdna
# ----------------------------------------------------------------------------
def test_ensembl_gene_knockout_cannot_rank_cdna():
    app = AppTest.from_file(str(APP_PATH), default_timeout=30).run()
    next(x for x in app.radio if x.label == 'Sequence database').set_value('Ensembl')
    app.button[0].click().run()
    assert not app.exception
    assert any('NCBI genomic CDS' in x.value for x in app.error)
    assert 'analysis' not in app.session_state
