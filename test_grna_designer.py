"""
Unit tests for grna_designer.py

These tests exercise the pure sequence-processing logic only
(PAM scanning, scoring, validation, and now off-target analysis).
They do not hit the network, so they run anywhere without an internet connection or NCBI/Ensembl access.

Run with:
    pytest -v
"""

import sys

import pytest

from grna_designer import (
    GuideRNA,
    analyze_offtargets,
    clean_dna_sequence,
    design_guides,
    validate_guide,
    _gc_content,
    _has_homopolymer,
    _doench_like_score,
    calculate_offtarget_score,
    score_breakdown,
)


# A synthetic 200 bp sequence with several NGG / CCN sites on both strands.
SAMPLE_SEQUENCE = (
    "ATGGCTAGCTAGCTAGGCTAGCATCGATCGATCGGATCGATCGATCGATCGGCTAGCTAGCTAGCTAGG"
    "CTAGCATCGATCGATCGGATCGATCGATCGATCGGCTAGCTAGCTAGCTAGGCTAGCATCGATCGATCG"
    "GATCGATCGATCGATCGGCTAGCTAGCTAGCTAGGCTAGCATCGATCGATCGGATCGATCGATCGATCG"
)


class TestSequenceHelpers:
    def test_clean_fasta_and_rna(self):
        raw = ">example\nAUGC RYSW\n12"
        assert clean_dna_sequence(raw) == "ATGCNNNN"

    def test_clean_sequence_rejects_unexpected_punctuation(self):
        with pytest.raises(ValueError, match="unsupported character"):
            clean_dna_sequence("ATGC!")

    def test_gc_content_all_gc(self):
        assert _gc_content("GGGGCCCC") == 100.0

    def test_gc_content_all_at(self):
        assert _gc_content("AAAATTTT") == 0.0

    def test_gc_content_mixed(self):
        assert _gc_content("GCGCGCGCGCATATATATAT") == 50.0

    def test_gc_content_empty_string(self):
        assert _gc_content("") == 0.0

    def test_homopolymer_detected(self):
        assert _has_homopolymer("ACGTGGGGACGT", max_run=4) is True

    def test_homopolymer_not_detected(self):
        assert _has_homopolymer("ACGTGGGACGT", max_run=4) is False


class TestScoring:
    def test_score_in_valid_range(self):
        score = _doench_like_score("A" * 20, "AGG")
        assert 0.0 <= score <= 100.0

    def test_score_rejects_wrong_length(self):
        assert _doench_like_score("ACGT", "AGG") == 0.0

    def test_good_gc_and_terminal_g_scores_higher_than_poor_spacer(self):
        good_spacer = "ACGTACGTACGTACGTACGG"[:20]
        poor_spacer = "TTTTTTTTTTTTTTTTTTTT"
        good_score = _doench_like_score(good_spacer, "AGG")
        poor_score = _doench_like_score(poor_spacer, "AGG")
        assert good_score > poor_score

    def test_application_adjustment_is_clamped_to_100(self):
        components = score_breakdown(
            "GCGCGCGCGCGCGCGCGCGG",
            "CGG",
            application="knockdown",
            start=5,
            sequence_length=1000,
        )
        assert 0.0 <= components["Final score"] <= 100.0


class TestDesignGuides:
    def test_returns_guide_rna_objects(self):
        guides = design_guides(SAMPLE_SEQUENCE, application="knockout", min_score=0)
        assert len(guides) > 0
        assert all(isinstance(g, GuideRNA) for g in guides)

    def test_guides_sorted_descending_by_score(self):
        guides = design_guides(SAMPLE_SEQUENCE, application="knockout", min_score=0)
        scores = [g.score for g in guides]
        assert scores == sorted(scores, reverse=True)

    def test_respects_max_guides(self):
        guides = design_guides(SAMPLE_SEQUENCE, min_score=0, max_guides=3)
        assert len(guides) <= 3

    def test_all_spacers_are_20nt(self):
        guides = design_guides(SAMPLE_SEQUENCE, min_score=0)
        assert all(len(g.sequence) == 20 for g in guides)

    def test_all_pams_end_in_gg(self):
        guides = design_guides(SAMPLE_SEQUENCE, min_score=0)
        assert all(g.pam.upper().endswith("GG") for g in guides)

    def test_knockdown_biases_toward_5_prime(self):
        long_seq = SAMPLE_SEQUENCE * 4
        guides = design_guides(long_seq, application="knockdown", min_score=0, max_guides=5)
        assert len(guides) > 0
        assert guides[0].start < len(long_seq) * 0.6

    def test_empty_sequence_returns_no_guides(self):
        assert design_guides("", min_score=0) == []

    def test_sequence_with_no_pam_returns_no_guides(self):
        assert design_guides("A" * 100, min_score=0) == []


class TestValidateGuide:
    def test_well_formed_guide_passes(self):
        guide = GuideRNA(sequence="ACGTACGTACGTACGTACGG", pam="AGG", strand="+", start=0, end=23, gc_content=55.0, score=75.0)
        guide.sequence = "ACGTACGTACGTACGTACGG"[:20]
        checks = validate_guide(guide)
        assert checks["length_ok"] is True
        assert checks["pam_ok"] is True
        assert checks["overall_pass"] is True

    def test_bad_pam_fails(self):
        guide = GuideRNA(sequence="A" * 20, pam="ATT", strand="+", start=0, end=23, gc_content=50.0, score=75.0)
        checks = validate_guide(guide)
        assert checks["pam_ok"] is False
        assert checks["overall_pass"] is False

    def test_low_score_fails_overall(self):
        guide = GuideRNA(sequence="ACGTACGTACGTACGTACGG"[:20], pam="AGG", strand="+", start=0, end=23, gc_content=50.0, score=10.0)
        checks = validate_guide(guide)
        assert checks["score_above_threshold"] is False
        assert checks["overall_pass"] is False


# ======================== NEW OFF-TARGET TESTS ========================
class TestOffTarget:
    def test_offtarget_score_calculation(self):
        spacer = "GCTAGCTAGCTAGCTAGCTA"
        pam = "AGG"
        genome_with_many = (spacer + pam) * 10 + "X" * 50
        count, score, _ = calculate_offtarget_score(spacer, pam, genome_with_many)
        # One exact match is treated as the intended target and excluded.
        assert count == 9
        assert 0 <= score <= 100

    def test_offtarget_no_matches(self):
        spacer = "GCTAGCTAGCTAGCTAGCTA"
        pam = "AGG"
        genome_short = "X" * 50
        count, score, _ = calculate_offtarget_score(spacer, pam, genome_short)
        assert count == 0
        assert score == 100.0

    def test_offtarget_with_mismatches(self):
        spacer = "GCTAGCTAGCTAGCTAGCTA"
        pam = "AGG"
        full_target = spacer + pam
        # Flip 2 bases in the spacer portion -> within the default
        # max_mismatches=3 tolerance, so these repeats should still be
        # flagged as near-match off-targets.
        near_match = full_target[:9] + "A" + full_target[10:14] + "A" + full_target[15:]
        assert sum(a != b for a, b in zip(near_match, full_target)) == 2
        genome_mixed = near_match * 4 + "X" * 50
        count, score, notes = calculate_offtarget_score(spacer, pam, genome_mixed)
        assert count == 4
        assert 0.0 < score < 100.0
        assert all("2 mismatch" in n for n in notes)

    def test_offtarget_empty_genome(self):
        spacer = "GCTAGCTAGCTAGCTAGCTA"
        pam = "AGG"
        count, score, notes = calculate_offtarget_score(spacer, pam, None)
        assert count == 0
        assert score == 0.0
        assert "Genome not provided" in notes

    def test_exact_intended_target_is_excluded(self):
        spacer = "GCTAGCTAGCTAGCTAGCTA"
        report = analyze_offtargets(spacer, spacer + "AGG")
        assert report.on_target_excluded is True
        assert report.hits == ()
        assert report.specificity_score == 100.0

    def test_additional_exact_copy_is_critical(self):
        spacer = "GCTAGCTAGCTAGCTAGCTA"
        report = analyze_offtargets(spacer, (spacer + "AGG") * 2)
        assert len(report.hits) == 1
        assert report.hits[0].risk == "Critical"
        assert report.specificity_score < 50.0

    def test_design_guides_attaches_screening_results(self):
        guides = design_guides(
            SAMPLE_SEQUENCE,
            min_score=0,
            max_guides=2,
            genome_context=SAMPLE_SEQUENCE,
        )
        assert guides
        assert all(guide.specificity_score is not None for guide in guides)
        assert all(guide.off_target_count is not None for guide in guides)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))

