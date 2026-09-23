# ============================================================================
# REPORTING
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Turns analysis results into clear tables or downloadable report content.
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
# - function: sequence_hash
# - function: make_metadata
# - function: _finite
# - function: export_json
# ============================================================================

"""Strict JSON export and reproducible sequence/reference fingerprints."""
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import math



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: sequence_hash
# ----------------------------------------------------------------------------
def sequence_hash(sequence):
    return hashlib.sha256(sequence.encode('ascii')).hexdigest()



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: make_metadata
# ----------------------------------------------------------------------------
def make_metadata(sequence, reference, **settings):
    return {
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'target_length': len(sequence), 'target_sha256': sequence_hash(sequence),
        'reference_records': [{'id': k, 'length': len(v), 'sha256': sequence_hash(v)} for k, v in (reference or {}).items()],
        'coordinate_system': 'Tables use 1-based inclusive intervals; Start/End cover spacer plus PAM; Cut after base is a between-base position. Internal intervals are 0-based half-open.',
        'local_search_scope': 'Supplied linear reference only; NGG PAM; substitutions only; ambiguous sites skipped; no variants, bulges or circular junctions.',
        'ranking': 'Sequence heuristic with local specificity as tie-breaker; CRISPRi uses TSS bands first. RS2 is not used to rank CRISPRi.',
        **settings,
    }



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: _finite
# ----------------------------------------------------------------------------
def _finite(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(k): _finite(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_finite(v) for v in value]
    return value



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: export_json
# ----------------------------------------------------------------------------
def export_json(metadata, guides, rows=None, genome_rows=None):
    data = {
        'metadata': metadata,
        'guides': rows if rows is not None else [g.to_dict() for g in guides],
        'local_hit_details': [{
            'guide_index': i, 'screened_mismatch_radius': g.screened_mismatch_radius,
            'total_hits': g.off_target_count,
            'risk_counts': None if g.screened_mismatch_radius is None else g.risk_counts,
            'maximum_per_site_mit_risk': g.max_mit_risk,
            'maximum_per_site_cfd_risk': g.max_cfd_risk,
            'displayed_hits': len(g.off_target_details),
            'details_truncated': (g.off_target_count or 0) > len(g.off_target_details),
            'hits': g.off_target_details,
        } for i, g in enumerate(guides, 1)],
        'guidescan_results': genome_rows or [],
    }
    return json.dumps(_finite(data), indent=2, allow_nan=False).encode('utf-8')
