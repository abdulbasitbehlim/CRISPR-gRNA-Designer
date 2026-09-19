"""CFD weights and optional Rule Set 2 interface with explicit provider status."""
from functools import lru_cache
import json
import math
from pathlib import Path
import warnings


@lru_cache(maxsize=1)
def cfd_weights():
    return json.loads((Path(__file__).parent / 'data' / 'cfd_data.json').read_text())


def cfd_score(spacer, off_target, pam='NGG'):
    """Doench CFD pair score; spacer and off-target are in guide orientation."""
    if len(spacer) != 20 or len(off_target) != 20 or set(spacer + off_target) - set('ACGT'):
        raise ValueError('CFD requires two unambiguous 20 nt DNA sequences.')
    if len(pam) != 3 or pam[1:] not in cfd_weights()['pam']:
        raise ValueError('Unsupported CFD PAM.')
    weight = float(cfd_weights()['pam'][pam[1:]])
    complement = dict(A='T', C='G', G='C', T='A')
    for position, (guide, target) in enumerate(zip(spacer, off_target), 1):
        if guide != target:
            key = f"r{guide.replace('T', 'U')}:d{complement[target]},{position}"
            weight *= cfd_weights()['mm'][key]
    return weight


@lru_cache(maxsize=1)
def rs2_provider():
    try:
        from guidemaker.doench_predict import predict
        return predict
    except ImportError:
        return None
    except Exception as exc:
        warnings.warn(f'Rule Set 2 provider failed to load: {exc}', RuntimeWarning)
        return None


def rs2_status():
    return 'GuideMaker Rule Set 2 available' if rs2_provider() else 'Rule Set 2 unavailable; sequence heuristic only'


def rs2_score(context):
    if not context or len(context) != 30 or set(context) - set('ACGT') or context[25:27] != 'GG':
        return None
    provider = rs2_provider()
    if provider is None:
        return None
    import numpy as np
    try:
        value = float(np.asarray(provider(np.array([context]), num_threads=1)).reshape(-1)[0])
        if not math.isfinite(value):
            raise ValueError('Provider returned a nonfinite score')
        # Report regression output directly, on the existing 100x display scale.
        # Do not clip a regression prediction or call it an editing probability.
        return round(value * 100, 2)
    except Exception as exc:
        warnings.warn(f'Rule Set 2 prediction failed: {exc}', RuntimeWarning)
        return None
