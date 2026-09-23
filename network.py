# ============================================================================
# NETWORK
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Contains network-related helper functions used when the program communicates with external sequence services.
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
# - function: get
# ============================================================================

"""Bounded public-data GETs with transient retries and NCBI request pacing."""
import threading
import time
import requests

_ncbi_lock = threading.Lock()
_last_ncbi = 0.0



# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: get
# ----------------------------------------------------------------------------
def get(url, **kwargs):
    global _last_ncbi
    kwargs.setdefault('timeout', 30)
    for attempt in range(3):
        if 'eutils.ncbi.nlm.nih.gov/' in url:
            with _ncbi_lock:
                wait = max(0, 0.35 - (time.monotonic() - _last_ncbi))
                if wait:
                    time.sleep(wait)
                _last_ncbi = time.monotonic()
        try:
            response = requests.get(url, **kwargs)
        except (requests.Timeout, requests.ConnectionError):
            if attempt == 2:
                raise
        else:
            if response.status_code not in {429, 500, 502, 503, 504} or attempt == 2:
                return response
        time.sleep(0.5 * (attempt + 1))
    raise RuntimeError('Unreachable request state')
