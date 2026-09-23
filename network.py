"""Bounded public-data GETs with transient retries and NCBI request pacing."""
import threading
import time
import requests

_ncbi_lock = threading.Lock()
_last_ncbi = 0.0


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
