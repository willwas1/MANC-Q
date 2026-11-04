"""Checks of the default number of parallel fits of MANC-Q 2.0.6: one per CPU core but one, and no more than the
free memory allows (seconds, no fitting).
Run with:  python -m pytest tests   (or:  python tests/test_workers.py)"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mancq import settings as st

GB = 2 ** 30


def _workers(cpus, free_gb, n_samples=None):
    keep = os.cpu_count, st.available_memory
    os.cpu_count = lambda: cpus
    st.available_memory = lambda: None if free_gb is None else free_gb * GB
    try:
        return st.default_workers(n_samples)
    finally:
        os.cpu_count, st.available_memory = keep


def test_free_memory_is_read():
    free = st.available_memory()
    assert free is None or free > 0


def test_memory_limits_workers():
    # a 12-core PC with 16 GB, about 10 GB free: 11 fits at once ran out of memory
    assert _workers(12, 10) == int((10 - st.MEMORY_HEADROOM_GB) // st.WORKER_MEMORY_GB) < 11
    assert _workers(12, 64) == 11                  # plenty of memory: cores - 1
    assert _workers(12, None) == 11                # free memory not known: cores - 1, as before
    assert _workers(12, 1) == 1                    # always at least one
    assert _workers(4, 64, n_samples=2) == 2       # never more than the spectra


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
    print("ok")
