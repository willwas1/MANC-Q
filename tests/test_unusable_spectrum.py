"""Two failures found by the stress test on synthetic spectra (MANC-Q 2.0.8):
1. baseline correction with closely spaced knots gave a runaway baseline (values near 10^98) where a stretch of
   signal, or the water gap, was wider than the knot spacing;
2. a spectrum with nothing usable at 0 ppm stopped with "Each lower bound must be strictly less than each upper
   bound" instead of a message a user can act on.
Synthetic demo spectrum only.
Run with:  python -m pytest tests   (or:  python tests/test_unusable_spectrum.py)"""
import os
import sys
import tempfile
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mancq import engine as q
from mancq import demo as dm
from mancq import fidproc as fp

SF = 600.13


def demo_spectrum(tmp):
    """The synthetic demo spectrum: ppm (descending, as in a Bruker file) and intensities."""
    ppm0, sw_hz, y = dm.make_spectrum(SF, cache_dir=os.path.join(tmp, "cache"))
    ppm = ppm0 - np.arange(y.size) * sw_hz / SF / y.size
    return ppm, np.asarray(y, float)


def test_baseline_correction_stays_bounded():
    with tempfile.TemporaryDirectory() as tmp:
        ppm, y = demo_spectrum(tmp)
        keep = np.abs(ppm - 4.75) > 0.25
        tsp = np.abs(ppm) < 0.01
        for knot_ppm in (0.5, 0.25, 0.1, 0.05):
            out = fp.baseline_correct(ppm, y, keep, knot_ppm=knot_ppm)
            assert np.isfinite(out).all(), knot_ppm
            assert np.abs(out).max() <= 2 * np.abs(y).max(), (knot_ppm, np.abs(out).max() / np.abs(y).max())
            assert 0.8 < out[tsp].max() / y[tsp].max() < 1.05, (knot_ppm, out[tsp].max() / y[tsp].max())


def test_unusable_spectrum_message():
    with tempfile.TemporaryDirectory() as tmp:
        ppm, y = demo_spectrum(tmp)
        S = q.Settings()
        nothing = y.copy()
        nothing[np.abs(ppm) < 3.0] = 0.0                       # like the damaged stress spectra: zeros around 0 ppm
        nothing[ppm > 9.5] = 0.0
        for bad, words in ((nothing, "peak near 0 ppm"), (np.where(np.abs(ppm - 2) < 0.1, np.nan, y), "not numbers")):
            try:
                q.fit_tsp(ppm[::-1], bad[::-1], SF, S=S)
                raise AssertionError("an unusable spectrum was accepted")
            except ValueError as exc:
                assert words in str(exc), str(exc)
        q.fit_tsp(ppm[::-1], y[::-1], SF, S=S)                 # the good spectrum is still accepted


if __name__ == "__main__":
    test_baseline_correction_stays_bounded()
    test_unusable_spectrum_message()
    print("ok")
