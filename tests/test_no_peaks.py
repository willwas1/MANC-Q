"""A spectrum with no peak tall enough to list (very dilute or very noisy) stays in the results when a run is
carried on, regraded or adjusted by hand (MANC-Q 2.0.7). Up to 2.0.6 its table of observed peaks was saved as an
empty file that could not be read back, so the spectrum was reported as failed and left out of every table.
Synthetic spectra only: the demo spectrum and one with the reference line and noise alone.
Run with:  python -m pytest tests   (or:  python tests/test_no_peaks.py)"""
import os
import sys
import tempfile
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mancq import engine as q
from mancq import demo as dm
from mancq import regrade as rg

FAST = dict(FIT_RANGE=(0.80, 1.60), OVERLAYS=False, N_WORKERS=2)


def samples_in(res):
    g = pd.read_csv(os.path.join(res, "concentrations_long.csv"), dtype={"sample": str})
    return sorted(g["sample"].unique())


def test_no_peaks():
    with tempfile.TemporaryDirectory() as tmp:
        spec, res, cache = os.path.join(tmp, "spectra"), os.path.join(tmp, "results"), os.path.join(tmp, "cache")
        ppm0, sw_hz, y = dm.make_spectrum(cache_dir=cache)
        dm.write_bruker(os.path.join(spec, "10"), ppm0, sw_hz, 600.13, y, "MANC-Q synthetic demo (processed)")
        ppm0, sw_hz, y = dm.make_spectrum(seed=7, cache_dir=cache, truth={})          # reference line and noise only
        dm.write_bruker(os.path.join(spec, "20"), ppm0, sw_hz, 600.13, y, "MANC-Q synthetic test: no peaks")
        q.run(spec, res, dm.REF_MM, **FAST)
        assert samples_in(res) == ["10", "20"]
        obs = os.path.join(res, "_per_sample", "20_obs.csv")
        saved = pd.read_csv(obs)                                # readable: the column headings are there
        assert len(saved) == 0 and list(saved.columns) == q.OBS_COLUMNS, saved
        g = pd.read_csv(os.path.join(res, "concentrations_long.csv"), dtype={"sample": str})
        assert not (g[g["sample"] == "20"].tier == "Quantified").any()

        # a results folder from 2.0.6 or earlier has an empty file there: it must read as "no peaks"
        with open(obs, "w", encoding="utf-8") as fh:
            fh.write("\n")
        assert len(q.read_observed(obs)) == 0
        # carrying on the run re-uses both spectra and keeps both in the tables
        q.run(spec, res, dm.REF_MM, existing="resume", **FAST)
        assert samples_in(res) == ["10", "20"]
        with open(os.path.join(res, "run_log.txt"), encoding="utf-8") as fh:
            log = fh.read()
        assert "Sample 20: re-used earlier result" in log and "FAILED" not in log, log[-2000:]
        # and so does regrading
        rg.regrade(res)
        assert samples_in(res) == ["10", "20"]
        peaks = pd.read_csv(os.path.join(res, "observed_peaks_assigned.csv"), dtype={"sample": str})
        assert set(peaks["sample"]) == {"10"}


if __name__ == "__main__":
    test_no_peaks()
    print("ok")
