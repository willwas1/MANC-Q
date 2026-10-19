"""Rerunning only some spectra into a results folder never leaves a mix of old and new results (MANC-Q 1.8.1),
on the synthetic demo spectra only. Fits a narrow ppm range so it takes about a minute.
Run with:  python -m pytest tests   (or:  python tests/test_subset_rerun.py)"""
import os
import sys
import tempfile
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mancq import engine as q
from mancq import demo as dm

FAST = dict(FIT_RANGE=(0.80, 1.60), OVERLAYS=False, N_WORKERS=2)


def samples_in(res):
    g = pd.read_csv(os.path.join(res, "concentrations_long.csv"), dtype={"sample": str})
    return sorted(g["sample"].unique())


def test_subset_rerun():
    with tempfile.TemporaryDirectory() as tmp:
        spec, res = os.path.join(tmp, "spectra"), os.path.join(tmp, "results")
        dm.write_demo(spec, cache_dir=os.path.join(tmp, "cache"), extra=True)
        q.run(spec, res, dm.REF_MM, SAMPLES=["10", "13"], **FAST)
        assert samples_in(res) == ["10", "13"]
        # into the same folder without saying what to do: refused
        try:
            q.run(spec, res, dm.REF_MM, SAMPLES=["10"], **FAST)
            raise AssertionError("a run into a folder with results was not refused")
        except q.ResultsExist:
            pass
        # carrying on a run with other spectra than it had: refused
        try:
            q.run(spec, res, dm.REF_MM, SAMPLES=["10"], existing="resume", **FAST)
            raise AssertionError("carrying on with other spectra was not refused")
        except ValueError:
            pass
        assert samples_in(res) == ["10", "13"]
        # overwrite with a subset: nothing of EXPNO 13 is left anywhere
        q.run(spec, res, dm.REF_MM, SAMPLES=["10"], existing="overwrite", **FAST)
        assert samples_in(res) == ["10"]
        left = [f for d, _, fs in os.walk(res) for f in fs if f.startswith("13_") or "sample_13" in f]
        assert not left, left
        wb = pd.read_excel(os.path.join(res, "metabolite_concentrations.xlsx"), "QC", dtype={"sample": str})
        assert list(wb["sample"]) == ["10"]
        # a new subfolder keeps the results beside it as they are
        sub = q.new_run_folder(res)
        q.run(spec, sub, dm.REF_MM, SAMPLES=["13"], **FAST)
        assert samples_in(sub) == ["13"] and samples_in(res) == ["10"]


if __name__ == "__main__":
    test_subset_rerun()
    print("ok")
