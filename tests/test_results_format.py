"""Checks of the settings object, the results-folder format and the library check (MANC-Q 1.10.0), on synthetic
files only (seconds, no fitting).
Run with:  python -m pytest tests   (or:  python tests/test_results_format.py)"""
import json
import os
import pickle
import re
import sys
import tempfile
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mancq import engine as q
from mancq import results as rs

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def folder(tmp, **saved):
    os.makedirs(os.path.join(tmp, "_per_sample"), exist_ok=True)
    with open(os.path.join(tmp, "_per_sample", "run_settings.json"), "w", encoding="utf-8") as fh:
        json.dump(saved, fh)
    pd.DataFrame([dict(sample="1", metabolite="L-Lactate", tier="Quantified", conc_mM=2.0, upper_bound_mM=2.1),
                  dict(sample="1", metabolite="Citrate", tier="Upper bound only", conc_mM=0.5, upper_bound_mM=np.nan),
                  dict(sample="1", metabolite="Acetone", tier="Not detected", conc_mM=0.0, upper_bound_mM=0.01)]
                 ).to_csv(os.path.join(tmp, "concentrations_long.csv"), index=False)
    pd.DataFrame([dict(sample="1", TSP_fwhm_Hz=1.0, chi2_per_point=1.0, signal_explained_by_library=0.97,
                       unassigned_singlets=0)]).to_csv(os.path.join(tmp, "qc_per_sample.csv"), index=False)
    return tmp


def test_settings_object():
    S = q.Settings(LOQ_SNR=8)
    assert S["LOQ_SNR"] == 8 and q.SETTINGS["LOQ_SNR"] == 10          # the defaults are never changed
    for bad in (lambda: q.Settings(LOQ_SNRR=8), lambda: S.update(FIT_RANGEE=(1, 2))):
        try:
            bad()
            raise AssertionError("a misspelt setting was accepted")
        except KeyError:
            pass
    S2 = pickle.loads(pickle.dumps(S))                                  # as sent to the worker processes
    assert isinstance(S2, q.Settings) and S2 == S
    with tempfile.TemporaryDirectory() as tmp:
        os.makedirs(os.path.join(tmp, "_per_sample"))
        S["EXCLUDE"] = [(4.6, 5.0), (5.5, 6.0)]
        S.save(os.path.join(tmp, "_per_sample", "run_settings.json"), RUN_SAMPLES=["1"])
        L = q.Settings.load(tmp)
        assert L["LOQ_SNR"] == 8 and L["EXCLUDE"] == [(4.6, 5.0), (5.5, 6.0)] and L["FIT_RANGE"] == (0.6, 9.6)
        assert L["RESULTS_FORMAT"] == q.RESULTS_FORMAT and L["MANCQ_VERSION"] == q.__version__


def test_version_from_pyproject():
    with open(os.path.join(ROOT, "pyproject.toml"), encoding="utf-8") as fh:
        v = re.search(r'^version\s*=\s*"([^"]+)"', fh.read(), re.M).group(1)
    assert q.__version__ == v


def test_newer_folder_refused():
    with tempfile.TemporaryDirectory() as tmp:
        folder(tmp, RESULTS_FORMAT=q.RESULTS_FORMAT + 1, MANCQ_VERSION="9.0.0")
        for f in (rs.check_readable, rs.load_results):
            try:
                f(tmp)
                raise AssertionError("a folder from a newer MANC-Q was read")
            except rs.ResultsTooNew as exc:
                assert "9.0.0" in str(exc)


def test_other_library_refused():
    with tempfile.TemporaryDirectory() as tmp:
        folder(tmp, RESULTS_FORMAT=q.RESULTS_FORMAT, MANCQ_VERSION="1.10.0", LIBRARY_HASH="0" * 32)
        try:
            rs.check_library(tmp)
            raise AssertionError("a folder fitted with another library was accepted")
        except rs.LibraryChanged:
            pass
        folder(tmp, RESULTS_FORMAT=q.RESULTS_FORMAT, LIBRARY_HASH=q.library_hash(os.path.join(q.HERE, "library")))
        rs.check_library(tmp)


def test_old_folder_brought_up_to_date():
    with tempfile.TemporaryDirectory() as tmp:
        folder(tmp, MANCQ_VERSION="1.8.0")               # format 1; upper bounds were stored in conc_mM
        d = rs.load_results(tmp)
        G = d["long"].set_index("metabolite")
        assert np.isnan(G.at["Citrate", "conc_mM"]) and G.at["Citrate", "upper_bound_mM"] == 0.5
        assert G.at["L-Lactate", "conc_mM"] == 2.0 and G.at["Acetone", "conc_mM"] == 0.0
        assert d["qc"].qc_flag.iloc[0] == "ok" and d["info"]["format"] == 1


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
    print("ok")
