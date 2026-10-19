"""Checks of the warnings and quality flags of MANC-Q 1.9.2 on synthetic data only (seconds, no fitting).
Run with:  python -m pytest tests   (or:  python tests/test_qc.py)"""
import os
import sys
import tempfile
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mancq import engine as q
from mancq import demo as dm


def tsp_spectrum(phase_deg=0.0, lo=-1.0, hi=12.0, n=60000, seed=0):
    sf = 600.13
    x = np.linspace(lo, hi, n)
    y = 1e9 * q.pv(x, 0.0, 1.0 / sf, 0.7)
    if phase_deg:
        yi = dm.imaginary_part(y[::-1])[::-1]
        y = y * np.cos(np.deg2rad(phase_deg)) + yi * np.sin(np.deg2rad(phase_deg))
    y += np.random.default_rng(seed).normal(0, 1, n) * y.max() / 5000
    return x, y, sf


def test_tsp_phase_error_is_measured():
    for ph in (0.0, 1.5, -3.0):
        x, y, sf = tsp_spectrum(ph)
        t = q.fit_tsp(x, y, sf)
        assert abs(t["phase_deg"] - ph) < 0.4, (ph, t["phase_deg"])


def test_no_tsp_message():
    x, y, sf = tsp_spectrum()
    y = np.random.default_rng(1).normal(0, 1, len(x))
    try:
        q.fit_tsp(x, y, sf)
        raise AssertionError("no error without a TSP peak")
    except ValueError as exc:
        assert "near 0 ppm" in str(exc)


def test_noise_without_free_region():
    x, y, sf = tsp_spectrum(lo=-0.4, hi=9.5)        # narrow spectral width: no signal-free ends
    why = []
    s = q.noise_sigma(x, y, why)
    assert why, "no warning for the fallback noise estimate"
    true = y.max() / 5000
    assert 0.7 < s / true < 1.3, s / true


def test_dilutions_file_and_warnings():
    with tempfile.TemporaryDirectory() as tmp:
        f = os.path.join(tmp, "dil.csv")
        open(f, "w").write("Sample;Dilution\nEXPNO 10;2\n11.0;3\n99;4\n")
        msgs = []
        d = q.read_dilutions(f, ["10", "11", "12"], msgs.append)
        assert d["10"] == 2 and d["11"] == 3
        text = " ".join(msgs)
        assert "no dilution for EXPNO 12" in text and "99" in text


def test_plain_messages_when_reading():
    with tempfile.TemporaryDirectory() as tmp:
        e = os.path.join(tmp, "20"); os.makedirs(os.path.join(e, "pdata", "1"))
        for what, setup in (("no 'procs' file", lambda: None),
                            ("choose that processing number", lambda: None),
                            ("2D experiment", lambda: open(os.path.join(e, "acqu2s"), "w").close())):
            setup()
            procno = 2 if "processing number" in what else 1
            try:
                q.read_spectrum(e, procno)
                raise AssertionError(f"no error ({what})")
            except ValueError as exc:
                assert what in str(exc), (what, str(exc))


def test_qc_flags():
    Q = pd.DataFrame([
        dict(sample="1", TSP_fwhm_Hz=1.0, metabolite_fwhm_Hz=1.0, chi2_per_point=1.2,
             signal_explained_by_library=0.97, unassigned_singlets=2, TSP_phase_deg=0.1, D1_s=8.0, AQ_s=2.7),
        dict(sample="2", TSP_fwhm_Hz=1.9, metabolite_fwhm_Hz=1.0, chi2_per_point=1.2,
             signal_explained_by_library=0.97, unassigned_singlets=2, TSP_phase_deg=1.5, D1_s=8.0, AQ_s=2.7),
        dict(sample="3", TSP_fwhm_Hz=1.0, metabolite_fwhm_Hz=1.0, chi2_per_point=12.0,
             signal_explained_by_library=0.6, unassigned_singlets=50, TSP_phase_deg=0.1, D1_s=1.0, AQ_s=2.0)])
    f = q.qc_flags(Q).set_index("sample")
    assert f.at["1", "qc_flag"] == "ok", f.at["1", "qc_reason"]
    assert f.at["2", "qc_flag"] == "check"
    r2 = f.at["2", "qc_reason"]
    assert "phase" in r2 and "broader than the metabolites" in r2
    assert f.at["3", "qc_flag"] == "poor" and "D1 + AQ" in f.at["3", "qc_reason"]


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
    print("ok")
