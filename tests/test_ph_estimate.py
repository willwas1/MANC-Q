"""Checks of the histidine pH estimate on synthetic spectra only.
Run with:  python -m pytest tests   (or:  python tests/test_ph_estimate.py)"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mancq import engine as q

SF = 600.13


def lorentz(x, x0, fwhm_hz, height):
    g = fwhm_hz / SF / 2
    return height * g ** 2 / ((x - x0) ** 2 + g ** 2)


def spectrum(ph, extra=(), seed=0, height=500.0):
    """TSP-referenced ppm axis (ascending, 0.25 Hz per point) with His H2 and H4 singlets at the shifts predicted
    for this pH, plus optional extra peaks (ppm, height); noise sigma 1."""
    rng = np.random.default_rng(seed)
    x = np.arange(6.5, 8.6, 0.25 / SF)
    y = rng.normal(0, 1, x.size)
    if ph is not None:
        for lim in (q.SETTINGS["PH_HIS_H2_LIMITS"], q.SETTINGS["PH_HIS_H4_LIMITS"]):
            y += lorentz(x, q.shift_from_ph(ph, lim, q.SETTINGS["PH_HIS_PKA"]), 1.2, height)
    for pos, h in extra:
        y += lorentz(x, pos, 1.2, h)
    return x, y


def test_known_ph():
    for ph in (6.6, 7.0, 7.4):
        x, y = spectrum(ph)
        r = q.estimate_ph(x, y, 1.0, S=q.Settings())
        assert r["pH_flag"] == "ok", r
        for k in ("pH_est", "pH_H2", "pH_H4"):
            assert abs(r[k] - ph) < 0.01, (ph, k, r)


def test_round_trip():
    for lim in (q.SETTINGS["PH_HIS_H2_LIMITS"], q.SETTINGS["PH_HIS_H4_LIMITS"]):
        for ph in (6.0, 6.8, 7.5):
            d = q.shift_from_ph(ph, lim, q.SETTINGS["PH_HIS_PKA"])
            assert abs(q.ph_from_shift(d, lim, q.SETTINGS["PH_HIS_PKA"]) - ph) < 1e-9


def test_taller_neighbour_is_not_taken():
    # a taller doublet near 7.74 ppm (as tryptophan) inside the H2 window: the H2 peak that agrees with H4 wins
    x, y = spectrum(7.0, extra=[(7.7405, 2000.0), (7.7270, 2000.0)])
    r = q.estimate_ph(x, y, 1.0, S=q.Settings())
    assert abs(r["pH_H2"] - 7.0) < 0.01 and r["pH_flag"] == "ok", r


def test_weak_and_absent():
    x, y = spectrum(7.0, height=5.0)            # below 10 x noise
    r = q.estimate_ph(x, y, 1.0, S=q.Settings())
    assert r["pH_flag"] == "no_histidine" and np.isnan(r["pH_est"]), r
    x, y = spectrum(None, extra=[(q.shift_from_ph(7.2, q.SETTINGS["PH_HIS_H2_LIMITS"], q.SETTINGS["PH_HIS_PKA"]), 300.0)])
    r = q.estimate_ph(x, y, 1.0, S=q.Settings())
    assert r["pH_flag"] == "one_proton_only" and abs(r["pH_est"] - 7.2) < 0.01, r


def test_placement_targets():
    x, y = spectrum(7.0)
    r = q.estimate_ph(x, y, 1.0, S=q.Settings())
    t = q.histidine_targets(r, S=q.Settings())
    assert abs(t["H2"][0] - r["his_H2_ppm"]) < 1e-9 and t["H2"][1] == q.SETTINGS["PH_HIS_PLACE_HZ"]
    r["his_H4_ppm"] = np.nan                    # not seen: predicted from the pH instead, wider window
    t = q.histidine_targets(r, S=q.Settings())
    assert t["H4"][1] == q.SETTINGS["PH_HIS_PLACE_PREDICTED_HZ"]


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f()
            print("ok ", name)
