"""The reported upper bound has a safety margin (MANC-Q 2.0.10). In the stress test on synthetic spectra the bound
taken from the data alone was below the true amount for 23 to 26 % of present compounds graded Upper bound only.
upper_bound_mM is now 1.8 x the data bound + 1 x LOD (setting UB_MARGIN); the grading decisions still use the data
bound, so no value and no tier changes. Synthetic demo spectrum only, narrow ppm range (about a minute).
Run with:  python -m pytest tests   (or:  python tests/test_upper_bound_margin.py)"""
import os
import sys
import tempfile
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mancq import engine as q
from mancq import demo as dm
from mancq import regrade as rg
from mancq.grade import reported_bound

FAST = dict(FIT_RANGE=(0.80, 1.60), OVERLAYS=False, N_WORKERS=1)


def test_reported_bound_rule():
    S = q.Settings()
    a, b = S["UB_MARGIN"]
    assert (a, b) == (1.8, 1.0)
    assert np.isclose(reported_bound(0.5, 0.02, S=S), a * 0.5 + b * 0.02)
    assert reported_bound(0.5, 0.02, S=S) > 0.5                                     # never below the data bound
    assert np.isclose(reported_bound(0.5, np.nan, S=S), a * 0.5)                  # no LOD: the factor only
    assert np.isnan(reported_bound(np.nan, 0.02, S=S))                            # no bound stays no bound
    assert reported_bound(0.5, 0.02, S=q.Settings(UB_MARGIN=(1.0, 0.0))) == 0.5   # (1, 0) switches it off


def test_bound_margin_in_results():
    with tempfile.TemporaryDirectory() as tmp:
        spec, res = os.path.join(tmp, "spectra"), os.path.join(tmp, "results")
        ppm0, sw_hz, y = dm.make_spectrum(cache_dir=os.path.join(tmp, "cache"))
        dm.write_bruker(os.path.join(spec, "10"), ppm0, sw_hz, 600.13, y, "MANC-Q synthetic demo (processed)")
        q.run(spec, res, dm.REF_MM, **FAST)
        g = pd.read_csv(os.path.join(res, "concentrations_long.csv"), dtype={"sample": str})
        assert (g.tier == "Upper bound only").any(), g.tier.value_counts()
        # regrading with the margin switched off changes the bound only: same tiers, values and SE
        rg.regrade(res, UB_MARGIN=(1.0, 0.0))
        g0 = pd.read_csv(os.path.join(res, "concentrations_long.csv"), dtype={"sample": str})
        assert list(g0.tier) == list(g.tier)
        assert np.allclose(g0.conc_mM.fillna(-1), g.conc_mM.fillna(-1))
        assert np.allclose(g0.se_mM.fillna(-1), g.se_mM.fillna(-1))
        ok = np.isfinite(g0.upper_bound_mM)
        assert ok.any()
        assert np.allclose(g.upper_bound_mM[ok], 1.8 * g0.upper_bound_mM[ok] + g0.LOD_mM[ok].fillna(0), rtol=1e-6)
        wb = pd.read_excel(os.path.join(res, "metabolite_concentrations.xlsx"), "README")
        assert wb.text.str.contains("x LOD, a margin").any()


if __name__ == "__main__":
    test_reported_bound_rule()
    test_bound_margin_in_results()
    print("ok")
