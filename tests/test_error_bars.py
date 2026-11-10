"""The reported standard error includes model error, not only the noise (MANC-Q 2.0.9). The stress test on
synthetic spectra showed the noise-only SE held the truth for 59 % of Quantified values instead of 95 %.
se_mM is now the noise SE with (a + b x best_misfit) x value added in quadrature; the noise SE stays in
se_noise_mM, and no value or tier changes. Synthetic demo spectrum only, narrow ppm range (about a minute).
Run with:  python -m pytest tests   (or:  python tests/test_error_bars.py)"""
import os
import sys
import tempfile
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mancq import engine as q
from mancq import demo as dm
from mancq import regrade as rg
from mancq.grade import total_se

FAST = dict(FIT_RANGE=(0.80, 1.60), OVERLAYS=False, N_WORKERS=1)
Q, O = "Quantified", "Overlapped (semi-quantitative)"


def test_total_se_rule():
    S = q.Settings()
    a, b = S["SE_MODEL_QUANTIFIED"]
    assert np.isclose(total_se(Q, 2.0, 0.01, 0.04, S=S), np.hypot(0.01, (a + b * 0.04) * 2.0))
    assert total_se(Q, 2.0, 0.01, 0.2, S=S) > total_se(Q, 2.0, 0.01, 0.02, S=S)          # worse fit, wider error bar
    assert total_se(O, 2.0, 0.01, 0.04, S=S) > total_se(Q, 2.0, 0.01, 0.04, S=S)          # Overlapped wider than Quantified
    assert total_se(Q, 2.0, 0.01, np.nan, S=S) == np.hypot(0.01, a * 2.0)                # no misfit: the constant part only
    assert total_se(Q, 1.0, 0.01, 50.0, S=S) <= np.hypot(0.01, 1.0) + 1e-12              # the term never exceeds the value
    assert np.isnan(total_se("Upper bound only", np.nan, np.nan, 0.1, S=S))              # a bound has no SE
    assert total_se("Not detected", 0.0, 0.003, 0.1, S=S) == 0.003                        # no value: unchanged
    off = q.Settings(SE_MODEL_QUANTIFIED=(0, 0))
    assert total_se(Q, 2.0, 0.01, 0.04, S=off) == 0.01                                    # (0, 0) switches it off


def test_error_bars_in_results():
    with tempfile.TemporaryDirectory() as tmp:
        spec, res = os.path.join(tmp, "spectra"), os.path.join(tmp, "results")
        ppm0, sw_hz, y = dm.make_spectrum(cache_dir=os.path.join(tmp, "cache"))
        dm.write_bruker(os.path.join(spec, "10"), ppm0, sw_hz, 600.13, y, "MANC-Q synthetic demo (processed)")
        q.run(spec, res, dm.REF_MM, **FAST)
        g = pd.read_csv(os.path.join(res, "concentrations_long.csv"), dtype={"sample": str})
        v = g[g.tier.isin([Q, O])]
        assert len(v) >= 3, g.tier.value_counts()
        assert (v.se_mM > v.se_noise_mM).all()
        S = q.Settings()
        for r in v.itertuples():
            a, b = S["SE_MODEL_QUANTIFIED"] if r.tier == Q else S["SE_MODEL_OVERLAPPED"]
            lo, hi = (np.hypot(r.se_noise_mM, min(a + b * max(r.best_misfit + d, 0.0), 1.0) * r.conc_mM)
                      for d in (-0.0006, 0.0006))                                        # best_misfit is saved to 3 decimals
            assert lo * 0.999 <= r.se_mM <= hi * 1.001, (r.metabolite, lo, r.se_mM, hi)
        ub = g[g.tier == "Upper bound only"]
        assert ub.se_mM.isna().all() and ub.se_noise_mM.isna().all()
        # the truth is inside 1.96 SE for the Quantified values of the demo
        truth = pd.Series(dm.TRUTH_MM)
        qv = g[(g.tier == Q) & g.metabolite.isin(truth.index)]
        inside = (qv.conc_mM - qv.metabolite.map(truth)).abs() <= 1.96 * qv.se_mM
        assert inside.mean() >= 0.8, qv.assign(truth=qv.metabolite.map(truth))[["metabolite", "conc_mM", "se_mM", "truth"]]
        wb = pd.read_excel(os.path.join(res, "metabolite_concentrations.xlsx"), "README")
        assert wb.text.str.contains("se_noise_mM").any()
        # regrading with the term switched off changes the SE only: same tiers, same values, SE = noise SE
        rg.regrade(res, SE_MODEL_QUANTIFIED=(0, 0), SE_MODEL_OVERLAPPED=(0, 0), SE_MODEL_DECONV=(0, 0))
        g0 = pd.read_csv(os.path.join(res, "concentrations_long.csv"), dtype={"sample": str})
        assert list(g0.tier) == list(g.tier)
        assert np.allclose(g0.conc_mM.fillna(-1), g.conc_mM.fillna(-1))
        assert np.allclose(g0.se_mM.fillna(-1), g.se_noise_mM.fillna(-1))


if __name__ == "__main__":
    test_total_se_rule()
    test_error_bars_in_results()
    print("ok")
