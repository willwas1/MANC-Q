"""A value that rests on one singlet needs stronger evidence (MANC-Q 2.0.11). In the stress test on synthetic
spectra, compounds whose only evidence is one singlet (glycolate, succinate, oxaloacetate, acetone, ...) gave 70 % of
the values reported for compounds that were not there, mostly from another singlet near their position. Such a value
is now given only if the singlet has at least half of the fitted signal at its peak and sits within 5 Hz of its
library position (or clearly owns its peak, dominance 0.9 or more); otherwise the compound is Upper bound only.
Synthetic demo spectrum only (about a minute).
Run with:  python -m pytest tests   (or:  python tests/test_lone_singlet.py)"""
import os
import sys
import tempfile
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mancq import engine as q
from mancq import demo as dm
from mancq import regrade as rg
from mancq.grade import _lone_singlet_doubt

FAST = dict(FIT_RANGE=(1.80, 2.20), OVERLAYS=False, N_WORKERS=1)
OFF = dict(LONE_SINGLET_MIN_DOMINANCE=0.0, LONE_SINGLET_MAX_SHIFT_HZ=1e9, LONE_SINGLET_SHIFT_OK_DOMINANCE=0.0)


def singlet(dom, shift_ppm, lines=1):
    it = np.ones(lines)
    return dict(p=np.linspace(0, 0.01, lines), it=it, d=1.0 + shift_ppm, base=1.0, loc_dom=dom)


def test_lone_singlet_rule():
    S, sf = q.Settings(), 600.0
    assert _lone_singlet_doubt([singlet(0.95, 0.0)], sf, S=S) == ""                 # clean and in place
    assert _lone_singlet_doubt([singlet(0.6, 3 / sf)], sf, S=S) == ""                # 3 Hz away: fine
    assert "only 40%" in _lone_singlet_doubt([singlet(0.4, 0.0)], sf, S=S)           # shared peak
    assert "Hz from its library" in _lone_singlet_doubt([singlet(0.7, 8 / sf)], sf, S=S)   # moved and shared
    assert _lone_singlet_doubt([singlet(0.95, 8 / sf)], sf, S=S) == ""               # moved but clearly its own
    assert _lone_singlet_doubt([singlet(0.3, 0.0, lines=3)], sf, S=S) == ""          # not a singlet: rule not used
    assert _lone_singlet_doubt([singlet(0.3, 0.0), singlet(0.3, 0.0)], sf, S=S) == ""     # two reporters
    assert _lone_singlet_doubt([singlet(0.3, 0.0)], sf, S=q.Settings(**OFF)) == ""        # switched off


def test_lone_singlet_in_results():
    with tempfile.TemporaryDirectory() as tmp:
        spec, res = os.path.join(tmp, "spectra"), os.path.join(tmp, "results")
        ppm0, sw_hz, y = dm.make_spectrum(cache_dir=os.path.join(tmp, "cache"))
        dm.write_bruker(os.path.join(spec, "10"), ppm0, sw_hz, 600.13, y, "MANC-Q synthetic demo (processed)")
        q.run(spec, res, dm.REF_MM, **FAST)
        g = pd.read_csv(os.path.join(res, "concentrations_long.csv"), dtype={"sample": str}).set_index("metabolite")
        rg.regrade(res, **OFF)
        g0 = pd.read_csv(os.path.join(res, "concentrations_long.csv"), dtype={"sample": str}).set_index("metabolite")
        # the rule only ever takes a value away, and says why
        changed = g.index[g.tier != g0.tier]
        assert (g.loc[changed, "tier"] == "Upper bound only").all(), g.loc[changed, ["tier"]].join(g0.tier, rsuffix="_off")
        assert g.loc[changed, "method"].str.contains("rests on one singlet").all()
        same = g.index[g.tier == g0.tier]
        assert np.allclose(g.loc[same, "conc_mM"].fillna(-1), g0.loc[same, "conc_mM"].fillna(-1))
        # acetate is in the demo and its singlet is its own: it keeps its value
        assert g.at["Acetate", "tier"] == "Quantified", g.loc["Acetate"]
        assert abs(g.at["Acetate", "conc_mM"] / dm.TRUTH_MM["Acetate"] - 1) < 0.10


if __name__ == "__main__":
    test_lone_singlet_rule()
    test_lone_singlet_in_results()
    print("ok")
