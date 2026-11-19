"""A compound is not left out of the fit because a neighbour inflated one of its multiplets (MANC-Q 2.0.12). Up to
2.0.11 a compound was left out as a "ghost" if any one of its well fitted multiplets, taken alone, implied more of it
than its other multiplets show. On real 800 MHz culture spectra this dropped glucose at 20-50 mM: its multiplet at
3.54 ppm shares a peak with other compounds and implied four times more glucose than the rest of the molecule, so
the alpha-H1 doublet at 5.24 ppm looked too small. The amount is now judged from the multiplet that implies the least.
The rule is checked with the numbers from that spectrum; the whole-spectrum check (synthetic: the demo mixture with
a foreign singlet, not in the library, on glucose's 3.54 ppm multiplet) only shows that glucose keeps its value.
Run with:  python -m pytest tests   (or:  python tests/test_ghost_least_amount.py)"""
import os
import sys
import tempfile
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mancq import engine as q
from mancq import demo as dm
from mancq.fit import _ghost_test_multiplet

FAST = dict(OVERLAYS=False, N_WORKERS=1)
TRUTH = dm.TRUTH_MM


def glucose(tmp, **settings):
    spec, res = os.path.join(tmp, "spectra"), os.path.join(tmp, "results_" + str(len(settings)))
    if not os.path.isdir(spec):
        ppm0, sw_hz, y = dm.make_spectrum(cache_dir=os.path.join(tmp, "cache"), truth=TRUTH,
                                          singlets=[(3.5425, 4.0)])           # a foreign peak on glucose's alpha-H2
        dm.write_bruker(os.path.join(spec, "10"), ppm0, sw_hz, 600.13, y, "MANC-Q synthetic test: inflated multiplet")
    q.run(spec, res, dm.REF_MM, **FAST, **settings)
    g = pd.read_csv(os.path.join(res, "concentrations_long.csv"), dtype={"sample": str}).set_index("metabolite")
    return g.loc["D-Glucose"]


def test_least_amount_rule():
    # glucose in a real 800 MHz culture spectrum (sample 49): multiplet 6 at 3.54 ppm shares its peak and implies
    # four times more glucose than the others, including the alpha-H1 doublet at 5.24 ppm that glucose owns
    S = q.Settings()
    amps = {1: 2.066e4, 2: 2.362e4, 4: 2.068e4, 5: 2.387e4, 6: 9.857e4, 13: 2.363e4}
    doms = {1: 0.65, 2: 0.61, 4: 0.67, 5: 0.63, 6: 0.69, 13: 1.00}
    carry = [dict(j=j - 1, local_amp=a, loc_dom=doms[j]) for j, a in amps.items()]
    carry.append(dict(j=7, local_amp=5.0e3, loc_dom=0.23))                  # not dominant: not used
    assert _ghost_test_multiplet(carry, True, S=S)["j"] == 0                # multiplet 1, the least of the dominant ones
    assert _ghost_test_multiplet(carry, False, S=S)["j"] == 7               # no dominant multiplet: the least of all
    # a look-alike with one dominant multiplet (creatinine on the creatine CH3 singlet) is still judged from it
    one = [dict(j=0, local_amp=5.0e4, loc_dom=0.9), dict(j=1, local_amp=1.0e3, loc_dom=0.1)]
    assert _ghost_test_multiplet(one, True, S=S)["j"] == 0
    assert _ghost_test_multiplet(carry, True, S=q.Settings(GHOST_LEAST_AMOUNT=False)) is None   # switched off


def test_inflated_multiplet_does_not_drop_glucose():
    with tempfile.TemporaryDirectory() as tmp:
        new = glucose(tmp)
        assert not str(new.method).startswith("left out of the fit"), new.method
        assert new.tier in ("Quantified", "Overlapped (semi-quantitative)"), (new.tier, new.method)
        assert abs(new.conc_mM / TRUTH["D-Glucose"] - 1) < 0.10, new.conc_mM
        old = glucose(tmp, GHOST_LEAST_AMOUNT=False)                              # the rule as up to 2.0.11
        print("up to 2.0.11:", old.tier, "|", old.method)


if __name__ == "__main__":
    test_least_amount_rule()
    test_inflated_multiplet_does_not_drop_glucose()
    print("ok")
