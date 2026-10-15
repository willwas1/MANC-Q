"""Checks of the chemical shift assignment table (mancq.shifttable) on synthetic results only.
Run with:  python -m pytest tests   (or:  python tests/test_shift_table.py)"""
import os
import sys
import tempfile
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mancq import shifttable as T


def synthetic():
    G = pd.DataFrame([
        dict(sample="1", metabolite="L-Lactate", tier="Quantified", conc_mM=5.0),
        dict(sample="2", metabolite="L-Lactate", tier="Quantified", conc_mM=6.0),
        dict(sample="1", metabolite="L-Threonine", tier="Overlapped (semi-quantitative)", conc_mM=0.5),
        dict(sample="2", metabolite="L-Threonine", tier="Upper bound only", conc_mM=0.6),
        dict(sample="1", metabolite="Citrate", tier="Deconvolution estimate (low confidence)", conc_mM=0.2),
        dict(sample="2", metabolite="Citrate", tier="Not detected", conc_mM=0.0),
        dict(sample="1", metabolite="Acetone", tier="Not detected", conc_mM=0.0),
        dict(sample="2", metabolite="Acetone", tier="Upper bound only", conc_mM=0.1)])
    P = pd.DataFrame([
        # lactate: doublet used in both samples, quartet never used
        dict(sample="1", metabolite="L-Lactate", multiplet=1, centre_ppm=1.3300, height_snr=900, used_as_reporter=True),
        dict(sample="2", metabolite="L-Lactate", multiplet=1, centre_ppm=1.3320, height_snr=950, used_as_reporter=True),
        dict(sample="1", metabolite="L-Lactate", multiplet=2, centre_ppm=4.1100, height_snr=200, used_as_reporter=False),
        dict(sample="2", metabolite="L-Lactate", multiplet=2, centre_ppm=4.1120, height_snr=210, used_as_reporter=False),
        # threonine: one shared reporter; a multiplet in an ignored region (never fitted) is left out
        dict(sample="1", metabolite="L-Threonine", multiplet=1, centre_ppm=1.3250, height_snr=40, used_as_reporter=True),
        dict(sample="1", metabolite="L-Threonine", multiplet=2, centre_ppm=4.2500, height_snr=12, used_as_reporter=False),
        dict(sample="1", metabolite="L-Threonine", multiplet=3, centre_ppm=4.8000, height_snr=0, used_as_reporter=False),
        dict(sample="1", metabolite="Citrate", multiplet=1, centre_ppm=2.6600, height_snr=8, used_as_reporter=False),
        dict(sample="1", metabolite="Acetone", multiplet=1, centre_ppm=2.2300, height_snr=5, used_as_reporter=False)])
    lib = [dict(name="L-Lactate", source="GISSMO", provenance="GISSMO bmse000208/1 (L-lactic acid), QM-simulated"),
           dict(name="L-Threonine", source="GISSMO", provenance="GISSMO bmse000049/1 (L-threonine), QM-simulated"),
           dict(name="Citrate", source="CASMDB", provenance="CASMDB citrate peak list (origin HMDB), spacing rescaled"),
           dict(name="Acetone", source="GISSMO", provenance="GISSMO bmse000032/1")]
    return G, P, lib


def test_table():
    G, P, lib = synthetic()
    # default: Quantified / Overlapped / Manual in >= 25 % of samples (citrate is only an estimate, acetone absent)
    assert list(T.build(G, P, lib).Metabolite) == ["Lactate", "Threonine"]
    assert list(T.build(G, P, lib, min_fraction=0.75).Metabolite) == ["Lactate"]
    tiers = T.REPORTED + ("Deconvolution estimate (low confidence)",)
    df = T.build(G, P, lib, tiers=tiers).set_index("Metabolite")
    assert list(df.index) == ["Citrate", "Lactate", "Threonine"]
    assert df.loc["Lactate", "1H chemical shift / ppm"] == "1.331*, 4.111"
    assert df.loc["Threonine", "1H chemical shift / ppm"] == "1.325**, 4.250"
    assert df.loc["Citrate", "1H chemical shift / ppm"] == "2.660"
    assert df.loc["Threonine", "Reported in (samples)"] == "1 of 2"
    assert df.loc["Lactate", "Reference spectrum"] == "GISSMO bmse000208"
    assert df.loc["Citrate", "Reference spectrum"] == "CASMDB peak list (HMDB)"


def test_file():
    G, P, lib = synthetic()
    Q = pd.DataFrame(dict(sample=["1", "2"], SF_MHz=[800.31, 800.31]))
    with tempfile.TemporaryDirectory() as d:
        df, note = T.make(d, G, P, Q, lib, "test", "TSP")
        x = pd.read_excel(os.path.join(d, T.FILE), header=None)
        assert x.iloc[0, 0].startswith("Chemical shift assignments")
        assert list(x.iloc[2, [0, 2, 3]]) == ["Metabolite", "Reported in (samples)", "Reference spectrum"]
        assert x.iloc[3:5, 0].tolist() == ["Lactate", "Threonine"]
        assert "at least 25% of the spectra" in note and "800 MHz" in note and "*" in note and "\u2014" not in note


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok ", name)
