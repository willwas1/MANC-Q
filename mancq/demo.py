"""Synthetic demonstration: builds a spectrum of known composition from the library, writes it in Bruker
format, runs the full quantification on it and compares the result with the truth. Useful as an installation
check and as a template for the expected input layout."""
import os
import numpy as np
import pandas as pd
import nmrglue as ng
from . import engine as q

TRUTH_MM = {  # concentrations in the NMR tube
    "D-Glucose": 5.0, "L-Lactate": 3.0, "L-Glutamine": 2.0, "Glycine": 1.0, "L-Alanine": 1.0, "L-Valine": 0.8,
    "L-Leucine": 0.8, "L-Threonine": 0.8, "L-Glutamate": 0.6, "L-Isoleucine": 0.5, "Acetate": 0.5,
    "Formate": 0.3, "L-Tyrosine": 0.3, "L-Phenylalanine": 0.3, "Creatine": 0.3, "Choline": 0.2,
    "L-Methionine": 0.2, "L-Histidine": 0.2,
}
REF_MM = 0.5


def write_bruker(folder, ppm0, sw_hz, sf, y, title):
    pdir = os.path.join(folder, "pdata", "1")
    os.makedirs(pdir, exist_ok=True)
    nc = int(np.ceil(np.log2(np.abs(y).max() / 2 ** 30))) if np.abs(y).max() > 2 ** 30 else 0
    np.round(y / 2.0 ** nc).astype("<i4").tofile(os.path.join(pdir, "1r"))
    procs = {"SF": sf, "SW_p": sw_hz, "OFFSET": ppm0, "SI": len(y), "FTSIZE": len(y), "STSI": len(y), "STSR": 0,
             "NC_proc": nc, "BYTORDP": 0, "DTYPP": 0, "AXNUC": "<1H>"}
    text = "##TITLE= Parameter file\n##JCAMPDX= 5.0\n" + "".join(f"##${k}= {v}\n" for k, v in procs.items()) + "##END=\n"
    for fn in ("procs", "proc"):
        open(os.path.join(pdir, fn), "w").write(text)
    open(os.path.join(pdir, "title"), "w").write(title)


def make_spectrum(sf=600.13, seed=1, lw_hz=1.0, eta=0.7, cache_dir=".mancq_cache"):
    rng = np.random.default_rng(seed)
    si, sw_ppm, ppm0 = 65536, 14.0, 11.5
    ppm = ppm0 - np.arange(si) * sw_ppm / si
    x = ppm[::-1]                                     # ascending for the model functions
    lib = q.prepare_library(q.build_library(sf, os.path.join(q.HERE, "library"), cache_dir), sf)
    names = [c["name"] for c in lib]
    tsp_area = 1.0e9
    y = np.zeros(si)
    # reference singlet at 0 ppm with 29Si satellites (4.67 %)
    fsi, jsi = q.S["TSP_SI29_FRACTION"], 6.6
    fw = lw_hz / sf
    y += tsp_area * (1 - fsi) * q.pv(x, 0.0, fw, eta)
    y += tsp_area * fsi / 2 * (q.pv(x, -jsi / 2 / sf, fw, eta) + q.pv(x, jsi / 2 / sf, fw, eta))
    per_proton = tsp_area / 9 / REF_MM                # area of one proton at 1 mM
    for met, mm in TRUTH_MM.items():
        c = lib[names.index(met)]
        d = rng.uniform(-2, 2) / sf                   # small compound-level shift, as in real samples
        for m in c["mults"]:
            y += mm * per_proton * q.mcol(x, m["ppm"], m["inten"], d, fw * rng.uniform(0.95, 1.2), eta)
    y += 2e4 * np.sin(x / 3.0) + 1e4                  # gentle baseline
    ref_height = y[np.abs(x) < 0.01].max()
    y += rng.normal(0, 1, si) * ref_height / 5000   # reference signal-to-noise about 5000
    return ppm0, sw_ppm * sf, y[::-1]


def demo(out="mancq_demo", sf=600.13):
    spec_dir = os.path.join(out, "spectra")
    res_dir = os.path.join(out, "results")
    print(f"Writing a synthetic {sf:.2f} MHz spectrum of {len(TRUTH_MM)} metabolites to {spec_dir} ...")
    ppm0, sw_hz, y = make_spectrum(sf, cache_dir=os.path.join(out, "results", "_cache"))
    write_bruker(os.path.join(spec_dir, "10"), ppm0, sw_hz, sf, y, "MANC-Q synthetic demo")
    q.run(spec_dir, res_dir, REF_MM, dilution=1.0, N_WORKERS=1, RESUME=False)
    g = pd.read_csv(os.path.join(res_dir, "concentrations_long.csv"))
    g = g[g.metabolite.isin(TRUTH_MM)].set_index("metabolite")
    chk = pd.DataFrame({"true_mM": pd.Series(TRUTH_MM), "fitted_mM": g.conc_mM, "tier": g.tier})
    chk["ratio"] = chk.fitted_mM / chk.true_mM
    chk.round(3).to_csv(os.path.join(out, "demo_check.csv"))
    print("\nRecovered vs true (tube mM):")
    print(chk.round(3).to_string())
    good = chk[chk.tier.isin(["Quantified", "Overlapped (semi-quantitative)"])]
    print(f"\nMedian ratio over graded compounds: {good.ratio.median():.3f}   "
          f"(results in {res_dir}; comparison in {os.path.join(out, 'demo_check.csv')})")
    return chk
