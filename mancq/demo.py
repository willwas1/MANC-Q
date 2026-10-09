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


def imaginary_part(y):
    """The imaginary spectrum that goes with the real spectrum y (stored high ppm first, as TopSpin does), as
    TopSpin writes it to 1i: the transform of a time signal that starts at t = 0."""
    N = y.size
    s = np.fft.ifft(y)
    s2 = np.zeros_like(s)
    s2[0], s2[N // 2] = s[0], s[N // 2]
    s2[N // 2 + 1:] = 2 * s[N // 2 + 1:]
    return np.fft.fft(s2).imag


def write_bruker(folder, ppm0, sw_hz, sf, y, title):
    pdir = os.path.join(folder, "pdata", "1")
    os.makedirs(pdir, exist_ok=True)
    yi = imaginary_part(y)
    big = max(np.abs(y).max(), np.abs(yi).max())
    nc = int(np.ceil(np.log2(big / 2 ** 30))) if big > 2 ** 30 else 0
    np.round(y / 2.0 ** nc).astype("<i4").tofile(os.path.join(pdir, "1r"))
    np.round(yi / 2.0 ** nc).astype("<i4").tofile(os.path.join(pdir, "1i"))
    procs = {"SF": sf, "SW_p": sw_hz, "OFFSET": ppm0, "SI": len(y), "FTSIZE": len(y), "STSI": len(y), "STSR": 0,
             "NC_proc": nc, "BYTORDP": 0, "DTYPP": 0, "AXNUC": "<1H>", "WDW": 1, "LB": 0.0, "PHC0": 0.0,
             "PHC1": 0.0}
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


def make_fid(sf=600.13, seed=2, lw_hz=1.0, td=32768, ph0_deg=40.0, cache_dir=".mancq_cache"):
    """A raw FID of the same mixture (Bruker layout, no digital filter), deliberately out of phase by ph0_deg, so
    the FID-processing path can be demonstrated and tested."""
    rng = np.random.default_rng(seed)
    sw_ppm = 14.0
    sw_hz = sw_ppm * sf
    o1_ppm = 4.75                                     # carrier on water, as for presaturation
    lib = q.prepare_library(q.build_library(sf, os.path.join(q.HERE, "library"), cache_dir), sf)
    names = [c["name"] for c in lib]
    t = np.arange(td // 2) / sw_hz
    per_proton = 1.0 / 9 / REF_MM
    pos, amp = [0.0, -6.6 / 2 / sf, 6.6 / 2 / sf], [1 - 0.0467, 0.0467 / 2, 0.0467 / 2]   # TSP + 29Si satellites
    for met, mm in TRUTH_MM.items():
        c = lib[names.index(met)]
        d = rng.uniform(-2, 2) / sf
        for m in c["mults"]:
            pos.extend(list(np.asarray(m["ppm"]) + d)); amp.extend(list(mm * per_proton * np.asarray(m["inten"])))
    pos, amp = np.array(pos), np.array(amp)
    fid = np.zeros(t.size, complex)
    for k in range(0, len(pos), 200):
        nu = (pos[k:k + 200] - o1_ppm) * sf                 # Hz from the carrier
        fid += (amp[k:k + 200, None] * np.exp(2j * np.pi * nu[:, None] * t[None, :])).sum(0)
    fid *= np.exp(-np.pi * lw_hz * t) * np.exp(-1j * np.deg2rad(ph0_deg))
    fid += (rng.normal(0, 1, t.size) + 1j * rng.normal(0, 1, t.size)) * np.abs(fid[0]) / 3e4
    return fid, sw_hz, o1_ppm


def write_fid(folder, fid, sw_hz, sf, o1_ppm, title=""):
    os.makedirs(folder, exist_ok=True)
    scale = 2e8 / np.abs(fid).max()
    inter = np.empty(fid.size * 2); inter[0::2] = fid.real * scale; inter[1::2] = fid.imag * scale
    np.round(inter).astype("<i4").tofile(os.path.join(folder, "fid"))
    acqus = {"TD": fid.size * 2, "SW_h": sw_hz, "SW": sw_hz / sf, "SFO1": sf + o1_ppm * sf * 1e-6, "BF1": sf,
             "O1": o1_ppm * sf, "DTYPA": 0, "BYTORDA": 0, "GRPDLY": 0, "DECIM": 1, "DSPFVS": 20, "DIGMOD": 1,
             "AQ_mod": 3, "NS": 64, "NUC1": "<1H>", "PULPROG": "<zgpr>"}
    text = "##TITLE= Parameter file\n##JCAMPDX= 5.0\n" + "".join(f"##${k}= {v}\n" for k, v in acqus.items()) + "##END=\n"
    open(os.path.join(folder, "acqus"), "w").write(text)
    open(os.path.join(folder, "acqu"), "w").write(text)


def make_tocsy(sf=600.13, seed=3, si2=1024, si1=512, sw_ppm=12.0, ppm0=10.5):
    """A synthetic 1H-1H TOCSY of the same mixture: a diagonal peak for every proton and a cross peak for every
    pair of protons in one coupled spin system, at the library shifts (as the 2D check in the Results step
    expects). Returns the 2D array [F1, F2]."""
    from . import twod
    rng = np.random.default_rng(seed)
    f2 = ppm0 - np.arange(si2) * sw_ppm / si2
    f1 = ppm0 - np.arange(si1) * sw_ppm / si1
    z = np.zeros((si1, si2))
    w2, w1 = 0.006, 0.012                             # peak widths (ppm), broader in F1 as in a real TOCSY

    def put(x, y, h):
        gx = np.exp(-0.5 * ((f2 - x) / w2) ** 2)
        gy = np.exp(-0.5 * ((f1 - y) / w1) ** 2)
        z[:] += h * np.outer(gy, gx)
    put(0.0, 0.0, 9 * REF_MM)
    for met, mm in TRUTH_MM.items():
        e = twod.spin_system(met)
        if e is None:
            continue
        cs = np.asarray(e["shifts"], float)
        for c in cs:
            put(c, c, mm)
        for i, j in twod._groups(len(cs), e["J"], cosy=False):
            if abs(cs[i] - cs[j]) >= twod.MIN_SEP_PPM:
                put(cs[i], cs[j], 0.3 * mm); put(cs[j], cs[i], 0.3 * mm)
    z += rng.normal(0, 1, z.shape) * z.max() / 3000
    return z, f2, f1


def write_tocsy(folder, z, sf, sw_ppm=12.0, ppm0=10.5, title=""):
    """Write a 2D spectrum as TopSpin does (2rr in one submatrix, procs + proc2s, acqus + acqu2s)."""
    pdir = os.path.join(folder, "pdata", "1")
    os.makedirs(pdir, exist_ok=True)
    si1, si2 = z.shape
    scale = 2e8 / np.abs(z).max()
    np.round(z * scale).astype("<i4").tofile(os.path.join(pdir, "2rr"))
    def jcamp(p):
        return "##TITLE= Parameter file\n##JCAMPDX= 5.0\n" + "".join(f"##${k}= {v}\n" for k, v in p.items()) + "##END=\n"
    for fn, si in (("procs", si2), ("proc2s", si1)):
        p = {"SF": sf, "SW_p": sw_ppm * sf, "OFFSET": ppm0, "SI": si, "XDIM": si, "NC_proc": 0, "BYTORDP": 0,
             "DTYPP": 0, "AXNUC": "<1H>"}
        open(os.path.join(pdir, fn), "w").write(jcamp(p))
    for fn in ("acqus", "acqu2s"):
        a = {"PULPROG": "<mlevphpr.2>", "NUC1": "<1H>", "SFO1": sf, "TD": 2048, "SW": sw_ppm}
        open(os.path.join(folder, fn), "w").write(jcamp(a))
    open(os.path.join(pdir, "title"), "w").write(title)


def demo(out="mancq_demo", sf=600.13):
    spec_dir = os.path.join(out, "spectra")
    res_dir = os.path.join(out, "results")
    print(f"Writing two synthetic {sf:.2f} MHz spectra of {len(TRUTH_MM)} metabolites to {spec_dir}:\n"
          f"  EXPNO 10 processed (1r), EXPNO 11 raw FID only (processed by MANC-Q), EXPNO 12 a 2D TOCSY of 10 ...")
    write_demo(spec_dir, sf, cache_dir=os.path.join(out, "results", "_cache"))
    q.run(spec_dir, res_dir, REF_MM, dilution=1.0, N_WORKERS=2, RESUME=False)
    g = pd.read_csv(os.path.join(res_dir, "concentrations_long.csv"), dtype={"sample": str})
    rows = []
    for smp in ("10", "11"):
        h = g[(g["sample"] == smp) & g.metabolite.isin(TRUTH_MM)].set_index("metabolite")
        chk = pd.DataFrame({"true_mM": pd.Series(TRUTH_MM), "fitted_mM": h.conc_mM, "tier": h.tier})
        chk["ratio"] = chk.fitted_mM / chk.true_mM
        chk.insert(0, "sample", smp)
        rows.append(chk)
    chk = pd.concat(rows)
    chk.round(3).to_csv(os.path.join(out, "demo_check.csv"), index_label="metabolite")
    good_tiers = ["Quantified", "Overlapped (semi-quantitative)"]
    for smp, label in (("10", "processed spectrum"), ("11", "raw FID")):
        c = chk[chk["sample"] == smp].drop(columns="sample")
        print(f"\nEXPNO {smp} ({label}): recovered vs true (tube mM)")
        print(c.round(3).to_string())
        print(f"Median ratio over graded compounds: {c[c.tier.isin(good_tiers)].ratio.median():.3f}")
    print(f"\nResults in {res_dir}; comparison in {os.path.join(out, 'demo_check.csv')}")
    return chk


def write_demo(spec_dir, sf=600.13, cache_dir=".mancq_cache"):
    """Write the demo experiments: 10 (processed), 11 (raw FID only) and 12 (a TOCSY of the same sample as 10)."""
    ppm0, sw_hz, y = make_spectrum(sf, cache_dir=cache_dir)
    write_bruker(os.path.join(spec_dir, "10"), ppm0, sw_hz, sf, y, "MANC-Q synthetic demo (processed)")
    fid, fsw, o1 = make_fid(sf, cache_dir=cache_dir)
    write_fid(os.path.join(spec_dir, "11"), fid, fsw, sf, o1, "MANC-Q synthetic demo (raw FID)")
    z, _, _ = make_tocsy(sf)
    write_tocsy(os.path.join(spec_dir, "12"), z, sf, title="MANC-Q synthetic demo (processed)")
