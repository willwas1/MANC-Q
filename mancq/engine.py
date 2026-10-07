"""
mancq.engine - 1H NMR metabolite quantification by whole-spectrum library fitting.

For each processed 1D 1H spectrum (Bruker pdata/<procno>/1r) the engine:
  1. fits the reference singlet (TSP or DSS, set to 0.000 ppm) for chemical-shift referencing,
     linewidth and the area used for absolute quantification;
  2. builds the compound library at the spectrometer frequency of the data: GISSMO spin systems
     (chemical shifts and J couplings, BMRB) are simulated exactly, so second-order multiplets are
     correct at any field; a few compounds use fixed peak lists instead;
  3. fits the whole spectrum as  sum(amount x compound spectrum) + smooth baseline + free
     "unassigned" singlets, with small per-compound and per-multiplet shift and linewidth freedom;
  4. converts amounts to mM with the reference area and concentration, applies dilution;
  5. grades every compound (Quantified / Overlapped / Deconvolution estimate / Upper bound only /
     Not detected) from objective measures, and writes tables and overlay plots.

Command line: see `python -m mancq --help`.
"""
# ----------------------------- SETTINGS --------------------------------------
SETTINGS = dict(
    DATA_DIR=None,                # folder of EXPNO folders, or a single EXPNO folder
    OUT_DIR=None,
    DILUTION_FILE=None,           # optional CSV: sample,dilution
    PROCNO=1,
    TARGET_POINT_HZ=0.25,         # spectra with finer digital resolution are block-averaged to this
    SAMPLES=None,                 # None = every EXPNO folder found, or a list such as ["10", "11"]
    N_WORKERS=None,               # parallel processes (None = number of CPU cores - 1)
    RESUME=True,                  # re-use per-sample results already in <out>/_per_sample
    GALLERY_MAX_SAMPLES=8,        # the per-metabolite gallery is only drawn for small runs
    LIBRARY_DIR=None,             # None = the library folder shipped with the package
    OVERLAYS=True,                # write one overlay PDF per spectrum
    PLURONIC="auto",              # Pluronic F-68 (poloxamer): "auto" = model it if its PEO line is present, "yes" = always
                                  # model it, "ignore" = leave its regions out of the fit, "no" = not in the samples
    PLURONIC_REGIONS=[(3.69, 3.74), (1.10, 1.20)],   # PEO line and PPO methyl, used when PLURONIC = "ignore"
    SAMPLE_DIRS=None,             # optional {sample: (EXPNO folder, procno)}, e.g. FIDs processed into another folder
    FID_MODE="missing",           # raw FIDs: "missing" = process an FID only where there is no processed spectrum,
                                  # "all" = always process from the FID, "never" = only use TopSpin-processed spectra
    FID_LB=0.3,                   # line broadening (Hz) when MANC-Q processes FIDs
    PEO_DETECT_SNR=500,           # "auto": PEO line at 3.71 ppm taller than this x noise
    PEAKLIST_FIELD_MHZ=800.0,     # field at which the fixed peak lists in the library are defined

    # internal standard (the only reference used): TSP or DSS singlet at 0.000 ppm
    REFERENCE="TSP",              # name used in reports only; TSP-d4 and DSS-d6 are handled identically
    TSP_MM_IN_TUBE=None,          # reference concentration in the NMR tube (mM) - REQUIRED
    TSP_PROTONS=9,
    TSP_SI29_FRACTION=0.0467,     # 29Si satellites are outside the fitted main line
    DEFAULT_DILUTION=1.0,         # sample -> tube dilution factor used when none is given

    # fitting
    FIT_RANGE=(0.60, 9.60),
    EXCLUDE=[(4.60, 5.00)],       # regions left out of the fit (ppm); default: residual water
    PEO_CORE_HZ=8.0,              # Pluronic PEO line centre excluded when taller than PEO_CORE_MIN_SNR x noise
    PEO_CORE_MIN_SNR=20000,
    MAX_SHIFT_HZ=8.0,             # compound-level shift freedom around library + global offset
    MAX_SHIFT_TITRATABLE_HZ=25.0, # pH-sensitive compounds (citrate, histidine, ...)
    MULTIPLET_SHIFT_HZ=8.0,       # extra independent freedom of each proton group (grid searched)
    SHIFT_PRIOR=0.02,             # relative SSR penalty per (10 Hz)^2 of multiplet movement
    MULTIPLET_SHIFT_TITRATABLE_HZ=6.0,
    COMPOUND_GRID_HZ=0.5,
    LINEWIDTH_BASE_HZ=(0.8, 1.5),  # base linewidth = TSP linewidth clipped to this range; each multiplet scales it
    ETA_SEARCH=[0.4, 0.6, 0.8, 1.0],  # Lorentzian fraction of the lineshape (searched once)
    GLOBAL_OFFSET_SEARCH_PPM=0.025,   # library-vs-TSP referencing offset search
    REL_MODEL_ERROR=0.03,         # weighting: sigma_i^2 = noise^2 + (0.03*|y_i|)^2 + (0.005*local max)^2
    TAIL_ERROR=0.005,
    TAIL_WINDOW_PPM=0.02,
    BASELINE_KNOT_PPM=0.04,
    BASELINE_PENALTY=3.0,
    # sharp, common compounds used to measure the library-vs-reference offset; an anchor is only used
    # if its pattern is actually present (correlation >= OFFSET_MIN_CORR)
    OFFSET_ANCHORS=["L-Lactate", "L-Alanine", "L-Valine", "Acetate", "Formate", "L-Isoleucine",
                    "D-Glucose", "L-Tyrosine", "L-Phenylalanine", "Choline", "L-Leucine"],
    OFFSET_MIN_CORR=0.3,
    OFFSET_MIN_ANCHORS=3,
    N_ITER=3,
    ROI_MIN_SNR=5,                # multiplets taller than this take part in local refinement
    ROI_PAD_HZ=6.0,
    ROI_MAX_NFEV=60,
    UNKNOWN_MIN_SNR=15,           # residual peaks above this become free "unassigned" singlets
    UNKNOWN_MIN_FRACTION=0.30,    # ... and only where the residual is >=30% of the local signal
    UNKNOWN_BLOCK_SNR=10,         # ... and not within UNKNOWN_EXCLUSION_HZ of a fitted library line
    UNKNOWN_EXCLUSION_HZ=3.0,
    COMPOUND_TIE=0.3,             # soft constraint: equal per-proton amplitude of multiplets in a compound
    MAX_UNKNOWNS=300,

    # grading
    LOD_SNR=3, LOQ_SNR=10,
    Q_MIN_DOMINANCE=0.60, Q_MAX_MISFIT=0.25, MAX_REPORTER_RATIO=2.0,
    MAX_CONC_OVER_BOUND=1.5,      # reported value may not exceed 1.5x the data upper bound
    SINGLE_REPORTER_MIN_SHARE=0.3,  # a single clean reporter must carry >=30% of the molecule's protons
    CONTRADICT_MIN_SNR=10, CONTRADICT_FRACTION=0.3,   # predicted multiplet >10 sigma but <30% observed -> reject
    O_MIN_DOMINANCE=0.25, O_MAX_MISFIT=0.50,
    DECONV_MAX_REL_SE=0.35, DECONV_MIN_MULTIPLETS=2,

    # plots
    OVERLAY_REGIONS=[(0.75, 1.10), (1.10, 1.55), (1.55, 2.00), (2.00, 2.35), (2.35, 2.75),
                     (2.75, 3.10), (3.10, 3.35), (3.35, 3.62), (3.62, 3.85), (3.85, 4.15),
                     (4.15, 4.60), (5.00, 5.50), (5.50, 6.60), (6.60, 7.10), (7.10, 7.50),
                     (7.50, 8.00), (8.00, 9.40)],
)
# -----------------------------------------------------------------------------

import os, sys, json, time, hashlib, warnings
import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear, least_squares, nnls
from scipy.interpolate import BSpline
from scipy.signal import find_peaks
import nmrglue as ng
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages


class _plt:
    """Backend-independent stand-in for the two pyplot calls used here, so importing the engine never changes
    the matplotlib backend (the GUI draws with Qt in the same process)."""
    @staticmethod
    def subplots(nrows=1, ncols=1, figsize=None, squeeze=True, **kw):
        fig = Figure(figsize=figsize)
        ax = fig.subplots(nrows, ncols, squeeze=squeeze, **kw)
        return fig, ax

    @staticmethod
    def close(fig=None):
        pass


plt = _plt

from .spinsim import simulate_spin_system, merge_sticks

__version__ = "1.1.2"
HERE = os.path.dirname(os.path.abspath(__file__))

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
S = SETTINGS
import copy as _copy
DEFAULT_SETTINGS = _copy.deepcopy(SETTINGS)
TIMING = {}


# =============================================================================
# 1. Reading
# =============================================================================
DECONV = "Deconvolution estimate (low confidence)"
NOTMEAS = "Not measurable (region ignored)"


def _json_default(v):
    """Lets json.dumps write numpy numbers and arrays (used for the saved fit state)."""
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, np.generic):
        return v.item()
    raise TypeError(f"cannot save {type(v)}")


def read_spectrum(expdir, procno):
    pdir = os.path.join(expdir, "pdata", str(procno))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dic, data = ng.bruker.read_pdata(pdir, scale_data=True, read_acqus=False)
    p = dic["procs"]
    if np.ndim(data) != 1:
        raise ValueError(f"{pdir}: not a 1D spectrum")
    si = len(data)
    # TopSpin convention: point i is at OFFSET - i * SW_p / SF / SI  (ppm)
    ppm = float(p["OFFSET"]) - np.arange(si) * float(p["SW_p"]) / float(p["SF"]) / si
    y = np.asarray(data, float)
    o = np.argsort(ppm)
    ppm, y = ppm[o], y[o]
    sf = float(dic["procs"]["SF"])
    # spectra processed with a large SI are block-averaged to ~TARGET_POINT_HZ per point:
    # linewidths are >= 1 Hz, so nothing is lost, and memory/time stay the same for any SI
    dx_hz = abs(ppm[1] - ppm[0]) * sf
    k = max(1, int(round(S["TARGET_POINT_HZ"] / dx_hz)))
    if k >= 2:
        n = (len(ppm) // k) * k
        ppm = ppm[:n].reshape(-1, k).mean(1)
        y = y[:n].reshape(-1, k).mean(1)
    o = np.arange(len(ppm))
    title = ""
    tf = os.path.join(pdir, "title")
    if os.path.exists(tf):
        title = open(tf, errors="ignore").read().strip()
    return ppm, y, sf, title


def _usable(d, procno, fids):
    return os.path.exists(os.path.join(d, "pdata", str(procno), "1r")) or \
        (fids and os.path.exists(os.path.join(d, "fid")) and os.path.exists(os.path.join(d, "acqus")))


def list_samples(data_dir, procno, fids=False):
    """EXPNO folders under data_dir with pdata/<procno>/1r (or, if fids, a raw fid + acqus). If data_dir is itself
    an EXPNO folder, returns [""]."""
    if os.path.exists(os.path.join(data_dir, "acqus")) or os.path.isdir(os.path.join(data_dir, "pdata")):
        return [""] if _usable(data_dir, procno, fids) else []
    out = []
    for d in sorted(os.listdir(data_dir), key=lambda x: (len(x), x)):
        if os.path.isdir(os.path.join(data_dir, d)) and _usable(os.path.join(data_dir, d), procno, fids):
            out.append(d)
    return out


def read_dilutions(path):
    if not path:
        return {}
    if not os.path.exists(path):
        raise FileNotFoundError(f"dilution file not found: {path}")
    df = pd.read_csv(path)
    return {str(k).strip(): float(v) for k, v in zip(df.iloc[:, 0], df.iloc[:, 1])}


# =============================================================================
# 2. Lineshape
# =============================================================================
def pv(x, x0, fwhm, eta):
    """Area-normalised pseudo-Voigt (x, x0, fwhm in the same units)."""
    g = fwhm / 2.0
    lor = (g / np.pi) / ((x - x0) ** 2 + g ** 2)
    sig = fwhm / (2 * np.sqrt(2 * np.log(2)))
    gau = np.exp(-0.5 * ((x - x0) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
    return eta * lor + (1 - eta) * gau


def noise_sigma(ppm, y):
    parts = []
    for lo, hi in [(10.0, 12.0), (-2.5, -0.5)]:
        m = (ppm > lo) & (ppm < hi)
        if m.sum() > 200:
            yy = y[m] - np.polyval(np.polyfit(ppm[m], y[m], 2), ppm[m])
            parts.append(1.4826 * np.median(np.abs(yy - np.median(yy))))
    return float(np.min(parts)) if parts else float(np.std(y[:2000]))


def fit_tsp(ppm, y, sf):
    """Pseudo-Voigt main line + 29Si doublet (J 6.6 Hz) + linear baseline."""
    m = (ppm > -0.2) & (ppm < 0.2)
    x0 = ppm[m][np.argmax(y[m])]
    w = (ppm > x0 - 0.03) & (ppm < x0 + 0.03)
    xx, yy = ppm[w], y[w]
    h = yy.max()
    fsi = S["TSP_SI29_FRACTION"]

    def model(p):
        c, f, eta, A, b0, b1, jsi = p
        fp = f / sf
        main = A * pv(xx, c, fp, eta)
        sat = A * fsi / (1 - fsi) / 2 * (pv(xx, c - jsi / 2 / sf, fp, eta) + pv(xx, c + jsi / 2 / sf, fp, eta))
        return main + sat + b0 + b1 * (xx - c)

    A0 = h * 1.2 / sf * np.pi / 2
    p0 = [x0, 1.2, 0.6, A0, 0, 0, 6.6]
    lb = [x0 - 0.002, 0.3, 0.0, 0, -h, -np.inf, 5.5]
    ub = [x0 + 0.002, 6.0, 1.0, np.inf, h, np.inf, 7.5]
    r = least_squares(lambda p: model(p) - yy, p0, bounds=(lb, ub), x_scale="jac")
    c, f, eta, A, b0, b1, jsi = r.x
    resid = np.sqrt(np.mean(r.fun ** 2)) / h
    area_total = A / (1 - fsi)          # main + 29Si satellites = all TSP protons
    return dict(centre=c, fwhm_hz=f, eta=eta, area_main=A, area=area_total, height=h, rel_rms=resid, jsi=jsi)


# =============================================================================
# 3. Library
# =============================================================================
def parse_gissmo(path):
    lib = {}
    for line in open(path):
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        gid, name, rest = line.split("|", 2)
        p = rest.split("#")
        cs = [float(x) for x in p[1].split(",")]
        J = {}
        for x in p[2].split(";"):
            if x:
                a, v = x.split(":")
                i, j = a.split("-")
                J[(int(i) - 1, int(j) - 1)] = float(v)
        het = None
        if len(p) > 3 and p[3].startswith("ADD"):
            het = {int(k) - 1: float(v) for k, v in (e.split("=") for e in p[3][3:].split(","))}
        lib[gid] = dict(name=name, shifts=cs, J=J, hetero=het)
    return lib


ANOMER_WEIGHTS = {"D-Glucose": {5.2223: 0.36, 4.6360: 0.64}}


def build_library(sf, libdir, cache_dir):
    key = f"library_{sf:.3f}MHz.json"
    cpath = os.path.join(cache_dir, key)
    comp = pd.read_csv(os.path.join(libdir, "compounds.csv"))
    sig = hashlib.md5((open(os.path.join(libdir, "compounds.csv")).read() +
                       open(os.path.join(libdir, "gissmo_spin_systems.txt")).read() +
                       open(os.path.join(libdir, "other_spin_systems.txt")).read() +
                       open(os.path.join(libdir, "peaklists.csv")).read() + __version__).encode()).hexdigest()
    if os.path.exists(cpath):
        d = json.load(open(cpath))
        if d.get("sig") == sig:
            return d["compounds"]
    g = parse_gissmo(os.path.join(libdir, "gissmo_spin_systems.txt"))
    other = parse_gissmo(os.path.join(libdir, "other_spin_systems.txt"))
    pk = pd.read_csv(os.path.join(libdir, "peaklists.csv"))
    man = pd.read_csv(os.path.join(libdir, "peaklists_manifest.csv")).set_index("metabolite")
    out = []
    for _, r in comp.iterrows():
        if r.source in ("GISSMO", "SPIN"):
            e = g[r.source_id] if r.source == "GISSMO" else other[r.source_id]
            ppm, it, grp = simulate_spin_system(e["shifts"], e["J"], sf, e["hetero"], ANOMER_WEIGHTS.get(r.metabolite))
            ppm, it, grp = merge_sticks(ppm, it, grp)
            keep = it > 2e-3 * it.max()
            ppm, it, grp = ppm[keep], it[keep], grp[keep]
            nprot = len(e["shifts"]) if r.metabolite not in ANOMER_WEIGHTS else len(e["shifts"]) / 2
            it = it * nprot / it.sum()
            prov = (f"GISSMO bmse{r.source_id} ({e['name']}), QM-simulated at {sf:.2f} MHz" if r.source == "GISSMO" else
                    f"spin system {r.source_id} ({e['name']}, see other_spin_systems.txt), QM-simulated at {sf:.2f} MHz")
        else:
            sub = pk[pk.metabolite == r.source_id]
            ppm = sub.ppm.values * 1.0
            it = sub.height.values * 1.0
            grp = sub.cluster.values.astype(int)
            # peak lists are defined at PEAKLIST_FIELD_MHZ: rescale the line spacing within each cluster
            # (first-order approximation: J in Hz is field independent)
            scale = S["PEAKLIST_FIELD_MHZ"] / sf
            for g_ in np.unique(grp):
                m_ = grp == g_
                c_ = (ppm[m_] * it[m_]).sum() / it[m_].sum()
                ppm[m_] = c_ + (ppm[m_] - c_) * scale
            nprot = float(man.loc[r.source_id, "n_protons"])
            it = it * nprot / it.sum()
            prov = (f"CASMDB {man.loc[r.source_id,'casmdb_name']} peak list (origin {man.loc[r.source_id,'origin']}), "
                    f"spacing rescaled to {sf:.2f} MHz")
        ov = {}
        for col in ("max_shift_hz", "multiplet_shift_hz"):
            if col in r and pd.notna(r[col]) and str(r[col]).strip() != "":
                ov[col] = float(r[col])
        out.append(dict(name=r.metabolite, cls=r["class"], source=r.source, titratable=int(r.titratable), **ov,
                        nprot=float(nprot), ppm=[float(v) for v in ppm], inten=[float(v) for v in it], group=[int(v) for v in grp],
                        provenance=prov))
    os.makedirs(cache_dir, exist_ok=True)
    json.dump(dict(sig=sig, compounds=out), open(cpath, "w"))
    return out


def prepare_library(lib, sf, pluronic=False):
    """Attach multiplet groupings and, if requested, add the Pluronic F-68 polymer terms
    (a PEO/PPO surfactant common in cell-culture media)."""
    lib = [dict(c) for c in lib]
    jch = 142.0 / sf
    for tag, wf0, wlo, whi in ([("sharp", 1.3, 0.8, 5.0), ("broad", 25.0, 10.0, 60.0)] if pluronic else []):
        lib.append(dict(name=f"Pluronic F-68 PEO ({tag})", cls="polymer", source="polymer", titratable=0, nprot=1.0,
                        ppm=[3.7128 - jch / 2, 3.7128, 3.7128 + jch / 2], inten=[0.0055, 1.0, 0.0055], group=[1, 1, 1],
                        provenance="free polymer term (PEO -CH2CH2O-), position/width fitted, with 13C satellites",
                        wf0=wf0, wlo=wlo, whi=whi, dmax_hz=3.0))
    j = 6.3 / sf
    if pluronic:
        lib.append(dict(name="Pluronic F-68 PPO (CH3)", cls="polymer", source="polymer", titratable=0, nprot=1.0,
                        ppm=[1.150 - j / 2, 1.150 + j / 2], inten=[0.5, 0.5], group=[1, 1],
                        provenance="free polymer term (PPO methyl doublet, broad), position/width fitted",
                        wf0=6.0, wlo=2.0, whi=25.0, dmax_hz=30.0))
    for c in lib:
        if c["source"] == "polymer":
            c["mults"] = [dict(ppm=np.array(c["ppm"]), inten=np.array(c["inten"]), lo=min(c["ppm"]), hi=max(c["ppm"]),
                               centre=float(np.mean(c["ppm"])), nprot=float(sum(c["inten"])))]
        else:
            c["mults"] = group_multiplets(c["ppm"], c["inten"], c["group"], sf)
    return lib


def group_multiplets(ppm, inten, group, sf):
    """One multiplet per proton group (from the spin simulation / CASMDB cluster);
    a group wider than 0.10 ppm is split at its largest gap."""
    ppm, inten, group = np.asarray(ppm, float), np.asarray(inten, float), np.asarray(group)
    out = []
    for g in np.unique(group):
        sel = group == g
        if inten[sel].sum() <= 0:
            continue
        out.extend(split_multiplets(ppm[sel], inten[sel], sf, gap_hz=60.0, max_span_ppm=0.10))
    out.sort(key=lambda m: m["centre"])
    return out


def split_multiplets(ppm, inten, sf, gap_hz=20.0, max_span_ppm=0.08):
    """Group stick lines into multiplets: lines closer than gap_hz belong together;
    clusters wider than max_span_ppm are split at their largest internal gap
    (so each multiplet gets its own small shift/width freedom)."""
    o = np.argsort(ppm)
    p, it = np.asarray(ppm, float)[o], np.asarray(inten, float)[o]
    gap = gap_hz / sf
    groups, cur = [], [0]
    for k in range(1, len(p)):
        if p[k] - p[k - 1] > gap:
            groups.append(cur)
            cur = []
        cur.append(k)
    groups.append(cur)

    def split(g):
        if len(g) < 2 or p[g[-1]] - p[g[0]] <= max_span_ppm:
            return [g]
        d = np.diff(p[g])
        cut = int(np.argmax(d)) + 1
        return split(g[:cut]) + split(g[cut:])

    out = []
    for g in groups:
        for gg in split(g):
            pp, ii = p[gg], it[gg]
            out.append(dict(ppm=pp, inten=ii, lo=float(pp.min()), hi=float(pp.max()),
                            centre=float((pp * ii).sum() / ii.sum()), nprot=float(ii.sum())))
    return out


# =============================================================================
# 4. Model pieces
# =============================================================================
HALF_WIN = 0.06      # ppm evaluated either side of a line (48 Hz at 800 MHz)


def mcol(x, p, it, d, fw, eta):
    """Column (on sorted x) of one multiplet: sticks p (ppm), intensities it,
    shifted by d (ppm), pseudo-Voigt fwhm fw (ppm)."""
    out = np.zeros(len(x))
    c = p + d
    lo = np.searchsorted(x, c.min() - HALF_WIN)
    hi = np.searchsorted(x, c.max() + HALF_WIN)
    if hi > lo:
        xx = x[lo:hi]
        if len(c) * (hi - lo) > 4e6:
            for cc, ii in zip(c, it):
                out[lo:hi] += ii * pv(xx, cc, fw, eta)
        else:
            out[lo:hi] = (it[:, None] * pv(xx[None, :], c[:, None], fw, eta)).sum(0)
    return out


def mcol_d(x, p, it, d, fw, eta):
    """Multiplet column and its analytic derivatives w.r.t. shift d and fwhm fw."""
    n = len(x)
    col, dd, df = np.zeros(n), np.zeros(n), np.zeros(n)
    c = p + d
    lo = np.searchsorted(x, c.min() - HALF_WIN)
    hi = np.searchsorted(x, c.max() + HALF_WIN)
    if hi <= lo:
        return col, dd, df
    xx = x[lo:hi][None, :]
    u = xx - c[:, None]
    g = fw / 2.0
    den = u ** 2 + g ** 2
    L = (g / np.pi) / den
    dL_dx0 = (g / np.pi) * 2 * u / den ** 2
    dL_df = 0.5 * (1 / np.pi) * (u ** 2 - g ** 2) / den ** 2
    s = fw / (2 * np.sqrt(2 * np.log(2)))
    z = u / s
    Gs = np.exp(-0.5 * z ** 2) / (s * np.sqrt(2 * np.pi))
    dG_dx0 = Gs * u / s ** 2
    dG_df = Gs * (z ** 2 - 1) / fw
    w_ = it[:, None]
    col[lo:hi] = (w_ * (eta * L + (1 - eta) * Gs)).sum(0)
    dd[lo:hi] = (w_ * (eta * dL_dx0 + (1 - eta) * dG_dx0)).sum(0)
    df[lo:hi] = (w_ * (eta * dL_df + (1 - eta) * dG_df)).sum(0)
    return col, dd, df


def spline_basis(x_all, lo, hi, knot):
    k = 3
    inner = np.arange(lo, hi + knot, knot)
    t = np.r_[[inner[0]] * k, inner, [inner[-1]] * k]
    n = len(t) - k - 1
    B = BSpline.design_matrix(np.clip(x_all, inner[0], inner[-1]), t, k).toarray()
    D = np.diff(np.eye(n), 2, axis=0)
    return B, D


class Grid:
    def __init__(self, ppm, mask):
        self.idx = np.where(mask)[0]
        self.x = ppm[self.idx]
        self.dx = abs(ppm[1] - ppm[0])


# =============================================================================
# 5. Fit one spectrum
# =============================================================================
def fit_spectrum(ppm_raw, y, sf, lib, log):
    """Fit one spectrum. Model:
         y(x) = sum_m a_m * M_m(x; d_m, w_m)  +  sum_u b_u * singlet_u(x)  +  spline(x)
       M_m = simulated multiplet m (library intensities, locked internal geometry)
       a_m >= 0, softly tied to be equal (per proton) within a compound
       d_m = library position + global offset + compound shift D_k + small delta_m
       w_m = linewidth factor (pseudo-Voigt, shared eta)"""
    t0 = time.time()
    tsp = fit_tsp(ppm_raw, y, sf)
    ppm = ppm_raw - tsp["centre"]
    sigma = noise_sigma(ppm, y)
    fw0 = float(np.clip(tsp["fwhm_hz"], S["LINEWIDTH_BASE_HZ"][0], S["LINEWIDTH_BASE_HZ"][1])) / sf
    eta = float(np.clip(tsp["eta"], 0.3, 1.0))
    log(f"   {S['REFERENCE']}: area={tsp['area']:.4g}  fwhm={tsp['fwhm_hz']:.2f} Hz  eta={tsp['eta']:.2f}  "
        f"rel.rms={tsp['rel_rms']:.3f}  noise sigma={sigma:.4g}  SNR(TSP)={tsp['height']/sigma:.0f}")

    lo, hi = S["FIT_RANGE"]
    mask = (ppm > lo) & (ppm < hi)
    for a_, b_ in S["EXCLUDE"]:
        mask &= ~((ppm > a_) & (ppm < b_))
    # core of the Pluronic PEO line (often >100x any metabolite): its exact lineshape is
    # not a pseudo-Voigt, so the central +/-PEO_CORE_HZ is left out of the fit
    pm = (ppm > 3.69) & (ppm < 3.74)
    peo_apex = float(ppm[pm][np.argmax(y[pm])]) if pm.any() else None
    if peo_apex is not None and y[pm].max() > S["PEO_CORE_MIN_SNR"] * sigma:
        mask &= ~((ppm > peo_apex - S["PEO_CORE_HZ"] / sf) & (ppm < peo_apex + S["PEO_CORE_HZ"] / sf))
        log(f"   Pluronic PEO line at {peo_apex:.4f} ppm ({y[pm].max()/sigma:.0f} x noise): core +/-{S['PEO_CORE_HZ']:.0f} Hz excluded")
    grid = Grid(ppm, mask)
    x, yv = grid.x, y[grid.idx]
    from scipy.ndimage import maximum_filter1d
    local_max = maximum_filter1d(np.abs(yv), size=max(3, int(S["TAIL_WINDOW_PPM"] / grid.dx)))
    wts = 1.0 / np.sqrt(sigma ** 2 + (S["REL_MODEL_ERROR"] * np.abs(yv)) ** 2 + (S["TAIL_ERROR"] * local_max) ** 2)
    chi = lambda f: float((((yv - f) * wts) ** 2).sum() / len(yv))

    # ---------- global referencing offset of the GISSMO library ----------
    offs = []
    for c in lib:
        if c["name"] in S["OFFSET_ANCHORS"] and c["source"] == "GISSMO":
            best = None
            for d in np.arange(-S["GLOBAL_OFFSET_SEARCH_PPM"], S["GLOBAL_OFFSET_SEARCH_PPM"] + 1e-9, 0.25 / sf):
                col = sum(mcol(x, m["ppm"], m["inten"], d, fw0, eta) for m in c["mults"])
                nrm = np.sqrt((col * col).sum())
                if nrm > 0:
                    sc = (col * yv).sum() / nrm
                    if best is None or sc > best[0]:
                        best = (sc, d, col)
            if best is None:
                continue
            # use the anchor only if its pattern is clearly present at the best position
            col = best[2]
            sup = col > 0.02 * col.max()
            cc = np.corrcoef(col[sup], yv[sup])[0, 1] if sup.sum() >= 5 else 0.0
            if np.isfinite(cc) and cc >= S["OFFSET_MIN_CORR"]:
                offs.append(best[1])
    if len(offs) >= S["OFFSET_MIN_ANCHORS"]:
        g_off = float(np.median(offs))
        log(f"   GISSMO library offset: {g_off*sf:+.1f} Hz ({g_off:+.4f} ppm) from {len(offs)} anchors")
    else:
        g_off = 0.0
        log(f"   WARNING: only {len(offs)} offset anchors found; GISSMO library offset set to 0 "
            f"(per-compound shift freedom still applies)")

    # ---------- multiplet state ----------
    ncomp = len(lib)
    MP = []
    DMAX, DLT, Dk = np.zeros(ncomp), np.zeros(ncomp), np.zeros(ncomp)
    for k, c in enumerate(lib):
        base = g_off if c["source"] in ("GISSMO", "SPIN") else 0.0
        DMAX[k] = c.get("max_shift_hz", S["MAX_SHIFT_TITRATABLE_HZ"] if c["titratable"] else S["MAX_SHIFT_HZ"]) / sf
        DLT[k] = c.get("multiplet_shift_hz", S["MULTIPLET_SHIFT_TITRATABLE_HZ"] if c["titratable"] else S["MULTIPLET_SHIFT_HZ"]) / sf
        if c["source"] == "polymer":
            DMAX[k], DLT[k] = c.get("dmax_hz", 3.0) / sf, 0.5 / sf
        for j, m in enumerate(c["mults"]):
            MP.append(dict(k=k, j=j, p=m["ppm"], it=m["inten"], nprot=m["nprot"], base=base,
                           d=base, dlo=base - DLT[k], dhi=base + DLT[k],
                           wf=c.get("wf0", 1.0), wlo=c.get("wlo", 0.5), whi=c.get("whi", 4.0),
                           titr=c["titratable"], local_amp=np.nan, local_se=np.nan,
                           loc_dom=0.0, loc_mis=np.inf, loc_h=0.0))
    nM = len(MP)
    MPk = {k: [i for i, m in enumerate(MP) if m["k"] == k] for k in range(ncomp)}
    B, Dpen = spline_basis(x, lo, hi, S["BASELINE_KNOT_PPM"])
    unknown = []

    def mult_columns():
        M = np.zeros((len(x), nM))
        for i, m in enumerate(MP):
            M[:, i] = mcol(x, m["p"], m["it"], m["d"], fw0 * m["wf"], eta)
        return M

    def unk_columns():
        if not unknown:
            return np.zeros((len(x), 0))
        return np.array([mcol(x, np.array([u["pos"]]), np.array([1.0]), 0.0, fw0 * u["wf"], eta) for u in unknown]).T

    def global_solve(M, U):
        A = np.hstack([M, U, B])
        nnon = nM + U.shape[1]
        Aw = A * wts[:, None]
        cn = np.sqrt((Aw ** 2).sum(0))
        rows = []
        # baseline curvature penalty
        P = np.zeros((Dpen.shape[0], A.shape[1]))
        P[:, nnon:] = Dpen * S["BASELINE_PENALTY"] * np.median(cn[nnon:][cn[nnon:] > 0])
        rows.append(P)
        # soft tie: multiplets of one compound should have equal per-proton amplitude
        lam = S["COMPOUND_TIE"]
        for k, ii in MPk.items():
            if len(ii) < 2:
                continue
            for i1, i2 in zip(ii[:-1], ii[1:]):
                s_ = lam * np.sqrt(cn[i1] * cn[i2])
                if s_ <= 0:
                    continue
                r_ = np.zeros(A.shape[1])
                r_[i1], r_[i2] = s_, -s_
                rows.append(r_[None, :])
        Pr = np.vstack(rows)
        cs = np.sqrt((Aw ** 2).sum(0) + (Pr ** 2).sum(0))
        cs[cs == 0] = 1
        # normal equations -> Cholesky -> small bounded least squares (same solution, much faster)
        Aw /= cs
        Pr = Pr / cs
        G = Aw.T @ Aw + Pr.T @ Pr
        h = Aw.T @ (yv * wts)
        G[np.diag_indices_from(G)] += 1e-10 * np.trace(G) / len(G)
        try:
            R = np.linalg.cholesky(G).T
        except np.linalg.LinAlgError:
            G[np.diag_indices_from(G)] += 1e-6 * np.trace(G) / len(G)
            R = np.linalg.cholesky(G).T
        z = np.linalg.solve(R.T, h)
        lb = np.r_[np.zeros(nnon), -np.inf * np.ones(B.shape[1])]
        _t = time.time()
        r = lsq_linear(R, z, bounds=(lb, np.inf), method="bvls", max_iter=5000)
        TIMING["solve"] = TIMING.get("solve", 0) + time.time() - _t
        coef = r.x / cs
        return coef, A @ coef, A

    def roi_refine(coef, M, U, fit):
        am = coef[:nM]
        ua = coef[nM:nM + len(unknown)]
        base_fit = B @ coef[nM + len(unknown):]
        for m in MP:
            m.update(local_amp=np.nan, local_se=np.nan, loc_dom=0.0, loc_mis=np.inf, loc_h=0.0)
        act = [i for i in range(nM) if am[i] * M[:, i].max() > S["ROI_MIN_SNR"] * sigma]
        items = [("m", i, MP[i]["p"].min() + MP[i]["d"], MP[i]["p"].max() + MP[i]["d"]) for i in act]
        items += [("u", i, u["pos"], u["pos"]) for i, u in enumerate(unknown) if ua[i] > 0]
        pad = S["ROI_PAD_HZ"] / sf
        items.sort(key=lambda t: t[2])
        rois, cur, cur_hi = [], [], -1e9
        for it_ in items:
            if cur and it_[2] - pad > cur_hi:
                rois.append(cur)
                cur, cur_hi = [], -1e9
            cur.append(it_)
            cur_hi = max(cur_hi, it_[3] + pad)
        if cur:
            rois.append(cur)
        signal = fit - base_fit
        for roi in rois:
            rlo = min(t[2] for t in roi) - pad
            rhi = max(t[3] for t in roi) + pad
            w = (x >= rlo) & (x <= rhi)
            xs = x[w]
            if len(xs) < 5:
                continue
            ms = [t[1] for t in roi if t[0] == "m"]
            us = [t[1] for t in roi if t[0] == "u"]
            own = M[w][:, ms] @ am[ms] + (U[w][:, us] @ ua[us] if us else 0)
            target = yv[w] - base_fit[w] - (signal[w] - own)
            ww = wts[w]
            nm, nu = len(ms), len(us)
            v = np.array([MP[i]["d"] for i in ms] + [MP[i]["wf"] for i in ms] +
                         [unknown[i]["pos"] for i in us] + [unknown[i]["wf"] for i in us], float)
            lbv = np.array([MP[i]["dlo"] for i in ms] + [MP[i]["wlo"] for i in ms] +
                           [unknown[i]["pos"] - 1.0 / sf for i in us] + [0.7] * nu, float)
            ubv = np.array([MP[i]["dhi"] for i in ms] + [MP[i]["whi"] for i in ms] +
                           [unknown[i]["pos"] + 1.0 / sf for i in us] + [4.0] * nu, float)
            v = np.clip(v, lbv, ubv)
            tt = (xs - xs.mean()) / max(np.ptp(xs), 1e-9)
            blin = np.array([np.ones(len(xs)), -np.ones(len(xs)), tt, -tt]).T

            def lin(v):
                cols_, dd_, dw_ = [], [], []
                for a_, i in enumerate(ms):
                    c0, c1, c2 = mcol_d(xs, MP[i]["p"], MP[i]["it"], v[a_], fw0 * v[nm + a_], eta)
                    cols_.append(c0); dd_.append(c1); dw_.append(c2 * fw0)
                for a_, i in enumerate(us):
                    c0, c1, c2 = mcol_d(xs, np.array([0.0]), np.array([1.0]), v[2 * nm + a_], fw0 * v[2 * nm + nu + a_], eta)
                    cols_.append(c0); dd_.append(c1); dw_.append(c2 * fw0)
                C = np.array(cols_).T
                A_ = np.hstack([C, blin]) * ww[:, None]
                c_, _ = nnls(A_, target * ww, maxiter=50 * A_.shape[1])
                return c_, A_ @ c_ - target * ww, np.array(dd_).T, np.array(dw_).T, C, A_

            c_, r_, Dd, Dw, C, A_ = lin(v)
            ssr = (r_ ** 2).sum()
            lam = 1e-3
            ntot = nm + nu
            idx_d = list(range(nm)) + list(range(2 * nm, 2 * nm + nu))
            idx_w = list(range(nm, 2 * nm)) + list(range(2 * nm + nu, 2 * nm + 2 * nu))
            for _ in range(S["ROI_MAX_NFEV"]):
                amps = c_[:ntot]
                J = np.hstack([Dd * amps[None, :], Dw * amps[None, :]]) * ww[:, None]
                JTJ = J.T @ J
                g_ = J.T @ r_
                diag = np.diag(JTJ).copy()
                diag[diag <= 0] = 1.0
                improved, rel = False, 0.0
                for _t in range(6):
                    try:
                        step = np.linalg.solve(JTJ + lam * np.diag(diag), g_)
                    except np.linalg.LinAlgError:
                        lam *= 10
                        continue
                    vn = v.copy()
                    vn[idx_d] -= step[:ntot]
                    vn[idx_w] -= step[ntot:]
                    vn = np.clip(vn, lbv, ubv)
                    out = lin(vn)
                    sn = (out[1] ** 2).sum()
                    if sn < ssr:
                        rel = (ssr - sn) / ssr
                        v, (c_, r_, Dd, Dw, C, A_), ssr = vn, out, sn
                        lam = max(lam / 3, 1e-7)
                        improved = True
                        break
                    lam *= 8
                if not improved or rel < 1e-5:
                    break
            for a_, i in enumerate(ms):
                MP[i]["d"], MP[i]["wf"] = float(v[a_]), float(v[nm + a_])
            for a_, i in enumerate(us):
                unknown[i]["pos"], unknown[i]["wf"] = float(v[2 * nm + a_]), float(v[2 * nm + nu + a_])
            s2 = ssr / max(len(r_) - A_.shape[1], 1)
            try:
                se_ = np.sqrt(np.clip(np.diag(np.linalg.pinv(A_.T @ A_)) * s2, 0, None))
            except Exception:
                se_ = np.full(A_.shape[1], np.nan)
            local_model = C @ c_[:ntot] + blin @ c_[ntot:]
            for a_, i in enumerate(ms):
                fwp = fw0 * v[nm + a_]
                lw_ = (xs >= MP[i]["p"].min() + v[a_] - 2.5 * fwp) & (xs <= MP[i]["p"].max() + v[a_] + 2.5 * fwp)
                own_ = c_[a_] * C[lw_, a_]
                osum = own_.sum()
                tot = np.clip(C[lw_] @ c_[:ntot], 0, None).sum()
                MP[i].update(local_amp=float(c_[a_]), local_se=float(se_[a_]),
                             loc_dom=float(osum / tot) if tot > 0 else 0.0,
                             loc_mis=float(np.abs(target[lw_] - local_model[lw_]).sum() / osum) if osum > 0 else np.inf,
                             loc_h=float(own_.max()) if own_.size else 0.0)
        return len(rois)

    def compound_stage(coef, fit, M):
        """Grid search of the compound-level shift D_k (all multiplets of a compound
        move together, keeping their relative amplitudes and internal offsets)."""
        am = coef[:nM]
        fit = fit.copy()
        step = S["COMPOUND_GRID_HZ"] / sf
        strength = np.array([sum(am[i] * M[:, i].max() for i in MPk[k]) for k in range(ncomp)])
        moved = 0
        for k in np.argsort(-strength):
            if (lib[k]["source"] == "polymer" and lib[k].get("dmax_hz", 3.0) <= 3.0) or not MPk[k]:
                continue
            ii = MPk[k]
            spans = [(MP[i]["p"].min() + MP[i]["base"] - DMAX[k] - DLT[k] - HALF_WIN / 2,
                      MP[i]["p"].max() + MP[i]["base"] + DMAX[k] + DLT[k] + HALF_WIN / 2) for i in ii]
            sup = np.unique(np.concatenate([np.arange(np.searchsorted(x, a_), np.searchsorted(x, b_)) for a_, b_ in spans]))
            if sup.size == 0:
                continue
            if strength[k] <= 0 and yv[sup].max() < S["ROI_MIN_SNR"] * sigma:
                continue
            xs = x[sup]
            ww2 = wts[sup] ** 2
            own = M[sup][:, ii] @ am[ii]
            r = yv[sup] - fit[sup] + own
            base_ssr = (ww2 * r * r).sum()
            rel_amp = np.array([am[i] for i in ii])
            if rel_amp.sum() <= 0:
                rel_amp = np.ones(len(ii))
            rel_amp = rel_amp / rel_amp.max()
            rel_amp = np.maximum(rel_amp, 0.2)          # keep weak multiplets in the pattern
            deltas = [MP[i]["d"] - MP[i]["base"] - Dk[k] for i in ii]
            best = (np.inf, Dk[k], 0.0, None)
            for Dc in np.arange(-DMAX[k], DMAX[k] + 1e-12, step):
                col = np.zeros(len(xs))
                for i, dl, ra in zip(ii, deltas, rel_amp):
                    col += ra * mcol(xs, MP[i]["p"], MP[i]["it"], MP[i]["base"] + Dc + dl, fw0 * MP[i]["wf"], eta)
                den = (ww2 * col * col).sum()
                if den <= 0:
                    continue
                aa = max(0.0, (ww2 * col * r).sum() / den)
                ssr = (ww2 * (r - aa * col) ** 2).sum()
                if ssr < best[0]:
                    best = (ssr, Dc, aa, col)
            if best[3] is not None and best[0] < base_ssr:
                if abs(best[1] - Dk[k]) > 1e-9:
                    moved += 1
                Dk[k] = best[1]
                for i, dl in zip(ii, deltas):
                    MP[i]["d"] = MP[i]["base"] + Dk[k] + dl
                fit[sup] = fit[sup] - own + best[2] * best[3]
            # pH-sensitive compounds: each multiplet may move further on its own
            if DLT[k] > 2.0 / sf:
                for i in ii:
                    m = MP[i]
                    lo_ = np.searchsorted(x, m["p"].min() + m["base"] + Dk[k] - DLT[k] - HALF_WIN / 2)
                    hi_ = np.searchsorted(x, m["p"].max() + m["base"] + Dk[k] + DLT[k] + HALF_WIN / 2)
                    if hi_ - lo_ < 5:
                        continue
                    xs2 = x[lo_:hi_]
                    w2 = wts[lo_:hi_] ** 2
                    a_i = max(am[i], 0.0)
                    own2 = a_i * mcol(xs2, m["p"], m["it"], m["d"], fw0 * m["wf"], eta)
                    r2 = yv[lo_:hi_] - fit[lo_:hi_] + own2
                    bb = (np.inf, m["d"], None, 0.0)
                    for dd in np.arange(-DLT[k], DLT[k] + 1e-12, step):
                        col = mcol(xs2, m["p"], m["it"], m["base"] + Dk[k] + dd, fw0 * m["wf"], eta)
                        den = (w2 * col * col).sum()
                        if den <= 0:
                            continue
                        aa = max(0.0, (w2 * col * r2).sum() / den)
                        # mild preference for the library position (breaks ties between neighbours)
                        ss = (w2 * (r2 - aa * col) ** 2).sum() * (1.0 + S["SHIFT_PRIOR"] * (dd * sf / 10.0) ** 2)
                        if ss < bb[0]:
                            bb = (ss, m["base"] + Dk[k] + dd, col, aa)
                    if bb[2] is not None:
                        m["d"] = bb[1]
                        fit[lo_:hi_] = fit[lo_:hi_] - own2 + bb[3] * bb[2]
        for m in MP:
            k = m["k"]
            m["dlo"] = m["base"] + Dk[k] - DLT[k]
            m["dhi"] = m["base"] + Dk[k] + DLT[k]
            m["d"] = float(np.clip(m["d"], m["dlo"], m["dhi"]))
        return moved

    # ---------- iterate ----------
    M = mult_columns()
    U = unk_columns()
    coef, fit, A = global_solve(M, U)
    log(f"   start: chi2/n={chi(fit):.2f}  ({time.time()-t0:.0f}s)")
    for it in range(S["N_ITER"]):
        if it == 1 and S["ETA_SEARCH"]:
            best = (chi(fit), eta)
            eta_keep = eta
            for et in S["ETA_SEARCH"]:
                eta = et
                _, f2, _ = global_solve(mult_columns(), U)
                if chi(f2) < best[0]:
                    best = (chi(f2), et)
            eta = best[1] if best[1] is not None else eta_keep
            M = mult_columns()
            coef, fit, A = global_solve(M, U)
            log(f"   lineshape eta -> {eta:.2f}  chi2/n={chi(fit):.2f}")
        moved = compound_stage(coef, fit, M)
        M = mult_columns()
        coef, fit, A = global_solve(M, U)
        c1 = chi(fit)
        nroi = roi_refine(coef, M, U, fit)
        M = mult_columns()
        coef, fit, A = global_solve(M, U)
        log(f"   iter {it+1}: {moved} compounds re-positioned (chi2/n={c1:.2f}); {nroi} ROIs refined -> chi2/n={chi(fit):.2f}  ({time.time()-t0:.0f}s)")

    # ---------- unassigned singlets: only where no fitted library line is nearby ----------
    am = coef[:nM]
    occupied = []
    for i, m in enumerate(MP):
        if am[i] * M[:, i].max() > S["UNKNOWN_BLOCK_SNR"] * sigma:
            occupied.extend((m["p"][m["it"] > 0.05 * m["it"].max()] + m["d"]).tolist())
    occupied = np.sort(np.array(occupied)) if occupied else np.array([])
    r = yv - fit
    pk, _ = find_peaks(r, height=S["UNKNOWN_MIN_SNR"] * sigma, distance=max(1, int(1.5 * fw0 / grid.dx)))
    pk = sorted([p for p in pk if r[p] > S["UNKNOWN_MIN_FRACTION"] * yv[p]], key=lambda p: -r[p])
    for p in pk:
        if len(unknown) >= S["MAX_UNKNOWNS"]:
            break
        if occupied.size:
            j = np.searchsorted(occupied, x[p])
            near = min(abs(occupied[max(j - 1, 0)] - x[p]), abs(occupied[min(j, occupied.size - 1)] - x[p]))
            if near < S["UNKNOWN_EXCLUSION_HZ"] / sf:
                continue
        if unknown and min(abs(u["pos"] - x[p]) for u in unknown) < 1.5 * fw0:
            continue
        unknown.append(dict(pos=float(x[p]), wf=1.2))
    U = unk_columns()
    coef, fit, A = global_solve(M, U)
    log(f"   +{len(unknown)} unassigned singlets -> chi2/n={chi(fit):.2f}")
    roi_refine(coef, M, U, fit)
    M = mult_columns()
    U = unk_columns()
    coef, fit, A = global_solve(M, U)
    # local amplitudes consistent with the final positions (no further movement)
    keep_nfev = S["ROI_MAX_NFEV"]; S["ROI_MAX_NFEV"] = 0
    roi_refine(coef, M, U, fit)
    S["ROI_MAX_NFEV"] = keep_nfev
    log(f"   final chi2/n={chi(fit):.2f}  ({time.time()-t0:.0f}s)")

    am = coef[:nM]
    ua = coef[nM:nM + len(unknown)]
    base_fit = B @ coef[nM + len(unknown):]
    comp_fit = np.zeros((len(x), ncomp))
    amp = np.zeros(ncomp)
    se = np.full(ncomp, np.nan)
    try:
        actv = np.where(coef[:nM + len(unknown)] > 0)[0]
        Aact = np.hstack([A[:, actv], B]) * wts[:, None]
        s2 = chi(fit) * len(yv) / max(len(yv) - Aact.shape[1], 1)
        cov_diag = np.clip(np.diag(np.linalg.pinv(Aact.T @ Aact)) * s2, 0, None)
        se_m = np.full(nM, np.nan)
        for jj, kk in enumerate(actv):
            if kk < nM:
                se_m[kk] = np.sqrt(cov_diag[jj])
    except Exception:
        se_m = np.full(nM, np.nan)
    for k, ii in MPk.items():
        if not ii:
            continue
        comp_fit[:, k] = M[:, ii] @ am[ii]
        npr = np.array([MP[i]["nprot"] for i in ii])
        amp[k] = float((am[ii] * npr).sum() / npr.sum())
        if np.isfinite(se_m[ii]).any():
            se[k] = float(np.sqrt(np.nansum((se_m[ii] * npr) ** 2)) / npr.sum())
    for i, m in enumerate(MP):
        m["global_amp"] = float(am[i])
        m["unit_peak"] = float(M[:, i].max())
        m["apex_ppm"] = float(x[int(np.argmax(M[:, i]))])
        m["pinned"] = bool(abs(Dk[m["k"]]) >= DMAX[m["k"]] - 0.3 / sf) and lib[m["k"]]["source"] != "polymer"
        m["compound_shift_Hz"] = float(Dk[m["k"]] * sf)
    U = unk_columns()
    return dict(tsp=tsp, sigma=sigma, fw0=fw0, eta=eta, g_off=g_off, ppm=ppm, grid=grid, yv=yv, wts=wts,
                fit=fit, base_fit=base_fit, amp=amp, se=se, MP=MP, comp_fit=comp_fit, unknown=unknown, uamp=ua,
                unk_fit=U * ua[None, :] if len(unknown) else np.zeros((len(x), 0)),
                chi2n=chi(fit), sf=sf, mcols=None)


# =============================================================================
# 6. Grading, peak lists
# =============================================================================
def _excluded_fraction(m):
    """Fraction of a multiplet's intensity whose lines sit in an ignored region or outside the fitted range."""
    pos = np.asarray(m["p"], float) + float(m["d"])
    it = np.asarray(m["it"], float)
    lo, hi = S["FIT_RANGE"]
    out = (pos < lo) | (pos > hi)
    for a_, b_ in S["EXCLUDE"]:
        out |= (pos > a_) & (pos < b_)
    return float(it[out].sum() / it.sum()) if it.sum() > 0 else 1.0


def grade(res, lib, sf, cf):
    """cf converts a per-proton area to mM in the medium. Per-multiplet quality
    measures come from the local (ROI) fits, where each multiplet has its own
    amplitude, so a library intensity error elsewhere in the molecule cannot
    bias a clean reporter."""
    x, yv, fit = res["grid"].x, res["yv"], res["fit"]
    sig = res["sigma"]
    rows, prow = [], []
    for k, c in enumerate(lib):
        a = res["amp"][k]
        mps = [m for m in res["MP"] if m["k"] == k]
        if mps and all(m["unit_peak"] <= 0 or _excluded_fraction(m) >= 0.9 for m in mps):
            # every multiplet of this compound lies in an excluded region: nothing to fit
            rows.append(dict(metabolite=c["name"], cls=c["cls"], tier=NOTMEAS, conc_mM=np.nan, se_mM=np.nan,
                             whole_compound_fit_mM=np.nan, method="all multiplets inside ignored regions",
                             LOD_mM=np.nan, upper_bound_mM=np.nan, snr=0.0, best_dominance=np.nan,
                             best_misfit=np.nan, n_reporters=0, reporter_ratio=np.nan,
                             source=c["source"], provenance=c["provenance"]))
            continue
        unit_peak = max(m["unit_peak"] for m in mps) if mps else 0
        snr = a * unit_peak / sig if a > 0 else 0.0
        clean, usable = [], []
        best_dom, best_mis = 0.0, np.nan
        for m in mps:
            share = m["nprot"] / c["nprot"]
            dom, mis, h = m.get("loc_dom", 0.0), m.get("loc_mis", np.inf), m.get("loc_h", 0.0)
            la = m["local_amp"]
            pinned = m.get("pinned", False)
            ok_amp = np.isfinite(la) and la > 0
            is_clean = (ok_amp and not pinned and (share >= 0.1 or m["nprot"] >= 0.3) and dom >= S["Q_MIN_DOMINANCE"] and mis <= S["Q_MAX_MISFIT"]
                        and h / sig >= S["LOQ_SNR"])
            is_usable = (ok_amp and not pinned and (share >= 0.1 or m["nprot"] >= 0.3) and dom >= S["O_MIN_DOMINANCE"] and mis <= S["O_MAX_MISFIT"]
                         and h / sig >= S["LOQ_SNR"])
            if is_clean:
                clean.append(m)
            if is_usable:
                usable.append(m)
            if (share >= 0.1 or m["nprot"] >= 0.3) and ok_amp and dom > best_dom:
                best_dom, best_mis = dom, mis
            keep = m["it"] > 0.03 * m["it"].max()
            fwp = res["fw0"] * m["wf"]
            prow.append(dict(
                metabolite=c["name"], multiplet=m["j"] + 1,
                centre_ppm=round(float((m["p"] * m["it"]).sum() / m["it"].sum() + m["d"]), 4),
                from_ppm=round(float(m["p"].min() + m["d"]), 4), to_ppm=round(float(m["p"].max() + m["d"]), 4),
                n_protons=round(m["nprot"], 2), n_lines=int(keep.sum()),
                line_ppm=";".join(f"{v:.4f}" for v in m["p"][keep] + m["d"]),
                line_rel_intensity=";".join(f"{v:.2f}" for v in m["it"][keep] / m["it"].max()),
                shift_from_library_Hz=round((m["d"] - m["base"]) * sf, 2), compound_shift_at_bound=pinned,
                linewidth_Hz=round(fwp * sf, 2),
                fitted_height=round(h, 1), height_snr=round(h / sig, 1),
                dominance=round(dom, 3), misfit=round(mis, 3) if np.isfinite(mis) else np.nan,
                conc_this_multiplet_mM=round(la * cf, 5) if np.isfinite(la) else np.nan,
                se_this_multiplet_mM=round(m["local_se"] * cf, 5) if np.isfinite(m["local_se"]) else np.nan,
                used_as_reporter=bool(is_clean or (not clean and is_usable))))
        # a common linewidth for the whole molecule, so one badly fitted narrow multiplet
        # cannot produce a too-small upper bound
        wfs = [m["wf"] for m in mps if m.get("loc_h", 0) > S["LOQ_SNR"] * sig]
        wf_med = float(np.median(wfs)) if wfs else 1.0
        lod = S["LOD_SNR"] * sig / unit_peak * cf if unit_peak > 0 else np.nan
        cross = np.nan
        if clean:
            use = clean
        elif usable:
            use = usable
        else:
            use = []
        if use:
            key = "local_amp" if clean else "global_amp"     # shared reporters: use the compound-tied global fit
            vals = [m[key] for m in use]
            if len(vals) >= 3:
                amp_use = float(np.median(vals))             # robust to one contaminated multiplet
                how = "median of"
            elif len(vals) == 2:
                amp_use = float(use[int(np.argmax([m["loc_dom"] for m in use]))][key])
                how = "cleanest of"
            else:
                amp_use = float(vals[0])
                how = ""
            conc = amp_use * cf
            wsum = sum(m["nprot"] for m in use)
            se = np.sqrt(sum((m["local_se"] * m["nprot"]) ** 2 for m in use)) / wsum * cf
            cross = max(vals) / min(vals) if len(vals) > 1 and min(vals) > 0 else np.nan
            method = (how + " " if how else "") + ("clean" if clean else "shared") + \
                " reporter multiplet(s) " + ",".join(str(m["j"] + 1) for m in use)
        else:
            conc = a * cf if a > 0 else 0.0
            se = res["se"][k] * cf if np.isfinite(res["se"][k]) else np.nan
            method = "whole-compound fit"
        # ---- consistency with the data: a multiplet predicted clearly above the noise
        # must actually be there; the data also give an upper bound for every compound
        amp_rep = conc / cf if cf > 0 else 0.0
        xg, base = res["grid"].x, res["base_fit"]
        contradicted, ubs = [], []
        for m in mps:
            if m["unit_peak"] <= 0 or (m["nprot"] < 0.3 and m["nprot"] / c["nprot"] < 0.1):
                continue
            fwp = res["fw0"] * m["wf"]
            half = 1.5 * fwp + (S["MULTIPLET_SHIFT_HZ"] + 1.0) / sf
            j0 = np.searchsorted(xg, m["apex_ppm"] - half)
            j1 = np.searchsorted(xg, m["apex_ppm"] + half)
            seg = xg[max(j0 - 1, 0):min(j1 + 1, len(xg))]
            if seg.size < 2 or np.any(np.diff(seg) > 2 * res["grid"].dx):
                continue          # window touches an excluded region (water, PEO core)
            if j1 <= j0:
                continue          # apex in an excluded region (water, PEO core)
            hobs = max(float((yv[j0:j1] - base[j0:j1]).max()), 0.0)
            up_common = m["unit_peak"] * m["wf"] / wf_med          # unit peak at the common linewidth
            ubs.append(hobs / up_common)
            hpred = amp_rep * m["unit_peak"]
            if hpred > S["CONTRADICT_MIN_SNR"] * sig and hobs < S["CONTRADICT_FRACTION"] * hpred:
                contradicted.append(m["j"] + 1)
        # the most constraining multiplet, but not a single outlier when several are available
        ub = (np.percentile(ubs, 20) if len(ubs) >= 4 else min(ubs)) * cf if ubs else np.nan
        if np.isfinite(ub) and ub < lod:
            tier, conc = "Not detected", 0.0
        elif contradicted:
            tier, conc = "Upper bound only", ub
            method = f"multiplet(s) {','.join(map(str, contradicted))} predicted but absent -> data upper bound"
        elif np.isfinite(ub) and conc > S["MAX_CONC_OVER_BOUND"] * ub:
            tier, conc = "Upper bound only", ub
            method = "fitted value exceeds what the data allow at another multiplet -> data upper bound"
        elif clean and (not np.isfinite(cross) or cross <= S["MAX_REPORTER_RATIO"]) and \
                (len(clean) >= 2 or clean[0]["nprot"] / c["nprot"] >= S["SINGLE_REPORTER_MIN_SHARE"]):
            tier = "Quantified"
        elif use:
            tier = "Overlapped (semi-quantitative)"
        else:
            wc = a * cf if a > 0 else 0.0
            wc_se = res["se"][k] * cf if np.isfinite(res["se"][k]) else np.nan
            n_fit = sum(1 for m in mps if m.get("loc_h", 0) > S["LOQ_SNR"] * sig)
            if (np.isfinite(ub) and np.isfinite(wc_se) and wc > lod and wc <= ub
                    and n_fit >= S["DECONV_MIN_MULTIPLETS"]
                    and wc_se <= S["DECONV_MAX_REL_SE"] * wc):
                tier, conc, se = DECONV, wc, wc_se
                method = (f"whole-compound deconvolution over {n_fit} fitted multiplets, "
                          f"no resolved reporter (data upper bound {ub:.3g} mM)")
            else:
                tier, conc = "Upper bound only", ub if np.isfinite(ub) else conc
                method = "no usable reporter -> data upper bound"
        rows.append(dict(metabolite=c["name"], cls=c["cls"], tier=tier, conc_mM=conc, se_mM=se,
                         whole_compound_fit_mM=a * cf if a > 0 else 0.0, method=method,
                         LOD_mM=lod, upper_bound_mM=ub, snr=round(snr, 1), best_dominance=round(best_dom, 3),
                         best_misfit=round(best_mis, 3) if np.isfinite(best_mis) else np.nan,
                         n_reporters=len(use), reporter_ratio=round(cross, 3) if np.isfinite(cross) else np.nan,
                         source=c["source"], provenance=c["provenance"]))
    return pd.DataFrame(rows), pd.DataFrame(prow)


def observed_peaks(res, lib, sample, sf):
    x = res["grid"].x
    ys = res["yv"] - res["base_fit"]
    pk, _ = find_peaks(ys, height=10 * res["sigma"], prominence=5 * res["sigma"], distance=max(1, int(0.8 / sf / res["grid"].dx)))
    names = [c["name"] for c in lib] + [f"unassigned singlet {u['pos']:.4f}" for u in res["unknown"]]
    contrib = np.hstack([res["comp_fit"], res["unk_fit"]])
    rows = []
    for p in pk:
        v = contrib[p].clip(min=0)
        tot = v.sum()
        order = np.argsort(-v)[:4]
        assign = [f"{names[j]} ({v[j]/tot:.0%})" for j in order if tot > 0 and v[j] / tot >= 0.10]
        top = names[order[0]] if tot > 0 else ""
        rows.append(dict(sample=sample, ppm=round(float(x[p]), 4), height=float(ys[p]), snr=round(float(ys[p] / res["sigma"]), 1),
                         fitted_height=float(res["fit"][p] - res["base_fit"][p]),
                         assignment="; ".join(assign) if assign else "unassigned",
                         top_contributor=top, top_fraction=round(float(v[order[0]] / tot), 3) if tot > 0 else 0.0,
                         library_assigned=bool(tot > 0 and not top.startswith("unassigned"))))
    return pd.DataFrame(rows)


# =============================================================================
# 7. Plots
# =============================================================================
PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#17becf",
           "#bcbd22", "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#e6550d"]
TAG = {"Quantified": "Q", "Overlapped (semi-quantitative)": "O", DECONV: "D",
       "Upper bound only": "U", "Not detected": "", NOTMEAS: ""}


def gapped(x, *ys, dx=None):
    """Insert NaN where the ppm axis has a gap (excluded regions) so lines are not drawn across it."""
    if len(x) < 2:
        return (x,) + ys
    dx = dx or np.median(np.diff(x))
    cut = np.where(np.diff(x) > 3 * dx)[0] + 1
    if not len(cut):
        return (x,) + ys
    xo = np.insert(x.astype(float), cut, np.nan)
    return (xo,) + tuple(np.insert(np.asarray(y, float), cut, np.nan) for y in ys)


def overlay_pdf(res, lib, grades, sample, path):
    x = res["grid"].x
    names = [c["name"] for c in lib]
    tiers = dict(zip(grades.metabolite, grades.tier))
    conc = dict(zip(grades.metabolite, grades.conc_mM))
    with PdfPages(path) as pdf:
        fig, ax = plt.subplots(2, 1, figsize=(16, 9), gridspec_kw=dict(height_ratios=[3, 1]), sharex=True)
        gx, gy, gf, gr_ = gapped(x, res["yv"], res["fit"], res["yv"] - res["fit"])
        ax[0].plot(gx, gy, color="k", lw=0.5, label="observed")
        ax[0].plot(gx, gf, color="red", lw=0.5, alpha=0.8, label="simulated (library fit + baseline)")
        ax[0].set_ylim(-0.02 * np.percentile(res["yv"], 99.9), np.percentile(res["yv"], 99.8) * 1.1)
        ax[0].legend(loc="upper left")
        expl = (res["comp_fit"].sum() / max((res["yv"] - res["base_fit"]).clip(min=0).sum(), 1e-9))
        ax[0].set_title(f"Sample {sample}: observed vs simulated spectrum   chi2/n={res['chi2n']:.1f}   "
                        f"library explains {expl:.0%} of signal area   TSP linewidth {res['fw0']*res['sf']:.2f} Hz")
        ax[1].plot(gx, gr_, color="grey", lw=0.5)
        ax[1].set_ylabel("residual")
        ax[1].set_xlim(S["FIT_RANGE"][1], S["FIT_RANGE"][0])
        ax[1].set_xlabel("ppm (TSP = 0)")
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        for lo, hi in S["OVERLAY_REGIONS"]:
            w = (x >= lo) & (x <= hi)
            if w.sum() < 10:
                continue
            xs = x[w]
            b = res["base_fit"][w]
            fig, ax = plt.subplots(3, 1, figsize=(16, 11), gridspec_kw=dict(height_ratios=[3, 3, 1]), sharex=True)
            gx, gy, gf, gb, grr = gapped(xs, res["yv"][w], res["fit"][w], b, res["yv"][w] - res["fit"][w])
            ax[0].plot(gx, gy, "k", lw=0.8, label="observed")
            ax[0].plot(gx, gf, "r", lw=0.8, alpha=0.85, label="simulated total")
            ax[0].plot(gx, gb, color="tab:brown", lw=0.8, ls=":", label="baseline")
            top = max(res["yv"][w].max(), res["fit"][w].max())
            ax[0].set_ylim(min(0, b.min()) - 0.03 * top, top * 1.08)
            ax[0].legend(loc="upper right", fontsize=8)
            contrib = res["comp_fit"][w]
            order = np.argsort(-contrib.max(0))
            ci = 0
            for k in order:
                if contrib[:, k].max() < max(8 * res["sigma"], 0.003 * top) or ci >= 24:
                    break
                col = PALETTE[ci % len(PALETTE)]; ci += 1
                _, cb = gapped(xs, b + contrib[:, k])
                ax[1].fill_between(gx, gb, cb, color=col, alpha=0.40, lw=0)
                ax[1].plot(gx, cb, color=col, lw=0.7)
                j = np.argmax(contrib[:, k])
                lab = f"{names[k]} [{TAG.get(tiers.get(names[k]), '')}]"
                ax[1].annotate(lab, (xs[j], b[j] + contrib[j, k]), fontsize=7, color=col, rotation=55,
                               ha="left", va="bottom")
            if res["unk_fit"].shape[1]:
                uf = res["unk_fit"][w]
                for j in range(uf.shape[1]):
                    if uf[:, j].max() > 8 * res["sigma"]:
                        ax[1].plot(gx, gapped(xs, b + uf[:, j])[1], color="grey", lw=0.8, ls="--")
            ax[1].plot(gx, gy, color="k", lw=0.4, alpha=0.4)
            ax[1].set_ylim(ax[0].get_ylim())
            ax[1].set_ylabel("simulated metabolite spectra\n(grey dashed = unassigned peaks)", fontsize=8)
            ax[2].plot(gx, grr, color="grey", lw=0.6)
            ax[2].axhline(0, color="k", lw=0.4)
            ax[2].set_ylabel("residual")
            ax[2].set_xlim(hi, lo)
            ax[2].set_xlabel("ppm (TSP = 0)")
            ax[0].set_title(f"Sample {sample}   {lo:.2f}-{hi:.2f} ppm      [Q] quantified  [O] overlapped/semi-quantitative  [U] upper bound only")
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def gallery_pdf(results, lib, G, path):
    names = [c["name"] for c in lib]
    samples = list(results.keys())
    show = [m for m in G.sort_values("conc_mM", ascending=False).metabolite.unique()
            if (G[G.metabolite == m].tier != "Not detected").any() and not m.startswith("Pluronic")]
    per_page = 5
    with PdfPages(path) as pdf:
        for i0 in range(0, len(show), per_page):
            ch = show[i0:i0 + per_page]
            fig, axes = plt.subplots(len(ch), len(samples), figsize=(4.2 * len(samples), 2.8 * len(ch)), squeeze=False)
            for r_, met in enumerate(ch):
                k = names.index(met)
                for s_, smp in enumerate(samples):
                    res = results[smp]
                    ax = axes[r_, s_]
                    gr = G[(G["sample"] == smp) & (G.metabolite == met)].iloc[0]
                    mps = [m for m in res["MP"] if m["k"] == k]
                    # show the multiplet with the largest fitted contribution
                    best = max(mps, key=lambda m: m["nprot"] * (1 if not np.isfinite(m["local_amp"]) else 1 + m["local_amp"]))
                    cen = float((best["p"] * best["it"]).sum() / best["it"].sum() + best["d"])
                    span = max(0.02, (best["p"].max() - best["p"].min()) / 2 + 0.012)
                    x = res["grid"].x
                    w = (x > cen - span) & (x < cen + span)
                    if w.sum() < 5:
                        ax.axis("off"); continue
                    b = res["base_fit"][w]
                    gx, gy, gf, gb, gc = gapped(x[w], res["yv"][w], res["fit"][w], b, b + res["comp_fit"][w, k])
                    ax.plot(gx, gy, "k", lw=0.7)
                    ax.plot(gx, gf, "r", lw=0.7, alpha=0.8)
                    ax.fill_between(gx, gb, gc, color="tab:blue", alpha=0.5, lw=0)
                    ax.set_xlim(x[w].max(), x[w].min())
                    ax.set_title(f"{met} | {smp}\n{gr.tier.split(' (')[0]}  {gr.conc_mM:.3g} mM", fontsize=8)
                    ax.tick_params(labelsize=6)
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def gallery_from_csv(out, G, P, samples, path, max_metabolites=200):
    """Observed vs fitted around each metabolite's reporter multiplet, all samples side by
    side. Built from the per-sample simulated-spectrum CSVs, so it needs no fit in memory."""
    show = [m for m in G.groupby("metabolite").conc_mM.mean().sort_values(ascending=False).index
            if (G[G.metabolite == m].tier != "Not detected").any() and not m.startswith("Pluronic")][:max_metabolites]
    win = {}
    for smp in samples:
        f = os.path.join(out, f"simulated_spectrum_sample_{smp}.csv.gz")
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f)
        ps = P[P["sample"].astype(str) == str(smp)]
        for m in show:
            sub = ps[ps.metabolite == m]
            if not len(sub):
                continue
            row = sub[sub.used_as_reporter] if sub.used_as_reporter.any() else sub
            row = row.sort_values("n_protons", ascending=False).iloc[0]
            half = max(0.018, (row.to_ppm - row.from_ppm) / 2 + 0.012)
            w = (d.ppm > row.centre_ppm - half) & (d.ppm < row.centre_ppm + half)
            if w.sum() < 5:
                continue
            dd = d[w]
            win[(m, smp)] = (dd.ppm.values, dd.observed.values, dd.simulated_total.values,
                             dd.baseline.values, (dd[m].values if m in dd else np.zeros(w.sum())))
    ncol = len(samples)
    with PdfPages(path) as pdf:
        per = max(1, min(6, int(30 / max(ncol, 1)) + 1))
        for i0 in range(0, len(show), per):
            ch = show[i0:i0 + per]
            fig, axes = plt.subplots(len(ch), ncol, figsize=(max(2.0, 16 / ncol) * ncol, 2.4 * len(ch)), squeeze=False)
            for r, m in enumerate(ch):
                for c, smp in enumerate(samples):
                    ax = axes[r, c]
                    if (m, smp) not in win:
                        ax.axis("off")
                        continue
                    x, obs, fit, base, own = win[(m, smp)]
                    gx, go, gf, gb, gc = gapped(x, obs, fit, base, base + own)
                    ax.plot(gx, go, "k", lw=0.7)
                    ax.plot(gx, gf, "r", lw=0.7, alpha=0.8)
                    ax.fill_between(gx, gb, gc, color="tab:blue", alpha=0.5, lw=0)
                    ax.set_xlim(x.max(), x.min())
                    g = G[(G["sample"].astype(str) == str(smp)) & (G.metabolite == m)]
                    lab = f"{m} | {smp}"
                    if len(g):
                        lab += f"\n{g.iloc[0].tier.split(' (')[0]}  {g.iloc[0].conc_mM:.3g} mM"
                    ax.set_title(lab, fontsize=7.5)
                    ax.set_yticks([]); ax.tick_params(labelsize=6)
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)


def profiles_pdf(conc, tier, path, title="Concentration across samples"):
    """Small multiples: concentration vs sample for every metabolite seen in the run."""
    samples = list(conc.columns)
    x = np.arange(len(samples))
    ok = tier.isin(["Quantified", "Overlapped (semi-quantitative)", DECONV])
    keep = [m for m in conc.index if ok.loc[m].sum() >= max(2, 0.25 * len(samples))]
    keep = sorted(keep, key=lambda m: -np.nanmean(np.where(ok.loc[m].values,
                                                           conc.loc[m].values.astype(float), np.nan)))
    per, ncol = 12, 4
    with PdfPages(path) as pdf:
        for i0 in range(0, len(keep), per):
            ch = keep[i0:i0 + per]
            nrow = int(np.ceil(len(ch) / ncol))
            fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 2.6 * nrow), squeeze=False)
            for ax in axes.ravel():
                ax.axis("off")
            for j, m in enumerate(ch):
                ax = axes[j // ncol, j % ncol]
                ax.axis("on")
                v = conc.loc[m].values.astype(float)
                t = tier.loc[m].values
                for style, mk, fc, lab in [("Quantified", "o", None, "quantified"),
                                           ("Overlapped (semi-quantitative)", "s", "none", "overlapped"),
                                           (DECONV, "D", "none", "deconvolution estimate"),
                                           ("Upper bound only", "v", "none", "upper bound")]:
                    sel = t == style
                    if sel.any():
                        ax.scatter(x[sel], v[sel], marker=mk, s=18, facecolors=fc,
                                   color=("#2563eb" if style == "Quantified" else
                                          "#d97706" if style.startswith("Overlapped") else "#8a8a8a"),
                                   label=lab, zorder=3)
                good = t != "Not detected"
                ax.plot(x[good], v[good], color="#c9c9c7", lw=0.8, zorder=1)
                ax.set_title(m, fontsize=9)
                ax.set_ylabel("mM", fontsize=8)
                ax.set_ylim(bottom=0)
                step = max(1, len(samples) // 8)
                ax.set_xticks(x[::step])
                ax.set_xticklabels([samples[k] for k in range(0, len(samples), step)], fontsize=6, rotation=90)
                ax.tick_params(labelsize=7)
                for sp in ("top", "right"):
                    ax.spines[sp].set_visible(False)
                if j == 0:
                    ax.legend(frameon=False, fontsize=6.5, loc="upper right")
            fig.suptitle(title, x=0.02, ha="left", fontsize=11)
            fig.tight_layout(rect=[0, 0, 1, 0.97])
            pdf.savefig(fig)
            plt.close(fig)


# =============================================================================
# 8. Main
# =============================================================================
def process_sample(job):
    """Fit, grade and plot one spectrum (runs in a worker process). Results are written
    to <out>/_per_sample as they are produced, so a crash or a stop loses nothing: the
    next run re-uses them and only fits what is missing."""
    smp, d, libdir, out, settings = job
    S.update(settings)
    per = os.path.join(out, "_per_sample")
    done = [os.path.join(per, f"{smp}_{k}.csv") for k in ("grades", "peaks", "obs", "qc")]
    if S["RESUME"] and all(os.path.exists(f) for f in done):
        gr, pk, ob, qdf = (pd.read_csv(f) for f in done)
        return gr, pk, ob, qdf.iloc[0].to_dict(), None, [f"Sample {smp}: re-used earlier result"]
    try:
        return _fit_one(smp, d, libdir, out)
    except Exception as exc:
        import traceback
        return None, None, None, dict(sample=smp, error=f"{type(exc).__name__}: {exc}"), None, \
            [f"Sample {smp} FAILED: {type(exc).__name__}: {exc}", traceback.format_exc()]


def sample_source(smp):
    """(EXPNO folder, procno) for a sample: SAMPLE_DIRS if given, else DATA_DIR/<sample>, PROCNO."""
    sd = S.get("SAMPLE_DIRS") or {}
    if smp in sd:
        return sd[smp]
    return os.path.join(S["DATA_DIR"], smp), S["PROCNO"]


def detect_pluronic(ppm, y):
    """True if the spectrum has the tall PEO line of Pluronic F-68 near 3.71 ppm."""
    m = (ppm > 3.69) & (ppm < 3.74)
    return bool(m.any() and y[m].max() > S["PEO_DETECT_SNR"] * noise_sigma(ppm, y))


def _fit_one(smp, d, libdir, out):
    lines = []
    log = lambda msg: lines.append(msg)
    expdir, procno = sample_source(smp)
    ppm, y, sf, title = read_spectrum(expdir, procno)
    pluronic = detect_pluronic(ppm, y) if S["PLURONIC"] == "auto" else S["PLURONIC"] == "yes"
    lib = prepare_library(build_library(sf, libdir, os.path.join(out, "_cache")), sf, pluronic)
    log(f"Sample {smp}  (title: {title}, dilution x{d}, {sf:.2f} MHz"
        f"{', Pluronic terms added' if pluronic else ''})")
    t0 = time.time()
    res = fit_spectrum(ppm, y, sf, lib, log)
    cf = S["TSP_MM_IN_TUBE"] * S["TSP_PROTONS"] / res["tsp"]["area"] * d
    gr, pk = grade(res, lib, sf, cf)
    gr.insert(0, "sample", smp); pk.insert(0, "sample", smp)
    ob = observed_peaks(res, lib, smp, sf)
    sig_area = max((res["yv"] - res["base_fit"]).clip(min=0).sum(), 1e-9)
    q = dict(sample=smp, title=title, SF_MHz=sf, dilution=d, TSP_area=res["tsp"]["area"],
             TSP_fwhm_Hz=res["tsp"]["fwhm_hz"], TSP_eta=res["tsp"]["eta"], TSP_fit_rel_rms=res["tsp"]["rel_rms"],
             noise_sigma=res["sigma"], TSP_SNR=res["tsp"]["height"] / res["sigma"],
             GISSMO_offset_Hz=res["g_off"] * sf, lineshape_eta=res["eta"], chi2_per_point=res["chi2n"],
             unassigned_singlets=len(res["unknown"]),
             signal_explained_by_library=res["comp_fit"].sum() / sig_area,
             signal_in_unassigned_singlets=res["unk_fit"].sum() / sig_area if res["unk_fit"].size else 0.0,
             abs_residual_fraction=np.abs(res["yv"] - res["fit"]).sum() / sig_area,
             fit_time_s=round(time.time() - t0, 1))
    if S["OVERLAYS"]:
        overlay_pdf(res, lib, gr, smp, os.path.join(out, f"overlay_sample_{smp}.pdf"))
    per = os.path.join(out, "_per_sample")
    os.makedirs(per, exist_ok=True)
    x = res["grid"].x
    cols = {"ppm": x, "observed": res["yv"], "simulated_total": res["fit"], "baseline": res["base_fit"]}
    for k, c in enumerate(lib):
        if res["amp"][k] > 0:
            cols[c["name"]] = res["comp_fit"][:, k]
    if res["unk_fit"].size:
        cols["unassigned_singlets"] = res["unk_fit"].sum(1)
    pd.DataFrame(cols).to_csv(os.path.join(out, f"simulated_spectrum_sample_{smp}.csv.gz"), index=False, float_format="%.6g")
    np.savez_compressed(os.path.join(per, f"{smp}_fit.npz"), cf=cf, sigma=res["sigma"], fw0=res["fw0"],
                        x=res["grid"].x.astype(np.float32), yv=res["yv"].astype(np.float32),
                        fit=res["fit"].astype(np.float32), base=res["base_fit"].astype(np.float32),
                        dx=res["grid"].dx, sf=sf, amp=res["amp"], se=res["se"], pluronic=pluronic,
                        mp_json=np.array(json.dumps([{k: (v.tolist() if isinstance(v, np.ndarray) else v)
                                                      for k, v in m.items()} for m in res["MP"]],
                                                     default=_json_default)))
    gr.to_csv(os.path.join(per, f"{smp}_grades.csv"), index=False)
    pk.to_csv(os.path.join(per, f"{smp}_peaks.csv"), index=False)
    ob.to_csv(os.path.join(per, f"{smp}_obs.csv"), index=False)
    pd.DataFrame([q]).to_csv(os.path.join(per, f"{smp}_qc.csv"), index=False)
    open(os.path.join(per, f"{smp}_log.txt"), "w").write("\n".join(lines))
    return gr, pk, ob, q, None, lines


TIER_ORDER = {"Quantified": 0, "Overlapped (semi-quantitative)": 1, DECONV: 2,
              "Upper bound only": 3, "Not detected": 4, NOTMEAS: 5}


def aggregate_outputs(out, grades, peaks, obs, qc, lib, samples, log):
    """Build every table and figure from the per-sample results (also used by regrade.py)."""
    G = pd.concat(grades, ignore_index=True)
    P = pd.concat(peaks, ignore_index=True)
    O = pd.concat(obs, ignore_index=True)
    Q = pd.DataFrame(qc)

    # ---------- summary tables ----------
    G = G[~G.metabolite.str.startswith("Pluronic")].copy()
    best_tier = G.groupby("metabolite").tier.agg(lambda t: min(t, key=lambda v: TIER_ORDER[v]))
    conc = G.pivot(index="metabolite", columns="sample", values="conc_mM")
    tier = G.pivot(index="metabolite", columns="sample", values="tier")
    se = G.pivot(index="metabolite", columns="sample", values="se_mM")
    lod = G.pivot(index="metabolite", columns="sample", values="LOD_mM")
    report = conc.copy().astype(object)
    for m in report.index:
        for s_ in report.columns:
            t = tier.loc[m, s_]
            if t == NOTMEAS:
                report.loc[m, s_] = "n/m"
            elif t == "Not detected":
                report.loc[m, s_] = f"<{lod.loc[m, s_]:.2g}"
            elif t == "Upper bound only":
                report.loc[m, s_] = f"<={conc.loc[m, s_]:.2g}"
            elif t == DECONV:
                report.loc[m, s_] = f"~{conc.loc[m, s_]:.3g}"
            else:
                report.loc[m, s_] = round(float(conc.loc[m, s_]), 4)
    summ = pd.DataFrame(dict(best_tier=best_tier))
    summ["class"] = G.groupby("metabolite").cls.first()
    for lab, cond in [("n_quantified", G.tier == "Quantified"), ("n_overlapped", G.tier.str.startswith("Overlapped")),
                      ("n_detected_any", ~G.tier.isin(["Not detected", NOTMEAS]))]:
        summ[lab] = G[cond].groupby("metabolite").size()
    summ = summ.fillna({"n_quantified": 0, "n_overlapped": 0, "n_detected_any": 0})
    summ["mean_conc_mM"] = conc.mean(axis=1)
    summ["mean_rel_SE_%"] = (se / conc.replace(0, np.nan) * 100).mean(axis=1).round(1)
    summ["library_source"] = G.groupby("metabolite").source.first()
    summ["_o"] = summ.best_tier.map(TIER_ORDER)
    summ = summ.sort_values(["_o", "mean_conc_mM"], ascending=[True, False]).drop(columns="_o")
    idx = summ.index
    conc, report, tier, se, lod = conc.loc[idx], report.loc[idx], tier.loc[idx], se.loc[idx], lod.loc[idx]

    readme = pd.DataFrame(dict(item=[
        "Software", "Units", "Reference", "How numbers are produced", "Quantified", "Overlapped (semi-quantitative)",
        DECONV, "Upper bound only", "Not detected", "SE", "Library", "Peak_list_by_metabolite",
        "Observed_peaks_assigned", "QC", "Settings"], text=[
        f"MANC-Q {__version__} (https://github.com/willwas1/MANC-Q), results written {time.strftime('%Y-%m-%d %H:%M')}.",
        "mM in the original medium = tube mM x dilution factor (from the dilution file).",
        f"{S['REFERENCE']} internal standard only: {S['TSP_MM_IN_TUBE']} mM in tube, 9 H, area from a pseudo-Voigt fit incl. 29Si satellites. No calibration against any sample or standard.",
        "Whole spectrum fitted as a sum of simulated metabolite spectra + smooth baseline + free unassigned singlets; each multiplet has a small fitted shift and linewidth. Concentration = mean of the cleanly resolved 'reporter' multiplets (each fitted on its own); if none, the whole-compound fit.",
        f"Clean reporter multiplet(s): SNR>={S['LOQ_SNR']}, this metabolite >={S['Q_MIN_DOMINANCE']:.0%} of the fitted signal under it, |observed-fit| <= {S['Q_MAX_MISFIT']} x its own signal, shift not at bound; either two or more reporters agreeing within x{S['MAX_REPORTER_RATIO']}, or one reporter carrying >={S['SINGLE_REPORTER_MIN_SHARE']:.0%} of the protons; no multiplet of the molecule contradicted by the data.",
        f"SNR>={S['LOQ_SNR']} but the best multiplet is shared (dominance {S['O_MIN_DOMINANCE']:.0%}-{S['Q_MIN_DOMINANCE']:.0%}) or misfit <= {S['O_MAX_MISFIT']}. Number depends on the overlap partners: use for trends.",
        f"No multiplet of this metabolite is resolved enough to quantify on its own, but the whole-compound least-squares fit over >={S['DECONV_MIN_MULTIPLETS']} fitted multiplets is well determined (relative SE <= {S['DECONV_MAX_REL_SE']:.0%}) and sits below the data upper bound. The value is that deconvolution estimate, shown with a ~ in Report_mM. It depends on the library explaining the overlapping signals correctly, so treat it as an estimate, not a measurement.",
        f"Signal is present where the metabolite would appear but it cannot be attributed to it with confidence (no usable reporter, or a multiplet the library predicts >{S['CONTRADICT_MIN_SNR']}x noise is missing). The value is the largest concentration the data allow (upper bound).",
        "Below the detection limit; Report_mM shows <LOD in mM. 'n/m' = not measurable: every peak of the metabolite lies in a region left out of the fit.",
        "Standard error from the least-squares fit only (noise + local misfit). Excludes library, relaxation (T1) and TSP weighing errors.",
        f"{sum(c['source']=='GISSMO' for c in lib)} GISSMO spin systems (BMRB, QM-simulated at the spectrometer frequency) + {sum(c['source']=='SPIN' for c in lib)} other spin system(s) simulated the same way + {sum(c['source']=='CASMDB' for c in lib)} fixed peak lists (CASMDB) + optional Pluronic F-68 terms. See Library sheet.",
        "Every multiplet of every metabolite in every sample: fitted centre, range, line positions and relative intensities, protons, linewidth, shift from library, dominance, misfit, concentration from that multiplet alone, and whether it was used as a reporter.",
        "Peaks picked in the observed spectrum (>10 x noise) with the fitted contributors at each apex (the peak assignment list).",
        "Per-sample referencing, linewidth, noise, fit quality and fraction of signal explained.",
        json.dumps({k: v for k, v in S.items() if k != "OVERLAY_REGIONS"})]))
    libdf = pd.DataFrame([dict(metabolite=c["name"], cls=c["cls"], source=c["source"], n_protons=c["nprot"],
                               titratable=c["titratable"], n_lines=len(c["ppm"]), n_multiplets=len(c["mults"]),
                               provenance=c["provenance"]) for c in lib])
    xl = os.path.join(out, "metabolite_concentrations.xlsx")
    with pd.ExcelWriter(xl, engine="openpyxl") as xw:
        readme.to_excel(xw, sheet_name="README", index=False)
        summ.to_excel(xw, sheet_name="Summary")
        report.to_excel(xw, sheet_name="Report_mM")
        conc.to_excel(xw, sheet_name="Concentration_mM_all")
        se.to_excel(xw, sheet_name="SE_mM")
        tier.to_excel(xw, sheet_name="Tier")
        lod.to_excel(xw, sheet_name="LOD_mM")
        G.to_excel(xw, sheet_name="Long_table", index=False)
        P.to_excel(xw, sheet_name="Peak_list_by_metabolite", index=False)
        O.to_excel(xw, sheet_name="Observed_peaks_assigned", index=False)
        Q.to_excel(xw, sheet_name="QC", index=False)
        libdf.to_excel(xw, sheet_name="Library", index=False)
        for ws in xw.book.worksheets:
            ws.freeze_panes = "B2"
            for col in ws.columns:
                width = min(60, max(10, max(len(str(c.value)) if c.value is not None else 0 for c in col[:50]) + 2))
                ws.column_dimensions[col[0].column_letter].width = width
    if len(samples) <= S["GALLERY_MAX_SAMPLES"]:
        gallery_from_csv(out, G, P, samples, os.path.join(out, "metabolite_gallery.pdf"))
    profiles_pdf(conc, tier, os.path.join(out, "concentration_profiles.pdf"),
                 f"Concentration across samples ({os.path.basename(out)})")
    G.to_csv(os.path.join(out, "concentrations_long.csv"), index=False)
    report.to_csv(os.path.join(out, "concentrations_report_mM.csv"))
    P.to_csv(os.path.join(out, "peak_list_by_metabolite.csv"), index=False)
    O.to_csv(os.path.join(out, "observed_peaks_assigned.csv"), index=False)
    Q.to_csv(os.path.join(out, "qc_per_sample.csv"), index=False)
    log("\nBest tier per metabolite across samples:")
    log(summ.best_tier.value_counts().to_string())
    log(Q.round(3).T.to_string())
    log(f"\nWrote {xl}")
    return G


def run(data_dir, out_dir, ref_mm, dilution=1.0, dilution_file=None, progress=None, cancel=None, **settings):
    """Quantify every spectrum under data_dir (or the single EXPNO folder data_dir) into out_dir.

    progress(event, sample, info): optional callback; events are "run_started", "sample_started",
    "sample_done", "sample_failed", "cancelled", "aggregating", "finished".
    cancel(): optional function returning True to stop; spectra already running are finished first."""
    if ref_mm is None or ref_mm <= 0:
        raise ValueError("the reference (TSP/DSS) concentration in the tube, in mM, is required")
    say = progress or (lambda *a, **k: None)
    stop = cancel or (lambda: False)
    S.clear(); S.update(_copy.deepcopy(DEFAULT_SETTINGS))   # each run starts from the defaults
    S.update(settings)
    S.update(DATA_DIR=os.path.abspath(data_dir) if data_dir else None, OUT_DIR=os.path.abspath(out_dir),
             TSP_MM_IN_TUBE=float(ref_mm), DEFAULT_DILUTION=float(dilution), DILUTION_FILE=dilution_file)
    if S["PLURONIC"] == "ignore":
        S["EXCLUDE"] = list(S["EXCLUDE"]) + [r for r in S["PLURONIC_REGIONS"] if r not in S["EXCLUDE"]]
    if S["SAMPLES"]:
        samples = list(S["SAMPLES"])
    else:
        samples = list_samples(S["DATA_DIR"], S["PROCNO"], fids=S["FID_MODE"] != "never")
        if samples == [""]:                       # data_dir is itself one EXPNO folder
            S["DATA_DIR"], one = os.path.split(S["DATA_DIR"].rstrip("/\\"))
            samples = [one]
    if not samples:
        raise FileNotFoundError(f"no spectra found: expected <EXPNO>/pdata/{S['PROCNO']}/1r under {data_dir}")
    out = S["OUT_DIR"]
    os.makedirs(out, exist_ok=True)
    logf = open(os.path.join(out, "run_log.txt"), "w")

    def log(msg):
        print(msg, flush=True); logf.write(msg + "\n"); logf.flush()

    libdir = S["LIBRARY_DIR"] or os.path.join(HERE, "library")
    dil = read_dilutions(S["DILUTION_FILE"])
    # raw FIDs: process them (automatic phasing) into <out>/_fid_processed unless already given in SAMPLE_DIRS
    sd = dict(S.get("SAMPLE_DIRS") or {})
    if S["FID_MODE"] != "never":
        from . import fidproc
        for smp in samples:
            if smp in sd:
                continue
            ed = os.path.join(S["DATA_DIR"], smp)
            if fidproc.has_fid(ed) and (S["FID_MODE"] == "all" or not fidproc.has_processed(ed, S["PROCNO"])):
                dest = os.path.join(out, "_fid_processed", smp)
                say("fid_processing", smp, dict(k=len(sd) + 1, n=len(samples)))
                r = fidproc.process_fid(ed, lb=S["FID_LB"])
                fidproc.write_processed(r, dest, title=f"EXPNO {smp} (from FID)")
                open(os.path.join(dest, "processing.txt"), "w").write(
                    f"Processed from {ed} by MANC-Q: LB {r['lb']} Hz, SI {r['si']}, automatic phase "
                    f"ph0 {r['ph0']:.2f}, ph1 {r['ph1']:.2f}\n")
                sd[smp] = (dest, 1)
                log(f"EXPNO {smp}: processed from the FID (automatic phase ph0 {r['ph0']:.1f}, ph1 {r['ph1']:.1f}); "
                    f"check its overlay")
    S["SAMPLE_DIRS"] = sd
    os.makedirs(os.path.join(out, "_per_sample"), exist_ok=True)
    json.dump(dict({k: v for k, v in S.items() if k != "OVERLAY_REGIONS"}, MANCQ_VERSION=__version__),
              open(os.path.join(out, "_per_sample", "run_settings.json"), "w"), indent=1, default=str)
    log(f"MANC-Q {__version__}   {time.strftime('%Y-%m-%d %H:%M')}\nData: {S['DATA_DIR']}\n"
        f"Samples ({len(samples)}): {samples}\nReference: {S['REFERENCE']} {S['TSP_MM_IN_TUBE']} mM in tube\n"
        f"Ignored regions (ppm): {', '.join(f'{a:.2f}-{b:.2f}' for a, b in S['EXCLUDE'])}\n"
        f"Pluronic F-68: {S['PLURONIC']}")
    d0, p0 = sample_source(samples[0])
    _, _, sf, _ = read_spectrum(d0, p0)
    log(f"Library at {sf:.3f} MHz ...")
    lib = prepare_library(build_library(sf, libdir, os.path.join(out, "_cache")), sf, S["PLURONIC"] in ("auto", "yes"))
    log(f"  {sum(c['source']=='GISSMO' for c in lib)} GISSMO spin systems, {sum(c['source']=='SPIN' for c in lib)} other spin systems, "
        f"{sum(c['source']=='CASMDB' for c in lib)} peak lists, {sum(c['source']=='polymer' for c in lib)} polymer terms")
    nw = S["N_WORKERS"] or max(1, min(len(samples), (os.cpu_count() or 2) - 1))
    nw = max(1, min(nw, len(samples)))
    jobs = {smp: (smp, dil.get(smp, S["DEFAULT_DILUTION"]), libdir, out, dict(S)) for smp in samples}
    log(f"Fitting {len(samples)} spectra with {nw} worker(s); a few minutes per spectrum ...")
    say("run_started", None, dict(n=len(samples), workers=nw, sf=sf))
    results, cancelled = {}, []

    def finished(smp, output, t0):
        gr, pk, ob, q, res, lines = output
        for ln in lines:
            log(ln)
        results[smp] = output
        if gr is None:
            say("sample_failed", smp, dict(error=q.get("error", "unknown error"), seconds=time.time() - t0))
        else:
            say("sample_done", smp, dict(seconds=time.time() - t0, n_quantified=int((gr.tier == "Quantified").sum()),
                                         explained=float(q.get("signal_explained_by_library", np.nan)),
                                         reused=any("re-used" in ln for ln in lines)))

    queue = list(samples)
    # Spectra are fitted in worker processes with single-threaded linear algebra. Multi-threaded BLAS
    # gives tiny run-to-run rounding differences that can flip a grid-search choice; with one thread per worker the
    # results are identical every time, from the command line, the GUI or Python.
    for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                "NUMEXPR_NUM_THREADS"):
        os.environ[var] = "1"
    from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
    import multiprocessing
    # "spawn" everywhere: the Windows behaviour, and safe when called from a GUI
    with ProcessPoolExecutor(max_workers=nw, mp_context=multiprocessing.get_context("spawn")) as ex:
        running = {}
        while queue or running:
            while queue and len(running) < nw and not stop():
                smp = queue.pop(0)
                say("sample_started", smp, {})
                running[ex.submit(process_sample, jobs[smp])] = (smp, time.time())
            if stop() and queue:
                cancelled, queue = queue, []
            if not running:
                break
            done, _ = wait(list(running), timeout=1.0, return_when=FIRST_COMPLETED)
            for fut in done:
                smp, t0 = running.pop(fut)
                try:
                    output = fut.result()
                except Exception as exc:          # worker crashed
                    output = (None, None, None, dict(sample=smp, error=f"{type(exc).__name__}: {exc}"), None,
                              [f"Sample {smp} FAILED: {exc}"])
                finished(smp, output, t0)
    if cancelled:
        log(f"\nCancelled: {len(cancelled)} spectra not fitted: {cancelled}")
        say("cancelled", None, dict(not_fitted=cancelled))
    grades, peaks, obs, qc, failed = [], [], [], [], []
    for smp in samples:
        if smp not in results:
            continue
        gr, pk, ob, q, res, lines = results[smp]
        if gr is None:
            failed.append(smp); qc.append(q); continue
        grades.append(gr); peaks.append(pk); obs.append(ob); qc.append(q)
    if failed:
        log(f"\n*** {len(failed)} spectra failed and are not in the tables: {failed}")
    if not grades:
        logf.close()
        say("finished", None, dict(ok=False, failed=failed, cancelled=cancelled))
        raise RuntimeError("no spectrum could be fitted; see run_log.txt")
    say("aggregating", None, {})
    done_samples = [s_ for s_ in samples if s_ in results and results[s_][0] is not None]
    G = aggregate_outputs(out, grades, peaks, obs, qc, lib, done_samples, log)
    logf.close()
    say("finished", None, dict(ok=True, failed=failed, cancelled=cancelled, out=out))
    return G
