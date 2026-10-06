"""
mancq.fidproc - process raw Bruker 1D 1H FIDs into spectra MANC-Q can fit.

Steps, matching standard TopSpin processing (efp; apk; abs):
  1. read acqus + fid, remove the Bruker digital filter (group delay)
  2. exponential apodisation (line broadening, default 0.3 Hz)
  3. zero filling to SI points, Fourier transform
  4. phasing: automatic (zero order from the absorption integral, then a joint zero/first order refinement that
     minimises negative intensity), or user-supplied ph0/ph1
  5. optional baseline correction (off by default: MANC-Q fits a smooth baseline itself, and a separate
     correction here was found to distort broad signals such as the poloxamer PEO hump)
  6. referencing: the TSP/DSS singlet nearest 0 ppm is set to 0.000 ppm
The processed spectrum is written in Bruker format (pdata/1/1r + procs) to a folder of the caller's choosing, so
the raw data are never modified.
"""
import os
import warnings
import numpy as np
import nmrglue as ng

DEFAULTS = dict(lb=0.3, si=None, ph0=None, ph1=None, baseline=True, water_ppm=None, water_halfwidth=0.25)


def has_fid(expdir):
    return os.path.exists(os.path.join(expdir, "fid")) and os.path.exists(os.path.join(expdir, "acqus"))


def has_processed(expdir, procno=1):
    return os.path.exists(os.path.join(expdir, "pdata", str(procno), "1r"))


def _read(expdir):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dic, fid = ng.bruker.read(expdir, read_pulseprogram=False)
    if fid.ndim != 1:
        raise ValueError(f"{expdir}: not a 1D FID")
    fid = ng.bruker.remove_digital_filter(dic, fid)
    return dic, fid


def _transform(dic, fid, lb, si):
    a = dic["acqus"]
    t = np.arange(fid.size) / float(a["SW_h"])
    F = np.fft.fftshift(np.fft.fft(fid * np.exp(-np.pi * lb * t), n=si))
    return np.roll(F[::-1], 1)                       # high ppm first, as TopSpin stores 1r


def _axis(dic, si):
    a = dic["acqus"]
    bf1, sw_h, o1 = float(a["BF1"]), float(a["SW_h"]), float(a["O1"])
    first = (o1 + sw_h / 2) / bf1                    # ppm of point 0 before referencing (no SR)
    return first - np.arange(si) * sw_h / bf1 / si, bf1, sw_h


def _phase(F, p0, p1):
    x = np.linspace(0, 1, F.size, endpoint=False)
    return (F * np.exp(1j * np.deg2rad(p0 + p1 * x))).real


def autophase(F, keep, fit_ph1=False):
    """Zero-order phase: coarse search maximising the absorption integral, then a fine search minimising negative
    intensity (0.1 degree steps). After digital-filter removal Bruker spectra need almost no first-order phase, so
    ph1 is fixed at 0 unless fit_ph1=True. `keep` masks out the water region."""
    Fk = F[keep]
    p0s = np.arange(0, 360, 0.5)
    integ = [(Fk * np.exp(1j * np.deg2rad(p))).real.sum() for p in p0s]
    p0 = float(p0s[int(np.argmax(integ))])

    def cost(p0_, p1_):
        y = _phase(F, p0_, p1_)[keep]
        noise = 1.4826 * np.median(np.abs(y - np.median(y)))
        return float((np.clip(y, None, -3 * noise) ** 2).sum())

    cands = np.arange(p0 - 15, p0 + 15 + 1e-9, 0.1)
    p0 = float(cands[int(np.argmin([cost(c, 0.0) for c in cands]))])
    p1 = 0.0
    if fit_ph1:
        best = (cost(p0, 0.0), p0, 0.0)
        for a_ in np.arange(p0 - 3, p0 + 3 + 1e-9, 0.25):
            for b_ in np.arange(-10, 10 + 1e-9, 0.5):
                c = cost(a_, b_)
                if c < best[0]:
                    best = (c, a_, b_)
        p0, p1 = best[1], best[2]
    return p0 % 360, p1


def baseline_correct(ppm, y, keep, knot_ppm=0.5):
    """Smooth spline through points that look like baseline (iteratively excluding signal)."""
    from scipy.interpolate import make_lsq_spline
    x = ppm[::-1]
    yy = y[::-1]
    kk = keep[::-1].copy()
    noise = 1.4826 * np.median(np.abs(np.diff(yy))) / np.sqrt(2)
    base_mask = kk.copy()
    knots = np.arange(x[0] + knot_ppm, x[-1] - knot_ppm, knot_ppm)
    t = np.r_[[x[0]] * 4, knots, [x[-1]] * 4]
    b = np.zeros_like(yy)
    for _ in range(6):
        xs, ys = x[base_mask], yy[base_mask]
        try:
            spl = make_lsq_spline(xs, ys, t, k=3)
        except Exception:
            break
        b = spl(x)
        r = yy - b
        new = kk & (np.abs(r) < 4 * noise)
        if new.sum() < 50 or (new == base_mask).all():
            break
        base_mask = new
    return (yy - b)[::-1]


def reference(ppm, y, window=0.4):
    """Shift the axis so the tallest sharp peak within +/-window ppm of 0 sits at 0.000 ppm."""
    m = np.abs(ppm) < window
    i = np.argmax(np.where(m, y, -np.inf))
    return ppm - ppm[i]


def prepare(expdir, lb=0.3, si=None, water_halfwidth=0.25):
    """Read, filter-correct, apodise, zero-fill and transform. The result can be phased repeatedly and quickly."""
    dic, fid = _read(expdir)
    si = int(si or 1 << int(np.ceil(np.log2(2 * fid.size))))   # zero-fill x2, as TopSpin SI = TD
    F = _transform(dic, fid, lb, si)
    ppm, bf1, sw_hz = _axis(dic, si)
    water = float(dic["acqus"]["O1"]) / bf1          # presaturation carrier = water
    keep = np.abs(ppm - water) > water_halfwidth
    return dict(F=F, ppm=ppm, keep=keep, sf=bf1, sw_hz=sw_hz, lb=lb, si=si)


def finish(prep, ph0, ph1, baseline=False):
    """Phase (degrees), baseline-correct and reference a prepared spectrum."""
    y = _phase(prep["F"], ph0, ph1)
    if baseline:
        y = baseline_correct(prep["ppm"], y, prep["keep"])
    ppm = reference(prep["ppm"], y)
    return dict(ppm=ppm, y=y, sf=prep["sf"], sw_hz=prep["sw_hz"], ph0=float(ph0) % 360, ph1=float(ph1),
                lb=prep["lb"], si=prep["si"])


def process_fid(expdir, lb=0.3, si=None, ph0=None, ph1=None, baseline=False, water_halfwidth=0.25):
    """Full processing. ph0/ph1 None = automatic phasing. Returns dict(ppm, y, sf, sw_hz, ph0, ph1, auto, ...)
    with ppm descending (point 0 = highest ppm)."""
    prep = prepare(expdir, lb, si, water_halfwidth)
    auto = ph0 is None or ph1 is None
    if auto:
        ph0, ph1 = autophase(prep["F"], prep["keep"])
    res = finish(prep, ph0, ph1, baseline)
    res["auto"] = auto
    return res


def write_processed(res, folder, title=""):
    """Write a processed spectrum as <folder>/pdata/1/{1r,procs,title} (Bruker format)."""
    pdir = os.path.join(folder, "pdata", "1")
    os.makedirs(pdir, exist_ok=True)
    y = np.asarray(res["y"], float)
    big = np.abs(y).max()
    nc = int(np.ceil(np.log2(big / 2 ** 29))) if big > 2 ** 29 else 0
    np.round(y / 2.0 ** nc).astype("<i4").tofile(os.path.join(pdir, "1r"))
    procs = {"SF": res["sf"], "SW_p": res["sw_hz"], "OFFSET": float(res["ppm"][0]), "SI": len(y), "FTSIZE": len(y),
             "STSI": len(y), "STSR": 0, "NC_proc": nc, "BYTORDP": 0, "DTYPP": 0, "AXNUC": "<1H>",
             "LB": res["lb"], "PHC0": res["ph0"], "PHC1": res["ph1"]}
    text = "##TITLE= MANC-Q processed\n##JCAMPDX= 5.0\n" + "".join(f"##${k}= {v}\n" for k, v in procs.items()) + "##END=\n"
    for fn in ("procs", "proc"):
        open(os.path.join(pdir, fn), "w").write(text)
    open(os.path.join(pdir, "title"), "w").write(title or "processed from FID by MANC-Q")
    return folder
