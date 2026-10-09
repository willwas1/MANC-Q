"""
mancq.fidproc - process raw Bruker 1D 1H FIDs into spectra MANC-Q can fit, and adjust spectra already processed
in TopSpin.

TopSpin-processed spectra (from version 1.2) are adjusted starting from TopSpin's own result: the real and
imaginary parts (pdata/<procno>/1r and 1i) are re-phased by the change the user asks for, and extra line
broadening is applied by going back to the time domain, so "no change" gives exactly TopSpin's spectrum. If 1i
has been deleted but the raw FID is there, the FID is processed with TopSpin's line broadening and phased to
match TopSpin's 1r, and changes are applied from there.

Raw FID steps, matching standard TopSpin processing (efp; apk; abs):
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
    """Phase (degrees), baseline-correct and reference a prepared spectrum. For a TopSpin spectrum, ph0 and ph1 are
    changes from TopSpin's phasing (0, 0 = unchanged)."""
    b0, b1 = prep.get("ph_base", (0.0, 0.0))
    y = _phase(prep["F"], b0 + ph0, b1 + ph1)
    if baseline:
        y = baseline_correct(prep["ppm"], y, prep["keep"])
    ppm = reference(prep["ppm"], y)
    rel = prep.get("kind") == "topspin"
    return dict(ppm=ppm, y=y, sf=prep["sf"], sw_hz=prep["sw_hz"],
                ph0=(float(ph0) + 180) % 360 - 180 if rel else float(ph0) % 360, ph1=float(ph1),
                lb=prep["lb"], si=prep["si"], baseline=bool(baseline), kind=prep.get("kind", "fid"),
                lb_topspin=prep.get("lb_topspin"), how=prep.get("how", "raw FID"))


# ----------------------------------------------------------------------------- TopSpin-processed spectra
def _pdir(expdir, procno):
    return os.path.join(expdir, "pdata", str(procno))


def has_imaginary(expdir, procno=1):
    return os.path.exists(os.path.join(_pdir(expdir, procno), "1i"))


def can_adjust(expdir, procno=1):
    """A TopSpin spectrum can be re-phased / re-broadened if its imaginary part (1i) or the raw FID is there."""
    return has_processed(expdir, procno) and (has_imaginary(expdir, procno) or has_fid(expdir))


def topspin_params(expdir, procno=1):
    """TopSpin's processing parameters: dict(lb, wdw, phc0, phc1, si). lb is the exponential line broadening
    (0 if TopSpin used another window function, wdw != 1)."""
    p = ng.bruker.read_jcamp(os.path.join(_pdir(expdir, procno), "procs"))
    wdw = int(float(p.get("WDW", 1)))
    lb = float(p.get("LB", 0.0)) if wdw == 1 else 0.0
    return dict(lb=lb, wdw=wdw, phc0=float(p.get("PHC0", 0.0)), phc1=float(p.get("PHC1", 0.0)),
                si=int(p.get("SI", 0)), procs=p)


def broaden(F, dlb, sw_hz, n_acq=None):
    """Extra exponential line broadening dlb (Hz; negative narrows) of a complex spectrum, by inverse Fourier
    transform, multiplication of the time signal and transform back. Works whichever end of the array the
    time signal starts at (TopSpin stores high ppm first, which reverses time). n_acq = acquired complex points,
    if known: anything later is zero filling and is set to exactly zero."""
    N = F.size
    s = np.fft.ifft(F)
    a = np.abs(s)
    k = max(1, N // 16)
    origin = 0 if a[:k].sum() + a[-k:].sum() >= a[N // 2 - k:N // 2 + k].sum() else N // 2
    rel = (np.arange(N) - origin) % N
    fwd = a[(rel > 0) & (rel < N // 4)].sum()
    bwd = a[rel > 3 * N // 4].sum()
    if fwd >= bwd:
        t = np.where(rel < N - N // 32, rel, N - rel)
    else:
        t = np.where(rel > N // 32, N - rel, rel)
    if n_acq:
        s = np.where(t < n_acq, s, 0)
    return np.fft.fft(s * np.exp(-np.pi * dlb * t / float(sw_hz)))


def _acqus(expdir):
    f = os.path.join(expdir, "acqus")
    return ng.bruker.read_jcamp(f) if os.path.exists(f) else None


def _match_phase(F, y, keep):
    """Phase (ph0, ph1 in MANC-Q's convention) that makes Re(F) best match TopSpin's real spectrum y, and the
    remaining squared misfit."""
    x = np.linspace(0, 1, F.size, endpoint=False)
    m = keep & (np.abs(y) > 20 * 1.4826 * np.median(np.abs(y - np.median(y))))
    if m.sum() < 50:
        m = keep
    best = None
    for p1 in np.arange(-30, 30 + 1e-9, 0.5):
        G = F[m] * np.exp(1j * np.deg2rad(p1 * x[m]))
        A = np.c_[G.real, G.imag, np.ones(m.sum())]
        c, res, *_ = np.linalg.lstsq(A, y[m], rcond=None)
        r = float(((A @ c - y[m]) ** 2).sum())
        if best is None or r < best[0]:
            best = (r, float(np.rad2deg(np.arctan2(-c[1], c[0]))), float(p1))
    return best[1] % 360, best[2], best[0]


def prepare_topspin(expdir, procno=1, lb=None, water_halfwidth=0.25):
    """A TopSpin-processed spectrum ready for adjustment, in the same form as prepare(). lb = total line
    broadening wanted (None = TopSpin's). Phase changes passed to finish() are relative to TopSpin's phasing."""
    pdir = _pdir(expdir, procno)
    tp = topspin_params(expdir, procno)
    p = tp["procs"]
    sf, sw_hz, off = float(p["SF"]), float(p["SW_p"]), float(p["OFFSET"])
    aq = _acqus(expdir)
    water = float(aq["O1"]) / float(aq["BF1"]) if aq and "O1" in aq else 4.75
    lb = tp["lb"] if lb is None else float(lb)
    if has_imaginary(expdir, procno):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, data = ng.bruker.read_pdata(pdir, all_components=True, read_acqus=False, scale_data=True)
        F = np.asarray(data[0], float) + 1j * np.asarray(data[1], float)
        si = F.size
        ppm = off - np.arange(si) * sw_hz / sf / si
        if abs(lb - tp["lb"]) > 1e-9:
            n_acq = int(aq["TD"]) // 2 if aq and "TD" in aq else None
            F = broaden(F, lb - tp["lb"], sw_hz, n_acq)
        keep = np.abs(ppm - water) > water_halfwidth
        return dict(F=F, ppm=ppm, keep=keep, sf=sf, sw_hz=sw_hz, lb=lb, si=si, ph_base=(0.0, 0.0), kind="topspin",
                    lb_topspin=tp["lb"], how="TopSpin spectrum (1r and 1i)")
    if not has_fid(expdir):
        raise ValueError(f"{expdir}: neither the imaginary part (1i) nor the raw FID is there, so this spectrum "
                         "can only be used as TopSpin processed it")
    # no 1i: process the FID with the same line broadening and size, then match TopSpin's phase
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, y_ts = ng.bruker.read_pdata(pdir, read_acqus=False, scale_data=True)
    y_ts = np.asarray(y_ts, float)
    si = y_ts.size
    ppm_ts = off - np.arange(si) * sw_hz / sf / si
    pr = prepare(expdir, lb=lb, si=si, water_halfwidth=water_halfwidth)
    # same size and spectral width, so the two axes differ by a whole number of points (TopSpin's referencing):
    # cross-correlating the magnitude spectra finds it to within a point or two; the phase match decides
    keep = np.abs(ppm_ts - water) > water_halfwidth
    Fm = pr["F"] if abs(lb - tp["lb"]) <= 1e-9 else prepare(expdir, lb=tp["lb"], si=si,
                                                            water_halfwidth=water_halfwidth)["F"]
    a, b = np.abs(Fm), np.abs(y_ts)
    lag0 = int(np.argmax(np.fft.ifft(np.fft.fft(b) * np.conj(np.fft.fft(a))).real))
    best = None
    for lag in range(lag0 - 3, lag0 + 4):
        b0, b1, r = _match_phase(np.roll(Fm, lag), y_ts, keep)
        if best is None or r < best[0]:
            best = (r, lag, b0, b1)
    _, lag, b0, b1 = best
    F = np.roll(pr["F"], lag)
    return dict(F=F, ppm=ppm_ts, keep=keep, sf=sf, sw_hz=sw_hz, lb=lb, si=si, ph_base=(b0, b1), kind="topspin",
                lb_topspin=tp["lb"], how="raw FID phased to match TopSpin's 1r (no 1i)")


def autophase_relative(prep):
    """Automatic phase for a prepared spectrum, as a change from its starting phase (ph_base): (ph0, ph1) with
    ph0 in -180..180."""
    b0, b1 = prep.get("ph_base", (0.0, 0.0))
    x = np.linspace(0, 1, prep["F"].size, endpoint=False)
    p0, p1 = autophase(prep["F"] * np.exp(1j * np.deg2rad(b0 + b1 * x)), prep["keep"])
    return (p0 + 180) % 360 - 180, p1


def describe(res):
    """One line saying how a spectrum was processed (for processing.txt, the log and the QC table)."""
    bl = "baseline correction on" if res.get("baseline") else "no baseline correction"
    if res.get("kind") == "topspin":
        return (f"{res.get('how')}: phase change ph0 {res['ph0']:+.2f}, ph1 {res['ph1']:+.2f} deg from TopSpin; "
                f"line broadening {res['lb']:.2f} Hz (TopSpin {res.get('lb_topspin') or 0:.2f} Hz); {bl}")
    return (f"raw FID: line broadening {res['lb']:.2f} Hz, SI {res['si']}, "
            f"{'automatic phase' if res.get('auto') else 'phase'} ph0 {res['ph0']:.2f}, ph1 {res['ph1']:.2f} deg; {bl}")


def process(expdir, procno=1, source="fid", lb=None, ph0=None, ph1=None, baseline=False):
    """Process one experiment. source "fid": the raw FID (lb default 0.3 Hz; ph0/ph1 None = automatic phase,
    absolute angles). source "topspin": TopSpin's spectrum with changes (lb = total line broadening, None =
    TopSpin's; ph0/ph1 = phase change from TopSpin, None = automatic). Returns the dict from finish() plus
    'auto' and 'description'."""
    if source == "fid":
        res = process_fid(expdir, lb=0.3 if lb is None else lb, ph0=ph0, ph1=ph1, baseline=baseline)
    elif source == "topspin":
        prep = prepare_topspin(expdir, procno, lb)
        auto = ph0 is None
        if auto:
            ph0, ph1 = autophase_relative(prep)
        res = finish(prep, ph0, ph1 or 0.0, baseline)
        res["auto"] = auto
    else:
        raise ValueError(f"unknown source {source!r} (use 'fid' or 'topspin')")
    res["description"] = describe(res)
    return res


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
