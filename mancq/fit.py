"""
mancq.fit - fitting one spectrum: the reference (TSP/DSS), noise, the pH estimate, the lineshape and the
whole-spectrum fit with shift and linewidth freedom, and the tests for multiplets missing from the data.
"""
import logging
import time
import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import least_squares
from scipy.optimize import lsq_linear
from scipy.optimize import nnls
from scipy.signal import find_peaks
_log = logging.getLogger("mancq")      # see mancq.log (window) or the console
from .settings import TIMING


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


def noise_sigma(ppm, y, why=None):
    """Noise standard deviation from the signal-free ends of the spectrum (10-12 and -2.5 to -0.5 ppm). If
    neither has 200 points (a narrow spectral width), a robust estimate from the point-to-point differences of
    the whole spectrum is used instead, and the reason is appended to the list `why`."""
    parts = []
    for lo, hi in [(10.0, 12.0), (-2.5, -0.5)]:
        m = (ppm > lo) & (ppm < hi)
        if m.sum() > 200:
            yy = y[m] - np.polyval(np.polyfit(ppm[m], y[m], 2), ppm[m])
            parts.append(1.4826 * np.median(np.abs(yy - np.median(yy))))
    if parts:
        return float(np.min(parts))
    if why is not None:
        why.append("the spectrum has no signal-free region at 10-12 or -2.5 to -0.5 ppm, so the noise was estimated "
                   "from point-to-point differences (less reliable: check LOD and SNR)")
    d2 = y[2:] - 2 * y[1:-1] + y[:-2]               # second difference: removes smooth lines and the baseline
    return float(1.4826 * np.median(np.abs(d2 - np.median(d2))) / np.sqrt(6.0))


def fit_tsp(ppm, y, sf, *, S):
    """Pseudo-Voigt main line + 29Si doublet (J 6.6 Hz) + linear baseline."""
    m = (ppm > -0.2) & (ppm < 0.2)
    if not m.any():
        raise ValueError("the spectrum does not cover 0 ppm, where the TSP/DSS reference must be")
    x0 = ppm[m][np.argmax(y[m])]
    if y[m].max() < 20 * noise_sigma(ppm, y):
        raise ValueError(f"no {S['REFERENCE']} peak near 0 ppm. MANC-Q needs the TSP/DSS internal standard, "
                         f"referenced to 0.000 ppm: check that it was added and that the spectrum is referenced "
                         f"in TopSpin")
    w = (ppm > x0 - 0.03) & (ppm < x0 + 0.03)
    xx, yy = ppm[w], y[w]
    h = yy.max()
    fsi = S["TSP_SI29_FRACTION"]

    def disp(x, c, fwhm):
        g = fwhm / 2.0
        return ((x - c) / np.pi) / ((x - c) ** 2 + g ** 2)

    def model(p):
        c, f, eta, A, b0, b1, jsi, D = p
        fp = f / sf
        main = A * pv(xx, c, fp, eta) + D * disp(xx, c, fp)       # D: dispersive part (residual phase error)
        sat = A * fsi / (1 - fsi) / 2 * (pv(xx, c - jsi / 2 / sf, fp, eta) + pv(xx, c + jsi / 2 / sf, fp, eta))
        return main + sat + b0 + b1 * (xx - c)

    A0 = h * 1.2 / sf * np.pi / 2
    p0 = [x0, 1.2, 0.6, A0, 0, 0, 6.6, 0.0]
    lb = [x0 - 0.002, 0.3, 0.0, 0, -h, -np.inf, 5.5, -0.5 * A0]
    ub = [x0 + 0.002, 6.0, 1.0, np.inf, h, np.inf, 7.5, 0.5 * A0]
    r = least_squares(lambda p: model(p) - yy, p0, bounds=(lb, ub), x_scale="jac")
    c, f, eta, A, b0, b1, jsi, D = r.x
    resid = np.sqrt(np.mean(r.fun ** 2)) / h
    area_total = A / (1 - fsi)          # main + 29Si satellites = all TSP protons
    return dict(centre=c, fwhm_hz=f, eta=eta, area_main=A, area=area_total, height=h, rel_rms=resid, jsi=jsi,
                phase_deg=float(np.degrees(np.arctan2(-D, A))) if A > 0 else np.nan)   # size matters; sign is not TopSpin's


# =============================================================================
# 2b. Sample pH estimate (histidine imidazole H2 and H4)
# =============================================================================
PH_COLUMNS = ("pH_est", "pH_H2", "pH_H4", "his_H2_ppm", "his_H4_ppm", "pH_flag")


def ph_from_shift(d_obs, limits, pka):
    """Henderson-Hasselbalch: pH from an observed shift, limits = (shift protonated, shift deprotonated)."""
    d_ha, d_a = limits
    f = np.clip((d_ha - d_obs) / (d_ha - d_a), 0.001, 0.999)       # fraction deprotonated
    return float(pka + np.log10(f / (1 - f)))


def shift_from_ph(ph, limits, pka):
    """The shift Henderson-Hasselbalch predicts at a given pH (inverse of ph_from_shift)."""
    d_ha, d_a = limits
    f = 1.0 / (1.0 + 10 ** (pka - ph))
    return float(d_ha - f * (d_ha - d_a))


def _peak_candidates(ppm, y, sigma, window, *, S):
    """Peaks in a ppm window at least PH_MIN_SNR x noise tall: [(position, height SNR)], position refined with a
    3-point parabola through the maximum. ppm must be ascending and evenly spaced."""
    lo, hi = window
    idx = np.where((ppm >= lo) & (ppm <= hi))[0]
    if idx.size < 5:
        return []
    base = float(np.percentile(y[idx], 10))          # local baseline under the window
    pk, _ = find_peaks(y[idx] - base, height=S["PH_MIN_SNR"] * sigma)
    out = []
    for p in pk:
        i = idx[p]
        y0, y1, y2 = y[i - 1], y[i], y[i + 1]
        den = y0 - 2 * y1 + y2
        delta = 0.5 * (y0 - y2) / den if den < 0 else 0.0
        delta = float(np.clip(delta, -0.5, 0.5))
        pos = ppm[i] + delta * (ppm[i + 1] - ppm[i])
        height = y1 - 0.25 * (y0 - y2) * delta - base
        out.append((float(pos), float(height / sigma)))
    return out


def estimate_ph(ppm, y, sigma, *, S):
    """Estimated pH in the NMR tube from the histidine H2 and H4 singlets (TSP-referenced ppm, ascending).
    Returns pH_est (mean of the protons found), pH_H2, pH_H4, his_H2_ppm, his_H4_ppm and pH_flag
    (ok / one_proton_only / protons_disagree / no_histidine). An estimate: see the PH_ settings."""
    pka, l2, l4 = S["PH_HIS_PKA"], tuple(S["PH_HIS_H2_LIMITS"]), tuple(S["PH_HIS_H4_LIMITS"])
    c2 = _peak_candidates(ppm, y, sigma, S["PH_HIS_H2_WINDOW"], S=S)
    c4 = _peak_candidates(ppm, y, sigma, S["PH_HIS_H4_WINDOW"], S=S)
    h2 = max(c2, key=lambda c: c[1]) if c2 else None        # tallest peak in each window
    h4 = max(c4, key=lambda c: c[1]) if c4 else None
    if h2 and h4 and abs(ph_from_shift(h2[0], l2, pka) - ph_from_shift(h4[0], l4, pka)) > S["PH_MAX_DISAGREE"]:
        # the tallest peaks disagree: a neighbour (e.g. the Trp doublet at 7.74 ppm) may have been taken.
        # Prefer the pair of peaks that agree, the stronger pair if several do
        ok = [(a, b) for a in c2 for b in c4
              if abs(ph_from_shift(a[0], l2, pka) - ph_from_shift(b[0], l4, pka)) <= S["PH_MAX_DISAGREE"]]
        if ok:
            h2, h4 = max(ok, key=lambda ab: min(ab[0][1], ab[1][1]))
    p2 = ph_from_shift(h2[0], l2, pka) if h2 else np.nan
    p4 = ph_from_shift(h4[0], l4, pka) if h4 else np.nan
    have = [v for v in (p2, p4) if np.isfinite(v)]
    if not have:
        flag = "no_histidine"
    elif len(have) == 1:
        flag = "one_proton_only"
    else:
        flag = "protons_disagree" if abs(p2 - p4) > S["PH_MAX_DISAGREE"] else "ok"
    return dict(pH_est=round(float(np.mean(have)), 3) if have else np.nan,
                pH_H2=round(p2, 3) if h2 else np.nan, pH_H4=round(p4, 3) if h4 else np.nan,
                his_H2_ppm=round(h2[0], 5) if h2 else np.nan, his_H4_ppm=round(h4[0], 5) if h4 else np.nan,
                pH_flag=flag)


def histidine_targets(ph, *, S):
    """{'H2'|'H4': (target ppm, half window Hz)} for placing histidine in the fit from a pH estimate: the
    measured peak where it was seen, else the position the other proton's pH predicts (wider window)."""
    out = {}
    if not ph or not np.isfinite(ph.get("pH_est", np.nan)):
        return out
    for tag, lim in (("H2", S["PH_HIS_H2_LIMITS"]), ("H4", S["PH_HIS_H4_LIMITS"])):
        meas = ph.get(f"his_{tag}_ppm", np.nan)
        if np.isfinite(meas):
            out[tag] = (float(meas), S["PH_HIS_PLACE_HZ"])
        else:
            out[tag] = (shift_from_ph(ph["pH_est"], tuple(lim), S["PH_HIS_PKA"]), S["PH_HIS_PLACE_PREDICTED_HZ"])
    return out


# =============================================================================
# 4. Model pieces
# =============================================================================


def mcol(x, p, it, d, fw, eta, *, half):
    """Column (on sorted x) of one multiplet: sticks p (ppm), intensities it,
    shifted by d (ppm), pseudo-Voigt fwhm fw (ppm); each line evaluated within +/-half (ppm) of its centre."""
    out = np.zeros(len(x))
    c = p + d
    lo = np.searchsorted(x, c.min() - half)
    hi = np.searchsorted(x, c.max() + half)
    if hi > lo:
        xx = x[lo:hi]
        if len(c) * (hi - lo) > 4e6:
            for cc, ii in zip(c, it):
                out[lo:hi] += ii * pv(xx, cc, fw, eta)
        else:
            out[lo:hi] = (it[:, None] * pv(xx[None, :], c[:, None], fw, eta)).sum(0)
    return out


def mcol_d(x, p, it, d, fw, eta, *, half):
    """Multiplet column and its analytic derivatives w.r.t. shift d and fwhm fw."""
    n = len(x)
    col, dd, df = np.zeros(n), np.zeros(n), np.zeros(n)
    c = p + d
    lo = np.searchsorted(x, c.min() - half)
    hi = np.searchsorted(x, c.max() + half)
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
def fit_spectrum(ppm_raw, y, sf, lib, log, *, S):
    """Fit one spectrum. Model:
         y(x) = sum_m a_m * M_m(x; d_m, w_m)  +  sum_u b_u * singlet_u(x)  +  spline(x)
       M_m = simulated multiplet m (library intensities, locked internal geometry)
       a_m >= 0, softly tied to be equal (per proton) within a compound
       d_m = library position + global offset + compound shift D_k + small delta_m
       w_m = linewidth factor (pseudo-Voigt, shared eta)"""
    t0 = time.time()
    HW = S["LINE_WINDOW_PPM"]                 # each line is computed within +/-HW ppm of its centre
    tsp = fit_tsp(ppm_raw, y, sf, S=S)
    ppm = ppm_raw - tsp["centre"]
    noise_why = []
    sigma = noise_sigma(ppm, y, noise_why)
    for w_ in noise_why:
        log(f"   WARNING: {w_}")
    fw0 = float(np.clip(tsp["fwhm_hz"], S["LINEWIDTH_BASE_HZ"][0], S["LINEWIDTH_BASE_HZ"][1])) / sf
    eta = float(np.clip(tsp["eta"], 0.3, 1.0))
    log(f"   {S['REFERENCE']}: area={tsp['area']:.4g}  fwhm={tsp['fwhm_hz']:.2f} Hz  eta={tsp['eta']:.2f}  "
        f"rel.rms={tsp['rel_rms']:.3f}  noise sigma={sigma:.4g}  SNR(TSP)={tsp['height']/sigma:.0f}")
    ph = estimate_ph(ppm, y, sigma, S=S) if S["PH_ESTIMATE"] else dict({k: np.nan for k in PH_COLUMNS}, pH_flag="not_estimated")
    if S["PH_ESTIMATE"]:
        log(f"   pH estimate (histidine, estimate only): {ph['pH_est']:.2f}  H2 {ph['his_H2_ppm']:.4f} ppm -> "
            f"{ph['pH_H2']:.2f}, H4 {ph['his_H4_ppm']:.4f} ppm -> {ph['pH_H4']:.2f}  ({ph['pH_flag']})")
    his_win = histidine_targets(ph, S=S) if S["PH_PLACE_HISTIDINE"] else {}

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
                col = sum(mcol(x, m["ppm"], m["inten"], d, fw0, eta, half=HW) for m in c["mults"])
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
            # histidine H2/H4 with a pH estimate: kept within a few Hz of the measured peak (absolute window on d),
            # so they cannot move onto a neighbour's peak; the compound shift is still chosen as for other compounds
            win = None
            if c["name"] == "L-Histidine" and his_win:
                for tag in ("H2", "H4"):
                    lo_, hi_ = S[f"PH_HIS_{tag}_WINDOW"]
                    if tag in his_win and lo_ <= m["centre"] <= hi_:
                        tgt, hw = his_win[tag]
                        win = (tgt - m["centre"] - hw / sf, tgt - m["centre"] + hw / sf)
            # with a pH estimate, the side-chain multiplets of histidine need only the usual pH-sensitive freedom
            # (the large freedom in the library is for H2/H4); otherwise a doublet of doublets can slip by one coupling
            dlt = DLT[k]
            if c["name"] == "L-Histidine" and his_win and win is None:
                dlt = min(dlt, S["MULTIPLET_SHIFT_TITRATABLE_HZ"] / sf)
            MP.append(dict(k=k, j=j, p=m["ppm"], it=m["inten"], nprot=m["nprot"], base=base, win=win, dlt=dlt,
                           d=base, dlo=base - dlt if win is None else win[0], dhi=base + dlt if win is None else win[1],
                           wf=c.get("wf0", 1.0), wlo=c.get("wlo", 0.5), whi=c.get("whi", 4.0),
                           titr=c["titratable"], local_amp=np.nan, local_se=np.nan,
                           loc_dom=0.0, loc_mis=np.inf, loc_h=0.0))
    nM = len(MP)
    MPk = {k: [i for i, m in enumerate(MP) if m["k"] == k] for k in range(ncomp)}
    B, Dpen = spline_basis(x, lo, hi, S["BASELINE_KNOT_PPM"])
    unknown = []
    dropped = {}                  # compounds left out of the fit (ghosts): {k: reason}

    def mult_columns(full=False):
        M = np.zeros((len(x), nM))
        for i, m in enumerate(MP):
            if full or m["k"] not in dropped:
                M[:, i] = mcol(x, m["p"], m["it"], m["d"], fw0 * m["wf"], eta, half=HW)
        return M

    def find_ghosts(coef, settled):
        """Compounds that only sit on other compounds' peaks (see DROP_GHOSTS); returns how many were dropped."""
        base_now = B @ coef[nM + len(unknown):]
        n = 0
        for k, ii in MPk.items():
            if k in dropped or len(ii) < 2 or lib[k]["source"] == "polymer":
                continue
            mps = [MP[i] for i in ii]
            carry = [m for m in mps if np.isfinite(m["local_amp"]) and m["local_amp"] > 0
                     and m["loc_h"] >= S["LOQ_SNR"] * sigma]
            if not carry or any(m["loc_dom"] >= S["GHOST_MAX_DOMINANCE"] for m in carry):
                continue
            cols = {id(m): mcol(x, m["p"], m["it"], m["d"], fw0 * m["wf"], eta, half=HW) for m in mps}
            nprot = sum(m["nprot"] for m in mps)
            for m in sorted(carry, key=lambda m: -m["loc_h"]):
                centre = float((m["p"] * m["it"]).sum() / m["it"].sum() + m["d"])
                where = (f"left out of the fit: its signal at {centre:.3f} ppm (multiplet {m['j'] + 1}, "
                         f"{m['loc_dom']:.0%} of the fitted signal there)")
                gone = _contradicted_by(m, mps, m["local_amp"], x, yv, base_now, fw0, sf, sigma, grid.dx, cols, S=S)
                if gone:
                    dropped[k] = f"{where} implies multiplet(s) {','.join(map(str, gone))}, which are not in the spectrum"
                    n += 1
                    break
                # moved well away from its library position, and the multiplets it implies get nothing in the fit
                # (they may sit under another compound's peaks, so the data alone cannot show they are missing)
                shift_hz = abs(m["d"] - m["base"]) * sf
                if not settled or shift_hz < S["GHOST_SHIFT_TITRATABLE_HZ" if lib[k]["titratable"] else "GHOST_SHIFT_HZ"]:
                    continue
                empty = [o["j"] + 1 for o in mps if o is not m and (o["nprot"] >= 0.3 or o["nprot"] / nprot >= 0.1)
                         and _excluded_fraction(o, S=S) < 0.5
                         and m["local_amp"] * cols[id(o)].max() > S["CONTRADICT_MIN_SNR"] * sigma
                         and not (np.isfinite(o["local_amp"])
                                  and o["local_amp"] >= S["CONTRADICT_FRACTION"] * m["local_amp"])]
                if empty:
                    dropped[k] = (f"{where} is {shift_hz:.1f} Hz from the library position, and its multiplet(s) "
                                  f"{','.join(map(str, empty))} get almost nothing in the fit")
                    n += 1
                    break
        if n:
            log(f"   left out of the fit (signal only on other compounds' peaks, other multiplets absent): "
                f"{', '.join(lib[k]['name'] for k in dropped)}")
        return n

    def unk_columns():
        if not unknown:
            return np.zeros((len(x), 0))
        return np.array([mcol(x, np.array([u["pos"]]), np.array([1.0]), 0.0, fw0 * u["wf"], eta, half=HW) for u in unknown]).T

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

    def roi_refine(coef, M, U, fit, max_nfev=None):
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
                           [unknown[i]["pos"] - S["UNKNOWN_MOVE_HZ"] / sf for i in us] + [S["UNKNOWN_WIDTH_RANGE"][0]] * nu, float)
            ubv = np.array([MP[i]["dhi"] for i in ms] + [MP[i]["whi"] for i in ms] +
                           [unknown[i]["pos"] + S["UNKNOWN_MOVE_HZ"] / sf for i in us] + [S["UNKNOWN_WIDTH_RANGE"][1]] * nu, float)
            v = np.clip(v, lbv, ubv)
            tt = (xs - xs.mean()) / max(np.ptp(xs), 1e-9)
            blin = np.array([np.ones(len(xs)), -np.ones(len(xs)), tt, -tt]).T

            def lin(v):
                cols_, dd_, dw_ = [], [], []
                for a_, i in enumerate(ms):
                    c0, c1, c2 = mcol_d(xs, MP[i]["p"], MP[i]["it"], v[a_], fw0 * v[nm + a_], eta, half=HW)
                    cols_.append(c0); dd_.append(c1); dw_.append(c2 * fw0)
                for a_, i in enumerate(us):
                    c0, c1, c2 = mcol_d(xs, np.array([0.0]), np.array([1.0]), v[2 * nm + a_], fw0 * v[2 * nm + nu + a_], eta, half=HW)
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
            for _ in range(S["ROI_MAX_NFEV"] if max_nfev is None else max_nfev):
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
                _log.debug("no standard errors for one fitted region", exc_info=True)
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
            if (lib[k]["source"] == "polymer" and lib[k].get("dmax_hz", 3.0) <= 3.0) or not MPk[k] or k in dropped:
                continue
            ii = MPk[k]
            spans = [(MP[i]["p"].min() + MP[i]["base"] - DMAX[k] - MP[i]["dlt"] - HW / 2,
                      MP[i]["p"].max() + MP[i]["base"] + DMAX[k] + MP[i]["dlt"] + HW / 2) for i in ii]
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
                    col += ra * mcol(xs, MP[i]["p"], MP[i]["it"], MP[i]["base"] + Dc + dl, fw0 * MP[i]["wf"], eta, half=HW)
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
                    if MP[i]["win"]:          # the compound shift is chosen as usual; held multiplets stay in their window
                        MP[i]["d"] = float(np.clip(MP[i]["d"], *MP[i]["win"]))
                col = best[3]
                if any(MP[i]["win"] for i in ii):
                    col = sum(ra * mcol(xs, MP[i]["p"], MP[i]["it"], MP[i]["d"], fw0 * MP[i]["wf"], eta, half=HW)
                              for i, ra in zip(ii, rel_amp))
                fit[sup] = fit[sup] - own + best[2] * col
            # pH-sensitive compounds: each multiplet may move further on its own (multiplets held near a
            # measured peak move only inside their own window)
            if any(MP[i]["dlt"] > 2.0 / sf or MP[i]["win"] for i in ii):
                for i in ii:
                    m = MP[i]
                    if m["win"]:
                        d_c, hw_ = 0.5 * (m["win"][0] + m["win"][1]), 0.5 * (m["win"][1] - m["win"][0])
                    elif m["dlt"] > 2.0 / sf:
                        d_c, hw_ = m["base"] + Dk[k], m["dlt"]
                    else:
                        continue
                    lo_ = np.searchsorted(x, m["p"].min() + d_c - hw_ - HW / 2)
                    hi_ = np.searchsorted(x, m["p"].max() + d_c + hw_ + HW / 2)
                    if hi_ - lo_ < 5:
                        continue
                    xs2 = x[lo_:hi_]
                    w2 = wts[lo_:hi_] ** 2
                    a_i = max(am[i], 0.0)
                    own2 = a_i * mcol(xs2, m["p"], m["it"], m["d"], fw0 * m["wf"], eta, half=HW)
                    r2 = yv[lo_:hi_] - fit[lo_:hi_] + own2
                    bb = (np.inf, m["d"], None, 0.0)
                    for dd in np.arange(-hw_, hw_ + 1e-12, step):
                        col = mcol(xs2, m["p"], m["it"], d_c + dd, fw0 * m["wf"], eta, half=HW)
                        den = (w2 * col * col).sum()
                        if den <= 0:
                            continue
                        aa = max(0.0, (w2 * col * r2).sum() / den)
                        # mild preference for the library position (breaks ties between neighbours)
                        ss = (w2 * (r2 - aa * col) ** 2).sum() * (1.0 + S["SHIFT_PRIOR"] * (dd * sf / 10.0) ** 2)
                        if ss < bb[0]:
                            bb = (ss, d_c + dd, col, aa)
                    if bb[2] is not None:
                        m["d"] = bb[1]
                        fit[lo_:hi_] = fit[lo_:hi_] - own2 + bb[3] * bb[2]
        for m in MP:
            k = m["k"]
            m["dlo"] = m["base"] + Dk[k] - m["dlt"] if not m["win"] else m["win"][0]
            m["dhi"] = m["base"] + Dk[k] + m["dlt"] if not m["win"] else m["win"][1]
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
        # missing multiplets in the data: from the first iteration; nothing in the fit: once shifts have settled
        if S["DROP_GHOSTS"] and find_ghosts(coef, settled=it >= 1):
            M = mult_columns()
            coef, fit, A = global_solve(M, U)

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
        unknown.append(dict(pos=float(x[p]), wf=S["UNKNOWN_WIDTH_START"]))
    U = unk_columns()
    coef, fit, A = global_solve(M, U)
    log(f"   +{len(unknown)} unassigned singlets -> chi2/n={chi(fit):.2f}")
    roi_refine(coef, M, U, fit)
    M = mult_columns()
    U = unk_columns()
    coef, fit, A = global_solve(M, U)
    # local amplitudes consistent with the final positions (no further movement)
    roi_refine(coef, M, U, fit, max_nfev=0)
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
        _log.warning("no standard errors for this spectrum's whole-compound fit", exc_info=True)
        se_m = np.full(nM, np.nan)
    for k, ii in MPk.items():
        if not ii:
            continue
        comp_fit[:, k] = M[:, ii] @ am[ii]
        npr = np.array([MP[i]["nprot"] for i in ii])
        amp[k] = float((am[ii] * npr).sum() / npr.sum())
        if np.isfinite(se_m[ii]).any():
            se[k] = float(np.sqrt(np.nansum((se_m[ii] * npr) ** 2)) / npr.sum())
    Mf = mult_columns(full=True) if dropped else M        # shapes of left-out compounds, for the upper bounds
    for i, m in enumerate(MP):
        m["global_amp"] = float(am[i])
        m["unit_peak"] = float(Mf[:, i].max())
        m["apex_ppm"] = float(x[int(np.argmax(Mf[:, i]))])
        m["dropped"] = dropped.get(m["k"], "")
        m["pinned"] = bool(abs(Dk[m["k"]]) >= DMAX[m["k"]] - 0.3 / sf) and lib[m["k"]]["source"] != "polymer"
        m["compound_shift_Hz"] = float(Dk[m["k"]] * sf)
    U = unk_columns()
    return dict(tsp=tsp, sigma=sigma, fw0=fw0, eta=eta, g_off=g_off, ppm=ppm, grid=grid, yv=yv, wts=wts,
                fit=fit, base_fit=base_fit, amp=amp, se=se, MP=MP, comp_fit=comp_fit, unknown=unknown, uamp=ua,
                unk_fit=U * ua[None, :] if len(unknown) else np.zeros((len(x), 0)),
                chi2n=chi(fit), sf=sf, mcols=None, ph=ph, noise_method="from point-to-point differences"
                if noise_why else "signal-free regions")


# =============================================================================
# 6. Grading, peak lists
# =============================================================================
def _observed_height(x, yv, base, apex, half, dx):
    """Largest observed signal (above the fitted baseline) within +/-half of apex, or None if that window
    touches an excluded region (water, PEO core) or the apex itself lies in one."""
    j0, j1 = np.searchsorted(x, apex - half), np.searchsorted(x, apex + half)
    seg = x[max(j0 - 1, 0):min(j1 + 1, len(x))]
    if seg.size < 2 or np.any(np.diff(seg) > 2 * dx) or j1 <= j0:
        return None
    return max(float((yv[j0:j1] - base[j0:j1]).max()), 0.0)


def _contradicted_by(i, mps, amp, x, yv, base, fw0, sf, sig, dx, cols=None, *, S):
    """Multiplets (numbers) of one molecule that `amp` (the amount implied by multiplet i, or the reported value)
    predicts above CONTRADICT_MIN_SNR x noise but that are (almost) absent from the data."""
    out = []
    nprot = sum(m["nprot"] for m in mps)
    for o in mps:
        if o is i or (o["nprot"] < 0.3 and o["nprot"] / nprot < 0.1):
            continue
        col = cols[id(o)] if cols is not None else None
        up = float(col.max()) if col is not None else o.get("unit_peak", 0.0)
        if up <= 0:
            continue
        apex = float(x[int(np.argmax(col))]) if col is not None else o["apex_ppm"]
        hobs = _observed_height(x, yv, base, apex, 1.5 * fw0 * o["wf"] + (S["MULTIPLET_SHIFT_HZ"] + 1.0) / sf, dx)
        if hobs is None:
            continue
        hpred = amp * up
        if hpred > S["CONTRADICT_MIN_SNR"] * sig and hobs < S["CONTRADICT_FRACTION"] * hpred:
            out.append(o["j"] + 1)
    return out


def _n_lines(m):
    return int((m["it"] > 0.03 * m["it"].max()).sum())


def _excluded_fraction(m, *, S):
    """Fraction of a multiplet's intensity whose lines sit in an ignored region or outside the fitted range."""
    pos = np.asarray(m["p"], float) + float(m["d"])
    it = np.asarray(m["it"], float)
    lo, hi = S["FIT_RANGE"]
    out = (pos < lo) | (pos > hi)
    for a_, b_ in S["EXCLUDE"]:
        out |= (pos > a_) & (pos < b_)
    return float(it[out].sum() / it.sum()) if it.sum() > 0 else 1.0


def detect_pluronic(ppm, y, *, S):
    """True if the spectrum has the tall PEO line of Pluronic F-68 near 3.71 ppm."""
    m = (ppm > 3.69) & (ppm < 3.74)
    return bool(m.any() and y[m].max() > S["PEO_DETECT_SNR"] * noise_sigma(ppm, y))
