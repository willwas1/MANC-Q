"""mancq.manual - manual (Chenomx-style) adjustment of fitted compounds, used by the Review step of the window and
by `python -m mancq apply-edits`.

A finished fit is rebuilt from _per_sample/<sample>_fit.npz: every multiplet of every compound with its fitted
amplitude, shift and linewidth, the smooth baseline, and the unassigned singlets (fit - baseline - compounds).
An adjusted compound is drawn as its whole simulated spectrum at one concentration (the library's intensity
ratios), moved by a shift and broadened by a linewidth factor chosen by the user, starting from the positions and
linewidths of the automatic fit.

The automatic fit itself is never changed. Adjustments are saved on their own in _per_sample/<sample>_manual.json
(plain JSON) and applied when the tables are written: the metabolite gets the tier "Manually adjusted", and its
automatic tier and value are kept alongside (engine.apply_manual). Fitting a sample again sets its adjustments
aside, because they refer to the earlier fit.
"""
import os
import json
import time
import getpass
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.ndimage import maximum_filter1d
from . import engine as q

WIDTH_GRID = (0.7, 0.85, 1.0, 1.15, 1.3, 1.6, 2.0)   # linewidth factors tried by best_fit
SHIFT_RANGE_HZ, SHIFT_STEP_HZ = 6.0, 0.25            # compound shifts tried by best_fit, around the current one


def stamp(edit, reason=None):
    """Record who made an adjustment and when."""
    edit = dict(edit)
    if reason:
        edit["reason"] = reason
    edit["user"] = getpass.getuser()
    edit["time"] = time.strftime("%Y-%m-%d %H:%M")
    return edit


def save_edits(out, smp, edits):
    """Write (or, with no edits, delete) _per_sample/<sample>_manual.json."""
    f = q.manual_file(out, smp)
    if not edits:
        if os.path.exists(f):
            os.remove(f)
        return
    json.dump(dict(mancq_version=q.__version__, sample=smp,
                   about="Manual adjustments made in the MANC-Q Review step. conc_mM: whole compound at this "
                         "concentration (mM in the original sample); shift_Hz: move of the whole compound from "
                         "the automatic fit; mult_shift_Hz: extra move of single multiplets; width: linewidth "
                         "factor relative to the automatic fit.",
                   edits=edits), open(f, "w"), indent=1)


class SampleFit:
    """One fitted spectrum rebuilt from its saved fit state, with manual adjustments applied on top."""

    def __init__(self, out, smp, library_dir=None):
        per = os.path.join(out, "_per_sample")
        f = os.path.join(per, f"{smp}_fit.npz")
        if not os.path.exists(f):
            raise FileNotFoundError(f"no saved fit for sample {smp} in {per}")
        d = np.load(f, allow_pickle=False)
        if "mp_json" not in d:
            raise ValueError(f"sample {smp} was fitted by MANC-Q 1.1.0 or earlier; fit it again with this version "
                             "(untick re-use or use a new results folder) to review it")
        self.out, self.smp = out, smp
        self.x = np.asarray(d["x"], float)
        self.yv = np.asarray(d["yv"], float)
        self.fit_auto = np.asarray(d["fit"], float)
        self.base = np.asarray(d["base"], float)
        self.sigma, self.fw0, self.sf, self.cf = float(d["sigma"]), float(d["fw0"]), float(d["sf"]), float(d["cf"])
        self.dx = float(d["dx"])
        pl = bool(d["pluronic"]) if "pluronic" in d else True
        libdir = library_dir or os.path.join(q.HERE, "library")
        self.lib = q.prepare_library(q.build_library(self.sf, libdir, os.path.join(out, "_cache")), self.sf, pl)
        self.names = [c["name"] for c in self.lib]
        self.index = {n: k for k, n in enumerate(self.names)}
        self.MP = []
        for m in json.loads(str(d["mp_json"])):
            m = dict(m)
            m["p"], m["it"] = np.asarray(m["p"], float), np.asarray(m["it"], float)
            self.MP.append(m)
        self.MPk = {k: [i for i, m in enumerate(self.MP) if m["k"] == k] for k in range(len(self.lib))}
        self.eta = float(d["eta"]) if "eta" in d else None
        lm = maximum_filter1d(np.abs(self.yv), size=max(3, int(q.S["TAIL_WINDOW_PPM"] / self.dx)))
        self.wts = 1.0 / np.sqrt(self.sigma ** 2 + (q.S["REL_MODEL_ERROR"] * np.abs(self.yv)) ** 2 +
                                 (q.S["TAIL_ERROR"] * lm) ** 2)
        if self.eta is None:
            self.eta = self._estimate_eta()
        self.auto = {}
        for k in range(len(self.lib)):
            c = self._auto_contribution(k)
            if c is not None:
                self.auto[k] = c
        self.unk = self.fit_auto - self.base - (sum(self.auto.values()) if self.auto else 0)
        g = pd.read_csv(os.path.join(per, f"{smp}_grades.csv"), dtype={"sample": str})
        self.grades = g.set_index("metabolite")
        self.edits = q.read_manual(out, smp)

    # ---------------------------------------------------------------- building blocks
    def _auto_contribution(self, k, eta=None):
        eta = self.eta if eta is None else eta
        out, any_ = np.zeros(len(self.x)), False
        for i in self.MPk[k]:
            m = self.MP[i]
            a = float(m.get("global_amp", 0.0))
            if a > 0:
                out += a * q.mcol(self.x, m["p"], m["it"], m["d"], self.fw0 * m["wf"], eta)
                any_ = True
        return out if any_ else None

    def _estimate_eta(self):
        """Results from MANC-Q 1.1.x do not store the lineshape: pick the one that leaves the unassigned part
        non-negative (as it is by construction)."""
        best = None
        for eta in np.arange(0.3, 1.0001, 0.05):
            comp = sum((c for c in (self._auto_contribution(k, eta) for k in range(len(self.lib))) if c is not None),
                       np.zeros(len(self.x)))
            unk = self.fit_auto - self.base - comp
            score = float((np.clip(unk, None, 0) ** 2).sum())
            if best is None or score < best[0]:
                best = (score, float(eta))
        return best[1]

    def contribution(self, name, edit=None, x=None):
        """Simulated spectrum of one compound: automatic (edit None) or adjusted (edit dict)."""
        k = self.index[name]
        if edit is None:
            if x is None:
                return self.auto.get(k, np.zeros(len(self.x)))
            c = np.zeros(len(x))
            for i in self.MPk[k]:
                m = self.MP[i]
                a = float(m.get("global_amp", 0.0))
                if a > 0:
                    c += a * q.mcol(x, m["p"], m["it"], m["d"], self.fw0 * m["wf"], self.eta)
            return c
        return float(edit.get("conc_mM", 0.0)) / self.cf * self.unit(name, edit, x)

    def unit(self, name, edit, x=None):
        """The adjusted shape of a compound at a per-proton amplitude of 1."""
        x = self.x if x is None else x
        k = self.index[name]
        out = np.zeros(len(x))
        sh = float(edit.get("shift_Hz", 0.0))
        wfac = float(edit.get("width", 1.0))
        ms = edit.get("mult_shift_Hz", {}) or {}
        for i in self.MPk[k]:
            m = self.MP[i]
            d = m["d"] + (sh + float(ms.get(str(m["j"] + 1), 0.0))) / self.sf
            out += q.mcol(x, m["p"], m["it"], d, self.fw0 * m["wf"] * wfac, self.eta)
        return out

    def current(self, name, edits=None):
        edits = self.edits if edits is None else edits
        return self.contribution(name, edits.get(name))

    def total(self, edits=None, without=()):
        """Baseline + unassigned singlets + every compound (adjusted where edited)."""
        edits = self.edits if edits is None else edits
        t = self.base + self.unk
        for k, c in self.auto.items():
            n = self.names[k]
            if n not in edits and n not in without:
                t = t + c
        for n, e in edits.items():
            if n in self.index and n not in without:
                t = t + self.contribution(n, e)
        return t

    def auto_value(self, name):
        """(tier, value in mM) from the automatic grading."""
        if name not in self.grades.index:
            return "", np.nan
        r = self.grades.loc[name]
        return str(r.tier), float(r.conc_mM) if np.isfinite(r.conc_mM) else np.nan

    def measurable(self, name):
        return self.auto_value(name)[0] != q.NOTMEAS

    def start_edit(self, name):
        """A new adjustment starts from the automatic value, at the automatic positions and linewidths."""
        if name in self.edits:
            return dict(self.edits[name])
        tier, v = self.auto_value(name)
        return dict(conc_mM=float(v) if np.isfinite(v) else 0.0, shift_Hz=0.0, width=1.0, mult_shift_Hz={},
                    note="", reason="set by hand")

    def multiplets(self, name, edit=None):
        """[(number, centre ppm, protons)] of a compound, at its current position."""
        k = self.index[name]
        e = edit or {}
        sh = float(e.get("shift_Hz", 0.0))
        ms = e.get("mult_shift_Hz", {}) or {}
        out = []
        for i in self.MPk[k]:
            m = self.MP[i]
            c = float((m["p"] * m["it"]).sum() / m["it"].sum() + m["d"])
            c += (sh + float(ms.get(str(m["j"] + 1), 0.0))) / self.sf
            out.append((m["j"] + 1, c, float(m["nprot"])))
        return sorted(out, key=lambda t: t[1])

    def window(self, name, edit=None, pad_hz=4.0):
        """Points near the compound's multiplets (those carrying a real share of its protons)."""
        k = self.index[name]
        e = edit or {}
        sh = float(e.get("shift_Hz", 0.0))
        wfac = float(e.get("width", 1.0))
        ms = e.get("mult_shift_Hz", {}) or {}
        mask = np.zeros(len(self.x), bool)
        nprot = self.lib[k]["nprot"]
        for i in self.MPk[k]:
            m = self.MP[i]
            if m["nprot"] < 0.3 and m["nprot"] / nprot < 0.1:
                continue
            d = m["d"] + (sh + float(ms.get(str(m["j"] + 1), 0.0))) / self.sf
            fw = self.fw0 * m["wf"] * wfac
            lo = m["p"].min() + d - 3 * fw - pad_hz / self.sf
            hi = m["p"].max() + d + 3 * fw + pad_hz / self.sf
            mask |= (self.x >= lo) & (self.x <= hi)
        return mask

    def misfit(self, name, total=None, edits=None):
        """How badly the fit matches the data around a compound's peaks: sum |measured - fit| there, as a share of
        the compound's own signal (0 = perfect; lower is better). None if the compound has no real signal."""
        edits = self.edits if edits is None else edits
        total = self.total(edits) if total is None else total
        mask = self.window(name, edits.get(name))
        own = self.current(name, edits)[mask]
        if own.size == 0 or own.max() < 3 * self.sigma:
            return None
        return float(np.abs(self.yv[mask] - total[mask]).sum() / own.sum())

    def chi2n(self, edits=None):
        f = self.total(edits)
        return float((((self.yv - f) * self.wts) ** 2).sum() / len(self.yv))

    # ---------------------------------------------------------------- fitting helpers
    def best_fit(self, name, edits=None):
        """Concentration, shift and linewidth of one compound that best explain the data around its peaks, with
        everything else as it is now (Chenomx's "fit this compound")."""
        edits = self.edits if edits is None else edits
        e0 = dict(edits.get(name) or self.start_edit(name))
        r_all = self.yv - self.total(edits, without=(name,))
        mask = self.window(name, e0, pad_hz=SHIFT_RANGE_HZ + 4.0)
        if mask.sum() < 5:
            return e0
        xs, r, w2 = self.x[mask], r_all[mask], self.wts[mask] ** 2
        best = None
        s0 = float(e0.get("shift_Hz", 0.0))
        for sh in np.arange(s0 - SHIFT_RANGE_HZ, s0 + SHIFT_RANGE_HZ + 1e-9, SHIFT_STEP_HZ):
            for wf in WIDTH_GRID:
                col = self.unit(name, dict(e0, shift_Hz=sh, width=wf), xs)
                den = (w2 * col * col).sum()
                if den <= 0:
                    continue
                a = max(0.0, float((w2 * col * r).sum() / den))
                ssr = float((w2 * (r - a * col) ** 2).sum()) * (1.0 + q.S["SHIFT_PRIOR"] * (sh / 10.0) ** 2)
                if best is None or ssr < best[0]:
                    best = (ssr, a, float(sh), float(wf))
        if best is None:
            return e0
        _, a, sh, wf = best
        return dict(e0, conc_mM=a * self.cf, shift_Hz=round(sh, 3), width=wf)

    def neighbours(self, name, edits=None, min_snr=3.0):
        """Other compounds with real signal under this compound's peaks."""
        edits = self.edits if edits is None else edits
        mask = self.window(name, edits.get(name))
        out = []
        for n in set(self.names[k] for k in self.auto) | set(edits):
            if n == name or n not in self.index:
                continue
            c = self.current(n, edits)[mask]
            if c.size and c.max() > min_snr * self.sigma:
                out.append(n)
        return sorted(out)

    def refit_neighbours(self, name, edits=None):
        """Keep this compound as it is and refit the concentrations of the compounds overlapping it (each as its
        whole simulated spectrum). Returns {metabolite: new edit} for the neighbours."""
        edits = self.edits if edits is None else edits
        nb = self.neighbours(name, edits)
        if not nb:
            return {}
        shapes = {n: dict(edits.get(n) or self.start_edit(n)) for n in nb}
        mask = self.window(name, edits.get(name))
        for n in nb:
            mask |= self.window(n, shapes[n])
        fixed = self.total(edits, without=tuple(nb))
        r = (self.yv - fixed)[mask]
        w = self.wts[mask]
        A = np.array([self.unit(n, shapes[n], self.x[mask]) for n in nb]).T
        a, _ = nnls(A * w[:, None], r * w, maxiter=50 * A.shape[1])
        return {n: dict(shapes[n], conc_mM=float(a[i]) * self.cf, reason=f"refitted around {name}")
                for i, n in enumerate(nb)}

    # ---------------------------------------------------------------- outputs
    def res_like(self, edits=None):
        """The pieces overlay_pdf needs, with the adjustments applied."""
        edits = self.edits if edits is None else edits
        comp = np.zeros((len(self.x), len(self.lib)))
        for k, c in self.auto.items():
            if self.names[k] not in edits:
                comp[:, k] = c
        for n, e in edits.items():
            if n in self.index:
                comp[:, self.index[n]] = self.contribution(n, e)

        class Grid:
            pass
        g = Grid()
        g.x, g.dx = self.x, self.dx
        fit = self.base + self.unk + comp.sum(1)
        return dict(grid=g, yv=self.yv, fit=fit, base_fit=self.base, comp_fit=comp, unk_fit=self.unk[:, None],
                    chi2n=float((((self.yv - fit) * self.wts) ** 2).sum() / len(self.yv)), sf=self.sf,
                    fw0=self.fw0, sigma=self.sigma)

    def write_outputs(self, overlay=True):
        """Overlay PDF and simulated-spectrum CSV of this sample, with its current adjustments."""
        res = self.res_like()
        cols = {"ppm": self.x, "observed": self.yv, "simulated_total": res["fit"], "baseline": self.base}
        for k, n in enumerate(self.names):
            if res["comp_fit"][:, k].any():
                cols[n] = res["comp_fit"][:, k]
        cols["unassigned_singlets"] = self.unk
        pd.DataFrame(cols).to_csv(os.path.join(self.out, f"simulated_spectrum_sample_{self.smp}.csv.gz"),
                                  index=False, float_format="%.6g")
        if overlay:
            gr = q.apply_manual(self.out, self.grades.reset_index())
            q.overlay_pdf(res, self.lib, gr, self.smp, os.path.join(self.out, f"overlay_sample_{self.smp}.pdf"))


def update_results(out, samples=None, log=print):
    """Rebuild the tables and plots of a finished run with the manual adjustments saved in _per_sample (nothing
    is refitted). samples: also redraw these samples' overlays (e.g. ones whose adjustments were just removed);
    samples with adjustments are always redrawn."""
    from .regrade import restore_settings
    keep = dict(q.S)
    try:
        restore_settings(out, grading_from_run=True)
        per = os.path.join(out, "_per_sample")
        qcf = os.path.join(out, "qc_per_sample.csv")
        order = list(pd.read_csv(qcf, dtype={"sample": str})["sample"]) if os.path.exists(qcf) else []
        order += sorted({f[:-len("_grades.csv")] for f in os.listdir(per) if f.endswith("_grades.csv")} - set(order),
                        key=lambda s: (len(s), s))
        done = [s for s in order if all(os.path.exists(os.path.join(per, f"{s}_{k}.csv"))
                                        for k in ("grades", "peaks", "obs", "qc"))]
        redraw = set(samples or []) | {s for s in done if os.path.exists(q.manual_file(out, s))}
        overlays = bool(q.S.get("OVERLAYS", True))
        for s in sorted(redraw, key=lambda s: (len(s), s)):
            if s in done and os.path.exists(os.path.join(per, f"{s}_fit.npz")):
                SampleFit(out, s).write_outputs(overlay=overlays)
                log(f"  {s}: overlay and simulated spectrum redrawn ({len(q.read_manual(out, s))} adjusted)")
        grades = [pd.read_csv(os.path.join(per, f"{s}_grades.csv"), dtype={"sample": str}) for s in done]
        peaks = [pd.read_csv(os.path.join(per, f"{s}_peaks.csv"), dtype={"sample": str}) for s in done]
        obs = [pd.read_csv(os.path.join(per, f"{s}_obs.csv"), dtype={"sample": str}) for s in done]
        qc = [pd.read_csv(os.path.join(per, f"{s}_qc.csv"), dtype={"sample": str}).iloc[0].to_dict() for s in done]
        sf = float(qc[0]["SF_MHz"])
        lib = q.prepare_library(q.build_library(sf, q.S["LIBRARY_DIR"] or os.path.join(q.HERE, "library"),
                                                os.path.join(out, "_cache")), sf, True)
        log(f"Applying manual adjustments to {out} ...")
        G = q.aggregate_outputs(out, grades, peaks, obs, qc, lib, done, log)
    finally:
        q.S.clear()
        q.S.update(keep)
    return G
