"""mancq.results - read a finished MANC-Q output folder back in: the tables (Results, Review and Figures steps)
and the saved fit of each sample (regrade, Review, apply-edits). Everything needed to read folders made by
older versions is here, and a folder made by a newer MANC-Q, or fitted with a different compound library, is
refused with a clear message instead of being misread."""
import json
import os
import numpy as np
import pandas as pd
from . import engine as q

TIERS = ["Quantified", "Overlapped (semi-quantitative)", "Deconvolution estimate (low confidence)",
         "Manually adjusted", "Upper bound only", "Not detected", "Not measurable (region ignored)"]


class ResultsTooNew(Exception):
    """The folder was written by a newer MANC-Q, in a layout this version does not know."""


class LibraryChanged(Exception):
    """The folder was fitted with a different compound library from the one installed now."""


def is_results_folder(out):
    return os.path.exists(os.path.join(out, "concentrations_long.csv"))


def _vtuple(v):
    try:
        return tuple(int(x) for x in str(v).split(".")[:3])
    except ValueError:
        return (0,)


def results_info(out):
    """dict(format, version, library) of a results folder: the layout version (1 = MANC-Q 1.9.2 or earlier), the
    MANC-Q version that wrote it and the fingerprint of its compound library (None before 1.8.1)."""
    f = os.path.join(out, "_per_sample", "run_settings.json")
    d = {}
    if os.path.exists(f):
        with open(f, encoding="utf-8") as fh:
            d = json.load(fh)
    return dict(format=int(d.get("RESULTS_FORMAT", 1)), version=str(d.get("MANCQ_VERSION", "1.1 or earlier")),
                library=d.get("LIBRARY_HASH"))


def check_readable(out):
    """results_info(out), or ResultsTooNew if a newer MANC-Q wrote the folder."""
    info = results_info(out)
    if info["format"] > q.RESULTS_FORMAT:
        raise ResultsTooNew(f"This results folder was made by MANC-Q {info['version']}, which is newer than this "
                            f"MANC-Q ({q.__version__}). Install MANC-Q {info['version']} or later to open it.")
    return info


def check_library(out, libdir=None):
    """results_info(out), or LibraryChanged if the folder was fitted with another compound library (regrading or
    adjusting it would then mix two libraries)."""
    info = check_readable(out)
    libdir = libdir or os.path.join(q.HERE, "library")
    if info["library"] and info["library"] != q.library_hash(libdir):
        raise LibraryChanged(
            f"This results folder was fitted with a different compound library (by MANC-Q {info['version']}) from "
            f"the one in this MANC-Q ({q.__version__}). Regrading or adjusting it now would mix the two. Fit the "
            f"spectra again into a new results folder, or use MANC-Q {info['version']} for this folder.")
    return info


def load_fit(out, smp, lib, log=print):
    """The saved fit state of one sample: (the arrays of _per_sample/<sample>_fit.npz, list of multiplet
    records). Each record's compound is checked against `lib` (results from 2.0 on store the name)."""
    f = os.path.join(out, "_per_sample", f"{smp}_fit.npz")
    if not os.path.exists(f):
        raise FileNotFoundError(f"no saved fit for sample {smp} in {os.path.dirname(f)}")
    d = np.load(f, allow_pickle=False)
    if "mp_json" in d:
        mp_list = json.loads(str(d["mp_json"]))
    else:
        # MANC-Q 1.1.0 or earlier stored this part with Python pickle, which can run code when loaded; only
        # accept those files from your own runs
        log(f"  {smp}: older results format; loading it as you created it yourself")
        mp_list = list(np.load(f, allow_pickle=True)["mp"])
    MP = []
    for m in mp_list:
        m = dict(m)
        m["p"], m["it"] = np.asarray(m["p"], float), np.asarray(m["it"], float)
        if "name" in m and (m["k"] >= len(lib) or lib[m["k"]]["name"] != m["name"]):
            raise LibraryChanged(f"sample {smp}: the saved fit has {m['name']} where the installed library has "
                                 f"{lib[m['k']]['name'] if m['k'] < len(lib) else 'no compound'}; this folder was "
                                 f"fitted with another compound library. Fit the spectra again.")
        MP.append(m)
    return d, MP


def load_results(out):
    """dict(long=DataFrame sample x metabolite rows, samples=[...], metabolites=[...] ordered by median value,
    qc=DataFrame, info=results_info). Results of older versions are brought up to date here: before 1.8.2 an
    upper bound was stored in conc_mM; before 1.9.2 there was no qc_flag."""
    info = check_readable(out)
    G = pd.read_csv(os.path.join(out, "concentrations_long.csv"), dtype={"sample": str})
    G = G[~G.metabolite.str.startswith("Pluronic")].copy()
    if _vtuple(info["version"]) < (1, 8, 2):
        ub = G.tier == "Upper bound only"
        if "upper_bound_mM" not in G:
            G["upper_bound_mM"] = np.nan
        G.loc[ub, "upper_bound_mM"] = G.loc[ub, "upper_bound_mM"].fillna(G.loc[ub, "conc_mM"])
        G.loc[ub, "conc_mM"] = np.nan
    Q = pd.read_csv(os.path.join(out, "qc_per_sample.csv"), dtype={"sample": str})
    if "qc_flag" not in Q:
        Q = q.qc_flags(Q, S=q.Settings.load(out))
    samples = list(dict.fromkeys(G["sample"]))
    good = G[G.tier.isin(TIERS[:4])]
    order = good.groupby("metabolite").conc_mM.median().sort_values(ascending=False).index.tolist()
    rest = [m for m in dict.fromkeys(G.metabolite) if m not in order]
    per_cell, per_met = reporter_peaks(out, G)
    return dict(long=G, qc=Q, samples=samples, metabolites=order + rest, reporter_cell=per_cell, reporter_met=per_met,
                info=info)


load = load_results


def reporter_peaks(out, G):
    """Where each value was measured: {(metabolite, sample): [ppm, ...]} for Quantified and Overlapped values (the
    fitted centre of the reporter multiplets), and {metabolite: "1.331" or "3.776, 3.752"}: the usual reporter
    peak(s), median over the samples (peaks used in at least a quarter of them). Results written before MANC-Q 1.5
    have no reporter_ppm column; it is then taken from the peak list."""
    use = G[G.tier.isin(TIERS[:2])]
    cell = {}
    if "reporter_ppm" in G:
        for r in use.itertuples():
            v = [float(t) for t in str(r.reporter_ppm).split(";") if t.strip() not in ("", "nan")]
            if v:
                cell[(r.metabolite, r.sample)] = v
    else:
        f = os.path.join(out, "peak_list_by_metabolite.csv")
        if os.path.exists(f):
            P = pd.read_csv(f, dtype={"sample": str})
            P = P[P.used_as_reporter.astype(bool)]
            ok = set(zip(use.metabolite, use["sample"]))
            for (m, s), sub in P.groupby(["metabolite", "sample"]):
                if (m, s) in ok:
                    cell[(m, s)] = sorted(sub.centre_ppm.astype(float), reverse=True)
    met = {}
    by_met = {}
    for (m, s), v in cell.items():
        by_met.setdefault(m, []).append(v)
    for m, lists in by_met.items():
        # peaks within 0.01 ppm of each other (the same peak in different samples) form one group
        groups, cur = [], []
        for p in sorted(p for v in lists for p in v):
            if cur and p - cur[-1] > 0.01:
                groups.append(cur); cur = []
            cur.append(p)
        if cur:
            groups.append(cur)
        keep = [np.median(v) for v in groups if len(v) >= max(1, 0.25 * len(lists))]
        met[m] = ", ".join(f"{p:.3f}" for p in sorted(keep, reverse=True))
    return cell, met


def metabolite_window(out, sample, metabolite, pad=0.012, min_half=0.018):
    """Arrays to plot one metabolite in one sample: the region around its main reporter multiplet.
    Returns dict(ppm, observed, fit, baseline, compound, lo, hi, row) or None if not available."""
    f = os.path.join(out, f"simulated_spectrum_sample_{sample}.csv.gz")
    pk = os.path.join(out, "_per_sample", f"{sample}_peaks.csv")
    if not (os.path.exists(f) and os.path.exists(pk)):
        return None
    P = pd.read_csv(pk)
    sub = P[P.metabolite == metabolite]
    if not len(sub):
        return None
    row = sub[sub.used_as_reporter] if sub.used_as_reporter.any() else sub
    row = row.sort_values(["height_snr", "n_protons"], ascending=False).iloc[0]
    half = max(min_half, (row.to_ppm - row.from_ppm) / 2 + pad)
    lo, hi = row.centre_ppm - half, row.centre_ppm + half
    d = pd.read_csv(f)
    w = (d.ppm > lo) & (d.ppm < hi)
    if w.sum() < 5:
        return None
    d = d[w]
    comp = d[metabolite].values if metabolite in d else np.zeros(len(d))
    return dict(ppm=d.ppm.values, observed=d.observed.values, fit=d.simulated_total.values,
                baseline=d.baseline.values, compound=comp, lo=lo, hi=hi, row=row)
