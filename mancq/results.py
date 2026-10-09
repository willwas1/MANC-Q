"""mancq.results - read a finished MANC-Q output folder back in (used by the GUI's Results tab)."""
import os
import numpy as np
import pandas as pd

TIERS = ["Quantified", "Overlapped (semi-quantitative)", "Deconvolution estimate (low confidence)",
         "Manually adjusted", "Upper bound only", "Not detected", "Not measurable (region ignored)"]


def is_results_folder(out):
    return os.path.exists(os.path.join(out, "concentrations_long.csv"))


def load(out):
    """dict(long=DataFrame sample x metabolite rows, samples=[...], metabolites=[...] ordered by median value,
    qc=DataFrame)."""
    G = pd.read_csv(os.path.join(out, "concentrations_long.csv"), dtype={"sample": str})
    G = G[~G.metabolite.str.startswith("Pluronic")]
    Q = pd.read_csv(os.path.join(out, "qc_per_sample.csv"), dtype={"sample": str})
    samples = list(dict.fromkeys(G["sample"]))
    good = G[G.tier.isin(TIERS[:4])]
    order = good.groupby("metabolite").conc_mM.median().sort_values(ascending=False).index.tolist()
    rest = [m for m in dict.fromkeys(G.metabolite) if m not in order]
    return dict(long=G, qc=Q, samples=samples, metabolites=order + rest)


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
