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
    per_cell, per_met = reporter_peaks(out, G)
    return dict(long=G, qc=Q, samples=samples, metabolites=order + rest, reporter_cell=per_cell, reporter_met=per_met)


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
