"""
mancq.report - the outputs: overlay, gallery and profile PDFs, the per-sample quality flag, and the
workbook and tables built from the per-sample results (aggregate_outputs).
"""
import json
import os
import time
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from .settings import DECONV, MANUAL, NOTMEAS, __version__
from .io import read_manual
from .grade import TIER_COLOURS, TIER_ORDER, apply_manual


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


# =============================================================================
# 7. Plots
# =============================================================================
PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#17becf",
           "#bcbd22", "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#e6550d"]
TAG = {"Quantified": "Q", "Overlapped (semi-quantitative)": "O", DECONV: "D", MANUAL: "M",
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


def overlay_pdf(res, lib, grades, sample, path, *, S):
    x = res["grid"].x
    names = [c["name"] for c in lib]
    tiers = dict(zip(grades.metabolite, grades.tier))
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
    ok = tier.isin(["Quantified", "Overlapped (semi-quantitative)", DECONV, MANUAL])
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
                                           (MANUAL, "P", None, "manually adjusted"),
                                           ("Upper bound only", "v", "none", "upper bound")]:
                    sel = t == style
                    if sel.any():
                        ax.scatter(x[sel], v[sel], marker=mk, s=18, facecolors=fc, color=TIER_COLOURS[style][0],
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


def _metabolite_fwhm(res, sf, *, S):
    """Median fitted linewidth (Hz) of the metabolite multiplets clearly above the noise."""
    w = [res["fw0"] * m["wf"] * sf for m in res["MP"] if m.get("loc_h", 0) > 10 * S["LOQ_SNR"] * res["sigma"]
         and m.get("loc_dom", 0) >= S["Q_MIN_DOMINANCE"]]
    return round(float(np.median(w)), 3) if w else np.nan


def qc_flags(Q, *, S):
    """Add qc_flag (ok / check / poor) and qc_reason (in words) to the per-sample QC table."""
    Q = Q.copy()
    flags, reasons = [], []
    num = lambda r, k: float(r[k]) if k in r and pd.notna(r[k]) and not isinstance(r[k], str) else np.nan
    for _, r in Q.iterrows():
        bad, chk = [], []

        def test(v, lim, text, low=False):
            if not np.isfinite(v):
                return
            c, p = lim
            if (v < p) if low else (v > p):
                bad.append(text)
            elif (v < c) if low else (v > c):
                chk.append(text)
        if "error" in r and isinstance(r["error"], str) and r["error"]:
            flags.append("poor"); reasons.append(f"not fitted: {r['error']}"); continue
        fw = num(r, "TSP_fwhm_Hz")
        test(fw, S["QC_TSP_FWHM_HZ"], f"TSP line {fw:.1f} Hz wide (shimming?)")
        mfw = num(r, "metabolite_fwhm_Hz")
        if np.isfinite(fw) and np.isfinite(mfw) and fw > S["QC_TSP_BROAD_RATIO"] * mfw:
            chk.append(f"TSP line ({fw:.1f} Hz) broader than the metabolites ({mfw:.1f} Hz): TSP may be bound to "
                       f"protein, which makes every concentration too high")
        test(num(r, "chi2_per_point"), S["QC_CHI2"], f"poor fit (chi2 per point {num(r, 'chi2_per_point'):.1f})")
        ex = num(r, "signal_explained_by_library")
        test(ex, S["QC_EXPLAINED"], f"only {ex:.0%} of the signal explained by the library", low=True)
        un = num(r, "unassigned_singlets")
        test(un, S["QC_UNASSIGNED"], f"{un:.0f} unassigned singlets (unknown compounds or a poor fit)")
        ph = num(r, "TSP_phase_deg")
        test(abs(ph), S["QC_PHASE_DEG"], f"residual phase error about {abs(ph):.1f} degrees at TSP: re-phase this "
                                         f"spectrum (step 2 of the window, or TopSpin)")
        if str(r.get("pH_flag", "")) == "protons_disagree":
            chk.append("the two histidine pH estimates disagree")
        d1, aq = num(r, "D1_s"), num(r, "AQ_s")
        if np.isfinite(d1) and np.isfinite(aq) and d1 + aq < S["QC_MIN_RELAXATION_S"]:
            chk.append(f"D1 + AQ = {d1 + aq:.1f} s (under {S['QC_MIN_RELAXATION_S']:g} s): signals with long T1 "
                       f"(e.g. formate, aromatic protons) may read low")
        if str(r.get("noise_method", "")).startswith("from point"):
            chk.append("noise estimated without a signal-free region (check LOD values)")
        flags.append("poor" if bad else "check" if chk else "ok")
        reasons.append("; ".join(bad + chk))
    Q["qc_flag"], Q["qc_reason"] = flags, reasons
    return Q


NOT_DETECTED_NOTE = ("Metabolites graded Not detected are recorded as 0 mM. A zero means the signal was below this "
                     "sample's detection limit (LOD_mM). It may reflect true absence, a concentration too low to "
                     "detect, or a peak MANC-Q could not resolve in that spectrum. It is not proof that the "
                     "metabolite is absent.")
UPPER_BOUND_NOTE = ("Upper bound only values are left blank here: they are not measurements. The bound is in "
                    "Report_mM (as <=value) and in the upper_bound_mM column of Long_table.")


def ph_readme(*, S):
    """The paragraph on the pH estimate for the results README."""
    return (f"pH_est in the QC sheet is an ESTIMATE of the pH in the NMR tube when the spectrum was recorded (not in the "
            f"culture), from the positions of the histidine imidazole singlets H2 (looked for in "
            f"{S['PH_HIS_H2_WINDOW'][0]:.2f}-{S['PH_HIS_H2_WINDOW'][1]:.2f} ppm) and H4 "
            f"({S['PH_HIS_H4_WINDOW'][0]:.2f}-{S['PH_HIS_H4_WINDOW'][1]:.2f} ppm), TSP-referenced. Each peak top is "
            f"refined with a 3-point parabola and used only if at least {S['PH_MIN_SNR']} x noise. Henderson-Hasselbalch: "
            f"f = (d_HA - d_obs)/(d_HA - d_A), pH = pKa + log10(f/(1-f)), with approximate literature constants for free "
            f"histidine in water at 25 C: pKa {S['PH_HIS_PKA']}, H2 d_HA {S['PH_HIS_H2_LIMITS'][0]} / d_A "
            f"{S['PH_HIS_H2_LIMITS'][1]} ppm, H4 d_HA {S['PH_HIS_H4_LIMITS'][0]} / d_A {S['PH_HIS_H4_LIMITS'][1]} ppm. "
            f"pH_est is the mean of pH_H2 and pH_H4. pH_flag: ok, one_proton_only, protons_disagree (more than "
            f"{S['PH_MAX_DISAGREE']} apart) or no_histidine. Changes in pH between samples are reliable; the absolute "
            f"value is only good to about +/-0.2 pH units, because the constants were not measured in your medium "
            f"(H2 and H4 can then read a steady amount apart). "
            + (f"Histidine H2/H4 were fitted within {S['PH_HIS_PLACE_HZ']:g} Hz of the peaks used for the estimate."
               if S['PH_PLACE_HISTIDINE'] else ""))


def aggregate_outputs(out, grades, peaks, obs, qc, lib, samples, log, *, S):
    """Build every table and figure from the per-sample results (also used by regrade.py and after manual
    adjustments). Manual adjustments saved in _per_sample are applied here."""
    grades = [apply_manual(out, g) for g in grades]
    G = pd.concat(grades, ignore_index=True)
    edits = []
    for g in grades:
        smp = str(g["sample"].iloc[0]) if len(g) else ""
        for met, e in read_manual(out, smp).items():
            row = g[g.metabolite == met]
            edits.append(dict(sample=smp, metabolite=met,
                              automatic_tier=row.auto_tier.iloc[0] if len(row) else "",
                              automatic_mM=row.auto_conc_mM.iloc[0] if len(row) else np.nan,
                              manual_mM=e.get("conc_mM"), shift_Hz=e.get("shift_Hz", 0.0),
                              multiplet_shifts_Hz=json.dumps(e.get("mult_shift_Hz", {})),
                              linewidth_factor=e.get("width", 1.0), how=e.get("reason", "set by hand"),
                              note=e.get("note", ""), user=e.get("user", ""), time=e.get("time", "")))
    if not edits:
        G = G.drop(columns=["auto_tier", "auto_conc_mM", "manual_note"], errors="ignore")
    P = pd.concat(peaks, ignore_index=True)
    O = pd.concat(obs, ignore_index=True)
    Q = qc_flags(pd.DataFrame(qc), S=S)

    # ---------- summary tables ----------
    G = G[~G.metabolite.str.startswith("Pluronic")].copy()
    best_tier = G.groupby("metabolite").tier.agg(lambda t: min(t, key=lambda v: TIER_ORDER[v]))
    conc = G.pivot(index="metabolite", columns="sample", values="conc_mM")
    tier = G.pivot(index="metabolite", columns="sample", values="tier")
    se = G.pivot(index="metabolite", columns="sample", values="se_mM")
    lod = G.pivot(index="metabolite", columns="sample", values="LOD_mM")
    ubd = G.pivot(index="metabolite", columns="sample", values="upper_bound_mM")
    report = conc.copy().astype(object)
    for m in report.index:
        for s_ in report.columns:
            t = tier.loc[m, s_]
            if t == NOTMEAS:
                report.loc[m, s_] = "n/m"
            elif t == "Not detected":
                report.loc[m, s_] = f"<{lod.loc[m, s_]:.2g}"
            elif t == "Upper bound only":
                report.loc[m, s_] = f"<={ubd.loc[m, s_]:.2g}"
            elif t == DECONV:
                report.loc[m, s_] = f"~{conc.loc[m, s_]:.3g}"
            elif t == MANUAL:
                report.loc[m, s_] = f"{conc.loc[m, s_]:.4g} (manual)"
            else:
                report.loc[m, s_] = round(float(conc.loc[m, s_]), 4)
    summ = pd.DataFrame(dict(best_tier=best_tier))
    summ["class"] = G.groupby("metabolite").cls.first()
    for lab, cond in [("n_quantified", G.tier == "Quantified"), ("n_overlapped", G.tier.str.startswith("Overlapped")),
                      ("n_manually_adjusted", G.tier == MANUAL),
                      ("n_detected_any", ~G.tier.isin(["Not detected", NOTMEAS])),
                      ("n_measured", G.conc_mM.notna() & (G.tier != "Not detected")),
                      ("n_not_detected", G.tier == "Not detected")]:
        summ[lab] = G[cond].groupby("metabolite").size()
    summ = summ.fillna({k: 0 for k in ("n_quantified", "n_overlapped", "n_manually_adjusted", "n_detected_any",
                                       "n_measured", "n_not_detected")})
    # over measured values and zeros (Not detected) only; upper bounds are blank in conc and never averaged
    summ["mean_conc_mM"] = conc.mean(axis=1)
    summ["mean_rel_SE_%"] = (se / conc.replace(0, np.nan) * 100).mean(axis=1).round(1)
    summ["library_source"] = G.groupby("metabolite").source.first()
    summ["_o"] = summ.best_tier.map(TIER_ORDER)
    summ = summ.sort_values(["_o", "mean_conc_mM"], ascending=[True, False]).drop(columns="_o")
    idx = summ.index
    conc, report, tier, se, lod = conc.loc[idx], report.loc[idx], tier.loc[idx], se.loc[idx], lod.loc[idx]
    ubd = ubd.loc[idx]

    readme = pd.DataFrame(dict(item=[
        "Software", "Units", "Reference", "How numbers are produced", "Quantified", "Overlapped (semi-quantitative)",
        DECONV, MANUAL, "Upper bound only", "Not detected", "SE", "Library", "Peak_list_by_metabolite",
        "Observed_peaks_assigned", "QC", "pH estimate", "Manual_edits", "Shift_assignments", "Settings"], text=[
        f"MANC-Q {__version__} (https://github.com/willwas1/MANC-Q), results written {time.strftime('%Y-%m-%d %H:%M')}.",
        "mM in the original medium = tube mM x dilution factor (from the dilution file).",
        f"{S['REFERENCE']} internal standard only: {S['TSP_MM_IN_TUBE']} mM in tube, 9 H, area from a pseudo-Voigt fit incl. 29Si satellites. No calibration against any sample or standard.",
        "Whole spectrum fitted as a sum of simulated metabolite spectra + smooth baseline + free unassigned singlets; each multiplet has a small fitted shift and linewidth. Concentration = mean of the cleanly resolved 'reporter' multiplets (each fitted on its own); if none, the whole-compound fit. The reporter_ppm column of Long_table gives the fitted centre (ppm) of the reporter multiplet(s) a Quantified or Overlapped value comes from. " + ("A compound that only sits on other compounds' peaks (the amount its signal implies would put another of its multiplets clearly above the noise where the data, or, for a compound with no dominant multiplet, the fit, show almost nothing) is left out of the fit, so it cannot take signal from the compounds that are there; its Long_table method says why. " if S["DROP_GHOSTS"] else "") + "A multiplet whose own fitted amount predicts the molecule's other multiplets clearly above the noise where they are missing is never used as a reporter. check_note in Long_table warns when a compound graded without a value is fitted on this compound's peak.",
        f"Clean reporter multiplet(s): SNR>={S['LOQ_SNR']}, this metabolite >={S['Q_MIN_DOMINANCE']:.0%} of the fitted signal under it, |observed-fit| <= {S['Q_MAX_MISFIT']} x its own signal, shift not at bound; either two or more reporters agreeing within x{S['MAX_REPORTER_RATIO']}, or one reporter carrying >={S['SINGLE_REPORTER_MIN_SHARE']:.0%} of the protons; no multiplet of the molecule contradicted by the data; the molecule's other multiplets above LOQ (whatever their dominance) give a median value within {S['MAX_REPORTER_DISAGREE']:.0%} of the reporter value" + ("; a lone singlet reporter must be confirmed by another multiplet of the molecule within that margin (molecules with no other multiplet, such as acetate or formate, are exempt)" if S["SINGLET_NEEDS_SECOND"] else "") + ". Otherwise the value is graded Overlapped and its method says why.",
        f"SNR>={S['LOQ_SNR']} but the best multiplet is shared (dominance {S['O_MIN_DOMINANCE']:.0%}-{S['Q_MIN_DOMINANCE']:.0%}) or misfit <= {S['O_MAX_MISFIT']}. Number depends on the overlap partners: use for trends.",
        f"No multiplet of this metabolite is resolved enough to quantify on its own, but the whole-compound least-squares fit over >={S['DECONV_MIN_MULTIPLETS']} fitted multiplets (each above LOQ, at least {S['DECONV_MIN_DOMINANCE']:.0%} of the fitted signal under it, and with its own signal at least {1 / S['DECONV_MAX_MISFIT']:g} x the local misfit) is well determined (relative SE <= {S['DECONV_MAX_REL_SE']:.0%}) and sits below the data upper bound. The value is that deconvolution estimate, shown with a ~ in Report_mM. It depends on the library explaining the overlapping signals correctly, so treat it as an estimate, not a measurement.",
        "The value was set by hand in the Review step of the window (Chenomx-style: the whole simulated compound at one concentration, with a shift and linewidth chosen by the user), or refitted around a compound that was. The automatic tier and value are in the auto_tier and auto_conc_mM columns of Long_table, and every change, with who made it and when, is in Manual_edits. Report_mM marks these values with (manual). No standard error is given.",
        f"Signal is present where the metabolite would appear but it cannot be attributed to it with confidence (no usable reporter, or a multiplet the library predicts >{S['CONTRADICT_MIN_SNR']}x noise is missing" + (f", or a value from a single reporter multiplet whose other multiplets, predicted >{S['CONTRADICT_MIN_SNR']}x noise, get less than {S['FIT_SUPPORT_FRACTION']:.0%} of it in the fit (their signal belongs to other compounds)" if S["NEEDS_FIT_SUPPORT"] else "") + f"). The largest concentration the data allow (upper bound) is in Report_mM as <=value and in upper_bound_mM of Long_table; conc_mM and Concentration_mM_all are blank, and the bound is never included in mean_conc_mM.",
        "Below the detection limit; Report_mM shows <LOD in mM. 'n/m' = not measurable: every peak of the metabolite lies in a region left out of the fit. " + NOT_DETECTED_NOTE,
        "Standard error from the least-squares fit only (noise + local misfit). Excludes library, relaxation (T1) and TSP weighing errors.",
        f"{sum(c['source']=='GISSMO' for c in lib)} GISSMO spin systems (BMRB, QM-simulated at the spectrometer frequency) + {sum(c['source']=='SPIN' for c in lib)} other spin system(s) simulated the same way + {sum(c['source']=='CASMDB' for c in lib)} fixed peak lists (CASMDB) + optional Pluronic F-68 terms. See Library sheet.",
        "Every multiplet of every metabolite in every sample: fitted centre, range, line positions and relative intensities, protons, linewidth, shift from library, dominance, misfit, concentration from that multiplet alone, and whether it was used as a reporter.",
        "Peaks picked in the observed spectrum (>10 x noise) with the fitted contributors at each apex (the peak assignment list).",
        "Per-sample referencing, linewidth, noise, fit quality, fraction of signal explained, how the spectrum was processed, the estimated pH, the acquisition parameters from acqus (D1, AQ, NS, RG, PULPROG, TE, P1), the residual phase error measured on TSP (TSP_phase_deg) and the median metabolite linewidth. qc_flag sums these up: ok, check or poor, with the reasons in qc_reason (TSP linewidth, TSP much broader than the metabolites, chi2, signal explained, unassigned singlets, phase error, histidine pH flag, D1 + AQ under " + f"{S['QC_MIN_RELAXATION_S']:g}" + " s, noise estimate).",
        ph_readme(S=S) if "pH_est" in Q else "Not estimated: these results were fitted by a MANC-Q version before 1.4.",
        "Every manual adjustment: sample, metabolite, automatic and manual value, shift and linewidth used, note, user and time. Empty if nothing was adjusted.",
        "Chemical shift assignment table for publication: every metabolite with a value in the tiers SHIFT_TABLE_TIERS in at least SHIFT_TABLE_MIN_FRACTION of the samples (default Quantified, Overlapped or Manually adjusted in 25 %), its 1H shifts (fitted multiplet centres, median over those samples), * = peaks used to quantify it, ** = used but shared with another metabolite (semi-quantitative), the number of samples and the reference spectrum. The same table, formatted for a supplement, is in chemical_shift_assignments.xlsx.",
        json.dumps({k: v for k, v in S.items() if k != "OVERLAY_REGIONS"})]))
    libdf = pd.DataFrame([dict(metabolite=c["name"], cls=c["cls"], source=c["source"], n_protons=c["nprot"],
                               titratable=c["titratable"], n_lines=len(c["ppm"]), n_multiplets=len(c["mults"]),
                               provenance=c["provenance"]) for c in lib])
    # chemical shift assignment table for publication (also written on its own as chemical_shift_assignments.xlsx)
    from . import shifttable
    shifts, shift_note = shifttable.make(out, G, P, Q, lib, __version__, S["REFERENCE"],
                                         tuple(S["SHIFT_TABLE_TIERS"]), float(S["SHIFT_TABLE_MIN_FRACTION"]))
    shift_sheet = pd.concat([shifts[shifts.columns[:4]],
                             pd.DataFrame({"Metabolite": ["", shift_note]})], ignore_index=True)
    xl = os.path.join(out, "metabolite_concentrations.xlsx")
    with pd.ExcelWriter(xl, engine="openpyxl") as xw:
        readme.to_excel(xw, sheet_name="README", index=False)
        summ.to_excel(xw, sheet_name="Summary")
        report.to_excel(xw, sheet_name="Report_mM")
        conc.to_excel(xw, sheet_name="Concentration_mM_all", startrow=3)
        ws = xw.sheets["Concentration_mM_all"]
        ws["A1"], ws["A2"] = NOT_DETECTED_NOTE, UPPER_BOUND_NOTE
        se.to_excel(xw, sheet_name="SE_mM")
        tier.to_excel(xw, sheet_name="Tier")
        lod.to_excel(xw, sheet_name="LOD_mM")
        G.to_excel(xw, sheet_name="Long_table", index=False)
        P.to_excel(xw, sheet_name="Peak_list_by_metabolite", index=False)
        O.to_excel(xw, sheet_name="Observed_peaks_assigned", index=False)
        Q.to_excel(xw, sheet_name="QC", index=False)
        (pd.DataFrame(edits) if edits else pd.DataFrame(dict(note=["No manual adjustments."]))).to_excel(
            xw, sheet_name="Manual_edits", index=False)
        libdf.to_excel(xw, sheet_name="Library", index=False)
        shift_sheet.to_excel(xw, sheet_name="Shift_assignments", index=False)
        for ws in xw.book.worksheets:
            notes = ws.title == "Concentration_mM_all"          # two notes above the table
            ws.freeze_panes = "B5" if notes else "B2"
            for col in ws.columns:
                width = min(60, max(10, max(len(str(c.value)) if c.value is not None else 0
                                            for c in col[2 if notes else 0:50]) + 2))
                ws.column_dimensions[col[0].column_letter].width = width
    if len(samples) <= S["GALLERY_MAX_SAMPLES"]:
        gallery_from_csv(out, G, P, samples, os.path.join(out, "metabolite_gallery.pdf"))
    profiles_pdf(conc.where(tier != "Upper bound only", ubd), tier, os.path.join(out, "concentration_profiles.pdf"),
                 f"Concentration across samples ({os.path.basename(out)})")
    G.to_csv(os.path.join(out, "concentrations_long.csv"), index=False)
    report.to_csv(os.path.join(out, "concentrations_report_mM.csv"))
    P.to_csv(os.path.join(out, "peak_list_by_metabolite.csv"), index=False)
    O.to_csv(os.path.join(out, "observed_peaks_assigned.csv"), index=False)
    Q.to_csv(os.path.join(out, "qc_per_sample.csv"), index=False)
    me = os.path.join(out, "manual_edits.csv")
    if edits:
        pd.DataFrame(edits).to_csv(me, index=False)
    elif os.path.exists(me):
        os.remove(me)
    log("\nBest tier per metabolite across samples:")
    log(summ.best_tier.value_counts().to_string())
    log(Q.round(3).T.to_string())
    if "pH_est" in Q:
        log("\npH estimate: " + ph_readme(S=S))
    log(f"\nWrote {xl}")
    return G
