"""
mancq.grade - grading every compound (Quantified ... Not detected) from the fit, the peak list, and
manual adjustments applied on top.
"""
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from .settings import DECONV, MANUAL, NOTMEAS
from .io import read_manual
from .fit import _contradicted_by, _excluded_fraction, _n_lines, mcol


def grade(res, lib, sf, cf, *, S):
    """cf converts a per-proton area to mM in the medium. Per-multiplet quality
    measures come from the local (ROI) fits, where each multiplet has its own
    amplitude, so a library intensity error elsewhere in the molecule cannot
    bias a clean reporter."""
    yv = res["yv"]
    sig = res["sigma"]
    xg, base, dx, fw0 = res["grid"].x, res["base_fit"], res["grid"].dx, res["fw0"]
    rows, prow = [], []
    notes = {}                    # {compound: [why another compound's fit on its peaks is suspect]}
    for k, c in enumerate(lib):
        a = res["amp"][k]
        mps = [m for m in res["MP"] if m["k"] == k]
        if mps and all(m["unit_peak"] <= 0 or _excluded_fraction(m, S=S) >= 0.9 for m in mps):
            # every multiplet of this compound lies in an excluded region: nothing to fit
            rows.append(dict(metabolite=c["name"], cls=c["cls"], tier=NOTMEAS, conc_mM=np.nan, se_mM=np.nan,
                             whole_compound_fit_mM=np.nan, method="all multiplets inside ignored regions", reporter_ppm="",
                             LOD_mM=np.nan, upper_bound_mM=np.nan, snr=0.0, best_dominance=np.nan,
                             best_misfit=np.nan, n_reporters=0, reporter_ratio=np.nan,
                             source=c["source"], provenance=c["provenance"]))
            continue
        unit_peak = max(m["unit_peak"] for m in mps) if mps else 0
        snr = a * unit_peak / sig if a > 0 else 0.0
        clean, usable, self_contra = [], [], {}
        best_dom, best_mis = 0.0, np.nan
        for m in mps:
            share = m["nprot"] / c["nprot"]
            dom, mis, h = m.get("loc_dom", 0.0), m.get("loc_mis", np.inf), m.get("loc_h", 0.0)
            la = m["local_amp"]
            pinned = m.get("pinned", False)
            ok_amp = np.isfinite(la) and la > 0
            # the amount this multiplet's own local fit implies must not predict the molecule's other multiplets
            # clearly above the noise where the data show nothing; if it does, the signal belongs to something else
            gone = _contradicted_by(m, mps, la, xg, yv, base, fw0, sf, sig, dx, S=S) \
                if ok_amp and h / sig >= S["LOQ_SNR"] else []
            if gone:
                self_contra[id(m)] = gone
                ok_amp = False
                centre = float((m["p"] * m["it"]).sum() / m["it"].sum() + m["d"])
                partner = _main_partner(res, k, centre, S["LINE_WINDOW_PPM"])
                if partner is not None and res["amp"][k] > 0:
                    notes.setdefault(lib[partner]["name"], []).append((c["name"],
                        f"{c['name']} (graded without a value) is fitted on its peak at {centre:.3f} ppm although "
                        f"its other multiplet(s) {','.join(map(str, gone))} are much smaller than that would need: "
                        f"this value may be too low by {c['name']}'s share, or {c['name']} may be wrongly assigned"))
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
                implies_missing_multiplets=",".join(map(str, gone)),
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
        why_not = ""
        if np.isfinite(ub) and ub < lod:
            tier, conc = "Not detected", 0.0
        elif contradicted:
            tier, conc = "Upper bound only", np.nan
            method = f"multiplet(s) {','.join(map(str, contradicted))} predicted but absent -> data upper bound"
        elif np.isfinite(ub) and conc > S["MAX_CONC_OVER_BOUND"] * ub:
            tier, conc = "Upper bound only", np.nan
            method = "fitted value exceeds what the data allow at another multiplet -> data upper bound"
        elif clean and (not np.isfinite(cross) or cross <= S["MAX_REPORTER_RATIO"]) and \
                (len(clean) >= 2 or clean[0]["nprot"] / c["nprot"] >= S["SINGLE_REPORTER_MIN_SHARE"]) and \
                not (why_not := _reporter_doubt(clean, mps, amp_use, self_contra, c, sig, cf, S=S)):
            tier = "Quantified"
        elif use:
            tier = "Overlapped (semi-quantitative)"
            if clean and why_not:
                method += f"; not Quantified: {why_not}"
        else:
            wc = a * cf if a > 0 else 0.0
            wc_se = res["se"][k] * cf if np.isfinite(res["se"][k]) else np.nan
            n_fit = sum(1 for m in mps if m.get("loc_h", 0) > S["LOQ_SNR"] * sig and id(m) not in self_contra
                        and m.get("loc_mis", np.inf) <= S["DECONV_MAX_MISFIT"]
                        and m.get("loc_dom", 0.0) >= S["DECONV_MIN_DOMINANCE"])
            if (np.isfinite(ub) and np.isfinite(wc_se) and wc > lod and wc <= ub
                    and n_fit >= S["DECONV_MIN_MULTIPLETS"]
                    and wc_se <= S["DECONV_MAX_REL_SE"] * wc):
                tier, conc, se = DECONV, wc, wc_se
                method = (f"whole-compound deconvolution over {n_fit} fitted multiplets, "
                          f"no resolved reporter (data upper bound {ub:.3g} mM)")
            else:
                if not np.isfinite(ub):
                    ub = conc                 # no multiplet gave a bound: the whole-compound fit is the bound
                tier, conc = "Upper bound only", np.nan
                method = "no usable reporter -> data upper bound"
        if tier == "Upper bound only":
            se = np.nan                       # a bound is not a measurement (the bound is in upper_bound_mM)
        left_out = next((m["dropped"] for m in mps if m.get("dropped")), "")
        if left_out:
            method = f"{left_out}; {method}"
        # centre (ppm, as fitted) of the multiplet(s) the value comes from, for the tiers that use reporters
        rep_ppm = ";".join(f"{float((m['p'] * m['it']).sum() / m['it'].sum() + m['d']):.4f}" for m in use) \
            if tier in ("Quantified", "Overlapped (semi-quantitative)") else ""
        rows.append(dict(metabolite=c["name"], cls=c["cls"], tier=tier, conc_mM=conc, se_mM=se,
                         whole_compound_fit_mM=a * cf if a > 0 else 0.0, method=method, reporter_ppm=rep_ppm,
                         LOD_mM=lod, upper_bound_mM=ub, snr=round(snr, 1), best_dominance=round(best_dom, 3),
                         best_misfit=round(best_mis, 3) if np.isfinite(best_mis) else np.nan,
                         n_reporters=len(use), reporter_ratio=round(cross, 3) if np.isfinite(cross) else np.nan,
                         source=c["source"], provenance=c["provenance"]))
    # a note on a compound only when the one fitted on its peak ends up without a value (likely not really there)
    no_value = {r["metabolite"] for r in rows if r["tier"] in ("Upper bound only", "Not detected")}
    for r in rows:
        r["check_note"] = "; ".join(t for who, t in notes.get(r["metabolite"], []) if who in no_value) \
            if r["tier"] != NOTMEAS else ""
    return pd.DataFrame(rows), pd.DataFrame(prow)


def _main_partner(res, k, ppm, half):
    """The compound (index) other than k with the largest fitted signal at ppm, or None."""
    best, who = 0.0, None
    xp = np.array([ppm])
    for o in res["MP"]:
        if o["k"] == k or o.get("global_amp", 0.0) <= 0:
            continue
        v = o["global_amp"] * float(mcol(xp, o["p"], o["it"], o["d"], res["fw0"] * o["wf"], res.get("eta", 0.7), half=half)[0])
        if v > best:
            best, who = v, o["k"]
    return who


def _reporter_doubt(clean, mps, amp_use, self_contra, c, sig, cf, *, S):
    """Why reporter value(s) that pass the per-multiplet tests should still not be graded Quantified (or ""):
    the molecule's other multiplets (above LOQ, outside ignored regions, any dominance) give a different value,
    or a lone singlet reporter is not confirmed by another multiplet."""
    def measurable(m):
        return (m["unit_peak"] > 0 and _excluded_fraction(m, S=S) < 0.5
                and (m["nprot"] >= 0.3 or m["nprot"] / c["nprot"] >= 0.1))
    rest = [m for m in mps if all(m is not u for u in clean) and measurable(m)]
    others = [m for m in rest if np.isfinite(m["local_amp"]) and m["local_amp"] > 0 and id(m) not in self_contra
              and m.get("loc_h", 0.0) >= S["LOQ_SNR"] * sig]
    lim = S["MAX_REPORTER_DISAGREE"]
    if others and amp_use > 0:
        med = float(np.median([m["local_amp"] for m in others]))
        rel = med / amp_use - 1
        if abs(rel) > lim:
            return (f"the other multiplet(s) {','.join(str(m['j'] + 1) for m in others)} give {med * cf:.3g} mM, "
                    f"{rel:+.0%} from the reporter (more than {lim:.0%})")
    if S["SINGLET_NEEDS_SECOND"] and len(clean) == 1 and _n_lines(clean[0]) == 1 and rest:
        agree = [m for m in others if abs(m["local_amp"] / amp_use - 1) <= lim] if amp_use > 0 else []
        if not agree:
            return "the only clean reporter is a singlet and no other multiplet of the molecule confirms it"
    return ""


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


def shown_value(tier, conc, upper_bound):
    """The number to show for one result: the upper bound for 'Upper bound only' (to be shown as <=value),
    otherwise the measured value. (Results from before 1.8.2 kept the bound in conc_mM.)"""
    if tier == "Upper bound only":
        return upper_bound if np.isfinite(upper_bound) or not np.isfinite(conc) else conc
    return conc


# one colour table for every tier, used by the PDF profiles and by the window's tables:
# tier: (marker colour in plots, table cell background, table cell text); cell colours keep text contrast >= 4.5:1
TIER_COLOURS = {"Quantified": ("#1d4ed8", "#1d4ed8", "#ffffff"),
                "Overlapped (semi-quantitative)": ("#b45309", "#b45309", "#ffffff"),
                DECONV: ("#6b7280", "#9ca3af", "#111827"),
                MANUAL: ("#0f766e", "#0f766e", "#ffffff"),
                "Upper bound only": ("#8a8a8a", "#e5e7eb", "#374151"),
                "Not detected": ("#8a8a8a", "#ffffff", "#4b5563"),
                NOTMEAS: ("#8a8a8a", "#f3f4f6", "#4b5563")}

TIER_ORDER = {"Quantified": 0, "Overlapped (semi-quantitative)": 1, DECONV: 2, MANUAL: 3,
              "Upper bound only": 4, "Not detected": 5, NOTMEAS: 6}


def apply_manual(out, gr):
    """Grades of one sample with its manual adjustments applied: an adjusted metabolite gets the tier
    'Manually adjusted' and the value set by hand; the automatic tier and value are kept alongside."""
    if not len(gr):
        return gr
    smp = str(gr["sample"].iloc[0])
    ed = read_manual(out, smp)
    gr = gr.copy()
    for c in ("auto_tier", "auto_conc_mM", "manual_note"):
        if c not in gr:
            gr[c] = np.nan if c == "auto_conc_mM" else ""
    gr["auto_tier"] = gr["auto_tier"].astype(object)
    if "reporter_ppm" in gr:
        gr["reporter_ppm"] = gr["reporter_ppm"].astype(object)
    gr["manual_note"] = gr["manual_note"].astype(object)
    for i in gr.index:
        e = ed.get(gr.at[i, "metabolite"])
        if e is None or gr.at[i, "tier"] == NOTMEAS:
            continue
        gr.at[i, "auto_tier"], gr.at[i, "auto_conc_mM"] = gr.at[i, "tier"], gr.at[i, "conc_mM"]
        gr.at[i, "tier"], gr.at[i, "conc_mM"], gr.at[i, "se_mM"] = MANUAL, float(e["conc_mM"]), np.nan
        if "reporter_ppm" in gr:
            gr.at[i, "reporter_ppm"] = ""        # a manual value is the whole compound, not one peak
        how = e.get("reason") or "set by hand"
        gr.at[i, "method"] = f"{how} in the Review step (shift {e.get('shift_Hz', 0):+.2f} Hz, " \
                             f"linewidth x{e.get('width', 1.0):.2f})"
        gr.at[i, "manual_note"] = f"{e.get('note', '')} [{e.get('user', '')}, {e.get('time', '')}]".strip()
    return gr
