"""Re-apply the grading rules to a finished run without refitting (seconds instead of minutes per spectrum).
Uses the fit state cached in <output>/_per_sample/*_fit.npz. Grading values come from SETTINGS in engine.py
(or regrade(out, LOQ_SNR=8, ...)); everything else from the run's own saved settings."""
import os
import time
import numpy as np
import pandas as pd
from . import engine as q
from . import results as rs


GRADING_KEYS = ("LOD_SNR", "LOQ_SNR", "Q_MIN_DOMINANCE", "Q_MAX_MISFIT", "MAX_REPORTER_RATIO", "MAX_CONC_OVER_BOUND",
                "SINGLE_REPORTER_MIN_SHARE", "CONTRADICT_MIN_SNR", "CONTRADICT_FRACTION", "O_MIN_DOMINANCE",
                "O_MAX_MISFIT", "DECONV_MAX_REL_SE", "DECONV_MIN_MULTIPLETS", "MAX_REPORTER_DISAGREE",
                "SINGLET_NEEDS_SECOND", "DECONV_MAX_MISFIT",
                "DECONV_MIN_DOMINANCE")


def run_settings(out, grading_from_run=False, **grading):
    """A finished run's own settings (ignored regions, reference, ...). The grading rules are the current ones
    (SETTINGS, changed by `grading`) unless grading_from_run is True."""
    S = q.Settings.load(out)
    if not grading_from_run:
        S.update({k: q.SETTINGS[k] for k in GRADING_KEYS})
    bad = [k for k in grading if k not in GRADING_KEYS]
    if bad:
        raise KeyError(f"not a grading setting: {', '.join(bad)} (regrade can only change {', '.join(GRADING_KEYS)})")
    S.update(grading)
    return S


def regrade(out, library_dir=None, **grading):
    per = os.path.join(out, "_per_sample")
    libdir = library_dir or os.path.join(q.HERE, "library")
    rs.check_library(out, libdir)               # same layout and compound library as this MANC-Q
    S = run_settings(out, **grading)
    files = sorted([f for f in os.listdir(per) if f.endswith("_fit.npz")], key=lambda f: (len(f), f))
    samples = [f[:-len("_fit.npz")] for f in files]
    libdir = library_dir or os.path.join(q.HERE, "library")
    with open(os.path.join(out, "regrade_log.txt"), "w", encoding="utf-8") as logf:

        def log(msg):
            print(msg, flush=True); logf.write(str(msg) + "\n"); logf.flush()

        log(f"regrade  {time.strftime('%Y-%m-%d %H:%M')}   {len(samples)} samples in {out}")
        grades, peaks, obs, qc, lib = [], [], [], [], None
        for smp in samples:
            with np.load(os.path.join(per, f"{smp}_fit.npz"), allow_pickle=False) as z:
                sf = float(z["sf"])
                pl = bool(z["pluronic"]) if "pluronic" in z else True
            slib = q.prepare_library(q.build_library(sf, libdir, os.path.join(out, "_cache"), S=S), sf, pl)
            lib = lib or q.prepare_library(q.build_library(sf, libdir, os.path.join(out, "_cache"), S=S), sf, True)
            d, MP = rs.load_fit(out, smp, slib, log)

            class Grid:
                x = np.asarray(d["x"], float)
                dx = float(d["dx"])
            res = dict(grid=Grid(), yv=np.asarray(d["yv"], float), fit=np.asarray(d["fit"], float),
                       base_fit=np.asarray(d["base"], float), sigma=float(d["sigma"]), fw0=float(d["fw0"]),
                       amp=d["amp"], se=d["se"], MP=MP, eta=float(d["eta"]) if "eta" in d else 0.7)
            gr, pk = q.grade(res, slib, sf, float(d["cf"]), S=S)
            gr.insert(0, "sample", smp); pk.insert(0, "sample", smp)
            gr.to_csv(os.path.join(per, f"{smp}_grades.csv"), index=False)
            pk.to_csv(os.path.join(per, f"{smp}_peaks.csv"), index=False)
            grades.append(gr); peaks.append(pk)
            obs.append(pd.read_csv(os.path.join(per, f"{smp}_obs.csv")))
            qc.append(pd.read_csv(os.path.join(per, f"{smp}_qc.csv"), dtype={"sample": str}).iloc[0].to_dict())
            log(f"  {smp}: {(gr.tier == 'Quantified').sum()} quantified")
        q.aggregate_outputs(out, grades, peaks, obs, qc, lib, samples, log, S=S)
