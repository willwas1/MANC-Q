"""Re-apply the grading rules to a finished run without refitting (seconds instead of minutes per spectrum).
Uses the fit state cached in <output>/_per_sample/*_fit.npz. Change grading values in SETTINGS first."""
import json
import os
import time
import numpy as np
import pandas as pd
from . import engine as q


GRADING_KEYS = ("LOD_SNR", "LOQ_SNR", "Q_MIN_DOMINANCE", "Q_MAX_MISFIT", "MAX_REPORTER_RATIO", "MAX_CONC_OVER_BOUND",
                "SINGLE_REPORTER_MIN_SHARE", "CONTRADICT_MIN_SNR", "CONTRADICT_FRACTION", "O_MIN_DOMINANCE",
                "O_MAX_MISFIT", "DECONV_MAX_REL_SE", "DECONV_MIN_MULTIPLETS")


def restore_settings(out, grading_from_run=False):
    """Put a finished run's own settings (ignored regions, reference, ...) back into the engine settings.
    Grading rules come from the current code unless grading_from_run is True."""
    sfile = os.path.join(out, "_per_sample", "run_settings.json")
    if not os.path.exists(sfile):
        return
    saved = json.load(open(sfile))
    for k in ("EXCLUDE", "FIT_RANGE", "PLURONIC_REGIONS"):
        if k in saved:
            saved[k] = [tuple(x) for x in saved[k]] if k != "FIT_RANGE" else tuple(saved[k])
    saved.pop("SAMPLE_DIRS", None)
    saved.pop("MANCQ_VERSION", None)
    keep = {} if grading_from_run else {k: q.S[k] for k in GRADING_KEYS}
    q.S.update(saved)
    q.S.update(keep)


def regrade(out, library_dir=None):
    per = os.path.join(out, "_per_sample")
    restore_settings(out)
    files = sorted([f for f in os.listdir(per) if f.endswith("_fit.npz")], key=lambda f: (len(f), f))
    samples = [f[:-len("_fit.npz")] for f in files]
    libdir = library_dir or os.path.join(q.HERE, "library")
    logf = open(os.path.join(out, "regrade_log.txt"), "w")

    def log(msg):
        print(msg, flush=True); logf.write(str(msg) + "\n"); logf.flush()

    log(f"regrade  {time.strftime('%Y-%m-%d %H:%M')}   {len(samples)} samples in {out}")
    grades, peaks, obs, qc, lib = [], [], [], [], None
    for smp in samples:
        d = np.load(os.path.join(per, f"{smp}_fit.npz"), allow_pickle=False)
        if "mp_json" in d:
            mp_list = json.loads(str(d["mp_json"]))
        else:
            # results written by MANC-Q 1.1.0 or earlier store this part with Python pickle, which can run code
            # when loaded; only accept those files from your own runs
            log(f"  {smp}: older results format; loading it as you created it yourself")
            mp_list = list(np.load(os.path.join(per, f"{smp}_fit.npz"), allow_pickle=True)["mp"])
        sf = float(d["sf"])
        pl = bool(d["pluronic"]) if "pluronic" in d else True
        slib = q.prepare_library(q.build_library(sf, libdir, os.path.join(out, "_cache")), sf, pl)
        lib = lib or q.prepare_library(q.build_library(sf, libdir, os.path.join(out, "_cache")), sf, True)
        MP = []
        for m in mp_list:
            m = dict(m)
            m["p"], m["it"] = np.asarray(m["p"], float), np.asarray(m["it"], float)
            MP.append(m)

        class Grid:
            x = np.asarray(d["x"], float)
            dx = float(d["dx"])
        res = dict(grid=Grid(), yv=np.asarray(d["yv"], float), fit=np.asarray(d["fit"], float),
                   base_fit=np.asarray(d["base"], float), sigma=float(d["sigma"]), fw0=float(d["fw0"]),
                   amp=d["amp"], se=d["se"], MP=MP)
        gr, pk = q.grade(res, slib, sf, float(d["cf"]))
        gr.insert(0, "sample", smp); pk.insert(0, "sample", smp)
        gr.to_csv(os.path.join(per, f"{smp}_grades.csv"), index=False)
        pk.to_csv(os.path.join(per, f"{smp}_peaks.csv"), index=False)
        grades.append(gr); peaks.append(pk)
        obs.append(pd.read_csv(os.path.join(per, f"{smp}_obs.csv")))
        qc.append(pd.read_csv(os.path.join(per, f"{smp}_qc.csv"), dtype={"sample": str}).iloc[0].to_dict())
        log(f"  {smp}: {(gr.tier == 'Quantified').sum()} quantified")
    q.aggregate_outputs(out, grades, peaks, obs, qc, lib, samples, log)
    logf.close()
