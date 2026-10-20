"""
mancq.engine - 1H NMR metabolite quantification by whole-spectrum library fitting.

For each processed 1D 1H spectrum (Bruker pdata/<procno>/1r) the engine:
  1. fits the reference singlet (TSP or DSS, set to 0.000 ppm) for chemical-shift referencing,
     linewidth and the area used for absolute quantification;
  2. builds the compound library at the spectrometer frequency of the data: GISSMO spin systems
     (chemical shifts and J couplings, BMRB) are simulated exactly, so second-order multiplets are
     correct at any field; a few compounds use fixed peak lists instead;
  3. fits the whole spectrum as  sum(amount x compound spectrum) + smooth baseline + free
     "unassigned" singlets, with small per-compound and per-multiplet shift and linewidth freedom;
  4. converts amounts to mM with the reference area and concentration, applies dilution;
  5. grades every compound (Quantified / Overlapped / Deconvolution estimate / Upper bound only /
     Not detected) from objective measures, and writes tables and overlay plots.

The parts: settings.py (defaults, Settings), io.py (reading spectra and files), library.py
(compound library), fit.py (fitting one spectrum), grade.py (tiers), report.py (tables, workbook,
PDFs). This module runs them for a folder of spectra (run) and re-exports all of them.

Command line: see `python -m mancq --help`.
"""
import json
import os
import time
import numpy as np
import pandas as pd
# everything is importable from mancq.engine, as before the split into modules
from .settings import (SETTINGS, HERE, _read_version, __version__, DEFAULT_SETTINGS, TIMING, RESULTS_FORMAT,
                       Settings, DECONV, NOTMEAS, MANUAL, _json_default)
from .io import (read_spectrum, _usable, list_samples, _num, _flag, read_processing, read_acquisition,
                 read_dilutions, EXISTING_MESSAGE, _RESULT_FILES, _RESULT_DIRS, ResultsExist, has_results,
                 previous_run, new_run_folder, clear_results, sample_source, _raw_folder, manual_file, read_manual)
from .library import parse_gissmo, build_library, prepare_library, group_multiplets, split_multiplets, library_hash
from .fit import (pv, noise_sigma, fit_tsp, PH_COLUMNS, ph_from_shift, shift_from_ph, _peak_candidates, estimate_ph,
                  histidine_targets, mcol, mcol_d, spline_basis, Grid, fit_spectrum, _observed_height,
                  _contradicted_by, _n_lines, _excluded_fraction, detect_pluronic)
from .grade import (grade, _main_partner, _reporter_doubt, observed_peaks, shown_value, TIER_COLOURS, TIER_ORDER,
                    apply_manual)
from .report import (_plt, plt, PALETTE, TAG, gapped, overlay_pdf, gallery_from_csv, profiles_pdf, _metabolite_fwhm,
                     qc_flags, NOT_DETECTED_NOTE, UPPER_BOUND_NOTE, ph_readme, aggregate_outputs)


# =============================================================================
# 8. Main
# =============================================================================
def process_sample(job):
    """Fit, grade and plot one spectrum (runs in a worker process). Results are written
    to <out>/_per_sample as they are produced, so a crash or a stop loses nothing: the
    next run re-uses them and only fits what is missing."""
    smp, d, libdir, out, settings = job
    S = Settings(settings)
    per = os.path.join(out, "_per_sample")
    done = [os.path.join(per, f"{smp}_{k}.csv") for k in ("grades", "peaks", "obs", "qc")]
    note = []
    if S["RESUME"] and all(os.path.exists(f) for f in done):
        why = _reuse_problem(per, smp, d, S=S)
        if why is None:
            gr, pk, ob, qdf = (pd.read_csv(f) for f in done)
            return gr, pk, ob, qdf.iloc[0].to_dict(), None, [f"Sample {smp}: re-used earlier result (same run, "
                                                             f"same settings)"]
        note = [f"Sample {smp}: its earlier result cannot be used ({why}); fitting it again"]
    try:
        gr, pk, ob, q, res, lines = _fit_one(smp, d, libdir, out, S=S)
        return gr, pk, ob, q, res, note + lines
    except Exception as exc:
        import traceback
        return None, None, None, dict(sample=smp, error=f"{type(exc).__name__}: {exc}"), None, \
            note + [f"Sample {smp} FAILED: {type(exc).__name__}: {exc}", traceback.format_exc()]


# settings that do not change a sample's fitted values (where files go, speed, plots); every other setting is
# written to <sample>_source.json and must match for a saved result to be re-used
_NOT_FIT_SETTINGS = {"DATA_DIR", "OUT_DIR", "DILUTION_FILE", "DEFAULT_DILUTION", "SAMPLES", "N_WORKERS", "RESUME",
                     "GALLERY_MAX_SAMPLES", "LIBRARY_DIR", "OVERLAYS", "SAMPLE_DIRS", "PROCESSING", "FID_MODE",
                     "FID_LB", "OVERLAY_REGIONS", "SHIFT_TABLE_TIERS", "SHIFT_TABLE_MIN_FRACTION"}


def _source_desc(smp, dilution, *, S):
    """Everything a sample's saved result depends on: the spectrum and its processing, the run it belongs to,
    the MANC-Q version, the library, the reference concentration, its dilution and every fitting and grading
    setting. A saved result is only re-used if all of these are unchanged."""
    d, pn = sample_source(smp, S=S)
    u = (S.get("PROCESSING_USED") or {}).get(smp)
    d = os.path.normcase(os.path.abspath(d))
    out = os.path.normcase(os.path.abspath(S["OUT_DIR"])) if S.get("OUT_DIR") else None
    if out and d.startswith(out + os.sep):        # processed by MANC-Q into the results folder: the folder may move
        d = "<results>" + d[len(out):].replace(os.sep, "/")
    desc = dict(spectrum=d, procno=pn, processing=u["requested"] if u else None, run_id=S.get("RUN_ID"),
                version=__version__, library=S.get("LIBRARY_HASH"), dilution=float(dilution),
                settings={k: S.get(k) for k in DEFAULT_SETTINGS if k not in _NOT_FIT_SETTINGS})
    return json.loads(json.dumps(desc, sort_keys=True, default=str))     # tuples become lists, as when read back


_SOURCE_WHY = dict(spectrum="a different spectrum", procno="a different procno", processing="different processing",
                   run_id="it belongs to an earlier run", version="made by another MANC-Q version",
                   library="the compound library has changed", dilution="the dilution has changed")


def _reuse_problem(per, smp, dilution, *, S):
    """None if the saved result of `smp` can be re-used, otherwise the reason (in words) why not."""
    f = os.path.join(per, f"{smp}_source.json")
    if not os.path.exists(f):
        return "no record of how it was made"
    try:
        with open(f, encoding="utf-8") as fh:
            old = json.load(fh)
    except ValueError:
        old = None
    if not isinstance(old, dict) or "run_id" not in old:
        return "made by MANC-Q 1.8.0 or earlier"
    new = _source_desc(smp, dilution, S=S)
    why = [_SOURCE_WHY.get(k, k) for k in new if k != "settings" and old.get(k) != new[k]]
    so, sn = old.get("settings") or {}, new["settings"]
    changed = sorted(k for k in set(so) | set(sn) if so.get(k) != sn.get(k))
    if changed:
        why.append("changed setting" + ("s " if len(changed) > 1 else " ") + ", ".join(changed))
    return "; ".join(why) or None


def _fit_one(smp, d, libdir, out, *, S):
    lines = []
    log = lambda msg: lines.append(msg)
    per = os.path.join(out, "_per_sample")
    os.makedirs(per, exist_ok=True)
    mf = os.path.join(per, f"{smp}_manual.json")
    if os.path.exists(mf):
        aside = os.path.join(per, f"{smp}_manual_set_aside_{time.strftime('%Y%m%d_%H%M%S')}.json")
        os.replace(mf, aside)
        lines.append(f"   NOTE: sample {smp} was fitted again, so the manual adjustments made on its earlier fit "
                     f"no longer apply; they were moved to _per_sample/{os.path.basename(aside)}")
    expdir, procno = sample_source(smp, S=S)
    ppm, y, sf, title = read_spectrum(expdir, procno, S=S)
    pluronic = detect_pluronic(ppm, y, S=S) if S["PLURONIC"] == "auto" else S["PLURONIC"] == "yes"
    lib = prepare_library(build_library(sf, libdir, os.path.join(out, "_cache"), S=S), sf, pluronic, S=S)
    log(f"Sample {smp}  (title: {title}, dilution x{d}, {sf:.2f} MHz"
        f"{', Pluronic terms added' if pluronic else ''})")
    t0 = time.time()
    res = fit_spectrum(ppm, y, sf, lib, log, S=S)
    cf = S["TSP_MM_IN_TUBE"] * S["TSP_PROTONS"] / res["tsp"]["area"] * d
    gr, pk = grade(res, lib, sf, cf, S=S)
    gr.insert(0, "sample", smp); pk.insert(0, "sample", smp)
    ob = observed_peaks(res, lib, smp, sf)
    sig_area = max((res["yv"] - res["base_fit"]).clip(min=0).sum(), 1e-9)
    q = dict(sample=smp, title=title, SF_MHz=sf, dilution=d, TSP_area=res["tsp"]["area"],
             TSP_fwhm_Hz=res["tsp"]["fwhm_hz"], TSP_eta=res["tsp"]["eta"], TSP_fit_rel_rms=res["tsp"]["rel_rms"],
             noise_sigma=res["sigma"], TSP_SNR=res["tsp"]["height"] / res["sigma"],
             GISSMO_offset_Hz=res["g_off"] * sf, lineshape_eta=res["eta"], chi2_per_point=res["chi2n"],
             unassigned_singlets=len(res["unknown"]),
             signal_explained_by_library=res["comp_fit"].sum() / sig_area,
             signal_in_unassigned_singlets=res["unk_fit"].sum() / sig_area if res["unk_fit"].size else 0.0,
             abs_residual_fraction=np.abs(res["yv"] - res["fit"]).sum() / sig_area,
             fit_time_s=round(time.time() - t0, 1), **res["ph"],
             TSP_phase_deg=round(res["tsp"].get("phase_deg", np.nan), 2), noise_method=res["noise_method"],
             metabolite_fwhm_Hz=_metabolite_fwhm(res, sf, S=S),
             **read_acquisition(_raw_folder(smp, expdir, S=S)),
             processing=((S.get("PROCESSING_USED") or {}).get(smp) or {}).get(
                 "description", f"as processed in TopSpin (pdata/{procno})"))
    if S["OVERLAYS"]:
        overlay_pdf(res, lib, gr, smp, os.path.join(out, f"overlay_sample_{smp}.pdf"), S=S)
    x = res["grid"].x
    cols = {"ppm": x, "observed": res["yv"], "simulated_total": res["fit"], "baseline": res["base_fit"]}
    for k, c in enumerate(lib):
        if res["amp"][k] > 0:
            cols[c["name"]] = res["comp_fit"][:, k]
    if res["unk_fit"].size:
        cols["unassigned_singlets"] = res["unk_fit"].sum(1)
    pd.DataFrame(cols).to_csv(os.path.join(out, f"simulated_spectrum_sample_{smp}.csv.gz"), index=False, float_format="%.6g")
    np.savez_compressed(os.path.join(per, f"{smp}_fit.npz"), cf=cf, sigma=res["sigma"], fw0=res["fw0"],
                        x=res["grid"].x.astype(np.float32), yv=res["yv"].astype(np.float32),
                        fit=res["fit"].astype(np.float32), base=res["base_fit"].astype(np.float32),
                        dx=res["grid"].dx, sf=sf, amp=res["amp"], se=res["se"], pluronic=pluronic, eta=res["eta"],
                        mp_json=np.array(json.dumps([dict({k: (v.tolist() if isinstance(v, np.ndarray) else v)
                                                           for k, v in m.items()}, name=lib[m["k"]]["name"])
                                                      for m in res["MP"]], default=_json_default)))
    gr.to_csv(os.path.join(per, f"{smp}_grades.csv"), index=False)
    pk.to_csv(os.path.join(per, f"{smp}_peaks.csv"), index=False)
    ob.to_csv(os.path.join(per, f"{smp}_obs.csv"), index=False)
    pd.DataFrame([q]).to_csv(os.path.join(per, f"{smp}_qc.csv"), index=False)
    with open(os.path.join(per, f"{smp}_log.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    with open(os.path.join(per, f"{smp}_source.json"), "w", encoding="utf-8") as fh:
        json.dump(_source_desc(smp, d, S=S), fh, sort_keys=True, indent=1)
    return gr, pk, ob, q, None, lines


def run(data_dir, out_dir, ref_mm, dilution=1.0, dilution_file=None, progress=None, cancel=None, existing="stop",
        **settings):
    """Quantify every spectrum under data_dir (or the single EXPNO folder data_dir) into out_dir.

    existing: what to do if out_dir already holds MANC-Q results. "stop" (default) raises ResultsExist;
    "overwrite" deletes them first (see clear_results); "resume" carries on the stopped run in out_dir, re-using
    the spectra it finished if nothing that affects them has changed.
    progress(event, sample, info): optional callback; events are "run_started", "sample_started",
    "sample_done", "sample_failed", "cancelled", "aggregating", "finished".
    cancel(): optional function returning True to stop; spectra already running are finished first."""
    if ref_mm is None or ref_mm <= 0:
        raise ValueError("the reference (TSP/DSS) concentration in the tube, in mM, is required")
    say = progress or (lambda *a, **k: None)
    stop = cancel or (lambda: False)
    S = Settings(settings)                     # this run's own settings: the defaults with these changes
    S.update(DATA_DIR=os.path.abspath(data_dir) if data_dir else None, OUT_DIR=os.path.abspath(out_dir),
             TSP_MM_IN_TUBE=float(ref_mm), DEFAULT_DILUTION=float(dilution), DILUTION_FILE=dilution_file)
    if S["PLURONIC"] == "ignore":
        S["EXCLUDE"] = list(S["EXCLUDE"]) + [r for r in S["PLURONIC_REGIONS"] if r not in S["EXCLUDE"]]
    if S["SAMPLES"]:
        samples = list(S["SAMPLES"])
    else:
        samples = list_samples(S["DATA_DIR"], S["PROCNO"], fids=S["FID_MODE"] != "never")
        if samples == [""]:                       # data_dir is itself one EXPNO folder
            S["DATA_DIR"], one = os.path.split(S["DATA_DIR"].rstrip("/\\"))
            samples = [one]
    if not samples:
        raise FileNotFoundError(f"no spectra found: expected <EXPNO>/pdata/{S['PROCNO']}/1r under {data_dir}")
    out = S["OUT_DIR"]
    if existing not in ("stop", "overwrite", "resume"):
        raise ValueError(f"existing must be 'stop', 'overwrite' or 'resume', not {existing!r}")
    prev = previous_run(out) if has_results(out) else None
    if has_results(out) and existing == "stop":
        raise ResultsExist(f"{EXISTING_MESSAGE}\nResults folder: {out}")
    cleared = False
    if has_results(out) and existing == "overwrite":
        clear_results(out, keep=[dilution_file])
        cleared = True
    S["RUN_ID"] = None
    if existing == "resume" and prev and prev["run_id"]:
        if prev["samples"] and sorted(map(str, prev["samples"])) != sorted(map(str, samples)):
            raise ValueError(f"To carry on the stopped run in {out}, choose the same spectra as before "
                             f"(EXPNO {', '.join(map(str, prev['samples']))}).")
        S["RUN_ID"] = prev["run_id"]                 # same run: its finished spectra may be re-used
    resumed = S["RUN_ID"] is not None
    if not resumed:
        S["RUN_ID"] = f"{time.strftime('%Y%m%d-%H%M%S')}-{os.urandom(4).hex()}"
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "run_log.txt"), "w", encoding="utf-8") as logf:

        def log(msg):
            print(msg, flush=True); logf.write(msg + "\n"); logf.flush()

        if cleared:
            log("Overwrite: the results of the earlier run in this folder were deleted.")
        elif existing == "resume" and not resumed:
            log("Nothing to carry on in this folder (no stopped run of MANC-Q 1.8.1 or later): every spectrum is fitted.")
        libdir = S["LIBRARY_DIR"] or os.path.join(HERE, "library")
        S["LIBRARY_HASH"] = library_hash(libdir)
        dil = read_dilutions(S["DILUTION_FILE"], samples, log, S=S)
        # spectra MANC-Q processes itself, into <out>/_processed/<sample> (the raw data are never changed): raw FIDs
        # (automatic phasing unless PROCESSING says otherwise) and TopSpin spectra the user adjusted (PROCESSING)
        sd = dict(S.get("SAMPLE_DIRS") or {})
        for smp in samples:
            if smp not in sd and not os.path.isdir(os.path.join(S["DATA_DIR"], smp)):
                raise FileNotFoundError(f"EXPNO {smp} is not in {S['DATA_DIR']}")
        req = {str(k): dict(v) for k, v in (S.get("PROCESSING") or {}).items()}
        used = {}
        from . import fidproc
        todo = []
        for smp in samples:
            if smp in sd:
                continue
            ed = os.path.join(S["DATA_DIR"], smp)
            p = req.get(smp)
            if p is None and S["FID_MODE"] != "never" and fidproc.has_fid(ed) and \
                    (S["FID_MODE"] == "all" or not fidproc.has_processed(ed, S["PROCNO"])):
                p = dict(source="fid", lb=S["FID_LB"], ph0=None, ph1=None, baseline=False)
            if p is not None:
                todo.append((smp, ed, p))
        for k_, (smp, ed, p) in enumerate(todo, 1):
            if stop():
                log("Cancelled while processing spectra: nothing was fitted.")
                say("cancelled", None, dict(not_fitted=list(samples)))
                return None
            say("fid_processing", smp, dict(k=k_, n=len(todo), source=p.get("source", "fid")))
            p = dict(source=p.get("source") or "fid", lb=_num(p.get("lb")), ph0=_num(p.get("ph0")),
                     ph1=_num(p.get("ph1")), baseline=_flag(p.get("baseline")))
            r = fidproc.process(ed, S["PROCNO"], **p)
            dest = os.path.join(out, "_processed", smp)
            fidproc.write_processed(r, dest, title=f"EXPNO {smp} (processed by MANC-Q)")
            with open(os.path.join(dest, "processing.txt"), "w", encoding="utf-8") as fh:
                fh.write(
                f"Processed from {ed} by MANC-Q {__version__}: {r['description']}\n")
            sd[smp] = (dest, 1)
            used[smp] = dict(requested=p, source=p["source"], lb=round(r["lb"], 4), ph0=round(r["ph0"], 3),
                             ph1=round(r["ph1"], 3), baseline=r["baseline"], automatic_phase=bool(r.get("auto")),
                             description=r["description"])
            log(f"EXPNO {smp}: {r['description']}; check its overlay")
        S["SAMPLE_DIRS"] = sd
        S["PROCESSING_USED"] = used
        # the processing actually used, in the form --processing reads (blank phase = automatic)
        rows = [dict(sample=s_, source=u["source"], lb=u["lb"],
                     ph0="" if u["requested"]["ph0"] is None else u["ph0"],
                     ph1="" if u["requested"]["ph0"] is None else u["ph1"],
                     baseline=u["baseline"], description=u["description"]) for s_, u in used.items()]
        if rows:
            pd.DataFrame(rows).to_csv(os.path.join(out, "processing_used.csv"), index=False)
        os.makedirs(os.path.join(out, "_per_sample"), exist_ok=True)

        def save_run_settings(complete):
            S.save(os.path.join(out, "_per_sample", "run_settings.json"), RUN_SAMPLES=list(samples),
                   RUN_COMPLETE=complete)

        save_run_settings(False)
        log(f"MANC-Q {__version__}   {time.strftime('%Y-%m-%d %H:%M')}\nData: {S['DATA_DIR']}\n"
            f"Samples ({len(samples)}): {samples}\nReference: {S['REFERENCE']} {S['TSP_MM_IN_TUBE']} mM in tube\n"
            f"Ignored regions (ppm): {', '.join(f'{a:.2f}-{b:.2f}' for a, b in S['EXCLUDE'])}\n"
            f"Pluronic F-68: {S['PLURONIC']}")
        d0, p0 = sample_source(samples[0], S=S)
        _, _, sf, _ = read_spectrum(d0, p0, S=S)
        log(f"Library at {sf:.3f} MHz ...")
        lib = prepare_library(build_library(sf, libdir, os.path.join(out, "_cache"), S=S), sf, S["PLURONIC"] in ("auto", "yes"), S=S)
        log(f"  {sum(c['source']=='GISSMO' for c in lib)} GISSMO spin systems, {sum(c['source']=='SPIN' for c in lib)} other spin systems, "
            f"{sum(c['source']=='CASMDB' for c in lib)} peak lists, {sum(c['source']=='polymer' for c in lib)} polymer terms")
        nw = S["N_WORKERS"] or max(1, min(len(samples), (os.cpu_count() or 2) - 1))
        nw = max(1, min(nw, len(samples)))
        jobs = {smp: (smp, dil.get(smp, S["DEFAULT_DILUTION"]), libdir, out, dict(S)) for smp in samples}
        log(f"Fitting {len(samples)} spectra with {nw} worker(s); a few minutes per spectrum ...")
        say("run_started", None, dict(n=len(samples), workers=nw, sf=sf))
        results, cancelled = {}, []

        def finished(smp, output, t0):
            gr, pk, ob, q, res, lines = output
            for ln in lines:
                log(ln)
            results[smp] = output
            if gr is None:
                say("sample_failed", smp, dict(error=q.get("error", "unknown error"), seconds=time.time() - t0))
            else:
                say("sample_done", smp, dict(seconds=time.time() - t0, n_quantified=int((gr.tier == "Quantified").sum()),
                                             explained=float(q.get("signal_explained_by_library", np.nan)),
                                             reused=any("re-used" in ln for ln in lines)))

        queue = list(samples)
        # Spectra are fitted in worker processes with single-threaded linear algebra. Multi-threaded BLAS
        # gives tiny run-to-run rounding differences that can flip a grid-search choice; with one thread per worker the
        # results are identical every time, from the command line, the GUI or Python.
        for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                    "NUMEXPR_NUM_THREADS"):
            os.environ[var] = "1"
        from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
        import multiprocessing
        # "spawn" everywhere: the Windows behaviour, and safe when called from a GUI
        with ProcessPoolExecutor(max_workers=nw, mp_context=multiprocessing.get_context("spawn")) as ex:
            running = {}
            while queue or running:
                while queue and len(running) < nw and not stop():
                    smp = queue.pop(0)
                    say("sample_started", smp, {})
                    running[ex.submit(process_sample, jobs[smp])] = (smp, time.time())
                if stop() and queue:
                    cancelled, queue = queue, []
                if not running:
                    break
                done, _ = wait(list(running), timeout=1.0, return_when=FIRST_COMPLETED)
                for fut in done:
                    smp, t0 = running.pop(fut)
                    try:
                        output = fut.result()
                    except Exception as exc:          # worker crashed
                        output = (None, None, None, dict(sample=smp, error=f"{type(exc).__name__}: {exc}"), None,
                                  [f"Sample {smp} FAILED: {exc}"])
                    finished(smp, output, t0)
        if cancelled:
            log(f"\nCancelled: {len(cancelled)} spectra not fitted: {cancelled}")
            say("cancelled", None, dict(not_fitted=cancelled))
        grades, peaks, obs, qc, failed = [], [], [], [], []
        for smp in samples:
            if smp not in results:
                continue
            gr, pk, ob, q, res, lines = results[smp]
            if gr is None:
                failed.append(smp); qc.append(q); continue
            grades.append(gr); peaks.append(pk); obs.append(ob); qc.append(q)
        if failed:
            log(f"\n*** {len(failed)} spectra failed and are not in the tables: {failed}")
        if not grades:
            say("finished", None, dict(ok=False, failed=failed, cancelled=cancelled))
            raise RuntimeError("no spectrum could be fitted; see run_log.txt")
        say("aggregating", None, {})
        done_samples = [s_ for s_ in samples if s_ in results and results[s_][0] is not None]
        G = aggregate_outputs(out, grades, peaks, obs, qc, lib, done_samples, log, S=S)
        save_run_settings(not cancelled and not failed)
        say("finished", None, dict(ok=True, failed=failed, cancelled=cancelled, out=out))
        return G
