"""mancq.gui.runner - runs a GUI job in a separate process (no Qt imports here)."""
import os
import traceback


def run_job(job, events, stop_flag):
    """Runs in a separate process: FID processing then the engine. Progress goes back through `events`.
    A separate process keeps the window responsive (no shared interpreter lock) and survives engine crashes."""
    try:
        from mancq import engine as q_, fidproc as fp_
        sample_dirs = {}
        fids = [s for s in job["spectra"] if s["use_fid"]]
        for k, s in enumerate(fids, 1):
            if stop_flag.is_set():
                break
            events.put(("fid_processing", s["expno"], dict(k=k, n=len(fids))))
            res = fp_.process_fid(s["path"], lb=s["fid_lb"], ph0=s["fid_ph0"], ph1=s["fid_ph1"])
            dest = os.path.join(job["out"], "_fid_processed", s["expno"])
            fp_.write_processed(res, dest, title=f"{s['title'] or 'EXPNO ' + s['expno']} (from FID)")
            open(os.path.join(dest, "processing.txt"), "w").write(
                f"Processed from {s['path']} by MANC-Q: LB {res['lb']} Hz, SI {res['si']}, "
                f"ph0 {res['ph0']:.2f}, ph1 {res['ph1']:.2f} ({'automatic' if s['fid_ph0'] is None else 'set by the user'})\n")
            sample_dirs[s["expno"]] = (dest, 1)
        if stop_flag.is_set():
            events.put(("cancelled", "", dict(not_fitted=[s["expno"] for s in job["spectra"]])))
            events.put(("__end__", "", {}))
            return
        settings = dict(job["settings"])
        settings["SAMPLE_DIRS"] = sample_dirs
        q_.run(job["folder"], job["out"], job["ref_mm"], dilution=job["dilution"], dilution_file=job["dilution_file"],
               progress=lambda e, s, info: events.put((e, s or "", info or {})),
               cancel=stop_flag.is_set, **settings)
    except Exception as exc:
        events.put(("__error__", "", dict(msg=f"{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}")))
    events.put(("__end__", "", {}))


