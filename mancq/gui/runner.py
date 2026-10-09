"""mancq.gui.runner - runs a GUI job in a separate process (no Qt imports here)."""
import traceback


def run_job(job, events, stop_flag):
    """Runs in a separate process: the engine, including any processing chosen in the Process spectra step
    (passed as the PROCESSING setting, so the window and the command line process spectra with the same code).
    Progress goes back through `events`. A separate process keeps the window responsive (no shared interpreter
    lock) and survives engine crashes."""
    try:
        from mancq import engine as q_
        q_.run(job["folder"], job["out"], job["ref_mm"], dilution=job["dilution"], dilution_file=job["dilution_file"],
               progress=lambda e, s, info: events.put((e, s or "", info or {})),
               cancel=stop_flag.is_set, **job["settings"])
    except Exception as exc:
        events.put(("__error__", "", dict(msg=f"{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}")))
    events.put(("__end__", "", {}))
