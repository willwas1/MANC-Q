"""Command line:

    python -m mancq          (no arguments) opens the window
    python -m mancq run  <spectra folder> <output folder> --ref-mm 0.5
    python -m mancq demo [<output folder>]
    python -m mancq regrade <output folder>
    python -m mancq gui
    python -m mancq shortcut      (Windows: put a MANC-Q shortcut on the desktop)
"""
import argparse
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(prog="mancq", description=(
        "1H NMR metabolite quantification by whole-spectrum fitting of simulated GISSMO spin systems. "
        "Input: Bruker processed 1D 1H spectra with a TSP or DSS internal standard."))
    sub = ap.add_subparsers(dest="cmd")

    r = sub.add_parser("run", help="quantify a folder of spectra (or one EXPNO folder)")
    r.add_argument("spectra", help="folder containing EXPNO folders (each with pdata/<procno>/1r), or one EXPNO folder")
    r.add_argument("output", help="output folder (created if needed)")
    r.add_argument("--ref-mm", type=float, required=True,
                   help="concentration of the TSP/DSS reference in the NMR tube, mM (e.g. 0.5)")
    r.add_argument("--reference", default="TSP", choices=["TSP", "DSS"], help="reference compound name (default TSP)")
    r.add_argument("--dilution", type=float, default=1.0,
                   help="sample-to-tube dilution factor applied to every spectrum (default 1 = report tube mM)")
    r.add_argument("--dilution-file", default=None,
                   help="optional CSV with columns sample,dilution to set the factor per EXPNO")
    r.add_argument("--procno", type=int, default=1, help="processing number to read (default 1)")
    r.add_argument("--samples", nargs="*", default=None, help="only these EXPNO folders")
    r.add_argument("--workers", type=int, default=None, help="parallel processes (default: CPU cores - 1)")
    r.add_argument("--pluronic", default="auto", choices=["auto", "yes", "ignore", "no"],
                   help="Pluronic F-68 / poloxamer (cell-culture media): auto = model it if present (default), "
                        "yes = always model it, ignore = leave its regions out of the fit, no = not in the samples")
    r.add_argument("--exclude", nargs=2, type=float, action="append", metavar=("FROM", "TO"),
                   help="extra ppm region to leave out of the fit (repeatable); water 4.60-5.00 is always excluded")
    r.add_argument("--fids", default="missing", choices=["missing", "all", "never"],
                   help="raw FIDs: process them where no TopSpin-processed spectrum exists (default), always, or never")
    r.add_argument("--lb", type=float, default=0.3, help="line broadening in Hz when processing FIDs (default 0.3)")
    r.add_argument("--no-overlays", action="store_true", help="skip the per-spectrum overlay PDFs (faster)")
    r.add_argument("--no-resume", action="store_true", help="refit spectra that already have results in the output folder")

    d = sub.add_parser("demo", help="make synthetic spectra (one processed, one raw FID) of known composition and quantify them")
    d.add_argument("output", nargs="?", default="mancq_demo", help="demo folder (default ./mancq_demo)")
    d.add_argument("--field", type=float, default=600.13, help="spectrometer frequency to simulate, MHz")

    sub.add_parser("gui", help="open the window (the same as running mancq with no arguments)")

    sub.add_parser("shortcut", help="Windows: put a MANC-Q shortcut on the desktop that uses this Python")

    g = sub.add_parser("regrade", help="re-apply grading rules to a finished run without refitting")
    g.add_argument("output", help="output folder of an earlier run")

    a = ap.parse_args(argv)
    if a.cmd is None:
        a.cmd = "gui"
    if a.cmd == "run":
        import os
        if not os.path.isdir(a.spectra):
            print(f"Cannot find the spectra folder: {a.spectra}\n"
                  "Give the TopSpin dataset folder (the one containing the numbered EXPNO folders). "
                  "Put the path in double quotes if it contains spaces.")
            return 2
        from .engine import run
        extra = dict(REFERENCE=a.reference, PROCNO=a.procno, SAMPLES=a.samples, N_WORKERS=a.workers,
                     PLURONIC=a.pluronic, OVERLAYS=not a.no_overlays, RESUME=not a.no_resume,
                     FID_MODE=a.fids, FID_LB=a.lb)
        if a.exclude:
            extra["EXCLUDE"] = [(4.60, 5.00)] + [tuple(sorted(e)) for e in a.exclude]
        run(a.spectra, a.output, a.ref_mm, dilution=a.dilution, dilution_file=a.dilution_file, **extra)
    elif a.cmd == "demo":
        from .demo import demo
        demo(a.output, a.field)
    elif a.cmd == "gui":
        try:
            import PySide6  # noqa: F401
        except ImportError:
            print("The window needs PySide6, which is not installed.\n"
                  "  On Windows, the easiest way is to double-click MANC-Q.bat in the MANC-Q folder.\n"
                  "  Or install it with:  python -m pip install PySide6\n"
                  "The command line works without it; see:  python -m mancq --help")
            return 1
        try:
            from .gui.app import main as gui_main
            return gui_main()
        except Exception:
            if sys.stdout is None and sys.platform == "win32":   # started from a shortcut: no console to show errors
                import ctypes
                import traceback
                ctypes.windll.user32.MessageBoxW(None, "MANC-Q could not start:\n\n" + traceback.format_exc()[-1500:],
                                                 "MANC-Q", 0x10)
            raise
    elif a.cmd == "shortcut":
        from .shortcut import make_shortcut
        return make_shortcut()
    elif a.cmd == "regrade":
        from .regrade import regrade
        regrade(a.output)


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()
    sys.exit(main())
