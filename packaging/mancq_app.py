"""Entry point for the packaged Windows application (PyInstaller).

    MANC-Q.exe              opens the window
    MANC-Q.exe --selftest   runs the synthetic demo without a window and writes selftest.txt next to the
                            results (checks a fresh install, including parallel fitting inside the .exe)
"""
import multiprocessing
import os
import sys
import tempfile

if __name__ == "__main__":
    multiprocessing.freeze_support()          # needed for parallel fitting inside a packaged .exe
    if "--selftest" in sys.argv:
        from mancq.demo import demo
        out = os.path.join(tempfile.gettempdir(), "MANC-Q selftest")
        chk = demo(out)
        good = chk[chk.tier.isin(["Quantified", "Overlapped (semi-quantitative)"])]
        msg = (f"MANC-Q self-test: median recovered/true = {good.ratio.median():.3f} over {len(good)} values "
               f"(processed spectrum and raw FID)")
        open(os.path.join(out, "selftest.txt"), "w").write(msg + "\n")
        print(msg)
        sys.exit(0 if abs(good.ratio.median() - 1) < 0.05 else 1)
    from mancq.gui.app import main
    sys.exit(main())
