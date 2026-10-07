# MANC-Q: notes for Claude Code

MANC-Q (Metabolite Abundance by NMR Curve-fitting, Quantified) quantifies metabolites (mM) in 1D 1H NMR spectra
by fitting the whole spectrum with GISSMO spin systems simulated at the spectrometer field, using only the
TSP/DSS internal standard. Author: Will Smith, University of Manchester. Repository: github.com/willwas1/MANC-Q.
Will is a biologist, not a programmer: explain things plainly, give exact commands to paste, and never use
em dashes in anything you write for him.

## Rules (agreed with Will)

* Every change to the tool, however small, gets a new version number (1.1.2 -> 1.1.3 ...). Raise it in all of:
  `pyproject.toml`, `mancq/engine.py` (`__version__`), `CITATION.cff`, `packaging/installer.iss`
  (`MyAppVersion`), and the first line of `constraints.txt`. Add an entry at the top of `CHANGES.md`.
* Commit each version separately with a clear message, tag it `vX.Y.Z`, and only push or publish a release
  when Will says so. Publishing a GitHub release runs `.github/workflows/build-windows.yml`, which builds and
  attaches `MANC-Q-setup.exe`, a portable zip and `SHA256SUMS.txt`.
* Never add Will's own spectra, results, sample names, dilution factors, media details, or anything from
  industry collaborations. Demo data must stay synthetic (`mancq/demo.py`).
* Do not reintroduce pickle for saved files, or any network access in the tool itself.

## Layout

* `mancq/engine.py` fitting and grading (SETTINGS block at the top); `spinsim.py` spin simulation;
  `fidproc.py` raw FID processing; `regrade.py`; `results.py`; `demo.py` synthetic demo (one processed
  spectrum, one raw FID); `shortcut.py` Windows desktop shortcut; `__main__.py` command line
  (`mancq` with no arguments opens the window).
* `mancq/gui/app.py` PySide6 window (6 steps; Process FIDs step only shown when raw FIDs exist);
  `gui/runner.py` runs the fit in a separate process.
* `mancq/library/`: `compounds.csv` (149 compounds: 123 GISSMO, 25 CASMDB peak lists, 1 literature spin
  system), `gissmo_spin_systems.txt`, `other_spin_systems.txt`, `peaklists.csv`, `peaklists_manifest.csv`.
* `MANC-Q.bat` double-click launcher (downloads uv + private Python into `.runtime`); `constraints.txt` exact
  tested package versions; `packaging/` PyInstaller and Inno Setup build; `docs/` quick start, command line,
  screenshots (from the synthetic demo only).

## Checking a change

* `python -m mancq demo <folder>`: fits the synthetic spectra; median recovered/true should be about 1.00 for
  both EXPNO 10 (processed) and EXPNO 11 (raw FID). Takes a few minutes.
* `python -m mancq regrade <results folder>` should reproduce the same tiers and values.
* For window changes, open it (`python -m mancq`) and look at every step, including with no raw FIDs.

## Will's setup

* University Windows PC. Unsigned programs (the installer, MANC-Q.bat) are blocked by IT; he runs MANC-Q in
  Anaconda (`C:\Anaconda`). Install there with `python -m pip install -e . --no-deps --force-reinstall`
  (no `--no-deps` breaks on Anaconda's own PySide6/shiboken6). Desktop shortcut: `python -m mancq shortcut`.
* Windows dark mode: the window forces a light theme on purpose; keep text contrast at 4.5:1 or better.

## Known open items

* GitHub tidy-up: there are two old releases (v1.1 with files, v1.1.0 empty).
* Code signing (university IT or SignPath Foundation) so the installer runs on managed PCs.
* Automatic FID phasing can differ from TopSpin by 1 to 3 degrees, moving a few values 10 to 30 %.
