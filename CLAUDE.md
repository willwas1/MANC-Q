# MANC-Q: notes for Claude Code

MANC-Q (Metabolite Abundance by NMR Curve-fitting, Quantified) quantifies metabolites (mM) in 1D 1H NMR spectra
by fitting the whole spectrum with GISSMO spin systems simulated at the spectrometer field, using only the
TSP/DSS internal standard. Author: Will Smith, University of Manchester. Repository: github.com/willwas1/MANC-Q.
Will is a biologist, not a programmer: explain things plainly, give exact commands to paste, and never use
em dashes in anything you write for him.

## Rules (agreed with Will)

* Every change to the tool, however small, gets a new version number (1.1.2 -> 1.1.3 ...). The number is kept
  only in `pyproject.toml` (MANC-Q reads it from there, or from the installed package's metadata). Set it with
  `python packaging/set_version.py X.Y.Z`, which also copies it into `CITATION.cff` (with today's date),
  `packaging/installer.iss` and the first line of `constraints.txt`. Add an entry at the top of `CHANGES.md`.
* Settings: `SETTINGS` at the top of `engine.py` are the defaults and are never changed while MANC-Q runs. Each
  run makes its own `Settings` object and passes it to every function that needs it (as `S=`); it is saved in
  `_per_sample/run_settings.json` and read back with `Settings.load(out)`. Reading results folders, including
  everything needed for older ones, goes through `mancq/results.py` (`load_results`, `load_fit`); raise
  `RESULTS_FORMAT` in `engine.py` when the folder layout changes.
* Line endings: keep each file's existing ending when editing (`mancq/engine.py` is LF, `gui/app.py` CRLF).
* Commit each version separately with a clear message, tag it `vX.Y.Z`, and only push or publish a release
  when Will says so. Publishing a GitHub release runs `.github/workflows/build-windows.yml`, which builds and
  attaches `MANC-Q-setup.exe`, a portable zip and `SHA256SUMS.txt`.
* Never add Will's own spectra, results, sample names, dilution factors, media details, or anything from
  industry collaborations. Demo data must stay synthetic (`mancq/demo.py`).
* Do not reintroduce pickle for saved files, or any network access in the tool itself.

## Layout

* `mancq/engine.py` fitting and grading (SETTINGS block at the top); `spinsim.py` spin simulation;
  `fidproc.py` raw FID processing; `regrade.py`; `results.py`; `demo.py` synthetic demo and its checks
  (EXPNO 10, 11, 13, 14); `shortcut.py` Windows desktop shortcut; `__main__.py` command line
  (`mancq` with no arguments opens the window).
* `mancq/gui/app.py` PySide6 window (8 steps: Spectra, Process spectra, Regions, Settings, Run, Results,
  Review, Figures); `gui/runner.py` runs the fit in a separate process; `gui/figures_tab.py` the Figures
  step, drawn by `mancq/figures.py` (templates, forgiving sample sheet reader, own TIFF writer because
  Anaconda's Pillow crashes on compressed TIFFs).
* `mancq/shifttable.py` chemical shift assignment table (chemical_shift_assignments.xlsx), written with
  every results update.
* `tests/`: checks on synthetic data only, run with `python -m pytest tests` (each file also runs on its own:
  pH estimate, figures, shift table, overwrite warning, subset rerun, results tables).
* `mancq/library/`: `compounds.csv` (149 compounds: 123 GISSMO, 25 CASMDB peak lists, 1 literature spin
  system), `gissmo_spin_systems.txt`, `other_spin_systems.txt`, `peaklists.csv`, `peaklists_manifest.csv`.
* `MANC-Q.bat` double-click launcher (downloads uv + private Python into `.runtime`); `constraints.txt` exact
  tested package versions; `packaging/` PyInstaller and Inno Setup build; `docs/` quick start, command line,
  screenshots (from the synthetic demo only).

## Checking a change

* `python -m mancq demo <folder>`: fits the synthetic spectra (10 processed, 11 raw FID, 13 look-alike neighbours,
  14 shifts near the limit + phase error + unknown singlets) and checks itself; it must end with "Demo check
  passed" (exit code 0). About 5 minutes. `--quick` fits only 10 and 13.
* `python -m pytest tests` (also run by `.github/workflows/tests.yml` on every push).
* `python -m mancq regrade <results folder>` should reproduce the same tiers and values.
* For window changes, open it (`python -m mancq`) and look at every step, including with no raw FIDs.

## Will's setup

* University Windows PC. Unsigned programs (the installer, MANC-Q.bat) are blocked by IT; he runs MANC-Q in
  Anaconda (`C:\Anaconda`). Install there with `python -m pip install -e . --no-deps --force-reinstall`
  (no `--no-deps` breaks on Anaconda's own PySide6/shiboken6). Desktop shortcut: `python -m mancq shortcut`.
* Look (agreed with Will, 1.7.0 / 1.8.0): light (white) by default; File > Dark mode switches to near black and
  is remembered (it does not follow Windows). Manchester purple accent, Manchester gold for highlights only (never as text on white);
  colours come from `mancq/gui/theme.py` tokens, not hex values in widgets. Plots inside the window stay white
  in both modes. Keep text contrast at 4.5:1 or better. `MANCQ_THEME=dark|light` overrides the choice.
* Logo: merge 2 chosen by Will (lactate quartet, white line over gold fitted area, purple tile). Drawn by
  `packaging/make_logo.py` (needs PySide6 + Pillow) into `mancq/gui/assets/` (SVG, PNGs, mancq.ico) and the
  docs/help `logo.png`; re-run it after changing the drawing. Never use the University of Manchester crest.

## Known open items

* GitHub tidy-up: there are two old releases (v1.1 with files, v1.1.0 empty).
* Code signing (university IT or SignPath Foundation) so the installer runs on managed PCs.
* Automatic FID phasing can differ from TopSpin by 1 to 3 degrees, moving a few values 10 to 30 %.
