# MANC-Q change log

Each version is kept in its own folder. Every change, however small, gets a new version number.

## 1.3.0 (9 October 2026)

The Results step can now show the whole spectrum and 2D spectra. The fitting and grading are unchanged.

* **Whole spectrum and zoom in Results.** New buttons above the plot: **This metabolite** / **Whole spectrum**
  (the whole spectrum with the selected metabolite filled in purple and arrows over its peaks), and a **Peak**
  list to zoom to the main peak, all peaks or any single multiplet. The mouse wheel zooms around the pointer and
  a toolbar gives box zoom, pan and back to the start.
* **2D spectra (TOCSY, COSY) in Results** (new `mancq/twod.py`). The **2D** button shows the 2D spectrum of the
  sample as contours, with the 1D and the fitted metabolite above it. The 2D is found automatically (same data
  folder, same title, nearest EXPNO; a nearby 2D with a different title is shown with a warning) or chosen with
  **Choose 2D spectrum...** (kept in `_per_sample/twod_links.json`). TopSpin submatrix `2rr` files are read
  directly, and the 2D is re-referenced to its TSP diagonal peak. For the selected metabolite MANC-Q predicts
  its cross peaks from the library spin system (couplings of 2 Hz or more, at the shifts fitted in the 1D,
  including the sample's library offset) and marks each as seen (a peak top at least 10 x the 2D noise within
  0.015 ppm) or missing. **2D check of all metabolites...** gives this for every metabolite of the sample and
  saves `twod_check_sample_<EXPNO>.csv`. This is an identity check only; no value is changed. On real 800 MHz
  TOCSY data, quantified metabolites had a median of all their cross peaks seen, upper bounds about a quarter.
* The synthetic demo has a third experiment, EXPNO 12: a TOCSY of the same mixture as EXPNO 10, so the 2D view
  can be tried. 2D experiments are never taken as 1D samples.

## 1.2.0 (9 October 2026)

Two new features and a set of improvements to the window. The fitting and grading themselves are unchanged.

**Processing controls for every spectrum, not only raw FIDs**

* Step 2 of the window is now **Process spectra** and is always shown. Spectra processed in TopSpin are still
  used exactly as TopSpin processed them by default ("Use TopSpin's processing as is"). Unticking it lets you
  change the phase (zero and first order), line broadening and baseline correction, starting from TopSpin's own
  result: phase changes are applied to TopSpin's real and imaginary parts (`1r`, `1i`), so "no change" gives
  exactly TopSpin's spectrum, and extra line broadening is applied through the time domain. If `1i` is missing
  but the raw FID is there, the FID is processed with TopSpin's line broadening and phased to match TopSpin's
  `1r` first (`mancq/fidproc.py`: `prepare_topspin`, `broaden`, `process`). The grey trace shows TopSpin's
  spectrum for comparison.
* Processing now runs inside the engine (setting `PROCESSING`), so the window and the command line use the same
  code. Processed spectra are written to `<results>/_processed/<EXPNO>` (was `_fid_processed`), with a
  `processing.txt`. Every run that processes spectra writes `processing_used.csv`; the new command line option
  `mancq run ... --processing processing_used.csv` repeats that processing exactly.
* The QC sheet has a `processing` column saying how each spectrum was processed.
* Re-using finished spectra (resume) now checks that the spectrum and its processing are unchanged; if the
  processing changed, the spectrum is fitted again.
* The synthetic demo's processed spectrum (EXPNO 10) now has an imaginary part (`1i`), as TopSpin writes, so
  the adjustment can be tried on it. **Try the demo** rebuilds a demo folder made by an older version (no `1i`).

**Manual review and adjustment of fitted compounds, as in Chenomx (new step 7, Review)**

* Pick a sample and a compound: it is drawn over the measured spectrum with the total fit, its neighbours and
  the residual. Drag it up or down to change its concentration and sideways to move it, or type the
  concentration, shift and linewidth; single multiplets can be moved on their own. **Best fit for this
  compound** fits it alone with everything else fixed; **Refit overlapping compounds** keeps it fixed and refits
  the compounds under its peaks (shown for approval first). **Set to zero**, **Back to automatic** and a note
  field are there too. Double-clicking a value in step 6 opens it here.
* Adjusted values get a new tier, **Manually adjusted** (teal; `number (manual)` in Report_mM), with the
  automatic tier and value kept in new `auto_tier` and `auto_conc_mM` columns. The workbook has a new
  `Manual_edits` sheet (also `manual_edits.csv`) listing every change with the note, user name and time.
* Adjustments are saved separately in `_per_sample/<EXPNO>_manual.json` (plain JSON; the automatic fit is never
  changed) and applied whenever the tables are written, including by `regrade`. The new command
  `mancq apply-edits <results folder>` rebuilds the tables and overlays from them. Fitting a spectrum again sets
  its adjustments aside (`<EXPNO>_manual_set_aside_<time>.json`), because they refer to the earlier fit.
* The saved fit state now also stores the lineshape (`eta`), which the Review step needs; results from 1.1.1 and
  1.1.2 can still be reviewed (the lineshape is then worked out from the saved fit).
* Settings step: "Spectra fitted at the same time" is now called **CPU cores to use**, with a note saying how
  many cores the PC has and when more cores help (the setting itself is unchanged; `--workers` on the command
  line).
**Window improvements**

* Remembers your settings between sessions (dataset folder, PROCNO, reference, dilution, cores, regions,
  notification choice) in `window_settings.json` in the AppData\MANC-Q folder, and lists recently opened
  results folders in step 6 (**Recent results...**).
* Review: new **Fit** column (|measured - fit| around a compound's peaks as a share of its own signal; red above
  25 %). Every column sorts when its heading is clicked, again to reverse, as in File Explorer; blanks stay at
  the bottom. **Undo** / **Redo** (Ctrl+Z, Ctrl+Y; a drag or a run of typing on one compound is one step) and
  Ctrl+Up / Ctrl+Down for the previous / next compound.
* Results: Ctrl+C copies the selected values with metabolite and sample names; right-click offers
  **Copy numbers only** (for Excel or Prism). **Save this plot...** saves the plot as PNG (300 dpi), SVG or PDF.
  Full sample titles on hover.
* Process spectra: **Auto-phase all raw FIDs**, **Mark all as checked** (asks first), and arrow keys for 0.1
  degree phase steps (Page Up / Page Down: 1 degree). The line-broadening button is greyed out with the boxes
  when a spectrum is used as TopSpin processed it.
* Spectra: select several rows to overlay them, each divided by its reference peak (can be switched off);
  scroll to zoom.
* Settings: **CPU cores to use** (was "Spectra fitted at the same time"), with how many cores the PC has;
  **Tell me when a run finishes** (Windows notification, sound and taskbar flash). Run progress is also shown in
  the window title.
* A **Help** button on every step shows its part of the quick start guide, which now ships inside the package
  (`mancq/help/`, a copy of `docs/QUICKSTART_GUI.md` and its pictures).
* Error messages have **Copy error details** (version, Python, system and the full error, ready to paste), and
  unexpected errors are now shown instead of disappearing silently.
* Wider columns and full text on hover in the tables; screenshots in `docs/images` renewed.
* New module `mancq/manual.py`; `mancq/regrade.py` has a shared `restore_settings`.
* README, `docs/QUICKSTART_GUI.md` and `docs/COMMAND_LINE.md` describe both features.
* Version number raised to 1.2.0 in `pyproject.toml`, `mancq/engine.py`, `CITATION.cff`,
  `packaging/installer.iss` and `constraints.txt`.

## 1.1.2 (7 October 2026)

Compared with the 1.1.1 files on GitHub (uploaded 6 October):

* New command `python -m mancq shortcut` puts a MANC-Q shortcut on the Windows desktop that starts MANC-Q with
  the same Python (Anaconda or normal Python), with no console window (`mancq/shortcut.py`, `mancq/__main__.py`).
* If MANC-Q fails to start from a shortcut, a message box now shows the error instead of nothing happening
  (`mancq/__main__.py`).
* Every results workbook now records the MANC-Q version that produced it (first line of the README sheet of
  `metabolite_concentrations.xlsx`), and `_per_sample/run_settings.json` stores it as `MANCQ_VERSION`
  (`mancq/engine.py`). The fitting itself is unchanged.
* README: desktop shortcut step for university PCs; update instruction now uses `--force-reinstall` (without it
  pip can skip the update); the library offset is described as "a few hertz" instead of a figure from our own
  samples. `docs/COMMAND_LINE.md`: shortcut command.
* `examples/dilution_example.csv`: example factor changed from 5.5 to 2.0 (5.5 was our own dilution).
* Version number raised to 1.1.2 in `pyproject.toml`, `mancq/engine.py`, `CITATION.cff`,
  `packaging/installer.iss` and `constraints.txt`.

## 1.1.1 (6 October 2026)

Window always uses a light theme (readable in Windows dark mode) with higher-contrast tier colours; steps are
numbered 1, 2, 3... with no gap when the Process FIDs step is hidden; DMSO library entry changed from an SDBS
peak list to a spin system from published shifts (`library/other_spin_systems.txt`); saved fit state uses JSON
instead of pickle; installer builds attach SHA256 checksums; README: "How it works" with GISSMO, university-PC
(Anaconda) route, "Is it safe?" section.

## 1.1.0 (6 October 2026)

Graphical interface, raw FID processing, regions to ignore, Windows launcher and installer build.
