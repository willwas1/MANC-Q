# MANC-Q change log

Each version is kept in its own folder. Every change, however small, gets a new version number.

## 2.0.1 (20 October 2026)

The quick start screenshots (mancq/help/images, also shown by the window's Help and the README) are retaken from MANC-Q 2.0.1 with the synthetic demo. They now show the File menu, the warning when the results folder already has results, the "Load from a file..." and "Paste from Excel" dilution buttons, all four demo spectra and the "QC: check" label above EXPNO 14 in Results. They are 1280 x 880 pixels (the window cannot be shorter any more). Nothing in the calculation has changed.

## 2.0.0 (20 October 2026)

MANC-Q 2.0 brings together the fixes from the review of 1.8.0 (versions 1.8.1 to 1.11.1 below). In short:

* **Results you can trust not to be stale:** a folder with results is never silently reused (1.8.1).
* **Bounds are not values:** "Upper bound only" is blank in the value tables, and "Not detected" is 0 mM with a
  plain note on what that means (1.8.2).
* **Better grading:** compounds that are not in the sample no longer take signal from those that are, and a
  reporter peak must agree with the rest of its molecule (1.9.0, 1.9.1). **Rerun data you plan to publish**,
  especially creatine, choline and the 3.19 to 3.27 ppm region.
* **Checks:** a demo that checks itself on four synthetic spectra, 25 automatic tests run on every change, a
  quality flag for every spectrum, and clearer warnings (1.9.1, 1.9.2).
* **Built to last:** settings per run, a results-folder format number, a library check, one version number,
  and the engine and window split into readable parts (1.10.0, 1.11.0).

Why 2.0 rather than 1.12: some changes are not backwards compatible.

* Results folders made by 2.0 (format 2) cannot be opened by MANC-Q 1.x; 2.0 still opens 1.x folders.
* Python 3.12 or later is needed (it was 3.9). Anaconda with Python 3.13 works.
* For your own Python scripts: the engine's functions now take the settings as a keyword argument (`S=`), and
  the defaults moved to `mancq/settings.py`. `mancq.run(...)` works as before.

The numbers in this version are those of 1.11.1.

## 1.11.1 (20 October 2026)

* **The window over which each line is computed is now set in Hz (`LINE_WINDOW_HZ`, 36 Hz) instead of ppm**,
  so it is the same physical width at any field. Before, it was 0.06 ppm: 36 Hz at 600 MHz but 48 Hz at
  800 MHz. At 600 MHz this moves the line tails by a hair (0.00001 ppm), which is enough to change some
  values slightly. Demo, before and after: medians unchanged; the largest change of a true compound is
  1.6 % (choline in the hardest spectrum, EXPNO 14), and every other is under 1.1 %. Two tiers changed:
  threonine in EXPNO 13 is now Quantified, phenylalanine in EXPNO 14 now Overlapped. The demo check passes.
  Spectra at other fields change a little more, because the window there was not 36 Hz before.

## 1.11.0 (20 October 2026)

Split and slim. No change in results: the demo gives exactly the same values, digit for digit, as 1.9.2
and 1.10.0.

* **The engine in parts.** `engine.py` (2,500 lines) is now `settings.py` (the defaults, with the SETTINGS
  block at the top, and the settings of a run), `io.py` (reading spectra and files), `library.py`,
  `fit.py`, `grade.py` and `report.py`. `engine.py` runs them for a folder of spectra and still offers every
  name, so scripts that use `mancq.engine` keep working.
* **The window in parts:** one file per step (`gui/spectra_step.py` ... `gui/review_step.py`), shared
  helpers in `gui/common.py`, and the main window in `gui/app.py`.
* **One copy of the quick start and its screenshots,** in `mancq/help/` (also the window's Help); the copies
  in `docs/` are gone and the README links to it.
* **Errors that were caught and ignored are now logged** (Python's `logging`, name "mancq"): in the window
  to `mancq.log` next to its remembered settings, on the command line to the screen.
* **One colour table for the tiers** (`TIER_COLOURS`), used by the window and the concentration profiles
  PDF (whose markers now use the same blue and orange as the window).
* **Fixed numbers moved into SETTINGS with comments:** the line window (`LINE_WINDOW_PPM`), multiplet grouping
  (`MULTIPLET_GAP_HZ`, `MULTIPLET_MAX_SPAN_PPM`), the glucose anomer ratio (`GLUCOSE_ALPHA_FRACTION`) and
  the limits of unassigned singlets (`UNKNOWN_WIDTH_START`, `UNKNOWN_WIDTH_RANGE`, `UNKNOWN_MOVE_HZ`). The
  library cache now notices a change to these.
* Removed the unused `gallery_pdf`; tidied a few lines.

## 1.10.0 (19 October 2026)

Structure for the future. No change in results: the demo gives exactly the same values as 1.9.2.

* **Settings per run.** The settings at the top of `engine.py` are now only the defaults. Each run makes its
  own settings, passes them to every step, saves them with the results, and regrade, the Review step and
  apply-edits read them back from the results folder. Before, there was one shared set that some steps
  changed and then put back, so a step could see another run's settings. A misspelt setting name is now an
  error instead of being ignored. For scripts: `regrade(out, LOQ_SNR=8)` changes a grading setting for one
  regrade.
* **Results folders carry a format number** (`RESULTS_FORMAT` in `run_settings.json`). A folder made by a
  newer MANC-Q is refused with a clear message ("Install MANC-Q x or later to open it") instead of being
  misread. Everything needed to read older folders is now in one place (`mancq/results.py`).
* **The compound library is checked.** Each saved fit now records the compound name of every multiplet, and
  regrade, the Review step and apply-edits refuse, with a clear message, a folder fitted with a different
  compound library from the installed one, instead of mixing the two.
* **One version number,** kept in `pyproject.toml` only. MANC-Q reads it from there (or from the installed
  package's metadata; the Windows app now includes it). `python packaging/set_version.py X.Y.Z` sets it and
  copies it into `CITATION.cff`, `installer.iss` and `constraints.txt`.
* **Dependencies:** the minimum versions in `pyproject.toml` are now the ones actually tested (Python 3.12,
  numpy 2.3, scipy 1.16, pandas 2.3, matplotlib 3.10, nmrglue 0.12, openpyxl 3.1, PySide6 6.9). The GitHub
  tests now run on Python 3.12 and 3.13.
* Every text file is opened with an explicit encoding (UTF-8 for MANC-Q's own files, Latin-1 for TopSpin's)
  and closed straight after, also when something goes wrong.
* New check `tests/test_results_format.py`.

## 1.9.2 (19 October 2026)

* **A quality flag for every spectrum.** `qc_flag` (ok / check / poor) and `qc_reason` in the QC sheet, and
  in the Results step as `QC: check` (amber) or `QC: poor` (red) above the column, with the reasons when you
  hover over it. It looks at the TSP linewidth, the fit (chi2), the share of signal explained, unassigned
  singlets, the histidine pH flag, the residual phase error, the relaxation delay and the noise estimate.
  Thresholds are settings (`QC_...`).
* **Acquisition parameters** from acqus in the QC sheet: D1, AQ, NS, RG, pulse programme, temperature and
  P1. A warning if D1 + AQ is under 10 s (signals with long T1 may read low).
* **TSP bound to protein:** a warning when the TSP line is much broader (1.5 x) than the metabolites
  (`metabolite_fwhm_Hz`), because every concentration is then too high.
* **Residual phase error:** the TSP fit now includes a small dispersive part, and `TSP_phase_deg` gives the
  size of the phase error, so you know which spectra to re-phase (the demo's 1.5 degree error reads 1.6).
* **Dilutions:** a warning for every spectrum without a dilution and every dilution with no spectrum. The
  dilution file may now be CSV or Excel with headers and "EXPNO 10" or "10.0" style names. On the Settings
  step, **Load from a file...** and **Paste from Excel** fill the per-sample table.
* **Noise:** if the spectrum has no signal-free region, MANC-Q warns and estimates the noise from
  point-to-point differences, instead of from the first 2000 points (which could include TSP).
* **Plain messages** for common problems: a 2D experiment, no processed spectrum with this processing
  number (with the ones that exist), a missing procs or 1r file, no TSP/DSS near 0 ppm, an EXPNO that is
  not in the folder.
* Deconvolution estimate: a multiplet now counts only if the compound is at least 10 % of the fitted signal
  under it (`DECONV_MIN_DOMINANCE`). A slight change in the fit had made a proline sliver under glucose
  into an estimate.
* New check `tests/test_qc.py`.

## 1.9.1 (19 October 2026)

* **The demo now checks itself and fails loudly.** `mancq demo` ends with "Demo check passed" or "DEMO CHECK
  FAILED" and a list: every Quantified value must be within 10 % of the truth, every Overlapped value within
  25 %, and no compound that is not in a spectrum may be graded with a value. The exit code is 1 on failure,
  and the Windows self-test uses the same checks. `mancq demo --quick` fits only EXPNO 10 and 13.
* **Two harder synthetic spectra.** EXPNO 13 has creatinine and O-acetylcarnitine without their look-alikes
  creatine and choline. EXPNO 14 has shifts 5.5 to 6.5 Hz from the library, a 1.5 degree phase error and
  three unknown singlets.
* **Fixes these harder cases found:**
  * Ghost compounds are now looked for in two stages. Missing multiplets in the data are checked from the
    first iteration. "Nothing in the fit" is checked once shifts have settled (from the second iteration),
    because otherwise methionine, displaced 6 Hz, was wrongly dropped. A pH-sensitive compound now counts
    as a ghost too when it has moved at least `GHOST_SHIFT_TITRATABLE_HZ` (15 Hz) and the fit gives its other
    multiplets nothing (cystine and cysteine had slid 26 to 33 Hz onto histidine). `GHOST_MAX_DOMINANCE` is
    now 0.6, the same as for a clean reporter.
  * With a pH estimate, histidine's side-chain multiplets get the usual pH-sensitive freedom (6 Hz) instead
    of the 60 Hz meant for its imidazole H2/H4 (which are held at their measured peaks anyway). Before, its
    CH2 could slip by one coupling (13 Hz) and leave two lines for dimethylsulfone to take.
  * Deconvolution estimate: a multiplet now counts only if its own signal is at least twice the local
    misfit (`DECONV_MAX_MISFIT`), so slivers under a neighbour's imperfect fit no longer make an estimate.
* **Automatic tests on GitHub** (`.github/workflows/tests.yml`) on every push and pull request: all of
  `tests/` (now run with `python -m pytest tests`) and the quick demo. New check
  `tests/test_subset_rerun.py`.

## 1.9.0 (19 October 2026)

**This version changes tiers and some values. Rerun any data you plan to publish**, in particular anything
with creatine, choline or the 3.19 to 3.27 ppm trimethylammonium region.

* **Compounds that are not there no longer take signal from those that are.** A library compound that only
  sits on other compounds' peaks is now left out of the fit (setting `DROP_GHOSTS`). That means none of its
  multiplets is dominant, and either the amount its signal implies would put another of its multiplets
  clearly above the noise where the spectrum shows almost nothing, or (not for pH-sensitive compounds) it
  has moved at least `GHOST_SHIFT_HZ` (5 Hz) from its library position and the fit gives its other
  multiplets almost nothing. Its Long_table `method` says why, and it is graded Upper bound only or Not
  detected. In the demo this was what made choline 18 % low (O-acetylcarnitine sat on its singlet) and
  creatine 5 % low (creatinine sat on its singlet).
* **A reporter is checked against the rest of its molecule.** Quantified now also needs the molecule's
  other multiplets above LOQ (whatever their overlap) to give a median value within `MAX_REPORTER_DISAGREE`
  (20 %) of the reporter value. A lone singlet reporter needs another multiplet of the molecule to agree
  (`SINGLET_NEEDS_SECOND`). Molecules with no other multiplet, such as acetate and formate, are exempt,
  because there is nothing to confirm them with. Otherwise the value is graded Overlapped and its `method`
  says why.
* **Each multiplet's own fitted amount is checked**, not only the final value: a multiplet whose own amount
  would put the molecule's other multiplets clearly above the noise where they are missing is never used as
  a reporter (new column `implies_missing_multiplets` in Peak_list_by_metabolite). If that compound ends up
  without a value, the compound whose peak it sits on gets a `check_note` (also shown in the Results step).
* Demo, EXPNO 10 (processed) before and after: choline 0.816 -> 0.999, creatine 0.947 -> 1.000, glutamine
  0.935 (Overlapped) -> 0.998 (Quantified), threonine now Quantified; three absent compounds that were
  graded with a value (oxaloacetate Quantified, succinate Overlapped, carnitine Deconvolution) are now
  Upper bound only. EXPNO 11 (raw FID): unchanged values, methionine and threonine now Quantified.
* The new grading settings are re-applied by `mancq regrade`.

## 1.8.2 (19 October 2026)

* **Upper bounds are no longer stored as measured values.** Up to 1.8.1, a metabolite graded "Upper bound
  only" had its bound written into `conc_mM`, so it appeared as a plain number in `concentrations_long.csv`,
  the workbook's Concentration_mM_all sheet and the mean on the Summary sheet. Now `conc_mM` is blank for
  these; the bound is in `upper_bound_mM` and in Report_mM as `<=value`, as before. The window, the Review
  step, the Figures labels and the concentration profiles still show the bound with `<=`. Its SE is blank too.
* **Not detected stays 0 mM, with a disclaimer** in the README sheet, at the top of Concentration_mM_all
  (the table now starts on row 4) and in the help: "Metabolites graded Not detected are recorded as 0 mM. A
  zero means the signal was below this sample's detection limit (LOD_mM). It may reflect true absence, a
  concentration too low to detect, or a peak MANC-Q could not resolve in that spectrum. It is not proof that
  the metabolite is absent."
* **Summary:** `mean_conc_mM` is over measured values and zeros only, never bounds; new columns
  `n_measured` and `n_not_detected`.
* Older results folders still open: their bounds are read from `conc_mM`, as before.
* New check `tests/test_results_tables.py` (synthetic demo).

## 1.8.1 (19 October 2026)

* **No more silent re-use of old results (please read).** Up to 1.8.0, running again into a results folder
  that already had results re-used the earlier fit of every spectrum, even if the TSP concentration, the
  dilution, the ignored regions or other settings had changed, so the new values could be wrong without any
  sign. A rerun of only some spectra also dropped the others from the tables. If you ever ran twice into the
  same folder with different settings, run those spectra again into a new folder.
* **The window warns when a results folder already has results**, on the Settings step as soon as you choose
  it and again when you click Start: "This folder already has MANC-Q results. Running here will overwrite
  them." Choose **New subfolder** (the default, e.g. `run_2026-10-19_1340` inside it), **Overwrite** (deletes
  the old results, including manual adjustments and spectra processed by MANC-Q) or **Cancel**. If the run in
  that folder was stopped or crashed, **Carry on the stopped run** fits only what is missing.
* **Command line:** `mancq run` stops with the same message unless the folder is new, `--overwrite` is given,
  or `--resume` is given to carry on a stopped run. See `docs/COMMAND_LINE.md`.
* **Carrying on a stopped run is checked.** A spectrum is only re-used if it belongs to the same run and
  nothing that affects its values has changed: the spectrum and its processing, the TSP concentration, its
  dilution, every fitting and grading setting, the compound library and the MANC-Q version (all written to
  `_per_sample/<EXPNO>_source.json`). Otherwise it is fitted again and the log says why.
* New check `tests/test_overwrite.py` (synthetic demo).

## 1.8.0 (16 October 2026)

* **File menu and a dark mode switch.** A **File** menu in the top left corner: Open a results folder,
  **Dark mode**, Help, About and Quit. The window is now light (white) by default; tick **Dark mode** for
  the black window. The choice is remembered, and MANC-Q no longer follows the Windows setting. Switching
  redraws the tables and the plot toolbars in the new colours.
* **The MANC-Q logo.** The lactate CH quartet as a white measured spectrum over its gold fitted area, on a
  Manchester-purple tile. It is the window and taskbar icon, sits beside the steps at the top of the window,
  heads the Help and About windows, the desktop shortcut (`python -m mancq shortcut`, make it again to
  get the new icon), the installer and the .exe, and the top of README and the quick start guide. All logo
  files are drawn by `packaging/make_logo.py` (in `mancq/gui/assets/`).
* Help pictures smaller than the help window are no longer stretched.

## 1.7.0 (15 October 2026)

A new look for the window. Nothing about the fitting, results or files changes.

* **Light and dark mode.** The window now follows the Windows light / dark setting (it used to force light),
  with a white or near-black background, and changes while it is open. The spectra plots stay white in
  both, like paper. `MANCQ_THEME=dark` or `light` overrides the detection.
* **Purple and gold.** Manchester purple is the accent (main buttons, headings of each box, the selected
  step, ticked boxes, the progress bar); Manchester gold marks the selected step and the line above the
  status bar. Rounded buttons, boxes and fields; flat tabs.
* All colours come from one place (`mancq/gui/theme.py`), with text contrast of at least 5 : 1 in both
  modes (status colours such as green, amber and red are lightened in dark mode).
* Help: a note on light and dark mode; screenshots of steps 6 and 8 updated.

## 1.6.0 (15 October 2026)

* **Chemical shift assignment table for publication.** Every run (and every regrade or saved Review
  adjustment) writes `chemical_shift_assignments.xlsx`: a formatted supplementary table (title, ruled
  header, footnote) with one row per metabolite, its 1H chemical shifts (fitted multiplet centres, median
  over the samples where it was reported), `*` on the peaks used to quantify it, `**` on peaks used but
  shared with another metabolite (Overlapped, semi-quantitative), the number of samples, and the reference
  spectrum (GISSMO bmse number or CASMDB peak list). The same table is the new `Shift_assignments`
  sheet of the workbook. It lists metabolites Quantified, Overlapped or Manually adjusted in at least 25 %
  of the samples (new settings `SHIFT_TABLE_TIERS`, `SHIFT_TABLE_MIN_FRACTION`; the footnote states the
  rule). No 13C shifts: MANC-Q fits 1H spectra only. New module `mancq/shifttable.py`.
* **Figures: colour of the fit.** **Fit colour...** sets the colour of the fitted spectrum (fit evidence, and
  **Draw the fit as well**); the default is now Manchester purple (#660099, was orange). The caption
  names the colours used.
* New test `tests/test_shift_table.py` (synthetic data only).

## 1.5.1 (15 October 2026)

Figures step: colours and offsets.

* **A colour for each spectrum.** Click a spectrum in the list and **Colour...** (or double-click it); a
  small square shows each spectrum's colour, and the caption names the colours. **Reset colours** goes
  back to the template's own (black; colour-blind safe colours in overlays). Works in every template (in
  fit evidence it colours the measured spectrum).
* **Offsets.** Time course and overlay have **Offset upwards** (a share of the plot height; the overlay
  now has it too, default 0) and **Offset sideways** (ppm per spectrum, for a staggered, angled stack;
  peak labels and cut-off marks move with the spectrum). The caption gives both offsets and, with a
  sideways offset, says the ppm scale belongs to the first spectrum.

## 1.5.0 (15 October 2026)

A new step for publication figures, and the Results table now says where each value is measured. The fitting
itself is unchanged.

**Figures (new step 8, `mancq/figures.py`, `mancq/gui/figures_tab.py`)**

* Four templates, following the figure forms seen in published 1H NMR papers on media and biofluids:
  representative spectrum (water region cut out as an axis break, aromatic region magnified xN, numbered or
  abbreviated peak labels with the key in the caption, tall peaks cut off and marked ~); fit evidence (measured
  spectrum and total fit, each chosen metabolite on its own row on the same scale, residual); time course
  (spectra stacked with a small offset on one scale); overlay in colour with a zoomed inset and a dashed box
  marking its region.
* Spectra are the ones MANC-Q fitted, each divided by its own TSP area so all share one vertical scale; the
  smooth fitted baseline is subtracted by default (said in the caption).
* Groups: replicates of one group and time point become one line, either their mean or one chosen replicate.
  Groups are typed in (or pasted from Excel) in **Samples and groups...**, or read from a sample sheet (CSV or
  Excel) that is forgiving: column names in many spellings (EXPNO or sample, group or condition or reactor,
  time or hours or day, replicate or rep), times as 12, 12 h, 2 d or day 3, any delimiter or encoding; samples
  listed twice, samples not in the results and unreadable times are reported and the rest is used; extra
  bioprocess columns are kept. Saved as `sample_sheet.csv` in the results folder.
* Labels come from MANC-Q's own fit (metabolites in the tiers you tick, at their fitted peak positions; one
  per metabolite by default), are moved up with leader lines where crowded, get brackets over wide
  multiplets, and can be edited: hide, change the text, add symbols (* contaminant, ~ cut-off peak,
  # 13C satellite, each explained in the caption) or your own labels.
* Export at journal column widths (85, 120, 180 mm or custom) as PDF, SVG, EPS, TIFF (1200 dpi, Deflate
  compressed, written by MANC-Q itself because the TIFF writer of some Pillow builds, including Anaconda's,
  crashes on compressed TIFFs) and PNG (600 dpi), optionally with the caption, the data behind the figure
  (CSV), a processing record and a Python script that redraws it. A preflight line checks the smallest text
  and thinnest line at the printed size.
* A draft caption writes itself from the settings (field, pulse sequence, reference, cut and magnified
  regions, label key) and can be edited and copied.

**Where each value is measured**

* New column `reporter_ppm` in `concentrations_long.csv` and the `Long_table` sheet: the fitted centre (ppm)
  of the reporter peak(s) each Quantified or Overlapped value comes from (blank for other tiers and for manual
  values).
* The Results step has a first column **Peak used (ppm)**: the usual reporter peak(s) of each metabolite
  (median over samples); hovering over a value shows the peak(s) used in that sample. Results from earlier
  versions take this from the peak list.

* New tests on synthetic data only: `tests/test_figures.py`.

## 1.4.0 (15 October 2026)

MANC-Q now estimates the pH of every sample from its spectrum, and uses it to keep histidine on the right peaks.

* **Estimated pH.** The histidine ring singlets H2 (looked for in 7.70 to 7.96 ppm) and H4 (7.00 to 7.12 ppm)
  move with pH. Each peak top is found on the TSP-referenced axis, refined with a 3-point parabola, used only if
  at least 10 x noise, and converted to pH with Henderson-Hasselbalch (pKa 6.05; H2 8.66 / 7.70 ppm and H4
  7.39 / 7.02 ppm protonated / deprotonated, approximate literature values for free histidine in water, 25 C).
  `pH_est` is the mean of the two. If a window holds more than one peak (for example the tryptophan doublet at
  7.74 ppm) and the tallest peaks disagree, the pair of peaks that agree is used. New columns in
  `qc_per_sample.csv` and the QC sheet: `pH_est`, `pH_H2`, `pH_H4`, `his_H2_ppm`, `his_H4_ppm`, `pH_flag`
  (`ok`, `one_proton_only`, `protons_disagree` when more than 0.25 apart, `no_histidine`). The constants are
  settings (`PH_HIS_PKA`, `PH_HIS_H2_LIMITS`, `PH_HIS_H4_LIMITS`, ...) so a titration in your own medium can
  replace them. It is the pH in the NMR tube at acquisition and an estimate: trends are reliable, absolute pH
  about plus or minus 0.2. The README sheet, the run log, the help and README.md say so.
* **Results and Review steps** show `pH ~7.05 (estimate)` under each sample name (hover for details).
* **Histidine placed from its measured peaks** (setting `PH_PLACE_HISTIDINE`, on by default). When the pH
  could be estimated, histidine H2 and H4 are kept within 3 Hz of the peaks the estimate used (25 Hz of the
  position predicted from the other proton if one peak was too weak), instead of the 60 Hz freedom per
  multiplet they had, which let them jump onto a neighbouring peak. The shift of the molecule as a whole (which
  places its CH2 and CH peaks) is still chosen as before. Without an estimate nothing changes. The fitting of
  other compounds is unchanged, but because histidine now sits differently in some spectra, nearby and weak
  compounds can move a little: Quantified values hardly at all, weak ones sometimes change tier.
* New test on synthetic spectra only: `python -m pytest tests` (or `python tests/test_ph_estimate.py`).

## 1.3.1 (9 October 2026)

* **Try the demo** now rebuilds a demo folder made before 1.3 (it had no 2D), so the demo TOCSY (EXPNO 12) is
  always there, and its message says how to see it (click **2D** in the Results step).

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
