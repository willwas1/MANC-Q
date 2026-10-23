# MANC-Q from the command line

For scripts, batch jobs, servers, or if you prefer typing to clicking. The window and the command line run the
same code and give the same numbers.

## Starting

* If you set MANC-Q up with **MANC-Q.bat** (Windows): open a Command Prompt in the MANC-Q folder (click the
  File Explorer address bar, type `cmd`, press Enter) and type `MANC-Q.bat` wherever this page says `mancq`,
  for example `MANC-Q.bat run "C:\NMR\Experiment 1" "C:\NMR\Results" --ref-mm 0.5`.
* If you installed it with pip or uv: type `mancq` (or `python -m mancq`).
* `mancq` with nothing after it opens the window. `mancq --help` lists the commands.

Put any path that contains spaces in double quotes, for example `"C:\NMR data\Experiment 1"`.

## Commands

Check that everything works: this makes synthetic 600 MHz spectra of known composition, quantifies them, and
prints the recovered values next to the true ones (about 7 minutes). EXPNO 10 is a processed spectrum, EXPNO 11
a raw FID of the same mixture, EXPNO 13 has creatinine and O-acetylcarnitine instead of their look-alikes
creatine and choline, EXPNO 14 has shifts near the limit, a small phase error and three unknown singlets, and
EXPNO 15 is a late culture (little glucose, much lactate).

```
mancq demo
```

The demo checks itself: every Quantified value must be within 10 % of the truth, every Overlapped value within
25 %, and no compound that is not in a spectrum may be graded with a value. It ends with "Demo check passed",
or "DEMO CHECK FAILED" and the list of what failed (then the command's exit code is 1, so scripts can tell).
`mancq demo --quick` fits only EXPNO 10 and 13 (used by the automatic tests on GitHub).

Quantify your own spectra:

```
mancq run  path/to/experiment_folder  path/to/results  --ref-mm 0.5
```

`experiment_folder` is the TopSpin dataset folder containing the numbered EXPNO folders; a single EXPNO folder
also works. `--ref-mm` is the TSP/DSS concentration in the tube.

Optional settings:

| Option | Meaning |
|---|---|
| `--dilution 5` | multiply every result by this factor (sample-to-tube dilution); default 1 = mM in the tube |
| `--dilution-file dil.csv` | per-sample factors, CSV with columns `sample,dilution` (see `examples/`) |
| `--reference DSS` | label the reference as DSS instead of TSP (the fit is the same) |
| `--procno 2` | read `pdata/2` instead of `pdata/1` |
| `--samples 10 11 12` | only these EXPNOs |
| `--exclude 5.5 6.5` | leave an extra region out of the fit (repeatable); water 4.60-5.00 ppm is always excluded |
| `--pluronic auto/yes/ignore/no` | Pluronic F-68 / poloxamer (cell-culture media): model it if present (default), always model it, leave its regions (3.69-3.74, 1.10-1.20 ppm) out of the fit, or not present |
| `--fids missing/all/never` | raw FIDs: process them where there is no TopSpin-processed spectrum (default), always, or never |
| `--lb 0.3` | line broadening (Hz) used when processing FIDs |
| `--processing proc.csv` | process these spectra in MANC-Q (see below); `processing_used.csv` from an earlier run repeats its processing exactly |
| `--no-overlays` | skip the overlay PDFs |
| `--workers 4` | number of spectra fitted in parallel |
| `--overwrite` | the results folder already has MANC-Q results: delete them and run again (see below) |
| `--resume` | carry on a run that was stopped or crashed in this results folder (see below) |

A fit takes roughly 1 to 5 minutes per spectrum on one core.

### A results folder that already has results

If the results folder already holds MANC-Q results (`concentrations_long.csv` or `_per_sample`), `mancq run`
stops with "This folder already has MANC-Q results. Running here will overwrite them." and runs nothing. Then:

* give a new results folder (the safest), or
* add `--overwrite`: the old results are deleted first, including manual adjustments made in the Review step
  and the spectra MANC-Q processed (`_processed`). Your own files in the folder, the Figures settings and
  `sample_sheet.csv` are kept, or
* add `--resume` to carry on a run that was stopped or crashed: give the same spectra and settings as before.
  Spectra it finished are re-used only if nothing that affects their values has changed (the spectrum and its
  processing, the TSP concentration, the dilution, the ignored regions, the fit range, Pluronic, any other
  fitting or grading setting, the compound library and the MANC-Q version); any that changed are fitted again,
  and `run_log.txt` says why.

Before version 1.8.1, MANC-Q silently re-used earlier results in the same folder, even when the TSP
concentration or other settings had changed.

### Changing the processing of a spectrum

Spectra processed in TopSpin are used as they are. To change some, give a CSV file with `--processing`:

```
sample,source,lb,ph0,ph1,baseline
10,topspin,0.5,1.5,0,False
11,fid,0.3,,,False
```

* `source topspin`: start from TopSpin's spectrum; `ph0` and `ph1` are the change from TopSpin's phase in
  degrees (0 = unchanged), `lb` is the total line broadening in Hz (blank = TopSpin's).
* `source fid`: process the raw FID; `ph0` and `ph1` are the phase itself, blank = automatic.
* `baseline True` adds a spline baseline correction (normally left off).

Every run that processes spectra writes `processing_used.csv` in the results folder in this format.

Put a desktop shortcut (Windows) that opens the window with the Python you ran this from:

```
mancq shortcut
```

Re-grade a finished run after changing grading settings (seconds rather than minutes per spectrum):

```
mancq regrade path/to/results
```

Rebuild the tables and overlays after manual adjustments (made in the Review step of the window and saved in
`_per_sample/<EXPNO>_manual.json`); nothing is refitted:

```
mancq apply-edits path/to/results
```

From Python or a notebook:

```python
import mancq
mancq.run("experiment_folder", "results", ref_mm=0.5, dilution=1.0)
# a results folder that already has results: existing="overwrite" or existing="resume" (as on the command line)
# any setting can be changed for one run by name, e.g. a stricter quantification limit:
mancq.run("experiment_folder", "results_loq20", ref_mm=0.5, LOQ_SNR=20)
```

From version 2.0 each run has its own settings (`mancq.Settings`); the defaults in `mancq/settings.py` are
never changed by a run, and a misspelt setting name is an error.
