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

Check that everything works: this makes two synthetic 600 MHz spectra of 18 metabolites at known
concentrations (one processed, one raw FID), quantifies them, and prints the recovered values next to the true
ones (a few minutes).

```
mancq demo
```

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
| `--no-overlays` | skip the overlay PDFs |
| `--workers 4` | number of spectra fitted in parallel |

A fit takes roughly 1 to 5 minutes per spectrum on one core. Runs resume: if you stop a run and start it again,
finished spectra are reused.

Re-grade a finished run after changing grading settings (seconds rather than minutes per spectrum):

```
mancq regrade path/to/results
```

From Python or a notebook:

```python
import mancq
mancq.run("experiment_folder", "results", ref_mm=0.5, dilution=1.0)
```
