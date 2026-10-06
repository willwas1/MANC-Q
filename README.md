# MANC-Q

**M**etabolite **A**bundance by **N**MR **C**urve-fitting, **Q**uantified. Developed at the University of Manchester.

Absolute quantification of metabolites in 1D ¹H NMR spectra by fitting the whole spectrum with simulated
compound spectra.

Every compound in the library is described by its GISSMO spin system (chemical shifts and J couplings, from
BMRB). The program simulates each one exactly at the field of your spectrometer, so multiplet shapes are
correct at 400, 600 or 800 MHz, including second-order patterns. It then finds the set of amounts that best
reproduces the measured spectrum, allowing each compound and each multiplet a few hertz of shift and some
linewidth freedom. Amounts are converted to mM using only the internal standard (TSP or DSS). No spiked
standards and no calibration against other samples are needed.

Every result carries a confidence tier, so it is clear which numbers are measurements and which are estimates.

## What you need

* Bruker processed 1D ¹H spectra (phased and baseline-corrected in TopSpin): `<EXPNO>/pdata/1/1r` and `procs`.
* A TSP or DSS internal standard of known concentration in the NMR tube.
* Optionally, the dilution from your original sample to the tube, if you want concentrations in the original sample.

That is all. Sample names are taken from the EXPNO folder names; nothing else is read from the acquisition.

## Install

Python 3.9 or later.

```
git clone https://github.com/willwas1/manc-q.git
cd manc-q
pip install -r requirements.txt
```

Or install it as a command: `pip install .` (this gives you the `mancq` command in place of `python -m mancq`).

## Quick start

Check that everything works: this makes a synthetic 600 MHz spectrum of 18 metabolites at known
concentrations, quantifies it, and prints the recovered values next to the true ones (about 2 minutes).

```
python -m mancq demo
```

Quantify your own spectra:

```
python -m mancq run  path/to/experiment_folder  path/to/results  --ref-mm 0.5
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
| `--pluronic yes/no` | model Pluronic F-68 (a surfactant in many cell-culture media); default: detected automatically |
| `--no-overlays` | skip the overlay PDFs |
| `--workers 4` | number of spectra fitted in parallel |

A fit takes roughly 1 to 5 minutes per spectrum on one core. Runs resume: if you stop a run and start it again,
finished spectra are reused.

## Output

In the results folder:

* `metabolite_concentrations.xlsx`: start here. The `README` sheet explains everything; `Report_mM` has one
  row per metabolite and one column per sample, each value written according to its tier (see below).
* `overlay_sample_<EXPNO>.pdf`: the measured spectrum, the fitted model and each compound's contribution, in
  17 windows. Look here first whenever a number is surprising.
* `concentrations_long.csv`, `peak_list_by_metabolite.csv`, `observed_peaks_assigned.csv`, `qc_per_sample.csv`:
  the same results as plain tables.
* `run_log.txt`: what was done to each spectrum.

### Confidence tiers

| Tier | Meaning | In `Report_mM` |
|---|---|---|
| Quantified | at least one of the compound's multiplets is resolved: strong, mostly this compound, fitted well | number |
| Overlapped (semi-quantitative) | the best multiplet is shared with other compounds | number |
| Deconvolution estimate | no multiplet is resolved; the value comes from the whole-compound fit and is well determined | `~number` |
| Upper bound only | signal is present where the compound would be, but cannot be attributed to it | `<=number` |
| Not detected | below 3 x noise | `<LOD` |

The tier is decided per spectrum. Read the tiers along with the numbers: a compound "fitted at 1 mM" in the
Upper bound tier has not been measured at 1 mM.

## Changing things

* **Fit and grading settings** are in the `SETTINGS` block at the top of `mancq/engine.py`, each with a comment.
* **Re-grade without refitting**: after changing grading settings, run
  `python -m mancq regrade path/to/results` (seconds rather than minutes per spectrum).
* **Adding a compound**:
  1. Find it at <https://gissmo.bmrb.io> and copy its spin-system matrix (shifts and couplings).
  2. Add one line to `mancq/library/gissmo_spin_systems.txt` in the same format as the others:
     `<entry id>|<name>|<note>#<shift1>,<shift2>,...#<i>-<j>:<J in Hz>;...` (protons numbered from 1).
  3. Add one row to `mancq/library/compounds.csv`. Set `titratable` to 1 for compounds whose shifts
     depend on pH; the last two columns widen the allowed shift for that compound, in Hz.
  4. The library is rebuilt automatically on the next run.
* **Using it from Python**:
  ```python
  import mancq
  mancq.run("experiment_folder", "results", ref_mm=0.5, dilution=1.0)
  ```

## How it works, briefly

1. The reference singlet at 0 ppm is fitted (pseudo-Voigt plus ²⁹Si satellites). It sets the ppm scale, the
   starting linewidth and the concentration scale:
   `c = amplitude per proton x (ref_mM x 9 / reference area) x dilution`.
2. The library is simulated at the spectrometer frequency and split into multiplets, one per proton group.
3. A small global offset between the library and the referenced spectrum is measured from sharp compounds that
   are clearly present (typically +10 to +12 Hz at 800 MHz in cell-culture media).
4. The whole spectrum (0.6 to 9.6 ppm) is fitted as non-negative compound amounts plus a smooth spline baseline.
   Each compound is then slid as a rigid pattern (a few Hz; more for pH-sensitive compounds), and overlapping
   multiplets are refined locally in position and width. This repeats three times. Peaks that no library compound
   explains become free singlets, so they cannot inflate a neighbour.
5. Each compound's concentration is taken from its cleanly resolved multiplets, fitted individually; the tier is
   assigned from signal-to-noise, how much of the local signal is the compound's own, how well it fits, and
   whether every multiplet the compound should produce is actually there.

## Known limitations

* Matrix effects on chemical shift are only partly absorbed by the shift freedom. Citrate in media containing
  Ca²⁺/Mg²⁺ is a known failure (its AB pattern changes shape); histidine and other pH-sensitive compounds move
  between samples.
* Compounds sharing an identical singlet (for example the N(CH₃)₃ group of choline, phosphocholine, betaine,
  acetylcholine and carnitine near 3.2 ppm, or creatine and creatinine near 3.03 ppm) cannot always be separated.
  In the synthetic demo one or two singlet compounds can read 10 to 30 % low while still graded Quantified,
  because a look-alike compound takes part of the singlet. Check the overlays for singlet-only compounds.
* Crowded regions (branched-chain amino acid methyls, 2.0 to 2.5 ppm) are handled, but values there are more
  often Overlapped than Quantified.
* Not corrected: incomplete T1 relaxation (use a long enough recycle delay), weighing error of the reference,
  and reference binding to proteins (serum and plasma need ultrafiltration or an ERETIC-type reference).
* The 26 compounds described by fixed peak lists (see `library/peaklists_manifest.csv`) are exact only near
  800 MHz; at other fields their line spacing is rescaled with a first-order approximation.
* Overlapped and Deconvolution-tier values are sensitive to small numerical differences and can shift by 10 % or
  more between runs on different machines; Quantified values are stable to about 3 %.
* The library was assembled for cell-culture media and similar aqueous samples. Other matrices will need
  compounds added.

## Library data

* `gissmo_spin_systems.txt`: spin-system parameters from GISSMO / BMRB (<https://gissmo.bmrb.io>). Please cite
  Dashti H. et al., *Anal. Chem.* 2017, 89, 12201-12208, and BMRB.
* `peaklists.csv`: fixed peak lists for 26 compounds taken from CASMDB, which itself draws on HMDB, BMRB and
  SDBS spectra (the origin of each is in `peaklists_manifest.csv`). These sources have their own terms of use.

## Citation

If you use this software, please cite it (see `CITATION.cff`) and GISSMO.

## Licence

MIT for the code (see `LICENSE`). The library data files remain subject to the terms of their sources.
