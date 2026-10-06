# MANC-Q

**M**etabolite **A**bundance by **N**MR **C**urve-fitting, **Q**uantified. Developed at the University of Manchester.

MANC-Q measures metabolite concentrations (mM) in 1D ¹H NMR spectra. It fits the whole spectrum with simulated
spectra of known compounds and uses the TSP or DSS peak as the concentration reference. No spiked standards
and no calibration curves are needed. Every result carries a confidence tier, so it is clear which numbers are
measurements and which are estimates.

![MANC-Q results screen](docs/images/6_results.png)

## Get started on Windows

You do not need Python or anything else installed. Pick one of these.

### Option 1: the installer (simplest)

1. Open the [Releases page](https://github.com/willwas1/MANC-Q/releases) and download **MANC-Q-setup.exe**.
2. Double-click it and follow the steps. No administrator rights are needed.
3. Start **MANC-Q** from the Start menu.

Windows may say "Windows protected your PC" because the installer is not signed. Click **More info**, then
**Run anyway**.

### Option 2: download the folder and double-click

Use this if there is no installer on the Releases page yet, or if you want the Python code as well.

1. On this page, click the green **Code** button, then **Download ZIP**.
2. In your Downloads folder, right-click the ZIP file and choose **Extract All**. Put the folder somewhere easy
   to find, such as Documents.
3. Open the extracted **MANC-Q** folder and double-click **MANC-Q.bat**.

The first time, a black window shows the set-up: it downloads a private copy of Python and the packages MANC-Q
needs into the MANC-Q folder (about 5 minutes and 650 MB; it needs internet). Nothing is installed anywhere else
on the computer and no administrator rights are needed. It then puts a MANC-Q shortcut on your desktop. After
that, MANC-Q opens in a few seconds. To remove it, delete the folder and the shortcut.

If Windows says "Windows protected your PC", click **More info**, then **Run anyway**.

## Your first analysis

MANC-Q opens on the first of six tabs. The full guide with pictures is in
[docs/QUICKSTART_GUI.md](docs/QUICKSTART_GUI.md).

1. **Spectra.** Click **Try the demo** the first time: it makes two synthetic spectra with known
   concentrations, so you can practise. For your own data, click **Browse** and choose the TopSpin dataset
   folder (the folder that contains the numbered experiment folders such as `10`, `11`, `12`).
2. **Process FIDs.** Only appears if some experiments have not been processed in TopSpin. Check the phasing.
3. **Regions.** Parts of the spectrum to leave out. The defaults suit most cell-culture samples.
4. **Settings.** Type the concentration of TSP (or DSS) in the NMR tube, and the dilution if you want
   concentrations in the original sample.
5. **Run.** Click **Start**. Each spectrum takes 1 to 5 minutes.
6. **Results.** A table of metabolites by samples. Click any value to see the fit it came from.
   **Open Excel workbook** gives every table.

## What your data needs to be

* **Bruker** 1D ¹H experiments. Either processed in TopSpin (each experiment folder contains `pdata/1/1r`), or
  raw FIDs (`fid` and `acqus`), which MANC-Q processes itself.
* A **TSP or DSS** internal standard of known concentration in the tube.

```
MyExperiment            <- choose this folder
├── 10
│   ├── fid, acqus ...
│   └── pdata/1/1r, procs
├── 11
└── 12
```

Other vendors' formats are not supported yet.

## Reading the results

The main file is `metabolite_concentrations.xlsx` in the results folder. Its `README` sheet explains every
sheet. Each value has a confidence tier:

| Tier | What it means | Shown as |
|---|---|---|
| Quantified | at least one of the compound's signals is clean and well fitted; treat it as a measurement | a number |
| Overlapped | its best signal is shared with other compounds; semi-quantitative, good for trends | a number |
| Deconvolution estimate | no clean signal; estimated from the fit of the whole compound | `~number` |
| Upper bound only | something is there, but it cannot be attributed to this compound | `<=number` |
| Not detected | below 3 x noise | `<LOD` |
| Not measurable | all its peaks are in a region you chose to leave out | `n/m` |

Always read the tier with the number. "1.2 mM, Quantified" is a measurement; "<=1.2 mM, Upper bound" only says
it is not more than 1.2 mM. When a value looks surprising, open `overlay_sample_<EXPNO>.pdf`, which shows the
measured spectrum, the fit and each compound's contribution.

## If something goes wrong

| Problem | What to do |
|---|---|
| MANC-Q.bat says the download failed | Check the internet connection. University and company firewalls sometimes block the download; try another network once, or use the installer (Option 1). Running MANC-Q.bat again carries on where it stopped. |
| "This file must stay inside the MANC-Q folder" | You opened MANC-Q.bat from inside the ZIP. Extract the ZIP first (right-click, Extract All). |
| No experiments are listed | Choose the folder that contains the numbered experiment folders, not an experiment folder or a single file. |
| Results look wrong | Check the TSP/DSS concentration, the phasing (tab 2, or in TopSpin) and the overlay PDF for that sample. |

To ask for help, open an [issue](https://github.com/willwas1/MANC-Q/issues) with a screenshot of the error and,
if there is one, the `run_log.txt` file from the results folder.

## For Python users

The same program can be used from the command line, from scripts and from notebooks; see
[docs/COMMAND_LINE.md](docs/COMMAND_LINE.md).

* With Python 3.9 or later: `pip install git+https://github.com/willwas1/MANC-Q.git`, then `mancq` (opens the
  window) or `mancq --help`.
* On a Mac or Linux without Python: install uv (`curl -LsSf https://astral.sh/uv/install.sh | sh`), then run
  `uvx --python 3.12 --from git+https://github.com/willwas1/MANC-Q mancq`. MANC-Q is developed and tested on
  Windows; other systems should work but are less tested.
* `constraints.txt` lists the exact package versions MANC-Q was tested with. Small numerical differences between
  package versions can move Overlapped and Deconvolution values; use
  `pip install . -c constraints.txt` if you need identical numbers to someone else.

### Building the Windows installer

You normally do not need to: publishing a release on GitHub (Releases > Draft a new release > choose a tag such
as `v1.1.0` > Publish) starts `.github/workflows/build-windows.yml`, which builds `MANC-Q-setup.exe` and a portable
zip on GitHub's computers, tests them, and attaches them to the release (about 15 minutes; progress is on the
Actions tab). To build on your own PC instead, run `packaging\build_windows.bat`.

## Changing the fitting (for developers)

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
* Phasing matters. The same FIDs processed by MANC-Q (automatic phasing) and by TopSpin gave the same
  concentrations on the median (ratio 1.00), and the MANC-Q-phased spectra fitted slightly better, but a 1 to 3
  degree difference in zero-order phase moved a few individual values by 10 to 30 % (for example choline, and
  glucose in one spectrum). Check the phasing (tab 2 of the window shows it before the run; from the command line, look at the overlays) for any
value you rely on, or process the FIDs in TopSpin.
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
