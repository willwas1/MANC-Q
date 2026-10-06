# MANC-Q

## Metabolite Abundance by NMR Curve-fitting, Quantified

**MANC-Q** is software developed at the **University of Manchester** for measuring metabolite concentrations from 1D proton (^1H) NMR spectra.

It is designed to automatically analyse NMR spectra and estimate the concentrations of metabolites such as amino acids and other compounds.

This guide is written for people who **are not familiar with Python, Git, or the Windows Command Prompt**.

---

# Before you start

MANC-Q currently runs using Python.

You will need:

* A Windows PC
* Python 3.9 or newer
* Bruker TopSpin processed 1D ^1H NMR data
* A TSP or DSS internal standard with a known concentration

You **do not need to know how to program** to use MANC-Q.

However, the first installation requires a few steps. Once MANC-Q has been installed, you can create a shortcut so that you do not need to type commands every time.

---

# What does MANC-Q do?

MANC-Q takes a processed NMR spectrum and compares it against a library of simulated metabolite spectra.

It then estimates how much of each metabolite is present.

The software produces results including:

* Metabolite concentrations
* Quality/confidence information
* Graphs showing how well the spectrum was fitted
* Excel spreadsheets containing the results
* CSV files containing the results
* A log showing what the software did

The main results file is:

**`metabolite_concentrations.xlsx`**

Start with this file when you want to look at your results.

---

# Step 1 — Install Python

If Python is not already installed on your computer:

1. Go to:

   https://www.python.org/downloads/

2. Download the latest Python 3 version.

3. Start the installer.

4. **IMPORTANT:** Before clicking Install, tick:

   **Add Python to PATH**

5. Complete the installation.

### Checking Python

You can check whether Python is installed by opening Command Prompt and typing:

```text
python --version
```

You should see something similar to:

```text
Python 3.12.5
```

If you see a Python version number, Python is installed correctly.

---

# Step 2 — Download MANC-Q

You have two choices.

## Option A — Download without using Git

This is recommended for beginners.

1. Go to:

   https://github.com/willwas1/MANC-Q

2. Click the green **Code** button.

3. Click:

   **Download ZIP**

4. Windows will download a ZIP file.

5. Open your **Downloads** folder.

6. Right-click the downloaded ZIP file.

7. Select:

   **Extract All...**

8. Extract the folder somewhere easy to find.

For example:

```text
Documents\MANC-Q
```

You should now have a folder containing files and folders such as:

```text
MANC-Q
├── examples
├── mancq
├── README.md
├── requirements.txt
└── pyproject.toml
```

---

# Step 3 — Install MANC-Q

This is the only part where you need to use Command Prompt.

Don't worry — you do not need to understand programming.

### Open the MANC-Q folder

Find the folder where you extracted MANC-Q.

For example:

```text
Documents\MANC-Q
```

Click in the **address bar** at the top of Windows File Explorer.

Type:

```text
cmd
```

and press **Enter**.

A black Command Prompt window will appear.

It should already be opened in the MANC-Q folder.

---

# Step 4 — Install the required software

Copy and paste this command into the black window:

```text
python -m pip install -r requirements.txt
```

Then press **Enter**.

Windows will install the software MANC-Q needs.

This may take a few minutes.

### If Windows says that `python` is not recognised

Try:

```text
py -m pip install -r requirements.txt
```

If that also doesn't work, Python has probably not been installed correctly.

---

# Step 5 — Test MANC-Q

Before using your own NMR data, it is recommended that you test the program.

In the Command Prompt window, type:

```text
python -m mancq demo
```

Then press **Enter**.

MANC-Q will run a demonstration using a simulated NMR spectrum.

The demonstration takes approximately **2 minutes**.

If it completes successfully, MANC-Q is installed correctly.

---

# Using your own NMR data

## What do I need?

MANC-Q expects **processed Bruker 1D ^1H NMR spectra**.

The spectra should already have been:

* Phased
* Baseline corrected
* Processed in Bruker TopSpin

MANC-Q expects the processed spectrum files to be in the normal Bruker structure.

For example:

```text
MyExperiment
│
├── 10
│   └── pdata
│       └── 1
│           ├── 1r
│           └── procs
│
├── 11
│   └── pdata
│       └── 1
│           ├── 1r
│           └── procs
│
└── 12
    └── pdata
        └── 1
            ├── 1r
            └── procs
```

The numbers (`10`, `11`, `12`, etc.) are the **EXPNO** numbers from TopSpin.

---

# Internal standard

MANC-Q uses an internal standard to calculate the concentrations.

You need either:

* **TSP**, or
* **DSS**

in your NMR tube.

You must know the concentration of the internal standard.

For example:

```text
TSP concentration = 0.5 mM
```

You will give this value to MANC-Q when starting an analysis.

---

# Running your first analysis

The basic command is:

```text
python -m mancq run YOUR_NMR_DATA YOUR_RESULTS_FOLDER --ref-mm 0.5
```

For example:

```text
python -m mancq run "C:\NMR\Experiment1" "C:\NMR\MANCQ Results" --ref-mm 0.5
```

Replace:

```text
C:\NMR\Experiment1
```

with the location of your NMR data.

Replace:

```text
0.5
```

with the concentration of your TSP or DSS internal standard.

---

# A simpler way to run MANC-Q

If you are uncomfortable typing commands, you can create a small Windows file that starts MANC-Q for you.

## Creating a MANC-Q launcher

Open **Notepad**.

Copy the following:

```bat
@echo off
title MANC-Q
echo.
echo ==============================
echo        MANC-Q
echo ==============================
echo.
echo Starting MANC-Q...
echo.
python -m mancq
pause
```

Save the file as:

```text
MANC-Q.bat
```

Make sure Windows does **not** save it as:

```text
MANC-Q.bat.txt
```

You can then double-click the `.bat` file.

---

# Important: running an analysis

MANC-Q needs to know:

1. Where your NMR data is
2. Where you want the results saved
3. The concentration of your TSP/DSS reference

For example:

```text
NMR data:
C:\NMR\Experiment1

Results:
C:\NMR\MANC-Q Results

TSP:
0.5 mM
```

The command would be:

```text
python -m mancq run "C:\NMR\Experiment1" "C:\NMR\MANC-Q Results" --ref-mm 0.5
```

If you regularly analyse data, this can be turned into a **one-click Windows launcher**.

---

# How long does an analysis take?

A typical spectrum takes approximately:

**1–5 minutes per spectrum on one CPU core.**

For example, 20 spectra could take considerably longer than a single spectrum.

MANC-Q can process multiple spectra in parallel if you increase the number of workers.

For example:

```text
--workers 4
```

means that up to four spectra can be processed at the same time.

If you are not comfortable with computer settings, you can simply leave this option alone.

---

# Where are my results?

After MANC-Q finishes, look inside the results folder you specified.

You should find files including:

### `metabolite_concentrations.xlsx`

This is the main results file.

Open it using:

* Microsoft Excel
* LibreOffice Calc
* Another program capable of opening `.xlsx` files

The workbook contains a `README` sheet explaining the results.

---

# Understanding the results

MANC-Q gives each result a confidence category.

## Quantified

This is the strongest result category.

The metabolite has at least one clearly resolved signal that fits well.

The reported concentration can generally be treated as a measured value.

---

## Overlapped

The metabolite signal overlaps with signals from other compounds.

The result is therefore **semi-quantitative**.

Use more caution when interpreting these values.

---

## Deconvolution estimate

The software could not find a clean, isolated signal for the compound.

Instead, the concentration was estimated from the overall spectrum.

These values should be treated as estimates rather than direct measurements.

---

## Upper bound only

The software sees a signal in the area where the metabolite could be present, but cannot confidently determine how much of that signal belongs to the metabolite.

The result is therefore an **upper limit**, rather than a measured concentration.

---

## Not detected

The metabolite was not detected above the software's detection threshold.

The results table will show:

```text
<LOD
```

---

# Very important: don't look only at the number

A result such as:

```text
1.2 mM
```

does not tell the whole story.

Always look at the confidence category as well.

For example:

```text
Quantified       1.2 mM
```

is very different from:

```text
Upper bound      <=1.2 mM
```

The second result does **not** mean that the sample contains exactly 1.2 mM.

---

# Checking a surprising result

MANC-Q creates PDF files called:

```text
overlay_sample_<EXPNO>.pdf
```

These show:

* The measured NMR spectrum
* The fitted spectrum
* Individual compound contributions

If a result looks surprising, these PDF files are a good place to investigate what happened.

---

# If you stop the program

Don't panic.

MANC-Q is designed so that completed spectra can be reused.

If you stop a run and start it again, finished spectra do not normally need to be processed again.

---

# Optional settings

MANC-Q has additional options for experienced users.

| Option               | What it does                                  |
| -------------------- | --------------------------------------------- |
| `--dilution 5`       | Corrects results for sample dilution          |
| `--reference DSS`    | Uses DSS as the reference label               |
| `--procno 2`         | Uses processed data in `pdata/2`              |
| `--samples 10 11 12` | Analyses only selected EXPNOs                 |
| `--exclude 5.5 6.5`  | Excludes an additional region of the spectrum |
| `--pluronic yes/no`  | Controls Pluronic F-68 modelling              |
| `--no-overlays`      | Does not create overlay PDFs                  |
| `--workers 4`        | Processes several spectra simultaneously      |

If you are new to MANC-Q, you can ignore these options initially.

---

# Common problems

## "python is not recognised"

Python is either not installed or was not added to Windows PATH.

Reinstall Python and make sure you select:

**Add Python to PATH**

during installation.

---

## "pip is not recognised"

Instead of using:

```text
pip install ...
```

use:

```text
python -m pip install ...
```

---

## The program cannot find my NMR data

Make sure you select the **Bruker experiment folder**, rather than an individual file.

MANC-Q expects the normal Bruker folder structure containing the EXPNO folders.

---

## My results don't look correct

First check:

1. The NMR spectrum has been correctly phased.
2. The baseline has been corrected.
3. The TSP/DSS concentration is correct.
4. The correct experiment folder was selected.
5. The overlay PDF for the sample.

Remember that some metabolites naturally have less reliable results because their NMR signals overlap with other compounds.

---

# Important limitations

MANC-Q is designed primarily for **cell-culture media and similar aqueous samples**.

Results can be affected by:

* Overlapping metabolite signals
* Changes in chemical shifts
* pH
* Matrix effects
* Incomplete T1 relaxation
* Errors in the internal-standard concentration
* Binding of the reference compound to proteins
* Differences between spectrometer frequencies

Some metabolites cannot always be separated reliably when they have very similar or identical signals.

Always consider the confidence tier and inspect the overlay when interpreting an unexpected result.

---

# I don't know anything about Python. Can I still use MANC-Q?

**Yes.**

You do not need to learn Python programming.

Python is simply the software environment that MANC-Q uses behind the scenes.

Once MANC-Q has been installed, you can set up a Windows shortcut or launcher so that the process is much easier.

If you are setting MANC-Q up for routine laboratory use, it is recommended to create a simple launcher specifically for your laboratory's folder structure.

---

# Getting help

If something does not work:

1. Take a screenshot of the error message.
2. Note which step you were doing.
3. Note which version of Windows and Python you are using.
4. If possible, include the `run_log.txt` file produced by MANC-Q.

The MANC-Q project is available here:

https://github.com/willwas1/MANC-Q

---

# Citation

If you use MANC-Q in research, please follow the citation information provided in the repository's `CITATION.cff` file.

MANC-Q was developed at the University of Manchester.

---

# Quick reference

### Install

```text
python -m pip install -r requirements.txt
```

### Test

```text
python -m mancq demo
```

### Analyse NMR data

```text
python -m mancq run "NMR_DATA_FOLDER" "RESULTS_FOLDER" --ref-mm 0.5
```

### Main result

```text
metabolite_concentrations.xlsx
```

### Remember

**Always check the confidence tier alongside the concentration.**
