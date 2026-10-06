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

**1–5**
