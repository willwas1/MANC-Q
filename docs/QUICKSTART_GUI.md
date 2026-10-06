# MANC-Q quick start (graphical interface)

You need: a folder of Bruker 1D ¹H experiments (processed in TopSpin, or raw FIDs), and the TSP (or DSS)
concentration in the NMR tube. A run takes 1 to 5 minutes per spectrum.

**Start MANC-Q** from the Start menu (if you used the installer), or double-click MANC-Q.bat in the MANC-Q folder (or its desktop shortcut).

## 1 Spectra

![Spectra](images/1_spectra.png)

Click **Browse...** and choose the TopSpin dataset folder (the one that contains the numbered experiment
folders). Every experiment is listed. Untick any you do not want. The **Check** column flags spectra with a
weak or broad reference peak; click a row to see the spectrum.

First time? **Try the demo** makes two synthetic spectra (one processed, one raw FID) with known concentrations to practise on.

## 2 Process FIDs (only if some experiments are raw FIDs)

![Process FIDs](images/2_process.png)

Each FID is processed and phased automatically. For each one, look at the zoomed view: the baseline either
side of the peaks should be flat, with no dips below the dashed line. If it is not, move the **ph0** slider
until it is. Tick **I have checked the phasing**. Phasing errors change concentrations, so do not skip this.

## 3 Regions

![Regions](images/3_regions.png)

Regions left out of the fit. Water is ignored by default. If your medium contains poloxamer (Pluronic F-68),
leave it on **Model it**. To ignore anything else, drag across the spectrum or use **Add a common region**.
The purple box tells you which metabolites lose signal in the selected region.

## 4 Settings

![Settings](images/4_settings.png)

Enter the **reference concentration in the tube** (the calculator converts % w/v), the **dilution** from your
original sample to the tube (or "No dilution" for tube concentrations), and check the results folder. Click
**Start**.

## 5 Run

![Run](images/5_run.png)

Progress for each spectrum and the time remaining. You can cancel; finished spectra are kept.

## 6 Results

![Results](images/6_results.png)

Each value is coloured by how far it can be trusted:

| Colour | Meaning | Use |
|---|---|---|
| Blue | Quantified | a measurement |
| Orange | Overlapped | semi-quantitative; good for trends |
| Grey (`~`) | Deconvolution estimate | an estimate only |
| Pale grey (`<=`) | Upper bound | the most it can be; not a measurement |
| `<LOD` / `n/m` | not detected / not measurable (region ignored) | |

Click any value to see the fit it came from. **Open Excel workbook** gives every table; **Export for
GraphPad** writes samples x metabolites with only the values you chose to show.
