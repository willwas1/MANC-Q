# MANC-Q change log

Each version is kept in its own folder. Every change, however small, gets a new version number.

## 1.1.3 (7 October 2026)

Repository tidy-up only; the fitting and the window are unchanged.

* Added `CLAUDE.md` (notes for Claude Code: rules, layout, how to check a change).
* `README_windows.md` removed from GitHub (the 1.1.2 folder never had it, so the 1.1.2 commit already drops it).
* Version number raised to 1.1.3 in `pyproject.toml`, `mancq/engine.py`, `CITATION.cff`,
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
