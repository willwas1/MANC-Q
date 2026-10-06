"""
MANC-Q desktop interface.

Six steps, one tab each: Spectra -> Process FIDs (only if some spectra are raw FIDs) -> Regions -> Settings ->
Run -> Results. The interface is a thin layer over the MANC-Q engine (package mancq): everything it does can also be done from the
command line, and both give identical numbers.
"""
import os
import sys
import json
import time
import datetime
import subprocess
import traceback

import numpy as np
import pandas as pd
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.widgets import SpanSelector

from PySide6.QtCore import Qt, QThread, Signal, QTimer
from PySide6.QtGui import QFont, QColor, QPalette, QBrush, QAction
from PySide6.QtWidgets import (QApplication, QMainWindow, QTabWidget, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
                               QLabel, QPushButton, QLineEdit, QTableWidget, QTableWidgetItem, QGroupBox, QComboBox,
                               QDoubleSpinBox, QRadioButton, QCheckBox, QProgressBar, QPlainTextEdit, QSlider,
                               QHeaderView, QSplitter, QSpinBox, QStatusBar, QFileDialog, QMessageBox, QButtonGroup,
                               QAbstractItemView, QSizePolicy)

from mancq import engine as q
from mancq import fidproc as fp
from mancq import results as rs


ACCENT = "#660099"
TIER_COL = {"Quantified": "#1d4ed8", "Overlapped (semi-quantitative)": "#d97706",
            "Deconvolution estimate (low confidence)": "#9ca3af", "Upper bound only": "#e5e7eb",
            "Not detected": "#ffffff", "Not measurable (region ignored)": "#f3f4f6"}
TIER_TXT = {"Quantified": "#ffffff", "Overlapped (semi-quantitative)": "#ffffff",
            "Deconvolution estimate (low confidence)": "#111827", "Upper bound only": "#6b7280",
            "Not detected": "#9ca3af", "Not measurable (region ignored)": "#9ca3af"}
MW = {"TSP-d4": 172.27, "DSS-d6": 224.36}
COMMON_REGIONS = [("Urea", 5.70, 5.85), ("DMSO", 2.69, 2.74), ("Methanol", 3.34, 3.38), ("Acetone", 2.21, 2.24)]
POLOX_MODES = ["Model it (recommended)", "Ignore", "Not in my samples"]
CACHE = os.path.join(os.path.expanduser("~"), ".mancq_cache")


# ============================================================================ helpers
def open_path(path):
    if sys.platform.startswith("win"):
        os.startfile(path)                                  # noqa: windows only
    elif sys.platform == "darwin":
        subprocess.Popen(["open", path])
    else:
        subprocess.Popen(["xdg-open", path])


def heading(text, sub=None):
    box = QWidget(); v = QVBoxLayout(box); v.setContentsMargins(0, 0, 0, 6)
    t = QLabel(text); t.setStyleSheet("font-size: 16px; font-weight: 600;"); v.addWidget(t)
    if sub:
        s = QLabel(sub); s.setStyleSheet("color: #4b5563;"); s.setWordWrap(True); v.addWidget(s)
    return box


def primary(text):
    b = QPushButton(text)
    b.setStyleSheet(f"QPushButton {{ background:{ACCENT}; color:white; padding:6px 18px; font-weight:600; }}"
                    "QPushButton:disabled { background:#c4b5d4; }")
    return b


def warn(parent, text, title="MANC-Q"):
    QMessageBox.warning(parent, title, text)


class Canvas(FigureCanvasQTAgg):
    def __init__(self, w=6, h=4):
        self.fig = Figure(figsize=(w, h), dpi=100)
        self.fig.patch.set_facecolor("white")
        super().__init__(self.fig)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)


def tidy(ax):
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.set_yticks([])


def scan_folder(folder, procno):
    """Every EXPNO folder with a processed spectrum and/or an FID."""
    out = []
    if not folder or not os.path.isdir(folder):
        return out
    cands = [folder] if os.path.exists(os.path.join(folder, "acqus")) else \
        [os.path.join(folder, d) for d in os.listdir(folder)]
    for d in sorted([c for c in cands if os.path.isdir(c)], key=lambda p: (len(os.path.basename(p)), os.path.basename(p))):
        proc, fid = fp.has_processed(d, procno), fp.has_fid(d)
        if not (proc or fid):
            continue
        title, field = "", None
        for tf in (os.path.join(d, "pdata", str(procno), "title"), os.path.join(d, "pdata", "1", "title")):
            if os.path.exists(tf):
                title = " ".join(open(tf, errors="ignore").read().split())[:60]
                break
        try:
            if os.path.exists(os.path.join(d, "acqus")):
                for ln in open(os.path.join(d, "acqus"), errors="ignore"):
                    if ln.startswith("##$BF1="):
                        field = float(ln.split("=")[1])
        except Exception:
            pass
        out.append(dict(expno=os.path.basename(d), path=d, title=title, field=field, processed=proc, fid=fid,
                        selected=True, snr=None, fwhm=None, check="", fid_lb=0.3, fid_ph0=None, fid_ph1=0.0,
                        fid_checked=False))
    return out


# ============================================================================ background workers
class InfoWorker(QThread):
    """Reference signal-to-noise and linewidth for each processed spectrum."""
    row_done = Signal(int, float, float, str)

    def __init__(self, spectra, procno):
        super().__init__()
        self.spectra, self.procno = spectra, procno

    def run(self):
        for i, s in enumerate(self.spectra):
            if not s["processed"]:
                continue
            try:
                ppm, y, sf, _ = q.read_spectrum(s["path"], self.procno)
                t = q.fit_tsp(ppm, y, sf)
                snr = t["height"] / q.noise_sigma(ppm, y)
                check = "OK"
                if snr < 200:
                    check = "Weak or missing reference peak"
                elif t["fwhm_hz"] > 2.0:
                    check = f"Broad reference line ({t['fwhm_hz']:.1f} Hz)"
                self.row_done.emit(i, float(snr), float(t["fwhm_hz"]), check)
            except Exception as exc:
                self.row_done.emit(i, float("nan"), float("nan"), f"Cannot read: {exc}")


class LibraryWorker(QThread):
    """Builds the compound library at a given field (cached on disk) for the Regions tab."""
    ready = Signal(object)

    def __init__(self, sf):
        super().__init__()
        self.sf = sf

    def run(self):
        try:
            lib = q.prepare_library(q.build_library(self.sf, os.path.join(q.HERE, "library"), CACHE), self.sf, False)
            self.ready.emit(lib)
        except Exception:
            self.ready.emit(None)


class RunWorker(QThread):
    """Starts the run in a child process and relays its progress messages to the window."""
    event = Signal(str, str, object)
    failed = Signal(str)

    def __init__(self, job):
        super().__init__()
        self.job = job
        import multiprocessing
        self.ctx = multiprocessing.get_context("spawn")
        self.events = self.ctx.Queue()
        self.stop_flag = self.ctx.Event()
        self.proc = None

    def stop(self):
        self.stop_flag.set()

    def kill(self):
        """Stop at once (used when the window is closed during a run); finished spectra are already saved."""
        self.stop_flag.set()
        if self.proc is not None and self.proc.is_alive():
            self.proc.terminate()

    def run(self):
        import queue as _queue
        from .runner import run_job
        self.proc = self.ctx.Process(target=run_job, args=(self.job, self.events, self.stop_flag), daemon=False)
        self.proc.start()
        while True:
            try:
                ev, smp, info = self.events.get(timeout=0.5)
            except _queue.Empty:
                if not self.proc.is_alive():
                    if self.proc.exitcode not in (0, None):
                        self.failed.emit(f"The calculation process stopped unexpectedly (exit code {self.proc.exitcode}).")
                    break
                continue
            if ev == "__end__":
                break
            if ev == "__error__":
                self.failed.emit(info["msg"]); continue
            self.event.emit(ev, smp, info)
        self.proc.join(timeout=10)


# ============================================================================ tab 1: spectra
class SpectraTab(QWidget):
    def __init__(self, win):
        super().__init__()
        self.win = win
        v = QVBoxLayout(self)
        v.addWidget(heading("1  Choose spectra", "Pick the TopSpin dataset folder (the one containing the numbered "
                                                  "experiment folders). Untick any spectra you do not want."))
        row = QHBoxLayout()
        self.path = QLineEdit(); self.path.setPlaceholderText("Folder of Bruker experiments, e.g. ...\\my_dataset")
        row.addWidget(self.path, 1)
        b = QPushButton("Browse..."); b.clicked.connect(self.browse); row.addWidget(b)
        row.addWidget(QLabel("PROCNO")); self.procno = QSpinBox(); self.procno.setRange(1, 999); self.procno.setValue(1)
        row.addWidget(self.procno)
        b2 = QPushButton("Scan"); b2.clicked.connect(self.scan); row.addWidget(b2)
        b3 = QPushButton("Try the demo"); b3.clicked.connect(self.demo); row.addWidget(b3)
        v.addLayout(row)
        split = QSplitter(Qt.Horizontal)
        self.table = QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(["", "EXPNO", "Title", "Field (MHz)", "Data", "Ref. S/N", "Check"])
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemChanged.connect(self.tick_changed)
        self.table.currentCellChanged.connect(lambda r, *a: self.preview(r))
        split.addWidget(self.table)
        self.canvas = Canvas(6.2, 3.8)
        split.addWidget(self.canvas); split.setSizes([560, 540])
        v.addWidget(split, 1)
        self.status = QLabel("No folder chosen yet."); v.addWidget(self.status)
        h = QHBoxLayout(); h.addStretch(); self.next = primary("Next"); self.next.clicked.connect(win.next_tab)
        h.addWidget(self.next); v.addLayout(h)
        self.info_worker = None

    def browse(self):
        d = QFileDialog.getExistingDirectory(self, "Choose the TopSpin dataset folder", self.path.text() or "")
        if d:
            self.path.setText(d); self.scan()

    def demo(self):
        from mancq import demo as dm
        folder = os.path.join(os.path.expanduser("~"), "MANC-Q demo")
        spec = os.path.join(folder, "spectra")
        if not (os.path.exists(os.path.join(spec, "10", "pdata", "1", "1r"))
                and os.path.exists(os.path.join(spec, "11", "fid"))):
            QApplication.setOverrideCursor(Qt.WaitCursor)
            try:
                dm.write_demo(spec, 600.13, cache_dir=CACHE)
            finally:
                QApplication.restoreOverrideCursor()
        self.path.setText(spec); self.scan()
        self.win.settings_tab.ref_mm.setValue(dm.REF_MM)
        self.win.settings_tab.r_none.setChecked(True)
        QMessageBox.information(self, "Demo spectrum", "Two synthetic 600 MHz spectra of 18 metabolites at known "
                                "concentrations have been made (reference 0.5 mM, no dilution): EXPNO 10 is processed, "
                                "EXPNO 11 is a raw FID so you can practise phasing. Click Next through the steps and "
                                "Start; the true concentrations are listed in mancq/demo.py.")

    def scan(self):
        folder = self.path.text().strip()
        spectra = scan_folder(folder, self.procno.value())
        self.win.spectra = spectra
        self.win.folder = folder
        self.fill()
        if not spectra:
            self.status.setText("No Bruker experiments found in this folder (expected EXPNO folders with "
                                "pdata/<procno>/1r or fid + acqus).")
            return
        n_fid = sum(1 for s in spectra if not s["processed"])
        self.status.setText(f"{len(spectra)} spectra found: {len(spectra) - n_fid} processed in TopSpin, "
                            f"{n_fid} raw FID only.")
        self.info_worker = InfoWorker(spectra, self.procno.value())
        self.info_worker.row_done.connect(self.info_row)
        self.info_worker.start()
        self.win.spectra_changed()
        if spectra:
            self.table.selectRow(0); self.preview(0)

    def fill(self):
        self.table.blockSignals(True)
        sp = self.win.spectra
        self.table.setRowCount(len(sp))
        for i, s in enumerate(sp):
            cb = QTableWidgetItem(); cb.setFlags(cb.flags() | Qt.ItemIsUserCheckable)
            cb.setCheckState(Qt.Checked if s["selected"] else Qt.Unchecked)
            self.table.setItem(i, 0, cb)
            data = "Processed" if s["processed"] else "FID only"
            vals = [s["expno"], s["title"], f"{s['field']:.1f}" if s["field"] else "", data,
                    "" if s["snr"] is None else f"{s['snr']:,.0f}", s["check"] or ("..." if s["processed"] else
                                                                                   "Processed in step 2")]
            for j, val in enumerate(vals, 1):
                it = QTableWidgetItem(val)
                if j == 4 and val == "FID only":
                    it.setForeground(QBrush(QColor("#b45309")))
                self.table.setItem(i, j, it)
        self.table.blockSignals(False)

    def info_row(self, i, snr, fwhm, check):
        sp = self.win.spectra
        if i >= len(sp):
            return
        sp[i].update(snr=snr, fwhm=fwhm, check=check)
        self.table.blockSignals(True)
        self.table.setItem(i, 5, QTableWidgetItem("" if not np.isfinite(snr) else f"{snr:,.0f}"))
        it = QTableWidgetItem(check)
        it.setForeground(QBrush(QColor("#15803d" if check == "OK" else "#b45309")))
        self.table.setItem(i, 6, it)
        self.table.blockSignals(False)

    def tick_changed(self, item):
        if item.column() == 0:
            self.win.spectra[item.row()]["selected"] = item.checkState() == Qt.Checked
            self.win.spectra_changed()

    def preview(self, row):
        sp = self.win.spectra
        self.canvas.fig.clear()
        ax = self.canvas.fig.add_subplot(111)
        if 0 <= row < len(sp):
            s = sp[row]
            if s["processed"]:
                try:
                    ppm, y, sf, _ = q.read_spectrum(s["path"], self.procno.value())
                    ax.plot(ppm, y, "k", lw=0.5)
                    ax.set_xlim(min(10, ppm.max()), max(-0.5, ppm.min()))
                    ax.set_ylim(-0.005 * y.max(), 0.08 * y.max())
                    ax.set_title(f"EXPNO {s['expno']}  {s['title']}", fontsize=9, loc="left")
                except Exception as exc:
                    ax.text(0.5, 0.5, f"Cannot read this spectrum:\n{exc}", ha="center", transform=ax.transAxes)
            else:
                ax.text(0.5, 0.5, "Raw FID only.\nIt is processed and checked in step 2.", ha="center",
                        va="center", transform=ax.transAxes, color="#6b7280")
        ax.set_xlabel("ppm"); tidy(ax)
        self.canvas.fig.tight_layout(); self.canvas.draw()


# ============================================================================ tab 2: process FIDs
class ProcessTab(QWidget):
    def __init__(self, win):
        super().__init__()
        self.win = win
        self.prep, self.cur = {}, None
        v = QVBoxLayout(self)
        v.addWidget(heading("2  Process raw FIDs", "Spectra without TopSpin processing are processed here: digital-"
                            "filter correction, line broadening, zero filling, Fourier transform, automatic phasing "
                            "and referencing to TSP. Check each one; a phasing error goes straight into "
                            "every concentration."))
        split = QSplitter(Qt.Horizontal)
        left = QWidget(); lv = QVBoxLayout(left); lv.setContentsMargins(0, 0, 0, 0)
        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["EXPNO", "Title", "Phase (ph0, ph1)", "Checked"])
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.currentCellChanged.connect(lambda r, *a: self.select(r))
        lv.addWidget(self.table)
        g = QGroupBox("Phase"); gl = QGridLayout(g)
        self.s0 = QSlider(Qt.Horizontal); self.s0.setRange(0, 3600)
        self.s1 = QSlider(Qt.Horizontal); self.s1.setRange(-900, 900)
        self.l0, self.l1 = QLabel(), QLabel()
        gl.addWidget(QLabel("Zero order (ph0)"), 0, 0); gl.addWidget(self.s0, 0, 1); gl.addWidget(self.l0, 0, 2)
        gl.addWidget(QLabel("First order (ph1)"), 1, 0); gl.addWidget(self.s1, 1, 1); gl.addWidget(self.l1, 1, 2)
        for sl in (self.s0, self.s1):
            sl.valueChanged.connect(self.phase_moved)
            sl.sliderReleased.connect(self.redraw_full)
        b = QPushButton("Auto-phase this spectrum again"); b.clicked.connect(self.auto_one); gl.addWidget(b, 2, 1)
        lv.addWidget(g)
        g2 = QGroupBox("Processing"); g2l = QGridLayout(g2)
        g2l.addWidget(QLabel("Line broadening"), 0, 0)
        self.lb = QDoubleSpinBox(); self.lb.setRange(0, 5); self.lb.setSingleStep(0.1); self.lb.setValue(0.3)
        self.lb.setSuffix(" Hz"); self.lb.editingFinished.connect(self.lb_changed); g2l.addWidget(self.lb, 0, 1)
        g2l.addWidget(QLabel("Baseline"), 1, 0); g2l.addWidget(QLabel("none needed: MANC-Q fits a smooth baseline itself"), 1, 1)
        g2l.addWidget(QLabel("Reference"), 2, 0); g2l.addWidget(QLabel("TSP/DSS set to 0.000 ppm automatically"), 2, 1)
        lv.addWidget(g2)
        b2 = QPushButton("Use this line broadening for all FIDs"); b2.clicked.connect(self.lb_all); lv.addWidget(b2)
        self.ok = QCheckBox("I have checked the phasing of this spectrum"); self.ok.toggled.connect(self.checked)
        lv.addWidget(self.ok)
        self.zoom = QComboBox()
        for lab in ["Zoom 0.8-1.6 ppm (methyl region)", "Zoom 3.0-4.2 ppm (sugar region)",
                    "Zoom 6.8-8.5 ppm (aromatic region)", "Zoom -0.1-0.1 ppm (reference)"]:
            self.zoom.addItem(lab)
        self.zoom.currentIndexChanged.connect(self.redraw_full)
        lv.addWidget(self.zoom)
        lv.addStretch()
        split.addWidget(left)
        self.canvas = Canvas(6.2, 4.6); split.addWidget(self.canvas); split.setSizes([440, 660])
        v.addWidget(split, 1)
        h = QHBoxLayout(); bk = QPushButton("Back"); bk.clicked.connect(win.prev_tab); h.addWidget(bk); h.addStretch()
        nx = primary("Next"); nx.clicked.connect(self.go_next); h.addWidget(nx); v.addLayout(h)

    def fids(self):
        return [s for s in self.win.spectra if s["selected"] and not s["processed"]]

    def refresh(self):
        f = self.fids()
        self.table.setRowCount(len(f))
        for i, s in enumerate(f):
            ph = "auto (not yet run)" if s["fid_ph0"] is None else f"{s['fid_ph0']:.1f}, {s['fid_ph1']:.1f}"
            for j, val in enumerate([s["expno"], s["title"], ph, "yes" if s["fid_checked"] else "no"]):
                self.table.setItem(i, j, QTableWidgetItem(val))
        if f and self.cur is None:
            self.table.selectRow(0); self.select(0)

    def _prep(self, s):
        key = (s["path"], s["fid_lb"])
        if key not in self.prep:
            QApplication.setOverrideCursor(Qt.WaitCursor)
            try:
                self.prep[key] = fp.prepare(s["path"], lb=s["fid_lb"])
            finally:
                QApplication.restoreOverrideCursor()
        return self.prep[key]

    def select(self, row):
        f = self.fids()
        if not (0 <= row < len(f)):
            return
        s = self.cur = f[row]
        try:
            p = self._prep(s)
        except Exception as exc:
            warn(self, f"Could not read the FID of EXPNO {s['expno']}:\n{exc}"); return
        if s["fid_ph0"] is None:
            QApplication.setOverrideCursor(Qt.WaitCursor)
            try:
                s["fid_ph0"], s["fid_ph1"] = fp.autophase(p["F"], p["keep"])
            finally:
                QApplication.restoreOverrideCursor()
        for sl, val in ((self.s0, s["fid_ph0"] * 10), (self.s1, s["fid_ph1"] * 10)):
            sl.blockSignals(True); sl.setValue(int(round(val))); sl.blockSignals(False)
        self.lb.setValue(s["fid_lb"])
        self.ok.blockSignals(True); self.ok.setChecked(s["fid_checked"]); self.ok.blockSignals(False)
        self.update_labels(); self.redraw_full(); self.refresh_row()

    def update_labels(self):
        self.l0.setText(f"{self.s0.value() / 10:.1f} deg"); self.l1.setText(f"{self.s1.value() / 10:.1f} deg")

    def phase_moved(self):
        if self.cur is None:
            return
        self.cur["fid_ph0"], self.cur["fid_ph1"] = self.s0.value() / 10, self.s1.value() / 10
        self.cur["fid_checked"] = False
        self.ok.blockSignals(True); self.ok.setChecked(False); self.ok.blockSignals(False)
        self.update_labels()
        if not (self.s0.isSliderDown() or self.s1.isSliderDown()):
            self.redraw_full()
        else:
            self.draw(baseline=False)
        self.refresh_row()

    def redraw_full(self):
        self.draw(baseline=False)  # same processing the run uses (no baseline correction)

    def draw(self, baseline=False):
        s = self.cur
        if s is None:
            return
        p = self._prep(s)
        r = fp.finish(p, s["fid_ph0"], s["fid_ph1"], baseline=baseline)
        zl = [(1.6, 0.8), (4.2, 3.0), (8.5, 6.8), (0.1, -0.1)][self.zoom.currentIndex()]
        self.canvas.fig.clear()
        ax1 = self.canvas.fig.add_subplot(211); ax2 = self.canvas.fig.add_subplot(212)
        for ax, (lo, hi), top, lab in [(ax1, (10, -0.5), 0.08, f"EXPNO {s['expno']}: whole spectrum"),
                                       (ax2, zl, 1.05, "Zoom: the baseline either side of each peak should be flat")]:
            m = (r["ppm"] < lo) & (r["ppm"] > hi)
            if m.sum() < 2:
                continue
            ax.plot(r["ppm"][m], r["y"][m], "k", lw=0.5)
            ax.axhline(0, color=ACCENT, lw=0.8, ls="--")
            mx = r["y"][m].max()
            ax.set_xlim(lo, hi); ax.set_ylim(-0.05 * mx, top * mx); tidy(ax)
            ax.set_title(lab, fontsize=9, loc="left")
        ax2.set_xlabel("ppm")
        self.canvas.fig.tight_layout(); self.canvas.draw()

    def auto_one(self):
        if self.cur is None:
            return
        p = self._prep(self.cur)
        self.cur["fid_ph0"], self.cur["fid_ph1"] = fp.autophase(p["F"], p["keep"])
        self.cur["fid_checked"] = False
        self.select(self.fids().index(self.cur))

    def lb_changed(self):
        if self.cur is not None and abs(self.cur["fid_lb"] - self.lb.value()) > 1e-9:
            self.cur["fid_lb"] = self.lb.value(); self.redraw_full()

    def lb_all(self):
        for s in self.fids():
            s["fid_lb"] = self.lb.value()
        self.redraw_full()

    def checked(self, on):
        if self.cur is not None:
            self.cur["fid_checked"] = on; self.refresh_row()

    def refresh_row(self):
        f = self.fids()
        if self.cur in f:
            i = f.index(self.cur); s = self.cur
            self.table.setItem(i, 2, QTableWidgetItem(f"{s['fid_ph0']:.1f}, {s['fid_ph1']:.1f}"))
            it = QTableWidgetItem("yes" if s["fid_checked"] else "no")
            it.setForeground(QBrush(QColor("#15803d" if s["fid_checked"] else "#b45309")))
            self.table.setItem(i, 3, it)

    def go_next(self):
        if any(not s["fid_checked"] for s in self.fids()):
            if QMessageBox.question(self, "Phasing not checked",
                                    "Not every FID has been marked as checked. Continue anyway?") != QMessageBox.Yes:
                return
        self.win.next_tab()


# ============================================================================ tab 3: regions
class RegionsTab(QWidget):
    def __init__(self, win):
        super().__init__()
        self.win = win
        self.lib, self.lib_worker = None, None
        self.regions = [dict(on=True, name="Water", lo=4.60, hi=5.00, kind="ignore"),
                        dict(on=True, name="Poloxamer PEO line", lo=3.69, hi=3.74, kind="polox"),
                        dict(on=True, name="Poloxamer PPO methyl", lo=1.10, hi=1.20, kind="polox")]
        self.polox_mode = POLOX_MODES[0]
        v = QVBoxLayout(self)
        v.addWidget(heading("3  Regions to ignore", "Ignored regions are left out of the fit in every spectrum. Water "
                            "is ignored by default. Drag across the spectrum to add your own region."))
        split = QSplitter(Qt.Horizontal)
        left = QWidget(); lv = QVBoxLayout(left); lv.setContentsMargins(0, 0, 0, 0)
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["", "Region", "From (ppm)", "To (ppm)", "How"])
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.table.itemChanged.connect(self.edited)
        self.table.currentCellChanged.connect(lambda r, *a: self.explain(r))
        lv.addWidget(self.table)
        hb = QHBoxLayout()
        b = QPushButton("Add region"); b.clicked.connect(lambda: self.add("My region", 6.00, 6.10)); hb.addWidget(b)
        self.pre = QComboBox(); self.pre.addItem("Add a common region...")
        for n, a, bb in COMMON_REGIONS:
            self.pre.addItem(f"{n} ({a:.2f}-{bb:.2f} ppm)")
        self.pre.activated.connect(self.add_common); hb.addWidget(self.pre)
        b2 = QPushButton("Remove"); b2.clicked.connect(self.remove); hb.addWidget(b2); hb.addStretch()
        lv.addLayout(hb)
        self.info = QLabel(); self.info.setWordWrap(True); self.info.setTextFormat(Qt.RichText)
        self.info.setStyleSheet("background:#f5f3ff; border:1px solid #ddd6fe; padding:8px;")
        self.info.setAlignment(Qt.AlignTop); self.info.setMinimumHeight(150)
        lv.addWidget(self.info)
        note = QLabel("Metabolites whose peaks all fall inside ignored regions are reported as \"n/m\" (not "
                      "measurable), not as zero.")
        note.setWordWrap(True); note.setStyleSheet("color:#4b5563;"); lv.addWidget(note)
        lv.addStretch()
        split.addWidget(left)
        self.canvas = Canvas(6.2, 4.4); split.addWidget(self.canvas); split.setSizes([520, 580])
        v.addWidget(split, 1)
        h = QHBoxLayout(); bk = QPushButton("Back"); bk.clicked.connect(win.prev_tab); h.addWidget(bk); h.addStretch()
        nx = primary("Next"); nx.clicked.connect(win.next_tab); h.addWidget(nx); v.addLayout(h)
        self.spec = None
        self.fill()

    # ---- table
    def fill(self):
        self.table.blockSignals(True)
        self.table.setRowCount(len(self.regions))
        for i, r in enumerate(self.regions):
            cb = QTableWidgetItem(); cb.setFlags(cb.flags() | Qt.ItemIsUserCheckable)
            cb.setCheckState(Qt.Checked if r["on"] else Qt.Unchecked); self.table.setItem(i, 0, cb)
            self.table.setItem(i, 1, QTableWidgetItem(r["name"]))
            self.table.setItem(i, 2, QTableWidgetItem(f"{r['lo']:.2f}"))
            self.table.setItem(i, 3, QTableWidgetItem(f"{r['hi']:.2f}"))
            combo = QComboBox()
            if r["kind"] == "polox":
                combo.addItems(POLOX_MODES); combo.setCurrentText(self.polox_mode)
                combo.currentTextChanged.connect(self.polox_changed)
            else:
                combo.addItems(["Ignore"])
            self.table.setCellWidget(i, 4, combo)
        self.table.blockSignals(False)
        self.draw()

    def edited(self, item):
        i, c = item.row(), item.column()
        r = self.regions[i]
        if c == 0:
            r["on"] = item.checkState() == Qt.Checked
        elif c == 1:
            r["name"] = item.text()
        elif c in (2, 3):
            try:
                val = float(item.text())
                r["lo" if c == 2 else "hi"] = val
                if r["lo"] > r["hi"]:
                    r["lo"], r["hi"] = r["hi"], r["lo"]
            except ValueError:
                pass
            self.fill(); return
        self.draw(); self.explain(i)

    def polox_changed(self, txt):
        self.polox_mode = txt
        self.fill()
        self.explain(self.table.currentRow())

    def add(self, name, lo, hi):
        self.regions.append(dict(on=True, name=name, lo=min(lo, hi), hi=max(lo, hi), kind="ignore"))
        self.fill(); self.table.selectRow(len(self.regions) - 1); self.explain(len(self.regions) - 1)

    def add_common(self, idx):
        if idx > 0:
            n, a, b = COMMON_REGIONS[idx - 1]
            self.add(n, a, b)
        self.pre.setCurrentIndex(0)

    def remove(self):
        i = self.table.currentRow()
        if 0 <= i < len(self.regions):
            if self.regions[i]["kind"] == "polox":
                warn(self, "The poloxamer rows cannot be removed. Set \"How\" to \"Not in my samples\" instead.")
                return
            self.regions.pop(i); self.fill()

    # ---- what each region costs
    def affected(self, lo, hi):
        if self.lib is None:
            return None
        out = []
        for c in self.lib:
            lost = sum(m["nprot"] for m in c["mults"] if m["hi"] > lo and m["lo"] < hi)
            if lost > 0:
                out.append((c["name"], lost / c["nprot"]))
        return sorted(out, key=lambda t: -t[1])

    def explain(self, i):
        if not (0 <= i < len(self.regions)):
            self.info.setText(""); return
        r = self.regions[i]
        txt = f"<b>{r['name']}, {r['lo']:.2f}-{r['hi']:.2f} ppm</b><br>"
        aff = self.affected(r["lo"], r["hi"])
        if r["kind"] == "polox":
            txt += ("<b>Model it</b> (recommended): the polymer is fitted, so metabolite peaks underneath are kept "
                    "(the centre of the tall PEO line is always left out).<br><b>Ignore</b>: both poloxamer "
                    "regions are left out of the fit; the metabolite peaks below are lost.<br><b>Not in my "
                    "samples</b>: nothing is fitted or removed.<br>")
        if aff is None:
            txt += "<i>Working out which metabolites have peaks here...</i>"
        elif not aff:
            txt += "No library metabolite has peaks in this region."
        else:
            whole = [n for n, f in aff if f >= 0.999]
            part = [f"{n} ({f:.0%} of its signal)" for n, f in aff if f < 0.999][:10]
            if whole:
                txt += f"<b>Not measurable if ignored:</b> {', '.join(whole)}<br>"
            if part:
                txt += f"<b>Lose part of their signal if ignored:</b> {', '.join(part)}"
                if len(aff) - len(whole) > 10:
                    txt += f" and {len(aff) - len(whole) - 10} more"
        self.info.setText(txt)

    # ---- spectrum with shaded regions
    def refresh(self):
        sel = [s for s in self.win.spectra if s["selected"] and s["processed"]]
        if sel and (self.spec is None or self.spec[0] != sel[0]["path"]):
            try:
                ppm, y, sf, _ = q.read_spectrum(sel[0]["path"], self.win.spectra_tab.procno.value())
                t = q.fit_tsp(ppm, y, sf)
                self.spec = (sel[0]["path"], ppm - t["centre"], y, sel[0]["expno"])
                self.start_library(sf)
            except Exception:
                self.spec = None
        elif not sel:
            fs = [s for s in self.win.spectra if s["selected"] and s["field"]]
            if fs:
                self.start_library(fs[0]["field"])
        self.draw(); self.explain(self.table.currentRow() if self.table.currentRow() >= 0 else 0)

    def start_library(self, sf):
        if self.lib is None and self.lib_worker is None:
            self.lib_worker = LibraryWorker(sf)
            self.lib_worker.ready.connect(self.library_ready)
            self.lib_worker.start()

    def library_ready(self, lib):
        self.lib = lib or []
        self.explain(max(0, self.table.currentRow()))

    def draw(self):
        self.canvas.fig.clear()
        ax = self.canvas.fig.add_subplot(111)
        if self.spec is not None:
            _, ppm, y, ex = self.spec
            m = (ppm < 10) & (ppm > -0.5)
            ax.plot(ppm[m], y[m], "k", lw=0.5)
            ax.set_ylim(-0.01 * y[m].max(), 0.08 * y[m].max())
            ax.set_title(f"EXPNO {ex}  (drag across the spectrum to add a region; scroll to zoom)", fontsize=9, loc="left")
        else:
            ax.set_title("No processed spectrum selected yet", fontsize=9, loc="left")
        ax.set_xlim(10, -0.5)
        for r in self.regions:
            if not r["on"]:
                continue
            if r["kind"] == "polox":
                if self.polox_mode == POLOX_MODES[2]:
                    continue
                col, a = (ACCENT, 0.18) if self.polox_mode == POLOX_MODES[0] else ("#9ca3af", 0.45)
            else:
                col, a = "#9ca3af", 0.45
            ax.axvspan(r["lo"], r["hi"], color=col, alpha=a, lw=0)
        ax.set_xlabel("ppm"); tidy(ax)
        self.canvas.fig.tight_layout()
        self.span = SpanSelector(ax, self.dragged, "horizontal", useblit=True, minspan=0.005,
                                 props=dict(alpha=0.25, facecolor=ACCENT))
        self.canvas.mpl_connect("scroll_event", self.scroll)
        self.canvas.draw()

    def dragged(self, a, b):
        if abs(b - a) >= 0.005:
            self.add("Custom region", round(min(a, b), 3), round(max(a, b), 3))

    def scroll(self, ev):
        ax = ev.inaxes
        if ax is None or ev.xdata is None:
            return
        f = 0.8 if ev.button == "up" else 1.25
        x0, x1 = ax.get_xlim()
        c = ev.xdata
        ax.set_xlim(c + (x0 - c) * f, c + (x1 - c) * f)
        if self.spec is not None:
            _, ppm, y, _ = self.spec
            lo, hi = sorted(ax.get_xlim())
            m = (ppm > lo) & (ppm < hi) & ~((ppm > 4.6) & (ppm < 5.0))
            if m.sum() > 2:
                ax.set_ylim(-0.02 * y[m].max(), 1.05 * y[m].max())
        self.canvas.draw_idle()

    # ---- what the engine needs
    def engine_settings(self):
        excl = [(r["lo"], r["hi"]) for r in self.regions if r["on"] and r["kind"] == "ignore"]
        polox_on = any(r["on"] for r in self.regions if r["kind"] == "polox")
        mode = {POLOX_MODES[0]: "auto", POLOX_MODES[1]: "ignore", POLOX_MODES[2]: "no"}[self.polox_mode]
        if not polox_on:
            mode = "no"
        pregions = [(r["lo"], r["hi"]) for r in self.regions if r["kind"] == "polox" and r["on"]]
        return dict(EXCLUDE=excl, PLURONIC=mode, PLURONIC_REGIONS=pregions or q.DEFAULT_SETTINGS["PLURONIC_REGIONS"])

    def to_json(self):
        return dict(regions=self.regions, polox_mode=self.polox_mode)

    def from_json(self, d):
        if "regions" in d:
            self.regions = d["regions"]
        self.polox_mode = d.get("polox_mode", self.polox_mode)
        self.fill()


# ============================================================================ tab 4: settings
class SettingsTab(QWidget):
    def __init__(self, win):
        super().__init__()
        self.win = win
        v = QVBoxLayout(self)
        v.addWidget(heading("4  Settings", "One number is required: the reference concentration in the NMR tube."))
        g = QGroupBox("Internal reference (required)"); gl = QGridLayout(g)
        gl.addWidget(QLabel("Compound"), 0, 0)
        self.ref = QComboBox(); self.ref.addItems(list(MW)); self.ref.currentTextChanged.connect(self.ref_changed)
        gl.addWidget(self.ref, 0, 1)
        gl.addWidget(QLabel("Concentration in the tube"), 1, 0)
        self.ref_mm = QDoubleSpinBox(); self.ref_mm.setDecimals(4); self.ref_mm.setRange(0, 1000)
        self.ref_mm.setSingleStep(0.01); self.ref_mm.setValue(0.0); self.ref_mm.setSuffix(" mM")
        gl.addWidget(self.ref_mm, 1, 1)
        calc = QGroupBox("Calculator: I know the % w/v in the tube"); cl = QHBoxLayout(calc)
        self.pct = QDoubleSpinBox(); self.pct.setDecimals(4); self.pct.setRange(0, 10); self.pct.setValue(0.01)
        self.pct.setSuffix(" % w/v"); cl.addWidget(self.pct)
        cl.addWidget(QLabel("molar mass")); self.mw = QDoubleSpinBox(); self.mw.setRange(1, 2000); self.mw.setDecimals(2)
        self.mw.setValue(MW["TSP-d4"]); self.mw.setSuffix(" g/mol"); cl.addWidget(self.mw)
        self.calc_out = QLabel(); cl.addWidget(self.calc_out)
        use = QPushButton("Use this value"); use.clicked.connect(lambda: self.ref_mm.setValue(self.calc_value()))
        cl.addWidget(use)
        for w_ in (self.pct, self.mw):
            w_.valueChanged.connect(self.calc_update)
        self.calc_update()
        gl.addWidget(calc, 2, 0, 1, 2)
        v.addWidget(g)
        g2 = QGroupBox("Dilution (original sample to NMR tube)"); g2l = QVBoxLayout(g2)
        self.r_none = QRadioButton("No dilution: report concentrations in the tube")
        self.r_same = QRadioButton("Same for every sample:")
        self.r_each = QRadioButton("Different per sample (edit the table)")
        grp = QButtonGroup(self)
        for r_ in (self.r_none, self.r_same, self.r_each):
            grp.addButton(r_)
        self.r_same.setChecked(True)
        g2l.addWidget(self.r_none)
        h = QHBoxLayout(); h.addWidget(self.r_same)
        self.dil = QDoubleSpinBox(); self.dil.setRange(0.001, 10000); self.dil.setDecimals(3); self.dil.setValue(1.0)
        self.dil.setPrefix("x "); h.addWidget(self.dil); h.addStretch(); g2l.addLayout(h)
        g2l.addWidget(self.r_each)
        self.dtable = QTableWidget(0, 3); self.dtable.setHorizontalHeaderLabels(["EXPNO", "Title", "Dilution"])
        self.dtable.verticalHeader().setVisible(False); self.dtable.setMaximumHeight(140)
        self.dtable.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        g2l.addWidget(self.dtable)
        self.r_each.toggled.connect(self.dtable.setEnabled); self.dtable.setEnabled(False)
        self.dil.valueChanged.connect(self.sync_dilutions)
        v.addWidget(g2)
        g3 = QGroupBox("Output and speed"); g3l = QGridLayout(g3)
        g3l.addWidget(QLabel("Results folder"), 0, 0)
        self.out = QLineEdit(); g3l.addWidget(self.out, 0, 1)
        bo = QPushButton("Browse..."); bo.clicked.connect(self.browse_out); g3l.addWidget(bo, 0, 2)
        g3l.addWidget(QLabel("Spectra fitted at the same time"), 1, 0)
        self.workers = QSpinBox(); self.workers.setRange(1, max(1, os.cpu_count() or 1))
        self.workers.setValue(max(1, (os.cpu_count() or 2) - 1)); g3l.addWidget(self.workers, 1, 1)
        self.overlays = QCheckBox("Write an overlay PDF for every spectrum (recommended)"); self.overlays.setChecked(True)
        g3l.addWidget(self.overlays, 2, 0, 1, 3)
        v.addWidget(g3)
        h2 = QHBoxLayout()
        bs = QPushButton("Save settings..."); bs.clicked.connect(self.save); h2.addWidget(bs)
        bl = QPushButton("Load settings..."); bl.clicked.connect(self.load); h2.addWidget(bl); h2.addStretch()
        v.addLayout(h2)
        v.addStretch()
        h = QHBoxLayout(); bk = QPushButton("Back"); bk.clicked.connect(win.prev_tab); h.addWidget(bk); h.addStretch()
        self.start = primary("Start"); self.start.clicked.connect(win.start_run); h.addWidget(self.start)
        v.addLayout(h)

    def ref_changed(self, name):
        self.mw.setValue(MW.get(name, self.mw.value()))

    def calc_value(self):
        return self.pct.value() * 10 / self.mw.value() * 1000

    def calc_update(self):
        self.calc_out.setText(f"=  {self.calc_value():.4f} mM")

    def browse_out(self):
        d = QFileDialog.getExistingDirectory(self, "Choose where to put the results", self.out.text() or "")
        if d:
            self.out.setText(d)

    def refresh(self):
        sel = [s for s in self.win.spectra if s["selected"]]
        old = {self.dtable.item(i, 0).text(): self.dtable.item(i, 2).text() for i in range(self.dtable.rowCount())
               if self.dtable.item(i, 0) and self.dtable.item(i, 2)}
        self.dtable.setRowCount(len(sel))
        for i, s in enumerate(sel):
            for j, val in enumerate([s["expno"], s["title"], old.get(s["expno"], f"{self.dil.value():g}")]):
                it = QTableWidgetItem(val)
                if j < 2:
                    it.setFlags(it.flags() & ~Qt.ItemIsEditable)
                self.dtable.setItem(i, j, it)
        if not self.out.text() and self.win.folder:
            base = os.path.join(os.path.expanduser("~"), "Documents", "MANC-Q results")
            name = os.path.basename(os.path.normpath(self.win.folder)) or "run"
            self.out.setText(os.path.join(base, f"{name}_{datetime.datetime.now():%Y-%m-%d_%H%M}"))

    def sync_dilutions(self):
        """While one factor applies to all, keep the per-sample table showing it (a starting point for editing)."""
        if not self.r_each.isChecked():
            for i in range(self.dtable.rowCount()):
                self.dtable.setItem(i, 2, QTableWidgetItem(f"{self.dil.value():g}"))

    def dilutions(self):
        """(single factor, per-sample dict or None)."""
        if self.r_none.isChecked():
            return 1.0, None
        if self.r_same.isChecked():
            return self.dil.value(), None
        per = {}
        for i in range(self.dtable.rowCount()):
            try:
                per[self.dtable.item(i, 0).text()] = float(self.dtable.item(i, 2).text())
            except (ValueError, AttributeError):
                raise ValueError(f"The dilution for EXPNO {self.dtable.item(i, 0).text()} is not a number.")
        return 1.0, per

    def to_json(self):
        return dict(reference=self.ref.currentText(), ref_mm=self.ref_mm.value(),
                    dilution_mode="none" if self.r_none.isChecked() else "same" if self.r_same.isChecked() else "each",
                    dilution=self.dil.value(), workers=self.workers.value(), overlays=self.overlays.isChecked())

    def from_json(self, d):
        self.ref.setCurrentText(d.get("reference", "TSP-d4")); self.ref_mm.setValue(d.get("ref_mm", 0.0))
        {"none": self.r_none, "same": self.r_same, "each": self.r_each}[d.get("dilution_mode", "same")].setChecked(True)
        self.dil.setValue(d.get("dilution", 1.0)); self.workers.setValue(d.get("workers", self.workers.value()))
        self.overlays.setChecked(d.get("overlays", True))

    def save(self):
        f, _ = QFileDialog.getSaveFileName(self, "Save settings", "MANC-Q settings.json", "Settings (*.json)")
        if f:
            json.dump(dict(settings=self.to_json(), regions=self.win.regions_tab.to_json()), open(f, "w"), indent=2)

    def load(self):
        f, _ = QFileDialog.getOpenFileName(self, "Load settings", "", "Settings (*.json)")
        if f:
            d = json.load(open(f))
            self.from_json(d.get("settings", {})); self.win.regions_tab.from_json(d.get("regions", {}))


# ============================================================================ tab 5: run
class RunTab(QWidget):
    def __init__(self, win):
        super().__init__()
        self.win = win
        self.worker, self.out, self.t_start, self.rows, self.log_pos = None, None, None, {}, 0
        self.durations, self.n_total, self.n_done = [], 0, 0
        v = QVBoxLayout(self)
        v.addWidget(heading("5  Run", "Each spectrum takes 1 to 5 minutes. You can leave this running. If it is "
                            "stopped, finished spectra are kept, and starting again with the same results folder "
                            "carries on from where it stopped."))
        self.overall = QLabel("Not started."); v.addWidget(self.overall)
        self.bar = QProgressBar(); self.bar.setValue(0); v.addWidget(self.bar)
        self.table = QTableWidget(0, 4); self.table.setHorizontalHeaderLabels(["EXPNO", "Title", "Status", "Result"])
        self.table.verticalHeader().setVisible(False); self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True); self.table.setMaximumHeight(230)
        v.addWidget(self.table)
        self.log = QPlainTextEdit(); self.log.setReadOnly(True)
        self.log.setStyleSheet("font-family: Consolas, 'DejaVu Sans Mono', monospace; font-size: 11px;")
        v.addWidget(self.log, 1)
        h = QHBoxLayout(); h.addStretch()
        self.cancel = QPushButton("Cancel after the spectra now running"); self.cancel.setEnabled(False)
        self.cancel.clicked.connect(self.do_cancel); h.addWidget(self.cancel)
        self.view = primary("View results"); self.view.setEnabled(False); self.view.clicked.connect(self.show_results)
        h.addWidget(self.view); v.addLayout(h)
        self.timer = QTimer(self); self.timer.timeout.connect(self.tick)

    def start(self, job):
        self.out = job["out"]
        self.table.setRowCount(len(job["spectra"])); self.rows = {}
        for i, s in enumerate(job["spectra"]):
            self.rows[s["expno"]] = i
            for j, val in enumerate([s["expno"], s["title"], "Waiting", ""]):
                self.table.setItem(i, j, QTableWidgetItem(val))
        self.n_total, self.n_done, self.durations = len(job["spectra"]), 0, []
        self.bar.setRange(0, self.n_total); self.bar.setValue(0)
        self.log.clear(); self.log_pos = 0; self.t_start = time.time()
        self.overall.setText("Starting...")
        self.cancel.setEnabled(True); self.view.setEnabled(False)
        self.worker = RunWorker(job)
        self.worker.event.connect(self.on_event)
        self.worker.failed.connect(self.on_failed)
        self.worker.finished.connect(self.on_thread_end)
        self.worker.start(); self.timer.start(1000)

    def set_status(self, smp, text, colour="#111827", result=None):
        i = self.rows.get(smp)
        if i is None:
            return
        it = QTableWidgetItem(text); it.setForeground(QBrush(QColor(colour))); self.table.setItem(i, 2, it)
        if result is not None:
            self.table.setItem(i, 3, QTableWidgetItem(result))

    def on_event(self, ev, smp, info):
        if ev == "fid_processing":
            self.set_status(smp, f"Processing FID ({info['k']} of {info['n']})", ACCENT)
        elif ev == "run_started":
            self.overall.setText(f"Fitting {info['n']} spectra, {info['workers']} at a time...")
        elif ev == "sample_started":
            self.set_status(smp, "Fitting...", ACCENT)
        elif ev == "sample_done":
            self.n_done += 1
            if not info.get("reused"):
                self.durations.append(info["seconds"])
            res = f"{info['n_quantified']} quantified; library explains {info['explained']:.0%} of the signal"
            self.set_status(smp, "Done (earlier result reused)" if info.get("reused") else
                            f"Done in {info['seconds'] / 60:.1f} min", "#15803d", res)
        elif ev == "sample_failed":
            self.n_done += 1
            self.set_status(smp, "Failed", "#b91c1c", info.get("error", ""))
        elif ev == "cancelled":
            for s in info.get("not_fitted", []):
                self.set_status(s, "Not fitted (cancelled)", "#6b7280")
        elif ev == "aggregating":
            self.overall.setText("Writing the results workbook and plots...")
        elif ev == "finished":
            if info.get("ok"):
                self.overall.setText(f"Finished. Results are in {info.get('out', self.out)}")
                self.view.setEnabled(True)
        self.bar.setValue(self.n_done)

    def tick(self):
        f = os.path.join(self.out or "", "run_log.txt")
        if os.path.exists(f):
            txt = open(f, errors="ignore").read()
            if len(txt) > self.log_pos:
                self.log.appendPlainText(txt[self.log_pos:].rstrip("\n")); self.log_pos = len(txt)
        if self.worker is not None and self.worker.isRunning() and self.n_done < self.n_total:
            el = time.time() - self.t_start
            msg = f"{self.n_done} of {self.n_total} spectra done, {el / 60:.0f} min elapsed"
            if self.durations:
                w = self.win.settings_tab.workers.value()
                left = (self.n_total - self.n_done) * np.mean(self.durations) / max(1, min(w, self.n_total))
                msg += f", about {max(1, round(left / 60))} min remaining"
            self.overall.setText(msg)

    def do_cancel(self):
        if self.worker is not None:
            self.worker.stop(); self.cancel.setEnabled(False)
            self.overall.setText("Cancelling: spectra already running will finish first...")

    def on_failed(self, msg):
        self.overall.setText("The run stopped with an error.")
        self.log.appendPlainText("\n" + msg)
        warn(self, "The run stopped with an error. The details are in the log on this page.")

    def on_thread_end(self):
        self.timer.stop(); self.tick(); self.cancel.setEnabled(False)
        self.win.running = False
        if rs.is_results_folder(self.out or ""):
            self.view.setEnabled(True)

    def show_results(self):
        self.win.results_tab.load(self.out)
        self.win.tabs.setCurrentWidget(self.win.results_tab)


# ============================================================================ tab 6: results
class ResultsTab(QWidget):
    FILTERS = [("All, including upper bounds", None),
               ("Quantified, Overlapped and Estimates", ["Quantified", "Overlapped (semi-quantitative)",
                                                         "Deconvolution estimate (low confidence)"]),
               ("Quantified and Overlapped", ["Quantified", "Overlapped (semi-quantitative)"]),
               ("Quantified only", ["Quantified"])]

    def __init__(self, win):
        super().__init__()
        self.win = win
        self.out, self.data = None, None
        v = QVBoxLayout(self)
        v.addWidget(heading("6  Results", "Colour shows how far each number can be trusted. Click a cell to see the "
                            "fit it came from."))
        bar = QHBoxLayout()
        bo = QPushButton("Open a results folder..."); bo.clicked.connect(self.choose); bar.addWidget(bo)
        bar.addWidget(QLabel("Show:"))
        self.filter = QComboBox(); self.filter.addItems([f[0] for f in self.FILTERS]); self.filter.setCurrentIndex(2)
        self.filter.currentIndexChanged.connect(self.fill); bar.addWidget(self.filter)
        bar.addSpacing(12)
        for lab, t in [("Quantified", "Quantified"), ("Overlapped", "Overlapped (semi-quantitative)"),
                       ("Estimate", "Deconvolution estimate (low confidence)"), ("Upper bound", "Upper bound only")]:
            sw = QLabel(f"  {lab}  ")
            sw.setStyleSheet(f"background:{TIER_COL[t]}; color:{TIER_TXT[t]}; border:1px solid #d1d5db; padding:2px;")
            bar.addWidget(sw)
        bar.addStretch(); v.addLayout(bar)
        self.where = QLabel(""); self.where.setStyleSheet("color:#4b5563;"); v.addWidget(self.where)
        split = QSplitter(Qt.Horizontal)
        self.table = QTableWidget(0, 0); self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.currentCellChanged.connect(lambda r, c, *a: self.cell(r, c))
        split.addWidget(self.table)
        self.canvas = Canvas(5.6, 4.4); split.addWidget(self.canvas); split.setSizes([560, 540])
        v.addWidget(split, 1)
        h = QHBoxLayout()
        for lab, fn in [("Open Excel workbook", self.open_xlsx), ("Open results folder", self.open_folder),
                        ("Export for GraphPad", self.export_graphpad), ("Open overlay PDF (this sample)", self.open_overlay)]:
            b = QPushButton(lab); b.clicked.connect(fn); h.addWidget(b)
        h.addStretch(); v.addLayout(h)

    def choose(self):
        d = QFileDialog.getExistingDirectory(self, "Choose a MANC-Q results folder", self.out or "")
        if d:
            if not rs.is_results_folder(d):
                warn(self, "That folder does not contain MANC-Q results (concentrations_long.csv)."); return
            self.load(d)

    def load(self, out):
        self.out = out
        self.data = rs.load(out)
        self.where.setText(f"Results folder: {out}")
        self.fill()

    def fill(self):
        if self.data is None:
            return
        G, samples = self.data["long"], self.data["samples"]
        allowed = self.FILTERS[self.filter.currentIndex()][1]
        titles = dict(zip(self.data["qc"]["sample"], self.data["qc"].get("title", pd.Series(dtype=str)).fillna("")))
        key = G.set_index(["metabolite", "sample"])
        mets = [m for m in self.data["metabolites"]
                if allowed is None or G[(G.metabolite == m) & G.tier.isin(allowed)].shape[0] > 0]
        if allowed is None:
            mets = [m for m in mets if (G[G.metabolite == m].tier != "Not detected").any()]
        self.mets, self.samples = mets, samples
        self.table.clear()
        self.table.setRowCount(len(mets)); self.table.setColumnCount(len(samples))
        self.table.setHorizontalHeaderLabels([f"{s}\n{str(titles.get(s, ''))[:14]}" for s in samples])
        self.table.setVerticalHeaderLabels(mets)
        for i, m in enumerate(mets):
            for j, s in enumerate(samples):
                if (m, s) not in key.index:
                    continue
                r = key.loc[(m, s)]
                t = r.tier
                if t == "Not measurable (region ignored)":
                    txt = "n/m"
                elif t == "Not detected":
                    txt = "<LOD"
                elif t == "Upper bound only":
                    txt = f"<={r.conc_mM:.2g}"
                elif t.startswith("Deconvolution"):
                    txt = f"~{r.conc_mM:.3g}"
                else:
                    txt = f"{r.conc_mM:.3g}"
                show = allowed is None or t in allowed
                it = QTableWidgetItem(txt if show else "")
                it.setTextAlignment(Qt.AlignCenter)
                if show:
                    it.setBackground(QBrush(QColor(TIER_COL.get(t, "#ffffff"))))
                    it.setForeground(QBrush(QColor(TIER_TXT.get(t, "#111827"))))
                it.setToolTip(f"{m}, EXPNO {s}: {t}, {r.method}")
                self.table.setItem(i, j, it)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        if mets and samples:
            self.table.setCurrentCell(0, 0)

    def cell(self, r, c):
        self.canvas.fig.clear()
        if self.data is None or r < 0 or c < 0 or r >= len(self.mets) or c >= len(self.samples):
            self.canvas.draw(); return
        m, s = self.mets[r], self.samples[c]
        G = self.data["long"]
        row = G[(G.metabolite == m) & (G["sample"] == s)]
        win_ = rs.metabolite_window(self.out, s, m)
        ax = self.canvas.fig.add_subplot(211); ax2 = self.canvas.fig.add_subplot(212, sharex=ax)
        if win_ is None or not len(row):
            ax.text(0.5, 0.5, "No fitted region to show for this metabolite.", ha="center", transform=ax.transAxes)
            tidy(ax); tidy(ax2); self.canvas.draw(); return
        x = win_["ppm"]
        ax.plot(x, win_["observed"], "k", lw=0.8, label="measured")
        ax.plot(x, win_["fit"], color="#dc2626", lw=0.8, label="fit")
        ax.fill_between(x, win_["baseline"], win_["baseline"] + win_["compound"], color=ACCENT, alpha=0.35, label=m)
        t = row.iloc[0]
        val = "n/m" if t.tier.startswith("Not measurable") else f"{t.conc_mM:.3g} mM"
        ax.set_title(f"{m} in EXPNO {s}: {val}, {t.tier}", fontsize=9, loc="left")
        ax.legend(fontsize=7, frameon=False); tidy(ax)
        res = win_["observed"] - win_["fit"]
        top = max(np.abs(win_["observed"]).max(), 1e-9)
        ax2.plot(x, res, color="grey", lw=0.7); ax2.axhline(0, color="k", lw=0.4)
        ax2.set_ylim(-0.3 * top, 0.3 * top); tidy(ax2); ax2.set_ylabel("residual", fontsize=8)
        ax2.set_xlim(x.max(), x.min()); ax2.set_xlabel("ppm")
        self.canvas.fig.tight_layout(); self.canvas.draw()

    def open_xlsx(self):
        if self.out:
            open_path(os.path.join(self.out, "metabolite_concentrations.xlsx"))

    def open_folder(self):
        if self.out:
            open_path(self.out)

    def open_overlay(self):
        c = self.table.currentColumn()
        if self.out and 0 <= c < len(self.samples):
            f = os.path.join(self.out, f"overlay_sample_{self.samples[c]}.pdf")
            if os.path.exists(f):
                open_path(f)
            else:
                warn(self, "No overlay PDF for this sample (overlays may have been switched off).")

    def export_graphpad(self):
        """Two tables, samples as rows and metabolites as columns: values (only tiers with a number) and tiers."""
        if self.data is None:
            return
        G = self.data["long"]
        allowed = self.FILTERS[self.filter.currentIndex()][1] or ["Quantified", "Overlapped (semi-quantitative)",
                                                                 "Deconvolution estimate (low confidence)"]
        val = G.assign(v=np.where(G.tier.isin(allowed), G.conc_mM, np.nan)).pivot(index="sample", columns="metabolite", values="v")
        tier = G.pivot(index="sample", columns="metabolite", values="tier")
        cols = [m for m in self.mets if m in val.columns]
        f = os.path.join(self.out, "GraphPad_export.xlsx")
        with pd.ExcelWriter(f) as xw:
            val.loc[self.samples, cols].to_excel(xw, sheet_name="Values_mM")
            tier.loc[self.samples, cols].to_excel(xw, sheet_name="Tiers")
            pd.DataFrame({"note": [f"Values shown only for: {', '.join(allowed)}. Blank = no reportable number.",
                                   "Rows = EXPNO; columns = metabolites; mM in the original sample (dilution applied)."]}
                         ).to_excel(xw, sheet_name="README", index=False)
        QMessageBox.information(self, "Exported", f"Saved {f}")


# ============================================================================ main window
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"MANC-Q {q.__version__}")
        self.spectra, self.folder, self.running = [], "", False
        self.tabs = QTabWidget()
        self.tabs.setStyleSheet(f"QTabBar::tab:selected {{ color:{ACCENT}; font-weight:600; }} "
                                "QTabBar::tab { padding: 6px 14px; }")
        self.spectra_tab = SpectraTab(self)
        self.process_tab = ProcessTab(self)
        self.regions_tab = RegionsTab(self)
        self.settings_tab = SettingsTab(self)
        self.run_tab = RunTab(self)
        self.results_tab = ResultsTab(self)
        for w, lab in [(self.spectra_tab, "1  Spectra"), (self.process_tab, "2  Process FIDs"),
                       (self.regions_tab, "3  Regions"), (self.settings_tab, "4  Settings"),
                       (self.run_tab, "5  Run"), (self.results_tab, "6  Results")]:
            self.tabs.addTab(w, lab)
        self.tabs.setTabVisible(1, False)
        self.tabs.currentChanged.connect(self.tab_changed)
        self.setCentralWidget(self.tabs)
        sb = QStatusBar()
        sb.showMessage(f"MANC-Q {q.__version__}  |  Metabolite Abundance by NMR Curve-fitting, Quantified  |  "
                       "University of Manchester")
        self.setStatusBar(sb)
        self.resize(1180, 760)

    def spectra_changed(self):
        has_fid = any(s["selected"] and not s["processed"] for s in self.spectra)
        self.tabs.setTabVisible(1, has_fid)

    def tab_changed(self, i):
        w = self.tabs.widget(i)
        if hasattr(w, "refresh"):
            w.refresh()

    def next_tab(self):
        i = self.tabs.currentIndex() + 1
        while i < self.tabs.count() and not self.tabs.isTabVisible(i):
            i += 1
        if i < self.tabs.count():
            self.tabs.setCurrentIndex(i)

    def prev_tab(self):
        i = self.tabs.currentIndex() - 1
        while i >= 0 and not self.tabs.isTabVisible(i):
            i -= 1
        if i >= 0:
            self.tabs.setCurrentIndex(i)

    def start_run(self):
        if self.running:
            warn(self, "A run is already in progress."); return
        sel = [s for s in self.spectra if s["selected"]]
        st = self.settings_tab
        if not sel:
            warn(self, "No spectra are selected (step 1)."); return
        if st.ref_mm.value() <= 0:
            warn(self, "Enter the reference (TSP/DSS) concentration in the tube (step 4)."); return
        try:
            single, per = st.dilutions()
        except ValueError as exc:
            warn(self, str(exc)); return
        out = st.out.text().strip()
        if not out:
            warn(self, "Choose a results folder (step 4)."); return
        try:
            os.makedirs(out, exist_ok=True)
        except OSError as exc:
            warn(self, f"Cannot create the results folder:\n{exc}"); return
        dil_file = None
        if per:
            dil_file = os.path.join(out, "dilutions_used.csv")
            pd.DataFrame(dict(sample=list(per), dilution=list(per.values()))).to_csv(dil_file, index=False)
        spectra = [dict(s, use_fid=not s["processed"]) for s in sel]
        settings = dict(REFERENCE=st.ref.currentText().split("-")[0], PROCNO=self.spectra_tab.procno.value(),
                        SAMPLES=[s["expno"] for s in sel], N_WORKERS=st.workers.value(),
                        OVERLAYS=st.overlays.isChecked(), RESUME=True)
        settings.update(self.regions_tab.engine_settings())
        json.dump(dict(folder=self.folder, settings=st.to_json(), regions=self.regions_tab.to_json(),
                       fid_processing={s["expno"]: dict(lb=s["fid_lb"], ph0=s["fid_ph0"], ph1=s["fid_ph1"])
                                       for s in spectra if s["use_fid"]}),
                  open(os.path.join(out, "MANC-Q settings used.json"), "w"), indent=2, default=str)
        job = dict(folder=self.folder if os.path.exists(os.path.join(self.folder, "acqus")) is False else
                   os.path.dirname(self.folder.rstrip("/\\")),
                   out=out, ref_mm=st.ref_mm.value(), dilution=single, dilution_file=dil_file,
                   spectra=spectra, settings=settings)
        self.running = True
        self.tabs.setCurrentWidget(self.run_tab)
        self.run_tab.start(job)

    def closeEvent(self, ev):
        if self.running and QMessageBox.question(self, "Run in progress",
                                                 "A run is in progress. Quit anyway? Finished spectra are kept.") \
                != QMessageBox.Yes:
            ev.ignore(); return
        if self.run_tab.worker is not None and self.run_tab.worker.isRunning():
            self.run_tab.worker.kill()
            self.run_tab.worker.wait(5000)
        ev.accept()


def make_app(argv=None):
    app = QApplication.instance() or QApplication(argv or sys.argv)
    app.setStyle("Fusion")
    f = QFont(); f.setFamilies(["Segoe UI", "Roboto", "Arial"]); f.setPointSize(10); app.setFont(f)
    pal = app.palette(); pal.setColor(QPalette.Highlight, QColor(ACCENT)); app.setPalette(pal)
    return app


def main():
    import multiprocessing
    multiprocessing.freeze_support()
    app = make_app()
    w = MainWindow(); w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
