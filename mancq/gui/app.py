"""
MANC-Q desktop interface.

Seven steps, one tab each: Spectra -> Process spectra -> Regions -> Settings -> Run -> Results -> Review (manual,
Chenomx-style adjustment of fitted compounds). The interface is a thin layer over the MANC-Q engine (package mancq):
everything it does can also be done from the command line (processing with --processing processing_used.csv,
adjustments with apply-edits), and both give identical numbers.
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

from PySide6.QtCore import Qt, QThread, Signal, QTimer, QUrl
from PySide6.QtGui import QFont, QColor, QPalette, QBrush, QAction, QKeySequence, QShortcut, QImage, QTextDocument
from PySide6.QtWidgets import (QApplication, QMainWindow, QTabWidget, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
                               QLabel, QPushButton, QLineEdit, QTableWidget, QTableWidgetItem, QGroupBox, QComboBox,
                               QDoubleSpinBox, QRadioButton, QCheckBox, QProgressBar, QPlainTextEdit, QSlider,
                               QHeaderView, QSplitter, QSpinBox, QStatusBar, QFileDialog, QMessageBox, QButtonGroup,
                               QAbstractItemView, QSizePolicy, QDialog, QTextBrowser, QSystemTrayIcon, QStyle)

from mancq import engine as q
from mancq import fidproc as fp
from mancq import results as rs
from mancq import twod


ACCENT = "#660099"
TIER_COL = {"Quantified": "#1d4ed8", "Overlapped (semi-quantitative)": "#b45309",
            "Deconvolution estimate (low confidence)": "#9ca3af", "Manually adjusted": "#0f766e",
            "Upper bound only": "#e5e7eb",
            "Not detected": "#ffffff", "Not measurable (region ignored)": "#f3f4f6"}
TIER_TXT = {"Quantified": "#ffffff", "Overlapped (semi-quantitative)": "#ffffff",
            "Deconvolution estimate (low confidence)": "#111827", "Manually adjusted": "#ffffff",
            "Upper bound only": "#374151",
            "Not detected": "#4b5563", "Not measurable (region ignored)": "#4b5563"}
MW = {"TSP-d4": 172.27, "DSS-d6": 224.36}
COMMON_REGIONS = [("Urea", 5.70, 5.85), ("DMSO", 2.69, 2.74), ("Methanol", 3.34, 3.38), ("Acetone", 2.21, 2.24)]
POLOX_MODES = ["Model it (recommended)", "Ignore", "Not in my samples"]
CACHE = os.path.join(os.path.expanduser("~"), ".mancq_cache")
# what the window remembers between sessions (plain JSON; delete the file to start afresh)
PREFS_FILE = os.environ.get("MANCQ_PREFS") or os.path.join(
    os.environ.get("APPDATA") or os.path.expanduser("~"), "MANC-Q", "window_settings.json")
HELP_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "help")


# ============================================================================ helpers
def open_path(path):
    if sys.platform.startswith("win"):
        os.startfile(path)                                  # noqa: windows only
    elif sys.platform == "darwin":
        subprocess.Popen(["open", path])
    else:
        subprocess.Popen(["xdg-open", path])


def load_prefs():
    try:
        return json.load(open(PREFS_FILE, encoding="utf-8"))
    except Exception:
        return {}


def save_prefs(d):
    try:
        os.makedirs(os.path.dirname(PREFS_FILE), exist_ok=True)
        json.dump(d, open(PREFS_FILE, "w", encoding="utf-8"), indent=1, default=str)
    except Exception:
        pass                                  # remembering settings is a convenience, never a reason to fail


class SortItem(QTableWidgetItem):
    """Table cell that sorts by a hidden key (numbers sort as numbers), as in File Explorer. Blank cells
    (blank=True) stay at the bottom whichever way the column is sorted."""
    def __init__(self, text, key=None, blank=False):
        super().__init__(text)
        self.key = text.lower() if key is None else key
        self.blank = blank

    def __lt__(self, other):
        ob = getattr(other, "blank", False)
        if self.blank or ob:
            t = self.tableWidget()
            desc = t is not None and t.horizontalHeader().sortIndicatorOrder() == Qt.DescendingOrder
            return self.blank and not ob if desc else ob and not self.blank
        try:
            return self.key < other.key
        except Exception:
            return str(self.key) < str(getattr(other, "key", ""))


def error_details(details=""):
    import platform
    from PySide6 import __version__ as pyside_version
    return (f"MANC-Q {q.__version__}, Python {platform.python_version()}, PySide6 {pyside_version}, "
            f"{platform.platform()}\n\n{details}").strip()


def error_box(parent, text, details=""):
    """An error message with the technical details behind "Show Details..." and a button that copies them,
    ready to paste into an e-mail or a GitHub issue."""
    box = QMessageBox(parent)
    box.setIcon(QMessageBox.Critical); box.setWindowTitle("MANC-Q")
    box.setText(text)
    full = error_details(details)
    box.setDetailedText(full)
    copy = box.addButton("Copy error details", QMessageBox.ActionRole)
    box.addButton(QMessageBox.Close)
    box.exec()
    if box.clickedButton() is copy:
        QApplication.clipboard().setText(full)
        QMessageBox.information(parent, "MANC-Q", "The error details are on the clipboard: paste them (Ctrl+V) "
                                "into an e-mail or a GitHub issue.")


def help_sections():
    """{step number: markdown} from the quick start guide shipped with MANC-Q."""
    f = os.path.join(HELP_DIR, "QUICKSTART_GUI.md")
    if not os.path.exists(f):
        return {}, ""
    text = open(f, encoding="utf-8").read()
    out = {}
    for part in text.split("\n## ")[1:]:
        num = part.split()[0]
        if num.isdigit():
            out[int(num)] = "## " + part
    return out, text


def show_help(parent, step):
    sections, whole = help_sections()
    dlg = QDialog(parent); dlg.setWindowTitle("MANC-Q help"); dlg.resize(860, 720)
    v = QVBoxLayout(dlg)
    tb = QTextBrowser(); tb.setOpenExternalLinks(True); tb.setSearchPaths([HELP_DIR])
    doc = tb.document()
    imgs = os.path.join(HELP_DIR, "images")
    if os.path.isdir(imgs):                   # screenshots scaled to fit the window
        for fn in os.listdir(imgs):
            im = QImage(os.path.join(imgs, fn))
            if not im.isNull():
                doc.addResource(QTextDocument.ImageResource, QUrl(f"images/{fn}"),
                                im.scaledToWidth(780, Qt.SmoothTransformation))
    if not whole:
        tb.setPlainText("The guide was not found. It is in docs/QUICKSTART_GUI.md on the MANC-Q GitHub page.")
    else:
        tb.setMarkdown(sections.get(step, whole))
    v.addWidget(tb, 1)
    h = QHBoxLayout()
    ball = QPushButton("Show the whole guide"); ball.clicked.connect(lambda: tb.setMarkdown(whole)); h.addWidget(ball)
    h.addStretch()
    bc = QPushButton("Close"); bc.clicked.connect(dlg.accept); h.addWidget(bc)
    v.addLayout(h)
    dlg.exec()


def heading(text, sub=None):
    box = QWidget(); v = QVBoxLayout(box); v.setContentsMargins(0, 0, 0, 6)
    top = QHBoxLayout()
    t = QLabel(text); t.setStyleSheet("font-size: 16px; font-weight: 600;"); t.setObjectName("step_heading")
    top.addWidget(t); top.addStretch()
    step = int(text.split()[0]) if text.split()[0].isdigit() else 0
    hb = QPushButton("Help"); hb.setToolTip("Open the guide for this step")
    hb.clicked.connect(lambda: show_help(box, step))
    top.addWidget(hb)
    v.addLayout(top)
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
        # processing: TopSpin spectra start "as is" (phase changes are relative to TopSpin's, lb read later from
        # procs); raw FIDs start with automatic phasing
        out.append(dict(expno=os.path.basename(d), path=d, title=title, field=field, processed=proc, fid=fid,
                        selected=True, snr=None, fwhm=None, check="", proc_asis=proc, proc_lb=None,
                        proc_ph0=0.0 if proc else None, proc_ph1=0.0, proc_baseline=False, proc_checked=False))
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
        self.table.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.table.itemChanged.connect(self.tick_changed)
        self.table.itemSelectionChanged.connect(self.preview_selected)
        split.addWidget(self.table)
        right = QWidget(); rv = QVBoxLayout(right); rv.setContentsMargins(0, 0, 0, 0)
        ch = QHBoxLayout()
        hint = QLabel("Select several rows (Ctrl or Shift + click) to compare them. Scroll over the plot to zoom.")
        hint.setStyleSheet("color:#4b5563;"); hint.setWordWrap(True); ch.addWidget(hint, 1)
        self.scale = QCheckBox("Scale to the reference peak"); self.scale.setChecked(True)
        self.scale.setToolTip("Divide each spectrum by the height of its TSP/DSS peak, so spectra with different "
                              "numbers of scans or receiver gain can be compared")
        self.scale.toggled.connect(self.preview_selected); ch.addWidget(self.scale)
        rv.addLayout(ch)
        self.canvas = Canvas(6.2, 3.8)
        self.canvas.mpl_connect("scroll_event", self.scroll)
        rv.addWidget(self.canvas, 1)
        split.addWidget(right); split.setSizes([560, 540])
        self.cache = {}
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
        if not (os.path.exists(os.path.join(spec, "10", "pdata", "1", "1i"))      # older demos had no 1i
                and os.path.exists(os.path.join(spec, "11", "fid"))
                and os.path.exists(os.path.join(spec, "12", "pdata", "1", "2rr"))):  # nor a 2D (before 1.3)
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
                                "EXPNO 11 is a raw FID so you can practise phasing (EXPNO 10 can be adjusted too: "
                                "untick \"Use TopSpin's processing as is\"). EXPNO 12 is a 2D TOCSY of EXPNO 10: after "
                                "the run, click 2D in the Results step to see it. Click Next through the steps and "
                                "Start; the true concentrations are listed in mancq/demo.py. Afterwards, step 7 lets "
                                "you review and adjust the fit of each compound.")

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
        self.win.remember()
        if spectra:
            self.table.selectRow(0)

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
                it.setToolTip(val)
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

    def spectrum(self, s):
        """(ppm referenced to the reference peak, intensity, reference height) of a processed spectrum, cached."""
        key = (s["path"], self.procno.value())
        if key not in self.cache:
            ppm, y, sf, _ = q.read_spectrum(s["path"], self.procno.value())
            ppm = fp.reference(ppm, y)
            ref = float(y[np.abs(ppm) < 0.02].max()) if (np.abs(ppm) < 0.02).any() else float(y.max())
            self.cache[key] = (ppm, y, ref)
        return self.cache[key]

    def preview(self, row):
        self.table.selectRow(row)

    def preview_selected(self):
        sp = self.win.spectra
        rows = sorted({i.row() for i in self.table.selectedIndexes() if 0 <= i.row() < len(sp)})
        self.canvas.fig.clear()
        ax = self.canvas.fig.add_subplot(111)
        shown, fids, top = [], [], 0.0
        for k, r in enumerate(rows[:12]):
            s = sp[r]
            if not s["processed"]:
                fids.append(s["expno"]); continue
            try:
                ppm, y, ref = self.spectrum(s)
            except Exception as exc:
                ax.text(0.5, 0.5, f"Cannot read EXPNO {s['expno']}:\n{exc}", ha="center", transform=ax.transAxes)
                continue
            yy = y / ref if (self.scale.isChecked() and len(rows) > 1 and ref > 0) else y
            col = "k" if len(rows) == 1 else q.PALETTE[len(shown) % len(q.PALETTE)]
            ax.plot(ppm, yy, color=col, lw=0.5, label=f"EXPNO {s['expno']}")
            m = (ppm > 0.5) & (ppm < 10)
            top = max(top, float(yy[m].max()) if m.any() else float(yy.max()))
            shown.append(s)
        if shown:
            ax.set_xlim(10, -0.5)
            ax.set_ylim(-0.01 * top, 1.05 * top if len(shown) > 1 else 0.08 * self.spectrum(shown[0])[1].max())
            if len(shown) == 1:
                ax.set_title(f"EXPNO {shown[0]['expno']}  {shown[0]['title']}", fontsize=9, loc="left")
            else:
                ax.legend(fontsize=7, frameon=False, loc="upper left")
                ax.set_title(f"{len(shown)} spectra" + (", each divided by its reference peak"
                                                         if self.scale.isChecked() else ""), fontsize=9, loc="left")
            self.top = top if len(shown) > 1 else 0.08 * self.spectrum(shown[0])[1].max()
        elif fids:
            ax.text(0.5, 0.5, "Raw FID only.\nIt is processed and checked in step 2.", ha="center",
                    va="center", transform=ax.transAxes, color="#6b7280")
        if shown and fids:
            ax.text(0.99, 0.98, f"Raw FIDs not shown: {', '.join(fids)}", ha="right", va="top", fontsize=7,
                    transform=ax.transAxes, color="#6b7280")
        ax.set_xlabel("ppm"); tidy(ax)
        self.canvas.fig.tight_layout(); self.canvas.draw()

    def scroll(self, ev):
        ax = ev.inaxes
        if ax is None or ev.xdata is None:
            return
        f = 0.8 if ev.button == "up" else 1.25
        x0, x1 = ax.get_xlim()
        c = ev.xdata
        ax.set_xlim(c + (x0 - c) * f, c + (x1 - c) * f)
        lo, hi = sorted(ax.get_xlim())
        tops = []
        for line in ax.get_lines():
            x, y = line.get_xdata(), line.get_ydata()
            m = (x > lo) & (x < hi) & ~((x > 4.6) & (x < 5.0))
            if m.sum() > 2:
                tops.append(np.nanmax(y[m]))
        if tops and hi - lo < 9:
            ax.set_ylim(-0.02 * max(tops), 1.05 * max(tops))
        self.canvas.draw_idle()


# ============================================================================ tab 2: process spectra
class ProcessTab(QWidget):
    """Every selected spectrum. TopSpin spectra are used as TopSpin processed them unless the user unticks 'as is';
    changes then start from TopSpin's own result. Raw FIDs are always processed here."""
    ZOOMS = [("Zoom 0.8-1.6 ppm (methyl region)", (1.6, 0.8)), ("Zoom 3.0-4.2 ppm (sugar region)", (4.2, 3.0)),
             ("Zoom 6.8-8.5 ppm (aromatic region)", (8.5, 6.8)), ("Zoom -0.1-0.1 ppm (reference)", (0.1, -0.1))]

    def __init__(self, win):
        super().__init__()
        self.win = win
        self.prep, self.cur = {}, None
        v = QVBoxLayout(self)
        v.addWidget(heading("2  Process spectra", "Spectra processed in TopSpin are used exactly as TopSpin processed "
                            "them, unless you untick \"Use TopSpin's processing as is\": phase changes then start from "
                            "TopSpin's phasing, and the line broadening you set replaces TopSpin's. Raw FIDs are "
                            "processed here (digital-filter correction, line broadening, Fourier transform, automatic "
                            "phasing). Check each spectrum you change; a phasing error goes straight into every "
                            "concentration."))
        split = QSplitter(Qt.Horizontal)
        left = QWidget(); lv = QVBoxLayout(left); lv.setContentsMargins(0, 0, 0, 0)
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["EXPNO", "Title", "Data", "Processing", "Checked"])
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.currentCellChanged.connect(lambda r, *a: self.select(r))
        lv.addWidget(self.table)
        self.asis = QCheckBox("Use TopSpin's processing as is")
        self.asis.toggled.connect(self.asis_toggled)
        lv.addWidget(self.asis)
        self.asis_note = QLabel(); self.asis_note.setWordWrap(True); self.asis_note.setStyleSheet("color:#4b5563;")
        lv.addWidget(self.asis_note)
        g = QGroupBox("Phase"); gl = QGridLayout(g)
        self.s0 = QSlider(Qt.Horizontal); self.s0.setRange(0, 3600)
        self.s1 = QSlider(Qt.Horizontal); self.s1.setRange(-900, 900)
        self.n0, self.n1 = QLabel("Zero order (ph0)"), QLabel("First order (ph1)")
        self.l0, self.l1 = QLabel(), QLabel()
        gl.addWidget(self.n0, 0, 0); gl.addWidget(self.s0, 0, 1); gl.addWidget(self.l0, 0, 2)
        gl.addWidget(self.n1, 1, 0); gl.addWidget(self.s1, 1, 1); gl.addWidget(self.l1, 1, 2)
        for sl in (self.s0, self.s1):
            sl.valueChanged.connect(self.phase_moved)
            sl.sliderReleased.connect(self.redraw_full)
        self.b_auto = QPushButton("Auto-phase this spectrum"); self.b_auto.clicked.connect(self.auto_one)
        gl.addWidget(self.b_auto, 2, 1)
        for sl in (self.s0, self.s1):
            sl.setSingleStep(1); sl.setPageStep(10)
        kh = QLabel("Click a slider, then use the arrow keys for 0.1 degree steps (Page Up / Page Down: 1 degree).")
        kh.setStyleSheet("color:#4b5563;"); kh.setWordWrap(True)
        gl.addWidget(kh, 3, 0, 1, 3)
        self.phase_box = g
        lv.addWidget(g)
        g2 = QGroupBox("Processing"); g2l = QGridLayout(g2)
        g2l.addWidget(QLabel("Line broadening"), 0, 0)
        self.lb = QDoubleSpinBox(); self.lb.setRange(0, 5); self.lb.setSingleStep(0.1); self.lb.setValue(0.3)
        self.lb.setSuffix(" Hz"); self.lb.editingFinished.connect(self.lb_changed); g2l.addWidget(self.lb, 0, 1)
        self.lb_note = QLabel(""); self.lb_note.setStyleSheet("color:#4b5563;"); g2l.addWidget(self.lb_note, 0, 2)
        self.base = QCheckBox("Baseline correction (usually leave off: MANC-Q fits a smooth baseline itself)")
        self.base.toggled.connect(self.baseline_toggled); g2l.addWidget(self.base, 1, 0, 1, 3)
        g2l.addWidget(QLabel("Reference"), 2, 0); g2l.addWidget(QLabel("TSP/DSS set to 0.000 ppm automatically"), 2, 1, 1, 2)
        # inside the box, so it is greyed out with it when this spectrum is used as TopSpin processed it
        b2 = QPushButton("Use this line broadening for every spectrum MANC-Q processes"); b2.clicked.connect(self.lb_all)
        g2l.addWidget(b2, 3, 0, 1, 3)
        self.proc_box = g2
        lv.addWidget(g2)
        self.ok = QCheckBox("I have checked the processing of this spectrum"); self.ok.toggled.connect(self.checked)
        lv.addWidget(self.ok)
        self.zoom = QComboBox()
        for lab, _ in self.ZOOMS:
            self.zoom.addItem(lab)
        self.zoom.currentIndexChanged.connect(self.redraw_full)
        lv.addWidget(self.zoom)
        hb = QHBoxLayout()
        ba = QPushButton("Auto-phase all raw FIDs"); ba.clicked.connect(self.auto_all)
        ba.setToolTip("Automatic phasing of every raw FID (TopSpin spectra are not changed)")
        hb.addWidget(ba)
        bm = QPushButton("Mark all as checked"); bm.clicked.connect(self.check_all)
        bm.setToolTip("Only after you have looked at every spectrum MANC-Q processes")
        hb.addWidget(bm)
        lv.addLayout(hb)
        lv.addStretch()
        split.addWidget(left)
        self.canvas = Canvas(6.2, 4.6); split.addWidget(self.canvas); split.setSizes([520, 600])
        v.addWidget(split, 1)
        h = QHBoxLayout(); bk = QPushButton("Back"); bk.clicked.connect(win.prev_tab); h.addWidget(bk); h.addStretch()
        nx = primary("Next"); nx.clicked.connect(self.go_next); h.addWidget(nx); v.addLayout(h)

    # ---- which spectra, and how each is processed
    def spectra(self):
        return [s for s in self.win.spectra if s["selected"]]

    @staticmethod
    def is_topspin(s):
        return s["processed"]

    def procno(self):
        return self.win.spectra_tab.procno.value()

    def summary(self, s):
        if self.is_topspin(s) and s["proc_asis"]:
            return "As in TopSpin"
        lb = s["proc_lb"]
        bl = ", baseline corr." if s["proc_baseline"] else ""
        if self.is_topspin(s):
            return f"Changed: ph0 {s['proc_ph0']:+.1f}, ph1 {s['proc_ph1']:+.1f}, LB {lb:.2f} Hz{bl}"
        ph = "auto (not yet run)" if s["proc_ph0"] is None else f"ph0 {s['proc_ph0']:.1f}, ph1 {s['proc_ph1']:.1f}"
        return f"{ph}, LB {lb:.2f} Hz{bl}"

    def needs_check(self, s):
        return not (self.is_topspin(s) and s["proc_asis"])

    def refresh(self):
        f = self.spectra()
        self.table.setRowCount(len(f))
        for i, s in enumerate(f):
            self._ensure_params(s)
            for j, val in enumerate([s["expno"], s["title"], "TopSpin" if self.is_topspin(s) else "Raw FID",
                                     self.summary(s), ""]):
                it = QTableWidgetItem(val); it.setToolTip(val)
                self.table.setItem(i, j, it)
            self._check_cell(i, s)
        if f and (self.cur is None or self.cur not in f):
            self.cur = None
            self.table.selectRow(0); self.select(0)
        elif not f:
            self.cur = None
            self.canvas.fig.clear(); self.canvas.draw()

    def _check_cell(self, i, s):
        if not self.needs_check(s):
            it = QTableWidgetItem("not needed")
            it.setForeground(QBrush(QColor("#4b5563")))
        else:
            it = QTableWidgetItem("yes" if s["proc_checked"] else "no")
            it.setForeground(QBrush(QColor("#15803d" if s["proc_checked"] else "#b45309")))
        self.table.setItem(i, 4, it)

    def _ensure_params(self, s):
        """TopSpin's own line broadening, read once per spectrum, is the starting value."""
        if s.get("proc_lb") is None:
            s["proc_lb"] = 0.3
            if self.is_topspin(s):
                try:
                    s["ts_lb"] = fp.topspin_params(s["path"], self.procno())["lb"]
                    s["proc_lb"] = s["ts_lb"]
                except Exception:
                    s["ts_lb"] = None
        if self.is_topspin(s) and "adjustable" not in s:
            s["adjustable"] = fp.can_adjust(s["path"], self.procno())

    def _prep(self, s):
        key = (s["path"], self.is_topspin(s), round(s["proc_lb"], 4), self.procno())
        if key not in self.prep:
            QApplication.setOverrideCursor(Qt.WaitCursor)
            try:
                if self.is_topspin(s):
                    self.prep[key] = fp.prepare_topspin(s["path"], self.procno(), lb=s["proc_lb"])
                else:
                    self.prep[key] = fp.prepare(s["path"], lb=s["proc_lb"])
            finally:
                QApplication.restoreOverrideCursor()
        return self.prep[key]

    def select(self, row):
        f = self.spectra()
        if not (0 <= row < len(f)):
            return
        s = self.cur = f[row]
        self._ensure_params(s)
        ts = self.is_topspin(s)
        self.asis.blockSignals(True)
        self.asis.setVisible(ts); self.asis.setChecked(ts and s["proc_asis"])
        self.asis.setEnabled(ts and s.get("adjustable", False))
        self.asis.blockSignals(False)
        if not ts:
            self.asis_note.setText("Raw FID: MANC-Q processes it with the settings below.")
        elif not s.get("adjustable", False):
            self.asis_note.setText("This spectrum can only be used as TopSpin processed it: its imaginary part (1i) "
                                   "and the raw FID are both missing.")
        elif fp.has_imaginary(s["path"], self.procno()):
            self.asis_note.setText("Untick to change the phase or line broadening, starting from TopSpin's result.")
        else:
            self.asis_note.setText("Untick to change the phase or line broadening. The imaginary part (1i) is missing, "
                                   "so MANC-Q processes the raw FID and phases it to match TopSpin first.")
        if ts:
            self.n0.setText("Zero-order change (ph0)"); self.n1.setText("First-order change (ph1)")
            self.s0.blockSignals(True); self.s0.setRange(-1800, 1800); self.s0.blockSignals(False)
            self.b_auto.setText("Auto-phase (change from TopSpin)")
            self.lb_note.setText(f"TopSpin used {s['ts_lb']:.2f} Hz" if s.get("ts_lb") is not None else "")
        else:
            self.n0.setText("Zero order (ph0)"); self.n1.setText("First order (ph1)")
            self.s0.blockSignals(True); self.s0.setRange(0, 3600); self.s0.blockSignals(False)
            self.b_auto.setText("Auto-phase this spectrum again")
            self.lb_note.setText("")
        if not ts and s["proc_ph0"] is None:
            try:
                p = self._prep(s)
            except Exception as exc:
                warn(self, f"Could not read the FID of EXPNO {s['expno']}:\n{exc}"); return
            QApplication.setOverrideCursor(Qt.WaitCursor)
            try:
                s["proc_ph0"], s["proc_ph1"] = fp.autophase(p["F"], p["keep"])
            finally:
                QApplication.restoreOverrideCursor()
        for sl, val in ((self.s0, s["proc_ph0"] * 10), (self.s1, s["proc_ph1"] * 10)):
            sl.blockSignals(True); sl.setValue(int(round(val))); sl.blockSignals(False)
        self.lb.blockSignals(True); self.lb.setValue(s["proc_lb"]); self.lb.blockSignals(False)
        self.base.blockSignals(True); self.base.setChecked(s["proc_baseline"]); self.base.blockSignals(False)
        self.ok.blockSignals(True); self.ok.setChecked(s["proc_checked"]); self.ok.blockSignals(False)
        self.set_enabled()
        self.update_labels(); self.redraw_full(); self.refresh_row()

    def set_enabled(self):
        s = self.cur
        on = s is not None and not (self.is_topspin(s) and s["proc_asis"])
        self.phase_box.setEnabled(on); self.proc_box.setEnabled(on)
        self.ok.setEnabled(on)

    def update_labels(self):
        self.l0.setText(f"{self.s0.value() / 10:+.1f} deg" if self.cur is not None and self.is_topspin(self.cur)
                        else f"{self.s0.value() / 10:.1f} deg")
        self.l1.setText(f"{self.s1.value() / 10:+.1f} deg")

    def changed(self):
        """Any change to a spectrum's processing means it has to be checked again."""
        self.cur["proc_checked"] = False
        self.ok.blockSignals(True); self.ok.setChecked(False); self.ok.blockSignals(False)

    def asis_toggled(self, on):
        if self.cur is None:
            return
        self.cur["proc_asis"] = on
        if not on:
            self.changed()
        self.set_enabled(); self.redraw_full(); self.refresh_row()

    def phase_moved(self):
        if self.cur is None:
            return
        self.cur["proc_ph0"], self.cur["proc_ph1"] = self.s0.value() / 10, self.s1.value() / 10
        self.changed()
        self.update_labels()
        self.draw()
        self.refresh_row()

    def redraw_full(self):
        self.draw()

    def draw(self):
        s = self.cur
        if s is None:
            return
        ts = self.is_topspin(s)
        before = None
        try:
            if ts and (s["proc_asis"] or not s.get("adjustable", False)):
                if s.get("adjustable", False):
                    r = fp.finish(self._prep(dict(s, proc_lb=s.get("ts_lb") if s.get("ts_lb") is not None
                                                   else s["proc_lb"])), 0.0, 0.0)
                else:
                    ppm, y, _, _ = q.read_spectrum(s["path"], self.procno())
                    r = dict(ppm=fp.reference(ppm, y), y=y)
            else:
                r = fp.finish(self._prep(s), s["proc_ph0"], s["proc_ph1"], baseline=s["proc_baseline"])
                if ts:
                    before = fp.finish(self._prep(dict(s, proc_lb=s.get("ts_lb") if s.get("ts_lb") is not None
                                                        else s["proc_lb"])), 0.0, 0.0)
        except Exception as exc:
            self.canvas.fig.clear()
            ax = self.canvas.fig.add_subplot(111)
            ax.text(0.5, 0.5, f"Cannot process EXPNO {s['expno']}:\n{exc}", ha="center", va="center",
                    transform=ax.transAxes, wrap=True)
            tidy(ax); self.canvas.draw(); return
        zl = self.ZOOMS[self.zoom.currentIndex()][1]
        self.canvas.fig.clear()
        ax1 = self.canvas.fig.add_subplot(211); ax2 = self.canvas.fig.add_subplot(212)
        title = f"EXPNO {s['expno']}: whole spectrum" + ("" if before is None else "   (grey = TopSpin as is)")
        for ax, (lo, hi), top, lab in [(ax1, (10, -0.5), 0.08, title),
                                       (ax2, zl, 1.05, "Zoom: the baseline either side of each peak should be flat")]:
            m = (r["ppm"] < lo) & (r["ppm"] > hi)
            if m.sum() < 2:
                continue
            if before is not None:
                mb = (before["ppm"] < lo) & (before["ppm"] > hi)
                ax.plot(before["ppm"][mb], before["y"][mb], color="#9ca3af", lw=0.6)
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
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            if self.is_topspin(self.cur):
                self.cur["proc_ph0"], self.cur["proc_ph1"] = fp.autophase_relative(p)
            else:
                self.cur["proc_ph0"], self.cur["proc_ph1"] = fp.autophase(p["F"], p["keep"])
        finally:
            QApplication.restoreOverrideCursor()
        self.changed()
        self.select(self.spectra().index(self.cur))

    def lb_changed(self):
        if self.cur is not None and abs(self.cur["proc_lb"] - self.lb.value()) > 1e-9:
            self.cur["proc_lb"] = self.lb.value(); self.changed(); self.redraw_full(); self.refresh_row()

    def baseline_toggled(self, on):
        if self.cur is not None:
            self.cur["proc_baseline"] = on; self.changed(); self.redraw_full(); self.refresh_row()

    def lb_all(self):
        if self.cur is None or not self.needs_check(self.cur):
            return
        for s in self.spectra():
            if self.needs_check(s) and abs(s["proc_lb"] - self.lb.value()) > 1e-9:
                s["proc_lb"] = self.lb.value(); s["proc_checked"] = False
        self.refresh()
        if self.cur is not None:
            self.select(self.spectra().index(self.cur))

    def checked(self, on):
        if self.cur is not None:
            self.cur["proc_checked"] = on; self.refresh_row()

    def refresh_row(self):
        f = self.spectra()
        if self.cur in f:
            i = f.index(self.cur)
            it = QTableWidgetItem(self.summary(self.cur)); it.setToolTip(it.text())
            self.table.setItem(i, 3, it)
            self._check_cell(i, self.cur)

    def auto_all(self):
        fids = [s for s in self.spectra() if not self.is_topspin(s)]
        if not fids:
            QMessageBox.information(self, "Auto-phase", "There are no raw FIDs among the selected spectra."); return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            for s in fids:
                self._ensure_params(s)
                p = self._prep(s)
                s["proc_ph0"], s["proc_ph1"] = fp.autophase(p["F"], p["keep"])
                s["proc_checked"] = False
        except Exception as exc:
            QApplication.restoreOverrideCursor()
            error_box(self, f"Automatic phasing failed for EXPNO {s['expno']}.", traceback.format_exc()); return
        QApplication.restoreOverrideCursor()
        self.refresh()
        if self.cur is not None:
            self.select(self.spectra().index(self.cur))

    def check_all(self):
        todo = [s for s in self.spectra() if self.needs_check(s) and not s["proc_checked"]]
        if not todo:
            return
        if QMessageBox.question(self, "Mark all as checked",
                                f"Mark {len(todo)} spectra as checked? Only do this if you have looked at the "
                                "phasing and baseline of each one.") != QMessageBox.Yes:
            return
        for s in todo:
            s["proc_checked"] = True
        self.refresh()
        if self.cur is not None:
            self.select(self.spectra().index(self.cur))

    def go_next(self):
        if any(self.needs_check(s) and not s["proc_checked"] for s in self.spectra()):
            if QMessageBox.question(self, "Processing not checked",
                                    "Not every spectrum MANC-Q processes has been marked as checked. "
                                    "Continue anyway?") != QMessageBox.Yes:
                return
        self.win.next_tab()

    def engine_processing(self):
        """The PROCESSING setting for the engine: every spectrum MANC-Q processes itself."""
        out = {}
        for s in self.spectra():
            self._ensure_params(s)
            if self.is_topspin(s):
                if s["proc_asis"] or not s.get("adjustable", False):
                    continue
                out[s["expno"]] = dict(source="topspin", lb=s["proc_lb"], ph0=s["proc_ph0"], ph1=s["proc_ph1"],
                                       baseline=s["proc_baseline"])
            else:
                out[s["expno"]] = dict(source="fid", lb=s["proc_lb"], ph0=s["proc_ph0"], ph1=s["proc_ph1"],
                                       baseline=s["proc_baseline"])
        return out


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
        self.info.setStyleSheet("background:#f5f3ff; color:#1f2937; border:1px solid #ddd6fe; padding:8px;")
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
        ncpu = os.cpu_count() or 1
        g3l.addWidget(QLabel("CPU cores to use"), 1, 0)
        self.workers = QSpinBox(); self.workers.setRange(1, max(1, ncpu))
        self.workers.setValue(max(1, ncpu - 1)); g3l.addWidget(self.workers, 1, 1)
        cn = QLabel(f"This PC has {ncpu}. Each core fits one spectrum at a time, so more cores only help with "
                    "several spectra. Use fewer to keep the PC responsive for other work.")
        cn.setWordWrap(True); cn.setStyleSheet("color:#4b5563;")
        g3l.addWidget(cn, 2, 1, 1, 2)
        self.workers.setToolTip("Number of spectra fitted at the same time (one per CPU core)")
        self.overlays = QCheckBox("Write an overlay PDF for every spectrum (recommended)"); self.overlays.setChecked(True)
        g3l.addWidget(self.overlays, 3, 0, 1, 3)
        self.notify = QCheckBox("Tell me when a run finishes (Windows notification and a sound)")
        self.notify.setChecked(True); g3l.addWidget(self.notify, 4, 0, 1, 3)
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
                    dilution=self.dil.value(), workers=self.workers.value(), overlays=self.overlays.isChecked(),
                    notify=self.notify.isChecked())

    def from_json(self, d):
        self.ref.setCurrentText(d.get("reference", "TSP-d4")); self.ref_mm.setValue(d.get("ref_mm", 0.0))
        {"none": self.r_none, "same": self.r_same, "each": self.r_each}[d.get("dilution_mode", "same")].setChecked(True)
        self.dil.setValue(d.get("dilution", 1.0)); self.workers.setValue(d.get("workers", self.workers.value()))
        self.overlays.setChecked(d.get("overlays", True))
        self.notify.setChecked(d.get("notify", True))

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
            what = "FID" if info.get("source", "fid") == "fid" else "spectrum"
            self.set_status(smp, f"Processing {what} ({info['k']} of {info['n']})", ACCENT)
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
                bad = len(info.get("failed", []))
                self.notify("MANC-Q run finished", f"{self.n_total - bad} of {self.n_total} spectra fitted"
                            + (f", {bad} failed" if bad else "") + ". Click View results.")
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
            self.win.setWindowTitle(f"MANC-Q {q.__version__}  ({self.n_done} of {self.n_total} done)")

    def do_cancel(self):
        if self.worker is not None:
            self.worker.stop(); self.cancel.setEnabled(False)
            self.overall.setText("Cancelling: spectra already running will finish first...")

    def notify(self, title, text):
        """Taskbar flash, a sound and a Windows notification (if switched on in Settings)."""
        if not self.win.settings_tab.notify.isChecked():
            return
        QApplication.alert(self.win, 0)
        QApplication.beep()
        try:
            if QSystemTrayIcon.isSystemTrayAvailable():
                self.tray = QSystemTrayIcon(self.win.windowIcon() if not self.win.windowIcon().isNull() else
                                            self.style().standardIcon(QStyle.SP_DialogApplyButton), self)
                self.tray.show()
                self.tray.showMessage(title, text, QSystemTrayIcon.Information, 15000)
                QTimer.singleShot(20000, self.tray.hide)
        except Exception:
            pass

    def on_failed(self, msg):
        self.overall.setText("The run stopped with an error.")
        self.log.appendPlainText("\n" + msg)
        self.notify("MANC-Q run stopped", "The run stopped with an error.")
        error_box(self, "The run stopped with an error. The log on this page shows where.", msg)

    def on_thread_end(self):
        self.timer.stop(); self.tick(); self.cancel.setEnabled(False)
        self.win.setWindowTitle(f"MANC-Q {q.__version__}")
        self.win.running = False
        if rs.is_results_folder(self.out or ""):
            self.view.setEnabled(True)

    def show_results(self):
        self.win.results_tab.load(self.out)
        self.win.tabs.setCurrentWidget(self.win.results_tab)


# ============================================================================ tab 6: results
class ResultsTab(QWidget):
    FILTERS = [("All, including upper bounds", None),
               ("Quantified, Overlapped, Estimates, Manual", ["Quantified", "Overlapped (semi-quantitative)",
                                                              "Deconvolution estimate (low confidence)",
                                                              "Manually adjusted"]),
               ("Quantified, Overlapped, Manual", ["Quantified", "Overlapped (semi-quantitative)",
                                                   "Manually adjusted"]),
               ("Quantified only", ["Quantified"])]

    def __init__(self, win):
        super().__init__()
        self.win = win
        self.out, self.data = None, None
        v = QVBoxLayout(self)
        v.addWidget(heading("6  Results", "Colour shows how far each number can be trusted. Click a cell to see the "
                            "fit it came from; double-click it to review and adjust that fit (next step)."))
        bar = QHBoxLayout()
        bo = QPushButton("Open a results folder..."); bo.clicked.connect(self.choose); bar.addWidget(bo)
        self.recent = QComboBox(); self.recent.setMinimumWidth(170); self.recent.activated.connect(self.open_recent)
        self.recent.setToolTip("Results folders you opened before")
        bar.addWidget(self.recent)
        bar.addWidget(QLabel("Show:"))
        self.filter = QComboBox(); self.filter.addItems([f[0] for f in self.FILTERS]); self.filter.setCurrentIndex(2)
        self.filter.currentIndexChanged.connect(self.fill); bar.addWidget(self.filter)
        bar.addSpacing(12)
        for lab, t in [("Quantified", "Quantified"), ("Overlapped", "Overlapped (semi-quantitative)"),
                       ("Estimate", "Deconvolution estimate (low confidence)"), ("Manual", "Manually adjusted"),
                       ("Upper bound", "Upper bound only")]:
            sw = QLabel(f"  {lab}  ")
            sw.setStyleSheet(f"background:{TIER_COL[t]}; color:{TIER_TXT[t]}; border:1px solid #d1d5db; padding:2px;")
            bar.addWidget(sw)
        bar.addStretch(); v.addLayout(bar)
        self.where = QLabel(""); self.where.setStyleSheet("color:#4b5563;"); v.addWidget(self.where)
        split = QSplitter(Qt.Horizontal)
        self.table = QTableWidget(0, 0); self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.currentCellChanged.connect(lambda r, c, *a: self.cell(r, c))
        self.table.cellDoubleClicked.connect(lambda r, c: self.review(r, c))
        self.table.setContextMenuPolicy(Qt.ActionsContextMenu)
        a1 = QAction("Copy (as shown)", self.table); a1.setShortcut(QKeySequence.Copy)
        a1.setShortcutContext(Qt.WidgetShortcut); a1.triggered.connect(lambda: self.copy_cells(False))
        a2 = QAction("Copy numbers only (for Excel, Prism)", self.table); a2.triggered.connect(lambda: self.copy_cells(True))
        a3 = QAction("Review and adjust this fit", self.table)
        a3.triggered.connect(lambda: self.review(self.table.currentRow(), self.table.currentColumn()))
        for a in (a1, a2, a3):
            self.table.addAction(a)
        split.addWidget(self.table)
        split.addWidget(self.plot_panel()); split.setSizes([540, 620])
        v.addWidget(split, 1)
        h = QHBoxLayout()
        for lab, fn in [("Open Excel workbook", self.open_xlsx), ("Open results folder", self.open_folder),
                        ("Export for GraphPad", self.export_graphpad), ("Open overlay PDF (this sample)", self.open_overlay),
                        ("Save this plot...", self.save_plot),
                        ("Review and adjust this fit", lambda: self.review(self.table.currentRow(),
                                                                           self.table.currentColumn()))]:
            b = QPushButton(lab); b.clicked.connect(fn); h.addWidget(b)
        h.addStretch(); v.addLayout(h)

    def plot_panel(self):
        """The plot on the right: 1D or 2D, one metabolite or the whole spectrum, with zoom and pan."""
        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
        self.spec_cache, self.twod_cache, self.peak_cache = {}, {}, {}
        self.cur, self.axes1d = None, None
        right = QWidget(); rv = QVBoxLayout(right); rv.setContentsMargins(0, 0, 0, 0)
        row = QHBoxLayout()

        def toggles(labels, tips):
            g = QButtonGroup(self); g.setExclusive(True)
            for i, (lab, tip) in enumerate(zip(labels, tips)):
                b = QPushButton(lab); b.setCheckable(True); b.setToolTip(tip); b.setChecked(i == 0)
                g.addButton(b, i); row.addWidget(b)
            g.idClicked.connect(lambda *a: self.view_changed())
            return g
        self.dim = toggles(["1D", "2D"], ["The 1D spectrum and its fit",
                                          "The 2D spectrum (TOCSY or COSY) of this sample, with the cross peaks this "
                                          "metabolite should give: green = seen, red = missing"])
        row.addSpacing(10)
        self.span = toggles(["This metabolite", "Whole spectrum"],
                            ["Zoom to the selected metabolite", "The whole spectrum, with the selected metabolite "
                             "filled in purple and its peaks marked with arrows"])
        row.addSpacing(10)
        row.addWidget(QLabel("Peak"))
        self.peak = QComboBox(); self.peak.setMinimumWidth(150)
        self.peak.setToolTip("Which of the metabolite's peaks (multiplets) to zoom to")
        self.peak.currentIndexChanged.connect(lambda *a: self.draw_plot())
        row.addWidget(self.peak); row.addStretch()
        rv.addLayout(row)
        self.twod_row = QWidget(); r2 = QHBoxLayout(self.twod_row); r2.setContentsMargins(0, 0, 0, 0)
        self.twod_label = QLabel(""); self.twod_label.setWordWrap(True); r2.addWidget(self.twod_label, 1)
        r2.addWidget(QLabel("Contours from"))
        self.levels = QDoubleSpinBox(); self.levels.setRange(2, 500); self.levels.setValue(30)
        self.levels.setSuffix(" x noise")
        self.levels.setToolTip("Lowest contour level, in multiples of the 2D noise. Lower shows weaker peaks.")
        self.levels.valueChanged.connect(lambda *a: self.draw_plot()); r2.addWidget(self.levels)
        b2 = QPushButton("Choose 2D spectrum..."); b2.clicked.connect(self.choose_2d)
        b2.setToolTip("Pick the 2D experiment (EXPNO folder) recorded on this sample"); r2.addWidget(b2)
        b3 = QPushButton("2D check of all metabolites..."); b3.clicked.connect(self.check_all_2d)
        b3.setToolTip("For every metabolite in this sample: how many of its expected cross peaks the 2D shows")
        r2.addWidget(b3)
        rv.addWidget(self.twod_row); self.twod_row.setVisible(False)
        self.canvas = Canvas(5.6, 4.4)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        rv.addWidget(self.toolbar); rv.addWidget(self.canvas, 1)
        hint = QLabel("Mouse wheel over the plot: zoom in and out. Magnifier: drag a box to zoom. Cross arrows: "
                      "move. House: back to the start.")
        hint.setStyleSheet("color:#4b5563;"); hint.setWordWrap(True); rv.addWidget(hint)
        self.canvas.mpl_connect("scroll_event", self.scroll)
        return right

    def choose(self):
        d = QFileDialog.getExistingDirectory(self, "Choose a MANC-Q results folder", self.out or "")
        if d:
            if not rs.is_results_folder(d):
                warn(self, "That folder does not contain MANC-Q results (concentrations_long.csv)."); return
            self.load(d)

    def load(self, out):
        self.out = out
        self.spec_cache, self.peak_cache = {}, {}
        self.data = rs.load(out)
        self.where.setText(f"Results folder: {out}  (Ctrl+C or right-click copies the selected values)")
        self.fill()
        self.win.add_recent(out)

    def fill_recent(self, folders):
        self.recent.blockSignals(True); self.recent.clear(); self.recent.addItem("Recent results...")
        for f in folders:
            self.recent.addItem(os.path.basename(os.path.normpath(f)) or f, f)
            self.recent.setItemData(self.recent.count() - 1, f, Qt.ToolTipRole)
        self.recent.blockSignals(False)

    def open_recent(self, i):
        f = self.recent.itemData(i)
        self.recent.setCurrentIndex(0)
        if not f:
            return
        if not rs.is_results_folder(f):
            warn(self, f"These results are no longer there:\n{f}"); return
        self.load(f)

    def copy_cells(self, numbers_only):
        """Selected cells as tab-separated text with metabolite and sample names, ready to paste into Excel."""
        if self.data is None:
            return
        idx = self.table.selectedIndexes()
        if not idx:
            return
        rows = sorted({i.row() for i in idx}); cols = sorted({i.column() for i in idx})
        G = self.data["long"].set_index(["metabolite", "sample"])
        allowed = self.FILTERS[self.filter.currentIndex()][1] or rs.TIERS[:4]
        lines = ["metabolite\t" + "\t".join(self.samples[c] for c in cols)]
        for r in rows:
            vals = []
            for c in cols:
                m, smp = self.mets[r], self.samples[c]
                if numbers_only:
                    ok = (m, smp) in G.index and G.loc[(m, smp)].tier in allowed
                    vals.append(f"{G.loc[(m, smp)].conc_mM:.6g}" if ok else "")
                else:
                    it = self.table.item(r, c)
                    vals.append(it.text() if it else "")
            lines.append(self.mets[r] + "\t" + "\t".join(vals))
        QApplication.clipboard().setText("\n".join(lines) + "\n")

    def save_plot(self):
        r, c = self.table.currentRow(), self.table.currentColumn()
        if self.data is None or r < 0 or c < 0:
            return
        name = f"{self.mets[r]}_EXPNO{self.samples[c]}" + ("_2D" if self.dim.checkedId() == 1 else "") + \
            ("_whole" if self.span.checkedId() == 1 else "")
        name = name.replace(" ", "_").replace("/", "-")
        f, _ = QFileDialog.getSaveFileName(self, "Save this plot", os.path.join(self.out, name + ".png"),
                                           "PNG image (*.png);;SVG for editing (*.svg);;PDF (*.pdf)")
        if f:
            self.canvas.fig.savefig(f, dpi=300, bbox_inches="tight")

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
        self.table.setHorizontalHeaderLabels([f"{s}\n{str(titles.get(s, ''))[:24]}" for s in samples])
        for j, s_ in enumerate(samples):
            self.table.horizontalHeaderItem(j).setToolTip(f"EXPNO {s_}: {titles.get(s_, '')}")
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
                tip = f"{m}, EXPNO {s}: {t}, {r.method}"
                if t == "Manually adjusted" and "auto_tier" in r:
                    tip += f"\nAutomatic: {r.auto_tier}, {r.auto_conc_mM:.3g} mM"
                it.setToolTip(tip)
                self.table.setItem(i, j, it)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        if mets and samples:
            self.table.setCurrentCell(0, 0)

    def cell(self, r, c):
        if self.data is None or r < 0 or c < 0 or r >= len(self.mets) or c >= len(self.samples):
            self.cur = None; self.canvas.fig.clear(); self.canvas.draw(); return
        self.cur = (self.mets[r], self.samples[c])
        self.fill_peaks()
        self.draw_plot()

    # ---- data for the plot (read once per sample)
    def spectrum_of(self, s):
        if s not in self.spec_cache:
            f = os.path.join(self.out, f"simulated_spectrum_sample_{s}.csv.gz")
            if len(self.spec_cache) > 4:
                self.spec_cache.pop(next(iter(self.spec_cache)))
            QApplication.setOverrideCursor(Qt.WaitCursor)
            try:
                self.spec_cache[s] = pd.read_csv(f).sort_values("ppm").reset_index(drop=True) if os.path.exists(f) else None
            finally:
                QApplication.restoreOverrideCursor()
        return self.spec_cache[s]

    def peaks_of(self, s):
        if s not in self.peak_cache:
            self.peak_cache[s] = twod.peak_table(self.out, s)
        return self.peak_cache[s]

    def twod_of(self, s):
        """(2D data or None, 2D folder or None, how it was found)."""
        path, how = twod.partner(self.out, s)
        if path is None:
            return None, None, how
        if path not in self.twod_cache:
            if len(self.twod_cache) > 2:
                self.twod_cache.pop(next(iter(self.twod_cache)))
            QApplication.setOverrideCursor(Qt.WaitCursor)
            try:
                self.twod_cache[path] = twod.load(path)
            except Exception as exc:
                self.twod_cache[path] = None
                QApplication.restoreOverrideCursor()
                error_box(self, f"Cannot read the 2D spectrum in {path}: {exc}", traceback.format_exc())
                return None, path, how
            QApplication.restoreOverrideCursor()
        return self.twod_cache[path], path, how

    def multiplets_of(self, m, s):
        """Rows of the peak table for metabolite m in sample s (high ppm first), or None."""
        P = self.peaks_of(s)
        if P is None:
            return None
        sub = P[P.metabolite == m]
        return sub.sort_values("centre_ppm", ascending=False) if len(sub) else None

    def fill_peaks(self):
        m, s = self.cur
        old = self.peak.currentIndex()
        self.peak.blockSignals(True); self.peak.clear()
        self.peak.addItem("Main peak", "main"); self.peak.addItem("All its peaks", "all")
        sub = self.multiplets_of(m, s)
        if sub is not None:
            for r in sub.itertuples():
                self.peak.addItem(f"{r.centre_ppm:.3f} ppm ({r.n_protons:g} H)", (r.from_ppm, r.to_ppm))
        self.peak.setCurrentIndex(old if 0 <= old < 2 else 0)
        self.peak.blockSignals(False)

    def view_changed(self):
        self.twod_row.setVisible(self.dim.checkedId() == 1)
        self.peak.setEnabled(self.span.checkedId() == 0)
        self.draw_plot()

    def peak_range(self, m, s):
        """(low, high) ppm to show for the peak chosen in the Peak list, or None."""
        sub = self.multiplets_of(m, s)
        if sub is None:
            return None
        sel = self.peak.currentData()
        if sel == "all":
            return float(sub.from_ppm.min()) - 0.03, float(sub.to_ppm.max()) + 0.03
        if sel == "main" or sel is None:
            rep = sub[sub.used_as_reporter] if sub.used_as_reporter.any() else sub
            r = rep.sort_values(["height_snr", "n_protons"], ascending=False).iloc[0]
            lo_, hi_ = float(r.from_ppm), float(r.to_ppm)
        else:
            lo_, hi_ = sel
        half = max(0.018, (hi_ - lo_) / 2 + 0.012)
        c_ = (lo_ + hi_) / 2
        return c_ - half, c_ + half

    # ---- drawing
    def draw_plot(self):
        self.canvas.fig.clear()
        self.axes1d = None
        if self.data is None or self.cur is None:
            self.canvas.draw(); return
        m, s = self.cur
        if self.dim.checkedId() == 1:
            self.draw_2d(m, s)
        else:
            self.draw_1d(m, s)
        self.canvas.draw_idle()
        self.toolbar.update()          # the house button goes back to this view

    def title_of(self, m, s):
        G = self.data["long"]
        row = G[(G.metabolite == m) & (G["sample"] == s)]
        if not len(row):
            return f"{m} in EXPNO {s}"
        t = row.iloc[0]
        val = "n/m" if t.tier.startswith("Not measurable") else value_text(t.tier, t.conc_mM) + " mM"
        return f"{m} in EXPNO {s}: {val}, {t.tier}"

    def draw_1d(self, m, s):
        d = self.spectrum_of(s)
        fig = self.canvas.fig
        gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.08)
        ax = fig.add_subplot(gs[0]); ax2 = fig.add_subplot(gs[1], sharex=ax)
        if d is None:
            ax.text(0.5, 0.5, "The fitted spectrum of this sample is missing from the results folder.", ha="center",
                    transform=ax.transAxes)
            tidy(ax); tidy(ax2); return
        x = d.ppm.values
        comp = d[m].values if m in d else np.zeros(len(d))
        gx, gy, gf, gb, gc, gr = q.gapped(x, d.observed.values, d.simulated_total.values, d.baseline.values,
                                          d.baseline.values + comp, d.observed.values - d.simulated_total.values)
        ax.plot(gx, gy, "k", lw=0.7, label="measured")
        ax.plot(gx, gf, color="#dc2626", lw=0.7, label="fit")
        ax.fill_between(gx, gb, gc, color=ACCENT, alpha=0.35, lw=0, label=m)
        ax2.plot(gx, gr, color="#6b7280", lw=0.6); ax2.axhline(0, color="k", lw=0.4)
        whole = self.span.checkedId() == 1
        rng_ = None if whole else self.peak_range(m, s)
        if rng_ is None and not whole:
            ax.text(0.5, 0.92, "No fitted peaks of this metabolite to zoom to: whole spectrum shown.",
                    ha="center", transform=ax.transAxes, fontsize=8)
        sub = self.multiplets_of(m, s)
        if sub is not None and (whole or rng_ is None or self.peak.currentData() == "all"):
            # arrows over the metabolite's peaks, so it can be found in a wide view
            for r in sub.itertuples():
                j = int(np.clip(np.searchsorted(x, r.centre_ppm), 0, len(x) - 1))
                ax.annotate("", (r.centre_ppm, d.observed.values[j]), xytext=(0, 16), textcoords="offset points",
                            arrowprops=dict(arrowstyle="-|>", color=ACCENT, lw=1.2))
        lo, hi = rng_ if rng_ is not None else (float(x.min()), float(x.max()))
        ax.set_title(self.title_of(m, s), fontsize=9, loc="left")
        ax.legend(fontsize=7, frameon=False, loc="upper left"); tidy(ax); ax.tick_params(labelbottom=False)
        tidy(ax2); ax2.set_ylabel("residual", fontsize=8); ax2.set_xlabel("ppm")
        self.axes1d = (ax, ax2, x, d.observed.values, d.simulated_total.values, d.baseline.values)
        ax2.set_xlim(hi, lo)
        self.rescale_1d()
        fig.subplots_adjust(left=0.05, right=0.98, top=0.92, bottom=0.1)

    def rescale_1d(self):
        """Fit the heights to what is visible (wide views leave out the reference peak at 0 ppm)."""
        ax, ax2, x, y, fit, base = self.axes1d
        hi, lo = ax.get_xlim()
        w = (x > min(lo, hi)) & (x < max(lo, hi))
        if abs(hi - lo) > 2 and (w & (np.abs(x) > 0.3)).any():
            w &= np.abs(x) > 0.3
        if not w.any():
            return
        top = max(y[w].max(), fit[w].max(), 1e-9)
        bot = min(0.0, base[w].min(), y[w].min())
        ax.set_ylim(bot - 0.03 * (top - bot), top + 0.15 * (top - bot))
        ax2.set_ylim(-0.3 * top, 0.3 * top)

    def twod_text(self, path, how, d):
        if path is None:
            return ("<b>No 2D spectrum found for this sample</b> (MANC-Q looks for a 2D experiment with the same "
                    "title next to it). Use \"Choose 2D spectrum...\" to pick one.")
        name = os.path.basename(os.path.normpath(path))
        kind = d["kind"] if d else "2D"
        why = {"chosen": "chosen by you", "automatic": "found automatically: same title",
               "automatic, titles differ": "<span style='color:#b91c1c'><b>found automatically, but its title is "
                                           "different: check it is the same sample</b></span>"}.get(how, how)
        extra = f"; moved {d['shift']:+.3f} ppm to put TSP at 0" if d and abs(d["shift"]) > 0.002 else ""
        return f"2D: <b>{kind}, EXPNO {name}</b> ({why}{extra}). Title: {d['title'] if d else ''}"

    def draw_2d(self, m, s):
        fig = self.canvas.fig
        d, path, how = self.twod_of(s)
        self.twod_label.setText(self.twod_text(path, how, d))
        if d is None:
            ax = fig.add_subplot(111)
            ax.text(0.5, 0.5, "No 2D spectrum for this sample.", ha="center", transform=ax.transAxes)
            ax.axis("off"); return
        homo = twod.homonuclear(d["path"])
        sf, off = self.field_of(s, d)
        exp = twod.expected_from_table(self.peaks_of(s), m, sf, d["kind"] == "COSY", off) if homo else []
        chk = twod.check(d, exp)
        # region: the metabolite's protons on both axes; a chosen peak narrows the horizontal (F2) axis to it
        sub = self.multiplets_of(m, s)
        if self.span.checkedId() == 1 or sub is None:
            x0, x1 = -0.3, 9.6
            y0, y1 = (x0, x1) if homo else (float(d["f1c"].min()), float(d["f1c"].max()))
        else:
            sh = [p for ab in exp for p in ab] + list(sub.centre_ppm)
            y0, y1 = min(sh) - 0.2, max(sh) + 0.2
            if y1 - y0 < 0.5:
                c_ = (y0 + y1) / 2; y0, y1 = c_ - 0.25, c_ + 0.25
            if self.peak.currentData() == "all":
                x0, x1 = y0, y1
            else:
                lo, hi = self.peak_range(m, s)
                c_, h_ = (lo + hi) / 2, max(0.15, hi - lo)
                x0, x1 = c_ - h_, c_ + h_
            if not homo:
                y0, y1 = float(d["f1c"].min()), float(d["f1c"].max())
        gs = fig.add_gridspec(2, 1, height_ratios=[1, 3.2], hspace=0.05)
        ax1 = fig.add_subplot(gs[0]); ax = fig.add_subplot(gs[1], sharex=ax1)
        # the 1D with the fitted metabolite on top, for orientation
        sp = self.spectrum_of(s)
        if sp is not None:
            w = ((sp.ppm > x0) & (sp.ppm < x1)).values
            if w.sum() > 3:
                xx, yy = sp.ppm.values[w], sp.observed.values[w]
                comp = sp[m].values[w] if m in sp else np.zeros(w.sum())
                ax1.plot(xx, yy, "k", lw=0.6)
                ax1.fill_between(xx, sp.baseline.values[w], sp.baseline.values[w] + comp, color=ACCENT, alpha=0.4, lw=0)
                ww = np.abs(xx) > 0.3 if (x1 - x0) > 2 else np.ones(len(xx), bool)
                if ww.any():
                    ax1.set_ylim(min(0, yy[ww].min()), yy[ww].max() * 1.08)
        tidy(ax1); ax1.tick_params(labelbottom=False)
        if sub is not None:
            for c_ in sub.centre_ppm:
                ax1.axvline(c_, color=ACCENT, lw=0.5, alpha=0.5)
        # contours of the region (every other point when a large area is drawn)
        c2 = np.where((d["f2c"] >= x0 - 0.05) & (d["f2c"] <= x1 + 0.05))[0]
        c1 = np.where((d["f1c"] >= y0 - 0.05) & (d["f1c"] <= y1 + 0.05))[0]
        step = 2 if c2.size * c1.size > 600000 else 1
        c2, c1 = c2[::step], c1[::step]
        if c2.size > 1 and c1.size > 1:
            z = d["z"][np.ix_(c1, c2)]
            base = self.levels.value() * d["noise"]
            lev = base * 1.6 ** np.arange(16)
            pos = lev[lev < z.max()]
            if len(pos):
                ax.contour(d["f2c"][c2], d["f1c"][c1], z, levels=pos, colors="#111827", linewidths=0.45)
            neg = lev[(lev >= 3 * base) & (lev < (-z).max())]     # only strong negative features (t1 noise)
            if len(neg):
                ax.contour(d["f2c"][c2], d["f1c"][c1], -z, levels=neg, colors="#f87171", linewidths=0.35)
        if homo:
            ax.plot([x0, x1], [x0, x1], color="#9ca3af", lw=0.5, ls="--")
        seen = [(a, b) for a, b, snr in chk if snr >= twod.SEEN_SNR]
        miss = [(a, b) for a, b, snr in chk if snr < twod.SEEN_SNR]
        for pts, col, lab in ((seen, "#15803d", "expected cross peak, seen"), (miss, "#b91c1c", "expected, missing")):
            if pts:
                xs = [p for a, b in pts for p in (a, b)]; ys = [p for a, b in pts for p in (b, a)]
                ax.scatter(xs, ys, s=110, facecolors="none", edgecolors=col, linewidths=1.4, label=lab, zorder=5)
        ax.set_xlim(x1, x0); ax.set_ylim(y1, y0)
        ax.set_xlabel("F2 ppm"); ax.set_ylabel("F1 ppm")
        if seen or miss:
            ax.legend(fontsize=7, loc="lower left", framealpha=0.85)
        txt = twod.summary(chk)[0] if homo else "the cross-peak check is for 1H-1H spectra (TOCSY, COSY) only"
        ax1.set_title(f"{self.title_of(m, s)}\n{d['kind']} EXPNO {os.path.basename(os.path.normpath(path))}: {txt}",
                      fontsize=8.5, loc="left")
        fig.subplots_adjust(left=0.1, right=0.97, top=0.88, bottom=0.1)

    def scroll(self, ev):
        """Mouse wheel: zoom in or out around the pointer (both axes in 2D, ppm only in 1D)."""
        ax = ev.inaxes
        if ax is None or ev.xdata is None:
            return
        f = 0.8 if ev.button == "up" else 1.25
        a, b = ax.get_xlim()
        ax.set_xlim(ev.xdata + (a - ev.xdata) * f, ev.xdata + (b - ev.xdata) * f)
        if self.axes1d is not None:
            self.rescale_1d()
        elif ev.ydata is not None and ax.get_ylabel().startswith("F1"):
            a, b = ax.get_ylim()
            ax.set_ylim(ev.ydata + (a - ev.ydata) * f, ev.ydata + (b - ev.ydata) * f)
        self.canvas.draw_idle()

    # ---- 2D: choose the experiment, check every metabolite
    def choose_2d(self):
        if self.data is None or self.cur is None:
            return
        s = self.cur[1]
        path, _ = twod.partner(self.out, s)
        start = os.path.dirname(path) if path else (twod.sample_dir(self.out, s) or self.out)
        d = QFileDialog.getExistingDirectory(self, f"Choose the 2D experiment (EXPNO folder) for sample {s}", start)
        if not d:
            return
        if not twod.is_2d(d):
            warn(self, "That folder has no processed 2D spectrum (pdata/1/2rr). Process it in TopSpin first, "
                       "then choose it again."); return
        twod.save_link(self.out, s, d)
        self.draw_plot()

    def check_all_2d(self):
        if self.data is None or self.cur is None:
            return
        s = self.cur[1]
        d, path, how = self.twod_of(s)
        if d is None:
            warn(self, "There is no 2D spectrum for this sample. Use \"Choose 2D spectrum...\" first."); return
        if not twod.homonuclear(d["path"]):
            warn(self, "The cross-peak check works for 1H-1H spectra (TOCSY, COSY) only."); return
        sf, off = self.field_of(s, d)
        df = twod.check_sample(d, self.peaks_of(s), self.data["long"], s, sf, off)
        TwoDCheck(self, s, os.path.basename(os.path.normpath(path)), d["kind"], how, df).exec()

    def field_of(self, s, d):
        """(spectrometer MHz, library offset Hz) of a sample, from the QC table."""
        qc = self.data["qc"].set_index("sample")
        sf = float(qc.SF_MHz.get(s, d["sf"])) if "SF_MHz" in qc else d["sf"]
        off = float(qc.GISSMO_offset_Hz.get(s, 0.0)) if "GISSMO_offset_Hz" in qc else 0.0
        return sf, (off if np.isfinite(off) else 0.0)

    def select_metabolite(self, m):
        if self.cur is None:
            return
        if m not in self.mets:
            self.filter.setCurrentIndex(0)
        if m in self.mets:
            self.table.setCurrentCell(self.mets.index(m), self.samples.index(self.cur[1]))

    def review(self, r, c):
        if self.data is None or not (0 <= r < len(self.mets)) or not (0 <= c < len(self.samples)):
            return
        self.win.tabs.setCurrentWidget(self.win.review_tab)
        self.win.review_tab.open(self.out, self.samples[c], self.mets[r])

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
                                                                 "Deconvolution estimate (low confidence)",
                                                                 "Manually adjusted"]
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


class TwoDCheck(QDialog):
    """Every metabolite of one sample against its 2D spectrum: expected cross peaks, and how many are seen."""

    def __init__(self, tab, smp, expno, kind, how, df):
        super().__init__(tab)
        self.tab, self.df, self.smp = tab, df, smp
        self.setWindowTitle(f"2D check, sample {smp}")
        self.resize(1100, 620)
        v = QVBoxLayout(self)
        warn_ = (" <span style='color:#b91c1c'><b>The 2D was found automatically but its title differs from the "
                 "1D: check it is the same sample.</b></span>") if how == "automatic, titles differ" else ""
        lab = QLabel(f"Sample {smp} against its {kind} (EXPNO {expno}).{warn_} For each metabolite MANC-Q predicts "
                     "where its cross peaks should be (protons coupled to each other in the library spin system, at "
                     f"the shifts fitted in the 1D) and looks for signal there: seen = above {twod.SEEN_SNR:g} x the "
                     "2D noise. Seen cross peaks support the assignment; for a strong metabolite, missing ones "
                     "deserve a look. Weak metabolites can simply be below what the 2D shows, and in a crowded "
                     "region another compound can make a cross peak look seen. This check does not change any "
                     "value. Double-click a row to see it.")
        lab.setWordWrap(True); v.addWidget(lab)
        t = QTableWidget(len(df), 5)
        t.setHorizontalHeaderLabels(["Metabolite", "mM", "Tier", "Seen / expected", "Cross peaks (ppm: x noise)"])
        t.verticalHeader().setVisible(False); t.setEditTriggers(QAbstractItemView.NoEditTriggers)
        t.setSelectionBehavior(QAbstractItemView.SelectRows)
        for i, r in enumerate(df.itertuples()):
            frac = r.seen / r.expected if r.expected else None
            items = [SortItem(r.metabolite, r.metabolite.lower()),
                     SortItem(value_text(r.tier, r.conc_mM), float(np.nan_to_num(r.conc_mM))),
                     SortItem(r.tier),
                     SortItem(f"{r.seen} / {r.expected}" if r.expected else "none expected",
                              -1.0 if frac is None else frac, blank=frac is None),
                     SortItem(r.cross_peaks)]
            if frac is not None:
                items[3].setBackground(QBrush(QColor("#dcfce7" if frac >= 0.5 else "#fee2e2")))
            for j, it in enumerate(items):
                t.setItem(i, j, it)
        t.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        t.horizontalHeader().setSectionResizeMode(4, QHeaderView.Stretch)
        t.setSortingEnabled(True)
        t.cellDoubleClicked.connect(lambda r, c: self.tab.select_metabolite(t.item(r, 0).text()))
        v.addWidget(t, 1)
        h = QHBoxLayout()
        bs = QPushButton("Save as CSV in the results folder"); bs.clicked.connect(self.save); h.addWidget(bs)
        h.addStretch()
        bc = QPushButton("Close"); bc.clicked.connect(self.accept); h.addWidget(bc)
        v.addLayout(h)

    def save(self):
        f = os.path.join(self.tab.out, f"twod_check_sample_{self.smp}.csv")
        try:
            self.df.to_csv(f, index=False)
        except OSError as exc:
            error_box(self, f"Could not save {f}: {exc}"); return
        QMessageBox.information(self, "Saved", f"Saved {f}")


# ============================================================================ tab 7: review and adjust
class UpdateWorker(QThread):
    """Rebuilds the results tables and plots after manual adjustments (nothing is refitted)."""
    done = Signal(str)

    def __init__(self, out, samples):
        super().__init__()
        self.out, self.samples = out, samples

    def run(self):
        try:
            from mancq import manual as mn
            mn.update_results(self.out, self.samples, log=lambda *a: None)
            self.done.emit("")
        except Exception as exc:
            self.done.emit(f"{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}")


def value_text(tier, v):
    if tier == q.NOTMEAS:
        return "n/m"
    if tier == "Not detected":
        return "<LOD"
    if not np.isfinite(v):
        return ""
    if tier == "Upper bound only":
        return f"<={v:.2g}"
    if tier.startswith("Deconvolution"):
        return f"~{v:.3g}"
    return f"{v:.3g}"


class ReviewTab(QWidget):
    """Chenomx-style review: each compound drawn over the measured spectrum; its concentration, shift and linewidth
    can be changed by hand (drag the compound, or type values). Adjusted values get the tier 'Manually adjusted'."""
    SHOW = ["Detected compounds", "All compounds", "Adjusted compounds"]

    def __init__(self, win):
        super().__init__()
        self.win = win
        self.out, self.fits, self.fit, self.name, self.dirty = None, {}, None, None, set()
        self.drag, self.worker, self.others = None, None, None
        self.undo_stack, self.redo_stack, self.last_push = {}, {}, (None, 0.0)
        self.names = []
        v = QVBoxLayout(self)
        v.addWidget(heading("7  Review and adjust", "See how each compound was fitted and, where the automatic fit is "
                            "wrong, adjust it by hand as in Chenomx: drag the purple compound up or down to change its "
                            "concentration and sideways to move it, or type the values. Adjusted values get their own "
                            "tier, \"Manually adjusted\", and the automatic value is kept beside them."))
        bar = QHBoxLayout()
        bar.addWidget(QLabel("Sample")); self.sample = QComboBox(); self.sample.setMinimumWidth(160)
        self.sample.currentIndexChanged.connect(self.sample_changed); bar.addWidget(self.sample)
        bar.addWidget(QLabel("Show")); self.show = QComboBox(); self.show.addItems(self.SHOW)
        self.show.currentIndexChanged.connect(self.fill); bar.addWidget(self.show)
        self.search = QLineEdit(); self.search.setPlaceholderText("Find a compound..."); self.search.setMaximumWidth(220)
        self.search.textChanged.connect(self.fill); bar.addWidget(self.search)
        bar.addStretch()
        bo = QPushButton("Open a results folder..."); bo.clicked.connect(self.choose); bar.addWidget(bo)
        v.addLayout(bar)
        self.where = QLabel("Open a results folder (or finish a run) to review it."); self.where.setStyleSheet("color:#4b5563;")
        v.addWidget(self.where)
        split = QSplitter(Qt.Horizontal)
        left = QWidget(); lv = QVBoxLayout(left); lv.setContentsMargins(0, 0, 0, 0)
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["Metabolite", "Automatic", "Now (mM)", "Fit", "Tier now"])
        fit_head = self.table.horizontalHeaderItem(3)
        fit_head.setToolTip("How far the measured spectrum differs from the fit around this compound's peaks, as a "
                            "share of the compound's own signal. Lower is better. Click to sort; click again to "
                            "put the worst fits first.")
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.table.setSortingEnabled(True)
        self.table.sortItems(1, Qt.AscendingOrder)
        self.table.currentCellChanged.connect(lambda r, *a: self.select(r))
        lv.addWidget(self.table, 1)
        sh = QLabel("Click a column heading to sort (again to reverse). Ctrl+Up / Ctrl+Down: previous / next "
                    "compound.")
        sh.setStyleSheet("color:#4b5563;"); sh.setWordWrap(True); lv.addWidget(sh)
        split.addWidget(left)
        right = QWidget(); rv = QVBoxLayout(right); rv.setContentsMargins(0, 0, 0, 0)
        self.canvas = Canvas(6.4, 4.4)
        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        rv.addWidget(self.toolbar); rv.addWidget(self.canvas, 1)
        g = QGroupBox("Adjust the selected compound"); gl = QGridLayout(g)
        gl.addWidget(QLabel("Concentration"), 0, 0)
        self.conc = QDoubleSpinBox(); self.conc.setDecimals(4); self.conc.setRange(0, 100000); self.conc.setSuffix(" mM")
        self.conc.setSingleStep(0.01); gl.addWidget(self.conc, 0, 1)
        gl.addWidget(QLabel("Shift"), 0, 2)
        self.shift = QDoubleSpinBox(); self.shift.setDecimals(2); self.shift.setRange(-60, 60); self.shift.setSuffix(" Hz")
        self.shift.setSingleStep(0.25); gl.addWidget(self.shift, 0, 3)
        gl.addWidget(QLabel("Linewidth"), 0, 4)
        self.width = QDoubleSpinBox(); self.width.setDecimals(2); self.width.setRange(0.3, 5.0); self.width.setPrefix("x ")
        self.width.setSingleStep(0.05); gl.addWidget(self.width, 0, 5)
        gl.addWidget(QLabel("Multiplet"), 1, 0)
        self.mult = QComboBox(); self.mult.currentIndexChanged.connect(self.mult_changed); gl.addWidget(self.mult, 1, 1, 1, 3)
        gl.addWidget(QLabel("Extra shift"), 1, 4)
        self.mshift = QDoubleSpinBox(); self.mshift.setDecimals(2); self.mshift.setRange(-30, 30); self.mshift.setSuffix(" Hz")
        self.mshift.setSingleStep(0.25); self.mshift.setToolTip("Moves only the multiplet chosen on the left")
        gl.addWidget(self.mshift, 1, 5)
        for sp in (self.conc, self.shift, self.width, self.mshift):
            sp.valueChanged.connect(self.spin_changed)
        hb = QHBoxLayout()
        for lab, fn, tip in [("Best fit for this compound", self.best_fit,
                              "Concentration, shift and linewidth that best match the data around this compound's "
                              "peaks, with everything else as it is now"),
                             ("Refit overlapping compounds", self.refit_neighbours,
                              "Keep this compound as it is and refit the compounds that overlap it"),
                             ("Set to zero", self.set_zero, "This compound is not in the sample"),
                             ("Back to automatic", self.reset_one, "Remove the adjustment of this compound")]:
            b = QPushButton(lab); b.setToolTip(tip); b.clicked.connect(fn); hb.addWidget(b)
        gl.addLayout(hb, 2, 0, 1, 6)
        gl.addWidget(QLabel("Note"), 3, 0)
        self.note = QLineEdit(); self.note.setPlaceholderText("Why it was adjusted (saved with the result)")
        self.note.editingFinished.connect(self.note_changed); gl.addWidget(self.note, 3, 1, 1, 5)
        self.adjust_box = g
        rv.addWidget(g)
        split.addWidget(right); split.setSizes([500, 700])
        v.addWidget(split, 1)
        self.status = QLabel(""); self.status.setWordWrap(True); v.addWidget(self.status)
        h = QHBoxLayout()
        br = QPushButton("Undo all adjustments in this sample"); br.clicked.connect(self.reset_sample); h.addWidget(br)
        self.b_undo = QPushButton("Undo"); self.b_undo.setToolTip("Undo the last adjustment (Ctrl+Z)")
        self.b_undo.clicked.connect(self.do_undo); h.addWidget(self.b_undo)
        self.b_redo = QPushButton("Redo"); self.b_redo.setToolTip("Redo (Ctrl+Y)")
        self.b_redo.clicked.connect(self.do_redo); h.addWidget(self.b_redo)
        h.addStretch()
        self.save_btn = primary("Save adjustments and update results"); self.save_btn.clicked.connect(self.save)
        h.addWidget(self.save_btn); v.addLayout(h)
        self.canvas.mpl_connect("button_press_event", self.press)
        self.canvas.mpl_connect("motion_notify_event", self.motion)
        self.canvas.mpl_connect("button_release_event", self.release)
        self.canvas.mpl_connect("scroll_event", self.scroll)
        self.adjust_box.setEnabled(False)
        for keys, fn in ((QKeySequence.Undo, self.do_undo), (QKeySequence("Ctrl+Y"), self.do_redo),
                         (QKeySequence("Ctrl+Shift+Z"), self.do_redo), (QKeySequence("Ctrl+Down"), lambda: self.step(1)),
                         (QKeySequence("Ctrl+Up"), lambda: self.step(-1))):
            sc = QShortcut(keys, self); sc.setContext(Qt.WidgetWithChildrenShortcut); sc.activated.connect(fn)
        self.update_undo_buttons()

    # ---- loading
    def refresh(self):
        out = self.win.results_tab.out
        if out and out != self.out:
            self.open(out)

    def choose(self):
        d = QFileDialog.getExistingDirectory(self, "Choose a MANC-Q results folder", self.out or "")
        if d:
            if not rs.is_results_folder(d):
                warn(self, "That folder does not contain MANC-Q results (concentrations_long.csv)."); return
            self.open(d)

    def confirm_discard(self):
        if not self.dirty:
            return True
        return QMessageBox.question(self, "Unsaved adjustments", "Some adjustments have not been saved. "
                                    "Discard them?") == QMessageBox.Yes

    def open(self, out, sample=None, metabolite=None):
        if out != self.out:
            if not self.confirm_discard():
                return
            self.out, self.fits, self.fit, self.name, self.dirty = out, {}, None, None, set()
            self.undo_stack, self.redo_stack = {}, {}
            data = rs.load(out)
            self.sample.blockSignals(True); self.sample.clear()
            titles = dict(zip(data["qc"]["sample"], data["qc"].get("title", pd.Series(dtype=str)).fillna("")))
            for s in data["samples"]:
                self.sample.addItem(f"{s}  {str(titles.get(s, ''))[:30]}", s)
            self.sample.blockSignals(False)
            self.where.setText(f"Results folder: {out}")
            self.sample_changed()
        if sample is not None:
            i = self.sample.findData(sample)
            if i >= 0 and i != self.sample.currentIndex():
                self.sample.setCurrentIndex(i)
        if metabolite is not None:
            self.goto(metabolite)

    def forget(self):
        """Drop the loaded fits (a new run may replace them)."""
        self.out, self.fits, self.fit, self.name, self.dirty = None, {}, None, None, set()
        self.undo_stack, self.redo_stack = {}, {}
        self.sample.blockSignals(True); self.sample.clear(); self.sample.blockSignals(False)
        self.table.setRowCount(0); self.adjust_box.setEnabled(False); self.status.setText("")
        self.where.setText("Open a results folder (or finish a run) to review it.")
        self.draw()

    def sample_changed(self):
        smp = self.sample.currentData()
        if smp is None:
            return
        if smp not in self.fits:
            QApplication.setOverrideCursor(Qt.WaitCursor)
            try:
                from mancq import manual as mn
                self.fits[smp] = mn.SampleFit(self.out, smp)
            except Exception as exc:
                QApplication.restoreOverrideCursor()
                self.fit = None; self.table.setRowCount(0)
                error_box(self, f"Cannot review sample {smp}: {exc}", traceback.format_exc()); return
            QApplication.restoreOverrideCursor()
        self.fit = self.fits[smp]
        self.name = None
        self.fill()
        self.update_undo_buttons()

    # ---- compound list
    def now(self, name):
        """(tier, value) of a compound with the current adjustments."""
        e = self.fit.edits.get(name)
        if e is not None:
            return q.MANUAL, float(e["conc_mM"])
        return self.fit.auto_value(name)

    def fill(self):
        f = self.fit
        if f is None:
            return
        mode, txt = self.show.currentIndex(), self.search.text().strip().lower()
        names = []
        for n in f.names:
            if n.startswith("Pluronic"):
                continue
            tier, _ = f.auto_value(n)
            if mode == 0 and tier in ("Not detected", q.NOTMEAS) and n not in f.edits:
                continue
            if mode == 2 and n not in f.edits:
                continue
            if txt and txt not in n.lower():
                continue
            names.append(n)
        self.names = names
        total = f.total()
        self.table.blockSignals(True)
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(names))
        for i, n in enumerate(names):
            self.set_row(i, n, total)
        self.table.setSortingEnabled(True)
        self.table.blockSignals(False)
        if names:
            keep = self.name in names
            name = self.name if keep else self.table.item(0, 0).text()
            self.table.blockSignals(True); self.table.setCurrentCell(self.row_of(name), 0); self.table.blockSignals(False)
            if not keep:
                self.select(self.row_of(name))
        else:
            self.name = None; self.adjust_box.setEnabled(False); self.draw()

    def row_of(self, name):
        for r in range(self.table.rowCount()):
            it = self.table.item(r, 0)
            if it is not None and it.text() == name:
                return r
        return -1

    def refresh_list(self):
        """Values and fits of every listed compound after a change (keeps the selection and the zoom)."""
        if self.fit is None:
            return
        total = self.fit.total()
        self.table.blockSignals(True)
        self.table.setSortingEnabled(False)
        for r in range(self.table.rowCount()):
            self.set_row(r, self.table.item(r, 0).text(), total)
        self.table.setSortingEnabled(True)
        if self.name is not None and self.row_of(self.name) >= 0:
            self.table.setCurrentCell(self.row_of(self.name), 0)
            self.table.scrollToItem(self.table.item(self.row_of(self.name), 0))
        self.table.blockSignals(False)

    def set_row(self, i, n, total=None):
        at, av = self.fit.auto_value(n)
        nt, nv = self.now(n)
        order = {t: k for k, t in enumerate(rs.TIERS)}
        mis = self.fit.misfit(n, total)
        it0 = SortItem(n, n.lower()); it0.setToolTip(n)
        it1 = SortItem(value_text(at, av), (order.get(at, 9), -float(np.nan_to_num(av)))); it1.setToolTip(at)
        it1.setBackground(QBrush(QColor(TIER_COL.get(at, "#ffffff")))); it1.setForeground(QBrush(QColor(TIER_TXT.get(at, "#111827"))))
        it2 = SortItem(value_text(nt, nv), float(np.nan_to_num(nv, nan=-1.0)))
        it3 = SortItem("" if mis is None else f"{mis:.0%}", 0.0 if mis is None else mis, blank=mis is None)
        it3.setToolTip("No real signal to judge" if mis is None else
                       f"|measured - fit| around its peaks = {mis:.0%} of its own signal (lower is better)")
        if mis is not None and mis > 0.25:
            it3.setForeground(QBrush(QColor("#b91c1c")))
        it4 = SortItem("Manually adjusted" if nt == q.MANUAL else "", 0 if nt == q.MANUAL else 1)
        if nt == q.MANUAL:
            for it in (it2, it4):
                it.setBackground(QBrush(QColor(TIER_COL[q.MANUAL]))); it.setForeground(QBrush(QColor(TIER_TXT[q.MANUAL])))
        for j, it in enumerate((it0, it1, it2, it3, it4)):
            self.table.setItem(i, j, it)

    def update_row(self):
        """Update the selected compound's row while it is being adjusted."""
        r = self.row_of(self.name) if self.name else -1
        if r < 0:
            return
        self.table.blockSignals(True)
        self.table.setSortingEnabled(False)
        self.set_row(r, self.name)
        self.table.setSortingEnabled(True)
        self.table.setCurrentCell(self.row_of(self.name), 0)
        self.table.scrollToItem(self.table.item(self.row_of(self.name), 0))
        self.table.blockSignals(False)

    def goto(self, name):
        if self.fit is None or name not in self.fit.index:
            return
        if name not in self.names:
            self.show.setCurrentIndex(1); self.search.clear()
        if name in self.names:
            self.table.setCurrentCell(self.row_of(name), 0)

    def step(self, d):
        """Previous / next compound in the list as it is sorted now."""
        r = self.table.currentRow() + d
        if 0 <= r < self.table.rowCount():
            self.table.setCurrentCell(r, 0)

    # ---- selected compound
    def select(self, row):
        if self.fit is None or not (0 <= row < self.table.rowCount()) or self.table.item(row, 0) is None:
            return
        self.name = self.table.item(row, 0).text()
        self.load_controls()
        self.zoom_to()
        self.draw()

    def edit_now(self):
        """The current adjustment of the selected compound, or where one would start."""
        return self.fit.edits.get(self.name) or self.fit.start_edit(self.name)

    def load_controls(self):
        e = self.edit_now()
        self.adjust_box.setEnabled(self.fit.measurable(self.name))
        for sp, val in ((self.conc, e["conc_mM"]), (self.shift, e.get("shift_Hz", 0.0)), (self.width, e.get("width", 1.0))):
            sp.blockSignals(True); sp.setValue(float(val)); sp.blockSignals(False)
        step = max(1e-4, 10 ** np.floor(np.log10(max(e["conc_mM"], 1e-4))) / 10)
        self.conc.setSingleStep(step)
        self.mult.blockSignals(True); self.mult.clear(); self.mult.addItem("All multiplets (whole compound)", 0)
        for j, c, n in self.fit.multiplets(self.name, e):
            self.mult.addItem(f"{j}: {c:.3f} ppm, {n:.1f} H", j)
        self.mult.blockSignals(False)
        self.mult_changed(redraw=False)
        self.note.setText(e.get("note", ""))

    def mult_changed(self, *a, redraw=True):
        j = self.mult.currentData() or 0
        e = self.edit_now()
        self.mshift.blockSignals(True)
        self.mshift.setEnabled(j != 0)
        self.mshift.setValue(float((e.get("mult_shift_Hz") or {}).get(str(j), 0.0)) if j else 0.0)
        self.mshift.blockSignals(False)
        if redraw:
            self.zoom_to(); self.draw()

    def push_undo(self, coalesce=False):
        """Remember the adjustments before a change. Typing or dragging on one compound counts as one change."""
        import copy as _cp
        import time as _t
        key = (self.fit.smp, self.name)
        if coalesce and self.last_push[0] == key and _t.time() - self.last_push[1] < 1.5:
            self.last_push = (key, _t.time()); return
        self.undo_stack.setdefault(self.fit.smp, []).append(_cp.deepcopy(self.fit.edits))
        self.redo_stack[self.fit.smp] = []
        self.last_push = (key, _t.time()) if coalesce else (None, 0.0)
        self.update_undo_buttons()

    def update_undo_buttons(self):
        smp = self.fit.smp if self.fit is not None else None
        self.b_undo.setEnabled(bool(self.undo_stack.get(smp)))
        self.b_redo.setEnabled(bool(self.redo_stack.get(smp)))

    def _swap(self, src, dst):
        import copy as _cp
        if self.fit is None or not src.get(self.fit.smp):
            return
        dst.setdefault(self.fit.smp, []).append(_cp.deepcopy(self.fit.edits))
        self.fit.edits = src[self.fit.smp].pop()
        self.dirty.add(self.fit.smp)
        self.last_push = (None, 0.0)
        self.refresh_list(); self.load_controls(); self.draw(); self.update_undo_buttons()

    def do_undo(self):
        self._swap(self.undo_stack, self.redo_stack)

    def do_redo(self):
        self._swap(self.redo_stack, self.undo_stack)

    def set_edit(self, e, reason=None):
        from mancq import manual as mn
        self.fit.edits[self.name] = mn.stamp(e, reason)
        self.dirty.add(self.fit.smp)
        self.update_row(); self.update_status()

    def spin_changed(self):
        if self.fit is None or self.name is None:
            return
        if self.drag is None:
            self.push_undo(coalesce=True)
        e = dict(self.edit_now())
        e.update(conc_mM=self.conc.value(), shift_Hz=self.shift.value(), width=self.width.value())
        j = self.mult.currentData() or 0
        if j:
            ms = dict(e.get("mult_shift_Hz") or {})
            ms[str(j)] = self.mshift.value()
            e["mult_shift_Hz"] = {k: v for k, v in ms.items() if abs(v) > 1e-9}
        self.set_edit(e, "set by hand")
        self.draw(fast=True)

    def note_changed(self):
        if self.fit is not None and self.name in self.fit.edits:
            e = dict(self.fit.edits[self.name]); e["note"] = self.note.text()
            self.fit.edits[self.name] = e; self.dirty.add(self.fit.smp)

    def best_fit(self):
        if self.name is None:
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            e = self.fit.best_fit(self.name)
        finally:
            QApplication.restoreOverrideCursor()
        e["note"] = self.note.text()
        self.push_undo()
        self.set_edit(e, "best fit of this compound alone")
        self.load_controls(); self.refresh_list(); self.draw()

    def refit_neighbours(self):
        if self.name is None:
            return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            new = self.fit.refit_neighbours(self.name)
        finally:
            QApplication.restoreOverrideCursor()
        if not new:
            QMessageBox.information(self, "Refit", f"No other compound has real signal under {self.name}'s peaks."); return
        lines = []
        for n, e in sorted(new.items()):
            t, v = self.now(n)
            lines.append(f"{n}: {value_text(t, v)} -> {e['conc_mM']:.3g} mM")
        if QMessageBox.question(self, "Refit overlapping compounds",
                                f"Keeping {self.name} as it is, these compounds would change and be marked "
                                f"\"Manually adjusted\":\n\n" + "\n".join(lines) + "\n\nApply?") != QMessageBox.Yes:
            return
        from mancq import manual as mn
        self.push_undo()
        for n, e in new.items():
            self.fit.edits[n] = mn.stamp(e)
        self.dirty.add(self.fit.smp)
        self.refresh_list(); self.draw(); self.update_status()

    def set_zero(self):
        if self.name is not None:
            self.push_undo()
            e = dict(self.edit_now()); e["conc_mM"] = 0.0
            self.set_edit(e, "removed by hand"); self.load_controls(); self.refresh_list(); self.draw()

    def reset_one(self):
        if self.name is not None and self.name in self.fit.edits:
            self.push_undo()
            del self.fit.edits[self.name]
            self.dirty.add(self.fit.smp)
            self.refresh_list(); self.load_controls(); self.draw(); self.update_status()

    def reset_sample(self):
        if self.fit is None or not self.fit.edits:
            return
        if QMessageBox.question(self, "Undo all", f"Remove every adjustment in sample {self.fit.smp}? "
                                "(Undo can bring them back.)") == QMessageBox.Yes:
            self.push_undo()
            self.fit.edits = {}
            self.dirty.add(self.fit.smp)
            self.refresh_list(); self.load_controls(); self.draw(); self.update_status()

    def update_status(self):
        if self.fit is None or self.name is None:
            self.status.setText(""); return
        at, av = self.fit.auto_value(self.name)
        nt, nv = self.now(self.name)
        txt = f"<b>{self.name}</b>, sample {self.fit.smp}: automatic {value_text(at, av)} mM ({at})"
        if nt == q.MANUAL:
            e = self.fit.edits[self.name]
            txt += f"; now <b>{nv:.4g} mM</b>, manually adjusted ({e.get('reason', 'set by hand')})"
        if not self.fit.measurable(self.name):
            txt += ". All its peaks are in ignored regions, so it cannot be adjusted."
        n = sum(len(f.edits) for f in self.fits.values())
        txt += f".   Adjusted compounds in this results folder: {n}" + (" (not saved yet)" if self.dirty else "")
        self.status.setText(txt)

    # ---- plot
    def zoom_to(self):
        if self.fit is None or self.name is None:
            return
        e = self.edit_now()
        ms = self.fit.multiplets(self.name, e)
        j = self.mult.currentData() or 0
        sel = [m for m in ms if m[0] == j] if j else ms
        if not sel:
            return
        if not j:
            # the strongest multiplet first: zoom to it, the whole compound may span several ppm
            contrib = self.fit.current(self.name)
            best = None
            for m in ms:
                w = np.abs(self.fit.x - m[1]) < 0.02
                h = contrib[w].max() if w.any() else 0
                if best is None or h > best[0]:
                    best = (h, m)
            sel = [best[1]] if best[0] > 0 else ms[:1]
        lo = min(m[1] for m in sel) - 0.04
        hi = max(m[1] for m in sel) + 0.04
        self.view = (hi, lo)

    def draw(self, fast=False):
        f = self.fit
        self.canvas.fig.clear()
        if f is None or self.name is None:
            self.canvas.draw(); return
        x = f.x
        e = f.edits.get(self.name)
        others = f.total(without=(self.name,))
        own = f.contribution(self.name, e)
        tot = others + own
        hi, lo = getattr(self, "view", (x.max(), x.min()))
        ax = self.canvas.fig.add_subplot(211); ax2 = self.canvas.fig.add_subplot(212, sharex=ax)
        self.ax = ax
        gx, gy, gt, gb, gc, gr = q.gapped(x, f.yv, tot, f.base, f.base + own, f.yv - tot)
        ax.plot(gx, gy, "k", lw=0.8, label="measured")
        ax.plot(gx, gt, color="#dc2626", lw=0.8, label="total fit")
        ax.fill_between(gx, gb, gc, color=ACCENT, alpha=0.35, lw=0)
        ax.plot(gx, gc, color=ACCENT, lw=0.9, label=self.name)
        if e is not None:
            _, ga = q.gapped(x, f.base + f.contribution(self.name, None))
            ax.plot(gx, ga, color="#6b7280", lw=0.8, ls="--", label="automatic fit of this compound")
        if not fast:
            w = (x > lo) & (x < hi)
            if w.any():
                cand = []
                for n in set(f.names[k] for k in f.auto) | set(f.edits):
                    if n == self.name:
                        continue
                    c = f.current(n)
                    if c[w].max() > 5 * f.sigma:
                        cand.append((c[w].max(), n, c))
                for ci, (_, n, c) in enumerate(sorted(cand, reverse=True)[:6]):
                    col = q.PALETTE[ci % len(q.PALETTE)]
                    _, gcn = q.gapped(x, f.base + c)
                    ax.plot(gx, gcn, color=col, lw=0.6, alpha=0.8)
                    j = np.argmax(np.where(w, c, -np.inf))
                    ax.annotate(n, (x[j], f.base[j] + c[j]), fontsize=7, color=col, ha="center", va="bottom")
        w = (x > lo) & (x < hi)
        top = max(f.yv[w].max(), tot[w].max(), 1e-9) if w.any() else 1.0
        bot = min(0.0, f.base[w].min()) if w.any() else 0.0
        ax.set_ylim(bot - 0.03 * top, top * 1.12)
        ax.legend(fontsize=7, frameon=False, loc="upper left")
        nt, nv = self.now(self.name)
        ax.set_title(f"{self.name} in sample {f.smp}: {value_text(nt, nv)} mM ({nt})", fontsize=9, loc="left")
        tidy(ax)
        ax2.plot(gx, gr, color="#4b5563", lw=0.7); ax2.axhline(0, color="k", lw=0.4)
        ax2.set_ylim(-0.3 * top, 0.3 * top); tidy(ax2); ax2.set_ylabel("residual", fontsize=8)
        ax2.set_xlim(hi, lo); ax2.set_xlabel("ppm")
        self.canvas.fig.tight_layout(); self.canvas.draw_idle()
        if not fast:
            self.update_status()

    # ---- mouse: drag the compound (Chenomx-style), scroll to zoom
    def press(self, ev):
        if (ev.button != 1 or ev.inaxes is not getattr(self, "ax", None) or self.toolbar.mode or self.fit is None
                or self.name is None or not self.fit.measurable(self.name) or ev.xdata is None):
            return
        i = int(np.clip(np.searchsorted(self.fit.x, ev.xdata), 0, len(self.fit.x) - 1))
        e = self.edit_now()
        h = ev.ydata - self.fit.base[i]
        lo, hi = self.ax.get_ylim()
        self.push_undo()
        self.drag = dict(x0=ev.xdata, h0=h if h > 0.05 * (hi - lo) else None, i=i,
                         conc0=float(e["conc_mM"]), shift0=float(e.get("shift_Hz", 0.0)))

    def motion(self, ev):
        d = self.drag
        if d is None or ev.xdata is None or ev.ydata is None:
            return
        conc = d["conc0"]
        if d["h0"] is not None:
            conc = max(0.0, d["conc0"] * (ev.ydata - self.fit.base[d["i"]]) / d["h0"])
        shift = d["shift0"] + (ev.xdata - d["x0"]) * self.fit.sf
        for sp, val in ((self.conc, conc), (self.shift, shift)):
            sp.blockSignals(True); sp.setValue(val); sp.blockSignals(False)
        self.spin_changed()

    def release(self, ev):
        if self.drag is not None:
            self.drag = None
            self.load_controls(); self.refresh_list(); self.draw()

    def scroll(self, ev):
        if ev.inaxes is None or ev.xdata is None:
            return
        f = 0.8 if ev.button == "up" else 1.25
        hi, lo = getattr(self, "view", ev.inaxes.get_xlim())
        c = ev.xdata
        self.view = (c + (hi - c) * f, c + (lo - c) * f)
        self.draw()

    # ---- save
    def save(self):
        if self.out is None:
            return
        if not self.dirty:
            QMessageBox.information(self, "Nothing to save", "There are no new adjustments to save."); return
        from mancq import manual as mn
        for smp in self.dirty:
            if smp in self.fits:
                mn.save_edits(self.out, smp, self.fits[smp].edits)
        self.save_btn.setEnabled(False); self.save_btn.setText("Updating the results...")
        QApplication.setOverrideCursor(Qt.WaitCursor)
        self.worker = UpdateWorker(self.out, sorted(self.dirty))
        self.worker.done.connect(self.saved)
        self.worker.start()

    def saved(self, err):
        QApplication.restoreOverrideCursor()
        self.save_btn.setEnabled(True); self.save_btn.setText("Save adjustments and update results")
        if err:
            error_box(self, "The adjustments were saved, but updating the results failed.", err); return
        self.dirty = set()
        self.update_status()
        self.win.results_tab.load(self.out)
        QMessageBox.information(self, "Saved", "Adjustments saved. The workbook, tables, overlay PDFs and the Results "
                                "step now show them (tier \"Manually adjusted\"; the automatic values are kept in "
                                "the Long_table and Manual_edits sheets).")


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
        self.review_tab = ReviewTab(self)
        for w, lab in [(self.spectra_tab, "1  Spectra"), (self.process_tab, "2  Process spectra"),
                       (self.regions_tab, "3  Regions"), (self.settings_tab, "4  Settings"),
                       (self.run_tab, "5  Run"), (self.results_tab, "6  Results"), (self.review_tab, "7  Review")]:
            self.tabs.addTab(w, lab)
        self.renumber_steps()
        self.tabs.currentChanged.connect(self.tab_changed)
        self.setCentralWidget(self.tabs)
        sb = QStatusBar()
        sb.showMessage(f"MANC-Q {q.__version__}  |  Metabolite Abundance by NMR Curve-fitting, Quantified  |  "
                       "University of Manchester")
        self.setStatusBar(sb)
        self.resize(1180, 760)
        self.prefs = load_prefs()
        self.restore_prefs()
        sys.excepthook = self.unexpected_error

    def unexpected_error(self, etype, value, tb):
        """Errors nobody expected: show them (with a copy button) instead of failing silently."""
        details = "".join(traceback.format_exception(etype, value, tb))
        sys.__stderr__ and sys.__stderr__.write(details)
        if getattr(self, "_in_error", False):
            return
        self._in_error = True
        try:
            QApplication.restoreOverrideCursor()
            error_box(self, f"Something went wrong: {etype.__name__}: {value}\n\nYour work so far is kept. If it "
                            "happens again, click Copy error details and send them to the author.", details)
        finally:
            self._in_error = False

    def restore_prefs(self):
        """Settings remembered from the last session (window_settings.json in the user's AppData folder)."""
        p = self.prefs
        try:
            if p.get("folder"):
                self.spectra_tab.path.setText(p["folder"])
            if p.get("procno"):
                self.spectra_tab.procno.setValue(int(p["procno"]))
            if p.get("settings"):
                self.settings_tab.from_json(p["settings"])
            if p.get("regions"):
                self.regions_tab.from_json(p["regions"])
        except Exception:
            pass
        self.results_tab.fill_recent([f for f in p.get("recent", []) if os.path.isdir(f)])

    def remember(self):
        self.prefs.update(folder=self.spectra_tab.path.text().strip(), procno=self.spectra_tab.procno.value(),
                          settings=self.settings_tab.to_json(), regions=self.regions_tab.to_json())
        save_prefs(self.prefs)

    def add_recent(self, out):
        out = os.path.abspath(out)
        rec = [out] + [f for f in self.prefs.get("recent", []) if os.path.abspath(f) != out]
        self.prefs["recent"] = rec[:10]
        save_prefs(self.prefs)
        self.results_tab.fill_recent([f for f in self.prefs["recent"] if os.path.isdir(f)])

    def spectra_changed(self):
        self.renumber_steps()

    def renumber_steps(self):
        """Number the visible steps 1, 2, 3..."""
        import re
        n = 0
        for i in range(self.tabs.count()):
            if not self.tabs.isTabVisible(i):
                continue
            n += 1
            self.tabs.setTabText(i, re.sub(r"^\d+", str(n), self.tabs.tabText(i)))
            lab = self.tabs.widget(i).findChild(QLabel, "step_heading")
            if lab is not None:
                lab.setText(re.sub(r"^\d+", str(n), lab.text()))

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
        if not self.review_tab.confirm_discard():
            return
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
        spectra = [dict(s) for s in sel]
        processing = self.process_tab.engine_processing()
        settings = dict(REFERENCE=st.ref.currentText().split("-")[0], PROCNO=self.spectra_tab.procno.value(),
                        SAMPLES=[s["expno"] for s in sel], N_WORKERS=st.workers.value(),
                        OVERLAYS=st.overlays.isChecked(), RESUME=True, PROCESSING=processing)
        settings.update(self.regions_tab.engine_settings())
        json.dump(dict(folder=self.folder, settings=st.to_json(), regions=self.regions_tab.to_json(),
                       processing=processing),
                  open(os.path.join(out, "MANC-Q settings used.json"), "w"), indent=2, default=str)
        job = dict(folder=self.folder if os.path.exists(os.path.join(self.folder, "acqus")) is False else
                   os.path.dirname(self.folder.rstrip("/\\")),
                   out=out, ref_mm=st.ref_mm.value(), dilution=single, dilution_file=dil_file,
                   spectra=spectra, settings=settings)
        self.review_tab.forget()
        self.remember()
        self.running = True
        self.tabs.setCurrentWidget(self.run_tab)
        self.run_tab.start(job)

    def closeEvent(self, ev):
        if self.running and QMessageBox.question(self, "Run in progress",
                                                 "A run is in progress. Quit anyway? Finished spectra are kept.") \
                != QMessageBox.Yes:
            ev.ignore(); return
        if not self.review_tab.confirm_discard():
            ev.ignore(); return
        self.remember()
        if self.run_tab.worker is not None and self.run_tab.worker.isRunning():
            self.run_tab.worker.kill()
            self.run_tab.worker.wait(5000)
        ev.accept()


def make_app(argv=None):
    app = QApplication.instance() or QApplication(argv or sys.argv)
    # Always use a light theme: the plots, tier colours and info boxes are designed for a light background, and
    # Windows dark mode would otherwise give white text on the light boxes.
    try:
        app.styleHints().setColorScheme(Qt.ColorScheme.Light)      # Qt 6.8 and later
    except Exception:
        pass
    app.setStyle("Fusion")
    f = QFont(); f.setFamilies(["Segoe UI", "Roboto", "Arial"]); f.setPointSize(10); app.setFont(f)
    pal = QPalette()
    light = {QPalette.Window: "#f3f4f6", QPalette.WindowText: "#111827", QPalette.Base: "#ffffff",
             QPalette.AlternateBase: "#f9fafb", QPalette.ToolTipBase: "#ffffff", QPalette.ToolTipText: "#111827",
             QPalette.PlaceholderText: "#6b7280", QPalette.Text: "#111827", QPalette.Button: "#f9fafb",
             QPalette.ButtonText: "#111827", QPalette.BrightText: "#b91c1c", QPalette.Link: ACCENT,
             QPalette.Highlight: ACCENT, QPalette.HighlightedText: "#ffffff", QPalette.Light: "#ffffff",
             QPalette.Midlight: "#e5e7eb", QPalette.Mid: "#9ca3af", QPalette.Dark: "#6b7280", QPalette.Shadow: "#374151"}
    for role, col in light.items():
        pal.setColor(QPalette.All, role, QColor(col))
    for role in (QPalette.WindowText, QPalette.Text, QPalette.ButtonText):
        pal.setColor(QPalette.Disabled, role, QColor("#6b7280"))
    app.setPalette(pal)
    return app


def main():
    import multiprocessing
    multiprocessing.freeze_support()
    app = make_app()
    w = MainWindow(); w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
