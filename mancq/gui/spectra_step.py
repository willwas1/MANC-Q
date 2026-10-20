"""
mancq.gui.spectra_step - Step 1, Spectra: choose the TopSpin dataset folder and the spectra to quantify.
"""
import os
import numpy as np
from PySide6.QtCore import QThread, Qt, Signal
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QCheckBox, QFileDialog, QHBoxLayout, QHeaderView,
                               QLabel, QLineEdit, QMessageBox, QPushButton, QSpinBox, QSplitter, QTableWidget,
                               QTableWidgetItem, QVBoxLayout, QWidget)
from mancq import fidproc as fp, engine as q
from mancq.gui import theme
from mancq.gui.common import CACHE, Canvas, DEFAULTS, heading, primary, scan_folder, tc, tidy


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
                ppm, y, sf, _ = q.read_spectrum(s["path"], self.procno, S=DEFAULTS)
                t = q.fit_tsp(ppm, y, sf, S=DEFAULTS)
                snr = t["height"] / q.noise_sigma(ppm, y)
                check = "OK"
                if snr < 200:
                    check = "Weak or missing reference peak"
                elif t["fwhm_hz"] > 2.0:
                    check = f"Broad reference line ({t['fwhm_hz']:.1f} Hz)"
                self.row_done.emit(i, float(snr), float(t["fwhm_hz"]), check)
            except Exception as exc:
                self.row_done.emit(i, float("nan"), float("nan"), f"Cannot read: {exc}")


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
        theme.role(hint, "muted"); hint.setWordWrap(True); ch.addWidget(hint, 1)
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
                    it.setForeground(QBrush(QColor(tc("#b45309"))))
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
        it.setForeground(QBrush(QColor(tc("#15803d" if check == "OK" else "#b45309"))))
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
            ppm, y, sf, _ = q.read_spectrum(s["path"], self.procno.value(), S=DEFAULTS)
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
