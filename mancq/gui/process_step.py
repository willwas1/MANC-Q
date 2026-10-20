"""
mancq.gui.process_step - Step 2, Process spectra: phase and line broadening of raw FIDs, or changes to TopSpin's processing.
"""
import logging
import traceback
from PySide6.QtCore import Qt
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QCheckBox, QComboBox, QDoubleSpinBox, QGridLayout,
                               QGroupBox, QHBoxLayout, QHeaderView, QLabel, QMessageBox, QPushButton, QSlider,
                               QSplitter, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)
from mancq import fidproc as fp, engine as q
from mancq.gui import theme
_log = logging.getLogger("mancq")      # see mancq.log (window) or the console
from mancq.gui.common import ACCENT, Canvas, DEFAULTS, error_box, heading, primary, tc, tidy, warn


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
        self.asis_note = QLabel(); self.asis_note.setWordWrap(True); theme.role(self.asis_note, "muted")
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
        theme.role(kh, "muted"); kh.setWordWrap(True)
        gl.addWidget(kh, 3, 0, 1, 3)
        self.phase_box = g
        lv.addWidget(g)
        g2 = QGroupBox("Processing"); g2l = QGridLayout(g2)
        g2l.addWidget(QLabel("Line broadening"), 0, 0)
        self.lb = QDoubleSpinBox(); self.lb.setRange(0, 5); self.lb.setSingleStep(0.1); self.lb.setValue(0.3)
        self.lb.setSuffix(" Hz"); self.lb.editingFinished.connect(self.lb_changed); g2l.addWidget(self.lb, 0, 1)
        self.lb_note = QLabel(""); theme.role(self.lb_note, "muted"); g2l.addWidget(self.lb_note, 0, 2)
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
            it.setForeground(QBrush(QColor(tc("#4b5563"))))
        else:
            it = QTableWidgetItem("yes" if s["proc_checked"] else "no")
            it.setForeground(QBrush(QColor(tc("#15803d" if s["proc_checked"] else "#b45309"))))
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
                    _log.debug("TopSpin line broadening not found for EXPNO %s", s["expno"], exc_info=True)
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
                    ppm, y, _, _ = q.read_spectrum(s["path"], self.procno(), S=DEFAULTS)
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
        except Exception:
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
