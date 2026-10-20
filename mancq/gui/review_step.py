"""
mancq.gui.review_step - Step 7, Review: adjust a compound's fit by hand (Chenomx-style); adjustments are saved apart.
"""
import traceback
import numpy as np
import pandas as pd
from PySide6.QtCore import QThread, Qt, Signal
from PySide6.QtGui import QBrush, QColor, QKeySequence, QShortcut
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QComboBox, QDoubleSpinBox, QFileDialog, QGridLayout,
                               QGroupBox, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMessageBox, QPushButton,
                               QSplitter, QTableWidget, QVBoxLayout, QWidget)
from mancq import engine as q, results as rs
from mancq.gui import theme
from mancq.gui.common import (ACCENT, Canvas, SortItem, TIER_COL, TIER_TXT, error_box, heading, ph_labels, primary,
                              tc, tidy, value_text, warn)


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
        self.where = QLabel("Open a results folder (or finish a run) to review it."); theme.role(self.where, "muted")
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
        theme.role(sh, "muted"); sh.setWordWrap(True); lv.addWidget(sh)
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
            try:
                data = rs.load(out)
            except rs.ResultsTooNew as exc:
                self.out = None
                warn(self, str(exc)); return
            self.sample.blockSignals(True); self.sample.clear()
            titles = dict(zip(data["qc"]["sample"], data["qc"].get("title", pd.Series(dtype=str)).fillna("")))
            ph = ph_labels(data["qc"])
            for s in data["samples"]:
                self.sample.addItem(f"{s}  {str(titles.get(s, ''))[:30]}" + (f"  {ph[s][0]}" if s in ph else ""), s)
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
            except (rs.LibraryChanged, rs.ResultsTooNew) as exc:
                QApplication.restoreOverrideCursor()
                self.fit = None; self.table.setRowCount(0)
                warn(self, f"Cannot review sample {smp}: {exc}"); return
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
            it3.setForeground(QBrush(QColor(tc("#b91c1c"))))
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
