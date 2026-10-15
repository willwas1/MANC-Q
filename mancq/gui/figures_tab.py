"""The Figures step of the MANC-Q window: publication figures of spectra from a results folder (mancq.figures does
the drawing; this file is only the controls, the preview and the dialogs)."""
import os
import io
import json
import traceback

import numpy as np
import pandas as pd

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QPixmap, QColor, QBrush, QIcon
from PySide6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QPushButton, QLineEdit,
                               QComboBox, QDoubleSpinBox, QSpinBox, QCheckBox, QRadioButton, QButtonGroup, QGroupBox,
                               QListWidget, QListWidgetItem, QPlainTextEdit, QScrollArea, QSplitter, QDialog,
                               QTableWidget, QTableWidgetItem, QHeaderView, QFileDialog, QMessageBox, QApplication,
                               QAbstractItemView, QSizePolicy, QFrame, QColorDialog)

from mancq import figures as F

LABEL_TIERS = [("Quantified", "Quantified"), ("Overlapped", "Overlapped (semi-quantitative)"),
               ("Estimate", "Deconvolution estimate (low confidence)"), ("Manual", "Manually adjusted")]
WHICH = [("Main peak of each metabolite", "main"), ("All its reporter peaks", "reporter"), ("All its peaks", "all")]
STYLES = [("Numbers (key in the caption)", "numbers"), ("Abbreviations", "abbrev"), ("No labels", "off")]


def _A():
    from mancq.gui import app
    return app


class FiguresTab(QWidget):
    def __init__(self, win):
        super().__init__()
        A = _A()
        self.win = win
        self.out, self.R, self.sheet, self.groups = None, None, None, {}
        self.spec = F.default_spec()
        self.loading = False
        self.caption_auto = True
        self.last_info, self.last_fig = None, None
        self.timer = QTimer(self); self.timer.setSingleShot(True); self.timer.timeout.connect(self.redraw)
        v = QVBoxLayout(self)
        v.addWidget(A.heading("8  Figures", "Figures of your spectra for papers and talks: pick a template and the "
                              "samples, check the preview (drawn at the printed size), then Export. Labels come from "
                              "MANC-Q's own fit, and the caption writes itself from the settings."))
        bar = QHBoxLayout()
        bo = QPushButton("Open a results folder..."); bo.clicked.connect(self.choose); bar.addWidget(bo)
        self.where = QLabel("Open a results folder (or finish a run) to make figures from it.")
        self.where.setStyleSheet("color:#4b5563;"); bar.addWidget(self.where, 1)
        v.addLayout(bar)
        split = QSplitter(Qt.Horizontal)
        split.addWidget(self.controls())
        split.addWidget(self.preview_panel())
        split.setSizes([380, 820])
        v.addWidget(split, 1)
        self.set_enabled(False)

    # ------------------------------------------------------------------ controls
    def controls(self):
        box = QWidget(); lay = QVBoxLayout(box); lay.setContentsMargins(0, 0, 6, 0)

        g1 = QGroupBox("1  What to show"); l1 = QGridLayout(g1)
        l1.addWidget(QLabel("Template"), 0, 0)
        self.template = QComboBox()
        for k, name in F.TEMPLATES.items():
            self.template.addItem(name, k)
        self.template.setToolTip("Representative spectrum: one spectrum with labelled peaks. Fit evidence: measured, "
                                 "fit, each chosen metabolite and the residual. Time course: spectra stacked on one "
                                 "scale. Overlay: spectra in colour with a zoomed inset.")
        l1.addWidget(self.template, 0, 1, 1, 2)
        l1.addWidget(QLabel("Samples / groups"), 1, 0, 1, 3)
        self.traces = QListWidget(); self.traces.setMinimumHeight(130)
        self.traces.setToolTip("Tick what to draw. With a sample sheet, replicates of one group and time point are "
                               "one line. Representative and fit evidence use the first ticked line.")
        l1.addWidget(self.traces, 2, 0, 1, 3)
        bs = QPushButton("Samples and groups..."); bs.clicked.connect(self.edit_sheet)
        bs.setToolTip("Say which samples belong together (group or condition, time, replicate): type it in, or "
                      "load a sample sheet (CSV or Excel)")
        ball = QPushButton("All"); ball.clicked.connect(lambda: self.tick_all(True))
        bnone = QPushButton("None"); bnone.clicked.connect(lambda: self.tick_all(False))
        l1.addWidget(bs, 3, 0); l1.addWidget(ball, 3, 1); l1.addWidget(bnone, 3, 2)
        bcol = QPushButton("Colour..."); bcol.clicked.connect(self.pick_colour)
        bcol.setToolTip("Choose the colour of the selected spectrum (or double-click it in the list)")
        bres = QPushButton("Reset colours"); bres.clicked.connect(self.reset_colours)
        bres.setToolTip("Back to the template's own colours (black, or colour-blind safe colours for overlays)")
        l1.addWidget(bcol, 6, 0); l1.addWidget(bres, 6, 1, 1, 2)
        self.traces.itemDoubleClicked.connect(lambda it: self.pick_colour())
        self.mode_mean = QRadioButton("Mean of replicates"); self.mode_one = QRadioButton("One replicate:")
        self.mode_mean.setToolTip("Average the replicate spectra of each group and time point (the caption says so)")
        self.mode_one.setToolTip("Show one real spectrum per group and time point: the replicate named on the right "
                                 "(the first one if a group has no replicate of that name)")
        self.mode_mean.setChecked(True)
        mg = QButtonGroup(self); mg.addButton(self.mode_mean); mg.addButton(self.mode_one)
        self.replicate = QComboBox(); self.replicate.setMinimumWidth(80)
        l1.addWidget(self.mode_mean, 4, 0, 1, 3)
        l1.addWidget(self.mode_one, 5, 0, 1, 2); l1.addWidget(self.replicate, 5, 2)
        lay.addWidget(g1)

        g2 = QGroupBox("2  Region"); l2 = QGridLayout(g2)
        l2.addWidget(QLabel("ppm shown"), 0, 0)
        self.region = QLineEdit(); self.region.setToolTip("For example 0.5-9.5, or several pieces: 0.8-1.6, 3.0-4.2")
        l2.addWidget(self.region, 0, 1, 1, 2)
        l2.addWidget(QLabel("Leave out"), 1, 0)
        self.leave = QLineEdit(); self.leave.setToolTip("Regions cut out of the axis (drawn as breaks), e.g. the water "
                                                        "region 4.6-5.0, or a contaminant: 4.6-5.0, 3.69-3.74")
        l2.addWidget(self.leave, 1, 1, 1, 2)
        self.magnify = QCheckBox("Magnify from"); self.magnify.setToolTip("Enlarge everything above this ppm (the "
                                                                          "weak aromatic region), marked xN")
        self.mag_from = QDoubleSpinBox(); self.mag_from.setRange(0, 12); self.mag_from.setDecimals(2)
        self.mag_from.setSuffix(" ppm")
        self.mag_factor = QSpinBox(); self.mag_factor.setRange(2, 50); self.mag_factor.setPrefix("x")
        l2.addWidget(self.magnify, 2, 0); l2.addWidget(self.mag_from, 2, 1); l2.addWidget(self.mag_factor, 2, 2)
        self.clip = QCheckBox("Cut tall peaks at"); self.clip.setToolTip("Peaks taller than this share of the tallest "
                                                                         "are cut off at the top and marked ~")
        self.clip_pct = QSpinBox(); self.clip_pct.setRange(5, 100); self.clip_pct.setSuffix(" %")
        l2.addWidget(self.clip, 3, 0); l2.addWidget(self.clip_pct, 3, 1)
        self.baseline = QCheckBox("Subtract the fitted baseline")
        self.baseline.setToolTip("Draw the spectra after removing the smooth baseline MANC-Q fitted (said in the "
                                 "caption)")
        l2.addWidget(self.baseline, 4, 0, 1, 3)
        self.show_fit = QCheckBox("Draw the fit as well (orange)")
        l2.addWidget(self.show_fit, 5, 0, 1, 3)
        self.stack_lab = QLabel("Offset upwards")
        self.stack_off = QSpinBox(); self.stack_off.setRange(0, 200); self.stack_off.setSuffix(" % of height")
        self.stack_off.setToolTip("How far each spectrum is raised above the one before, as a share of the plot "
                                  "height (0 = on top of each other)")
        l2.addWidget(self.stack_lab, 6, 0); l2.addWidget(self.stack_off, 6, 1, 1, 2)
        self.side_lab = QLabel("Offset sideways")
        self.side_off = QDoubleSpinBox(); self.side_off.setRange(0, 2); self.side_off.setDecimals(3)
        self.side_off.setSingleStep(0.01); self.side_off.setSuffix(" ppm")
        self.side_off.setToolTip("How far each spectrum is moved to the right of the one before (a staggered, "
                                 "angled stack). The ppm scale then applies to the first spectrum, and the caption "
                                 "says so.")
        l2.addWidget(self.side_lab, 11, 0); l2.addWidget(self.side_off, 11, 1, 1, 2)
        self.inset_lab = QLabel("Inset (ppm)")
        self.inset = QLineEdit(); self.inset.setToolTip("The region drawn enlarged in the inset, e.g. 0.8-1.1. Leave "
                                                        "empty for no inset.")
        l2.addWidget(self.inset_lab, 7, 0); l2.addWidget(self.inset, 7, 1, 1, 2)
        self.fit_lab = QLabel("Metabolites (one row each)")
        self.fit_mets = QListWidget(); self.fit_mets.setMinimumHeight(110)
        self.fit_mets.setToolTip("Tick the metabolites to draw on their own rows, e.g. lactate and threonine")
        self.fit_reg_lab = QLabel("ppm shown")
        self.fit_region = QLineEdit(); self.fit_region.setToolTip("Leave empty to show the chosen metabolites' main "
                                                                  "peaks, or type e.g. 1.25-1.40, 4.05-4.16")
        l2.addWidget(self.fit_lab, 8, 0, 1, 3); l2.addWidget(self.fit_mets, 9, 0, 1, 3)
        l2.addWidget(self.fit_reg_lab, 10, 0); l2.addWidget(self.fit_region, 10, 1, 1, 2)
        lay.addWidget(g2)

        g3 = QGroupBox("3  Labels"); l3 = QGridLayout(g3)
        self.style = QComboBox()
        for lab, k in STYLES:
            self.style.addItem(lab, k)
        l3.addWidget(QLabel("Style"), 0, 0); l3.addWidget(self.style, 0, 1, 1, 3)
        l3.addWidget(QLabel("Label"), 1, 0)
        self.tier_boxes = []
        for i, (lab, t) in enumerate(LABEL_TIERS):
            cb = QCheckBox(lab); cb.setToolTip(f"Label metabolites graded {t}")
            self.tier_boxes.append((cb, t)); l3.addWidget(cb, 1 + i // 2, 1 + i % 2)
        self.which = QComboBox()
        for lab, k in WHICH:
            self.which.addItem(lab, k)
        l3.addWidget(QLabel("Peaks"), 3, 0); l3.addWidget(self.which, 3, 1, 1, 3)
        be = QPushButton("Edit labels..."); be.clicked.connect(self.edit_labels)
        be.setToolTip("Change the text, hide a label, add symbols (* contaminant, ~ cut-off peak, # 13C satellite), "
                      "brackets over multiplets, or your own labels")
        self.n_labels = QLabel("")
        l3.addWidget(be, 4, 0, 1, 2); l3.addWidget(self.n_labels, 4, 2, 1, 2)
        lay.addWidget(g3)

        g4 = QGroupBox("4  Size and export"); l4 = QGridLayout(g4)
        self.width_pre = QComboBox()
        for k, mm in F.WIDTHS_MM.items():
            self.width_pre.addItem(k, mm)
        self.width_pre.addItem("Custom", None)
        self.width_pre.setToolTip("Column widths differ between journals (about 84 to 90 mm single, 174 to 190 mm "
                                  "double): check your journal and use Custom if needed")
        l4.addWidget(QLabel("Width"), 0, 0); l4.addWidget(self.width_pre, 0, 1, 1, 2)
        self.width = QDoubleSpinBox(); self.width.setRange(30, 400); self.width.setSuffix(" mm"); self.width.setDecimals(0)
        self.height = QDoubleSpinBox(); self.height.setRange(20, 400); self.height.setSuffix(" mm"); self.height.setDecimals(0)
        l4.addWidget(QLabel("Size"), 1, 0); l4.addWidget(self.width, 1, 1); l4.addWidget(self.height, 1, 2)
        self.font = QDoubleSpinBox(); self.font.setRange(4, 14); self.font.setSingleStep(0.5); self.font.setSuffix(" pt")
        self.line = QDoubleSpinBox(); self.line.setRange(0.1, 3); self.line.setSingleStep(0.1); self.line.setSuffix(" pt")
        l4.addWidget(QLabel("Text, lines"), 2, 0); l4.addWidget(self.font, 2, 1); l4.addWidget(self.line, 2, 2)
        bx = A_primary("Export...")
        bx.clicked.connect(self.export)
        l4.addWidget(bx, 3, 0, 1, 3)
        lay.addWidget(g4)
        lay.addStretch()
        for w in (self.template, self.style, self.which, self.width_pre, self.replicate):
            w.currentIndexChanged.connect(self.changed)
        for w in (self.region, self.leave, self.inset, self.fit_region):
            w.editingFinished.connect(self.changed)
        for w in (self.mag_from, self.mag_factor, self.clip_pct, self.stack_off, self.side_off, self.width,
                  self.height, self.font, self.line):
            w.valueChanged.connect(self.changed)
        for w in [self.magnify, self.clip, self.baseline, self.show_fit, self.mode_mean] + [b for b, _ in self.tier_boxes]:
            w.toggled.connect(self.changed)
        self.traces.itemChanged.connect(self.changed)
        self.fit_mets.itemChanged.connect(self.changed)
        self.width_pre.currentIndexChanged.connect(self.preset_width)
        scroll = QScrollArea(); scroll.setWidget(box); scroll.setWidgetResizable(True); scroll.setFrameShape(QFrame.NoFrame)
        scroll.setMinimumWidth(360)
        self.ctrl = box
        return scroll

    def preview_panel(self):
        right = QWidget(); v = QVBoxLayout(right); v.setContentsMargins(0, 0, 0, 0)
        self.preview = QLabel(); self.preview.setAlignment(Qt.AlignCenter)
        self.preview.setStyleSheet("background:#e5e7eb;")
        self.preview.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Ignored)
        self.preview.setMinimumHeight(260)
        v.addWidget(self.preview, 3)
        self.check = QLabel(""); self.check.setWordWrap(True); v.addWidget(self.check)
        h = QHBoxLayout()
        h.addWidget(QLabel("Caption (draft)"))
        self.cap_auto = QCheckBox("Rewrite from the settings"); self.cap_auto.setChecked(True)
        self.cap_auto.setToolTip("Untick to keep your own changes to the caption")
        self.cap_auto.toggled.connect(lambda on: on and self.redraw())
        h.addWidget(self.cap_auto); h.addStretch()
        bc = QPushButton("Copy caption"); bc.clicked.connect(lambda: QApplication.clipboard().setText(
            self.caption.toPlainText())); h.addWidget(bc)
        v.addLayout(h)
        self.caption = QPlainTextEdit(); self.caption.setMinimumHeight(90)
        self.caption.textChanged.connect(self.caption_edited)
        v.addWidget(self.caption, 1)
        return right

    def set_enabled(self, on):
        self.ctrl.setEnabled(on)

    # ------------------------------------------------------------------ loading
    def choose(self):
        A = _A()
        from mancq import results as rs
        d = QFileDialog.getExistingDirectory(self, "Choose a MANC-Q results folder", self.out or "")
        if d:
            if not rs.is_results_folder(d):
                A.warn(self, "That folder does not contain MANC-Q results (concentrations_long.csv)."); return
            self.open(d)

    def refresh(self):
        """When the step is shown: follow the results folder open in the Results step."""
        rt = getattr(self.win, "results_tab", None)
        if rt is not None and rt.out and rt.out != self.out:
            self.open(rt.out)

    def settings_file(self):
        return os.path.join(self.out, "figures", "figure_settings.json")

    def open(self, out):
        A = _A()
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            self.R = F.Results(out)
        except Exception as exc:
            QApplication.restoreOverrideCursor()
            A.error_box(self, f"Cannot read the results in {out}: {exc}", traceback.format_exc()); return
        QApplication.restoreOverrideCursor()
        self.out = out
        self.where.setText(f"Results folder: {out}")
        self.sheet = F.load_sheet(out, self.R.samples, self.R.titles)
        spec = F.default_spec(self.R)
        try:
            if os.path.exists(self.settings_file()):
                spec.update(json.load(open(self.settings_file())))
        except (ValueError, OSError):
            pass
        self.spec = spec
        self.fill_lists()
        self.to_controls()
        self.set_enabled(True)
        self.changed()

    def fill_lists(self):
        self.loading = True
        self.groups = F.trace_groups(self.R, self.sheet)
        ticked = set(self.spec.get("traces") or [])
        if not ticked & set(self.groups):
            ticked = set(list(self.groups)[:1])
        self.traces.clear()
        for k, e in self.groups.items():
            n = len(e["samples"])
            it = QListWidgetItem(e["label"] + (f"   ({n} replicates)" if n > 1 else
                                               f"   {self.R.titles.get(e['samples'][0], '')[:30]}" if "EXPNO" in e["label"] else ""))
            it.setData(Qt.UserRole, k)
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked if k in ticked else Qt.Unchecked)
            it.setToolTip("EXPNO " + ", ".join(e["samples"]))
            self.traces.addItem(it)
        reps = sorted({r for e in self.groups.values() if len(e["samples"]) > 1 for r in e["replicates"]})
        self.replicate.clear(); self.replicate.addItems(reps or ["(none)"])
        self.mode_one.setEnabled(bool(reps)); self.replicate.setEnabled(bool(reps))
        G = self.R.long
        good = G[G.tier.isin([t for _, t in LABEL_TIERS])]
        mets = list(good.groupby("metabolite").conc_mM.median().sort_values(ascending=False).index)
        mets = [m for m in mets if not m.startswith("Pluronic")]
        chosen = set(self.spec.get("fit_metabolites") or [])
        self.fit_mets.clear()
        for m in mets:
            it = QListWidgetItem(m); it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked if m in chosen else Qt.Unchecked)
            self.fit_mets.addItem(it)
        self.loading = False

    def to_controls(self):
        s = self.spec
        self.loading = True
        self.template.setCurrentIndex(max(0, self.template.findData(s["template"])))
        (self.mode_mean if s.get("mode", "mean") == "mean" else self.mode_one).setChecked(True)
        if s.get("replicate"):
            i = self.replicate.findText(str(s["replicate"]))
            if i >= 0:
                self.replicate.setCurrentIndex(i)
        self.region.setText(s["region"]); self.leave.setText(s["leave_out"]); self.inset.setText(s["inset_region"])
        self.fit_region.setText(s.get("fit_region", ""))
        self.magnify.setChecked(bool(s["magnify"])); self.mag_from.setValue(s["magnify_from"])
        self.mag_factor.setValue(int(s["magnify_factor"]))
        self.clip.setChecked(bool(s["clip"])); self.clip_pct.setValue(int(s["clip_percent"]))
        self.baseline.setChecked(bool(s.get("baseline_correct", True))); self.show_fit.setChecked(bool(s.get("show_fit")))
        self._tmpl = s["template"]
        self.stack_off.setValue(int(s.get("overlay_offset", 0) if s["template"] == "overlay" else s["stack_offset"]))
        self.side_off.setValue(float(s.get("offset_x_ppm") or 0.0))
        self.style.setCurrentIndex(max(0, self.style.findData(s["label_style"])))
        self.which.setCurrentIndex(max(0, self.which.findData(s.get("label_which", "main"))))
        for cb, t in self.tier_boxes:
            cb.setChecked(t in (s.get("label_tiers") or []))
        i = next((k for k in range(self.width_pre.count()) if self.width_pre.itemData(k) == s["width_mm"]),
                 self.width_pre.count() - 1)
        self.width_pre.setCurrentIndex(i)
        self.width.setValue(s["width_mm"]); self.height.setValue(s["height_mm"])
        self.font.setValue(s["font_pt"]); self.line.setValue(s["line_pt"])
        self.loading = False

    def from_controls(self):
        s = self.spec
        prev = getattr(self, "_tmpl", s["template"])           # the offset box belongs to the template shown
        key = {"stack": "stack_offset", "overlay": "overlay_offset"}.get(prev)
        if key:
            s[key] = self.stack_off.value()
        s["offset_x_ppm"] = self.side_off.value()
        s["template"] = self.template.currentData()
        s["traces"] = [self.traces.item(i).data(Qt.UserRole) for i in range(self.traces.count())
                       if self.traces.item(i).checkState() == Qt.Checked]
        s["mode"] = "mean" if self.mode_mean.isChecked() else "one"
        s["replicate"] = self.replicate.currentText() if self.mode_one.isChecked() else ""
        s["region"], s["leave_out"], s["inset_region"] = self.region.text(), self.leave.text(), self.inset.text()
        s["fit_region"] = self.fit_region.text()
        s["magnify"], s["magnify_from"], s["magnify_factor"] = self.magnify.isChecked(), self.mag_from.value(), \
            self.mag_factor.value()
        s["clip"], s["clip_percent"] = self.clip.isChecked(), self.clip_pct.value()
        s["baseline_correct"], s["show_fit"] = self.baseline.isChecked(), self.show_fit.isChecked()
        s["label_style"], s["label_which"] = self.style.currentData(), self.which.currentData()
        s["label_tiers"] = [t for cb, t in self.tier_boxes if cb.isChecked()]
        s["fit_metabolites"] = [self.fit_mets.item(i).text() for i in range(self.fit_mets.count())
                                if self.fit_mets.item(i).checkState() == Qt.Checked]
        s["width_mm"], s["height_mm"] = self.width.value(), self.height.value()
        s["font_pt"], s["line_pt"] = self.font.value(), self.line.value()

    # ------------------------------------------------------------------ updates
    def preset_width(self):
        mm = self.width_pre.currentData()
        if mm and not self.loading:
            self.width.setValue(mm)

    def changed(self, *a):
        if self.loading or self.R is None:
            return
        old = (self.spec.get("label_tiers"), self.spec.get("label_which"), tuple(self.spec.get("traces") or []),
               self.spec.get("mode"), self.spec.get("replicate"))
        self.from_controls()
        new = (self.spec.get("label_tiers"), self.spec.get("label_which"), tuple(self.spec.get("traces") or []),
               self.spec.get("mode"), self.spec.get("replicate"))
        if self.spec.get("labels") is None or (old != new and not self.spec.get("labels_edited")):
            self.make_labels()
        t = self.spec["template"]
        if t != getattr(self, "_tmpl", t):
            self.loading = True
            self.stack_off.setValue(int(self.spec.get("overlay_offset", 0) if t == "overlay"
                                        else self.spec.get("stack_offset", 25)))
            self.loading = False
        self._tmpl = t
        for w in (self.stack_lab, self.stack_off, self.side_lab, self.side_off):
            w.setVisible(t in ("stack", "overlay"))
        for w in (self.inset_lab, self.inset):
            w.setVisible(t == "overlay")
        for w in (self.fit_lab, self.fit_mets, self.fit_reg_lab, self.fit_region):
            w.setVisible(t == "fit")
        self.show_fit.setVisible(t == "representative")
        for w in (self.magnify, self.mag_from, self.mag_factor):
            w.setEnabled(t in ("representative", "stack"))
        self.width.setEnabled(self.width_pre.currentData() is None)
        self.timer.start(250)

    def make_labels(self):
        traces = self.current_traces()
        samples = [s for tr in traces for s in tr["samples"]] or self.R.samples
        self.spec["labels"] = F.auto_labels(self.R, samples, tuple(self.spec["label_tiers"]), self.spec["label_which"])
        self.spec["labels_edited"] = False

    def current_traces(self):
        return F.collect_traces(self.R, self.spec, self.sheet)

    def fit_default_region(self, traces):
        """Around the main peak of each chosen metabolite, when no region is typed for fit evidence."""
        if self.spec["template"] != "fit" or self.spec.get("fit_region") or not traces or self.R.peaks is None:
            return None
        s = traces[0]["samples"][0]
        P = self.R.peaks[(self.R.peaks["sample"] == s) & self.R.peaks.metabolite.isin(self.spec["fit_metabolites"])]
        if not len(P):
            return None
        rows = P[P.used_as_reporter.astype(bool)] if P.used_as_reporter.astype(bool).any() else P
        parts = [(float(r.from_ppm) - 0.03, float(r.to_ppm) + 0.03) for r in rows.itertuples()]
        parts.sort()
        merged = [list(parts[0])]
        for a, b in parts[1:]:
            if a <= merged[-1][1] + 0.05:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        return ", ".join(f"{a:.2f}-{b:.2f}" for a, b in merged[:3])

    def redraw(self):
        if self.R is None:
            return
        A = _A()
        try:
            traces = self.current_traces()
            spec = dict(self.spec)
            auto = self.fit_default_region(traces)
            if auto:
                spec["fit_region"] = auto
            fig, info = F.render(spec, self.R, traces, dpi=100)
            # preview at the printed size, scaled to the space available
            w_in, h_in = F.figure_size(spec)
            dpi = max(30.0, min(self.preview.width() / w_in, self.preview.height() / h_in, 220.0))
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=dpi, facecolor="white")
            pm = QPixmap(); pm.loadFromData(buf.getvalue())
            self.preview.setPixmap(pm)
            self.last_info, self.last_fig = info, fig
            txt, ok = F.preflight_text(fig, spec)
            warns = " ".join(info.get("warnings") or [])
            self.check.setText(("" if ok else "\u26a0 ") + txt + (f" {warns}" if warns else ""))
            self.check.setStyleSheet("color:#166534;" if ok and not warns else "color:#b45309;")
            self.show_colours(traces)
            n = sum(1 for r in self.spec.get("labels") or [] if r.get("show"))
            self.n_labels.setText(f"{n} label(s){' (edited)' if self.spec.get('labels_edited') else ''}")
            if self.cap_auto.isChecked():
                self.caption.blockSignals(True); self.caption.setPlainText(info.get("caption", ""))
                self.caption.blockSignals(False)
            self.save_settings()
        except Exception as exc:
            A.error_box(self, f"The figure could not be drawn: {exc}", traceback.format_exc())

    def caption_edited(self):
        if self.cap_auto.isChecked():
            self.cap_auto.blockSignals(True); self.cap_auto.setChecked(False); self.cap_auto.blockSignals(False)

    def resizeEvent(self, ev):
        super().resizeEvent(ev)
        if self.R is not None:
            self.timer.start(300)

    def save_settings(self):
        try:
            os.makedirs(os.path.dirname(self.settings_file()), exist_ok=True)
            json.dump(self.spec, open(self.settings_file(), "w"), indent=1, default=str)
        except OSError:
            pass

    def colour_of(self, key, k):
        """The colour a spectrum is drawn in: the one chosen, else the template's own (k = its place among the
        ticked spectra)."""
        own = F.COLOURS[k % len(F.COLOURS)] if self.spec["template"] == "overlay" else "#000000"
        return (self.spec.get("colours") or {}).get(key) or own

    def show_colours(self, traces):
        """A small swatch of each spectrum's colour in the list."""
        order = [tr.get("key") for tr in traces]
        self.traces.blockSignals(True)
        for i in range(self.traces.count()):
            it = self.traces.item(i)
            key = it.data(Qt.UserRole)
            if key in order or key in (self.spec.get("colours") or {}):
                pm = QPixmap(14, 14); pm.fill(QColor(self.colour_of(key, order.index(key) if key in order else 0)))
                it.setIcon(QIcon(pm))
            else:
                it.setIcon(QIcon())
        self.traces.blockSignals(False)

    def pick_colour(self):
        it = self.traces.currentItem()
        if it is None or self.R is None:
            _A().warn(self, "Click a spectrum in the list first, then Colour..."); return
        key = it.data(Qt.UserRole)
        ticked = [self.traces.item(i).data(Qt.UserRole) for i in range(self.traces.count())
                  if self.traces.item(i).checkState() == Qt.Checked]
        start = QColor(self.colour_of(key, ticked.index(key) if key in ticked else 0))
        c = QColorDialog.getColor(start, self, f"Colour of {it.text().split('   ')[0]}")
        if c.isValid():
            self.spec.setdefault("colours", {})[key] = c.name()
            self.changed()

    def reset_colours(self):
        self.spec["colours"] = {}
        self.changed()

    def tick_all(self, on):
        self.traces.blockSignals(True)
        for i in range(self.traces.count()):
            self.traces.item(i).setCheckState(Qt.Checked if on else Qt.Unchecked)
        self.traces.blockSignals(False)
        self.changed()

    # ------------------------------------------------------------------ dialogs
    def edit_sheet(self):
        if self.R is None:
            return
        dlg = SheetDialog(self, self.R, self.sheet)
        if dlg.exec() == QDialog.Accepted:
            self.sheet = dlg.result_sheet()
            F.save_sheet(self.out, self.sheet)
            self.spec["traces"] = []
            self.fill_lists()
            self.changed()

    def edit_labels(self):
        if self.R is None:
            return
        dlg = LabelDialog(self, self.spec.get("labels") or [], self.R)
        if dlg.exec() == QDialog.Accepted:
            self.spec["labels"] = dlg.labels()
            self.spec["labels_edited"] = True
            self.changed()
        elif dlg.reset:
            self.make_labels()
            self.changed()

    def export(self):
        if self.R is None:
            return
        traces = self.current_traces()
        if not traces:
            _A().warn(self, "Tick at least one sample or group to draw first."); return
        spec = dict(self.spec)
        auto = self.fit_default_region(traces)
        if auto:
            spec["fit_region"] = auto
        ExportDialog(self, spec, self.out, self.caption.toPlainText(), self.sheet).exec()


def A_primary(text):
    return _A().primary(text)


# ============================================================================ sample sheet
class SheetDialog(QDialog):
    """Which samples belong together: group (condition), time, replicate and an optional label, typed in or loaded
    from a CSV / Excel sample sheet. Saved in the results folder as sample_sheet.csv."""
    EDIT = ["group", "time_h", "replicate", "label", "notes"]
    HEAD = {"sample": "EXPNO", "nmr_title": "NMR title", "group": "Group / condition", "time_h": "Time (h)",
            "replicate": "Replicate", "label": "Label (optional)", "notes": "Notes"}

    def __init__(self, parent, R, sheet):
        super().__init__(parent)
        self.R = R
        self.setWindowTitle("Samples and groups")
        self.resize(980, 620)
        v = QVBoxLayout(self)
        v.addWidget(QLabel(
            "Replicates of the same <b>group</b> at the same <b>time</b> become one line in a figure (their mean, or "
            "the replicate you pick). Type in the table, paste a column from Excel (select the first cell, Ctrl+V), "
            "or load a sample sheet. A sample sheet is a CSV or Excel file with one row per sample and columns such "
            "as <i>sample</i> (or EXPNO), <i>group</i> (or condition, reactor), <i>time_h</i> (12, 12 h, 2 d and "
            "day 3 all work), <i>replicate</i> and <i>label</i>. Extra columns (feed, VCD, viability, ...) are kept. "
            "Small mistakes are reported, not refused."))
        self.msg = QLabel(""); self.msg.setWordWrap(True); self.msg.setStyleSheet("color:#b45309;")
        v.addWidget(self.msg)
        h = QHBoxLayout()
        bl = QPushButton("Load sample sheet (CSV or Excel)..."); bl.clicked.connect(self.load)
        bt = QPushButton("Save a sheet to fill in..."); bt.clicked.connect(self.save_template)
        bt.setToolTip("Writes a CSV with your samples and their NMR titles, ready to fill in Excel and load back")
        bc = QPushButton("Clear groups"); bc.clicked.connect(self.clear)
        h.addWidget(bl); h.addWidget(bt); h.addWidget(bc); h.addStretch()
        v.addLayout(h)
        self.cols = list(sheet.columns)
        self.table = QTableWidget(0, len(self.cols))
        self.table.setHorizontalHeaderLabels([self.HEAD.get(c, c) for c in self.cols])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.table.installEventFilter(self)
        v.addWidget(self.table, 1)
        self.show_sheet(sheet)
        h2 = QHBoxLayout(); h2.addStretch()
        ok = QPushButton("Use these groups"); ok.clicked.connect(self.accept)
        ca = QPushButton("Cancel"); ca.clicked.connect(self.reject)
        h2.addWidget(ok); h2.addWidget(ca); v.addLayout(h2)

    def show_sheet(self, df):
        self.cols = list(df.columns)
        self.table.setColumnCount(len(self.cols))
        self.table.setHorizontalHeaderLabels([self.HEAD.get(c, c) for c in self.cols])
        self.table.setRowCount(len(df))
        for i, r in enumerate(df.itertuples(index=False)):
            for j, c in enumerate(self.cols):
                val = r[j]
                txt = "" if val is None or (isinstance(val, float) and not np.isfinite(val)) else \
                    (f"{val:g}" if isinstance(val, float) else str(val))
                it = QTableWidgetItem(txt)
                if c in ("sample", "nmr_title"):
                    it.setFlags(it.flags() & ~Qt.ItemIsEditable)
                    it.setBackground(QBrush(QColor("#f3f4f6")))
                self.table.setItem(i, j, it)
        self.table.resizeColumnsToContents()

    def result_sheet(self):
        rows = []
        for i in range(self.table.rowCount()):
            rows.append({c: (self.table.item(i, j).text() if self.table.item(i, j) else "")
                         for j, c in enumerate(self.cols)})
        df = pd.DataFrame(rows, columns=self.cols)
        df["time_h"] = [F.parse_time_h(v) if str(v).strip() else np.nan for v in df["time_h"]]
        df["time_h"] = [np.nan if v is None else v for v in df["time_h"]]
        return df

    def load(self):
        f, _ = QFileDialog.getOpenFileName(self, "Load a sample sheet", "", "Sample sheet (*.csv *.txt *.tsv *.xlsx *.xls)")
        if not f:
            return
        new, warns = F.read_sample_sheet(f, self.R.samples)
        cur = self.result_sheet().set_index("sample")
        for c in new.columns:
            if c != "sample" and c not in cur.columns:
                cur[c] = ""
        n = 0
        for d in new.to_dict("records"):          # (column names such as "Viability (%)" are kept as written)
            s = d.pop("sample")
            if s in cur.index:
                n += 1
                for c, val in d.items():
                    if c == "nmr_title":
                        continue
                    if str(val).strip() not in ("", "nan"):
                        cur.at[s, c] = val
        self.show_sheet(cur.reset_index())
        self.msg.setStyleSheet("color:#166534;" if not warns else "color:#b45309;")
        self.msg.setText(f"Read {os.path.basename(f)}: {n} sample(s) matched. " + " ".join(warns))

    def save_template(self):
        f, _ = QFileDialog.getSaveFileName(self, "Save a sample sheet to fill in", "sample_sheet.csv", "CSV (*.csv)")
        if f:
            self.result_sheet().to_csv(f, index=False)
            self.msg.setStyleSheet("color:#166534;")
            self.msg.setText(f"Saved {f}. Fill in group, time_h and replicate (in Excel, keep it as CSV), then load it.")

    def clear(self):
        for i in range(self.table.rowCount()):
            for j, c in enumerate(self.cols):
                if c in self.EDIT:
                    self.table.setItem(i, j, QTableWidgetItem(""))

    def eventFilter(self, obj, ev):
        """Ctrl+V pastes a block copied from Excel, starting at the selected cell."""
        from PySide6.QtCore import QEvent
        from PySide6.QtGui import QKeySequence
        if obj is self.table and ev.type() == QEvent.KeyPress and ev.matches(QKeySequence.Paste):
            text = QApplication.clipboard().text()
            r0, c0 = max(self.table.currentRow(), 0), max(self.table.currentColumn(), 0)
            for i, line in enumerate(text.rstrip("\n").split("\n")):
                for j, val in enumerate(line.rstrip("\r").split("\t")):
                    r, c = r0 + i, c0 + j
                    if r < self.table.rowCount() and c < self.table.columnCount() and \
                            self.cols[c] not in ("sample", "nmr_title"):
                        self.table.setItem(r, c, QTableWidgetItem(val.strip()))
            return True
        return super().eventFilter(obj, ev)


# ============================================================================ labels
class LabelDialog(QDialog):
    COLS = ["Show", "Metabolite", "ppm", "Text", "Symbol", "Bracket"]

    def __init__(self, parent, labels, R):
        super().__init__(parent)
        self.reset = False
        self.setWindowTitle("Peak labels")
        self.resize(760, 560)
        v = QVBoxLayout(self)
        v.addWidget(QLabel("Labels come from MANC-Q's fit (metabolite at its fitted peak position). Untick to hide "
                           "one, change the text (used when the style is Abbreviations), add a symbol (* contaminant "
                           "or poloxamer, ~ peak cut off, # 13C satellite: each is explained in the caption) or a "
                           "bracket over a multiplet. <b>Add a label</b> puts your own text at any ppm."))
        self.table = QTableWidget(0, len(self.COLS))
        self.table.setHorizontalHeaderLabels(self.COLS)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        v.addWidget(self.table, 1)
        for r in labels:
            self.add_row(r)
        h = QHBoxLayout()
        ba = QPushButton("Add a label"); ba.clicked.connect(lambda: self.add_row(dict(show=True, metabolite="", ppm=1.0,
                                                                                         text="*", symbol="", bracket=False)))
        bd = QPushButton("Remove selected"); bd.clicked.connect(self.remove)
        br = QPushButton("Start again from the fit"); br.clicked.connect(self.do_reset)
        h.addWidget(ba); h.addWidget(bd); h.addWidget(br); h.addStretch()
        ok = QPushButton("OK"); ok.clicked.connect(self.accept)
        ca = QPushButton("Cancel"); ca.clicked.connect(self.reject)
        h.addWidget(ok); h.addWidget(ca)
        v.addLayout(h)

    def add_row(self, r):
        i = self.table.rowCount(); self.table.insertRow(i)
        sh = QTableWidgetItem(); sh.setFlags(sh.flags() | Qt.ItemIsUserCheckable)
        sh.setCheckState(Qt.Checked if r.get("show", True) else Qt.Unchecked)
        self.table.setItem(i, 0, sh)
        m = QTableWidgetItem(r.get("metabolite") or "(your own label)")
        m.setFlags(m.flags() & ~Qt.ItemIsEditable); m.setData(Qt.UserRole, dict(r))
        self.table.setItem(i, 1, m)
        p = QTableWidgetItem(f"{float(r.get('ppm', 0)):.4f}")
        if r.get("metabolite"):
            p.setFlags(p.flags() & ~Qt.ItemIsEditable)
        self.table.setItem(i, 2, p)
        self.table.setItem(i, 3, QTableWidgetItem(r.get("text", "")))
        sym = QComboBox(); sym.addItems(["", "*", "~", "#"]); sym.setCurrentText(r.get("symbol", ""))
        self.table.setCellWidget(i, 4, sym)
        br = QTableWidgetItem(); br.setFlags(br.flags() | Qt.ItemIsUserCheckable)
        br.setCheckState(Qt.Checked if r.get("bracket") else Qt.Unchecked)
        self.table.setItem(i, 5, br)

    def remove(self):
        for i in sorted({ix.row() for ix in self.table.selectedIndexes()}, reverse=True):
            self.table.removeRow(i)

    def do_reset(self):
        self.reset = True
        self.reject()

    def labels(self):
        out = []
        for i in range(self.table.rowCount()):
            r = dict(self.table.item(i, 1).data(Qt.UserRole) or {})
            try:
                ppm = float(self.table.item(i, 2).text())
            except ValueError:
                ppm = float(r.get("ppm", 0))
            if not r.get("metabolite"):
                r["lo"] = r["hi"] = ppm
            r.update(show=self.table.item(i, 0).checkState() == Qt.Checked, ppm=ppm,
                     text=self.table.item(i, 3).text(), symbol=self.table.cellWidget(i, 4).currentText(),
                     bracket=self.table.item(i, 5).checkState() == Qt.Checked)
            out.append(r)
        return out


# ============================================================================ export
class ExportDialog(QDialog):
    FORMATS = [("pdf", "PDF (vector)", True), ("svg", "SVG (vector, editable)", True), ("eps", "EPS (vector)", False),
               ("tiff", "TIFF, 1200 dpi", True), ("png", "PNG, 600 dpi", False)]
    EXTRAS = [("caption", "Caption text (.txt)", True), ("data", "The data behind the figure (.csv)", True),
              ("record", "Processing record (.txt)", True), ("script", "Python script that redraws it", False)]

    def __init__(self, parent, spec, out, caption, sheet):
        super().__init__(parent)
        self.spec, self.out, self.cap, self.sheet = spec, out, caption, sheet
        self.setWindowTitle("Export figure")
        self.resize(560, 420)
        v = QVBoxLayout(self)
        h = QHBoxLayout(); h.addWidget(QLabel("Save as"))
        name = f"figure_{spec['template']}_{int(spec['width_mm'])}mm"
        self.path = QLineEdit(os.path.join(out, "figures", name))
        bb = QPushButton("Browse..."); bb.clicked.connect(self.browse)
        h.addWidget(self.path, 1); h.addWidget(bb); v.addLayout(h)
        v.addWidget(QLabel(f"Printed size {spec['width_mm']:g} x {spec['height_mm']:g} mm; text {spec['font_pt']:g} pt."))
        g = QGroupBox("Formats"); gl = QVBoxLayout(g)
        self.fmt = []
        for k, lab, on in self.FORMATS:
            cb = QCheckBox(lab); cb.setChecked(on); gl.addWidget(cb); self.fmt.append((k, cb))
        v.addWidget(g)
        g2 = QGroupBox("Also save, next to the figure"); gl2 = QVBoxLayout(g2)
        self.ext = []
        for k, lab, on in self.EXTRAS:
            cb = QCheckBox(lab); cb.setChecked(on); gl2.addWidget(cb); self.ext.append((k, cb))
        v.addWidget(g2)
        self.msg = QLabel(""); self.msg.setWordWrap(True); v.addWidget(self.msg)
        h2 = QHBoxLayout(); h2.addStretch()
        self.go = _A().primary("Export"); self.go.clicked.connect(self.do_export)
        self.openb = QPushButton("Open folder"); self.openb.setEnabled(False)
        self.openb.clicked.connect(lambda: _A().open_path(os.path.dirname(self.path.text())))
        cl = QPushButton("Close"); cl.clicked.connect(self.accept)
        h2.addWidget(self.go); h2.addWidget(self.openb); h2.addWidget(cl); v.addLayout(h2)

    def browse(self):
        f, _ = QFileDialog.getSaveFileName(self, "Save the figure as", self.path.text(), "Figure name (*)")
        if f:
            self.path.setText(os.path.splitext(f)[0])

    def do_export(self):
        base = os.path.splitext(self.path.text().strip())[0]
        fmts = [k for k, cb in self.fmt if cb.isChecked()]
        extras = [k for k, cb in self.ext if cb.isChecked()]
        if not base or not (fmts or extras):
            return
        os.makedirs(os.path.dirname(base) or ".", exist_ok=True)
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            files = F.export(self.spec, self.out, base, fmts, extras, caption_text=self.cap, sheet=self.sheet)
        except Exception as exc:
            QApplication.restoreOverrideCursor()
            _A().error_box(self, f"Export failed: {exc}", traceback.format_exc()); return
        QApplication.restoreOverrideCursor()
        self.msg.setStyleSheet("color:#166534;")
        self.msg.setText("Saved:\n" + "\n".join(os.path.basename(f) for f in files) + f"\nin {os.path.dirname(base)}")
        self.openb.setEnabled(True)
