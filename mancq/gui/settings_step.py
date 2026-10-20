"""
mancq.gui.settings_step - Step 4, Settings: the reference concentration, dilutions, results folder and speed.
"""
import datetime
import json
import os
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QApplication, QButtonGroup, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog,
                               QGridLayout, QGroupBox, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMessageBox,
                               QPushButton, QRadioButton, QSpinBox, QTableWidget, QTableWidgetItem, QVBoxLayout,
                               QWidget)
from mancq import engine as q
from mancq.gui import theme
from mancq.gui.common import MW, heading, primary, warn


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
        hd = QHBoxLayout()
        bld = QPushButton("Load from a file..."); bld.clicked.connect(self.load_dilutions)
        bld.setToolTip("A CSV or Excel file: EXPNO (or sample) in the first column, dilution in the second")
        bpd = QPushButton("Paste from Excel"); bpd.clicked.connect(self.paste_dilutions)
        bpd.setToolTip("Copy two columns in Excel (EXPNO, dilution), then click here")
        hd.addWidget(bld); hd.addWidget(bpd); hd.addStretch(); g2l.addLayout(hd)
        self.r_each.toggled.connect(self.dtable.setEnabled); self.dtable.setEnabled(False)
        self.dil.valueChanged.connect(self.sync_dilutions)
        v.addWidget(g2)
        g3 = QGroupBox("Output and speed"); g3l = QGridLayout(g3)
        g3l.addWidget(QLabel("Results folder"), 0, 0)
        self.out = QLineEdit(); g3l.addWidget(self.out, 0, 1)
        bo = QPushButton("Browse..."); bo.clicked.connect(self.browse_out); g3l.addWidget(bo, 0, 2)
        self.out_warn = QLabel(); theme.role(self.out_warn, "warn"); self.out_warn.setWordWrap(True)
        self.out_warn.setVisible(False); g3l.addWidget(self.out_warn, 1, 1, 1, 2)
        self.out.textChanged.connect(self.check_out)
        ncpu = os.cpu_count() or 1
        g3l.addWidget(QLabel("CPU cores to use"), 2, 0)
        self.workers = QSpinBox(); self.workers.setRange(1, max(1, ncpu))
        self.workers.setValue(max(1, ncpu - 1)); g3l.addWidget(self.workers, 2, 1)
        cn = QLabel(f"This PC has {ncpu}. Each core fits one spectrum at a time, so more cores only help with "
                    "several spectra. Use fewer to keep the PC responsive for other work.")
        cn.setWordWrap(True); theme.role(cn, "muted")
        g3l.addWidget(cn, 3, 1, 1, 2)
        self.workers.setToolTip("Number of spectra fitted at the same time (one per CPU core)")
        self.overlays = QCheckBox("Write an overlay PDF for every spectrum (recommended)"); self.overlays.setChecked(True)
        g3l.addWidget(self.overlays, 4, 0, 1, 3)
        self.notify = QCheckBox("Tell me when a run finishes (Windows notification and a sound)")
        self.notify.setChecked(True); g3l.addWidget(self.notify, 5, 0, 1, 3)
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

    def check_out(self):
        """Warn as soon as the chosen results folder already holds MANC-Q results."""
        has = q.has_results(self.out.text().strip())
        self.out_warn.setText(f"{q.EXISTING_MESSAGE} When you click Start you can choose a new subfolder "
                              "inside it instead." if has else "")
        self.out_warn.setVisible(has)

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

    def load_dilutions(self):
        f, _ = QFileDialog.getOpenFileName(self, "Load dilutions", "", "Tables (*.csv *.tsv *.txt *.xlsx *.xls)")
        if f:
            from mancq.figures import _read_table
            try:
                self.apply_dilutions(_read_table(f).values.tolist(), os.path.basename(f))
            except Exception as exc:
                warn(self, f"Could not read {f}:\n{exc}")

    def paste_dilutions(self):
        import re
        text = QApplication.clipboard().text()
        rows = [re.split(r"\t|;|,(?=\s*[\d.]+\s*$)", ln) for ln in text.splitlines() if ln.strip()]
        if not rows:
            warn(self, "The clipboard is empty: copy two columns in Excel (EXPNO and dilution) first."); return
        self.apply_dilutions(rows, "the pasted table")

    def apply_dilutions(self, rows, where):
        """Fill the per-sample table from [EXPNO, dilution] rows (headers and unmatched rows are reported)."""
        from mancq.figures import _sample_id
        known = [self.dtable.item(i, 0).text() for i in range(self.dtable.rowCount()) if self.dtable.item(i, 0)]
        vals, extra = {}, []
        for r in rows:
            if len(r) < 2:
                continue
            s = _sample_id(r[0], set(known))
            try:
                v = float(str(r[1]).strip().replace(",", "."))
            except ValueError:
                continue                          # a header row
            (vals.__setitem__(s, v) if s in known else extra.append(s))
        if not vals:
            warn(self, f"No dilution in {where} matches the selected spectra (EXPNO {', '.join(known)})."); return
        self.r_each.setChecked(True)
        for i, s in enumerate(known):
            if s in vals:
                self.dtable.setItem(i, 2, QTableWidgetItem(f"{vals[s]:g}"))
        missing = [s for s in known if s not in vals]
        msg = f"Dilutions set for {len(vals)} of {len(known)} spectra from {where}."
        if missing:
            msg += (f"\n\nNo dilution for EXPNO {', '.join(missing)}: " +
                    ("it keeps" if len(missing) == 1 else "they keep") + " the value in the table.")
        if extra:
            msg += f"\n\n{where[0].upper() + where[1:]} also has {', '.join(extra)}, " + \
                ("which is not a selected spectrum." if len(extra) == 1 else "which are not selected spectra.")
        QMessageBox.information(self, "Dilutions", msg)

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
            with open(f, "w", encoding="utf-8") as fh:
                json.dump(dict(settings=self.to_json(), regions=self.win.regions_tab.to_json()), fh, indent=2)

    def load(self):
        f, _ = QFileDialog.getOpenFileName(self, "Load settings", "", "Settings (*.json)")
        if f:
            with open(f, encoding="utf-8") as fh:
                d = json.load(fh)
            self.from_json(d.get("settings", {})); self.win.regions_tab.from_json(d.get("regions", {}))
