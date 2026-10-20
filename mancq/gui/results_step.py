"""
mancq.gui.results_step - Step 6, Results: values coloured by tier, the fit behind each value, and the 2D check.
"""
import os
import traceback
import numpy as np
import pandas as pd
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QBrush, QColor, QKeySequence
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QButtonGroup, QComboBox, QDialog, QDoubleSpinBox,
                               QFileDialog, QHBoxLayout, QHeaderView, QLabel, QMessageBox, QPushButton, QSplitter,
                               QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)
from mancq import engine as q, results as rs, twod
from mancq.gui import theme
from mancq.gui.common import (ACCENT, Canvas, SortItem, TIER_COL, TIER_TXT, error_box, heading, open_path, ph_labels,
                              tc, tidy, value_text, warn)


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
        self.where = QLabel(""); theme.role(self.where, "muted"); v.addWidget(self.where)
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
        theme.role(hint, "muted"); hint.setWordWrap(True); rv.addWidget(hint)
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
        try:
            self.data = rs.load(out)
        except rs.ResultsTooNew as exc:
            self.out = None
            warn(self, str(exc)); return
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
        lines = ["metabolite\t" + "\t".join("peak_used_ppm" if c == 0 else self.sample_col(c) for c in cols)]
        for r in rows:
            vals = []
            for c in cols:
                m, smp = self.mets[r], self.sample_col(c)
                if c == 0:
                    it = self.table.item(r, 0)
                    vals.append(it.text() if it else "")
                elif numbers_only:
                    ok = (m, smp) in G.index and G.loc[(m, smp)].tier in allowed
                    ok = ok and np.isfinite(G.loc[(m, smp)].conc_mM)      # upper bounds are blank
                    vals.append(f"{G.loc[(m, smp)].conc_mM:.6g}" if ok else "")
                else:
                    it = self.table.item(r, c)
                    vals.append(it.text() if it else "")
            lines.append(self.mets[r] + "\t" + "\t".join(vals))
        QApplication.clipboard().setText("\n".join(lines) + "\n")

    def save_plot(self):
        if self.data is None or self.cur is None:
            return
        name = f"{self.cur[0]}_EXPNO{self.cur[1]}" + ("_2D" if self.dim.checkedId() == 1 else "") + \
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
        self.table.setRowCount(len(mets)); self.table.setColumnCount(len(samples) + 1)
        ph = ph_labels(self.data["qc"])
        qcd = self.data["qc"]
        qcf = dict(zip(qcd["sample"], zip(qcd["qc_flag"], qcd["qc_reason"].fillna("")))) if "qc_flag" in qcd else {}
        self.table.setHorizontalHeaderLabels(["Peak used\n(ppm)"] +
                                             [f"{s}\n{str(titles.get(s, ''))[:24]}" + (f"\n{ph[s][0]}" if s in ph else "")
                                              + (f"\nQC: {qcf[s][0]}" if qcf.get(s, ("ok",))[0] != "ok" else "")
                                              for s in samples])
        self.table.horizontalHeaderItem(0).setToolTip(
            "Where each metabolite's value is measured: the fitted centre (ppm) of the peak(s) MANC-Q used to "
            "quantify it (its reporter peaks), median over the samples where it was Quantified or Overlapped. "
            "Hover over a value for the peak used in that sample.")
        for j, s_ in enumerate(samples):
            hi = self.table.horizontalHeaderItem(j + 1)
            flag, why = qcf.get(s_, ("ok", ""))
            hi.setToolTip(f"EXPNO {s_}: {titles.get(s_, '')}" + (f"\n{ph[s_][1]}" if s_ in ph else "") +
                          (f"\nQuality: {flag}" + (f": {why}" if why else "") if qcf else ""))
            if flag in ("check", "poor"):
                hi.setForeground(QBrush(QColor(theme.T("warn" if flag == "check" else "bad"))))
        self.table.setVerticalHeaderLabels(mets)
        rep_met, rep_cell = self.data.get("reporter_met", {}), self.data.get("reporter_cell", {})
        for i, m in enumerate(mets):
            peaks_ = [p for p in rep_met.get(m, "").split(", ") if p]
            it = QTableWidgetItem(", ".join(peaks_[:2]) + (f" +{len(peaks_) - 2}" if len(peaks_) > 2 else ""))
            it.setTextAlignment(Qt.AlignCenter)
            it.setForeground(QBrush(QColor(tc("#111827"))))
            it.setToolTip(f"{m}: measured at {rep_met[m]} ppm (reporter peak centre, median over samples)"
                          if rep_met.get(m) else f"{m}: no Quantified or Overlapped value, so no reporter peak")
            self.table.setItem(i, 0, it)
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
                    txt = f"<={q.shown_value(t, r.conc_mM, r.upper_bound_mM):.2g}"
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
                if isinstance(getattr(r, "check_note", ""), str) and r.check_note:
                    tip += f"\nCheck: {r.check_note}"
                if (m, s) in rep_cell:
                    tip += f"\nMeasured at: {', '.join(f'{p:.3f}' for p in rep_cell[(m, s)])} ppm"
                if t == "Manually adjusted" and "auto_tier" in r:
                    tip += f"\nAutomatic: {r.auto_tier}, " + \
                        (f"{r.auto_conc_mM:.3g} mM" if np.isfinite(r.auto_conc_mM) else "upper bound")
                it.setToolTip(tip)
                self.table.setItem(i, j + 1, it)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        if mets and samples:
            self.table.setCurrentCell(0, 1)

    def sample_col(self, c):
        """The sample shown in table column c (column 0 is the peak position), or None."""
        return self.samples[c - 1] if 1 <= c <= len(self.samples) else None

    def cell(self, r, c):
        if self.data is not None and c == 0 and self.samples:
            c = 1 if self.cur is None or self.cur[1] not in self.samples else self.samples.index(self.cur[1]) + 1
        if self.data is None or r < 0 or r >= len(self.mets) or self.sample_col(c) is None:
            self.cur = None; self.canvas.fig.clear(); self.canvas.draw(); return
        self.cur = (self.mets[r], self.sample_col(c))
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
        val = "n/m" if t.tier.startswith("Not measurable") else \
            value_text(t.tier, q.shown_value(t.tier, t.conc_mM, t.upper_bound_mM)) + " mM"
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
               "automatic, titles differ": "<span style='color:" + theme.T("bad") + "'><b>found automatically, but its title is "
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
            self.table.setCurrentCell(self.mets.index(m), self.samples.index(self.cur[1]) + 1)

    def review(self, r, c):
        if c == 0 and self.cur is not None:
            c = self.samples.index(self.cur[1]) + 1
        if self.data is None or not (0 <= r < len(self.mets)) or self.sample_col(c) is None:
            return
        self.win.tabs.setCurrentWidget(self.win.review_tab)
        self.win.review_tab.open(self.out, self.sample_col(c), self.mets[r])

    def open_xlsx(self):
        if self.out:
            open_path(os.path.join(self.out, "metabolite_concentrations.xlsx"))

    def open_folder(self):
        if self.out:
            open_path(self.out)

    def open_overlay(self):
        if self.out and self.cur is not None:
            f = os.path.join(self.out, f"overlay_sample_{self.cur[1]}.pdf")
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
        shown = G.conc_mM.where(G.tier != "Upper bound only", G.upper_bound_mM.fillna(G.conc_mM))
        val = G.assign(v=np.where(G.tier.isin(allowed), shown, np.nan)).pivot(index="sample", columns="metabolite", values="v")
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
        warn_ = (" <span style='color:" + theme.T("bad") + "'><b>The 2D was found automatically but its title differs from the "
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
                     SortItem(value_text(r.tier, q.shown_value(r.tier, r.conc_mM, r.upper_bound_mM)),
                              float(np.nan_to_num(q.shown_value(r.tier, r.conc_mM, r.upper_bound_mM)))),
                     SortItem(r.tier),
                     SortItem(f"{r.seen} / {r.expected}" if r.expected else "none expected",
                              -1.0 if frac is None else frac, blank=frac is None),
                     SortItem(r.cross_peaks)]
            if frac is not None:
                items[3].setBackground(QBrush(QColor("#dcfce7" if frac >= 0.5 else "#fee2e2")))
                items[3].setForeground(QBrush(QColor("#14532d" if frac >= 0.5 else "#7f1d1d")))
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
