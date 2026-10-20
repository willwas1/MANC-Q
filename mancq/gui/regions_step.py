"""
mancq.gui.regions_step - Step 3, Regions: parts of the spectrum to leave out of the fit, and Pluronic.
"""
import logging
import os
from PySide6.QtCore import QThread, Qt, Signal
from PySide6.QtWidgets import (QAbstractItemView, QComboBox, QHBoxLayout, QHeaderView, QLabel, QPushButton,
                               QSplitter, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)
from matplotlib.widgets import SpanSelector
from mancq import engine as q
from mancq.gui import theme
_log = logging.getLogger("mancq")      # see mancq.log (window) or the console
from mancq.gui.common import (ACCENT, CACHE, COMMON_REGIONS, Canvas, DEFAULTS, POLOX_MODES, heading, primary, tidy,
                              warn)


class LibraryWorker(QThread):
    """Builds the compound library at a given field (cached on disk) for the Regions tab."""
    ready = Signal(object)

    def __init__(self, sf):
        super().__init__()
        self.sf = sf

    def run(self):
        try:
            lib = q.prepare_library(q.build_library(self.sf, os.path.join(q.HERE, "library"), CACHE, S=DEFAULTS), self.sf, False, S=DEFAULTS)
            self.ready.emit(lib)
        except Exception:
            _log.warning("the compound library could not be built for the Regions step", exc_info=True)
            self.ready.emit(None)


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
        theme.role(self.info, "info")
        self.info.setAlignment(Qt.AlignTop); self.info.setMinimumHeight(150)
        lv.addWidget(self.info)
        note = QLabel("Metabolites whose peaks all fall inside ignored regions are reported as \"n/m\" (not "
                      "measurable), not as zero.")
        note.setWordWrap(True); theme.role(note, "muted"); lv.addWidget(note)
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
                ppm, y, sf, _ = q.read_spectrum(sel[0]["path"], self.win.spectra_tab.procno.value(), S=DEFAULTS)
                t = q.fit_tsp(ppm, y, sf, S=DEFAULTS)
                self.spec = (sel[0]["path"], ppm - t["centre"], y, sel[0]["expno"])
                self.start_library(sf)
            except Exception:
                _log.warning("the spectrum for the Regions step could not be read", exc_info=True)
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
