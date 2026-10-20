"""
MANC-Q desktop interface.

Eight steps, one tab each: Spectra -> Process spectra -> Regions -> Settings -> Run -> Results -> Review (manual,
Chenomx-style adjustment of fitted compounds) -> Figures (publication figures of the spectra, gui/figures_tab.py). The interface is a thin layer over the MANC-Q engine (package mancq):
everything it does can also be done from the command line (processing with --processing processing_used.csv,
adjustments with apply-edits), and both give identical numbers.
"""
import json
import logging
import os
import sys
import traceback
import pandas as pd
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QKeySequence
from PySide6.QtWidgets import (QApplication, QHBoxLayout, QLabel, QMainWindow, QMenu, QMessageBox, QStatusBar,
                               QTabWidget, QWidget)
from mancq import engine as q
from mancq.gui import theme
_log = logging.getLogger("mancq")      # see mancq.log (window) or the console
# every part of the window is importable from mancq.gui.app, as before the split into files
from mancq.gui.common import (DEFAULTS, ACCENT, TIER_COL, TIER_TXT, MW, COMMON_REGIONS, POLOX_MODES, CACHE,
                              PREFS_FILE, HELP_DIR, ASSETS, logo_pixmap, logo_icon, tc, _tier_colours,
                              _toolbar_icons, open_path, load_prefs, save_prefs, SortItem, error_details, error_box,
                              help_sections, show_help, heading, primary, warn, Canvas, tidy, PH_FLAG_TEXT,
                              ph_labels, scan_folder, value_text)
from mancq.gui.spectra_step import InfoWorker, SpectraTab
from mancq.gui.process_step import ProcessTab
from mancq.gui.regions_step import LibraryWorker, RegionsTab
from mancq.gui.settings_step import SettingsTab
from mancq.gui.run_step import RunWorker, RunTab
from mancq.gui.results_step import ResultsTab, TwoDCheck
from mancq.gui.review_step import UpdateWorker, ReviewTab


# ============================================================================ main window
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"MANC-Q {q.__version__}")
        self.spectra, self.folder, self.running = [], "", False
        self.tabs = QTabWidget()
        self.spectra_tab = SpectraTab(self)
        self.process_tab = ProcessTab(self)
        self.regions_tab = RegionsTab(self)
        self.settings_tab = SettingsTab(self)
        self.run_tab = RunTab(self)
        self.results_tab = ResultsTab(self)
        self.review_tab = ReviewTab(self)
        from mancq.gui.figures_tab import FiguresTab
        self.figures_tab = FiguresTab(self)
        for w, lab in [(self.spectra_tab, "1  Spectra"), (self.process_tab, "2  Process spectra"),
                       (self.regions_tab, "3  Regions"), (self.settings_tab, "4  Settings"),
                       (self.run_tab, "5  Run"), (self.results_tab, "6  Results"), (self.review_tab, "7  Review"),
                       (self.figures_tab, "8  Figures")]:
            self.tabs.addTab(w, lab)
        self.renumber_steps()
        self.tabs.currentChanged.connect(self.tab_changed)
        self.setCentralWidget(self.tabs)
        self.make_menu()
        sb = QStatusBar()
        sb.showMessage(f"MANC-Q {q.__version__}  |  Metabolite Abundance by NMR Curve-fitting, Quantified  |  "
                       "University of Manchester")
        self.setStatusBar(sb)
        self.resize(1180, 760)
        self.prefs = load_prefs()
        self.restore_prefs()
        sys.excepthook = self.unexpected_error

    def make_menu(self):
        """File menu in the top left corner; the logo and name at the other end of the menu bar."""
        mb = self.menuBar()
        fm = self.file_menu = QMenu("&File", mb)          # kept by the window, so Python does not delete it
        mb.addMenu(fm)
        a = fm.addAction("Open a results folder...")
        a.triggered.connect(lambda: (self.tabs.setCurrentWidget(self.results_tab), self.results_tab.choose()))
        fm.addSeparator()
        self.dark_action = fm.addAction("Dark mode")
        self.dark_action.setCheckable(True)
        self.dark_action.setChecked(theme.dark())
        self.dark_action.setToolTip("A black window instead of a white one (the plots stay white). Remembered for next time.")
        self.dark_action.toggled.connect(self.set_dark)
        fm.addSeparator()
        fm.addAction("Help (quick start guide)").triggered.connect(lambda: show_help(self, 0))
        fm.addAction("About MANC-Q").triggered.connect(self.about)
        fm.addSeparator()
        qa = fm.addAction("Quit"); qa.setShortcut(QKeySequence.Quit); qa.triggered.connect(self.close)
        brand = QWidget(); bl = QHBoxLayout(brand); bl.setContentsMargins(0, 0, 6, 0); bl.setSpacing(6)
        lg = QLabel(); lg.setPixmap(logo_pixmap(20)); bl.addWidget(lg)
        nm = QLabel("MANC-Q"); nm.setObjectName("brand"); bl.addWidget(nm)
        self.tabs.setCornerWidget(brand, Qt.TopRightCorner)     # the logo and name, level with the steps

    def set_dark(self, on):
        theme.apply(QApplication.instance(), "dark" if on else "light")
        self.prefs["theme"] = "dark" if on else "light"
        save_prefs(self.prefs)
        # tables whose cells were coloured when they were filled: fill them again in the new colours
        cur = getattr(self.results_tab, "cur", None)
        if self.results_tab.data is not None:
            self.results_tab.fill()
            if cur is not None:
                self.results_tab.select_metabolite(cur[0])
        if getattr(self.review_tab, "fit", None) is not None:
            try:
                self.review_tab.fill()
            except Exception:
                _log.warning("the Review step could not be redrawn after the theme changed", exc_info=True)

    def about(self):
        box = QMessageBox(self)
        box.setWindowTitle("About MANC-Q")
        box.setIconPixmap(logo_pixmap(64))
        box.setText(f"<b>MANC-Q {q.__version__}</b><br>Metabolite Abundance by NMR Curve-fitting, Quantified<br>"
                    "Will Smith, University of Manchester<br>"
                    "<a href='https://github.com/willwas1/MANC-Q'>github.com/willwas1/MANC-Q</a>")
        box.exec()

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
            _log.warning("the remembered window settings could not be applied", exc_info=True)
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
        existing = "stop"
        if q.has_results(out):
            choice = self.ask_existing(out, [s["expno"] for s in sel])
            if choice is None:
                return
            if choice == "new":
                out = q.new_run_folder(out)
            elif choice == "overwrite":
                try:
                    q.clear_results(out)
                except OSError as exc:
                    warn(self, "Could not delete the old results (is a file from this folder open in Excel or "
                               f"a PDF viewer?):\n{exc}"); return
                st.check_out()
            else:
                existing = "resume"
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
        settings["existing"] = existing
        settings.update(self.regions_tab.engine_settings())
        with open(os.path.join(out, "MANC-Q settings used.json"), "w", encoding="utf-8") as fh:
            json.dump(dict(folder=self.folder, settings=st.to_json(), regions=self.regions_tab.to_json(),
                           processing=processing), fh, indent=2, default=str)
        job = dict(folder=self.folder if os.path.exists(os.path.join(self.folder, "acqus")) is False else
                   os.path.dirname(self.folder.rstrip("/\\")),
                   out=out, ref_mm=st.ref_mm.value(), dilution=single, dilution_file=dil_file,
                   spectra=spectra, settings=settings)
        self.review_tab.forget()
        self.remember()
        self.running = True
        self.tabs.setCurrentWidget(self.run_tab)
        self.run_tab.start(job)

    def ask_existing(self, out, expnos):
        """The results folder already has results: "new" (a subfolder), "overwrite", "resume" or None (cancel).
        Carry on is only offered when the run in the folder was stopped or crashed with the same spectra."""
        prev = q.previous_run(out)
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Warning); box.setWindowTitle("Results folder")
        box.setText(q.EXISTING_MESSAGE)
        info = (f"New subfolder: run into {os.path.basename(q.new_run_folder(out))} inside this folder and keep "
                "the old results.\nOverwrite: delete the old results, including manual adjustments and spectra "
                "processed by MANC-Q, then run.")
        resume = None
        if prev and prev["run_id"] and not prev["complete"] and \
                sorted(map(str, prev["samples"] or [])) == sorted(map(str, expnos)):
            resume = box.addButton("Carry on the stopped run", QMessageBox.AcceptRole)
            info += ("\nCarry on: the run in this folder did not finish. Spectra it finished are kept (if no "
                     "setting has changed) and only the rest are fitted.")
        box.setInformativeText(info)
        new = box.addButton("New subfolder", QMessageBox.AcceptRole)
        over = box.addButton("Overwrite", QMessageBox.DestructiveRole)
        box.addButton(QMessageBox.Cancel)
        box.setDefaultButton(new)
        box.exec()
        c = box.clickedButton()
        if c is new:
            return "new"
        if c is over:
            return "overwrite"
        return "resume" if resume is not None and c is resume else None

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
    app.setStyle("Fusion")
    f = QFont(); f.setFamilies(["Segoe UI", "Roboto", "Arial"]); f.setPointSize(10); app.setFont(f)
    if os.name == "nt":                     # its own taskbar icon, not Python's
        try:
            import ctypes
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("MANC-Q.MANC-Q")
        except Exception:
            _log.debug("Windows taskbar identity not set", exc_info=True)
    app.setWindowIcon(logo_icon())
    # light (white) unless dark mode was chosen in File > Dark mode; purple and gold accents (gui/theme.py); the
    # plots stay white in both modes
    theme.apply(app, "dark" if load_prefs().get("theme") == "dark" else "light")
    return app


def log_to_file():
    """Warnings, with the error behind them, go to MANC-Q/mancq.log next to the window settings (at most 1 MB)."""
    import logging.handlers
    try:
        f = os.path.join(os.path.dirname(PREFS_FILE), "mancq.log")
        os.makedirs(os.path.dirname(f), exist_ok=True)
        h = logging.handlers.RotatingFileHandler(f, maxBytes=1_000_000, backupCount=1, encoding="utf-8")
        h.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
        _log.addHandler(h)
        _log.setLevel(logging.INFO)
    except OSError:
        pass                                  # no log file is no reason not to start


def main():
    import multiprocessing
    multiprocessing.freeze_support()
    log_to_file()
    app = make_app()
    w = MainWindow(); w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
