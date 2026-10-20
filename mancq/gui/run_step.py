"""
mancq.gui.run_step - Step 5, Run: the fit in a separate process, with progress for each spectrum.
"""
import logging
import os
import time
import numpy as np
from PySide6.QtCore import QThread, QTimer, Signal
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QHBoxLayout, QHeaderView, QLabel, QPlainTextEdit,
                               QProgressBar, QPushButton, QStyle, QSystemTrayIcon, QTableWidget, QTableWidgetItem,
                               QVBoxLayout, QWidget)
from mancq import engine as q, results as rs
_log = logging.getLogger("mancq")      # see mancq.log (window) or the console
from mancq.gui.common import ACCENT, error_box, heading, primary, tc


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


# ============================================================================ tab 5: run
class RunTab(QWidget):
    def __init__(self, win):
        super().__init__()
        self.win = win
        self.worker, self.out, self.t_start, self.rows, self.log_pos = None, None, None, {}, 0
        self.durations, self.n_total, self.n_done = [], 0, 0
        v = QVBoxLayout(self)
        v.addWidget(heading("5  Run", "Each spectrum takes 1 to 5 minutes. You can leave this running. If it is "
                            "stopped, finished spectra are kept: click Start again with the same results folder and "
                            "choose Carry on to fit only what is missing."))
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
        it = QTableWidgetItem(text); it.setForeground(QBrush(QColor(tc(colour)))); self.table.setItem(i, 2, it)
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
            with open(f, encoding="utf-8", errors="replace") as fh:
                txt = fh.read()
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
            _log.debug("the Windows notification could not be shown", exc_info=True)

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
