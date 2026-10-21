"""
mancq.gui.common - what the steps of the window share: colours, helpers (messages, headings, help),
remembered settings, folder scanning, small widgets."""
import logging
import json
import os
import subprocess
import sys
import numpy as np
import pandas as pd
from PySide6.QtCore import QUrl, Qt
from PySide6.QtGui import QIcon, QImage, QPixmap, QTextDocument
from PySide6.QtWidgets import (QApplication, QDialog, QHBoxLayout, QLabel, QMessageBox, QPushButton, QSizePolicy,
                               QTableWidgetItem, QTextBrowser, QVBoxLayout, QWidget)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from mancq import fidproc as fp, engine as q
from mancq.gui import theme
_log = logging.getLogger("mancq")      # see mancq.log (window) or the console


DEFAULTS = q.Settings()        # engine defaults, for reading and checking spectra outside a run (never changed)


ACCENT = "#660099"
TIER_COL, TIER_TXT = theme.TIER_COL, theme.TIER_TXT        # from the one tier colour table (engine.TIER_COLOURS)
MW = {"TSP-d4": 172.27, "DSS-d6": 224.36}
COMMON_REGIONS = [("Urea", 5.70, 5.85), ("DMSO", 2.69, 2.74), ("Methanol", 3.34, 3.38), ("Acetone", 2.21, 2.24)]
POLOX_MODES = ["Model it (recommended)", "Ignore", "Not in my samples"]
CACHE = os.path.join(os.path.expanduser("~"), ".mancq_cache")
# what the window remembers between sessions (plain JSON; delete the file to start afresh)
PREFS_FILE = os.environ.get("MANCQ_PREFS") or os.path.join(
    os.environ.get("APPDATA") or os.path.expanduser("~"), "MANC-Q", "window_settings.json")
HELP_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "help")
ASSETS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets")      # the MANC-Q logo (packaging/make_logo.py)


def logo_pixmap(px):
    """The MANC-Q logo at px pixels (sharp on high-resolution screens)."""
    from PySide6.QtWidgets import QApplication as _QA
    ratio = _QA.instance().devicePixelRatio() if _QA.instance() else 1.0
    want = int(round(px * ratio))
    best = next((n for n in (16, 24, 32, 48, 64, 128, 256) if n >= want), 256)
    pm = QPixmap(os.path.join(ASSETS, f"mancq_logo_{best}.png"))
    if pm.isNull():
        return pm
    pm = pm.scaled(want, want, Qt.KeepAspectRatio, Qt.SmoothTransformation)
    pm.setDevicePixelRatio(ratio)
    return pm


def logo_icon():
    f = os.path.join(ASSETS, "mancq.ico")
    return QIcon(f) if os.path.exists(f) else QIcon()


def tc(colour):
    """A status colour written as a fixed hex value (green, amber, red, grey, purple) -> the matching colour of the
    current theme, so text stays readable in light and dark mode."""
    tok = {"#15803d": "ok", "#166534": "ok", "#b45309": "warn", "#b91c1c": "bad", "#4b5563": "muted",
           "#6b7280": "muted", "#111827": "text", ACCENT: "accent_fg"}.get(str(colour).lower())
    return theme.T(tok) if tok else colour


def _tier_colours():
    """Cells with no value (not detected, not measurable) take the table's own background in either mode."""
    TIER_COL["Not detected"], TIER_TXT["Not detected"] = theme.T("bg"), theme.T("muted")
    TIER_COL["Not measurable (region ignored)"] = theme.T("surface")
    TIER_TXT["Not measurable (region ignored)"] = theme.T("muted")


theme.on_change(_tier_colours)


def _toolbar_icons():
    """The zoom / pan toolbars above the plots pick light or dark icons from the window colour: redraw them."""
    from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
    for w in QApplication.allWidgets():
        if isinstance(w, NavigationToolbar2QT):
            for text, tip, image, callback in w.toolitems:
                if text and callback in getattr(w, "_actions", {}):
                    w._actions[callback].setIcon(w._icon(image + ".png"))


theme.on_change(_toolbar_icons)


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
        with open(PREFS_FILE, encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        _log.warning("the remembered window settings could not be read (%s)", PREFS_FILE, exc_info=True)
        return {}


def save_prefs(d):
    try:
        os.makedirs(os.path.dirname(PREFS_FILE), exist_ok=True)
        with open(PREFS_FILE, "w", encoding="utf-8") as fh:
            json.dump(d, fh, indent=1, default=str)
    except Exception:
        # remembering settings is a convenience, never a reason to fail
        _log.warning("the window settings could not be saved (%s)", PREFS_FILE, exc_info=True)


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
            _log.debug("values of different kinds compared as text", exc_info=True)
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
    with open(f, encoding="utf-8") as fh:
        text = fh.read()
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
    top = QHBoxLayout()
    lg = QLabel(); lg.setPixmap(logo_pixmap(40)); top.addWidget(lg)
    tl = QLabel(f"MANC-Q {q.__version__} help"); tl.setObjectName("step_heading"); top.addWidget(tl); top.addStretch()
    v.addLayout(top)
    tb = QTextBrowser(); tb.setOpenExternalLinks(True); tb.setSearchPaths([HELP_DIR])
    doc = tb.document()
    imgs = os.path.join(HELP_DIR, "images")
    if os.path.isdir(imgs):                   # screenshots scaled to fit the window
        for fn in os.listdir(imgs):
            im = QImage(os.path.join(imgs, fn))
            if not im.isNull():
                doc.addResource(QTextDocument.ImageResource, QUrl(f"images/{fn}"),
                                im.scaledToWidth(780, Qt.SmoothTransformation) if im.width() > 780 else im)
    if not whole:
        tb.setPlainText("The guide was not found. It is in mancq/help/QUICKSTART_GUI.md on the MANC-Q GitHub page.")
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
        s = QLabel(sub); theme.role(s, "muted"); s.setWordWrap(True); v.addWidget(s)
    return box


def fit_to_screen(dlg, w, h):
    """Open a dialog at w x h pixels, or smaller if the screen is (laptops, 125 % or 150 % Windows scaling)."""
    scr = dlg.screen() or QApplication.primaryScreen()
    if scr is not None:
        av = scr.availableGeometry()
        w, h = min(w, int(av.width() * 0.92)), min(h, int(av.height() * 0.88))
    dlg.resize(w, h)


def primary(text):
    b = QPushButton(text)
    b.setObjectName("primary")
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


PH_FLAG_TEXT = {"ok": "H2 and H4 agree", "one_proton_only": "from one histidine peak only",
                "protons_disagree": "H2 and H4 disagree, treat with care", "no_histidine": "no histidine peak found"}


def ph_labels(qc):
    """{sample: (short label, tooltip)} for the estimated pH in the QC table (results from 1.4.0 on)."""
    out = {}
    if qc is None or "pH_est" not in qc:
        return out
    for _, r in qc.iterrows():
        v, flag = r.get("pH_est"), str(r.get("pH_flag", ""))
        short = f"pH ~{v:.2f} (estimate)" if pd.notna(v) else "pH n/a"
        tip = (f"Estimated pH in the NMR tube: {v:.2f} " if pd.notna(v) else "pH could not be estimated: ") + \
            f"({PH_FLAG_TEXT.get(flag, flag)}). From the histidine H2/H4 peak positions; changes between samples " \
            f"are reliable, the absolute value only to about +/-0.2."
        out[str(r["sample"])] = (short, tip)
    return out


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
                with open(tf, encoding="latin-1") as fh:
                    title = " ".join(fh.read().split())[:60]
                break
        try:
            if os.path.exists(os.path.join(d, "acqus")):
                with open(os.path.join(d, "acqus"), encoding="latin-1") as fh:
                    acq = fh.read().splitlines()
                for ln in acq:
                    if ln.startswith("##$BF1="):
                        field = float(ln.split("=")[1])
        except Exception:
            _log.debug("could not read the title or field of %s", d, exc_info=True)
        # processing: TopSpin spectra start "as is" (phase changes are relative to TopSpin's, lb read later from
        # procs); raw FIDs start with automatic phasing
        out.append(dict(expno=os.path.basename(d), path=d, title=title, field=field, processed=proc, fid=fid,
                        selected=True, snr=None, fwhm=None, check="", proc_asis=proc, proc_lb=None,
                        proc_ph0=0.0 if proc else None, proc_ph1=0.0, proc_baseline=False, proc_checked=False))
    return out


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
