"""MANC-Q window theme: white or near-black background following the Windows light / dark setting, Manchester
purple as the accent, Manchester gold for highlights (the selected step, the status line). One set of colour
tokens per mode; widgets take their colours from here (by a "role" property in the style sheet, or from T()), so
the whole window changes when Windows switches between light and dark. Plots stay white in both modes, like
paper. Every text colour keeps a contrast of at least 4.5 : 1 against its background; gold is used only for
lines and fills, never for text on white (it would be too faint)."""
import os

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QPalette

PURPLE, GOLD = "#660099", "#FFCC33"          # Manchester purple and gold

LIGHT = dict(bg="#FFFFFF", surface="#F6F3F9", surface2="#ECE6F2", text="#17111D", muted="#574F61",
             line="#DCD4E5", accent=PURPLE, accent_hover="#7A1FB0", accent_fg=PURPLE, accent_text="#FFFFFF",
             gold=GOLD, ok="#166534", warn="#8A4B00", bad="#B91C1C", info_bg="#F5EEFB", info_line="#D6BFEA",
             disabled="#8C8496", sel=PURPLE, sel_text="#FFFFFF", edge="#857C91")
DARK = dict(bg="#0B0A0D", surface="#16131A", surface2="#211C27", text="#F2EEF6", muted="#B7AEC2",
            line="#352D3F", accent="#7A1FB0", accent_hover="#8E33C6", accent_fg="#C99DEB", accent_text="#FFFFFF",
            gold=GOLD, ok="#4ADE80", warn="#FBBF24", bad="#F87171", info_bg="#1D1326", info_line="#4B2F66",
            disabled="#6F6779", sel="#7A1FB0", sel_text="#FFFFFF", edge="#8C8398")

_current = dict(LIGHT)
_dark = False
_listeners = []


def is_dark(app=None):
    """True when Windows (or the desktop) is set to dark mode. The environment variable MANCQ_THEME=dark or
    MANCQ_THEME=light overrides it."""
    force = os.environ.get("MANCQ_THEME", "").strip().lower()
    if force in ("dark", "light"):
        return force == "dark"
    try:
        from PySide6.QtWidgets import QApplication
        app = app or QApplication.instance()
        return app.styleHints().colorScheme() == Qt.ColorScheme.Dark
    except Exception:
        return False


def T(name):
    """The colour token `name` for the current mode (e.g. T("muted"), T("ok"))."""
    return _current[name]


def dark():
    return _dark


def role(widget, name):
    """Give a widget a colour role from the style sheet (muted, ok, warn, bad, info, frame); follows the theme."""
    widget.setProperty("role", name)
    st = widget.style()
    st.unpolish(widget)
    st.polish(widget)
    return widget


def on_change(fn):
    """Call fn() after the theme changes (for colours set directly on table cells and the like)."""
    _listeners.append(fn)


def _tick_file():
    """A white tick (PNG, drawn here so no image plug-in is needed) for ticked boxes, kept with MANC-Q's other
    cached files."""
    from PySide6.QtGui import QImage, QPainter, QPen
    from PySide6.QtCore import QPointF
    d = os.path.join(os.path.expanduser("~"), ".mancq_cache")
    f = os.path.join(d, "tick_white_14.png")
    if not os.path.exists(f):
        try:
            os.makedirs(d, exist_ok=True)
            img = QImage(14, 14, QImage.Format_ARGB32)
            img.fill(Qt.transparent)
            p = QPainter(img)
            p.setRenderHint(QPainter.Antialiasing)
            pen = QPen(QColor("#ffffff"), 2.0)
            pen.setCapStyle(Qt.RoundCap); pen.setJoinStyle(Qt.RoundJoin)
            p.setPen(pen)
            p.drawPolyline([QPointF(3.0, 7.2), QPointF(5.8, 10.0), QPointF(11.0, 4.0)])
            p.end()
            if not img.save(f, "PNG"):
                return None
        except Exception:
            return None
    return f.replace(os.sep, "/")


def stylesheet(c):
    tick = _tick_file()
    tick_css = f"image: url({tick});" if tick else ""
    return f"""
QCheckBox::indicator, QListView::indicator, QTableView::indicator, QGroupBox::indicator {{
    width: 14px; height: 14px; border: 1px solid {c['edge']}; border-radius: 3px; background: {c['bg']}; }}
QCheckBox::indicator:hover, QListView::indicator:hover {{ border-color: {c['accent_fg']}; }}
QCheckBox::indicator:checked, QListView::indicator:checked, QTableView::indicator:checked,
QGroupBox::indicator:checked {{ background: {c['accent']}; border-color: {c['accent']}; {tick_css} }}
QCheckBox::indicator:disabled {{ border-color: {c['line']}; background: {c['surface']}; }}
QCheckBox::indicator:checked:disabled {{ background: {c['disabled']}; border-color: {c['disabled']}; }}
QRadioButton::indicator {{ width: 14px; height: 14px; border: 1px solid {c['edge']}; border-radius: 8px;
                          background: {c['bg']}; }}
QRadioButton::indicator:hover {{ border-color: {c['accent_fg']}; }}
QRadioButton::indicator:checked {{ border-color: {c['accent']};
    background: qradialgradient(cx:0.5, cy:0.5, radius:0.5, fx:0.5, fy:0.5, stop:0 #ffffff, stop:0.32 #ffffff,
                                stop:0.42 {c['accent']}, stop:1 {c['accent']}); }}
QRadioButton::indicator:disabled {{ border-color: {c['line']}; background: {c['surface']}; }}
QMainWindow, QDialog {{ background: {c['bg']}; }}
QLabel#step_heading {{ font-size: 18px; font-weight: 700; color: {c['text']}; }}
QLabel[role="muted"] {{ color: {c['muted']}; }}
QLabel[role="ok"] {{ color: {c['ok']}; }}
QLabel[role="warn"] {{ color: {c['warn']}; }}
QLabel[role="bad"] {{ color: {c['bad']}; }}
QLabel[role="info"] {{ background: {c['info_bg']}; color: {c['text']}; border: 1px solid {c['info_line']};
                      border-radius: 8px; padding: 8px; }}
QLabel[role="frame"] {{ background: {c['surface2']}; border-radius: 8px; }}
QPushButton {{ background: {c['surface']}; color: {c['text']}; border: 1px solid {c['line']}; border-radius: 6px;
              padding: 5px 12px; }}
QPushButton:hover {{ border-color: {c['accent_fg']}; }}
QPushButton:pressed {{ background: {c['surface2']}; }}
QPushButton:checked {{ background: {c['accent']}; color: {c['accent_text']}; border-color: {c['accent']}; }}
QPushButton:disabled {{ color: {c['disabled']}; }}
QPushButton#primary {{ background: {c['accent']}; color: {c['accent_text']}; border: 1px solid {c['accent']};
                      padding: 6px 18px; font-weight: 600; }}
QPushButton#primary:hover {{ background: {c['accent_hover']}; border-color: {c['gold']}; }}
QPushButton#primary:disabled {{ background: {c['surface2']}; color: {c['disabled']}; border-color: {c['line']}; }}
QTabWidget::pane {{ border: none; border-top: 1px solid {c['line']}; }}
QTabBar::tab {{ background: transparent; color: {c['muted']}; padding: 8px 14px; border: none;
               border-bottom: 3px solid transparent; }}
QTabBar::tab:hover {{ color: {c['text']}; }}
QTabBar::tab:selected {{ color: {c['accent_fg']}; font-weight: 600; border-bottom: 3px solid {c['gold']}; }}
QGroupBox {{ border: 1px solid {c['line']}; border-radius: 8px; margin-top: 16px; padding: 10px 6px 6px 6px; }}
QGroupBox::title {{ subcontrol-origin: margin; left: 10px; padding: 0 4px; color: {c['accent_fg']}; font-weight: 600; }}
QLineEdit, QPlainTextEdit, QTextEdit {{ background: {c['bg']}; color: {c['text']}; border: 1px solid {c['line']};
                                       border-radius: 5px; padding: 3px 6px; }}
QLineEdit:focus, QPlainTextEdit:focus, QTextEdit:focus {{ border: 1px solid {c['accent_fg']}; }}
QTableWidget, QTableView, QListWidget {{ background: {c['bg']}; color: {c['text']}; gridline-color: {c['line']};
                                       border: 1px solid {c['line']}; border-radius: 6px; }}
QHeaderView::section {{ background: {c['surface']}; color: {c['text']}; border: none;
                       border-right: 1px solid {c['line']}; border-bottom: 1px solid {c['line']}; padding: 4px; }}
QTableCornerButton::section {{ background: {c['surface']}; border: none; }}
QProgressBar {{ border: 1px solid {c['line']}; border-radius: 5px; background: {c['surface']}; text-align: center;
               color: {c['text']}; }}
QProgressBar::chunk {{ background: {c['accent']}; border-radius: 4px; }}
QStatusBar {{ background: {c['surface']}; color: {c['muted']}; border-top: 2px solid {c['gold']}; }}
QToolTip {{ background: {c['surface2']}; color: {c['text']}; border: 1px solid {c['line']}; padding: 4px; }}
QSplitter::handle {{ background: {c['line']}; }}
"""


def palette(c):
    pal = QPalette()
    roles = {QPalette.Window: c["bg"], QPalette.WindowText: c["text"], QPalette.Base: c["bg"],
             QPalette.AlternateBase: c["surface"], QPalette.ToolTipBase: c["surface2"], QPalette.ToolTipText: c["text"],
             QPalette.PlaceholderText: c["muted"], QPalette.Text: c["text"], QPalette.Button: c["surface"],
             QPalette.ButtonText: c["text"], QPalette.BrightText: c["bad"], QPalette.Link: c["accent_fg"],
             QPalette.LinkVisited: c["accent_fg"], QPalette.Highlight: c["sel"], QPalette.HighlightedText: c["sel_text"],
             QPalette.Light: c["surface2"], QPalette.Midlight: c["line"], QPalette.Mid: c["disabled"],
             QPalette.Dark: c["muted"], QPalette.Shadow: "#000000"}
    for r, col in roles.items():
        pal.setColor(QPalette.All, r, QColor(col))
    for r in (QPalette.WindowText, QPalette.Text, QPalette.ButtonText):
        pal.setColor(QPalette.Disabled, r, QColor(c["disabled"]))
    return pal


def apply(app):
    """Set the palette and style sheet for the current Windows mode, and follow later changes."""
    global _dark
    _dark = is_dark(app)
    _current.clear(); _current.update(DARK if _dark else LIGHT)
    app.setPalette(palette(_current))
    app.setStyleSheet(stylesheet(_current))
    for fn in list(_listeners):
        try:
            fn()
        except Exception:
            pass
    if not getattr(app, "_mancq_theme_hooked", False):
        try:
            app.styleHints().colorSchemeChanged.connect(lambda *a: apply(app))
            app._mancq_theme_hooked = True
        except Exception:
            pass
