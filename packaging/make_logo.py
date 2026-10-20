"""Draws the MANC-Q logo and writes every file made from it. Run from the repository folder:

    python packaging/make_logo.py

The mark: the lactate CH quartet (1 : 3 : 3 : 1) as a white measured spectrum over its fitted area in Manchester
gold, on a Manchester-purple rounded tile. The peaks are real Lorentzian lines computed here. Small icons (32 px
and below) get a thicker line so the four peaks stay visible. Writes:

    mancq/gui/assets/mancq_logo.svg          master drawing
    mancq/gui/assets/mancq_logo_<n>.png      16, 24, 32, 48, 64, 128, 256 px
    mancq/gui/assets/mancq.ico               Windows icon (all sizes; window, taskbar, shortcut, installer, .exe)
    mancq/help/images/logo.png               96 px, for the README and the quick start guide

Needs PySide6 (to draw the PNGs) and Pillow (to pack the .ico).
"""
import os
import sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ASSETS = os.path.join(ROOT, "mancq", "gui", "assets")
PURPLE, GOLD, WHITE = "#660099", "#FFCC33", "#FFFFFF"
SIZES = (16, 24, 32, 48, 64, 128, 256)


def lorentz_path(centres, heights, width, x0, x1, y_base, sx, sy, ox, n=400):
    x = np.linspace(x0, x1, n)
    y = np.zeros(n)
    for c, h in zip(centres, heights):
        y += h * (width / 2) ** 2 / ((x - c) ** 2 + (width / 2) ** 2)
    return "M" + " L".join(f"{ox + (a - x0) * sx:.2f},{y_base - b * sy:.2f}" for a, b in zip(x, y))


def svg(stroke=6.0):
    """The logo on a 128 x 128 canvas; stroke = width of the white spectrum line."""
    c, h, w, y0 = [-1.5, -0.5, 0.5, 1.5], [0.33, 1, 1, 0.33], 0.17, 96
    sx, sy, ox = 128 / 4.8 * 0.78, 64, 14
    meas = lorentz_path(c, h, w, -2.4, 2.4, y0, sx, sy, ox)
    area = lorentz_path(c, [v * 0.92 for v in h], w * 1.9, -2.4, 2.4, y0, sx, sy, ox) + f" L114,{y0} L14,{y0} Z"
    return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 128 128" width="128" height="128" role="img" '
            'aria-label="MANC-Q">'
            f'<rect x="1.5" y="1.5" width="125" height="125" rx="27" fill="{PURPLE}"/>'
            f'<path d="{area}" fill="{GOLD}"/>'
            f'<path d="{meas}" fill="none" stroke="{WHITE}" stroke-width="{stroke}" stroke-linecap="round" '
            'stroke-linejoin="round"/></svg>')


def stroke_for(px):
    return 11.0 if px <= 24 else 9.0 if px <= 32 else 6.0


def render_png(path, px):
    from PySide6.QtCore import QByteArray, QRectF, Qt
    from PySide6.QtGui import QImage, QPainter
    from PySide6.QtSvg import QSvgRenderer
    img = QImage(px, px, QImage.Format_ARGB32)
    img.fill(Qt.transparent)
    p = QPainter(img)
    p.setRenderHint(QPainter.Antialiasing)
    QSvgRenderer(QByteArray(svg(stroke_for(px)).encode())).render(p, QRectF(0, 0, px, px))
    p.end()
    assert img.save(path, "PNG"), path


def main():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtGui import QGuiApplication
    from PIL import Image
    app = QGuiApplication.instance() or QGuiApplication(sys.argv[:1])
    os.makedirs(ASSETS, exist_ok=True)
    with open(os.path.join(ASSETS, "mancq_logo.svg"), "w", encoding="utf-8") as fh:
        fh.write(svg())
    pngs = {}
    for px in SIZES:
        f = os.path.join(ASSETS, f"mancq_logo_{px}.png")
        render_png(f, px)
        pngs[px] = f
    for d in (os.path.join(ROOT, "mancq", "help", "images"),):
        render_png(os.path.join(d, "logo.png"), 96)
    big = Image.open(pngs[256]).convert("RGBA")
    big.save(os.path.join(ASSETS, "mancq.ico"), format="ICO", sizes=[(s, s) for s in SIZES],
             append_images=[Image.open(pngs[s]).convert("RGBA") for s in SIZES if s != 256])
    print("Logo files written to", ASSETS)
    del app


if __name__ == "__main__":
    main()
