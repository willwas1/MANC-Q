"""Builds the stand-alone application with PyInstaller: dist/MANC-Q/MANC-Q(.exe).

Used by packaging/build_windows.bat and by the GitHub Actions workflow. Run from the repository folder after
`pip install . pyinstaller`:   python packaging/build_app.py
"""
import os
import sys

import PyInstaller.__main__

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)          # the repository folder, so mancq is found even with an editable install
EXCLUDE = ["PySide6.QtWebEngineCore", "PySide6.QtWebEngineWidgets", "PySide6.QtQuick", "PySide6.QtQml",
           "PySide6.Qt3DCore", "PySide6.QtMultimedia", "PySide6.QtPdf", "PySide6.QtCharts", "PySide6.QtSql",
           "PySide6.QtNetwork", "tkinter"]

args = ["--noconfirm", "--clean", "--windowed", "--name", "MANC-Q", "--paths", ROOT,
        "--add-data", os.path.join(ROOT, "mancq", "library") + os.pathsep + "mancq/library",
        "--hidden-import", "mancq.gui.app", "--hidden-import", "mancq.gui.runner"]
for m in EXCLUDE:
    args += ["--exclude-module", m]
args.append(os.path.join(HERE, "mancq_app.py"))

if __name__ == "__main__":
    PyInstaller.__main__.run(args + sys.argv[1:])
