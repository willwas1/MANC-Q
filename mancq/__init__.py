"""MANC-Q (Metabolite Abundance by NMR Curve-fitting, Quantified): 1H NMR metabolite quantification by whole-spectrum fitting of simulated GISSMO spin systems."""
from .engine import __version__, run, read_spectrum, fit_spectrum, grade, build_library, prepare_library, SETTINGS

__all__ = ["__version__", "run", "read_spectrum", "fit_spectrum", "grade", "build_library", "prepare_library",
           "SETTINGS"]
