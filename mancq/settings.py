"""
mancq.settings - the default settings (SETTINGS), the settings object of a run (Settings), the version
number and the tier names. Change defaults here; a run's own changes are passed to run().
"""
# ----------------------------- SETTINGS --------------------------------------
SETTINGS = dict(
    DATA_DIR=None,                # folder of EXPNO folders, or a single EXPNO folder
    OUT_DIR=None,
    DILUTION_FILE=None,           # optional CSV: sample,dilution
    PROCNO=1,
    TARGET_POINT_HZ=0.25,         # spectra with finer digital resolution are block-averaged to this
    SAMPLES=None,                 # None = every EXPNO folder found, or a list such as ["10", "11"]
    N_WORKERS=None,               # parallel processes (None = number of CPU cores - 1)
    RESUME=True,                  # carrying on a stopped run: re-use the spectra it finished (same run, same settings)
    GALLERY_MAX_SAMPLES=8,        # the per-metabolite gallery is only drawn for small runs
    LIBRARY_DIR=None,             # None = the library folder shipped with the package
    OVERLAYS=True,                # write one overlay PDF per spectrum
    PLURONIC="auto",              # Pluronic F-68 (poloxamer): "auto" = model it if its PEO line is present, "yes" = always
                                  # model it, "ignore" = leave its regions out of the fit, "no" = not in the samples
    PLURONIC_REGIONS=[(3.69, 3.74), (1.10, 1.20)],   # PEO line and PPO methyl, used when PLURONIC = "ignore"
    SAMPLE_DIRS=None,             # optional {sample: (EXPNO folder, procno)}, e.g. FIDs processed into another folder
    FID_MODE="missing",           # raw FIDs: "missing" = process an FID only where there is no processed spectrum,
                                  # "all" = always process from the FID, "never" = only use TopSpin-processed spectra
    FID_LB=0.3,                   # line broadening (Hz) when MANC-Q processes FIDs
    PROCESSING=None,              # optional {sample: dict(source, lb, ph0, ph1, baseline)}: spectra MANC-Q processes
                                  # itself. source "topspin" = TopSpin's spectrum with these changes (ph0/ph1 = change
                                  # from TopSpin's phase, lb = total line broadening); source "fid" = the raw FID.
                                  # Spectra not listed are used as TopSpin processed them (FIDs: see FID_MODE).
    PEO_DETECT_SNR=500,           # "auto": PEO line at 3.71 ppm taller than this x noise
    PEAKLIST_FIELD_MHZ=800.0,     # field at which the fixed peak lists in the library are defined

    # internal standard (the only reference used): TSP or DSS singlet at 0.000 ppm
    REFERENCE="TSP",              # name used in reports only; TSP-d4 and DSS-d6 are handled identically
    TSP_MM_IN_TUBE=None,          # reference concentration in the NMR tube (mM) - REQUIRED
    TSP_PROTONS=9,
    TSP_SI29_FRACTION=0.0467,     # 29Si satellites are outside the fitted main line
    DEFAULT_DILUTION=1.0,         # sample -> tube dilution factor used when none is given

    # fitting
    FIT_RANGE=(0.60, 9.60),
    EXCLUDE=[(4.60, 5.00)],       # regions left out of the fit (ppm); default: residual water
    PEO_CORE_HZ=8.0,              # Pluronic PEO line centre excluded when taller than PEO_CORE_MIN_SNR x noise
    PEO_CORE_MIN_SNR=20000,
    MAX_SHIFT_HZ=8.0,             # compound-level shift freedom around library + global offset
    MAX_SHIFT_TITRATABLE_HZ=25.0, # pH-sensitive compounds (citrate, histidine, ...)
    MULTIPLET_SHIFT_HZ=8.0,       # extra independent freedom of each proton group (grid searched)
    SHIFT_PRIOR=0.02,             # relative SSR penalty per (10 Hz)^2 of multiplet movement
    MULTIPLET_SHIFT_TITRATABLE_HZ=6.0,
    COMPOUND_GRID_HZ=0.5,
    LINEWIDTH_BASE_HZ=(0.8, 1.5),  # base linewidth = TSP linewidth clipped to this range; each multiplet scales it
    ETA_SEARCH=[0.4, 0.6, 0.8, 1.0],  # Lorentzian fraction of the lineshape (searched once)
    GLOBAL_OFFSET_SEARCH_PPM=0.025,   # library-vs-TSP referencing offset search
    REL_MODEL_ERROR=0.03,         # weighting: sigma_i^2 = noise^2 + (0.03*|y_i|)^2 + (0.005*local max)^2
    TAIL_ERROR=0.005,
    TAIL_WINDOW_PPM=0.02,
    BASELINE_KNOT_PPM=0.04,
    BASELINE_PENALTY=3.0,
    # sharp, common compounds used to measure the library-vs-reference offset; an anchor is only used
    # if its pattern is actually present (correlation >= OFFSET_MIN_CORR)
    OFFSET_ANCHORS=["L-Lactate", "L-Alanine", "L-Valine", "Acetate", "Formate", "L-Isoleucine",
                    "D-Glucose", "L-Tyrosine", "L-Phenylalanine", "Choline", "L-Leucine"],
    OFFSET_MIN_CORR=0.3,
    OFFSET_MIN_ANCHORS=3,
    N_ITER=3,
    ROI_MIN_SNR=5,                # multiplets taller than this take part in local refinement
    ROI_PAD_HZ=6.0,
    ROI_MAX_NFEV=60,
    UNKNOWN_MIN_SNR=15,           # residual peaks above this become free "unassigned" singlets
    UNKNOWN_MIN_FRACTION=0.30,    # ... and only where the residual is >=30% of the local signal
    UNKNOWN_BLOCK_SNR=10,         # ... and not within UNKNOWN_EXCLUSION_HZ of a fitted library line
    UNKNOWN_EXCLUSION_HZ=3.0,
    COMPOUND_TIE=0.3,             # soft constraint: equal per-proton amplitude of multiplets in a compound
    LINE_WINDOW_PPM=0.06,         # each line is computed this far either side of its centre (36 Hz at 600 MHz)
    MULTIPLET_GAP_HZ=60.0,        # library lines of one proton group closer than this form one multiplet ...
    MULTIPLET_MAX_SPAN_PPM=0.10,  # ... and a group wider than this is split at its largest gap
    GLUCOSE_ALPHA_FRACTION=0.36,  # share of alpha-D-glucose at equilibrium in water (the rest is beta)
    UNKNOWN_WIDTH_START=1.2,      # unassigned singlets: starting linewidth (x the reference linewidth) ...
    UNKNOWN_WIDTH_RANGE=(0.7, 4.0),  # ... the range it may take ...
    UNKNOWN_MOVE_HZ=1.0,          # ... and how far the singlet may move from the residual peak it was found at
    MAX_UNKNOWNS=300,
    DROP_GHOSTS=True,             # leave out of the fit a compound that only sits on other compounds' peaks: none of its
                                  # multiplets is dominant, and the amount its signal implies would put another of its
                                  # multiplets above CONTRADICT_MIN_SNR x noise where the data show (almost) nothing
    GHOST_MAX_DOMINANCE=0.6,      # ... "dominant" = at least this share of the fitted signal under the multiplet
    GHOST_SHIFT_HZ=5.0,           # ... or its signal is this far from the library position (GHOST_SHIFT_TITRATABLE_HZ
    GHOST_SHIFT_TITRATABLE_HZ=15.0,  # for pH-sensitive compounds) and the other multiplets it implies get almost
                                  # nothing in the fit. Looked for from the second iteration on, once shifts have settled

    # sample pH estimate from the histidine imidazole singlets (H2 near 7.8 ppm, H4 near 7.05 ppm). It is the pH in
    # the NMR tube at acquisition, an ESTIMATE: trends are reliable, absolute pH about +/-0.2 with these constants
    PH_ESTIMATE=True,
    PH_HIS_PKA=6.05,                  # imidazole pKa of free histidine (approximate literature value, water, 25 C)
    PH_HIS_H2_LIMITS=(8.66, 7.70),    # H2 shift (ppm) fully protonated, fully deprotonated
    PH_HIS_H4_LIMITS=(7.39, 7.02),    # H4 shift (ppm) fully protonated, fully deprotonated
    PH_HIS_H2_WINDOW=(7.70, 7.96),    # where the H2 peak is looked for (ppm; covers about pH 6.5 and up)
    PH_HIS_H4_WINDOW=(7.00, 7.12),    # where the H4 peak is looked for (ppm)
    PH_MIN_SNR=10,                    # a proton is used only if its peak is at least this x noise
    PH_MAX_DISAGREE=0.25,             # H2 and H4 further apart than this (pH units) -> flag "protons_disagree"
    PH_PLACE_HISTIDINE=True,          # fit histidine H2/H4 close to the peaks the pH estimate used
    PH_HIS_PLACE_HZ=3.0,              # ... within this many Hz of the measured peak
    PH_HIS_PLACE_PREDICTED_HZ=25.0,
    # chemical shift assignment table (chemical_shift_assignments.xlsx): metabolites with a value in these tiers in
    # at least this share of the samples
    SHIFT_TABLE_TIERS=["Quantified", "Overlapped (semi-quantitative)", "Manually adjusted"],
    SHIFT_TABLE_MIN_FRACTION=0.25,   # ... or of the position predicted from the other proton, if one was not seen

    # per-sample quality flag (qc_flag in the QC sheet: ok / check / poor, with the reasons in qc_reason)
    QC_TSP_FWHM_HZ=(1.5, 3.0),    # TSP linewidth (Hz) above which a sample is "check" / "poor" (shimming)
    QC_TSP_BROAD_RATIO=1.5,       # TSP this many times broader than the metabolites -> "check" (TSP bound to protein?)
    QC_CHI2=(3.0, 10.0),          # fit chi2 per point -> check / poor
    QC_EXPLAINED=(0.90, 0.75),    # share of the signal explained by the library below this -> check / poor
    QC_UNASSIGNED=(10, 40),       # unassigned singlets -> check / poor
    QC_PHASE_DEG=(1.0, 3.0),      # residual phase error at TSP (degrees) -> check / poor
    QC_MIN_RELAXATION_S=10.0,     # D1 + AQ below this (s) -> check: signals with long T1 may be under-estimated

    # grading
    LOD_SNR=3, LOQ_SNR=10,
    Q_MIN_DOMINANCE=0.60, Q_MAX_MISFIT=0.25, MAX_REPORTER_RATIO=2.0,
    MAX_CONC_OVER_BOUND=1.5,      # reported value may not exceed 1.5x the data upper bound
    SINGLE_REPORTER_MIN_SHARE=0.3,  # a single clean reporter must carry >=30% of the molecule's protons
    CONTRADICT_MIN_SNR=10, CONTRADICT_FRACTION=0.3,   # predicted multiplet >10 sigma but <30% observed -> reject
    O_MIN_DOMINANCE=0.25, O_MAX_MISFIT=0.50,
    DECONV_MAX_REL_SE=0.35, DECONV_MIN_MULTIPLETS=2,
    DECONV_MAX_MISFIT=0.5,        # a multiplet counts towards DECONV_MIN_MULTIPLETS only if its own signal is at least
                                  # twice the local misfit (not a sliver under a neighbour's imperfect fit)
    DECONV_MIN_DOMINANCE=0.10,    # ... and this compound is at least 10% of the fitted signal under that multiplet
    MAX_REPORTER_DISAGREE=0.20,   # Quantified only if the molecule's other multiplets (above LOQ, any dominance) give a
                                  # median value within 20% of the reporter value; otherwise Overlapped
    SINGLET_NEEDS_SECOND=True,    # a lone singlet reporter is Quantified only if another multiplet of the molecule agrees
                                  # (molecules with no other multiplet outside the ignored regions, e.g. acetate, are exempt)

    # plots
    OVERLAY_REGIONS=[(0.75, 1.10), (1.10, 1.55), (1.55, 2.00), (2.00, 2.35), (2.35, 2.75),
                     (2.75, 3.10), (3.10, 3.35), (3.35, 3.62), (3.62, 3.85), (3.85, 4.15),
                     (4.15, 4.60), (5.00, 5.50), (5.50, 6.60), (6.60, 7.10), (7.10, 7.50),
                     (7.50, 8.00), (8.00, 9.40)],
)

import copy as _copy
import json
import logging
import os
import warnings
import numpy as np
import pandas as pd
_log = logging.getLogger("mancq")      # see mancq.log (window) or the console


HERE = os.path.dirname(os.path.abspath(__file__))


def _read_version():
    """The version number, kept in one place: pyproject.toml. Read from that file when MANC-Q runs from its own
    folder (so an editable install never reports an old number), otherwise from the installed package's
    metadata (importlib.metadata; also in the Windows app)."""
    pp = os.path.join(os.path.dirname(HERE), "pyproject.toml")
    if os.path.exists(pp):
        import re
        with open(pp, encoding="utf-8") as fh:
            text = fh.read()
        m = re.search(r'^version\s*=\s*"([^"]+)"', text, re.M)
        if m and re.search(r'^name\s*=\s*"mancq"', text, re.M):
            return m.group(1)
    try:
        from importlib.metadata import version
        return version("mancq")
    except Exception:
        _log.debug("version number not found in the package metadata", exc_info=True)
        return "unknown"


__version__ = _read_version()

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
DEFAULT_SETTINGS = SETTINGS       # the defaults above; never changed while MANC-Q runs
TIMING = {}
RESULTS_FORMAT = 2                # version of the results folder layout (1 = MANC-Q 1.9.2 and earlier)


class Settings(dict):
    """The settings of one run: the defaults in SETTINGS with this run's changes. Each run makes its own and
    passes it to every function that needs it (as S); it is saved with the results (run_settings.json) and
    read back by regrade, the Review step and apply-edits. A misspelt setting name is an error, not ignored."""
    RUNTIME = {"RUN_ID", "LIBRARY_HASH", "PROCESSING_USED", "MANCQ_VERSION", "RUN_SAMPLES", "RUN_COMPLETE",
               "RESULTS_FORMAT"}

    def __init__(self, *changes, **kw):
        super().__init__(_copy.deepcopy(SETTINGS))
        self.update(*changes, **kw)

    def __setitem__(self, k, v):
        if k not in SETTINGS and k not in self.RUNTIME:
            raise KeyError(f"unknown setting {k!r}")
        super().__setitem__(k, v)

    def update(self, *changes, **kw):
        for d in changes + (kw,):
            for k, v in dict(d).items():
                self[k] = v

    def setdefault(self, k, v=None):
        if k not in self:
            self[k] = v
        return self[k]

    def __reduce__(self):            # pickled to the worker processes as a plain mapping
        return (Settings, (dict(self),))

    def save(self, path, **extra):
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(dict({k: v for k, v in self.items() if k != "OVERLAY_REGIONS"}, MANCQ_VERSION=__version__,
                           RESULTS_FORMAT=RESULTS_FORMAT, **extra), fh, indent=1, default=str)

    @classmethod
    def load(cls, out):
        """The settings a results folder was made with (_per_sample/run_settings.json), or the defaults if it has
        none. Settings this version does not know are left out; JSON lists are turned back into tuples where the
        default is a tuple."""
        f = os.path.join(out, "_per_sample", "run_settings.json")
        saved = {}
        if os.path.exists(f):
            with open(f, encoding="utf-8") as fh:
                saved = json.load(fh)
        S = cls()
        for k, v in saved.items():
            if k not in SETTINGS and k not in cls.RUNTIME:
                continue
            dflt = SETTINGS.get(k)
            if isinstance(dflt, tuple) and isinstance(v, list):
                v = tuple(v)
            elif isinstance(dflt, list) and dflt and isinstance(dflt[0], tuple) and isinstance(v, list):
                v = [tuple(x) for x in v]
            S[k] = v
        S.pop("SAMPLE_DIRS", None)
        S["SAMPLE_DIRS"] = None
        return S


# =============================================================================
# 1. Reading
# =============================================================================
DECONV = "Deconvolution estimate (low confidence)"
NOTMEAS = "Not measurable (region ignored)"
MANUAL = "Manually adjusted"


def _json_default(v):
    """Lets json.dumps write numpy numbers and arrays (used for the saved fit state)."""
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, np.generic):
        return v.item()
    raise TypeError(f"cannot save {type(v)}")
