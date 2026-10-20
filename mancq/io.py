"""
mancq.io - reading spectra, acquisition parameters, processing and dilution files; results folders (the
overwrite check) and the files of manual adjustments.
"""
import json
import logging
import os
import time
import warnings
import nmrglue as ng
import numpy as np
import pandas as pd
_log = logging.getLogger("mancq")      # see mancq.log (window) or the console


def read_spectrum(expdir, procno, *, S):
    pdir = os.path.join(expdir, "pdata", str(procno))
    name = os.path.basename(os.path.normpath(expdir))
    if os.path.exists(os.path.join(expdir, "acqu2s")) or os.path.exists(os.path.join(pdir, "2rr")):
        raise ValueError(f"EXPNO {name} is a 2D experiment. MANC-Q quantifies 1D 1H spectra: leave it out "
                         f"(in the window, untick it in step 1; on the command line, use --samples).")
    if not os.path.isdir(pdir):
        have = sorted(os.listdir(os.path.join(expdir, "pdata"))) if os.path.isdir(os.path.join(expdir, "pdata")) else []
        raise ValueError(f"EXPNO {name} has no processed spectrum pdata/{procno}"
                         + (f" (it has pdata/{', pdata/'.join(have)}: choose that processing number)" if have else
                            " (process it in TopSpin, or let MANC-Q process its raw FID)") + ".")
    if not os.path.exists(os.path.join(pdir, "procs")):
        raise ValueError(f"EXPNO {name}: pdata/{procno} has no 'procs' file, so the spectrum's scale is unknown. It "
                         f"was probably copied without it: copy the whole EXPNO folder again.")
    if not os.path.exists(os.path.join(pdir, "1r")):
        raise ValueError(f"EXPNO {name}: pdata/{procno} has no 1r file (the processed 1D spectrum).")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dic, data = ng.bruker.read_pdata(pdir, scale_data=True, read_acqus=False)
    p = dic["procs"]
    if np.ndim(data) != 1:
        raise ValueError(f"{pdir}: not a 1D spectrum")
    si = len(data)
    # TopSpin convention: point i is at OFFSET - i * SW_p / SF / SI  (ppm)
    ppm = float(p["OFFSET"]) - np.arange(si) * float(p["SW_p"]) / float(p["SF"]) / si
    y = np.asarray(data, float)
    o = np.argsort(ppm)
    ppm, y = ppm[o], y[o]
    sf = float(dic["procs"]["SF"])
    # spectra processed with a large SI are block-averaged to ~TARGET_POINT_HZ per point:
    # linewidths are >= 1 Hz, so nothing is lost, and memory/time stay the same for any SI
    dx_hz = abs(ppm[1] - ppm[0]) * sf
    k = max(1, int(round(S["TARGET_POINT_HZ"] / dx_hz)))
    if k >= 2:
        n = (len(ppm) // k) * k
        ppm = ppm[:n].reshape(-1, k).mean(1)
        y = y[:n].reshape(-1, k).mean(1)
    o = np.arange(len(ppm))
    title = ""
    tf = os.path.join(pdir, "title")
    if os.path.exists(tf):
        with open(tf, encoding="latin-1") as fh:
            title = fh.read().strip()
    return ppm, y, sf, title


def _usable(d, procno, fids):
    return os.path.exists(os.path.join(d, "pdata", str(procno), "1r")) or \
        (fids and os.path.exists(os.path.join(d, "fid")) and os.path.exists(os.path.join(d, "acqus")))


def list_samples(data_dir, procno, fids=False):
    """EXPNO folders under data_dir with pdata/<procno>/1r (or, if fids, a raw fid + acqus). If data_dir is itself
    an EXPNO folder, returns [""]."""
    if os.path.exists(os.path.join(data_dir, "acqus")) or os.path.isdir(os.path.join(data_dir, "pdata")):
        return [""] if _usable(data_dir, procno, fids) else []
    out = []
    for d in sorted(os.listdir(data_dir), key=lambda x: (len(x), x)):
        if os.path.isdir(os.path.join(data_dir, d)) and _usable(os.path.join(data_dir, d), procno, fids):
            out.append(d)
    return out


def _num(v):
    """A number from a settings file or table; blank / None / NaN -> None."""
    if v is None or (isinstance(v, str) and not v.strip()):
        return None
    v = float(v)
    return None if not np.isfinite(v) else v


def _flag(v):
    if isinstance(v, str):
        return v.strip().lower() in ("1", "true", "yes", "on")
    return bool(v) if v is not None and not (isinstance(v, float) and np.isnan(v)) else False


def read_processing(path):
    """Processing table (CSV with columns sample, source, lb, ph0, ph1, baseline; blank ph0 = automatic phase),
    as MANC-Q writes it to processing_used.csv -> the PROCESSING setting."""
    df = pd.read_csv(path, dtype={"sample": str, "source": str})
    out = {}
    for _, r in df.iterrows():
        src = str(r.get("source", "topspin")).strip().lower()
        out[str(r["sample"]).strip()] = dict(source="fid" if src == "fid" else "topspin", lb=_num(r.get("lb")),
                                             ph0=_num(r.get("ph0")), ph1=_num(r.get("ph1")),
                                             baseline=_flag(r.get("baseline")))
    return out


def read_acquisition(expdir):
    """Acquisition parameters that matter for quantification, from the EXPNO folder's acqus (NaN if absent):
    D1 relaxation delay (s), AQ acquisition time (s), NS scans, RG receiver gain, PULPROG, TE temperature (K),
    P1 90 degree pulse (us)."""
    out = dict(D1_s=np.nan, AQ_s=np.nan, NS=np.nan, RG=np.nan, PULPROG="", TE_K=np.nan, P1_us=np.nan)
    f = os.path.join(expdir, "acqus")
    if not os.path.exists(f):
        return out
    try:
        a = ng.bruker.read_acqus_file(expdir)["acqus"]
    except Exception:
        _log.warning("could not read the acqus file in %s", expdir, exc_info=True)
        return out
    num = lambda v: float(v) if isinstance(v, (int, float)) else np.nan
    if isinstance(a.get("D"), list) and len(a["D"]) > 1:
        out["D1_s"] = num(a["D"][1])
    if isinstance(a.get("P"), list) and len(a["P"]) > 1:
        out["P1_us"] = num(a["P"][1])
    if num(a.get("TD")) > 0 and num(a.get("SW_h")) > 0:
        out["AQ_s"] = num(a["TD"]) / 2 / num(a["SW_h"])
    out.update(NS=num(a.get("NS")), RG=num(a.get("RG")), TE_K=num(a.get("TE")),
               PULPROG=str(a.get("PULPROG", "")).strip("<> "))
    return out


def read_dilutions(path, samples=None, log=None, *, S):
    """{sample: dilution} from a CSV or Excel file with the sample in the first column and the dilution in the
    second (headers, EXPNO prefixes and 10.0 for 10 are tolerated). With samples, every sample without an entry
    and every entry with no sample is reported through log."""
    if not path:
        return {}
    if not os.path.exists(path):
        raise FileNotFoundError(f"dilution file not found: {path}")
    from .figures import _read_table, _sample_id
    df = _read_table(path)
    known = set(map(str, samples or []))
    out = {}
    for k, v in zip(df.iloc[:, 0], df.iloc[:, 1]):
        s = _sample_id(k, known)
        try:
            out[s] = float(str(v).strip())
        except ValueError:
            if s:
                raise ValueError(f"dilution file {path}: the dilution for {s} is not a number ({v!r})")
    if samples is not None and log is not None:
        missing = [s for s in samples if s not in out]
        extra = [s for s in out if s not in known]
        if missing:
            log(f"WARNING: no dilution for EXPNO {', '.join(missing)} in {os.path.basename(path)}: "
                f"the default factor x{S['DEFAULT_DILUTION']:g} is used for them")
        if extra:
            log(f"WARNING: {os.path.basename(path)} has dilutions for {', '.join(extra)}, which are not among the "
                f"spectra of this run (typing error, or a different EXPNO?)")
    return out


# =============================================================================
# Results folders: refuse to mix a new run with the results of an earlier one
# =============================================================================
EXISTING_MESSAGE = "This folder already has MANC-Q results. Running here will overwrite them."
# what a run writes at the top of its results folder (overlay PDFs and simulated spectra are matched by pattern);
# _cache (the simulated library), figures/ and sample_sheet.csv (the user's Figures work) and other files are kept
_RESULT_FILES = ["concentrations_long.csv", "concentrations_report_mM.csv", "peak_list_by_metabolite.csv",
                 "observed_peaks_assigned.csv", "qc_per_sample.csv", "manual_edits.csv", "processing_used.csv",
                 "metabolite_concentrations.xlsx", "chemical_shift_assignments.xlsx", "metabolite_gallery.pdf",
                 "concentration_profiles.pdf", "run_log.txt", "regrade_log.txt", "MANC-Q settings used.json",
                 "dilutions_used.csv"]
_RESULT_DIRS = ["_per_sample", "_processed"]


class ResultsExist(Exception):
    """The results folder already holds MANC-Q results and the run was not told what to do with them."""


def has_results(out):
    """True if `out` already holds MANC-Q results (concentrations_long.csv or _per_sample)."""
    return bool(out) and any(os.path.exists(os.path.join(out, m)) for m in ("concentrations_long.csv", "_per_sample"))


def previous_run(out):
    """What run_settings.json says about the run in `out`: dict(run_id, samples, complete), or None."""
    f = os.path.join(out, "_per_sample", "run_settings.json")
    try:
        with open(f, encoding="utf-8") as fh:
            rs = json.load(fh)
    except (OSError, ValueError):
        return None
    return dict(run_id=rs.get("RUN_ID"), samples=rs.get("RUN_SAMPLES"), complete=bool(rs.get("RUN_COMPLETE")))


def new_run_folder(out):
    """A new, not yet existing subfolder of `out` for a run, e.g. <out>/run_2026-10-19_1340."""
    base = os.path.join(out, f"run_{time.strftime('%Y-%m-%d_%H%M')}")
    f, k = base, 2
    while os.path.exists(f):
        f, k = f"{base}_{k}", k + 1
    return f


def clear_results(out, keep=()):
    """Delete the results of earlier runs from `out`: the tables, workbooks, PDFs, logs, simulated spectra,
    _per_sample (including manual adjustments) and _processed. Files in `keep` are left alone."""
    import shutil
    import glob
    keep = {os.path.normcase(os.path.abspath(k)) for k in keep if k}
    files = [os.path.join(out, f) for f in _RESULT_FILES]
    files += glob.glob(os.path.join(out, "overlay_sample_*.pdf"))
    files += glob.glob(os.path.join(out, "simulated_spectrum_sample_*.csv.gz"))
    for f in files:
        if os.path.isfile(f) and os.path.normcase(os.path.abspath(f)) not in keep:
            os.remove(f)
    for d_ in _RESULT_DIRS:
        p = os.path.join(out, d_)
        if os.path.isdir(p):
            shutil.rmtree(p)


def sample_source(smp, *, S):
    """(EXPNO folder, procno) for a sample: SAMPLE_DIRS if given, else DATA_DIR/<sample>, PROCNO."""
    sd = S.get("SAMPLE_DIRS") or {}
    if smp in sd:
        return sd[smp]
    return os.path.join(S["DATA_DIR"], smp), S["PROCNO"]


def _raw_folder(smp, expdir, *, S):
    """The user's own EXPNO folder of a sample (with its acqus), also when MANC-Q processed it into _processed."""
    raw = os.path.join(S["DATA_DIR"], smp) if S.get("DATA_DIR") else expdir
    return raw if os.path.exists(os.path.join(raw, "acqus")) else expdir


# =============================================================================
# 9. Manual adjustments (made in the Review step of the window)
# =============================================================================
def manual_file(out, smp):
    return os.path.join(out, "_per_sample", f"{smp}_manual.json")


def read_manual(out, smp):
    """{metabolite: edit} for one sample, or {} (edits are kept in _per_sample/<sample>_manual.json)."""
    f = manual_file(out, smp)
    if not os.path.exists(f):
        return {}
    with open(f, encoding="utf-8") as fh:
        return json.load(fh).get("edits", {})
