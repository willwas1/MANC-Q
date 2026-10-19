"""mancq.figures - publication figures of 1H NMR spectra from a finished MANC-Q results folder.

Used by the Figures step of the window; everything here also works from Python without the window.
Four templates, following what published 1H NMR papers on media and biofluids show:

  representative  one spectrum, water region cut out, optional magnified aromatic region, numbered or abbreviated
                  peak labels (key in the caption), tall peaks cut off and marked ~
  fit             fit evidence: measured spectrum with the summed fit, each chosen metabolite's fitted spectrum on
                  its own row (same vertical scale) and the residual
  stack           time course or replicates, staggered on one vertical scale
  overlay         two or more spectra overlaid in colour, with a zoomed inset (dashed box marks its region)

A figure is described by a plain dictionary (a "spec", saved as JSON next to an exported figure), so the same
figure can be drawn again. Spectra are the ones MANC-Q fitted (simulated_spectrum_sample_<EXPNO>.csv.gz), each
divided by its own reference (TSP/DSS) area, so all traces share one intensity scale. Groups of samples
(condition, time, replicate) come from a sample sheet: a forgiving CSV or Excel file, see read_sample_sheet.
"""
import os
import re
import csv
import io
import json
import time
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib import rc_context
from scipy.signal import find_peaks

from . import engine as q

TEMPLATES = {"representative": "Representative spectrum", "fit": "Fit evidence",
             "stack": "Time course (stacked)", "overlay": "Overlay with inset"}
WIDTHS_MM = {"Single column (85 mm)": 85.0, "1.5 column (120 mm)": 120.0, "Double column (180 mm)": 180.0}
SYMBOLS = {"*": "contaminant or poloxamer", "~": "peak cut off at the top of the plot", "#": "13C satellite"}
MIN_FONT_PT, MIN_LINE_PT = 6.0, 0.25          # preflight: smallest text and thinnest line most journals accept
# colour-blind safe colours (Okabe and Ito), for overlays and fitted components
FIT_COLOUR = "#660099"       # Manchester purple: the default colour of the fitted (simulated) spectrum
COLOURS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000", "#F0E442"]
ABBREV = {
    "L-Lactate": "Lac", "L-Alanine": "Ala", "L-Glutamine": "Gln", "L-Glutamate": "Glu", "D-Glucose": "Glc",
    "L-Valine": "Val", "L-Leucine": "Leu", "L-Isoleucine": "Ile", "L-Threonine": "Thr", "Glycine": "Gly",
    "L-Serine": "Ser", "L-Proline": "Pro", "L-Lysine": "Lys", "L-Arginine": "Arg", "L-Asparagine": "Asn",
    "L-Aspartate": "Asp", "L-Methionine": "Met", "L-Cysteine": "Cys", "L-Cystine": "Cys2", "L-Histidine": "His",
    "L-Phenylalanine": "Phe", "L-Tyrosine": "Tyr", "L-Tryptophan": "Trp", "L-Ornithine": "Orn",
    "L-Citrulline": "Cit", "Pyruvate": "Pyr", "Citrate": "Cit", "Succinate": "Succ", "Fumarate": "Fum",
    "L-Malate": "Mal", "2-Oxoglutarate": "2-OG", "Acetate": "Ace", "Formate": "For", "Ethanol": "EtOH",
    "Choline": "Cho", "Phosphocholine": "PC", "Glycerophosphocholine": "GPC", "Creatine": "Cr",
    "Creatinine": "Crn", "Taurine": "Tau", "myo-Inositol": "m-Ino", "Glycerol": "Glyc", "Betaine": "Bet",
    "Nicotinamide": "NAM", "Pantothenate": "Pan", "HEPES": "HEPES", "Pyroglutamate": "pGlu", "Uridine": "Urd",
    "Hypoxanthine": "Hx", "Xanthine": "Xan", "Inosine": "Ino", "Adenosine": "Ado", "Guanosine": "Guo",
    "Thymidine": "Thd", "Uracil": "Ura", "Putrescine": "Put", "Glutathione_reduced": "GSH",
    "Glutathione_oxidized": "GSSG", "Isocitrate": "Icit", "cis-Aconitate": "Aco", "Isovalerate": "IVA",
    "Isobutyrate": "IBA", "2-Hydroxyisocaproate": "HICA", "2-Hydroxyisovalerate": "HIVA",
    "3-Hydroxybutyrate": "3-HB", "2-Hydroxybutyrate": "2-HB", "Acetoacetate": "AcAc", "Acetone": "Acn",
    "Methanol": "MeOH", "D-Mannose": "Man", "D-Galactose": "Gal", "D-Fructose": "Fru", "Sucrose": "Suc",
    "Folate": "Fol", "Thiamine": "Thi", "Riboflavin": "Rib", "Pyridoxine": "PN", "Biotin": "Bio",
    "Ascorbate": "Asc", "Dimethyl_sulfoxide": "DMSO", "Dimethylsulfone": "DMSO2", "Ethanolamine": "Etn",
    "Phosphoethanolamine": "PE", "Trimethylamine-N-oxide": "TMAO", "Trimethylamine": "TMA", "NAD+": "NAD+",
    "ATP": "ATP", "ADP": "ADP", "AMP": "AMP", "UMP": "UMP", "4-Aminobenzoate": "PABA",
    "L-Kynurenine": "Kyn", "Kynurenate": "KynA", "Anthranilate": "Ant", "Hippurate": "Hip",
    "Propionate": "Prop", "Butyrate": "But", "Glycolate": "Glyco", "beta-Alanine": "b-Ala", "GABA": "GABA",
    "Sarcosine": "Sar", "Dimethylglycine": "DMG", "L-Carnitine": "Car", "O-Acetylcarnitine": "ALCAR",
    "N-Acetyl-L-glutamine": "NAG", "N-Acetylaspartate": "NAA", "Allantoin": "All", "Orotate": "Oro",
    "Trigonelline": "Trig", "1-Methylnicotinamide": "MNA", "Nicotinate": "NA", "Carnosine": "Carn",
    "Anserine": "Ans", "Guanidoacetate": "GAA", "L-Homoserine": "Hse", "2-Aminobutyrate": "Abu",
}


def abbrev(name):
    return ABBREV.get(name, name.replace("L-", "").replace("D-", "")[:4])


def plain_name(name):
    """Metabolite name for a caption: 'L-Lactate' -> 'lactate', 'D-Glucose' -> 'glucose'."""
    n = name.replace("_", " ")
    n = re.sub(r"^(L|D)-", "", n)
    return n if any(ch.isupper() for ch in n[1:]) or "+" in n else n[:1].lower() + n[1:]


# =============================================================================
# Sample sheet: which samples belong together (condition, time point, replicate)
# =============================================================================
SHEET_FILE = "sample_sheet.csv"            # MANC-Q's own copy, in the results folder
SHEET_COLUMNS = ["sample", "nmr_title", "group", "time_h", "replicate", "label", "notes"]
_SYN = {
    "sample": ["sample", "expno", "exp_no", "exp", "experiment", "nmr_expno", "nmr_sample", "sample_no", "sample_number",
               "nmr", "spectrum", "id"],
    "group": ["group", "condition", "treatment", "culture", "reactor", "bioreactor", "vessel", "flask", "arm", "series",
              "medium", "media", "feed", "cell_line", "clone", "run"],
    "time_h": ["time_h", "time_hours", "time_hr", "time_hrs", "hours", "hour", "hrs", "hr", "h", "time", "timepoint",
               "time_point", "t", "sampling_time", "culture_time", "elapsed_time"],
    "time_d": ["time_d", "time_days", "day", "days", "culture_day", "d"],
    "replicate": ["replicate", "rep", "replicate_no", "bio_rep", "biological_replicate", "repeat", "n"],
    "label": ["label", "name", "sample_name", "description", "legend"],
    "notes": ["notes", "note", "comment", "comments", "remarks"],
}


def _norm_header(h):
    h = str(h).strip().lower().replace("\ufeff", "")
    h = re.sub(r"\((h|hr|hrs|hours?)\)", "_h", h)
    h = re.sub(r"\((d|days?)\)", "_d", h)
    h = re.sub(r"[^a-z0-9+]+", "_", h).strip("_")
    return h


def parse_time_h(v, default_unit="h"):
    """Hours from '12', '12 h', '12h', '0.5 hr', '2 d', 'day 3', 'D3', '48:00'. None if not a time."""
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return None
    if isinstance(v, (int, float, np.integer, np.floating)):
        return float(v) * (24.0 if default_unit == "d" else 1.0)
    s = str(v).strip().lower().replace(",", ".")
    if not s:
        return None
    m = re.match(r"^(\d+):(\d{2})$", s)
    if m:
        return int(m.group(1)) + int(m.group(2)) / 60.0
    m = re.match(r"^(day|d)\s*(-?\d+(\.\d+)?)$", s)
    if m:
        return float(m.group(2)) * 24.0
    m = re.match(r"^(-?\d+(\.\d+)?)\s*(h|hr|hrs|hour|hours|d|day|days|min|mins|minutes)?$", s)
    if not m:
        return None
    x, unit = float(m.group(1)), m.group(3) or default_unit
    return x * 24.0 if unit.startswith("d") else x / 60.0 if unit.startswith("m") else x


def _sample_id(v, known):
    """'10', '10.0', 'EXPNO 10', ' 10 ' -> '10' (matching the results' sample names where possible)."""
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return ""
    s = str(v).strip()
    if re.match(r"^\d+\.0+$", s):
        s = s.split(".")[0]
    if known and s not in known:
        m = re.findall(r"\d+", s)
        if len(m) == 1 and m[0] in known:
            s = m[0]
    return s


def _read_table(path):
    """Any CSV/TSV/semicolon file or Excel sheet, in any common encoding -> DataFrame of text."""
    if str(path).lower().endswith((".xlsx", ".xlsm", ".xls")):
        return pd.read_excel(path, dtype=str)
    raw = open(path, "rb").read()
    for enc in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            text = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    try:
        sep = csv.Sniffer().sniff(text[:5000], delimiters=",;\t").delimiter
    except csv.Error:
        sep = ","
    return pd.read_csv(io.StringIO(text), sep=sep, dtype=str, skip_blank_lines=True, engine="python",
                       on_bad_lines="skip")


def read_sample_sheet(path, known_samples=None):
    """Read a sample sheet without failing on small mistakes. Returns (DataFrame with columns SHEET_COLUMNS plus any
    extra columns as given, [warnings]). Columns are recognised by name in many spellings (EXPNO or sample;
    group, condition, reactor...; time, hours, time (h), day...). Times may be written 12, 12 h, 0.5 hr, 2 d or
    day 3. Problems are reported as warnings and the rest of the sheet is used."""
    warns = []
    known = set(map(str, known_samples or []))
    try:
        df = _read_table(path)
    except Exception as exc:
        return pd.DataFrame(columns=SHEET_COLUMNS), [f"The file could not be read ({exc}). Save it as CSV and try again."]
    df = df.dropna(how="all")
    df = df.loc[:, [not str(c).lower().startswith("unnamed") or df[c].notna().any() for c in df.columns]]
    if df.empty:
        return pd.DataFrame(columns=SHEET_COLUMNS), ["The sheet is empty."]
    norm = {c: _norm_header(c) for c in df.columns}
    found = {}
    for key, names in _SYN.items():
        for c, n in norm.items():
            if c not in found.values() and (n in names or (key != "sample" and any(n.startswith(x + "_") for x in names
                                                                                     if len(x) > 2))):
                found[key] = c
                break
    if "sample" not in found:
        # the column whose values best match the samples in the results, else the first column
        best = max(df.columns, key=lambda c: sum(_sample_id(v, known) in known for v in df[c])) if known else df.columns[0]
        found["sample"] = best
        warns.append(f"No column called 'sample' or 'EXPNO': used the column '{best}' for the sample numbers.")
    out = pd.DataFrame()
    out["sample"] = [_sample_id(v, known) for v in df[found["sample"]]]
    out["group"] = df[found["group"]].fillna("").astype(str).str.strip() if "group" in found else ""
    times, bad = [], []
    tcol, unit = (found.get("time_h"), "h") if "time_h" in found else (found.get("time_d"), "d")
    for v in (df[tcol] if tcol else [None] * len(df)):
        t = parse_time_h(v, unit)
        if t is None and v is not None and str(v).strip() not in ("", "nan"):
            bad.append(str(v))
        times.append(np.nan if t is None else t)
    out["time_h"] = times
    if bad:
        warns.append(f"{len(bad)} time value(s) not understood and left blank, e.g. '{bad[0]}'. Write times as "
                     f"hours (12 or 12 h) or days (2 d).")
    out["replicate"] = df[found["replicate"]].fillna("").astype(str).str.strip() if "replicate" in found else ""
    out["label"] = df[found["label"]].fillna("").astype(str).str.strip() if "label" in found else ""
    out["notes"] = df[found["notes"]].fillna("").astype(str) if "notes" in found else ""
    out["nmr_title"] = ""
    used = set(found.values())
    for c in df.columns:                      # anything else (feed, VCD, viability, ...) is kept as written
        if c not in used and _norm_header(c) not in ("nmr_title", "title"):
            out[str(c).strip()] = df[c].fillna("").astype(str).values
    blank = out["sample"] == ""
    if blank.any():
        warns.append(f"{int(blank.sum())} row(s) without a sample number were skipped.")
    out = out[~blank]
    dup = out["sample"].duplicated()
    if dup.any():
        warns.append(f"Sample(s) listed twice, first row used: {', '.join(sorted(set(out['sample'][dup])))}.")
        out = out[~dup]
    if known:
        unknown = sorted(set(out["sample"]) - known, key=lambda s: (len(s), s))
        if unknown:
            warns.append(f"{len(unknown)} sample(s) in the sheet are not in these results and are ignored: "
                         f"{', '.join(unknown[:8])}{' ...' if len(unknown) > 8 else ''}.")
        missing = sorted(known - set(out["sample"]), key=lambda s: (len(s), s))
        if missing:
            warns.append(f"{len(missing)} sample(s) in the results are not in the sheet (shown on their own): "
                         f"{', '.join(missing[:8])}{' ...' if len(missing) > 8 else ''}.")
    cols = SHEET_COLUMNS + [c for c in out.columns if c not in SHEET_COLUMNS]
    return out[cols].reset_index(drop=True), warns


def blank_sheet(samples, titles=None):
    """A sheet with one row per sample and the NMR titles filled in, for the user to complete."""
    titles = titles or {}
    return pd.DataFrame(dict(sample=list(samples), nmr_title=[titles.get(s, "") for s in samples], group="",
                             time_h=np.nan, replicate="", label="", notes=""))[SHEET_COLUMNS]


def load_sheet(out, samples=None, titles=None):
    """MANC-Q's saved sample sheet for a results folder, completed with any samples it does not list."""
    f = os.path.join(out, SHEET_FILE)
    if os.path.exists(f):
        df, _ = read_sample_sheet(f, samples)
    else:
        df = blank_sheet([])
    if samples is not None:
        extra = [s for s in samples if s not in set(df["sample"])]
        if extra:
            df = pd.concat([df, blank_sheet(extra, titles)], ignore_index=True)
        df = df[df["sample"].isin(list(samples))]
        order = {s: i for i, s in enumerate(samples)}
        df = df.sort_values("sample", key=lambda c: c.map(order)).reset_index(drop=True)
    if titles:
        df["nmr_title"] = [titles.get(s, t) for s, t in zip(df["sample"], df["nmr_title"])]
    return df


def save_sheet(out, df):
    df.to_csv(os.path.join(out, SHEET_FILE), index=False)


# =============================================================================
# Data: spectra, fits and peak positions from a results folder
# =============================================================================
class Results:
    """A finished results folder, read once: QC, values, peak list, and each sample's fitted spectrum on demand."""

    def __init__(self, out):
        from . import results
        results.check_readable(out)              # a folder from a newer MANC-Q is refused, not misread
        self.out = out
        self.long = pd.read_csv(os.path.join(out, "concentrations_long.csv"), dtype={"sample": str})
        self.qc = pd.read_csv(os.path.join(out, "qc_per_sample.csv"), dtype={"sample": str}).set_index("sample")
        f = os.path.join(out, "peak_list_by_metabolite.csv")
        self.peaks = pd.read_csv(f, dtype={"sample": str}) if os.path.exists(f) else None
        self.samples = list(dict.fromkeys(self.long["sample"]))
        self.titles = {s: str(self.qc.title.get(s, "")) if "title" in self.qc else "" for s in self.samples}
        rs = os.path.join(out, "_per_sample", "run_settings.json")
        self.run_settings = {}
        if os.path.exists(rs):
            with open(rs, encoding="utf-8") as fh:
                self.run_settings = json.load(fh)
        self._spec = {}

    def spectrum(self, s):
        """DataFrame ppm (ascending), observed, simulated_total, baseline, <compounds>..., each divided by the
        sample's reference area (so all samples share one intensity scale); None if missing."""
        if s not in self._spec:
            f = os.path.join(self.out, f"simulated_spectrum_sample_{s}.csv.gz")
            if not os.path.exists(f):
                self._spec[s] = None
            else:
                d = pd.read_csv(f).sort_values("ppm").reset_index(drop=True)
                a = float(self.qc.TSP_area.get(s, np.nan)) if "TSP_area" in self.qc else np.nan
                scale = 1.0 / a if np.isfinite(a) and a > 0 else 1.0 / max(np.nanmax(np.abs(d.observed.values)), 1e-12)
                for c in d.columns:
                    if c != "ppm":
                        d[c] = d[c] * scale
                self._spec[s] = d
                if len(self._spec) > 12:
                    self._spec.pop(next(iter(self._spec)))
        return self._spec[s]

    def sf(self):
        v = self.qc.SF_MHz.dropna() if "SF_MHz" in self.qc else []
        return float(np.median(v)) if len(v) else None

    def excluded(self):
        return [tuple(map(float, r)) for r in self.run_settings.get("EXCLUDE", [(4.6, 5.0)])]

    def reference(self):
        return self.run_settings.get("REFERENCE", "TSP")

    def pulse_program(self, s):
        try:
            from . import twod
            d = twod.sample_dir(self.out, s)
            if d and os.path.exists(os.path.join(d, "acqus")):
                p = twod._jcamp(os.path.join(d, "acqus")).get("PULPROG", "")
                return str(p).strip("<> ") or None
        except Exception:
            pass
        return None

    def tier(self, metabolite, s):
        r = self.long[(self.long.metabolite == metabolite) & (self.long["sample"] == s)]
        if not len(r):
            return "", np.nan
        t = str(r.tier.iloc[0])
        ub = float(r.upper_bound_mM.iloc[0]) if "upper_bound_mM" in r else np.nan
        return t, q.shown_value(t, float(r.conc_mM.iloc[0]), ub)      # the bound for Upper bound only


# ---------- traces (one line in a figure): a sample, or the mean of replicate samples
def trace_groups(R, sheet):
    """Ordered {key: dict(label, samples, replicates)}: one entry per (group, time) in the sheet; samples with no
    group or time are their own entry. Sorted by group, then time."""
    rows = sheet.set_index("sample") if len(sheet) else pd.DataFrame()
    groups = {}
    for s in R.samples:
        g = str(rows.group.get(s, "")) if len(rows) and "group" in rows else ""
        t = rows.time_h.get(s, np.nan) if len(rows) and "time_h" in rows else np.nan
        t = float(t) if t is not None and str(t) not in ("", "nan") else np.nan
        rep = str(rows.replicate.get(s, "")) if len(rows) and "replicate" in rows else ""
        lab = str(rows.label.get(s, "")) if len(rows) and "label" in rows else ""
        if not g and not np.isfinite(t):
            key = f"EXPNO {s}"
            groups[key] = dict(label=lab or f"EXPNO {s}", samples=[s], replicates=[rep or s], group="", time=np.nan,
                               sample_labels=[lab])
            continue
        key = f"{g}|{t:g}" if np.isfinite(t) else f"{g}|"
        e = groups.setdefault(key, dict(label=", ".join(x for x in (g, f"{t:g} h" if np.isfinite(t) else "") if x),
                                        samples=[], replicates=[], group=g, time=t, sample_labels=[]))
        e["samples"].append(s); e["replicates"].append(rep or s); e["sample_labels"].append(lab)
    for e in groups.values():               # a label column shared by every replicate names the trace
        if len(set(e["sample_labels"])) == 1 and e["sample_labels"][0] and len(e["samples"]) > 1:
            e["label"] = e["sample_labels"][0]
    order = sorted(groups, key=lambda k: (groups[k]["group"] == "", groups[k]["group"],
                                          groups[k]["time"] if np.isfinite(groups[k]["time"]) else 1e9,
                                          [int(x) if x.isdigit() else 0 for x in groups[k]["samples"][:1]]))
    return {k: groups[k] for k in order}


def build_trace(R, entry, mode="mean", replicate=None, components=()):
    """One line: mode 'mean' averages the replicate spectra (and their fits) on the first sample's ppm axis;
    'one' uses the replicate named `replicate` (else the first). Returns dict or None."""
    samples = list(entry["samples"])
    if mode != "mean" and len(samples) > 1:
        pick = [s for s, r in zip(samples, entry["replicates"]) if replicate and r == str(replicate)]
        samples = pick[:1] or samples[:1]
    ds = [(s, R.spectrum(s)) for s in samples]
    ds = [(s, d) for s, d in ds if d is not None]
    if not ds:
        return None
    x = ds[0][1].ppm.values
    cols = ["observed", "simulated_total", "baseline"] + list(components)

    def col(d, c):
        v = d[c].values if c in d else np.zeros(len(d))
        return v if len(d) == len(x) and np.allclose(d.ppm.values, x) else np.interp(x, d.ppm.values, v)
    arr = {c: np.mean([col(d, c) for _, d in ds], axis=0) for c in cols}
    used = [s for s, _ in ds]
    lab = entry["label"]
    if len(entry["samples"]) > 1 and mode != "mean":
        rep = entry["replicates"][entry["samples"].index(used[0])]
        lab = f"{lab}, replicate {rep}" if rep and rep != used[0] else f"{lab} (EXPNO {used[0]})"
    return dict(label=lab, samples=used, mean=len(used) > 1, ppm=x, **arr)


# =============================================================================
# Labels
# =============================================================================
def auto_labels(R, samples, tiers=("Quantified",), which="main"):
    """Labels from MANC-Q's own fit: for every metabolite with a value in one of `tiers` in these samples, its
    tallest reporter multiplet ('main'), all its reporter multiplets ('reporter') or all its multiplets ('all'),
    at the median fitted position.
    Returns a list of dicts (show, metabolite, ppm, lo, hi, text, symbol, bracket)."""
    if R.peaks is None:
        return []
    G = R.long[R.long["sample"].isin(samples) & R.long.tier.isin(tiers)]
    mets = list(dict.fromkeys(G.sort_values("conc_mM", ascending=False).metabolite))
    mets = [m for m in mets if not m.startswith("Pluronic")]
    P = R.peaks[R.peaks["sample"].isin(samples) & R.peaks.metabolite.isin(mets)]
    out = []
    for m in mets:
        sub = P[P.metabolite == m]
        if which in ("main", "reporter") and sub.used_as_reporter.astype(bool).any():
            sub = sub[sub.used_as_reporter.astype(bool)]
        elif which in ("main", "reporter"):
            continue
        if which == "main":
            best = sub.groupby("multiplet").height_snr.median().idxmax()
            sub = sub[sub.multiplet == best]
        for j, g in sub.groupby("multiplet"):
            lo, hi, c = float(g.from_ppm.median()), float(g.to_ppm.median()), float(g.centre_ppm.median())
            out.append(dict(show=True, metabolite=m, ppm=round(c, 4), lo=round(lo, 4), hi=round(hi, 4),
                            text=abbrev(m), symbol="", bracket=bool(hi - lo > 0.012)))
    return sorted(out, key=lambda r: -r["ppm"])


def label_numbers(labels):
    """{metabolite: number}, numbered from high to low ppm of each metabolite's first label."""
    order = []
    for r in sorted([r for r in labels if r.get("show") and r.get("metabolite")], key=lambda r: -r["ppm"]):
        if r["metabolite"] not in order:
            order.append(r["metabolite"])
    return {m: i + 1 for i, m in enumerate(order)}


# =============================================================================
# The figure description
# =============================================================================
def default_spec(R=None):
    excl = R.excluded() if R is not None else [(4.6, 5.0)]
    return dict(template="representative", traces=[], mode="mean", replicate="", region="0.5-9.5",
                leave_out=", ".join(f"{a:g}-{b:g}" for a, b in excl), magnify=True, magnify_from=5.5, magnify_factor=8,
                clip=False, clip_percent=30, label_style="numbers", label_tiers=["Quantified"], label_which="main",
                labels=None, fit_metabolites=[], fit_region="", stack_offset=25, inset_region="0.8-1.1",
                baseline_correct=True, overlay_offset=0, offset_x_ppm=0.0, colours={},
                show_fit=False, fit_colour=FIT_COLOUR, width_mm=180.0, height_mm=90.0, font_pt=7.0, line_pt=0.6, title="")


def parse_ranges(text):
    """'0.5-9.5, 4.6 to 5.0' -> [(0.5, 9.5), (4.6, 5.0)] (each low, high)."""
    out = []
    for part in re.split(r"[,;]", str(text or "")):
        m = re.search(r"(-?\d+(?:\.\d+)?)\s*(?:-|to|\.\.|\u2013)\s*(-?\d+(?:\.\d+)?)", part)
        if m:
            a, b = float(m.group(1)), float(m.group(2))
            if a != b:
                out.append((min(a, b), max(a, b)))
    return out


def segments(region, leave_out, magnify=None):
    """Pieces of the ppm axis to draw, high ppm first: [(lo, hi, factor)]. Regions to leave out become breaks;
    the magnified range gets its factor. Pieces narrower than 0.1 ppm are dropped."""
    pieces = []
    for lo, hi in sorted(region, key=lambda r: -r[1]):
        cuts = [(lo, hi)]
        for a, b in leave_out:
            nxt = []
            for c, d in cuts:
                if b <= c or a >= d:
                    nxt.append((c, d))
                else:
                    if a > c:
                        nxt.append((c, a))
                    if b < d:
                        nxt.append((b, d))
            cuts = nxt
        for c, d in sorted(cuts, key=lambda r: -r[1]):
            if magnify and magnify[1] > 1 and c < magnify[0] < d:
                pieces += [(magnify[0], d, magnify[1]), (c, magnify[0], 1.0)]
            else:
                pieces.append((c, d, magnify[1] if magnify and magnify[1] > 1 and c >= magnify[0] else 1.0))
    return [p for p in pieces if p[1] - p[0] >= 0.1]


# =============================================================================
# Drawing
# =============================================================================
_FAMILIES = []


def _families():
    """Arial (or a close equivalent) if installed, else matplotlib's own font; only installed fonts are listed, so
    matplotlib does not warn about the others."""
    if not _FAMILIES:
        from matplotlib import font_manager
        have = {f.name for f in font_manager.fontManager.ttflist}
        _FAMILIES.extend([f for f in ("Arial", "Helvetica", "Liberation Sans") if f in have][:1] + ["DejaVu Sans"])
    return list(_FAMILIES)


def _style(spec):
    return {"font.family": _families(), "font.size": spec["font_pt"],
            "axes.linewidth": spec["line_pt"], "xtick.major.width": spec["line_pt"],
            "xtick.major.size": 2.5, "xtick.minor.size": 1.5, "xtick.labelsize": spec["font_pt"],
            "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "mathtext.default": "regular"}


def _light(colour, f=0.65):
    """A pale solid version of a colour (fills without transparency, so EPS shows them too)."""
    import matplotlib.colors as mc
    r, g, b = mc.to_rgb(colour)
    return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)


COLOUR_NAMES = {"black": "#000000", "grey": "#808080", "blue": "#0072B2", "light blue": "#56B4E9",
                "green": "#009E73", "red": "#D62728", "orange": "#E69F00", "dark orange": "#D55E00",
                "purple": "#7B3294", "pink": "#CC79A7", "yellow": "#F0E442", "brown": "#8C564B"}


def trace_colour(spec, tr, k, default):
    """The colour chosen for a spectrum in the window (spec["colours"], by its key), else the template's own."""
    return (spec.get("colours") or {}).get(tr.get("key", ""), default) or default


def colour_name(colour):
    """The nearest plain colour name, for the caption."""
    import matplotlib.colors as mc
    r, g, b = mc.to_rgb(colour)
    return min(COLOUR_NAMES, key=lambda n: sum((x - y) ** 2 for x, y in zip((r, g, b), mc.to_rgb(COLOUR_NAMES[n]))))


def _shift_labels(labels, dx):
    """Labels moved sideways with an offset spectrum (dx ppm)."""
    if not dx:
        return labels
    return [dict(r, ppm=r["ppm"] - dx, lo=r.get("lo", r["ppm"]) - dx, hi=r.get("hi", r["ppm"]) - dx) for r in labels]


def _axes_row(fig, segs, rect):
    """One axes per segment, side by side, width proportional to ppm span, with break marks between them."""
    l, b, w, h = rect
    gap = 0.012
    total = sum(s[1] - s[0] for s in segs)
    avail = w - gap * (len(segs) - 1)
    axes, x0 = [], l
    for lo, hi, fac in segs:
        ww = avail * (hi - lo) / total
        ax = fig.add_axes([x0, b, ww, h])
        ax.set_xlim(hi, lo)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.set_yticks([])
        ax.patch.set_alpha(0)
        axes.append(ax)
        x0 += ww + gap
    for a1, a2 in zip(axes[:-1], axes[1:]):
        a1._right_break, a2._left_break = True, True
    for a1, a2 in zip(axes[:-1], axes[1:]):          # // marks on the axis line at each break
        for ax, xx in ((a1, 1.0), (a2, 0.0)):
            ax.plot([xx - 0.004 / ax.get_position().width * 0.5, xx + 0.004 / ax.get_position().width * 0.5],
                    [-0.025, 0.025], transform=ax.transAxes, color="k", lw=ax.spines["bottom"].get_linewidth(),
                    clip_on=False)
    return axes


def _ticks(ax, lo, hi, font_pt=7.0):
    """Tick marks at a round step, as many as fit the printed width of the axes without crowding."""
    fig = ax.figure
    width_pt = ax.get_position().width * fig.get_figwidth() * 72.0
    room = max(1, int(width_pt / (3.4 * font_pt)))
    span = hi - lo
    step = next((s for s in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0) if span / s <= room), 5.0)
    t = np.arange(np.ceil(lo / step - 1e-9) * step, hi + 1e-9, step)
    # no tick right at a break in the axis (the neighbouring piece shows that value)
    t = [v for v in t if hi - v > 0.03 * span or not getattr(ax, "_left_break", False)]
    t = [v for v in t if v - lo > 0.03 * span or not getattr(ax, "_right_break", False)]
    if not t:                                   # a narrow piece keeps one tick in its middle
        t = [round((lo + hi) / 2, 1 if span >= 0.3 else 2)]
    ax.set_xticks(t)
    ax.set_xticklabels([f"{v:.2f}".rstrip("0").rstrip(".") if step < 0.1 else f"{v:.1f}" for v in t])


def _visible(x, y, lo, hi, excl_ref=True):
    w = (x >= lo) & (x <= hi)
    if excl_ref:
        w &= np.abs(x) > 0.3
    return w


def _cut(x, *ys, leave_out=()):
    """NaN in the regions left out, and gaps where the data have holes."""
    m = np.zeros(len(x), bool)
    for a, b in leave_out:
        m |= (x > a) & (x < b)
    ys = [np.where(m, np.nan, y) for y in ys]
    return q.gapped(x, *ys)


def draw(fig, spec, R, traces):
    """Draw the figure described by spec on fig (any matplotlib Figure). traces: list of build_trace results.
    Returns info: dict(caption, data (DataFrame of what was plotted), warnings)."""
    fig.clear()
    fig.patch.set_facecolor("white")
    info = dict(warnings=[], plotted={}, clipped=False, magnified=None, inset=None, numbers={})
    if not traces:
        ax = fig.add_subplot(111); ax.axis("off")
        ax.text(0.5, 0.5, "Choose at least one sample or group to show.", ha="center", va="center", transform=ax.transAxes)
        info["caption"] = ""
        return info
    t = spec["template"]
    if spec.get("baseline_correct", True) and t != "fit":
        # the smooth baseline MANC-Q fitted is subtracted, so magnified regions and stacked traces sit on zero
        traces = [dict(tr, observed=tr["observed"] - tr["baseline"], simulated_total=tr["simulated_total"] - tr["baseline"],
                       baseline=np.zeros(len(tr["ppm"]))) for tr in traces]
    with rc_context(_style(spec)):
        if t == "fit":
            _draw_fit(fig, spec, R, traces[0], info)
        elif t == "stack":
            _draw_stack(fig, spec, R, traces, info)
        elif t == "overlay":
            _draw_overlay(fig, spec, R, traces, info)
        else:
            _draw_representative(fig, spec, R, traces[:1], info)
    info["caption"] = caption(spec, R, traces, info)
    return info


def _label_text(r, spec, numbers):
    sym = r.get("symbol") or ""
    if not r.get("metabolite"):                      # a free label or symbol typed by the user
        return (r.get("text") or "") + sym
    if spec["label_style"] == "numbers":
        return f"{numbers.get(r['metabolite'], '')}{sym}"
    return (r.get("text") or abbrev(r["metabolite"])) + sym


def _place_labels(fig, axes_segs, labels, spec, numbers, top_y, ytop_fn, info):
    """Labels above the peaks; in crowded places they are raised in steps and joined to the peak by a thin leader
    line, and multiplets get a bracket. ytop_fn(ppm_lo, ppm_hi, seg_index) -> data height of the spectrum there."""
    if spec["label_style"] == "off":
        return
    fs = spec["font_pt"]
    lw = max(spec["line_pt"] * 0.6, 0.3)
    dpi = fig.dpi
    placed = []                                     # display boxes already used
    rows = [r for r in labels if r.get("show")]
    for r in sorted(rows, key=lambda r: -r["ppm"]):
        seg = next((i for i, (ax, (lo, hi, fac)) in enumerate(axes_segs) if lo <= r["ppm"] <= hi), None)
        if seg is None:
            continue
        ax, (lo, hi, fac) = axes_segs[seg]
        txt = _label_text(r, spec, numbers)
        if not txt:
            continue
        tr = ax.transData
        wpx = (0.6 * len(txt) + 0.4) * fs * dpi / 72.0
        # the label must clear the spectrum under its whole width, not only its own peak
        x_l = tr.inverted().transform((tr.transform((r["ppm"], 0))[0] - wpx / 2, 0))[0]
        hw = abs(x_l - r["ppm"])
        a, b = (r.get("lo", r["ppm"]), r.get("hi", r["ppm"])) if r.get("bracket") else (r["ppm"], r["ppm"])
        a, b = min(a, b, r["ppm"] - hw), max(a, b, r["ppm"] + hw)
        ypk = min(ytop_fn(a, b, seg), top_y)
        xd, yd = tr.transform((r["ppm"], ypk))
        a, b = (r.get("lo", r["ppm"]), r.get("hi", r["ppm"])) if r.get("bracket") else (r["ppm"] - 0.004, r["ppm"] + 0.004)
        r = dict(r, bracket=bool(r.get("bracket")) and
                 abs(tr.transform((a, 0))[0] - tr.transform((b, 0))[0]) >= 1.2 * fs * dpi / 72.0)
        hpx = 1.15 * fs * dpi / 72.0
        base = yd + 0.9 * fs * dpi / 72.0 + (0.6 * fs * dpi / 72.0 if r.get("bracket") else 0)
        y = base
        for _ in range(12):
            box = (xd - wpx / 2, y, xd + wpx / 2, y + hpx)
            if not any(box[0] < p[2] and box[2] > p[0] and box[1] < p[3] and box[3] > p[1] for p in placed):
                break
            y += hpx
        placed.append((xd - wpx / 2, y, xd + wpx / 2, y + hpx))
        inv = tr.inverted()
        _, ylab = inv.transform((xd, y))
        _, ytip = inv.transform((xd, yd + 0.25 * fs * dpi / 72.0))
        if r.get("bracket"):
            _, yb = inv.transform((xd, yd + 0.55 * fs * dpi / 72.0))
            _, yb2 = inv.transform((xd, yd + 0.25 * fs * dpi / 72.0))
            ax.plot([max(a, b), max(a, b), min(a, b), min(a, b)], [yb2, yb, yb, yb2], color="k", lw=lw, clip_on=False)
            info["brackets_drawn"] = True
            ytip = yb
        if y > base + 1:
            _, ybot = inv.transform((xd, y))
            ax.plot([r["ppm"], r["ppm"]], [ytip, ybot], color="#555555", lw=lw * 0.8, clip_on=False)
        ax.text(r["ppm"], ylab, txt, ha="center", va="bottom", fontsize=fs, clip_on=False)
    info["labels_used"] = rows


def _mark_clipped(axes_segs, x, y, top, spec, info):
    """~ above every peak that is cut off by the top of the plot."""
    for ax, (lo, hi, fac) in axes_segs:
        w = (x >= lo) & (x <= hi)
        if not w.any():
            continue
        yy = np.where(np.isfinite(y[w]), y[w] * fac, -np.inf)
        pk, _ = find_peaks(yy, height=top, distance=max(1, int(0.03 / max(np.median(np.diff(x)), 1e-9))))
        for p in pk:
            ax.text(x[w][p], top, "~", ha="center", va="bottom", fontsize=spec["font_pt"] + 1, clip_on=False)
            info["clipped"] = True


def _top_of(traces, segs, clip, clip_percent, offset=0.0):
    tops = []
    for k, tr in enumerate(traces):
        x, y = tr["ppm"], tr["observed"]
        for lo, hi, fac in segs:
            w = _visible(x, y, lo, hi)
            if w.any():
                tops.append(np.nanmax(y[w]) * fac + k * offset)
    top = max(tops) if tops else 1.0
    if clip:
        peak = max(tops) if tops else 1.0
        top = max(min(top, peak * clip_percent / 100.0), 1e-12)
    return top


def _region(spec):
    reg = parse_ranges(spec["region"]) or [(0.5, 9.5)]
    return reg


def _draw_representative(fig, spec, R, traces, info):
    tr = traces[0]
    x = tr["ppm"]
    leave = parse_ranges(spec["leave_out"])
    mag = (spec["magnify_from"], spec["magnify_factor"]) if spec["magnify"] else None
    segs = segments(_region(spec), leave, mag)
    if not segs:
        info["warnings"].append("Nothing to draw: check the ppm region.")
        return
    top = _top_of([tr], segs, spec["clip"], spec["clip_percent"])
    lab_room = 0.30 if spec["label_style"] != "off" else 0.06
    axes = _axes_row(fig, segs, [0.03, 0.17, 0.94, 0.80])
    pairs = list(zip(axes, segs))
    gx, gy = _cut(x, tr["observed"], leave_out=leave)
    for ax, (lo, hi, fac) in pairs:
        ax.plot(gx, gy * fac, color=trace_colour(spec, tr, 0, "k"), lw=spec["line_pt"])
        if spec.get("show_fit"):
            gfx, gf = _cut(x, tr["simulated_total"], leave_out=leave)
            ax.plot(gfx, gf * fac, color=spec.get("fit_colour") or FIT_COLOUR, lw=spec["line_pt"] * 0.8)
        ax.set_ylim(-0.02 * top, top * (1 + lab_room))
        _ticks(ax, lo, hi, spec["font_pt"])
        if fac > 1:
            ax.text(0.5, 0.97, f"\u00d7{fac:g}", transform=ax.transAxes, ha="center", va="top", fontsize=spec["font_pt"])
            info["magnified"] = (lo, hi, fac)
    fig.text(0.5, 0.03, "\u03b4 \u00b9H (ppm)", ha="center", va="bottom", fontsize=spec["font_pt"])
    if spec["clip"]:
        for ax, s in pairs:
            ax.set_ylim(-0.02 * top, top * (1 + lab_room))
        _mark_clipped(pairs, x, tr["observed"], top, spec, info)
    labels = spec.get("labels") or []
    numbers = label_numbers(labels)
    info["numbers"] = numbers

    def ytop(a, b, i):
        lo, hi, fac = segs[i]
        w = (x >= a - 0.003) & (x <= b + 0.003)
        return float(np.nanmax(tr["observed"][w]) * fac) if w.any() else 0.0
    _place_labels(fig, pairs, labels, spec, numbers, top, ytop, info)
    info["plotted"] = {"ppm": x, f"{tr['label']} (observed)": tr["observed"]}
    if spec.get("show_fit"):
        info["plotted"][f"{tr['label']} (fit)"] = tr["simulated_total"]


def _draw_stack(fig, spec, R, traces, info):
    leave = parse_ranges(spec["leave_out"])
    mag = (spec["magnify_from"], spec["magnify_factor"]) if spec["magnify"] else None
    segs = segments(_region(spec), leave, mag)
    if not segs:
        info["warnings"].append("Nothing to draw: check the ppm region.")
        return
    top1 = _top_of(traces, segs, spec["clip"], spec["clip_percent"])
    off = top1 * spec["stack_offset"] / 100.0
    dx = float(spec.get("offset_x_ppm") or 0.0)          # sideways offset per spectrum (to the right)
    n = len(traces)
    top = top1 + off * (n - 1)
    lab_room = 0.18 if spec["label_style"] != "off" else 0.04
    axes = _axes_row(fig, segs, [0.03, 0.15, 0.82, 0.82])
    pairs = list(zip(axes, segs))
    for k, tr in enumerate(traces):
        gx, gy = _cut(tr["ppm"], tr["observed"], leave_out=leave)
        col = trace_colour(spec, tr, k, "k")
        for ax, (lo, hi, fac) in pairs:
            yy = np.minimum(gy * fac, top1) if spec["clip"] else gy * fac
            ax.plot(gx - k * dx, yy + k * off, color=col, lw=spec["line_pt"])
        last_ax, (lo, hi, fac) = pairs[-1]
        w = (tr["ppm"] >= lo) & (tr["ppm"] <= hi)
        ybase = float(np.nanmedian(tr["observed"][w])) * fac + k * off if w.any() else k * off
        last_ax.text(1.02, ybase, tr["label"], transform=last_ax.get_yaxis_transform(), ha="left", va="bottom",
                     fontsize=spec["font_pt"], clip_on=False, color=col)
        info["plotted"].setdefault("ppm", tr["ppm"])
        info["plotted"][f"{tr['label']} (observed)"] = np.interp(info["plotted"]["ppm"], tr["ppm"], tr["observed"])
    for ax, (lo, hi, fac) in pairs:
        ax.set_ylim(-0.02 * top1, top + top1 * lab_room)
        _ticks(ax, lo, hi, spec["font_pt"])
        if fac > 1:
            ax.text(0.5, 1.0, f"\u00d7{fac:g}", transform=ax.transAxes, ha="center", va="top", fontsize=spec["font_pt"])
            info["magnified"] = (lo, hi, fac)
    if spec["clip"]:
        info["clipped"] = True
        for ax, (lo, hi, fac) in pairs:
            for k, tr in enumerate(traces):
                w = (tr["ppm"] >= lo) & (tr["ppm"] <= hi)
                if w.any():
                    yy = tr["observed"][w] * fac
                    pk, _ = find_peaks(np.nan_to_num(yy, nan=-1), height=top1,
                                       distance=max(1, int(0.03 / max(np.median(np.diff(tr["ppm"])), 1e-9))))
                    for p in pk:
                        ax.text(tr["ppm"][w][p] - k * dx, top1 + k * off, "~", ha="center", va="bottom",
                                fontsize=spec["font_pt"], clip_on=False)
    fig.text(0.44, 0.03, "\u03b4 \u00b9H (ppm)", ha="center", va="bottom", fontsize=spec["font_pt"])
    labels = spec.get("labels") or []
    numbers = label_numbers(labels)
    info["numbers"] = numbers
    tr = traces[-1]
    sh = (n - 1) * dx                                   # labels go on the top spectrum, wherever it was moved
    x = tr["ppm"] - sh

    def ytop(a, b, i):
        lo, hi, fac = segs[i]
        w = (x >= a - 0.003) & (x <= b + 0.003)
        v = float(np.nanmax(tr["observed"][w]) * fac) if w.any() else 0.0
        return (min(v, top1) if spec["clip"] else v) + (n - 1) * off
    _place_labels(fig, pairs, _shift_labels(labels, sh), spec, numbers, top, ytop, info)
    info["offset"] = (spec["stack_offset"], dx)


def _draw_overlay(fig, spec, R, traces, info):
    leave = parse_ranges(spec["leave_out"])
    segs = segments(_region(spec), leave, None)
    if not segs:
        info["warnings"].append("Nothing to draw: check the ppm region.")
        return
    top1 = _top_of(traces, segs, spec["clip"], spec["clip_percent"])
    voff = top1 * float(spec.get("overlay_offset") or 0.0) / 100.0
    dx = float(spec.get("offset_x_ppm") or 0.0)
    top = top1 + voff * (len(traces) - 1)
    axes = _axes_row(fig, segs, [0.03, 0.15, 0.94, 0.80])
    pairs = list(zip(axes, segs))
    for k, tr in enumerate(traces):
        gx, gy = _cut(tr["ppm"], tr["observed"], leave_out=leave)
        for ax, _ in pairs:
            ax.plot(gx - k * dx, gy + k * voff, color=trace_colour(spec, tr, k, COLOURS[k % len(COLOURS)]),
                    lw=spec["line_pt"], label=tr["label"])
        info["plotted"].setdefault("ppm", tr["ppm"])
        info["plotted"][f"{tr['label']} (observed)"] = np.interp(info["plotted"]["ppm"], tr["ppm"], tr["observed"])
    for ax, (lo, hi, fac) in pairs:
        ax.set_ylim(-0.02 * top1, top + 0.05 * top1)
        _ticks(ax, lo, hi, spec["font_pt"])
    if spec["clip"]:
        info["clipped"] = True
    info["offset"] = (float(spec.get("overlay_offset") or 0.0), dx)
    axes[0].legend(loc="upper left", frameon=False, fontsize=spec["font_pt"], handlelength=1.5)
    fig.text(0.5, 0.03, "\u03b4 \u00b9H (ppm)", ha="center", va="bottom", fontsize=spec["font_pt"])
    ins = parse_ranges(spec["inset_region"])
    if ins:
        a, b = ins[0]
        host = next((ax for ax, (lo, hi, f) in pairs if lo <= a and b <= hi), None)
        tops = [np.nanmax(tr["observed"][(tr["ppm"] >= a) & (tr["ppm"] <= b)]) for tr in traces
                if ((tr["ppm"] >= a) & (tr["ppm"] <= b)).any()]
        if tops and host is not None:
            itop = max(max(tops), 1e-12)
            fac = top1 / itop
            ivoff = itop * float(spec.get("overlay_offset") or 0.0) / 100.0
            # the place (in any piece of the axis) where the inset hides the least signal
            best = None
            for ax_, (lo_, hi_, _f) in pairs:
                for x0 in (0.03, 0.29, 0.55):
                    p1, p0 = hi_ - x0 * (hi_ - lo_), hi_ - (x0 + 0.42) * (hi_ - lo_)
                    hidden = max([float(np.nanmax(np.where((tr["ppm"] >= p0) & (tr["ppm"] <= p1), tr["observed"], 0)))
                                  for tr in traces] + [0.0])
                    score = hidden / top - 0.02 * (hi_ - lo_)
                    if best is None or score < best[0]:
                        best = (score, ax_, x0)
            ia = best[1].inset_axes([best[2], 0.42, 0.42, 0.52])
            for k, tr in enumerate(traces):
                w = (tr["ppm"] >= a) & (tr["ppm"] <= b)
                ia.plot(tr["ppm"][w] - k * dx, tr["observed"][w] + k * ivoff,
                        color=trace_colour(spec, tr, k, COLOURS[k % len(COLOURS)]), lw=spec["line_pt"])
            ia.set_xlim(b, a); ia.set_ylim(-0.03 * itop, itop * 1.08 + ivoff * (len(traces) - 1))
            ia.set_yticks([])
            for sp in ("top", "right", "left"):
                ia.spines[sp].set_visible(False)
            _ticks(ia, a, b, spec["font_pt"])
            ia.tick_params(labelsize=spec["font_pt"])
            ia.text(0.02, 0.98, f"\u00d7{fac:.0f}" if fac >= 1.5 else "", transform=ia.transAxes, ha="left", va="top",
                    fontsize=spec["font_pt"])
            from matplotlib.patches import Rectangle
            host.add_patch(Rectangle((a, -0.01 * top1), b - a, itop * 1.1 + voff * (len(traces) - 1), fill=False, ls="--", lw=spec["line_pt"] * 0.8,
                                     edgecolor="#555555"))
            info["inset"] = (a, b, fac)


def _draw_fit(fig, spec, R, tr, info):
    reg = parse_ranges(spec["fit_region"]) or parse_ranges(spec["region"]) or [(0.5, 9.5)]
    leave = parse_ranges(spec["leave_out"])
    segs = segments(reg, leave, None)
    if not segs:
        info["warnings"].append("Nothing to draw: check the ppm region.")
        return
    x = tr["ppm"]
    mets = [m for m in spec.get("fit_metabolites") or [] if m in tr]
    base = tr["baseline"]
    obs = tr["observed"] - base
    fit = tr["simulated_total"] - base
    resid = tr["observed"] - tr["simulated_total"]
    w = np.zeros(len(x), bool)
    for lo, hi, _ in segs:
        w |= (x >= lo) & (x <= hi)
    top = float(np.nanmax(np.r_[obs[w], fit[w]])) if w.any() else 1.0
    gap = 0.12 * top
    rows = [("data", None)] + [(m, m) for m in mets] + [("residual", None)]
    comp_h = [float(np.nanmax(tr[m][w])) if w.any() else 0.0 for m in mets]
    heights = [top] + [max(h, 0.05 * top) for h in comp_h] + [max(float(np.nanmax(np.abs(resid[w]))) * 2, 0.05 * top) if w.any() else 0.1 * top]
    offsets, ycur = [], 0.0
    for h in reversed(heights):
        offsets.append(ycur)
        ycur += h + gap
    offsets = list(reversed(offsets))
    total = ycur
    axes = _axes_row(fig, segs, [0.03, 0.15, 0.76, 0.82])
    pairs = list(zip(axes, segs))
    lw = spec["line_pt"]
    names = []
    for ax, (lo, hi, fac) in pairs:
        o = offsets[0]
        gx, go, gf = _cut(x, obs, fit, leave_out=leave)
        ax.plot(gx, go + o, color=trace_colour(spec, tr, 0, "k"), lw=lw)
        ax.plot(gx, gf + o, color=spec.get("fit_colour") or FIT_COLOUR, lw=lw * 0.9)
        for k, m in enumerate(mets):
            o = offsets[k + 1]
            gx, gc = _cut(x, tr[m], leave_out=leave)
            col = COLOURS[k % len(COLOURS)] if COLOURS[k % len(COLOURS)] != "#000000" else "#555555"
            ax.fill_between(gx, o, gc + o, color=_light(col, 0.45), lw=0)
            ax.plot(gx, gc + o, color=col, lw=lw * 0.8)
        o = offsets[-1] + heights[-1] / 2
        gx, gr = _cut(x, resid, leave_out=leave)
        ax.plot(gx, gr + o, color="#777777", lw=lw * 0.8)
        ax.set_ylim(-0.02 * total, total)
        _ticks(ax, lo, hi, spec["font_pt"])
    last = pairs[-1][0]
    smp = tr["samples"][0] if len(tr["samples"]) == 1 else None
    labs = [("measured, fit", offsets[0] + heights[0] * 0.3)]
    for k, m in enumerate(mets):
        tier, conc = R.tier(m, smp) if smp else ("", np.nan)
        val = (f" <={conc:.2g} mM" if tier == "Upper bound only" else f" {conc:.3g} mM") \
            if smp and np.isfinite(conc) else ""
        labs.append((f"{plain_name(m)}{val}", offsets[k + 1] + heights[k + 1] * 0.2))
        names.append((m, tier, conc))
    labs.append(("residual", offsets[-1] + heights[-1] * 0.3))
    for txt, y in labs:
        last.text(1.02, y, txt, transform=last.get_yaxis_transform(), ha="left", va="bottom", fontsize=spec["font_pt"],
                  clip_on=False)
    fig.text(0.41, 0.03, "\u03b4 \u00b9H (ppm)", ha="center", va="bottom", fontsize=spec["font_pt"])
    info["fit_names"] = names
    info["fit_segments"] = segs
    info["plotted"] = {"ppm": x, "measured (baseline removed)": obs, "fit (baseline removed)": fit, "residual": resid,
                       **{f"{m} (fitted)": tr[m] for m in mets}}


# =============================================================================
# Caption
# =============================================================================
def _fmt_ppm(v):
    return f"{v:.2f}"


def _trace_desc(tr):
    if tr["mean"]:
        return f"{tr['label']} (mean of {len(tr['samples'])} spectra, EXPNO {', '.join(tr['samples'])})"
    return f"{tr['label']}" + ("" if f"EXPNO {tr['samples'][0]}" in tr["label"] else f" (EXPNO {tr['samples'][0]})")


def caption(spec, R, traces, info):
    """A draft figure legend from the settings (no em dashes; edit it freely)."""
    sf = R.sf()
    fld = f"{sf:.0f} MHz " if sf else ""
    pp = R.pulse_program(traces[0]["samples"][0]) if traces else None
    pp_txt = f", pulse sequence {pp}" if pp else ""
    ref = R.reference()
    parts = []
    t = spec["template"]
    if t == "representative":
        parts.append(f"Representative {fld}1H NMR spectrum of {_trace_desc(traces[0])}{pp_txt}." +
                     (f" The fitted spectrum (sum of all metabolites fitted by MANC-Q) is drawn in "
                      f"{colour_name(spec.get('fit_colour') or FIT_COLOUR)}." if spec.get("show_fit") else ""))
    elif t == "stack":
        cols = [trace_colour(spec, tr, k, "k") for k, tr in enumerate(traces)]
        named = len({c.lower() for c in cols}) > 1 or cols[0].lower() not in ("k", "#000000", "black")
        dx = float(spec.get("offset_x_ppm") or 0.0)
        parts.append(f"{fld}1H NMR spectra{pp_txt} of " +
                     "; ".join(_trace_desc(tr) + (f" ({colour_name(c)})" if named else "") for tr, c in zip(traces, cols)) +
                     f", offset vertically by {spec['stack_offset']:g} % of the plot height" +
                     (f" and sideways by {dx:g} ppm per spectrum" if dx else "") +
                     ", bottom to top in the order listed." +
                     (" The ppm scale applies to the bottom spectrum." if dx else ""))
    elif t == "overlay":
        names = "; ".join(f"{_trace_desc(tr)} ({colour_name(trace_colour(spec, tr, k, COLOURS[k % len(COLOURS)]))})"
                          for k, tr in enumerate(traces))
        parts.append(f"{'Overlaid ' if len(traces) > 1 else ''}{fld}1H NMR spectra{pp_txt} of {names}.")
        vo, dx = float(spec.get("overlay_offset") or 0.0), float(spec.get("offset_x_ppm") or 0.0)
        if len(traces) > 1 and (vo or dx):
            parts.append("Each spectrum is offset from the one before by " +
                         " and ".join(x for x in (f"{vo:g} % of the plot height upwards" if vo else "",
                                                  f"{dx:g} ppm sideways" if dx else "") if x) +
                         (" (the ppm scale applies to the first spectrum)." if dx else "."))
    else:
        tr = traces[0]
        names = info.get("fit_names", [])
        segs = info.get("fit_segments", [])
        rng = " and ".join(f"{_fmt_ppm(lo)} to {_fmt_ppm(hi)} ppm" for lo, hi, _ in sorted(segs))
        parts.append(f"Fit of the {fld}1H NMR spectrum of {_trace_desc(tr)} by MANC-Q {q.__version__}, {rng}. Top: "
                     f"measured spectrum ({colour_name(trace_colour(spec, tr, 0, 'k'))}) and the sum of all fitted metabolite "
                     f"spectra ({colour_name(spec.get('fit_colour') or FIT_COLOUR)}), baseline "
                     f"removed. Below: the fitted spectrum of each metabolite on the same vertical scale"
                     + (" (" + "; ".join(f"{plain_name(m)}{'' if not np.isfinite(c) else f', {c:.3g} mM, {tt.lower()}'}"
                                          for m, tt, c in names) + ")" if names else "") +
                     ". Bottom: residual (measured minus fit, grey).")
    parts.append(f"Chemical shifts are referenced to {ref} at 0.00 ppm.")
    if t != "fit" and spec.get("baseline_correct", True):
        parts.append("The smooth baseline fitted by MANC-Q is subtracted.")
    if t != "fit":
        parts.append("All spectra share one vertical scale (each divided by its own " + ref + " signal)."
                     if len(traces) > 1 else f"Intensity is scaled to the {ref} signal.")
    leave = parse_ranges(spec["leave_out"])
    if leave:
        parts.append("Not shown: " + ", ".join(f"{_fmt_ppm(a)} to {_fmt_ppm(b)} ppm" for a, b in leave) +
                     (" (residual water)" if any(a <= 4.8 <= b for a, b in leave) else "") + ".")
    if info.get("magnified"):
        lo, hi, f = info["magnified"]
        parts.append(f"The region {_fmt_ppm(lo)} to {_fmt_ppm(hi)} ppm is magnified \u00d7{f:g}.")
    if info.get("inset"):
        a, b, f = info["inset"]
        parts.append(f"Inset: {_fmt_ppm(a)} to {_fmt_ppm(b)} ppm" + (f", magnified \u00d7{f:.0f}" if f >= 1.5 else "") +
                     "; the dashed box marks this region in the main plot.")
    if info.get("clipped"):
        parts.append("Peaks taller than the plot are cut off" + (" and marked ~." if t != "overlay" else "."))
    rows = [r for r in (info.get("labels_used") or []) if r.get("show")]
    syms = sorted({r.get("symbol") for r in rows if r.get("symbol")} - {""})
    if rows and spec["label_style"] != "off":
        mets = [m for m in dict.fromkeys(r["metabolite"] for r in sorted(rows, key=lambda r: -r["ppm"]) if r.get("metabolite"))]
        if spec["label_style"] == "numbers":
            nums = info.get("numbers") or {}
            key = "; ".join(f"{nums[m]}, {plain_name(m)}" for m in sorted(mets, key=lambda m: nums.get(m, 0)))
            parts.append(f"Peak assignments: {key}.")
        else:
            seen, key = set(), []
            for r in rows:
                if r.get("metabolite") and r["metabolite"] not in seen:
                    seen.add(r["metabolite"]); key.append(f"{(r.get('text') or abbrev(r['metabolite']))}, {plain_name(r['metabolite'])}")
            parts.append(f"Abbreviations: {'; '.join(key)}.")
        if spec["label_style"] != "off" and info.get("brackets_drawn"):
            parts.append("Brackets span multiplets.")
    for s in syms:
        if s in SYMBOLS:
            parts.append(f"{s} {SYMBOLS[s]}.")
    return " ".join(parts)


# =============================================================================
# Export and preflight
# =============================================================================
def figure_size(spec):
    return spec["width_mm"] / 25.4, spec["height_mm"] / 25.4


def preflight(fig):
    """(smallest text in pt, thinnest line in pt) actually drawn, so they can be checked against journal rules."""
    from matplotlib.text import Text
    from matplotlib.lines import Line2D
    from matplotlib.collections import Collection
    from matplotlib.patches import Patch
    fonts = [t.get_fontsize() for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    lines = [l.get_linewidth() for l in fig.findobj(Line2D) if l.get_visible() and l.get_linewidth() > 0
             and len(l.get_xdata()) > 1]
    for c in fig.findobj(Collection):
        lw = c.get_linewidths()
        lines += [float(v) for v in np.atleast_1d(lw) if v > 0 and c.get_visible()]
    for p in fig.findobj(Patch):
        if p.get_visible() and p.get_linewidth() > 0 and p.get_edgecolor()[3] > 0 and p is not fig.patch \
                and not getattr(p, "axes", None) is None and p is not p.axes.patch:
            lines.append(p.get_linewidth())
    return (min(fonts) if fonts else np.nan), (min(lines) if lines else np.nan)


def preflight_text(fig, spec):
    f, l = preflight(fig)
    ok_f, ok_l = (not np.isfinite(f) or f >= MIN_FONT_PT), (not np.isfinite(l) or l >= MIN_LINE_PT)
    msg = f"At {spec['width_mm']:g} x {spec['height_mm']:g} mm: smallest text {f:.1f} pt, thinnest line {l:.2f} pt."
    if ok_f and ok_l:
        return msg + " OK for most journals.", True
    adv = []
    if not ok_f:
        adv.append(f"text below {MIN_FONT_PT:g} pt is too small for most journals: raise Font size")
    if not ok_l:
        adv.append(f"lines thinner than {MIN_LINE_PT:g} pt may vanish in print: raise Line width")
    return msg + " " + "; ".join(adv) + ".", False


def write_tiff(path, fig, dpi=1200):
    """Save fig as a Deflate-compressed TIFF at dpi, written here rather than by Pillow: some Pillow builds (for
    example in Anaconda) crash on any compressed TIFF. Greyscale when the figure has no colour."""
    import struct
    import zlib
    old = fig.dpi
    try:
        fig.set_dpi(dpi)
        cv = FigureCanvasAgg(fig)
        cv.draw()
        buf, (w, h) = cv.print_to_buffer()
    finally:
        fig.set_dpi(old)
    rgba = np.frombuffer(buf, np.uint8).reshape(h, w, 4)
    a = rgba[..., 3:4].astype(np.uint16)
    img = ((rgba[..., :3].astype(np.uint16) * a + 255 * (255 - a)) // 255).astype(np.uint8)   # onto white
    grey = bool(np.array_equal(img[..., 0], img[..., 1]) and np.array_equal(img[..., 1], img[..., 2]))
    data = np.ascontiguousarray(img[..., 0] if grey else img)
    spp = 1 if grey else 3
    rows = max(1, 262144 // (w * spp))
    strips = [zlib.compress(data[r:r + rows].tobytes(), 6) for r in range(0, h, rows)]
    extra = bytearray()                      # values that do not fit in a tag: placed after the IFD
    tags = []
    n_tags = 13
    ifd_at = 8 + sum(len(x) for x in strips)
    ifd_at += ifd_at % 2
    extra_at = ifd_at + 2 + 12 * n_tags + 4

    def tag(code, typ, values):
        fmt = {3: "H", 4: "I", 5: "II"}[typ]
        size = {3: 2, 4: 4, 5: 8}[typ]
        packed = b"".join(struct.pack("<" + fmt, *(v if isinstance(v, tuple) else (v,))) for v in values)
        if len(packed) <= 4:
            tags.append(struct.pack("<HHI", code, typ, len(values)) + packed.ljust(4, b"\0"))
        else:
            off = extra_at + len(extra)
            extra.extend(packed + (b"\0" if len(packed) % 2 else b""))
            tags.append(struct.pack("<HHII", code, typ, len(values), off))
    offsets, pos = [], 8
    for x in strips:
        offsets.append(pos); pos += len(x)
    tag(256, 4, [w]); tag(257, 4, [h]); tag(258, 3, [8] * spp); tag(259, 3, [8])
    tag(262, 3, [1 if grey else 2]); tag(273, 4, offsets); tag(277, 3, [spp]); tag(278, 4, [rows])
    tag(279, 4, [len(x) for x in strips]); tag(282, 5, [(int(dpi), 1)]); tag(283, 5, [(int(dpi), 1)])
    tag(284, 3, [1]); tag(296, 3, [2])
    with open(path, "wb") as f:
        f.write(b"II*\0" + struct.pack("<I", ifd_at))
        for x in strips:
            f.write(x)
        if f.tell() % 2:
            f.write(b"\0")
        f.write(struct.pack("<H", n_tags) + b"".join(tags) + struct.pack("<I", 0))
        f.write(bytes(extra))


def render(spec, R, traces, dpi=100):
    """A Figure at the final printed size, drawn."""
    fig = Figure(figsize=figure_size(spec), dpi=dpi)
    FigureCanvasAgg(fig)
    info = draw(fig, spec, R, traces)
    return fig, info


def collect_traces(R, spec, sheet):
    groups = trace_groups(R, sheet)
    comps = spec.get("fit_metabolites") or []
    out = []
    for key in spec.get("traces") or []:
        if key in groups:
            tr = build_trace(R, groups[key], spec.get("mode", "mean"), spec.get("replicate"), comps)
            if tr is not None:
                tr["key"] = key
                out.append(tr)
    return out


def export(spec, out, base, formats=("pdf",), extras=(), caption_text=None, sheet=None):
    """Save the figure described by spec (results folder out) as base.<format> for each format ('pdf', 'svg', 'eps',
    'tiff' at 1200 dpi, 'png' at 600 dpi), plus extras: 'caption' (base_caption.txt), 'data' (base_data.csv),
    'record' (base_processing.txt), 'script' (base_redraw.py with base_figure.json). Returns the files written."""
    R = Results(out)
    sheet = sheet if sheet is not None else load_sheet(out, R.samples, R.titles)
    traces = collect_traces(R, spec, sheet)
    written = []
    os.makedirs(os.path.dirname(os.path.abspath(base)), exist_ok=True)
    with rc_context(_style(spec)):
        fig, info = render(spec, R, traces, dpi=100)
        for fmt in formats:
            f = f"{base}.{'tif' if fmt == 'tiff' else fmt}"
            if fmt == "tiff":
                write_tiff(f, fig, 1200)
            else:
                fig.savefig(f, format=fmt, facecolor="white", **(dict(dpi=600) if fmt == "png" else {}))
            written.append(f)
    if "caption" in extras:
        f = f"{base}_caption.txt"
        with open(f, "w", encoding="utf-8") as fh:
            fh.write((caption_text or info.get("caption", "")) + "\n")
        written.append(f)
    if "data" in extras and info.get("plotted"):
        f = f"{base}_data.csv"
        pl = info["plotted"]
        x = pl["ppm"]
        pd.DataFrame({k: (np.interp(x, x, v) if len(v) == len(x) else v) for k, v in pl.items()}).to_csv(
            f, index=False, float_format="%.6g")
        written.append(f)
    if "record" in extras:
        f = f"{base}_processing.txt"
        lines = [f"Figure made with MANC-Q {q.__version__} on {time.strftime('%Y-%m-%d %H:%M')}",
                 f"Results folder: {out}", f"Template: {TEMPLATES.get(spec['template'], spec['template'])}",
                 f"Intensity: each spectrum divided by its own {R.reference()} area (from the fit)", ""]
        for tr in traces:
            lines.append(f"Trace '{tr['label']}': EXPNO {', '.join(tr['samples'])}" + (" (mean)" if tr["mean"] else ""))
            for s in tr["samples"]:
                qrow = R.qc.loc[s] if s in R.qc.index else {}
                lines.append(f"  EXPNO {s}: title '{R.titles.get(s, '')}'; "
                             f"{qrow.get('processing', 'processing not recorded')}; "
                             f"{qrow.get('SF_MHz', float('nan')):.2f} MHz; {R.reference()} linewidth "
                             f"{qrow.get('TSP_fwhm_Hz', float('nan')):.2f} Hz; pulse sequence {R.pulse_program(s) or 'unknown'}")
        lines += ["", "Figure settings:", json.dumps(spec, indent=1, default=str)]
        with open(f, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
        written.append(f)
    if "script" in extras:
        fj, fp = f"{base}_figure.json", f"{base}_redraw.py"
        with open(fj, "w", encoding="utf-8") as fh:
            json.dump(dict(spec=spec, results=os.path.abspath(out),
                           sheet=sheet.to_dict(orient="list")), fh, indent=1, default=str)
        with open(fp, "w", encoding="utf-8") as fh:
            fh.write(
            '"""Draws this figure again with MANC-Q (run it with the Python MANC-Q is installed in):\n'
            f'    python "{os.path.basename(fp)}"\n'
            'Change the settings in the .json file next to it (width_mm, font_pt, labels, ...) and run it again.\n'
            'It needs the MANC-Q results folder named in the .json file."""\n'
            "import json, os\nimport pandas as pd\nfrom mancq import figures\n\n"
            "here = os.path.dirname(os.path.abspath(__file__))\n"
            f"with open(os.path.join(here, {os.path.basename(fj)!r}), encoding='utf-8') as fh:\n"
            f"    cfg = json.load(fh)\n"
            f"base = os.path.join(here, {os.path.basename(base)!r})\n"
            "figures.export(cfg['spec'], cfg['results'], base, formats=('pdf', 'png'),\n"
            "               sheet=pd.DataFrame(cfg['sheet']).fillna(''))\n"
            "print('Saved', base + '.pdf and .png')\n")
        written += [fj, fp]
    return written
