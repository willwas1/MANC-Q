"""mancq.twod - 2D 1H-1H spectra (TOCSY, COSY) shown next to a 1D fit, used by the Results step of the window.

A 2D spectrum does not give concentrations here (TOCSY and COSY peak volumes depend on mixing, relaxation and
couplings), so it is used as a check of identity: for a fitted compound, MANC-Q predicts where its cross peaks
should be (every pair of protons in one coupled spin system of the library, at the shifts fitted in the 1D) and
looks for signal there. Cross peaks that are present support the assignment; strong compounds whose cross peaks
are missing deserve a second look. Compounds from fixed peak lists (no spin system) and singlets have none.

Each 1D sample is linked to a 2D experiment in the same data folder: the nearest EXPNO with the same title, or one
the user chooses. Links chosen by hand are kept in <results>/_per_sample/twod_links.json (plain JSON).
"""
import os
import json
import numpy as np
import pandas as pd
from . import engine as q

MIN_J_HZ = 2.0          # couplings smaller than this do not carry TOCSY transfer in practice
MIN_SEP_PPM = 0.06      # protons closer than this give cross peaks lost in the diagonal
BOX_PPM = 0.015         # half width of the box searched for each cross peak (both axes)
SEEN_SNR = 10.0         # a cross peak counts as seen when a peak top there is this many times the noise


def _jcamp(f):
    out = {}
    for line in open(f, errors="replace"):
        if line.startswith("##$"):
            k, _, v = line[3:].partition("=")
            v = v.strip()
            try:
                out[k.strip()] = float(v)
            except ValueError:
                out[k.strip()] = v.strip("<>")
    return out


def is_2d(expdir, procno=1):
    return os.path.exists(os.path.join(expdir, "pdata", str(procno), "2rr"))


def title(expdir, procno=1):
    for p in (os.path.join(expdir, "pdata", str(procno), "title"), os.path.join(expdir, "pdata", "1", "title")):
        if os.path.exists(p):
            return open(p, errors="replace").read().strip()
    return ""


def experiment(expdir):
    """Pulse programme name of an experiment, or ''."""
    f = os.path.join(expdir, "acqus")
    return str(_jcamp(f).get("PULPROG", "")) if os.path.exists(f) else ""


def kind(expdir):
    """'TOCSY', 'COSY', 'NOESY', 'HSQC' ... from the pulse programme and nuclei, for labels."""
    p = experiment(expdir).lower()
    nuc = ""
    f = os.path.join(expdir, "acqu2s")
    if os.path.exists(f):
        nuc = str(_jcamp(f).get("NUC1", ""))
    if nuc and nuc != "1H" and os.path.exists(os.path.join(expdir, "acqus")):
        n1 = str(_jcamp(os.path.join(expdir, "acqus")).get("NUC1", ""))
        if n1 != nuc:
            return "HSQC" if "hsqc" in p else f"{n1}-{nuc}"
    for key, name in (("mlev", "TOCSY"), ("dipsi", "TOCSY"), ("tocsy", "TOCSY"), ("cosy", "COSY"),
                      ("noesy", "NOESY"), ("roesy", "ROESY"), ("jres", "J-resolved")):
        if key in p:
            return name
    return "2D"


def homonuclear(expdir):
    k = kind(expdir)
    return k in ("TOCSY", "COSY", "NOESY", "ROESY", "2D")


def find_2d(datadir, procno=1):
    """EXPNO folders under datadir holding a processed 2D spectrum."""
    if not datadir or not os.path.isdir(datadir):
        return []
    return [d for d in sorted(os.listdir(datadir), key=lambda x: (len(x), x))
            if os.path.isdir(os.path.join(datadir, d)) and is_2d(os.path.join(datadir, d), procno)]


def guess_partner(expdir, procno=1):
    """The 2D experiment that belongs to a 1D one: in the same folder, homonuclear, same title, nearest EXPNO
    (later ones first, as 2D spectra are usually recorded after the 1D). None if there is none."""
    parent, name = os.path.split(os.path.normpath(expdir))
    cands = [d for d in find_2d(parent, procno) if homonuclear(os.path.join(parent, d))]
    if not cands:
        return None
    t = title(expdir, procno)
    same = [d for d in cands if t and title(os.path.join(parent, d), procno) == t]
    if not same and name.isdigit():
        # titles differ (often a typing slip): the next 2D within a few EXPNOs, which the window flags
        same = [d for d in cands if d.isdigit() and 0 < int(d) - int(name) <= 5]
    if not same:
        return None
    if name.isdigit():
        e = int(name)
        same.sort(key=lambda d: (abs(int(d) - e), int(d) < e) if d.isdigit() else (1e9, True))
    return os.path.join(parent, same[0])


def read_2rr(expdir, procno=1):
    """dict(z=2D array [F1 rows, F2 columns], f1=ppm, f2=ppm (both descending, as stored), path, kind, title)."""
    pdir = os.path.join(expdir, "pdata", str(procno))
    p2, p1 = _jcamp(os.path.join(pdir, "procs")), _jcamp(os.path.join(pdir, "proc2s"))
    si2, si1 = int(p2["SI"]), int(p1["SI"])
    xd2, xd1 = int(p2.get("XDIM", si2) or si2), int(p1.get("XDIM", si1) or si1)
    dt = ("<" if int(p2.get("BYTORDP", 0)) == 0 else ">") + ("i4" if int(p2.get("DTYPP", 0)) == 0 else "f8")
    raw = np.fromfile(os.path.join(pdir, "2rr"), dtype=dt).astype(float)
    if raw.size != si1 * si2:
        raise ValueError(f"{pdir}/2rr has {raw.size} points, expected {si1} x {si2}")
    # TopSpin stores 2D data in submatrices of XDIM(F1) x XDIM(F2) points
    n1, n2 = si1 // xd1, si2 // xd2
    z = raw.reshape(n1, n2, xd1, xd2).transpose(0, 2, 1, 3).reshape(si1, si2)
    z *= 2.0 ** float(p2.get("NC_proc", 0))

    def axis(p, n):
        return float(p["OFFSET"]) - np.arange(n) * float(p["SW_p"]) / float(p["SF"]) / n
    return dict(z=z, f2=axis(p2, si2), f1=axis(p1, si1), path=expdir, kind=kind(expdir),
                title=title(expdir, procno), sf=float(p2["SF"]))


def prepare(d):
    """Noise level and reference correction of a 2D spectrum read by read_2rr (adds noise, shift)."""
    z = d["z"]
    # noise: robust spread of a region with no signal (above 9.5 ppm in both dimensions, or the quietest strip)
    sel2, sel1 = d["f2"] > 9.5, d["f1"] > 9.5
    blk = z[np.ix_(sel1, sel2)] if sel1.sum() > 8 and sel2.sum() > 8 else z[: max(8, z.shape[0] // 16)]
    d["noise"] = float(1.4826 * np.median(np.abs(blk - np.median(blk)))) or float(np.std(z)) or 1.0
    # reference: the TSP/DSS diagonal peak should be at 0 ppm (the 1D fits are referenced to it)
    d["shift"] = 0.0
    w1 = np.where(np.abs(d["f1"]) < 0.3)[0]
    if w1.size > 2:
        cols = [int(np.argmin(np.abs(d["f2"] - d["f1"][i]))) for i in w1]
        diag = np.array([z[i, j] for i, j in zip(w1, cols)])
        k = int(np.argmax(diag))
        if diag[k] > 20 * d["noise"]:
            i, j = w1[k], cols[k]
            sub = z[max(0, i - 2):i + 3, max(0, j - 2):j + 3]
            a, b = np.unravel_index(np.argmax(sub), sub.shape)
            d["shift"] = float(-(d["f1"][max(0, i - 2) + a] + d["f2"][max(0, j - 2) + b]) / 2)
    d["f1c"], d["f2c"] = d["f1"] + d["shift"], d["f2"] + d["shift"]
    return d


def load(expdir, procno=1):
    return prepare(read_2rr(expdir, procno))


# ---------------------------------------------------------------- which 2D belongs to which sample
def _links_file(out):
    return os.path.join(out, "_per_sample", "twod_links.json")


def read_links(out):
    f = _links_file(out)
    if os.path.exists(f):
        try:
            return dict(json.load(open(f)).get("links", {}))
        except (ValueError, OSError):
            return {}
    return {}


def save_link(out, smp, path):
    """Remember the 2D chosen for a sample: a folder, "" for none, or None to go back to the automatic choice."""
    links = read_links(out)
    if path is None:
        links.pop(str(smp), None)
    else:
        links[str(smp)] = path
    os.makedirs(os.path.dirname(_links_file(out)), exist_ok=True)
    json.dump(dict(about="2D spectrum shown for each sample in the MANC-Q Results step (chosen by hand); "
                         "\"\" = none", links=links), open(_links_file(out), "w"), indent=1)


def sample_dir(out, smp):
    """The EXPNO folder a sample was read from, from the run settings saved with the results."""
    f = os.path.join(out, "_per_sample", "run_settings.json")
    if not os.path.exists(f):
        return None
    try:
        s = json.load(open(f))
    except (ValueError, OSError):
        return None
    # the original folder first: spectra MANC-Q processed itself (from a FID) sit in <results>/_processed,
    # away from the 2D experiments
    base = s.get("DATA_DIR")
    if base:
        d = os.path.join(base, str(smp))
        if os.path.isdir(d):
            return d
        if os.path.basename(os.path.normpath(base)) == str(smp) and os.path.isdir(base):
            return base
    d = (s.get("SAMPLE_DIRS") or {}).get(str(smp))
    if isinstance(d, (list, tuple)):
        d = d[0] if d else None
    return d if isinstance(d, str) and os.path.isdir(d) else None


def partner(out, smp, procno=1):
    """(2D folder or None, how it was found: 'chosen', 'automatic', 'automatic, titles differ', 'none' or
    'no data')."""
    links = read_links(out)
    if str(smp) in links:
        p = links[str(smp)]
        return (p, "chosen") if p and is_2d(p, procno) else (None, "none")
    d = sample_dir(out, smp)
    if d is None:
        return None, "no data"
    p = guess_partner(d, procno)
    if p is None:
        return None, "none"
    return p, "automatic" if title(p, procno) == title(d, procno) else "automatic, titles differ"


# ---------------------------------------------------------------- expected cross peaks
def _spin_systems(libdir=None):
    libdir = libdir or os.path.join(q.HERE, "library")
    g = q.parse_gissmo(os.path.join(libdir, "gissmo_spin_systems.txt"))
    o = q.parse_gissmo(os.path.join(libdir, "other_spin_systems.txt"))
    comp = pd.read_csv(os.path.join(libdir, "compounds.csv"))
    found = {}
    for _, r in comp.iterrows():
        e = g.get(r.source_id) if r.source == "GISSMO" else o.get(r.source_id) if r.source == "SPIN" else None
        if e is not None:
            found[r.metabolite] = e
    return found


_SPIN = None


def spin_system(name):
    global _SPIN
    if _SPIN is None:
        _SPIN = _spin_systems()
    return _SPIN.get(name)


def _groups(n, J, cosy):
    """Pairs of protons that give cross peaks: directly coupled (COSY) or in one coupled network (TOCSY)."""
    adj = {i: set() for i in range(n)}
    for (i, j), v in J.items():
        if abs(v) >= MIN_J_HZ and i < n and j < n:
            adj[i].add(j); adj[j].add(i)
    if cosy:
        return {(min(i, j), max(i, j)) for i in adj for j in adj[i]}
    seen, pairs = set(), set()
    for i in range(n):
        if i in seen:
            continue
        comp, stack = set(), [i]
        while stack:
            k = stack.pop()
            if k not in comp:
                comp.add(k); stack.extend(adj[k] - comp)
        seen |= comp
        comp = sorted(comp)
        pairs |= {(a, b) for ai, a in enumerate(comp) for b in comp[ai + 1:]}
    return pairs


def expected_peaks(name, mults, sf, cosy=False, offset_hz=0.0):
    """Cross peaks a fitted compound should give: [(ppm a, ppm b)] with a > b, merged when closer than BOX_PPM.
    mults: [(fitted centre ppm, shift from library Hz)] of its multiplets, as in the peak table of the results;
    offset_hz: the library-vs-reference offset of the sample (GISSMO_offset_Hz in the QC table). Each proton is
    moved by the shift of the multiplet nearest to it in the library. [] for compounds without a spin system or
    couplings."""
    e = spin_system(name)
    if e is None:
        return []
    cs = np.asarray(e["shifts"], float)
    if mults:
        now = np.array([float(c) for c, _ in mults])
        lib_c = now - (np.array([float(h) for _, h in mults]) + float(offset_hz)) / sf
        cs = np.array([c + (now - lib_c)[int(np.argmin(np.abs(lib_c - c)))] for c in cs])
    out = []
    for i, j in sorted(_groups(len(cs), e["J"], cosy)):
        a, b = max(cs[i], cs[j]), min(cs[i], cs[j])
        if a - b < MIN_SEP_PPM:
            continue
        if not any(abs(a - p[0]) < BOX_PPM and abs(b - p[1]) < BOX_PPM for p in out):
            out.append((float(a), float(b)))
    return out


def peak_table(out, smp):
    """The per-multiplet table of one sample (fitted centres and shifts), or None."""
    f = os.path.join(out, "_per_sample", f"{smp}_peaks.csv")
    return pd.read_csv(f) if os.path.exists(f) else None


def expected_from_table(P, name, sf, cosy=False, offset_hz=0.0):
    sub = P[P.metabolite == name] if P is not None else []
    mults = [(r.centre_ppm, r.shift_from_library_Hz) for r in sub.itertuples()] if len(sub) else []
    return expected_peaks(name, mults, sf, cosy, offset_hz)


def box_max(d, x, y, half=BOX_PPM):
    """Largest 2D intensity within +-half ppm of (F2 = x, F1 = y)."""
    w2 = np.abs(d["f2c"] - x) <= half
    w1 = np.abs(d["f1c"] - y) <= half
    if not w2.any() or not w1.any():
        return np.nan
    return float(d["z"][np.ix_(w1, w2)].max())


def apex(d, x, y, half=BOX_PPM):
    """Height (in noise units) of a real peak top within +-half ppm of (F2 = x, F1 = y), or 0. A peak top is a
    point higher than its 8 neighbours; the slope of a bigger peak next to the box does not count."""
    z = d["z"]
    j = np.where(np.abs(d["f2c"] - x) <= half)[0]
    i = np.where(np.abs(d["f1c"] - y) <= half)[0]
    if not i.size or not j.size:
        return 0.0
    i0, i1 = max(i[0] - 1, 0), min(i[-1] + 2, z.shape[0])
    j0, j1 = max(j[0] - 1, 0), min(j[-1] + 2, z.shape[1])
    sub = z[i0:i1, j0:j1]
    best = 0.0
    for a in range(1, sub.shape[0] - 1):
        for b in range(1, sub.shape[1] - 1):
            v = sub[a, b]
            if v > best and v >= sub[a - 1:a + 2, b - 1:b + 2].max():
                best = v
    return float(best / d["noise"])


def check(d, peaks):
    """[(a, b, height)] for each expected cross peak: the height of the peak top found at the stronger of its two
    symmetric positions, in noise units (0 = no peak top there). Seen when height >= SEEN_SNR."""
    return [(a, b, max(apex(d, a, b), apex(d, b, a))) for a, b in peaks]


def summary(checked):
    """'4 of 5 expected cross peaks seen' style text, and the fraction seen (None if nothing was expected)."""
    if not checked:
        return "no cross peaks expected (singlets, or no spin system in the library)", None
    n = sum(s >= SEEN_SNR for _, _, s in checked)
    return f"{n} of {len(checked)} expected cross peak{'s' if len(checked) != 1 else ''} seen", n / len(checked)


def check_sample(d, P, G, smp, sf, offset_hz=0.0):
    """Table of every metabolite with a value in one sample against its 2D spectrum d: expected and seen cross
    peaks. P: the sample's peak table; G: the long results table; offset_hz as for expected_peaks."""
    G = G[G["sample"] == str(smp)]
    cosy = d["kind"] == "COSY"
    rows = []
    for r in G.itertuples():
        if r.tier in ("Not detected", "Not measurable (region ignored)") or r.metabolite.startswith("Pluronic"):
            continue
        chk = check(d, expected_from_table(P, r.metabolite, sf, cosy, offset_hz))
        rows.append(dict(sample=str(smp), metabolite=r.metabolite, conc_mM=r.conc_mM, tier=r.tier,
                         expected=len(chk), seen=sum(s >= SEEN_SNR for *_, s in chk),
                         cross_peaks="; ".join(f"{a:.3f}/{b:.3f}: {s:.0f}" for a, b, s in chk)))
    return pd.DataFrame(rows, columns=["sample", "metabolite", "conc_mM", "tier", "expected", "seen", "cross_peaks"])
