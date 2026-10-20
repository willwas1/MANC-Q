"""
mancq.library - the compound library: GISSMO and other spin systems simulated at the spectrometer
frequency, CASMDB peak lists, multiplet grouping, Pluronic terms, and the library fingerprint.
"""
import hashlib
import json
import os
import numpy as np
import pandas as pd
from .spinsim import merge_sticks, simulate_spin_system
from .settings import __version__


# =============================================================================
# 3. Library
# =============================================================================
def parse_gissmo(path):
    lib = {}
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    for line in text.splitlines():
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        gid, name, rest = line.split("|", 2)
        p = rest.split("#")
        cs = [float(x) for x in p[1].split(",")]
        J = {}
        for x in p[2].split(";"):
            if x:
                a, v = x.split(":")
                i, j = a.split("-")
                J[(int(i) - 1, int(j) - 1)] = float(v)
        het = None
        if len(p) > 3 and p[3].startswith("ADD"):
            het = {int(k) - 1: float(v) for k, v in (e.split("=") for e in p[3][3:].split(","))}
        lib[gid] = dict(name=name, shifts=cs, J=J, hetero=het)
    return lib


def build_library(sf, libdir, cache_dir, *, S):
    key = f"library_{sf:.3f}MHz.json"
    cpath = os.path.join(cache_dir, key)
    comp = pd.read_csv(os.path.join(libdir, "compounds.csv"))
    text = ""
    for fn in ("compounds.csv", "gissmo_spin_systems.txt", "other_spin_systems.txt", "peaklists.csv"):
        with open(os.path.join(libdir, fn), encoding="utf-8") as fh:
            text += fh.read()
    # the H1 shift of each glucose anomer in the spin-system file, and its share
    anomers = {"D-Glucose": {5.2223: S["GLUCOSE_ALPHA_FRACTION"], 4.6360: 1.0 - S["GLUCOSE_ALPHA_FRACTION"]}}
    used = {k: S[k] for k in ("PEAKLIST_FIELD_MHZ", "GLUCOSE_ALPHA_FRACTION")}
    sig = hashlib.md5((text + __version__ + json.dumps(used, sort_keys=True)).encode()).hexdigest()
    if os.path.exists(cpath):
        with open(cpath, encoding="utf-8") as fh:
            d = json.load(fh)
        if d.get("sig") == sig:
            return d["compounds"]
    g = parse_gissmo(os.path.join(libdir, "gissmo_spin_systems.txt"))
    other = parse_gissmo(os.path.join(libdir, "other_spin_systems.txt"))
    pk = pd.read_csv(os.path.join(libdir, "peaklists.csv"))
    man = pd.read_csv(os.path.join(libdir, "peaklists_manifest.csv")).set_index("metabolite")
    out = []
    for _, r in comp.iterrows():
        if r.source in ("GISSMO", "SPIN"):
            e = g[r.source_id] if r.source == "GISSMO" else other[r.source_id]
            ppm, it, grp = simulate_spin_system(e["shifts"], e["J"], sf, e["hetero"], anomers.get(r.metabolite))
            ppm, it, grp = merge_sticks(ppm, it, grp)
            keep = it > 2e-3 * it.max()
            ppm, it, grp = ppm[keep], it[keep], grp[keep]
            nprot = len(e["shifts"]) if r.metabolite not in anomers else len(e["shifts"]) / 2
            it = it * nprot / it.sum()
            prov = (f"GISSMO bmse{r.source_id} ({e['name']}), QM-simulated at {sf:.2f} MHz" if r.source == "GISSMO" else
                    f"spin system {r.source_id} ({e['name']}, see other_spin_systems.txt), QM-simulated at {sf:.2f} MHz")
        else:
            sub = pk[pk.metabolite == r.source_id]
            ppm = sub.ppm.values * 1.0
            it = sub.height.values * 1.0
            grp = sub.cluster.values.astype(int)
            # peak lists are defined at PEAKLIST_FIELD_MHZ: rescale the line spacing within each cluster
            # (first-order approximation: J in Hz is field independent)
            scale = S["PEAKLIST_FIELD_MHZ"] / sf
            for g_ in np.unique(grp):
                m_ = grp == g_
                c_ = (ppm[m_] * it[m_]).sum() / it[m_].sum()
                ppm[m_] = c_ + (ppm[m_] - c_) * scale
            nprot = float(man.loc[r.source_id, "n_protons"])
            it = it * nprot / it.sum()
            prov = (f"CASMDB {man.loc[r.source_id,'casmdb_name']} peak list (origin {man.loc[r.source_id,'origin']}), "
                    f"spacing rescaled to {sf:.2f} MHz")
        ov = {}
        for col in ("max_shift_hz", "multiplet_shift_hz"):
            if col in r and pd.notna(r[col]) and str(r[col]).strip() != "":
                ov[col] = float(r[col])
        out.append(dict(name=r.metabolite, cls=r["class"], source=r.source, titratable=int(r.titratable), **ov,
                        nprot=float(nprot), ppm=[float(v) for v in ppm], inten=[float(v) for v in it], group=[int(v) for v in grp],
                        provenance=prov))
    os.makedirs(cache_dir, exist_ok=True)
    with open(cpath, "w", encoding="utf-8") as fh:
        json.dump(dict(sig=sig, compounds=out), fh)
    return out


def prepare_library(lib, sf, pluronic=False, *, S):
    """Attach multiplet groupings and, if requested, add the Pluronic F-68 polymer terms
    (a PEO/PPO surfactant common in cell-culture media)."""
    lib = [dict(c) for c in lib]
    jch = 142.0 / sf
    for tag, wf0, wlo, whi in ([("sharp", 1.3, 0.8, 5.0), ("broad", 25.0, 10.0, 60.0)] if pluronic else []):
        lib.append(dict(name=f"Pluronic F-68 PEO ({tag})", cls="polymer", source="polymer", titratable=0, nprot=1.0,
                        ppm=[3.7128 - jch / 2, 3.7128, 3.7128 + jch / 2], inten=[0.0055, 1.0, 0.0055], group=[1, 1, 1],
                        provenance="free polymer term (PEO -CH2CH2O-), position/width fitted, with 13C satellites",
                        wf0=wf0, wlo=wlo, whi=whi, dmax_hz=3.0))
    j = 6.3 / sf
    if pluronic:
        lib.append(dict(name="Pluronic F-68 PPO (CH3)", cls="polymer", source="polymer", titratable=0, nprot=1.0,
                        ppm=[1.150 - j / 2, 1.150 + j / 2], inten=[0.5, 0.5], group=[1, 1],
                        provenance="free polymer term (PPO methyl doublet, broad), position/width fitted",
                        wf0=6.0, wlo=2.0, whi=25.0, dmax_hz=30.0))
    for c in lib:
        if c["source"] == "polymer":
            c["mults"] = [dict(ppm=np.array(c["ppm"]), inten=np.array(c["inten"]), lo=min(c["ppm"]), hi=max(c["ppm"]),
                               centre=float(np.mean(c["ppm"])), nprot=float(sum(c["inten"])))]
        else:
            c["mults"] = group_multiplets(c["ppm"], c["inten"], c["group"], sf, S["MULTIPLET_GAP_HZ"],
                                          S["MULTIPLET_MAX_SPAN_PPM"])
    return lib


def group_multiplets(ppm, inten, group, sf, gap_hz, max_span_ppm):
    """One multiplet per proton group (from the spin simulation / CASMDB cluster);
    a group wider than max_span_ppm is split at its largest gap."""
    ppm, inten, group = np.asarray(ppm, float), np.asarray(inten, float), np.asarray(group)
    out = []
    for g in np.unique(group):
        sel = group == g
        if inten[sel].sum() <= 0:
            continue
        out.extend(split_multiplets(ppm[sel], inten[sel], sf, gap_hz=gap_hz, max_span_ppm=max_span_ppm))
    out.sort(key=lambda m: m["centre"])
    return out


def split_multiplets(ppm, inten, sf, gap_hz, max_span_ppm):
    """Group stick lines into multiplets: lines closer than gap_hz belong together;
    clusters wider than max_span_ppm are split at their largest internal gap
    (so each multiplet gets its own small shift/width freedom)."""
    o = np.argsort(ppm)
    p, it = np.asarray(ppm, float)[o], np.asarray(inten, float)[o]
    gap = gap_hz / sf
    groups, cur = [], [0]
    for k in range(1, len(p)):
        if p[k] - p[k - 1] > gap:
            groups.append(cur)
            cur = []
        cur.append(k)
    groups.append(cur)

    def split(g):
        if len(g) < 2 or p[g[-1]] - p[g[0]] <= max_span_ppm:
            return [g]
        d = np.diff(p[g])
        cut = int(np.argmax(d)) + 1
        return split(g[:cut]) + split(g[cut:])

    out = []
    for g in groups:
        for gg in split(g):
            pp, ii = p[gg], it[gg]
            out.append(dict(ppm=pp, inten=ii, lo=float(pp.min()), hi=float(pp.max()),
                            centre=float((pp * ii).sum() / ii.sum()), nprot=float(ii.sum())))
    return out


def library_hash(libdir):
    """Fingerprint of every file in the library folder (a changed library means a saved fit is out of date)."""
    h = hashlib.md5()
    for fn in sorted(os.listdir(libdir)):
        f = os.path.join(libdir, fn)
        if os.path.isfile(f):
            h.update(fn.encode())
            with open(f, "rb") as fh:
                h.update(fh.read())
    return h.hexdigest()
