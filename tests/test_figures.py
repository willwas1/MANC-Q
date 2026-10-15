"""Checks of the figure tool (mancq.figures) on a small synthetic results folder only.
Run with:  python -m pytest tests   (or:  python tests/test_figures.py)"""
import os
import sys
import struct
import tempfile
import zlib
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mancq import figures as F


def fake_results(folder):
    """Three synthetic spectra (a 'lactate' doublet and quartet, a 'formate' singlet) as MANC-Q would write them."""
    x = np.arange(0.6, 9.6, 0.0005)
    x = x[(x < 4.6) | (x > 5.0)]

    def lor(c, h, w=0.0015):
        return h * w ** 2 / ((x - c) ** 2 + w ** 2)
    rows, peaks, qc = [], [], []
    for k, s in enumerate(["10", "11", "12"]):
        lac = (lor(1.325, 1.0) + lor(1.313, 1.0) + sum(lor(c, 0.33) for c in (4.087, 4.099, 4.111, 4.123))) * (1 + k)
        frm = lor(8.45, 0.2)
        d = pd.DataFrame({"ppm": x, "observed": lac + frm + 0.01, "simulated_total": lac + frm + 0.01,
                          "baseline": np.full(len(x), 0.01), "L-Lactate": lac, "Formate": frm})
        d.to_csv(os.path.join(folder, f"simulated_spectrum_sample_{s}.csv.gz"), index=False)
        rows += [dict(sample=s, metabolite="L-Lactate", tier="Quantified", conc_mM=1.0 + k, reporter_ppm="4.1050"),
                 dict(sample=s, metabolite="Formate", tier="Quantified", conc_mM=0.1, reporter_ppm="8.4500")]
        peaks += [dict(sample=s, metabolite="L-Lactate", multiplet=1, centre_ppm=1.319, from_ppm=1.313, to_ppm=1.325,
                       height_snr=500, used_as_reporter=False),
                  dict(sample=s, metabolite="L-Lactate", multiplet=2, centre_ppm=4.105, from_ppm=4.087, to_ppm=4.123,
                       height_snr=150, used_as_reporter=True),
                  dict(sample=s, metabolite="Formate", multiplet=1, centre_ppm=8.45, from_ppm=8.45, to_ppm=8.45,
                       height_snr=90, used_as_reporter=True)]
        qc.append(dict(sample=s, title=f"synthetic {s}", SF_MHz=600.13, TSP_area=1.0, TSP_fwhm_Hz=1.0))
    pd.DataFrame(rows).to_csv(os.path.join(folder, "concentrations_long.csv"), index=False)
    pd.DataFrame(peaks).to_csv(os.path.join(folder, "peak_list_by_metabolite.csv"), index=False)
    pd.DataFrame(qc).to_csv(os.path.join(folder, "qc_per_sample.csv"), index=False)


def read_tiff(path):
    """Minimal reader for the TIFFs write_tiff makes: (width, height, samples per pixel, dpi, pixels)."""
    b = open(path, "rb").read()
    assert b[:4] == b"II*\0"
    ifd = struct.unpack("<I", b[4:8])[0]
    n = struct.unpack("<H", b[ifd:ifd + 2])[0]
    tags = {}
    for i in range(n):
        code, typ, cnt, val = struct.unpack("<HHII", b[ifd + 2 + 12 * i: ifd + 14 + 12 * i])
        size = {3: 2, 4: 4, 5: 8}[typ] * cnt
        raw = b[ifd + 10 + 12 * i: ifd + 14 + 12 * i] if size <= 4 else b[val:val + size]
        fmt = {3: "H", 4: "I", 5: "I"}[typ]
        tags[code] = struct.unpack("<" + fmt * (cnt * (2 if typ == 5 else 1)), raw[:size])
    w, h, spp = tags[256][0], tags[257][0], tags[277][0]
    data = b"".join(zlib.decompress(b[o:o + c]) for o, c in zip(tags[273], tags[279]))
    assert tags[259][0] == 8
    return w, h, spp, tags[282][0] / tags[282][1], np.frombuffer(data, np.uint8)


def test_sample_sheet_forgiving():
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "sheet.csv")
        open(f, "w", encoding="cp1252").write(
            "EXPNO;Condition;Time (h);Rep;Viability (%)\n10;Fed batch;0 h;1;98\n11.0;Fed batch;2 d;2;97\n"
            "12;Fed batch;day 3;1;95\n99;Ghost;5;1;90\n;;;;\n10;dup;5;1;1\n13;Batch;soon;1;80\n")
        df, warns = F.read_sample_sheet(f, ["10", "11", "12", "13"])
        assert list(df["sample"]) == ["10", "11", "12", "99", "13"]
        assert list(df["time_h"].iloc[:3]) == [0.0, 48.0, 72.0] and np.isnan(df["time_h"].iloc[4])
        assert list(df["group"].iloc[:3]) == ["Fed batch"] * 3
        assert "Viability (%)" in df.columns
        text = " ".join(warns)
        assert "listed twice" in text and "not in these results" in text and "not understood" in text


def test_parsers():
    assert F.parse_ranges("0.5-9.5, 4.6 to 5.0; 1.2\u20131.4") == [(0.5, 9.5), (4.6, 5.0), (1.2, 1.4)]
    assert F.parse_time_h("12 h") == 12 and F.parse_time_h("2 d") == 48 and F.parse_time_h("30 min") == 0.5
    assert F.parse_time_h("D3") == 72 and F.parse_time_h("48:00") == 48 and F.parse_time_h("abc") is None
    segs = F.segments([(0.5, 9.5)], [(4.6, 5.0)], (5.5, 8))
    assert segs == [(5.5, 9.5, 8), (5.0, 5.5, 1.0), (0.5, 4.6, 1.0)]


def test_templates_and_export():
    with tempfile.TemporaryDirectory() as d:
        fake_results(d)
        R = F.Results(d)
        sheet = F.blank_sheet(R.samples, R.titles)
        sheet["group"] = "Culture"
        sheet["time_h"] = [0.0, 0.0, 24.0]
        sheet["replicate"] = ["1", "2", "1"]
        groups = F.trace_groups(R, sheet)
        assert [g["samples"] for g in groups.values()] == [["10", "11"], ["12"]]
        tr = F.build_trace(R, list(groups.values())[0], "mean")
        assert tr["mean"] and np.allclose(np.nanmax(tr["observed"]), np.nanmax(R.spectrum("10").observed) * 1.5,
                                          rtol=0.01)
        labels = F.auto_labels(R, R.samples)
        assert {r["metabolite"] for r in labels} == {"L-Lactate", "Formate"}
        assert all(abs(r["ppm"] - (4.105 if r["metabolite"] == "L-Lactate" else 8.45)) < 1e-6 for r in labels)
        for t in F.TEMPLATES:
            spec = F.default_spec(R)
            spec.update(template=t, traces=list(groups), labels=labels, fit_metabolites=["L-Lactate"],
                        fit_region="1.25-1.40")
            fig, info = F.render(spec, R, F.collect_traces(R, spec, sheet))
            assert info["caption"] and "TSP" in info["caption"] and "\u2014" not in info["caption"]
            f, lw = F.preflight(fig)
            assert f >= F.MIN_FONT_PT and lw >= F.MIN_LINE_PT
        # colours chosen per spectrum, and offsets upwards and sideways
        import matplotlib.colors as mc
        keys = list(groups)
        spec = F.default_spec(R)
        spec.update(template="overlay", traces=keys, colours={keys[0]: "#d62728", keys[1]: "#0072b2"},
                    overlay_offset=10, offset_x_ppm=0.05, inset_region="")
        trs = F.collect_traces(R, spec, sheet)
        fig, info = F.render(spec, R, trs)
        lines = [l for ax in fig.axes for l in ax.get_lines() if len(l.get_xdata()) > 100]
        cols = {mc.to_hex(l.get_color()) for l in lines}
        assert cols == {"#d62728", "#0072b2"}
        red = [l for l in lines if mc.to_hex(l.get_color()) == "#d62728"]
        blue = [l for l in lines if mc.to_hex(l.get_color()) == "#0072b2"]
        assert np.allclose(np.nanmin(red[0].get_xdata()) - np.nanmin(blue[0].get_xdata()), 0.05, atol=1e-6)
        assert "(red)" in info["caption"] and "(blue)" in info["caption"] and "0.05 ppm sideways" in info["caption"]
        spec = F.default_spec(R)
        spec.update(template="representative", traces=list(groups)[:1], labels=labels, width_mm=85, height_mm=50)
        files = F.export(spec, d, os.path.join(d, "figs", "rep"), ["pdf", "svg", "tiff"],
                         ["caption", "data", "record", "script"], sheet=sheet)
        assert all(os.path.exists(x) for x in files) and len(files) == 8
        assert "1, formate; 2, lactate" in open(os.path.join(d, "figs", "rep_caption.txt")).read()
        w, h, spp, dpi, px = read_tiff(os.path.join(d, "figs", "rep.tif"))
        assert dpi == 1200 and abs(w - 85 / 25.4 * 1200) < 2 and abs(h - 50 / 25.4 * 1200) < 2
        assert px.size == w * h * spp and px.min() == 0 and px.max() == 255


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok ", name)
