"""Checks that a results folder is never silently re-used (MANC-Q 1.8.1), on the synthetic demo spectrum only.
Fits EXPNO 10 of the demo twice (a few minutes).
Run with:  python -m pytest tests   (or:  python tests/test_overwrite.py)"""
from pathlib import Path
import os
import sys
import json
import tempfile
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mancq import engine as q
from mancq import demo as dm
from mancq.__main__ import main


def _read(p):
    return Path(p).read_text(encoding="utf-8")


def _write(p, text):
    Path(p).write_text(text, encoding="utf-8")


def lactate(res):
    g = pd.read_csv(os.path.join(res, "concentrations_long.csv"), dtype={"sample": str})
    return float(g[(g["sample"] == "10") & (g.metabolite == "L-Lactate")].conc_mM.iloc[0])


def test_folder_helpers():
    with tempfile.TemporaryDirectory() as out:
        assert not q.has_results(out)
        os.makedirs(os.path.join(out, "_per_sample"))
        os.makedirs(os.path.join(out, "_cache"))
        os.makedirs(os.path.join(out, "_processed", "11"))
        for f in ("concentrations_long.csv", "overlay_sample_10.pdf", "simulated_spectrum_sample_10.csv.gz",
                  "my_notes.txt", "sample_sheet.csv", os.path.join("_per_sample", "10_manual.json")):
            _write(os.path.join(out, f), "")
        assert q.has_results(out)
        sub = q.new_run_folder(out)
        assert os.path.dirname(sub) == out and os.path.basename(sub).startswith("run_") and not os.path.exists(sub)
        q.clear_results(out)
        assert not q.has_results(out)
        assert sorted(os.listdir(out)) == ["_cache", "my_notes.txt", "sample_sheet.csv"]


def test_rerun_with_other_tsp():
    with tempfile.TemporaryDirectory() as tmp:
        spec, res = os.path.join(tmp, "spectra"), os.path.join(tmp, "results")
        dm.write_demo(spec, cache_dir=os.path.join(tmp, "cache"))
        args = ["run", spec, res, "--samples", "10", "--no-overlays", "--workers", "1"]
        assert main(args + ["--ref-mm", "0.5"]) in (None, 0)
        first = lactate(res)
        assert abs(first - 3.0) < 0.15, first
        # the same folder again: refused, nothing changed
        assert main(args + ["--ref-mm", "1.0"]) == 2
        assert lactate(res) == first
        # carrying on a run that finished with a different reference: refitted, not re-used
        prev = q.previous_run(res)
        assert prev["complete"] and prev["samples"] == ["10"]
        rs = os.path.join(res, "_per_sample", "run_settings.json")
        d = json.loads(_read(rs)); d["RUN_COMPLETE"] = False; _write(rs, json.dumps(d))   # as if it was stopped
        assert main(args + ["--ref-mm", "1.0", "--resume"]) in (None, 0)
        log = _read(os.path.join(res, "run_log.txt"))
        assert "cannot be used" in log and "TSP_MM_IN_TUBE" in log and "re-used" not in log, log
        doubled = lactate(res)
        assert abs(doubled / first - 2.0) < 0.01, (first, doubled)
        # carrying on with the same settings re-uses it
        d = json.loads(_read(rs)); d["RUN_COMPLETE"] = False; _write(rs, json.dumps(d))
        assert main(args + ["--ref-mm", "1.0", "--resume"]) in (None, 0)
        assert "re-used earlier result" in _read(os.path.join(res, "run_log.txt"))
        assert lactate(res) == doubled
        # overwrite: manual adjustments go, value from the new reference
        _write(os.path.join(res, "_per_sample", "10_manual.json"), '{"edits": {}}')
        assert main(args + ["--ref-mm", "0.5", "--overwrite"]) in (None, 0)
        assert not os.path.exists(os.path.join(res, "_per_sample", "10_manual.json"))
        assert abs(lactate(res) - first) < 1e-6


if __name__ == "__main__":
    test_folder_helpers()
    test_rerun_with_other_tsp()
    print("ok")
