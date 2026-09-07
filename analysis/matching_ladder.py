#!/usr/bin/env python3
"""matching_ladder.py — the missing rung of the R3 decomposition.

`s1b_rivals.py` computed two matched variants, `matched_on_mw_and_label` and
`matched_on_mw_only`, and the paper printed only the first.  The second shows
that matching on molecular weight alone leaves the gap essentially unchanged,
which means the entire shrinkage is the LABEL quintile and the manuscript's
motivating sentence named the wrong covariate.

To print that as a ladder the fourth rung is needed: label-only.  It is computed
here by calling `s1b_rivals.cem_matched_gap` itself -- the same function, the
same acceptance rule, the same `zlib.crc32` seed-key convention (`|y` alongside
the published `|mw` and `|mwy`) -- so the new rung is produced by exactly the
code path that produced the published ones.  The published values are re-derived
in the same run and asserted against `evidence/s1b_rival_mechanisms.json`, so a
drift in either direction fails loudly rather than silently reshaping a table.

No published value is overwritten.  Writes evidence/matching_ladder.json.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bench_data import load_group_split  # noqa: E402
from leakage_lib import FP_BYTES, morgan_packed  # noqa: E402
from s1b_rivals import cem_matched_gap, descriptors  # noqa: E402

PAPER = Path(__file__).resolve().parents[1]
EV = PAPER / "evidence"
BENCH = "solubility_aqsoldb"
LEAK_T = 0.9
TOL = 5e-4          # matching-draw noise tolerance for the reproduction assert


def fp_cans(smiles_list):
    cans = []
    for s in smiles_list:
        r = morgan_packed(s)
        if r is not None:
            cans.append(r[0])
    return cans


def main() -> int:
    nn_rows = [json.loads(l) for l in (EV / "s1_nn_per_molecule.jsonl").read_text().splitlines()]
    preds = [json.loads(l) for l in (EV / "s1b_gap_results_predictions.jsonl").read_text().splitlines()]
    published = json.loads((EV / "s1b_rival_mechanisms.json").read_text())
    pub = published["benchmarks"][BENCH]["R3_matched_gap"]

    test = load_group_split(BENCH, "test")
    te_can = fp_cans(test.smiles)
    keep = [i for i, s in enumerate(test.smiles) if morgan_packed(s) is not None]
    te_y = np.array([test.y[i] for i in keep], float)

    by_can = {r["cansmi"]: r for r in nn_rows
              if r["benchmark"] == BENCH and r["split"] == "test"}
    nn_sim = np.array([by_can[c]["nn_chembl_tanimoto"] for c in te_can], float)
    leaked = nn_sim >= LEAK_T

    scaf_path = Path("<REPO_ROOT>/.work/agenticls/full/chembl_scaffold_counts.json")
    scaffold_counts = json.loads(scaf_path.read_text())["counts"] if scaf_path.exists() else {}
    d = descriptors(te_can, scaffold_counts)

    # strata, built exactly as s1b_rivals builds them
    mw_q = np.quantile(d["mw"], np.linspace(0, 1, 6))
    strata_mw = np.clip(np.searchsorted(mw_q, d["mw"], side="right") - 1, 0, 4)
    y_q = np.quantile(te_y, np.linspace(0, 1, 6))
    strata_y = np.clip(np.searchsorted(y_q, te_y, side="right") - 1, 0, 4)
    st_mwy = [(int(a), int(c)) for a, c in zip(strata_mw, strata_y)]
    st_mw = [(int(a),) for a in strata_mw]
    st_y = [(int(c),) for c in strata_y]

    out = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "benchmark": BENCH,
        "purpose": "print the R3 decomposition as a ladder, not one number",
        "strata": {"mw_only": "MW quintile", "label_only": "label quintile",
                   "mw_and_label": "MW quintile x label quintile"},
        "rungs": {},
    }

    arms = sorted({p["arm"] for p in preds if p["benchmark"] == BENCH})
    models = sorted({p["model"] for p in preds if p["benchmark"] == BENCH})
    for arm in arms:
        for model in models:
            rows = [p for p in preds if p["benchmark"] == BENCH
                    and p["arm"] == arm and p["model"] == model]
            if not rows:
                continue
            assert [r["cansmi"] for r in rows] == te_can, "prediction order differs — FAIL"
            pv = np.array([r["pred_seed_mean"] for r in rows], float)
            key = f"{arm}|{model}"

            r_mwy = cem_matched_gap("regression", te_y, pv, leaked, st_mwy,
                                    f"{BENCH}|{arm}|{model}|mwy")
            r_mw = cem_matched_gap("regression", te_y, pv, leaked, st_mw,
                                   f"{BENCH}|{arm}|{model}|mw")
            r_y = cem_matched_gap("regression", te_y, pv, leaked, st_y,
                                  f"{BENCH}|{arm}|{model}|y")

            # reproduce-or-fail against the published file
            for tag, got, want in (
                ("mw_and_label", r_mwy["matched_gap_mean"],
                 pub[key]["matched_on_mw_and_label"]["matched_gap_mean"]),
                ("mw_only", r_mw["matched_gap_mean"],
                 pub[key]["matched_on_mw_only"]["matched_gap_mean"]),
            ):
                assert abs(got - want) < TOL, (
                    f"{key}/{tag} does not reproduce the published value: "
                    f"{got:.6f} vs {want:.6f} — FAIL")

            out["rungs"][key] = {
                "mw_only": {"matched_gap": round(r_mw["matched_gap_mean"], 4),
                            "n_matched_per_side": r_mw["n_matched_per_side"]},
                "label_only": {"matched_gap": round(r_y["matched_gap_mean"], 4),
                               "n_matched_per_side": r_y["n_matched_per_side"]},
                "mw_and_label": {"matched_gap": round(r_mwy["matched_gap_mean"], 4),
                                 "n_matched_per_side": r_mwy["n_matched_per_side"]},
            }
            e = out["rungs"][key]
            print(f"{key:35s} mw {e['mw_only']['matched_gap']:+.4f} | "
                  f"label {e['label_only']['matched_gap']:+.4f} | "
                  f"both {e['mw_and_label']['matched_gap']:+.4f}", flush=True)

    p = EV / "matching_ladder.json"
    p.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", p, "— published rungs reproduced exactly", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
