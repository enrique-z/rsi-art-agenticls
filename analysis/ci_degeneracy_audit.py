#!/usr/bin/env python3
"""ci_degeneracy_audit.py — audit EVERY interval this paper reports.

An n-unit percentile bootstrap of a mean is *degenerate* when its 2.5/97.5
quantiles coincide with the [min, max] of the n raw unit values, or when the
interval has zero width. At that point "the interval excludes zero" is a
restatement of "all n values share a sign" -- a sign test at n units, not an
interval. The coincidence is forced for n <= 3 and can happen by accident above.

This walks every interval that reaches the PDF, resolves the resampling unit and
the unit count behind it from the evidence files, checks degeneracy, and adds an
independent Student-t interval on the same units wherever the unit count is
small. It reads only; it writes evidence/ci_degeneracy_audit.json.

Interval families in this paper:
  W  Wilson score interval on a proportion            unit = test molecule
  M  molecule-level percentile bootstrap of a gap     unit = test molecule
  S  seed-level percentile bootstrap of a gap         unit = fitting seed
       (in the ledger, never quoted in the paper -- amendment A6)
  P  paired bootstrap over near-duplicate pairs       unit = pair
  N  nested bootstrap, 1000 outer molecule draws      unit = test molecule
  U  unit-level percentile bootstrap over the loop    unit = curator x seed
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np

PAPER = Path(__file__).resolve().parents[1]
EV = PAPER / "evidence"
LEDGER = EV / "agenticls-rsi-leakage_ledger.jsonl"
SMALL_N = 30
TOL = 5e-5


def t_crit(df: int) -> float:
    try:
        from scipy import stats
        return float(stats.t.ppf(0.975, df))
    except Exception:  # pragma: no cover
        return {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447,
                7: 2.365, 8: 2.306, 9: 2.262, 17: 2.110}.get(df, 1.96)


def t_interval(vals):
    v = np.asarray(vals, float)
    if len(v) < 2:
        return None
    m, se = float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v)))
    if se == 0.0:
        return [round(m, 4), round(m, 4)]
    h = t_crit(len(v) - 1) * se
    return [round(m - h, 4), round(m + h, 4)]


def loop_units():
    """Every interval whose unit is a curator x seed (or a curator call)."""
    out = {}
    s2 = json.loads((EV / "s2_rsi_results.json").read_text())
    llm, ctrl = [], []
    for name, arm in s2["arms"].items():
        vals = [s["gain_full_test"] - s["gain_nonleaked_test"]
                for s in arm["per_seed"]]
        (ctrl if "random" in name else llm).extend(vals)
    pooled = s2["pooled_llm_overstatement"]
    assert len(llm) == pooled["n_units"]
    assert abs(float(np.mean(llm)) - pooled["mean"]) < 1e-12
    out["s2_pooled::overstatement"] = ("LLM curator x seed", llm)
    out["s2_random_curator_control::full_minus_clean"] = ("seed", ctrl)

    c7 = json.loads((EV / "s2_shuffled_feedback.json").read_text())
    for k, row in c7["paired_real_minus_shuffled"].items():
        out[f"c7_paired::{k}"] = (row["unit"], row["per_unit"])
    out["c7_selection::jaccard"] = ("curator x seed x round call",
                                    c7["selection_agreement"]["per_call"])
    return out


def pair_counts():
    """r2_* claim prefix -> number of near-duplicate pairs behind the interval."""
    riv = json.loads((EV / "s1b_rival_mechanisms.json").read_text())
    out = {}
    for bench, b in riv["benchmarks"].items():
        for ref, r in b.get("R2_label_agreement", {}).items():
            if not r.get("n_near_pairs"):
                continue
            out[f"r2_{bench}_{ref}"] = r["n_near_pairs"]
    return out


def gap_cells():
    """s1b_* claim prefix -> (n_leaked, n_clean, molecule CI, seed CI)."""
    g = json.loads((EV / "s1b_gap_results.json").read_text())
    out = {}
    for bench, b in g["benchmarks"].items():
        for arm, a in b["arms"].items():
            for fam, m in a["models"].items():
                out[f"s1b_{bench}_{arm}_{fam}"] = {
                    "n_leaked": b["n_leaked"], "n_clean": b["n_nonleaked"],
                    "n_seeds": m["n_seeds"],
                    "mol_ci": [round(x, 4) for x in m["gap_ci_molecule_bootstrap"]],
                    "seed_ci": [round(x, 4) for x in m["gap_ci_seed_bootstrap"]],
                    "per_seed_gap": [s["gap"] for s in m["per_seed"]],
                }
    return out


def main() -> int:
    rows = [json.loads(l) for l in LEDGER.read_text().splitlines() if l.strip()]
    lu, pc, gc = loop_units(), pair_counts(), gap_cells()

    report, degenerate, small = [], [], []

    def add(cid, value, ci, kind, unit, n, vals=None, note=""):
        rec = {"claim_id": cid, "value": value, "ci": ci, "kind": kind,
               "unit": unit, "n_units": n, "degenerate": False,
               "t_ci": None, "note": note}
        if ci is not None and ci[0] is not None:
            if abs(ci[1] - ci[0]) < TOL:
                rec["degenerate"] = True
                rec["note"] = (note + " " if note else "") + \
                    "ZERO-WIDTH interval: every resample gives the same value."
            if vals is not None and len(vals):
                v = np.asarray(vals, float)
                rec["unit_min_max"] = [round(float(v.min()), 4),
                                       round(float(v.max()), 4)]
                if abs(ci[0] - v.min()) < TOL and abs(ci[1] - v.max()) < TOL:
                    rec["degenerate"] = True
                    rec["note"] = (rec["note"] + " " if rec["note"] else "") + \
                        "CI == [min,max] of the units: a sign test, not an interval."
        if vals is not None and 2 <= len(vals) <= SMALL_N:
            rec["t_ci"] = t_interval(vals)
        if rec["degenerate"]:
            degenerate.append(cid)
        if n is not None and n <= SMALL_N:
            small.append(cid)
        report.append(rec)
        return rec

    # ---- W: Wilson intervals on proportions -----------------------------
    wilson_marks = ("leak_rate", "dup_rate", "own_trainval_rate",
                    "union_rate", "drugbank_rate", "endpoint_pool_rate")
    for r in rows:
        cid = r["claim_id"]
        if r["ci_low"] is None:
            continue
        if any(w in cid for w in wilson_marks):
            add(cid, r["value"], [r["ci_low"], r["ci_high"]],
                "Wilson score interval (closed form)", "test molecule", None,
                note="closed form on n>=406 molecules; not a resample")

    # ---- M / S: the leaked-minus-clean gap ------------------------------
    for prefix, c in sorted(gc.items()):
        add(f"{prefix}::gap [molecule-level, QUOTED]",
            None, c["mol_ci"], "molecule-level percentile bootstrap",
            "test molecule", min(c["n_leaked"], c["n_clean"]),
            note=(f"leaked {c['n_leaked']} / clean {c['n_clean']}; the smaller "
                  "side bounds the resolution of the interval"))
        add(f"{prefix}::gap [seed-level, NOT QUOTED]",
            None, c["seed_ci"], "seed-level percentile bootstrap",
            "fitting seed", c["n_seeds"], vals=c["per_seed_gap"],
            note="ledger-only; amendment A6 forbids quoting it in the paper")

    # ---- P: paired bootstrap over near-duplicate pairs ------------------
    for r in rows:
        cid = r["claim_id"]
        if r["ci_low"] is None or not cid.startswith("r2_"):
            continue
        n = next((v for k, v in pc.items() if cid.startswith(k)), None)
        add(cid, r["value"], [r["ci_low"], r["ci_high"]],
            "paired bootstrap over near-duplicate pairs",
            "near-duplicate pair", n)

    # ---- N: nested bootstrap --------------------------------------------
    for r in rows:
        cid = r["claim_id"]
        if r["ci_low"] is None:
            continue
        if "matched_gap_nested" in cid or "paired_difference" in cid:
            add(cid, r["value"], [r["ci_low"], r["ci_high"]],
                "nested bootstrap, 1000 outer draws over test molecules",
                "test molecule (re-stratified and re-matched inside each draw)",
                1000)
        elif "matched_gap" in cid:
            add(cid, r["value"], [r["ci_low"], r["ci_high"]],
                "bootstrap over matched molecules", "matched test molecule",
                None)

    # ---- U: the loop's curator x seed units -----------------------------
    for cid, (unit, vals) in lu.items():
        led = next((r for r in rows if r["claim_id"] == cid), None)
        ci = [led["ci_low"], led["ci_high"]] if led and led["ci_low"] is not None else None
        add(cid, led["value"] if led else round(float(np.mean(vals)), 4), ci,
            "unit-level percentile bootstrap", unit, len(vals), vals=vals,
            note="" if ci else "no interval is printed for this quantity")

    out = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "small_n_threshold": SMALL_N,
        "n_intervals_audited": len(report),
        "degenerate_intervals": sorted(set(degenerate)),
        "n_degenerate": len(set(degenerate)),
        "rows": report,
    }
    (EV / "ci_degeneracy_audit.json").write_text(json.dumps(out, indent=1) + "\n")

    print(f"audited {len(report)} intervals across {len(set(r['kind'] for r in report))} families")
    print(f"DEGENERATE: {out['n_degenerate']}")
    for cid in out["degenerate_intervals"]:
        rec = next(r for r in report if r["claim_id"] == cid)
        print(f"  {cid:<52} n={rec['n_units']} ci={rec['ci']} :: {rec['note']}")
    print("\nsmall-unit intervals with an independent t-interval:")
    for rec in report:
        if rec["t_ci"] is not None:
            print(f"  {rec['claim_id']:<52} n={rec['n_units']:>3} "
                  f"boot {rec['ci']} t {rec['t_ci']} units{rec.get('unit_min_max')}")
    print("\nsmall-unit intervals with no per-unit vector on disk "
          "(resample is over molecules, not a listed unit):")
    for rec in report:
        if rec["t_ci"] is None and rec["n_units"] is not None \
                and rec["n_units"] <= SMALL_N:
            print(f"  {rec['claim_id']:<52} n={rec['n_units']:>3} ci={rec['ci']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
