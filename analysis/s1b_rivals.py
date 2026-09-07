#!/usr/bin/env python3
"""s1b_rivals.py — the three pre-declared rival mechanisms for the leakage claim.

R1  chemistry-overlap-is-normal
    Handled inside s1b_gap.py by arm B (train on the benchmark's own train_val
    split).  This script re-reports the arm-A-vs-arm-B contrast from the
    predictions file so it sits in one place.

R2  endpoint-value mismatch (activity cliffs)
    Do structural near-duplicates actually agree on the LABEL?  Two pairings:
      (a) test molecule  <->  its nearest neighbour in the benchmark's own
          train_val split, near-duplicate pairs (Tanimoto >= 0.9) vs random pairs;
      (b) test molecule  <->  its nearest neighbour in the ChEMBL endpoint pool
          (lipophilicity / solubility only -- ChEMBL carries no BBB labels).
    If near-duplicate pairs agree no better than random pairs, the duplicates
    carry no exploitable signal and the leakage mechanism is FALSIFIED.

R3  subset-composition confound
    Covariate distributions of the leaked and non-leaked test subsets (molecular
    weight, heavy-atom count, TPSA, Bemis-Murcko scaffold frequency in ChEMBL,
    label balance / label mean), plus a coarsened-exact-matched re-estimate of
    the leaked/non-leaked performance gap.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import zlib
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bench_data import BENCHMARKS, load_group_split  # noqa: E402
from leakage_lib import (  # noqa: E402
    FP_BYTES,
    bemis_murcko,
    bootstrap_ci,
    exact_nn_tanimoto,
    morgan_packed,
)

PAPER = Path(__file__).resolve().parents[1]
LEAK_T = 0.9
N_BOOT = 10000

ENDPOINT_FILE = {
    "lipophilicity_astrazeneca": "chembl_endpoint_lipophilicity.jsonl",
    "solubility_aqsoldb": "chembl_endpoint_solubility.jsonl",
    "bbb_martins": None,
}


def fp_matrix(smiles_list):
    packed, keep, cans = [], [], []
    for i, s in enumerate(smiles_list):
        r = morgan_packed(s)
        if r is None:
            continue
        cans.append(r[0])
        packed.append(r[1])
        keep.append(i)
    arr = np.frombuffer(b"".join(packed), dtype=np.uint8).reshape(len(packed), FP_BYTES)
    return arr, keep, cans


def pair_agreement(task, y_q, y_ref, sim, argmax, key) -> dict:
    """Near-duplicate vs random pair agreement on the LABEL.

    regression      -> mean |dy| (LOWER is more agreement)
    classification  -> fraction of pairs with equal label (HIGHER is more agreement)
    """
    near = np.flatnonzero(sim >= LEAK_T)
    if len(near) == 0:
        return {"n_near_pairs": 0, "note": "no near-duplicate pairs at t=0.9"}
    rng = np.random.default_rng(zlib.crc32(key.encode()))
    rand_ref = rng.integers(0, len(y_ref), size=len(near))
    if task == "regression":
        obs = np.abs(y_q[near] - y_ref[argmax[near]])
        rnd = np.abs(y_q[near] - y_ref[rand_ref])
        stat_name = "mean_abs_label_difference"
        better = "lower"
    else:
        obs = (y_q[near].astype(int) == y_ref[argmax[near]].astype(int)).astype(float)
        rnd = (y_q[near].astype(int) == y_ref[rand_ref].astype(int)).astype(float)
        stat_name = "label_agreement_rate"
        better = "higher"
    diff = obs - rnd
    lo, hi = bootstrap_ci(diff, n_boot=N_BOOT, seed_key=key + "|diff")
    olo, ohi = bootstrap_ci(obs, n_boot=N_BOOT, seed_key=key + "|obs")
    rlo, rhi = bootstrap_ci(rnd, n_boot=N_BOOT, seed_key=key + "|rnd")
    return {
        "statistic": stat_name,
        "agreement_direction": better,
        "n_near_pairs": int(len(near)),
        "near_duplicate_pairs": {"mean": float(obs.mean()), "ci": [olo, ohi]},
        "random_pairs": {"mean": float(rnd.mean()), "ci": [rlo, rhi]},
        "paired_difference_near_minus_random": {"mean": float(diff.mean()), "ci": [lo, hi]},
        "difference_ci_excludes_zero": bool(lo > 0 or hi < 0),
    }


def descriptors(cansmi_list, scaffold_counts) -> dict:
    mw, ha, tpsa, scf = [], [], [], []
    for s in cansmi_list:
        m = Chem.MolFromSmiles(s)
        mw.append(Descriptors.MolWt(m))
        ha.append(m.GetNumHeavyAtoms())
        tpsa.append(Descriptors.TPSA(m))
        sc = bemis_murcko(s)
        scf.append(scaffold_counts.get(sc or "", 0))
    return {"mw": np.array(mw), "heavy_atoms": np.array(ha, float),
            "tpsa": np.array(tpsa), "chembl_scaffold_count": np.array(scf, float)}


def smd(a: np.ndarray, b: np.ndarray) -> float:
    """Standardised mean difference (leaked minus non-leaked)."""
    s = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2) if len(a) > 1 and len(b) > 1 else np.nan
    return float((a.mean() - b.mean()) / s) if s and s > 0 else float("nan")


def cem_matched_gap(task, y, pred, leaked, strata_keys, key, n_rep=200) -> dict:
    """Coarsened-exact-matched re-estimate of the leaked/non-leaked gap.

    Molecules are binned by the covariate strata; within each stratum the larger
    side is down-sampled to the size of the smaller.  Repeated `n_rep` times with
    different draws; the gap distribution over repeats gives the CI.
    """
    rng = np.random.default_rng(zlib.crc32(key.encode()))
    by_stratum: dict[tuple, dict[bool, list[int]]] = defaultdict(lambda: {True: [], False: []})
    for i, s in enumerate(strata_keys):
        by_stratum[s][bool(leaked[i])].append(i)
    usable = {s: v for s, v in by_stratum.items() if v[True] and v[False]}
    n_matched_pairs = sum(min(len(v[True]), len(v[False])) for v in usable.values())
    if n_matched_pairs < 10:
        return {"n_matched_per_side": int(n_matched_pairs),
                "note": "too few matched molecules for a stable estimate"}
    gaps = []
    for _ in range(n_rep):
        li, ni = [], []
        for v in usable.values():
            k = min(len(v[True]), len(v[False]))
            li.extend(rng.choice(v[True], k, replace=False).tolist())
            ni.extend(rng.choice(v[False], k, replace=False).tolist())
        li = np.array(li)
        ni = np.array(ni)
        if task == "regression":
            g = np.abs(pred[ni] - y[ni]).mean() - np.abs(pred[li] - y[li]).mean()
        else:
            if len(set(y[li].tolist())) < 2 or len(set(y[ni].tolist())) < 2:
                continue
            g = roc_auc_score(y[li], pred[li]) - roc_auc_score(y[ni], pred[ni])
        gaps.append(g)
    if not gaps:
        return {"n_matched_per_side": int(n_matched_pairs),
                "note": "matched strata never contained both classes"}
    gaps = np.array(gaps)
    return {
        "n_strata_used": len(usable),
        "n_matched_per_side": int(n_matched_pairs),
        "matched_gap_mean": float(gaps.mean()),
        "matched_gap_ci": [float(np.quantile(gaps, 0.025)), float(np.quantile(gaps, 0.975))],
        "matched_ci_excludes_zero": bool(np.quantile(gaps, 0.025) > 0
                                         or np.quantile(gaps, 0.975) < 0),
        "n_repeats": len(gaps),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="<REPO_ROOT>/.work/agenticls/full")
    ap.add_argument("--nn", default=str(PAPER / "evidence" / "s1_nn_per_molecule.jsonl"))
    ap.add_argument("--preds", default=str(PAPER / "evidence" / "s1b_gap_results_predictions.jsonl"))
    ap.add_argument("--out", default=str(PAPER / "evidence" / "s1b_rival_mechanisms.json"))
    args = ap.parse_args()

    work = Path(args.work)
    nn_rows = [json.loads(l) for l in Path(args.nn).read_text().splitlines()]
    assert nn_rows, "empty NN ledger — FAIL"
    scaf_path = work / "chembl_scaffold_counts.json"
    scaffold_counts = json.loads(scaf_path.read_text())["counts"] if scaf_path.exists() else {}
    preds = [json.loads(l) for l in Path(args.preds).read_text().splitlines()]
    assert preds, "empty predictions file — FAIL"

    out = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "leak_threshold": LEAK_T,
        "chembl_scaffold_table": {"path": str(scaf_path), "available": bool(scaffold_counts)},
        "benchmarks": {},
    }

    for b, cfg in BENCHMARKS.items():
        task = cfg["task"]
        test = load_group_split(b, "test")
        tv = load_group_split(b, "train_val")
        te_fp, te_keep, te_can = fp_matrix(test.smiles)
        te_y = np.array([test.y[i] for i in te_keep], float)
        tv_fp, tv_keep, tv_can = fp_matrix(tv.smiles)
        tv_y = np.array([tv.y[i] for i in tv_keep], float)

        by_can = {r["cansmi"]: r for r in nn_rows
                  if r["benchmark"] == b and r["split"] == "test"}
        nn_sim = np.array([by_can[c]["nn_chembl_tanimoto"] for c in te_can], float)
        leaked = nn_sim >= LEAK_T

        entry: dict = {
            "display": cfg["display"], "task": task,
            "n_test": len(te_can),
            "n_leaked": int(leaked.sum()), "n_nonleaked": int((~leaked).sum()),
            "R2_label_agreement": {},
            "R3_covariates": {},
            "R3_matched_gap": {},
        }

        # ---- R2 (a) test <-> own train_val -------------------------------
        r = exact_nn_tanimoto(te_fp, tv_fp, chunk=8192)
        entry["R2_label_agreement"]["vs_own_trainval"] = pair_agreement(
            task, te_y, tv_y, r.max_sim, r.argmax, f"{b}|tv")

        # ---- R2 (b) test <-> ChEMBL endpoint pool ------------------------
        ep = ENDPOINT_FILE[b]
        if ep:
            pool = [json.loads(l) for l in (work / ep).read_text().splitlines()]
            p_fp, p_keep, p_can = fp_matrix([x["smiles"] for x in pool])
            p_y = np.array([pool[i]["y"] for i in p_keep], float)
            rp = exact_nn_tanimoto(te_fp, p_fp, chunk=8192)
            entry["R2_label_agreement"]["vs_chembl_endpoint_pool"] = pair_agreement(
                task, te_y, p_y, rp.max_sim, rp.argmax, f"{b}|pool")
            entry["R2_label_agreement"]["vs_chembl_endpoint_pool"]["pool_size"] = len(p_y)
        else:
            entry["R2_label_agreement"]["vs_chembl_endpoint_pool"] = {
                "available": False,
                "reason": "ChEMBL 36 carries no usable BBB label pool (see chembl_endpoints.json)",
            }

        # ---- R3 covariates -------------------------------------------------
        d = descriptors(te_can, scaffold_counts)
        cov = {}
        for name, v in d.items():
            lo_l, hi_l = bootstrap_ci(v[leaked], seed_key=f"{b}|{name}|L")
            lo_n, hi_n = bootstrap_ci(v[~leaked], seed_key=f"{b}|{name}|N")
            cov[name] = {
                "leaked": {"mean": float(v[leaked].mean()),
                           "median": float(np.median(v[leaked])), "ci": [lo_l, hi_l]},
                "nonleaked": {"mean": float(v[~leaked].mean()),
                              "median": float(np.median(v[~leaked])), "ci": [lo_n, hi_n]},
                "standardised_mean_difference": smd(v[leaked], v[~leaked]),
            }
        if task == "classification":
            cov["label_balance"] = {
                "leaked_positive_rate": float(te_y[leaked].mean()),
                "nonleaked_positive_rate": float(te_y[~leaked].mean()),
            }
        else:
            cov["label"] = {
                "leaked_mean": float(te_y[leaked].mean()),
                "nonleaked_mean": float(te_y[~leaked].mean()),
                "leaked_sd": float(te_y[leaked].std(ddof=1)),
                "nonleaked_sd": float(te_y[~leaked].std(ddof=1)),
                "standardised_mean_difference": smd(te_y[leaked], te_y[~leaked]),
            }
        scaf_counter = Counter(bemis_murcko(c) for c in te_can)
        cov["within_test_scaffold_frequency"] = {
            "leaked_mean": float(np.mean([scaf_counter[bemis_murcko(c)]
                                          for c, L in zip(te_can, leaked) if L])),
            "nonleaked_mean": float(np.mean([scaf_counter[bemis_murcko(c)]
                                             for c, L in zip(te_can, leaked) if not L])),
        }
        entry["R3_covariates"] = cov

        # ---- R3 matched gap, per arm x model --------------------------------
        mw_q = np.quantile(d["mw"], np.linspace(0, 1, 6))
        strata_mw = np.clip(np.searchsorted(mw_q, d["mw"], side="right") - 1, 0, 4)
        if task == "classification":
            strata = [(int(a), int(y)) for a, y in zip(strata_mw, te_y)]
            strata_desc = "MW quintile x label class"
        else:
            y_q = np.quantile(te_y, np.linspace(0, 1, 6))
            strata_y = np.clip(np.searchsorted(y_q, te_y, side="right") - 1, 0, 4)
            strata = [(int(a), int(c)) for a, c in zip(strata_mw, strata_y)]
            strata_desc = "MW quintile x label quintile"
        strata_mw_only = [(int(a),) for a in strata_mw]

        # Match predictions POSITIONALLY, not by canonical SMILES: the
        # BBB-Martins test split contains 12 duplicate canonical structures, so a
        # SMILES-keyed join silently drops rows.  s1b_gap.py writes one row per
        # test molecule in the same fp_matrix order this script reproduces, and
        # the SMILES sequence is asserted to match.
        for arm in sorted({p["arm"] for p in preds if p["benchmark"] == b}):
            for model in sorted({p["model"] for p in preds
                                 if p["benchmark"] == b and p["arm"] == arm}):
                rows = [p for p in preds if p["benchmark"] == b
                        and p["arm"] == arm and p["model"] == model]
                assert len(rows) == len(te_can), (
                    f"prediction row count {len(rows)} != {len(te_can)} for "
                    f"{b}/{arm}/{model} — FAIL")
                assert [r["cansmi"] for r in rows] == te_can, (
                    f"prediction row ORDER differs for {b}/{arm}/{model} — FAIL")
                pv = np.array([r["pred_seed_mean"] for r in rows], dtype=float)
                assert not np.isnan(pv).any(), (
                    f"NaN prediction for {b}/{arm}/{model} — FAIL")
                entry["R3_matched_gap"][f"{arm}|{model}"] = {
                    "strata_definition": strata_desc,
                    "matched_on_mw_and_label": cem_matched_gap(
                        task, te_y, pv, leaked, strata, f"{b}|{arm}|{model}|mwy"),
                    "matched_on_mw_only": cem_matched_gap(
                        task, te_y, pv, leaked, strata_mw_only, f"{b}|{arm}|{model}|mw"),
                }
        out["benchmarks"][b] = entry
        print(f"{b}: done", flush=True)

    Path(args.out).write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
