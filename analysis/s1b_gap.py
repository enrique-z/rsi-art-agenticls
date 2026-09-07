#!/usr/bin/env python3
"""s1b_gap.py — Stage 1b: does ChEMBL near-duplicate leakage inflate test scores?

Trains ADMET predictors and evaluates them SEPARATELY on the leaked
(NN-to-ChEMBL Tanimoto >= 0.9) and non-leaked (< 0.9) halves of each TDC ADMET
Benchmark Group test split.

Arms
----
A2  PRIMARY, pre-registered.  Train on the ChEMBL-derived pool for the endpoint
    with every molecule whose RDKit canonical SMILES equals a TDC *test* molecule
    removed -- the pre-declaration's "never on TDC test".
A1  Naive curator, supplementary.  Same pool WITHOUT that removal: what an LLM
    curator that does not deduplicate against the benchmark actually trains on.
B   Rival-mechanism control.  Train on the benchmark's OWN train_val split.  If
    the leaked/non-leaked gap also appears here, the gap is at least partly
    "molecules similar to ChEMBL are generically easier", not ChEMBL-derived
    training data leaking the answer.

Model families: RandomForest and HistGradientBoosting, both on the same
Morgan r=2 / 2048-bit fingerprints.  >= 5 seeds (default 10).

Every number written here comes from a real fit on real data.  Nothing is
simulated.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import zlib
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import LogisticRegression, Ridge
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bench_data import BENCHMARKS, load_group_split  # noqa: E402
from leakage_lib import FP_BYTES, exact_nn_tanimoto, morgan_packed, unpack_f32  # noqa: E402

os.environ.setdefault("PYTHONHASHSEED", "0")
PAPER = Path(__file__).resolve().parents[1]
LEAK_T = 0.9

ENDPOINT_FILE = {
    "lipophilicity_astrazeneca": "chembl_endpoint_lipophilicity.jsonl",
    "solubility_aqsoldb": "chembl_endpoint_solubility.jsonl",
    "bbb_martins": None,   # ChEMBL 36 has no usable BBB pool -- see DESIGN_NOTES
}


# --------------------------------------------------------------------------- #
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


MODEL_NAMES = ("random_forest", "linear_ecfp", "knn_tanimoto")


def make_models(task: str, seed: int, small: bool):
    """Three model families on identical Morgan r=2 / 2048-bit fingerprints.

    * `random_forest`  -- the pre-declared primary (max_features='sqrt' is the
      standard ECFP4 QSAR setting; sklearn's regressor default scans all 2048
      columns per split and is ~6x slower for no accuracy gain here).
    * `linear_ecfp`    -- the second model family: L2-regularised linear model on
      the same bits, a standard cheminformatics baseline and genuinely not a
      tree ensemble.
    * `knn_tanimoto`   -- k=5 similarity-weighted nearest neighbours under the
      SAME Tanimoto metric the leakage census uses.  This is the mechanism probe:
      it is the model that can exploit a near-duplicate most directly, so it
      upper-bounds how much a near-duplicate is worth.  Reported as an audit row
      (must_print=False), never as the headline.

    HistGradientBoosting was the original second family and was DROPPED after
    measurement: one 200-iteration fit on the 1,462-molecule BBB training pool
    took 1,259 s on this machine (sklearn 1.8.0, 2048 binary columns), which puts
    10 seeds x 7 arm/model combinations at ~24 h.  Recorded as amendment A5.
    """
    n_est = 60 if small else 200
    if task == "regression":
        return {
            "random_forest": RandomForestRegressor(
                n_estimators=n_est, n_jobs=-1, random_state=seed,
                min_samples_leaf=1, max_features="sqrt"),
            "linear_ecfp": Ridge(alpha=1.0, random_state=seed),
            "knn_tanimoto": TanimotoKNN(k=5, task=task),
        }
    return {
        "random_forest": RandomForestClassifier(
            n_estimators=n_est, n_jobs=-1, random_state=seed,
            min_samples_leaf=1, max_features="sqrt"),
        "linear_ecfp": LogisticRegression(C=1.0, max_iter=3000, random_state=seed),
        "knn_tanimoto": TanimotoKNN(k=5, task=task),
    }


class TanimotoKNN:
    """k-NN regressor/classifier under the exact Tanimoto metric on ECFP4 bits."""

    def __init__(self, k: int = 5, task: str = "regression"):
        self.k = k
        self.task = task

    def fit(self, X, y):
        self.X_ = np.ascontiguousarray(X, dtype=np.float32)
        self.pop_ = self.X_.sum(axis=1)
        self.y_ = np.asarray(y, dtype=np.float64)
        return self

    def _scores(self, Xq, chunk=4096):
        Xq = np.ascontiguousarray(Xq, dtype=np.float32)
        qa = Xq.sum(axis=1).astype(np.float64)
        out = np.empty(len(Xq))
        for s0 in range(0, len(Xq), chunk):
            q = Xq[s0:s0 + chunk]
            inter = (q @ self.X_.T).astype(np.float64)
            den = qa[s0:s0 + chunk, None] + self.pop_[None, :] - inter
            np.maximum(den, 1.0, out=den)
            sim = inter / den
            k = min(self.k, sim.shape[1])
            top = np.argpartition(-sim, k - 1, axis=1)[:, :k]
            w = np.take_along_axis(sim, top, axis=1)
            w = np.where(w.sum(axis=1, keepdims=True) > 0, w, 1.0)
            out[s0:s0 + chunk] = (w * self.y_[top]).sum(axis=1) / w.sum(axis=1)
        return out

    def predict(self, Xq):
        return self._scores(Xq)

    def predict_proba(self, Xq):
        p = np.clip(self._scores(Xq), 0.0, 1.0)
        return np.column_stack([1.0 - p, p])


MIN_SUBSET_N = 20   # pre-set: below this a subset metric is reported but not
                    # treated as an estimate (see EXPERIMENT_REPORT)

SIM_BINS = [0.0, 0.5, 0.7, 0.8, 0.9, 0.95, 0.999999, 1.0000001]
SIM_BIN_LABELS = ["<0.5", "0.5-0.7", "0.7-0.8", "0.8-0.9", "0.9-0.95",
                  "0.95-<1.0", "=1.0"]


def score(task, y_true, y_pred) -> float | None:
    if len(y_true) == 0:
        return None
    if task == "regression":
        return float(mean_absolute_error(y_true, y_pred))
    if len(set(np.asarray(y_true).tolist())) < 2:
        return None          # AUC undefined on a single-class subset
    return float(roc_auc_score(y_true, y_pred))


def similarity_response(task, te_y, pred_mean, nn_sim) -> dict:
    """Continuous dose-response: error vs nearest-neighbour similarity.

    The binary leaked/non-leaked split is degenerate when almost every test
    molecule is a near-duplicate (Lipophilicity-AZ has 4 non-leaked molecules).
    This per-molecule analysis keeps power in that regime and is reported
    alongside — never instead of — the pre-registered binary gap.
    """
    err = np.abs(pred_mean - te_y) if task == "regression" else np.abs(
        pred_mean - te_y.astype(float))
    rho = spearmanr(nn_sim, err)
    idx = np.digitize(nn_sim, SIM_BINS) - 1
    bins = []
    for k, lab in enumerate(SIM_BIN_LABELS):
        m = idx == k
        bins.append({
            "bin": lab, "n": int(m.sum()),
            "mean_abs_error": float(err[m].mean()) if m.sum() else None,
            "auc": (score(task, te_y[m], pred_mean[m])
                    if task == "classification" and m.sum() else None),
        })
    return {
        "statistic": "Spearman rho between NN-to-ChEMBL Tanimoto and "
                     "per-molecule absolute error (seed-averaged predictions)",
        "spearman_rho": float(rho.statistic),
        "spearman_p": float(rho.pvalue),
        "n_molecules": int(len(err)),
        "bins": bins,
    }


def _f(v) -> str:
    return "n/a" if v is None else f"{v:+.4f}"


def _mean_or_none(rows, key):
    vals = [r[key] for r in rows if r[key] is not None]
    return float(np.mean(vals)) if vals else None


def gap(task, s_leaked, s_nonleaked):
    """Positive gap == the leaked subset looks BETTER than the non-leaked one."""
    if s_leaked is None or s_nonleaked is None:
        return None
    return (s_leaked - s_nonleaked) if task == "classification" else (s_nonleaked - s_leaked)


def boot_ci(vals: np.ndarray, key: str, n_boot=10000, alpha=0.05):
    rng = np.random.default_rng(zlib.crc32(key.encode()))
    vals = np.asarray(vals, float)
    draws = rng.integers(0, len(vals), size=(n_boot, len(vals)))
    b = vals[draws].mean(axis=1)
    return float(np.quantile(b, alpha / 2)), float(np.quantile(b, 1 - alpha / 2))


def molecule_boot_ci(task, y_l, e_l, y_n, e_n, key, n_boot=10000, alpha=0.05):
    """Bootstrap the gap by resampling TEST MOLECULES within each subset.

    `e_*` are seed-averaged per-molecule absolute errors (regression) or
    seed-averaged predicted scores (classification).
    """
    rng = np.random.default_rng(zlib.crc32(key.encode()))
    nl, nn = len(e_l), len(e_n)
    out = np.empty(n_boot)
    for b in range(n_boot):
        il = rng.integers(0, nl, nl)
        inn = rng.integers(0, nn, nn)
        if task == "regression":
            out[b] = e_n[inn].mean() - e_l[il].mean()
        else:
            yl, yn = y_l[il], y_n[inn]
            if len(set(yl.tolist())) < 2 or len(set(yn.tolist())) < 2:
                out[b] = np.nan
                continue
            out[b] = roc_auc_score(yl, e_l[il]) - roc_auc_score(yn, e_n[inn])
    out = out[~np.isnan(out)]
    return float(np.quantile(out, alpha / 2)), float(np.quantile(out, 1 - alpha / 2)), len(out)


# --------------------------------------------------------------------------- #
def run_benchmark(b: str, nn_rows: list[dict], work: Path, seeds: list[int],
                  small: bool, log, pred_sink=None) -> dict:
    cfg = BENCHMARKS[b]
    task = cfg["task"]
    test = load_group_split(b, "test")
    trainval = load_group_split(b, "train_val")

    te_fp, te_keep, te_can = fp_matrix(test.smiles)
    te_y = np.array([test.y[i] for i in te_keep], float)
    by_can = {r["cansmi"]: r for r in nn_rows if r["benchmark"] == b and r["split"] == "test"}
    assert len(by_can) > 0, f"no S1 NN rows for {b} — FAIL"
    nn_sim = np.array([by_can[c]["nn_chembl_tanimoto"] for c in te_can], float)
    leaked = nn_sim >= LEAK_T
    test_can_set = set(te_can)

    log(f"{b}: test n={len(te_can)}  leaked={int(leaked.sum())}  "
        f"non-leaked={int((~leaked).sum())}")

    arms: dict[str, dict] = {}

    # ---- arm A: ChEMBL-derived training pools -----------------------------
    ep_file = ENDPOINT_FILE[b]
    if ep_file is not None:
        pool = [json.loads(l) for l in (work / ep_file).read_text().splitlines()]
        assert pool, f"empty ChEMBL pool {ep_file} — FAIL"
        p_fp, p_keep, p_can = fp_matrix([r["smiles"] for r in pool])
        p_y = np.array([pool[i]["y"] for i in p_keep], float)
        dup_mask = np.array([c in test_can_set for c in p_can])
        arms["A1_chembl_naive"] = {
            "fps": p_fp, "y": p_y, "task": task,
            "desc": "ChEMBL-derived pool, exact TDC-test duplicates NOT removed",
            "n_exact_test_duplicates_in_pool": int(dup_mask.sum()),
        }
        arms["A2_chembl_dedup"] = {
            "fps": p_fp[~dup_mask], "y": p_y[~dup_mask], "task": task,
            "desc": "ChEMBL-derived pool with exact TDC-test canonical-SMILES duplicates removed (PRIMARY)",
            "n_exact_test_duplicates_in_pool": int(dup_mask.sum()),
        }
        # sharper leakage definition: NN against the endpoint pool actually used
        pool_nn = exact_nn_tanimoto(te_fp, p_fp[~dup_mask], chunk=8192).max_sim
    else:
        pool_nn = None

    # ---- arm B: the benchmark's own train_val split -----------------------
    tv_fp, tv_keep, tv_can = fp_matrix(trainval.smiles)
    tv_y = np.array([trainval.y[i] for i in tv_keep], float)
    arms["B_tdc_trainval"] = {
        "fps": tv_fp, "y": tv_y, "task": task,
        "desc": "control: trained on the benchmark's own train_val split",
        "n_exact_test_duplicates_in_pool": int(sum(c in test_can_set for c in tv_can)),
    }

    te_X = unpack_f32(te_fp)
    results: dict = {}
    for arm_name, arm in arms.items():
        X = unpack_f32(arm["fps"])
        y = arm["y"]
        if task == "classification":
            y = y.astype(int)
        per_model: dict = {}
        for mname in MODEL_NAMES:
            per_seed = []
            pred_acc = np.zeros(len(te_y))
            for sd in seeds:
                rng = np.random.default_rng(zlib.crc32(f"{b}|{arm_name}|{mname}|{sd}".encode()))
                idx = rng.choice(len(y), size=int(0.9 * len(y)), replace=False)
                if task == "classification" and len(set(y[idx].tolist())) < 2:
                    continue
                m = make_models(task, sd, small)[mname]
                t0 = time.time()
                # Measured 2026-09-07 on the real logD pool (17,895 x 2048,
                # 200 trees): dense 13.7 s vs CSC 65.4 s.  Dense it is.
                Xtr, Xte = X[idx], te_X
                m.fit(Xtr, y[idx])
                pred = (m.predict(Xte) if task == "regression"
                        else m.predict_proba(Xte)[:, 1])
                pred_acc += pred
                s_l = score(task, te_y[leaked], pred[leaked])
                s_n = score(task, te_y[~leaked], pred[~leaked])
                s_all = score(task, te_y, pred)
                per_seed.append({
                    "seed": sd, "n_train": int(len(idx)),
                    "score_leaked": s_l, "score_nonleaked": s_n, "score_all": s_all,
                    "gap": gap(task, s_l, s_n), "fit_seconds": round(time.time() - t0, 1),
                })
                log(f"    {arm_name}/{mname} seed={sd} "
                    f"leaked={_f(s_l)} non-leaked={_f(s_n)} "
                    f"gap={_f(per_seed[-1]['gap'])} "
                    f"({per_seed[-1]['fit_seconds']}s)")
            assert per_seed, f"no successful fits for {b}/{arm_name}/{mname} — FAIL"
            gaps = np.array([r["gap"] for r in per_seed
                             if r["gap"] is not None], dtype=float)
            lo, hi = boot_ci(gaps, f"{b}|{arm_name}|{mname}|seedboot") if len(gaps) \
                else (None, None)
            pred_mean = pred_acc / len(per_seed)
            if task == "regression":
                err = np.abs(pred_mean - te_y)
                mlo, mhi, nb = molecule_boot_ci(
                    task, te_y[leaked], err[leaked], te_y[~leaked], err[~leaked],
                    f"{b}|{arm_name}|{mname}|molboot")
            else:
                mlo, mhi, nb = molecule_boot_ci(
                    task, te_y[leaked].astype(int), pred_mean[leaked],
                    te_y[~leaked].astype(int), pred_mean[~leaked],
                    f"{b}|{arm_name}|{mname}|molboot")
            if pred_sink is not None:
                for i, c in enumerate(te_can):
                    pred_sink.write(json.dumps({
                        "benchmark": b, "arm": arm_name, "model": mname,
                        "cansmi": c, "y_true": float(te_y[i]),
                        "pred_seed_mean": float(pred_mean[i]),
                        "nn_chembl_tanimoto": float(nn_sim[i]),
                        "leaked": bool(leaked[i]),
                        "n_seeds": len(per_seed),
                    }) + "\n")
            per_model[mname] = {
                "per_seed": per_seed,
                "n_seeds": len(per_seed),
                "mean_score_leaked": _mean_or_none(per_seed, "score_leaked"),
                "mean_score_nonleaked": _mean_or_none(per_seed, "score_nonleaked"),
                "mean_score_all": _mean_or_none(per_seed, "score_all"),
                "mean_gap": float(gaps.mean()) if len(gaps) else None,
                "sd_gap": float(gaps.std(ddof=1)) if len(gaps) > 1 else None,
                "gap_ci_seed_bootstrap": [lo, hi],
                "gap_ci_molecule_bootstrap": [mlo, mhi],
                "molecule_bootstrap_draws_used": nb,
                "seed_ci_excludes_zero": bool(lo is not None and (lo > 0 or hi < 0)),
                "molecule_ci_excludes_zero": bool(mlo > 0 or mhi < 0),
                "underpowered": bool(min(int(leaked.sum()), int((~leaked).sum()))
                                     < MIN_SUBSET_N),
                "similarity_response": similarity_response(
                    task, te_y, pred_mean, nn_sim),
            }
        results[arm_name] = {
            "description": arms[arm_name]["desc"],
            "n_train_pool": int(len(arm["y"])),
            "n_exact_test_duplicates_in_pool": arm["n_exact_test_duplicates_in_pool"],
            "models": per_model,
        }

    out = {
        "benchmark": b,
        "display": cfg["display"],
        "task": task,
        "metric": cfg["metric"],
        "gap_definition": ("AUC(leaked) - AUC(non-leaked)" if task == "classification"
                           else "MAE(non-leaked) - MAE(leaked)"),
        "leak_threshold": LEAK_T,
        "n_test": int(len(te_y)),
        "n_leaked": int(leaked.sum()),
        "n_nonleaked": int((~leaked).sum()),
        "seeds": seeds,
        "arms": results,
    }
    if pool_nn is not None:
        out["endpoint_pool_nn"] = {
            "note": "NN Tanimoto of each test molecule against the A2 ChEMBL "
                    "endpoint pool actually used for training (sharper than the "
                    "whole-corpus definition used for the leaked/non-leaked split)",
            "rate_ge_0p9": float((pool_nn >= 0.9).mean()),
            "rate_ge_0p95": float((pool_nn >= 0.95).mean()),
            "median": float(np.median(pool_nn)),
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="<REPO_ROOT>/.work/agenticls/full")
    ap.add_argument("--nn", default=str(PAPER / "evidence" / "s1_nn_per_molecule.jsonl"))
    ap.add_argument("--seeds", type=int, default=10)
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--out", default=str(PAPER / "evidence" / "s1b_gap_results.json"))
    args = ap.parse_args()

    work = Path(args.work)
    nn_rows = [json.loads(l) for l in Path(args.nn).read_text().splitlines()]
    assert nn_rows, "empty NN ledger — FAIL"

    # The pilot uses small forests but exercises EVERY code path: all three
    # benchmarks, both tasks (regression + classification), all arms, both model
    # families and both bootstrap routines.
    seeds = list(range(2 if args.pilot else args.seeds))
    benches = list(BENCHMARKS)

    t0 = time.time()
    def log(m):
        print(f"[{time.time()-t0:6.1f}s] {m}", flush=True)

    out = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "pilot": bool(args.pilot),
        "chembl_endpoints": json.loads((work / "chembl_endpoints.json").read_text()),
        "benchmarks": {},
    }
    pred_path = Path(args.out).with_name(Path(args.out).stem + "_predictions.jsonl")
    with open(pred_path, "w") as sink:
        for b in benches:
            out["benchmarks"][b] = run_benchmark(
                b, nn_rows, work, seeds, args.pilot, log, pred_sink=sink)
    out["predictions_file"] = str(pred_path)

    Path(args.out).write_text(json.dumps(out, indent=2) + "\n")
    log(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
