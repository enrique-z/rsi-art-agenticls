#!/usr/bin/env python3
"""review_round.py — the four analyses run in the 2026-09-07 analysis round.

Every one of these reuses artefacts that already exist on disk; none re-runs a
model-training sweep and none touches an existing ledger row.  Methodology is
imported from the original analysis modules wherever possible so that the new
numbers are produced by exactly the same code paths as the published ones.

A. ENDPOINT-POOL CENSUS  (C1)
   The published census defines a test molecule as leaked if it has a Tanimoto
   >= 0.9 neighbour ANYWHERE in ChEMBL 36.  Most such neighbours carry no
   measurement for the benchmark's endpoint and so cannot transfer a label.
   This recomputes the same census against the ENDPOINT-LABELLED pools only --
   the sets a curator for that endpoint can actually draw from -- giving the
   label-exploitable exposure rate beside the structural one.

B. MATCHING LADDER
   The published R3 matches on molecular-weight quintile x label quintile.  The
   MW-only variant was computed but never printed, and it turns out to leave the
   gap essentially unchanged: the entire shrinkage is the label quintile.  This
   computes the missing fourth rung -- label-only -- so the decomposition can be
   printed as a ladder (raw / MW-only / label-only / both) instead of collapsed
   into a single matched number.

C. NESTED BOOTSTRAP FOR THE MATCHED GAP  (C2)
   `s1b_rivals.cem_matched_gap` quantiles over redraws of the within-stratum
   down-sample with the matched population held fixed, so its interval covers
   matching variability and NOT sampling uncertainty over molecules.  This adds
   an OUTER bootstrap over test molecules -- resample, re-stratify, re-match,
   re-estimate -- and, because both arms are estimated on the same resampled
   molecules inside each outer draw, also yields a PAIRED interval for
   G(A2,matched) - G(B,matched), the quantity the paper previously (and wrongly)
   inferred from two non-overlapping marginal intervals.

D. CLUSTER-LEVEL LEAK-RATE READOUT  (C5)
   Section 6 asserts that benchmark error aggregated by chemical neighbourhood
   is itself partly a readout of which neighbourhoods the benchmark shares with
   ChEMBL.  That was a reading, not a measurement.  This reproduces the S2
   curation universe and its MiniBatchKMeans clustering exactly (same seed, same
   inputs), measures each cluster's test-reproduction rate, and reports (i) the
   rank correlation between that rate and the per-cluster benchmark error the
   curator was actually shown, and (ii) the reproduction rate of the clusters
   the LLM curators selected against the universe average.

Writes evidence/review_round_2026_09_07.json.  Nothing else is modified.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
import zlib
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from scipy.stats import spearmanr
from sklearn.cluster import MiniBatchKMeans
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bench_data import BENCHMARKS, load_group_split  # noqa: E402
from leakage_lib import (  # noqa: E402
    FP_BYTES,
    bemis_murcko,
    exact_nn_tanimoto,
    morgan_packed,
    unpack_f32,
    wilson_ci,
)

PAPER = Path(__file__).resolve().parents[1]
EV = PAPER / "evidence"
WORK = Path("<REPO_ROOT>/.work/agenticls/full")
LEAK_T = 0.9
N_OUTER = 1000          # outer bootstrap draws over molecules (C2)
N_INNER = 10            # matching redraws inside each outer draw (C2)
N_REP_POINT = 200       # matching redraws for the point estimate (matches s1b_rivals)

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


# --------------------------------------------------------------------------- #
# A. Endpoint-pool census
# --------------------------------------------------------------------------- #
def endpoint_pool_census() -> dict:
    out = {}
    for b, cfg in BENCHMARKS.items():
        ep = ENDPOINT_FILE[b]
        test = load_group_split(b, "test")
        te_fp, _, te_can = fp_matrix(test.smiles)
        if ep is None:
            out[b] = {
                "display": cfg["display"], "n_test": len(te_can),
                "available": False,
                "reason": "ChEMBL 36 carries no usable BBB endpoint pool (amendment A2)",
            }
            continue
        pool = [json.loads(l) for l in (WORK / ep).read_text().splitlines()]
        p_fp, _, _ = fp_matrix([x["smiles"] for x in pool])
        r = exact_nn_tanimoto(te_fp, p_fp, chunk=8192)
        n = len(te_can)
        k = int((r.max_sim >= LEAK_T).sum())
        lo, hi = wilson_ci(k, n)
        out[b] = {
            "display": cfg["display"], "n_test": n, "available": True,
            "pool_size": len(pool),
            "n_leaked_vs_endpoint_pool": k,
            "rate_0p9_vs_endpoint_pool": round(100.0 * k / n, 2),
            "ci_low": round(100.0 * lo, 2), "ci_high": round(100.0 * hi, 2),
            "median_nn_vs_endpoint_pool": round(float(np.median(r.max_sim)), 4),
        }
        print(f"[A] {b}: {k}/{n} = {out[b]['rate_0p9_vs_endpoint_pool']}% "
              f"vs endpoint pool of {len(pool)}", flush=True)
    return out


# --------------------------------------------------------------------------- #
# B/C. Matching ladder and nested bootstrap
# --------------------------------------------------------------------------- #
def _strata(mw, y, mode):
    """Coarsened-exact strata.  Quintile edges are recomputed on whatever sample
    is passed in, which is what makes the outer bootstrap in (C) valid."""
    def quint(v):
        q = np.quantile(v, np.linspace(0, 1, 6))
        return np.clip(np.searchsorted(q, v, side="right") - 1, 0, 4)
    if mode == "mw_only":
        return [(int(a),) for a in quint(mw)]
    if mode == "label_only":
        return [(int(c),) for c in quint(y)]
    if mode == "mw_and_label":
        return [(int(a), int(c)) for a, c in zip(quint(mw), quint(y))]
    raise ValueError(mode)


def _matched_gap_once(y, pred, leaked, strata, rng):
    """One matched draw: within each stratum down-sample the larger side to the
    smaller, then gap = MAE(clean) - MAE(leaked).  Same rule as
    s1b_rivals.cem_matched_gap."""
    by = defaultdict(lambda: {True: [], False: []})
    for i, s in enumerate(strata):
        by[s][bool(leaked[i])].append(i)
    usable = {s: v for s, v in by.items() if v[True] and v[False]}
    n_pairs = sum(min(len(v[True]), len(v[False])) for v in usable.values())
    if n_pairs < 10:
        return None, n_pairs
    li, ni = [], []
    for v in usable.values():
        k = min(len(v[True]), len(v[False]))
        li.extend(rng.choice(v[True], k, replace=False).tolist())
        ni.extend(rng.choice(v[False], k, replace=False).tolist())
    li, ni = np.array(li), np.array(ni)
    g = np.abs(pred[ni] - y[ni]).mean() - np.abs(pred[li] - y[li]).mean()
    return float(g), n_pairs


def matching_ladder_and_nested_bootstrap() -> dict:
    preds = [json.loads(l) for l in (EV / "s1b_gap_results_predictions.jsonl").read_text().splitlines()]
    b = "solubility_aqsoldb"
    rows_b = [p for p in preds if p["benchmark"] == b]
    arms = sorted({p["arm"] for p in rows_b})
    models = sorted({p["model"] for p in rows_b})

    # molecular weight, computed once per canonical SMILES
    cans = sorted({p["cansmi"] for p in rows_b})
    mw_of = {c: float(Descriptors.MolWt(Chem.MolFromSmiles(c))) for c in cans}

    ladder, nested, paired = {}, {}, {}
    per_arm_vectors = {}

    for arm in arms:
        for model in models:
            rs = [p for p in rows_b if p["arm"] == arm and p["model"] == model]
            if not rs:
                continue
            key = f"{arm}|{model}"
            y = np.array([p["y_true"] for p in rs], float)
            pred = np.array([p["pred_seed_mean"] for p in rs], float)
            leaked = np.array([p["leaked"] for p in rs], bool)
            mw = np.array([mw_of[p["cansmi"]] for p in rs], float)
            per_arm_vectors[key] = (y, pred, leaked, mw)

            raw = float(np.abs(pred[~leaked] - y[~leaked]).mean()
                        - np.abs(pred[leaked] - y[leaked]).mean())
            entry = {"raw_gap": round(raw, 4)}
            for mode in ("mw_only", "label_only", "mw_and_label"):
                rng = np.random.default_rng(zlib.crc32(f"{b}|{key}|{mode}".encode()))
                st = _strata(mw, y, mode)
                gs = []
                npair = 0
                for _ in range(N_REP_POINT):
                    g, npair = _matched_gap_once(y, pred, leaked, st, rng)
                    if g is not None:
                        gs.append(g)
                entry[mode] = {
                    "matched_gap": round(float(np.mean(gs)), 4),
                    "n_matched_per_side": int(npair),
                    "matching_draw_ci": [round(float(np.quantile(gs, 0.025)), 4),
                                         round(float(np.quantile(gs, 0.975)), 4)],
                }
            ladder[key] = entry
            print(f"[B] {key}: raw {entry['raw_gap']:+.4f} | mw {entry['mw_only']['matched_gap']:+.4f} "
                  f"| label {entry['label_only']['matched_gap']:+.4f} "
                  f"| both {entry['mw_and_label']['matched_gap']:+.4f}", flush=True)

    # ---- C. nested bootstrap over molecules, mw_and_label strata -----------
    # All arm/model cells share the same test molecules in the same order, so a
    # single outer index draw can be applied to every cell -- which is exactly
    # what makes the A2-minus-B difference PAIRED.
    order_key = f"{arms[0]}|{models[0]}"
    n_mol = len(per_arm_vectors[order_key][0])
    # sanity: every cell must be the same molecules in the same order
    ref_can = [p["cansmi"] for p in rows_b
               if p["arm"] == arms[0] and p["model"] == models[0]]
    for arm in arms:
        for model in models:
            c = [p["cansmi"] for p in rows_b if p["arm"] == arm and p["model"] == model]
            if c:
                assert c == ref_can, f"molecule order differs for {arm}|{model} — FAIL"

    rng = np.random.default_rng(zlib.crc32(f"{b}|nested_outer".encode()))
    boot = {k: [] for k in per_arm_vectors}
    diff_a2_b = {m: [] for m in models}
    t0 = time.time()
    for it in range(N_OUTER):
        idx = rng.integers(0, n_mol, size=n_mol)
        per_key_g = {}
        for key, (y, pred, leaked, mw) in per_arm_vectors.items():
            ys, ps, ls, ms = y[idx], pred[idx], leaked[idx], mw[idx]
            if ls.sum() < 20 or (~ls).sum() < 20:
                continue
            st = _strata(ms, ys, "mw_and_label")
            gs = [g for g in (_matched_gap_once(ys, ps, ls, st, rng)[0]
                              for _ in range(N_INNER)) if g is not None]
            if gs:
                per_key_g[key] = float(np.mean(gs))
                boot[key].append(per_key_g[key])
        for m in models:
            ka, kb = f"A2_chembl_dedup|{m}", f"B_tdc_trainval|{m}"
            if ka in per_key_g and kb in per_key_g:
                diff_a2_b[m].append(per_key_g[ka] - per_key_g[kb])
        if (it + 1) % 200 == 0:
            print(f"[C] outer draw {it+1}/{N_OUTER} ({time.time()-t0:.0f}s)", flush=True)

    for key, v in boot.items():
        v = np.array(v)
        nested[key] = {
            "n_outer_draws": int(len(v)),
            "nested_bootstrap_mean": round(float(v.mean()), 4),
            "nested_bootstrap_ci": [round(float(np.quantile(v, 0.025)), 4),
                                    round(float(np.quantile(v, 0.975)), 4)],
            "excludes_zero": bool(np.quantile(v, 0.025) > 0 or np.quantile(v, 0.975) < 0),
        }
    for m, v in diff_a2_b.items():
        v = np.array(v)
        paired[m] = {
            "n_outer_draws": int(len(v)),
            "paired_difference_mean": round(float(v.mean()), 4),
            "paired_difference_ci": [round(float(np.quantile(v, 0.025)), 4),
                                     round(float(np.quantile(v, 0.975)), 4)],
            "excludes_zero": bool(np.quantile(v, 0.025) > 0 or np.quantile(v, 0.975) < 0),
        }
        print(f"[C] A2-minus-B ({m}): {paired[m]['paired_difference_mean']:+.4f} "
              f"{paired[m]['paired_difference_ci']} excl0={paired[m]['excludes_zero']}", flush=True)

    return {"ladder": ladder, "nested_bootstrap": nested,
            "paired_a2_minus_b": paired,
            "config": {"n_outer": N_OUTER, "n_inner": N_INNER,
                       "n_rep_point_estimate": N_REP_POINT}}


# --------------------------------------------------------------------------- #
# D. Cluster-level leak-rate readout
# --------------------------------------------------------------------------- #
def cluster_leak_readout(n_clusters=40, budget=6000, trees=200, seeds=3) -> dict:
    BENCH = "solubility_aqsoldb"
    pool = [json.loads(l) for l in (WORK / ENDPOINT_FILE[BENCH]).read_text().splitlines()]
    p_fp, p_keep, p_can = fp_matrix([x["smiles"] for x in pool])
    p_y = np.array([pool[i]["y"] for i in p_keep], float)

    test = load_group_split(BENCH, "test")
    te_fp, te_keep, te_can = fp_matrix(test.smiles)
    te_y = np.array([test.y[i] for i in te_keep], float)
    nn_rows = {r["cansmi"]: r for l in (EV / "s1_nn_per_molecule.jsonl").read_text().splitlines()
               for r in [json.loads(l)] if r["benchmark"] == BENCH and r["split"] == "test"}
    nn_sim = np.array([nn_rows[c]["nn_chembl_tanimoto"] for c in te_can], float)
    leaked = nn_sim >= LEAK_T
    te_can_set = set(te_can)

    # reproduce the pre-committed hold-out and the curation universe exactly
    scaffolds = [bemis_murcko(c) or "" for c in p_can]
    uniq = sorted(set(scaffolds))
    held = {s for s in uniq
            if int(hashlib.sha256(("agenticls-rsi-leakage|" + s).encode())
                   .hexdigest()[:8], 16) % 5 == 0}
    is_held = np.array([s in held for s in scaffolds])
    is_test_dup = np.array([c in te_can_set for c in p_can])
    universe_idx = np.flatnonzero(~is_held & ~is_test_dup)

    commit = json.loads((EV / "s2_scaffold_holdout_commit.json").read_text())
    assert len(universe_idx) == commit["n_universe_compounds"], \
        "universe reconstruction differs from the committed S2 run — FAIL"

    U_dense = unpack_f32(p_fp[universe_idx])
    km = MiniBatchKMeans(n_clusters=n_clusters, random_state=0, n_init=5,
                         batch_size=2048).fit(U_dense)
    labels = km.labels_
    u_y = p_y[universe_idx]
    te_X = unpack_f32(te_fp)
    te_cluster = km.predict(te_X)

    # per-cluster TEST-REPRODUCTION RATE: fraction of the cluster's own compounds
    # that are a Tanimoto >= 0.9 near-duplicate of some TDC test molecule.  This
    # is "how much of the benchmark this cluster reproduces", i.e. what selecting
    # the cluster pulls into training.
    r = exact_nn_tanimoto(p_fp[universe_idx], te_fp, chunk=8192)
    u_is_dup = r.max_sim >= LEAK_T
    repro = np.array([float(u_is_dup[labels == c].mean()) for c in range(n_clusters)])
    csize = np.array([int((labels == c).sum()) for c in range(n_clusters)])
    universe_mean_repro = float(u_is_dup.mean())

    # per-cluster benchmark MAE actually shown to the curator at round 1, i.e.
    # computed from the round-0 model.  Round 0 is a fixed-seed random draw and
    # is reproduced here exactly.
    def draw(cids, seed_key):
        rng = np.random.default_rng(zlib.crc32(seed_key.encode()))
        cand = (np.flatnonzero(np.isin(labels, cids)) if cids is not None
                else np.arange(len(labels)))
        return rng.choice(cand, min(budget, len(cand)), replace=False)

    rhos, maes = [], []
    for sd in range(seeds):
        sel = draw(None, f"qwen2.5-coder:14b|{sd}|round0")
        m = RandomForestRegressor(n_estimators=trees, n_jobs=-1,
                                  random_state=sd, max_features="sqrt")
        m.fit(U_dense[sel], u_y[sel])
        p_te = m.predict(te_X)
        cl_mae = np.array([
            float(mean_absolute_error(te_y[te_cluster == c], p_te[te_cluster == c]))
            if (te_cluster == c).sum() else np.nan
            for c in range(n_clusters)])
        ok = ~np.isnan(cl_mae)
        rho = float(spearmanr(repro[ok], cl_mae[ok]).statistic)
        rhos.append(rho)
        maes.append(cl_mae)
        print(f"[D] seed {sd}: spearman(repro, cluster benchmark MAE) = {rho:+.4f} "
              f"over {int(ok.sum())} clusters", flush=True)

    # what the curators actually selected
    calls = [json.loads(l) for l in (EV / "s2_model_calls.jsonl").read_text().splitlines()]
    def weighted_repro(rows):
        sel_counts = np.zeros(n_clusters)
        for c in rows:
            for cid in c["selected"]:
                sel_counts[cid] += 1
        if sel_counts.sum() == 0:
            return None, None
        # weight each cluster by how often it was selected AND by its size, since
        # the draw is uniform over the union of selected clusters
        w = sel_counts * csize
        return float((w * repro).sum() / w.sum()), sel_counts

    llm_rows = [c for c in calls if c["curator"] != "random_control"]
    rnd_rows = [c for c in calls if c["curator"] == "random_control"]
    llm_repro, llm_counts = weighted_repro(llm_rows)
    rnd_repro, rnd_counts = weighted_repro(rnd_rows)
    print(f"[D] selection-weighted reproduction rate — LLM {llm_repro:.4f} | "
          f"random {rnd_repro:.4f} | universe mean {universe_mean_repro:.4f}", flush=True)

    return {
        "n_clusters": n_clusters,
        "n_universe_compounds": int(len(universe_idx)),
        "universe_mean_test_reproduction_rate": round(universe_mean_repro, 4),
        "per_cluster_test_reproduction_rate": [round(float(x), 4) for x in repro],
        "per_cluster_size": csize.tolist(),
        "spearman_repro_vs_cluster_benchmark_mae": {
            "per_seed": [round(x, 4) for x in rhos],
            "mean": round(float(np.mean(rhos)), 4),
        },
        "selection_weighted_reproduction_rate": {
            "llm_curators": round(llm_repro, 4),
            "random_control": round(rnd_repro, 4),
            "universe_mean": round(universe_mean_repro, 4),
            "note": ("the random control's value equals the universe mean by "
                     "construction; the informative contrast is LLM vs universe mean"),
        },
    }


def main() -> int:
    out = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "leak_threshold": LEAK_T,
        "purpose": "2026-09-07 analysis round: C1 endpoint-pool census, "
                   "matching ladder, C2 nested bootstrap, C5 cluster leak readout",
    }
    print("=== A. endpoint-pool census ===", flush=True)
    out["A_endpoint_pool_census"] = endpoint_pool_census()
    print("=== B/C. matching ladder + nested bootstrap ===", flush=True)
    out.update({"BC_matching": matching_ladder_and_nested_bootstrap()})
    print("=== D. cluster leak readout ===", flush=True)
    out["D_cluster_leak_readout"] = cluster_leak_readout()
    p = EV / "review_round_2026_09_07.json"
    p.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", p, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
