#!/usr/bin/env python3
"""s2_shuffled_feedback.py — C7, the feedback-shuffled control arm.

Design, primary quantities and the reporting rule for EVERY outcome were fixed
in `C7_PRECOMMITMENT_2026-09-07.md` before this file was written; its SHA-256 is
recorded in the output so the rule is verifiable against the version that
existed before the first model call.

The arm is identical to the S2 LLM arms in every respect -- same curators, same
temperature, same seeds and therefore the same round-0 draw, same universe, same
clustering, same budget, same specialized model -- except that the per-cluster
`benchmark MAE` column is RANDOMLY PERMUTED ACROSS CLUSTERS before the prompt is
built.  Cluster ids are not permuted, so selections still address real clusters
and the drawn pool is real chemistry; only the signal about which clusters need
work is destroyed.

Everything structural is imported from or reproduced exactly against
`s2_rsi_loop.py`, and the pre-committed hold-out SHA-256 is re-asserted at run
time, so a drift in the universe fails loudly instead of quietly changing what
"the same arm" means.

Writes evidence/s2_shuffled_feedback.json and
evidence/s2_shuffled_model_calls.jsonl.  Modifies nothing else.
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
import zlib
from pathlib import Path

import numpy as np
from sklearn.cluster import MiniBatchKMeans
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bench_data import load_group_split  # noqa: E402
from leakage_lib import FP_BYTES, bemis_murcko, morgan_packed, unpack_f32  # noqa: E402
from s2_rsi_loop import build_prompt, ollama, parse_clusters  # noqa: E402

PAPER = Path(__file__).resolve().parents[1]
EV = PAPER / "evidence"
WORK = Path("<REPO_ROOT>/.work/agenticls/full")
BENCH = "solubility_aqsoldb"
POOL_FILE = "chembl_endpoint_solubility.jsonl"
LEAK_T = 0.9
CURATORS = ["qwen2.5-coder:14b", "qwen3:8b"]
SEEDS, ROUNDS, CLUSTERS, BUDGET, TREES = 3, 3, 40, 6000, 200
PRECOMMIT = PAPER / "C7_PRECOMMITMENT_2026-09-07.md"


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


def jaccard(a, b) -> float:
    A, B = set(a), set(b)
    return len(A & B) / len(A | B) if (A | B) else 1.0


def main() -> int:
    assert PRECOMMIT.exists(), "pre-commitment file missing — refusing to run"
    precommit_sha = hashlib.sha256(PRECOMMIT.read_bytes()).hexdigest()
    print(f"pre-commitment SHA-256 {precommit_sha}", flush=True)

    # ---- rebuild the S2 universe exactly ----------------------------------
    pool = [json.loads(l) for l in (WORK / POOL_FILE).read_text().splitlines()]
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

    scaffolds = [bemis_murcko(c) or "" for c in p_can]
    uniq = sorted(set(scaffolds))
    held = {s for s in uniq
            if int(hashlib.sha256(("agenticls-rsi-leakage|" + s).encode())
                   .hexdigest()[:8], 16) % 5 == 0}
    is_held = np.array([s in held for s in scaffolds])
    is_test_dup = np.array([c in te_can_set for c in p_can])
    holdout_idx = np.flatnonzero(is_held & ~is_test_dup)
    universe_idx = np.flatnonzero(~is_held & ~is_test_dup)

    commit = json.loads((EV / "s2_scaffold_holdout_commit.json").read_text())
    ho_sha = hashlib.sha256(
        "\n".join(sorted(p_can[i] for i in holdout_idx)).encode()).hexdigest()
    assert ho_sha == commit["holdout_sha256"], (
        "pre-committed hold-out does NOT reproduce — refusing to run C7")
    assert len(universe_idx) == commit["n_universe_compounds"], "universe drift — FAIL"
    print(f"hold-out SHA-256 re-asserted: {ho_sha[:16]}…", flush=True)

    U_dense = unpack_f32(p_fp[universe_idx])
    km = MiniBatchKMeans(n_clusters=CLUSTERS, random_state=0, n_init=5,
                         batch_size=2048).fit(U_dense)
    labels = km.labels_
    u_y = p_y[universe_idx]
    u_can = [p_can[i] for i in universe_idx]
    te_X = unpack_f32(te_fp)
    te_cluster = km.predict(te_X)
    ho_X = unpack_f32(p_fp[holdout_idx])
    ho_y = p_y[holdout_idx]

    def fit_eval(sel_idx, seed):
        m = RandomForestRegressor(n_estimators=TREES, n_jobs=-1,
                                  random_state=seed, max_features="sqrt")
        m.fit(U_dense[sel_idx], u_y[sel_idx])
        p_te, p_ho = m.predict(te_X), m.predict(ho_X)
        per_cluster = {}
        for c in range(CLUSTERS):
            mask = te_cluster == c
            per_cluster[c] = (float(mean_absolute_error(te_y[mask], p_te[mask]))
                              if mask.sum() else None)
        return {
            "mae_full_test": float(mean_absolute_error(te_y, p_te)),
            "mae_nonleaked_test": float(mean_absolute_error(te_y[~leaked], p_te[~leaked])),
            "mae_leaked_test": float(mean_absolute_error(te_y[leaked], p_te[leaked])),
            "mae_scaffold_holdout": float(mean_absolute_error(ho_y, p_ho)),
        }, per_cluster

    def draw_from_clusters(cids, seed_key):
        rng = np.random.default_rng(zlib.crc32(seed_key.encode()))
        cand = (np.flatnonzero(np.isin(labels, cids)) if cids is not None
                else np.arange(len(labels)))
        return rng.choice(cand, min(BUDGET, len(cand)), replace=False)

    cluster_meta = []
    for c in range(CLUSTERS):
        idx = np.flatnonzero(labels == c)
        cluster_meta.append({
            "id": c, "n": int(len(idx)),
            "mean_y": float(u_y[idx].mean()), "sd_y": float(u_y[idx].std()),
            "example": u_can[int(idx[0])],
        })

    # real-feedback selections, for the paired selection-agreement statistic
    real_calls = [json.loads(l) for l in (EV / "s2_model_calls.jsonl").read_text().splitlines()]
    real_sel = {(c["arm"], c["seed"], c["round"]): c["selected"]
                for c in real_calls if c["curator"] != "random_control"}

    out = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "benchmark": BENCH,
        "arm": "shuffled_feedback",
        "precommitment_file": PRECOMMIT.name,
        "precommitment_sha256": precommit_sha,
        "manipulation": ("per-cluster benchmark MAE permuted across clusters before the "
                         "prompt is built; every other field byte-identical to the "
                         "real-feedback arm; cluster ids not permuted"),
        "config": {"rounds": ROUNDS, "seeds": SEEDS, "clusters": CLUSTERS,
                   "budget": BUDGET, "trees": TREES, "curators": CURATORS},
        "scaffold_holdout_sha256": ho_sha,
        "n_test": len(te_can), "n_leaked": int(leaked.sum()),
        "n_nonleaked": int((~leaked).sum()),
        "arms": {},
    }
    calls_out = []
    n_fallback = 0

    for arm_id in CURATORS:
        per_seed = []
        for sd in range(SEEDS):
            # IDENTICAL round-0 draw to the real-feedback arm: same seed key
            sel = draw_from_clusters(None, f"{arm_id}|{sd}|round0")
            m0, pc = fit_eval(sel, sd)
            rounds = [{"round": 0, "selected_clusters": None, **m0}]
            print(f"  [shuf] {arm_id} seed={sd} r0 full={m0['mae_full_test']:.4f} "
                  f"clean={m0['mae_nonleaked_test']:.4f} ho={m0['mae_scaffold_holdout']:.4f}",
                  flush=True)
            for r in range(1, ROUNDS + 1):
                rows = [dict(cm, bench_mae=pc[cm["id"]]) for cm in cluster_meta]
                # ---- THE MANIPULATION -------------------------------------
                vals = [x["bench_mae"] for x in rows]
                rng = np.random.default_rng(
                    zlib.crc32(f"shuffle|{arm_id}|{sd}|{r}".encode()))
                perm = rng.permutation(len(vals))
                for k, x in enumerate(rows):
                    x["bench_mae"] = vals[perm[k]]
                # -----------------------------------------------------------
                fb = {"benchmark_mae": rounds[-1]["mae_full_test"],
                      "prev_mae": rounds[-2]["mae_full_test"] if len(rounds) > 1 else None}
                prompt = build_prompt(r, rows, BUDGET, fb)
                resp = ollama(arm_id, prompt)
                cids = parse_clusters(resp.get("response", ""), CLUSTERS)
                ok = bool(cids)
                if not ok:
                    n_fallback += 1
                    cids = rounds[-1]["selected_clusters"] or list(range(CLUSTERS))
                calls_out.append({
                    "arm": arm_id, "seed": sd, "round": r, "curator": arm_id,
                    "feedback": "shuffled", "ok": ok, "selected": cids,
                    "latency_s": resp.get("latency_s"),
                    "eval_count": resp.get("eval_count"),
                    "prompt_sha256": hashlib.sha256(prompt.encode()).hexdigest(),
                    "prompt": prompt, "response": resp.get("response", "")[:4000],
                    "error": resp.get("error"),
                    "shuffle_permutation": perm.tolist(),
                    "jaccard_vs_real": jaccard(cids, real_sel.get((arm_id, sd, r), [])),
                })
                sel = draw_from_clusters(cids, f"{arm_id}|{sd}|round{r}|shuf")
                mr, pc = fit_eval(sel, sd)
                rounds.append({"round": r, "selected_clusters": cids, **mr})
                print(f"  [shuf] {arm_id} seed={sd} r{r} full={mr['mae_full_test']:.4f} "
                      f"clean={mr['mae_nonleaked_test']:.4f} ho={mr['mae_scaffold_holdout']:.4f} "
                      f"J={calls_out[-1]['jaccard_vs_real']:.2f} n_sel={len(cids)}", flush=True)
            per_seed.append({
                "seed": sd, "rounds": rounds,
                "gain_full_test": rounds[0]["mae_full_test"] - rounds[-1]["mae_full_test"],
                "gain_nonleaked_test": rounds[0]["mae_nonleaked_test"] - rounds[-1]["mae_nonleaked_test"],
                "gain_scaffold_holdout": rounds[0]["mae_scaffold_holdout"] - rounds[-1]["mae_scaffold_holdout"],
            })
        gf = np.array([r["gain_full_test"] for r in per_seed])
        gn = np.array([r["gain_nonleaked_test"] for r in per_seed])
        gh = np.array([r["gain_scaffold_holdout"] for r in per_seed])
        out["arms"][arm_id] = {
            "curator": arm_id, "per_seed": per_seed,
            "mean_gain_full_test": float(gf.mean()),
            "mean_gain_nonleaked_test": float(gn.mean()),
            "mean_gain_scaffold_holdout": float(gh.mean()),
            "mean_overstatement_full_minus_nonleaked": float((gf - gn).mean()),
        }

    out["n_fallbacks"] = n_fallback
    out["n_calls"] = len(calls_out)

    # ---- paired contrasts against the real-feedback arm --------------------
    real = json.loads((EV / "s2_rsi_results.json").read_text())

    def boot(d, key):
        d = np.array(d, float)
        rng = np.random.default_rng(zlib.crc32(key.encode()))
        m = np.array([rng.choice(d, len(d), replace=True).mean() for _ in range(10000)])
        return float(d.mean()), [float(np.quantile(m, .025)), float(np.quantile(m, .975))]

    pairs = {"gain_scaffold_holdout": [], "gain_full_test": [],
             "gain_nonleaked_test": [], "overstatement": []}
    for arm_id in CURATORS:
        rs = {s["seed"]: s for s in real["arms"][arm_id]["per_seed"]}
        for s in out["arms"][arm_id]["per_seed"]:
            r0 = rs[s["seed"]]
            pairs["gain_scaffold_holdout"].append(
                r0["gain_scaffold_holdout"] - s["gain_scaffold_holdout"])
            pairs["gain_full_test"].append(r0["gain_full_test"] - s["gain_full_test"])
            pairs["gain_nonleaked_test"].append(
                r0["gain_nonleaked_test"] - s["gain_nonleaked_test"])
            pairs["overstatement"].append(
                (r0["gain_full_test"] - r0["gain_nonleaked_test"])
                - (s["gain_full_test"] - s["gain_nonleaked_test"]))
    out["paired_real_minus_shuffled"] = {}
    for k, v in pairs.items():
        mean, ci = boot(v, f"c7|{k}")
        out["paired_real_minus_shuffled"][k] = {
            "n_units": len(v), "unit": "curator x seed", "per_unit": v,
            "mean": mean, "ci": ci,
            "ci_excludes_zero": bool(ci[0] > 0 or ci[1] < 0),
        }
        print(f"[C7] real-minus-shuffled {k}: {mean:+.4f} "
              f"[{ci[0]:+.4f}, {ci[1]:+.4f}] excl0={out['paired_real_minus_shuffled'][k]['ci_excludes_zero']}",
              flush=True)

    js = [c["jaccard_vs_real"] for c in calls_out]
    jm, jci = boot(js, "c7|jaccard")
    out["selection_agreement"] = {
        "statistic": "Jaccard(real-feedback selection, shuffled-feedback selection)",
        "n_calls": len(js), "mean": jm, "ci": jci, "per_call": js,
    }
    print(f"[C7] selection agreement Jaccard = {jm:.4f} [{jci[0]:.4f}, {jci[1]:.4f}] "
          f"over {len(js)} matched calls", flush=True)

    (EV / "s2_shuffled_feedback.json").write_text(json.dumps(out, indent=2) + "\n")
    with open(EV / "s2_shuffled_model_calls.jsonl", "w") as fh:
        for c in calls_out:
            fh.write(json.dumps(c) + "\n")
    print(f"wrote evidence/s2_shuffled_feedback.json ({len(calls_out)} calls, "
          f"{n_fallback} fallbacks)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
