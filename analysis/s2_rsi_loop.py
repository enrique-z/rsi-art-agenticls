#!/usr/bin/env python3
"""s2_rsi_loop.py — Stage 2: are LLM-curator RSI gains leakage-driven?

A 3-round recursive-self-improvement loop in which a general-purpose LLM curates
the ChEMBL training pool for a specialised ADMET model, guided by feedback from
the TDC benchmark, and the resulting improvement is measured three ways:

  (i)   the FULL TDC test set                 -- the number an RSI paper reports
  (ii)  the NON-LEAKED TDC test subset        -- NN-to-ChEMBL Tanimoto < 0.9
  (iii) a PRE-COMMITTED ChEMBL scaffold hold-out, whose Bemis-Murcko scaffolds
        are excluded from the curation universe entirely

Pre-registered claim: the apparent gain on (i) overstates the gain on (ii).

Design
------
* Curation universe: the ChEMBL solubility pool minus the scaffold hold-out,
  partitioned into `--clusters` chemically coherent clusters by MiniBatchKMeans
  on the same Morgan r=2 / 2048-bit fingerprints.
* Round 0 is a fixed-seed RANDOM selection of `--budget` compounds (baseline).
* Rounds 1..3: the curator is shown, per cluster, its size, its label statistics,
  example SMILES, and **the current model's mean absolute error on the TDC test
  molecules nearest to that cluster** -- i.e. real benchmark feedback -- and
  returns the cluster ids to draw the next training pool from.
* Budget (compounds per round) and number of curator calls are IDENTICAL across
  every arm, including the random-curator control arm.
* Every LLM call (prompt, raw response, latency, parse outcome) is logged.

Nothing is simulated: every metric is a real fit evaluated on real molecules.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
import zlib
from pathlib import Path

import numpy as np
from sklearn.cluster import MiniBatchKMeans
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bench_data import load_group_split  # noqa: E402
from leakage_lib import (  # noqa: E402
    FP_BYTES,
    bemis_murcko,
    exact_nn_tanimoto,
    morgan_packed,
    unpack_f32,
)

os.environ.setdefault("PYTHONHASHSEED", "0")
PAPER = Path(__file__).resolve().parents[1]
OLLAMA = "http://localhost:11434/api/generate"
LEAK_T = 0.9
# S1 measured only 4 non-leaked molecules in the Lipophilicity-AZ test split and
# 25 in BBB-Martins.  Solubility-AqSolDB (1,013 leaked / 982 non-leaked) is the
# ONLY benchmark on which "gain on the non-leaked subset" is estimable at all, so
# S2 runs there.  Recorded as amendment A7.
BENCH = "solubility_aqsoldb"
ENDPOINT_POOL = "chembl_endpoint_solubility.jsonl"
ENDPOINT_LABEL = "log S (aqueous solubility, mol/L)"


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


# Measured 2026-09-07 on this machine with a 12-cluster probe prompt:
#   qwen2.5-coder:14b            22.7 s, returned a well-formed selection
#   qwen3:8b   (think disabled)   1.0 s
#   qwen3:8b   (thinking on)     86.9 s, spent the whole budget inside <think>
#   gpt-oss:20b(think disabled) 277.6 s  -> too slow for 27 calls, not used
# Reasoning models are therefore called with thinking OFF and a hard token cap.
THINKING_MODELS = ("qwen3", "gpt-oss", "deepseek-r1")


def ollama(model: str, prompt: str, timeout: int = 420) -> dict:
    payload = {
        "model": model, "prompt": prompt, "stream": False,
        "options": {"temperature": 0.0, "num_ctx": 8192, "num_predict": 800},
    }
    if any(model.startswith(t) for t in THINKING_MODELS):
        payload["think"] = False
    body = json.dumps(payload).encode()
    req = urllib.request.Request(OLLAMA, data=body,
                                 headers={"Content-Type": "application/json"})
    t0 = time.time()
    try:
        with urllib.request.urlopen(req, timeout=timeout) as fh:
            raw = json.loads(fh.read().decode())
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        return {"ok": False, "error": repr(e), "latency_s": round(time.time() - t0, 1),
                "response": ""}
    txt = raw.get("response", "")
    assert isinstance(txt, str), "ollama returned a non-string response — FAIL"
    return {"ok": True, "response": txt, "latency_s": round(time.time() - t0, 1),
            "eval_count": raw.get("eval_count"), "model": raw.get("model")}


def parse_clusters(text: str, n_clusters: int) -> list[int]:
    """Extract cluster ids from the curator's reply.  Empty list == parse failure."""
    m = re.findall(r"\[[^\[\]]*\]", text.replace("\n", " "))
    for chunk in reversed(m):
        try:
            v = json.loads(chunk)
        except json.JSONDecodeError:
            v = [int(x) for x in re.findall(r"-?\d+", chunk)]
        ids = [int(x) for x in v if isinstance(x, (int, float))
               and 0 <= int(x) < n_clusters]
        if ids:
            return sorted(set(ids))
    ids = [int(x) for x in re.findall(r"\b\d+\b", text) if 0 <= int(x) < n_clusters]
    return sorted(set(ids))


def build_prompt(round_i, cluster_rows, budget, feedback) -> str:
    lines = [
        "You are curating a training set for a machine-learning model that "
        f"predicts {ENDPOINT_LABEL} of small molecules.",
        "",
        f"You may select AT MOST {budget} compounds in total for the next training "
        "round. Compounds are grouped into chemically similar clusters. Pick the "
        "clusters whose compounds will most improve the model's score on the "
        "held-out benchmark.",
        "",
        f"Current benchmark mean absolute error: {feedback['benchmark_mae']:.4f} "
        f"(lower is better). Previous round MAE: "
        + (f"{feedback['prev_mae']:.4f}" if feedback["prev_mae"] is not None else "n/a"),
        "",
        "CLUSTERS (id | n_compounds | mean_label | sd_label | benchmark MAE on the "
        "benchmark molecules closest to this cluster | example SMILES):",
    ]
    for r in cluster_rows:
        lines.append(
            f"{r['id']} | {r['n']} | {r['mean_y']:.2f} | {r['sd_y']:.2f} | "
            + (f"{r['bench_mae']:.3f}" if r["bench_mae"] is not None else "n/a")
            + f" | {r['example']}"
        )
    lines += [
        "",
        f"Reply with ONLY a JSON array of the cluster ids to include, e.g. [0, 3, 7]. "
        f"Selecting clusters totalling more than {budget} compounds is allowed; the "
        "pool will be truncated to the budget.",
    ]
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="<REPO_ROOT>/.work/agenticls/full")
    ap.add_argument("--nn", default=str(PAPER / "evidence" / "s1_nn_per_molecule.jsonl"))
    ap.add_argument("--curators", default="qwen2.5-coder:14b,qwen3:8b")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--clusters", type=int, default=40)
    ap.add_argument("--budget", type=int, default=6000)
    ap.add_argument("--trees", type=int, default=200)
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--out", default=str(PAPER / "evidence" / "s2_rsi_results.json"))
    ap.add_argument("--calls", default=str(PAPER / "evidence" / "s2_model_calls.jsonl"))
    args = ap.parse_args()

    if args.pilot:
        args.seeds, args.rounds, args.clusters, args.budget, args.trees = 1, 1, 12, 1500, 40
        args.curators = args.curators.split(",")[0]

    work = Path(args.work)
    t_start = time.time()

    # ---- data -------------------------------------------------------------
    pool = [json.loads(l) for l in
            (work / ENDPOINT_POOL).read_text().splitlines()]
    assert pool, f"empty ChEMBL pool {ENDPOINT_POOL} — FAIL"
    p_fp, p_keep, p_can = fp_matrix([r["smiles"] for r in pool])
    p_y = np.array([pool[i]["y"] for i in p_keep], float)

    test = load_group_split(BENCH, "test")
    te_fp, te_keep, te_can = fp_matrix(test.smiles)
    te_y = np.array([test.y[i] for i in te_keep], float)
    nn_rows = {r["cansmi"]: r for l in Path(args.nn).read_text().splitlines()
               for r in [json.loads(l)] if r["benchmark"] == BENCH and r["split"] == "test"}
    nn_sim = np.array([nn_rows[c]["nn_chembl_tanimoto"] for c in te_can], float)
    leaked = nn_sim >= LEAK_T
    te_can_set = set(te_can)

    # ---- PRE-COMMITTED scaffold hold-out (deterministic, hashed) ----------
    scaffolds = [bemis_murcko(c) or "" for c in p_can]
    uniq = sorted(set(scaffolds))
    held = {s for s in uniq
            if int(hashlib.sha256(("agenticls-rsi-leakage|" + s).encode())
                   .hexdigest()[:8], 16) % 5 == 0}          # deterministic 20%
    is_held = np.array([s in held for s in scaffolds])
    # a compound that is also an exact TDC test molecule never enters either side
    is_test_dup = np.array([c in te_can_set for c in p_can])
    holdout_idx = np.flatnonzero(is_held & ~is_test_dup)
    universe_idx = np.flatnonzero(~is_held & ~is_test_dup)
    assert len(holdout_idx) > 200 and len(universe_idx) > 2000, "degenerate split — FAIL"
    commit = {
        "rule": "sha256('agenticls-rsi-leakage|' + bemis_murcko_smiles)[:8] % 5 == 0",
        "n_scaffolds_total": len(uniq),
        "n_scaffolds_heldout": len(held),
        "n_holdout_compounds": int(len(holdout_idx)),
        "n_universe_compounds": int(len(universe_idx)),
        "n_excluded_exact_tdc_test_duplicates": int(is_test_dup.sum()),
        "holdout_sha256": hashlib.sha256(
            "\n".join(sorted(p_can[i] for i in holdout_idx)).encode()).hexdigest(),
    }
    (PAPER / "evidence" / "s2_scaffold_holdout_commit.json").write_text(
        json.dumps(commit, indent=2) + "\n")
    print("scaffold hold-out committed:", json.dumps(commit), flush=True)

    # ---- clustering of the curation universe ------------------------------
    U_dense = unpack_f32(p_fp[universe_idx])
    km = MiniBatchKMeans(n_clusters=args.clusters, random_state=0, n_init=5,
                         batch_size=2048).fit(U_dense)
    labels = km.labels_
    u_y = p_y[universe_idx]
    u_can = [p_can[i] for i in universe_idx]

    te_X = unpack_f32(te_fp)
    te_cluster = km.predict(te_X)
    ho_X = unpack_f32(p_fp[holdout_idx])
    ho_y = p_y[holdout_idx]

    def fit_eval(sel_idx, seed):
        m = RandomForestRegressor(n_estimators=args.trees, n_jobs=-1,
                                  random_state=seed, max_features="sqrt")
        m.fit(U_dense[sel_idx], u_y[sel_idx])
        p_te = m.predict(te_X)
        p_ho = m.predict(ho_X)
        per_cluster = {}
        for c in range(args.clusters):
            mask = te_cluster == c
            per_cluster[c] = (float(mean_absolute_error(te_y[mask], p_te[mask]))
                              if mask.sum() >= 3 else None)
        return {
            "mae_full_test": float(mean_absolute_error(te_y, p_te)),
            "mae_nonleaked_test": float(mean_absolute_error(te_y[~leaked], p_te[~leaked])),
            "mae_leaked_test": float(mean_absolute_error(te_y[leaked], p_te[leaked])),
            "mae_scaffold_holdout": float(mean_absolute_error(ho_y, p_ho)),
            "n_train": int(len(sel_idx)),
        }, per_cluster

    def draw_from_clusters(cids, seed_key):
        rng = np.random.default_rng(zlib.crc32(seed_key.encode()))
        pool_ids = np.flatnonzero(np.isin(labels, list(cids))) if cids else np.arange(len(labels))
        if len(pool_ids) > args.budget:
            pool_ids = rng.choice(pool_ids, args.budget, replace=False)
        return np.sort(pool_ids)

    cluster_meta = []
    for c in range(args.clusters):
        idx = np.flatnonzero(labels == c)
        cluster_meta.append({
            "id": c, "n": int(len(idx)),
            "mean_y": float(u_y[idx].mean()), "sd_y": float(u_y[idx].std()),
            "example": u_can[int(idx[0])][:90],
        })

    calls_fh = open(args.calls, "w")
    arms = [("random", "random_curator_control")]
    arms += [(m, m) for m in args.curators.split(",") if m.strip()]
    results = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "benchmark": BENCH, "pilot": bool(args.pilot),
        "config": {"rounds": args.rounds, "seeds": args.seeds,
                   "clusters": args.clusters, "budget": args.budget,
                   "trees": args.trees, "curators": [a[0] for a in arms]},
        "scaffold_holdout_commit": commit,
        "n_test": int(len(te_y)), "n_leaked": int(leaked.sum()),
        "n_nonleaked": int((~leaked).sum()),
        "arms": {},
    }

    for arm_id, arm_name in arms:
        per_seed = []
        for sd in range(args.seeds):
            sel = draw_from_clusters(None, f"{arm_id}|{sd}|round0")
            m0, pc = fit_eval(sel, sd)
            rounds = [{"round": 0, "selected_clusters": None,
                       "curator": "random baseline", **m0}]
            print(f"  {arm_name} seed={sd} r0 full={m0['mae_full_test']:.4f} "
                  f"nonleaked={m0['mae_nonleaked_test']:.4f} "
                  f"holdout={m0['mae_scaffold_holdout']:.4f}", flush=True)
            prev = None
            for r in range(1, args.rounds + 1):
                rows = [dict(cm, bench_mae=pc[cm["id"]]) for cm in cluster_meta]
                fb = {"benchmark_mae": rounds[-1]["mae_full_test"], "prev_mae": prev}
                if arm_id == "random":
                    rng = np.random.default_rng(zlib.crc32(f"rand|{sd}|{r}".encode()))
                    k = max(1, args.clusters // 2)
                    cids = sorted(rng.choice(args.clusters, k, replace=False).tolist())
                    call = {"arm": arm_id, "seed": sd, "round": r,
                            "curator": "random_control", "ok": True,
                            "selected": cids, "latency_s": 0.0,
                            "note": "no LLM call; identical budget and call count"}
                else:
                    prompt = build_prompt(r, rows, args.budget, fb)
                    resp = ollama(arm_id, prompt)
                    # An explicit empty selection is a curation failure, not a
                    # valid "select nothing" — it is counted as such below.
                    cids = parse_clusters(resp.get("response", ""), args.clusters)
                    call = {"arm": arm_id, "seed": sd, "round": r, "curator": arm_id,
                            "ok": bool(resp.get("ok")) and bool(cids),
                            "latency_s": resp.get("latency_s"),
                            "eval_count": resp.get("eval_count"),
                            "prompt_sha256": hashlib.sha256(prompt.encode()).hexdigest(),
                            "prompt": prompt, "response": resp.get("response", "")[:4000],
                            "selected": cids, "error": resp.get("error")}
                    if not cids:
                        # parse failure: recorded, arm continues with the previous
                        # selection so the budget stays identical across arms
                        cids = rounds[-1]["selected_clusters"] or list(range(args.clusters))
                        call["fallback_used"] = True
                calls_fh.write(json.dumps(call) + "\n")
                calls_fh.flush()
                sel = draw_from_clusters(cids, f"{arm_id}|{sd}|round{r}")
                prev = rounds[-1]["mae_full_test"]
                mr, pc = fit_eval(sel, sd)
                rounds.append({"round": r, "selected_clusters": cids,
                               "curator": arm_id, "llm_ok": call["ok"], **mr})
                print(f"  {arm_name} seed={sd} r{r} full={mr['mae_full_test']:.4f} "
                      f"nonleaked={mr['mae_nonleaked_test']:.4f} "
                      f"holdout={mr['mae_scaffold_holdout']:.4f} "
                      f"({len(cids)} clusters)", flush=True)
            last = rounds[-1]
            base = rounds[0]
            per_seed.append({
                "seed": sd, "rounds": rounds,
                "gain_full_test": base["mae_full_test"] - last["mae_full_test"],
                "gain_nonleaked_test": base["mae_nonleaked_test"] - last["mae_nonleaked_test"],
                "gain_scaffold_holdout": base["mae_scaffold_holdout"] - last["mae_scaffold_holdout"],
            })
        gf = np.array([r["gain_full_test"] for r in per_seed])
        gn = np.array([r["gain_nonleaked_test"] for r in per_seed])
        gh = np.array([r["gain_scaffold_holdout"] for r in per_seed])
        results["arms"][arm_name] = {
            "curator": arm_id,
            "per_seed": per_seed,
            "mean_gain_full_test": float(gf.mean()),
            "mean_gain_nonleaked_test": float(gn.mean()),
            "mean_gain_scaffold_holdout": float(gh.mean()),
            "mean_overstatement_full_minus_nonleaked": float((gf - gn).mean()),
            "per_seed_overstatement": (gf - gn).tolist(),
        }
    calls_fh.close()

    # ---- pooled paired test over (arm x seed), LLM arms only --------------
    d = []
    for a, v in results["arms"].items():
        if v["curator"] == "random":
            continue
        d += [s["gain_full_test"] - s["gain_nonleaked_test"] for s in v["per_seed"]]
    if d:
        d = np.array(d)
        rng = np.random.default_rng(zlib.crc32(b"s2|overstatement"))
        boots = d[rng.integers(0, len(d), size=(10000, len(d)))].mean(axis=1)
        results["pooled_llm_overstatement"] = {
            "n_units": len(d), "unit": "curator x seed",
            "mean": float(d.mean()),
            "ci": [float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))],
            "ci_excludes_zero": bool(np.quantile(boots, 0.025) > 0
                                     or np.quantile(boots, 0.975) < 0),
        }

    Path(args.out).write_text(json.dumps(results, indent=2) + "\n")
    print(f"wrote {args.out}  ({time.time()-t_start:.0f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
