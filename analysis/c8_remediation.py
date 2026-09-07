#!/usr/bin/env python3
"""c8_remediation.py — C8, does the R6 remedy actually work?

Design, primary quantity and the reporting rule for every outcome -- including
the one that undercuts our own design rule R6 -- were fixed in
`C8_PRECOMMITMENT_2026-09-07.md` before this file was written; its SHA-256 is
recorded in the output.

Re-scores the ALREADY-RUN S2 loop on progressively stricter decontaminations of
the evaluation split.  No new model calls: every curator selection is already in
`evidence/s2_rsi_results.json`, so each arm x seed x round is refitted on the
pool implied by its recorded selection, with the same `zlib.crc32` seed key, and
is therefore byte-identical to the original pool.

Scored on five sets:
    full            the full TDC test split
    clean@0.9 chembl   the paper's current leakage-clean subset
    clean@0.9 union    decontaminated against ChEMBL 36 AND DrugBank 5.1
    clean@0.8 union    progressively more conservative
    clean@0.7 union

Reproduction check (outcome D in the pre-commitment): the recomputed hold-out
gains must match the published S2 gains to printing precision, because the
hold-out does not depend on the test split at all.  A mismatch voids the run.

Writes evidence/c8_remediation.json.  Modifies nothing else.
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
from leakage_lib import (  # noqa: E402
    FP_BYTES,
    bemis_murcko,
    exact_nn_tanimoto,
    morgan_packed,
    unpack_f32,
)

PAPER = Path(__file__).resolve().parents[1]
EV = PAPER / "evidence"
WORK = Path("<REPO_ROOT>/.work/agenticls/full")
DRUGBANK_JSON = EV / "c6_drugbank_census.json"
DRUGBANK_XML = Path("<DATA_ROOT>/drugbank/5.1/full database.xml")
BENCH = "solubility_aqsoldb"
POOL_FILE = "chembl_endpoint_solubility.jsonl"
CLUSTERS, BUDGET, TREES = 40, 6000, 200
MIN_SUBSET_N = 20            # fixed in code before the S1b run; unchanged
REPRO_TOL = 5e-4             # hold-out gains must match to printing precision
PRECOMMIT = PAPER / "C8_PRECOMMITMENT_2026-09-07.md"
NS = "{http://www.drugbank.ca}"


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


def drugbank_fps():
    """Re-extract DrugBank SMILES exactly as c6_drugbank_census.py does."""
    import xml.etree.ElementTree as ET
    smiles = []
    for _, el in ET.iterparse(str(DRUGBANK_XML), events=("end",)):
        if el.tag != f"{NS}drug":
            continue
        cp = el.find(f"{NS}calculated-properties")
        if cp is not None:
            for prop in cp.findall(f"{NS}property"):
                if prop.findtext(f"{NS}kind") == "SMILES":
                    v = prop.findtext(f"{NS}value")
                    if v:
                        smiles.append(v.strip())
                    break
        el.clear()
    arr, _, cans = fp_matrix(smiles)
    keep = {}
    for i, c in enumerate(cans):
        keep.setdefault(c, i)
    idx = np.array(sorted(keep.values()))
    return arr[idx], len(keep)


def boot(d, key):
    d = np.array(d, float)
    rng = np.random.default_rng(zlib.crc32(key.encode()))
    m = np.array([rng.choice(d, len(d), replace=True).mean() for _ in range(10000)])
    return float(d.mean()), [float(np.quantile(m, .025)), float(np.quantile(m, .975))]


def main() -> int:
    assert PRECOMMIT.exists(), "pre-commitment missing — refusing to run"
    precommit_sha = hashlib.sha256(PRECOMMIT.read_bytes()).hexdigest()
    print(f"pre-commitment SHA-256 {precommit_sha}", flush=True)

    real = json.loads((EV / "s2_rsi_results.json").read_text())

    # ---- rebuild the S2 universe exactly ----------------------------------
    pool = [json.loads(l) for l in (WORK / POOL_FILE).read_text().splitlines()]
    p_fp, p_keep, p_can = fp_matrix([x["smiles"] for x in pool])
    p_y = np.array([pool[i]["y"] for i in p_keep], float)

    test = load_group_split(BENCH, "test")
    te_fp, te_keep, te_can = fp_matrix(test.smiles)
    te_y = np.array([test.y[i] for i in te_keep], float)
    nn_rows = {r["cansmi"]: r for l in (EV / "s1_nn_per_molecule.jsonl").read_text().splitlines()
               for r in [json.loads(l)] if r["benchmark"] == BENCH and r["split"] == "test"}
    chembl_sim = np.array([nn_rows[c]["nn_chembl_tanimoto"] for c in te_can], float)
    te_can_set = set(te_can)

    scaffolds = [bemis_murcko(c) or "" for c in p_can]
    held = {s for s in sorted(set(scaffolds))
            if int(hashlib.sha256(("agenticls-rsi-leakage|" + s).encode())
                   .hexdigest()[:8], 16) % 5 == 0}
    is_held = np.array([s in held for s in scaffolds])
    is_dup = np.array([c in te_can_set for c in p_can])
    holdout_idx = np.flatnonzero(is_held & ~is_dup)
    universe_idx = np.flatnonzero(~is_held & ~is_dup)
    commit = json.loads((EV / "s2_scaffold_holdout_commit.json").read_text())
    ho_sha = hashlib.sha256(
        "\n".join(sorted(p_can[i] for i in holdout_idx)).encode()).hexdigest()
    assert ho_sha == commit["holdout_sha256"], "hold-out does not reproduce — FAIL"
    print(f"hold-out SHA-256 re-asserted: {ho_sha[:16]}…", flush=True)

    U = unpack_f32(p_fp[universe_idx])
    km = MiniBatchKMeans(n_clusters=CLUSTERS, random_state=0, n_init=5,
                         batch_size=2048).fit(U)
    labels = km.labels_
    u_y = p_y[universe_idx]
    te_X = unpack_f32(te_fp)
    ho_X = unpack_f32(p_fp[holdout_idx])
    ho_y = p_y[holdout_idx]

    # ---- DrugBank neighbours for the test split ---------------------------
    print("re-extracting DrugBank …", flush=True)
    db_fp, db_n = drugbank_fps()
    db_sim = exact_nn_tanimoto(te_fp, db_fp, chunk=8192).max_sim
    print(f"  {db_n} distinct DrugBank structures", flush=True)
    union_sim = np.maximum(chembl_sim, db_sim)

    subsets = {
        "clean_0p9_chembl": chembl_sim < 0.9,
        "clean_0p9_union": union_sim < 0.9,
        "clean_0p8_union": union_sim < 0.8,
        "clean_0p7_union": union_sim < 0.7,
    }
    sizes = {k: int(v.sum()) for k, v in subsets.items()}
    print("subset sizes:", sizes, flush=True)

    def draw(cids, seed_key):
        rng = np.random.default_rng(zlib.crc32(seed_key.encode()))
        cand = (np.flatnonzero(np.isin(labels, cids)) if cids is not None
                else np.arange(len(labels)))
        return rng.choice(cand, min(BUDGET, len(cand)), replace=False)

    def score(sel, seed):
        m = RandomForestRegressor(n_estimators=TREES, n_jobs=-1,
                                  random_state=seed, max_features="sqrt")
        m.fit(U[sel], u_y[sel])
        p_te, p_ho = m.predict(te_X), m.predict(ho_X)
        out = {"full": float(mean_absolute_error(te_y, p_te)),
               "holdout": float(mean_absolute_error(ho_y, p_ho))}
        for k, mask in subsets.items():
            out[k] = (float(mean_absolute_error(te_y[mask], p_te[mask]))
                      if mask.sum() >= MIN_SUBSET_N else None)
        return out

    per_unit, repro_err = [], []
    for arm, v in real["arms"].items():
        if v["curator"] == "random":
            continue                      # LLM arms only, as pre-committed
        for s in v["per_seed"]:
            sd = s["seed"]
            r0 = score(draw(None, f"{arm}|{sd}|round0"), sd)
            last = s["rounds"][-1]
            r3 = score(draw(last["selected_clusters"], f"{arm}|{sd}|round{last['round']}"), sd)
            gains = {k: r0[k] - r3[k] for k in r0 if r0[k] is not None and r3[k] is not None}
            repro_err.append(abs(gains["holdout"] - s["gain_scaffold_holdout"]))
            per_unit.append({"arm": arm, "seed": sd, "gains": gains})
            print(f"  {arm} seed={sd} holdout {gains['holdout']:+.4f} "
                  f"(published {s['gain_scaffold_holdout']:+.4f})", flush=True)

    max_repro_err = float(max(repro_err))
    reproduced = max_repro_err < REPRO_TOL
    print(f"[C8] reproduction check: max hold-out deviation {max_repro_err:.2e} "
          f"-> {'OK' if reproduced else 'VOID (outcome D)'}", flush=True)

    out = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "precommitment_file": PRECOMMIT.name,
        "precommitment_sha256": precommit_sha,
        "benchmark": BENCH,
        "n_units": len(per_unit), "unit": "LLM curator x seed",
        "drugbank_distinct_structures": db_n,
        "subset_sizes": {"full": len(te_can), **sizes},
        "min_subset_n": MIN_SUBSET_N,
        "reproduction_check": {
            "max_holdout_gain_deviation": max_repro_err,
            "tolerance": REPRO_TOL, "reproduced": reproduced,
            "note": ("the hold-out does not depend on the test split, so a mismatch "
                     "would void the re-scoring (pre-commitment outcome D)"),
        },
        "per_unit": per_unit,
        "overstatement": {},
    }
    if not reproduced:
        out["verdict"] = "VOID — reproduction failed; no re-scored number is published"
        (EV / "c8_remediation.json").write_text(json.dumps(out, indent=2) + "\n")
        print("VOID", flush=True)
        return 1

    for k in subsets:
        if sizes[k] < MIN_SUBSET_N:
            out["overstatement"][k] = {"n_subset": sizes[k], "underpowered": True}
            continue
        d = [u["gains"]["full"] - u["gains"][k] for u in per_unit if k in u["gains"]]
        mean, ci = boot(d, f"c8|{k}")
        out["overstatement"][k] = {
            "n_subset": sizes[k], "n_units": len(d), "per_unit": d,
            "mean": mean, "ci": ci,
            "ci_excludes_zero": bool(ci[0] > 0 or ci[1] < 0),
            "underpowered": False,
        }
        print(f"[C8] overstatement vs {k} (n={sizes[k]}): {mean:+.4f} "
              f"[{ci[0]:+.4f}, {ci[1]:+.4f}] excl0="
              f"{out['overstatement'][k]['ci_excludes_zero']}", flush=True)

    hg = [u["gains"]["holdout"] for u in per_unit]
    out["holdout_all_negative"] = bool(all(g < 0 for g in hg))
    print(f"[C8] hold-out gains all negative: {out['holdout_all_negative']}", flush=True)

    (EV / "c8_remediation.json").write_text(json.dumps(out, indent=2) + "\n")
    print("wrote evidence/c8_remediation.json", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
