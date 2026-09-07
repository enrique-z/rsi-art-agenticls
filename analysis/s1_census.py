#!/usr/bin/env python3
"""s1_census.py — Stage 1: near-duplicate leakage census, TDC ADMET vs ChEMBL 36.

For every molecule of the canonical TDC ADMET Benchmark Group splits (test and
train_val) of BBB-Martins, Lipophilicity-AstraZeneca and Solubility-AqSolDB, this
script computes the EXACT nearest-neighbour Tanimoto similarity (Morgan r=2,
2048 bits) against the whole ChEMBL 36 compound corpus (2,854,800 parsed
structures), plus the same quantity against the benchmark's own train_val split
as a reference scale.

Outputs (paths relative to the paper root):
    evidence/s1_nn_per_molecule.jsonl     one row per benchmark molecule
    evidence/s1_census_aggregates.json    per-benchmark rates + Wilson CIs
    evidence/s1_exactness_audit.json      NN exactness audit

Pilot mode (`--pilot`) runs one benchmark x 100 test molecules against a small
ChEMBL subset and must finish inside ~60 s before the full census is launched.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bench_data import BENCHMARKS, load_group_split  # noqa: E402
from leakage_lib import (  # noqa: E402
    FP_BYTES,
    audit_exactness,
    bemis_murcko,
    exact_nn_tanimoto,
    morgan_packed,
    wilson_ci,
)

os.environ.setdefault("PYTHONHASHSEED", "0")

PAPER = Path(__file__).resolve().parents[1]
THRESHOLDS = [0.7, 0.8, 0.9, 0.95, 0.99, 1.0]
HIST_EDGES = np.round(np.arange(0.0, 1.0001, 0.05), 4)


def fingerprint_split(split) -> dict:
    """Canonicalise + fingerprint one split.  Unparseable molecules are counted."""
    packed, cans, keep = [], [], []
    n_fail = 0
    for i, smi in enumerate(split.smiles):
        r = morgan_packed(smi)
        if r is None:
            n_fail += 1
            continue
        cans.append(r[0])
        packed.append(r[1])
        keep.append(i)
    assert packed, f"no parseable molecules in {split.benchmark}/{split.split} — FAIL"
    arr = np.frombuffer(b"".join(packed), dtype=np.uint8).reshape(len(packed), FP_BYTES)
    return {
        "fps": arr,
        "cansmi": cans,
        "keep": keep,
        "n_parse_failures": n_fail,
        "split": split,
    }


def rates(sims: np.ndarray) -> dict:
    n = len(sims)
    out = {"n": int(n)}
    for t in THRESHOLDS:
        k = int((sims >= t - 1e-12).sum())
        lo, hi = wilson_ci(k, n)
        key = f"ge_{str(t).replace('.', 'p')}"
        out[key] = {"k": k, "rate": k / n, "ci_low": lo, "ci_high": hi}
    out["mean"] = float(sims.mean())
    out["median"] = float(np.median(sims))
    out["quantiles"] = {
        q: float(np.quantile(sims, q)) for q in (0.05, 0.25, 0.5, 0.75, 0.95)
    }
    hist, _ = np.histogram(sims, bins=HIST_EDGES)
    out["histogram_edges"] = [float(x) for x in HIST_EDGES]
    out["histogram_counts"] = [int(x) for x in hist]
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="<REPO_ROOT>/.work/agenticls/full")
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--audit-frac", type=float, default=0.05)
    ap.add_argument("--out-prefix", default="s1")
    args = ap.parse_args()

    work = Path(args.work)
    t_start = time.time()

    chembl_fps = np.load(work / "chembl_fps.npy", mmap_mode=None)
    chembl_regno = np.load(work / "chembl_molregno.npy")
    chembl_can = (work / "chembl_cansmi.txt").read_text().splitlines()
    assert chembl_fps.shape[0] == len(chembl_regno) == len(chembl_can), "corpus arrays disagree — FAIL"
    build_meta = json.loads((work / "chembl_build.json").read_text())
    print(f"ChEMBL corpus: {chembl_fps.shape[0]:,} fingerprints  ({time.time()-t_start:.1f}s)", flush=True)

    bench_names = list(BENCHMARKS)
    if args.pilot:
        bench_names = ["bbb_martins"]

    # ---- fingerprint every benchmark split ---------------------------------
    fp_data: dict[tuple[str, str], dict] = {}
    for b in bench_names:
        for s in ("test", "train_val"):
            sp = load_group_split(b, s)
            fp_data[(b, s)] = fingerprint_split(sp)
            if args.pilot and s == "test":
                d = fp_data[(b, s)]
                d["fps"] = d["fps"][:100]
                d["cansmi"] = d["cansmi"][:100]
                d["keep"] = d["keep"][:100]
            print(f"  {b}/{s}: {len(fp_data[(b, s)]['cansmi']):,} fp "
                  f"({fp_data[(b, s)]['n_parse_failures']} parse failures)", flush=True)

    # ---- one combined query matrix over the ChEMBL corpus -------------------
    order: list[tuple[str, str, int]] = []
    blocks = []
    for (b, s), d in fp_data.items():
        for i in range(d["fps"].shape[0]):
            order.append((b, s, i))
        blocks.append(d["fps"])
    q = np.concatenate(blocks, axis=0)
    print(f"query molecules: {q.shape[0]:,}", flush=True)

    corpus = chembl_fps[:100_000] if args.pilot else chembl_fps
    t0 = time.time()
    last = [0.0]

    def prog(done, total):
        now = time.time()
        if now - last[0] > 20:
            last[0] = now
            el = now - t0
            print(f"    NN {done:,}/{total:,}  {el:.0f}s  eta {el*(total-done)/max(done,1):.0f}s",
                  flush=True)

    res = exact_nn_tanimoto(q, corpus, chunk=4096, progress=prog)
    nn_seconds = time.time() - t0
    print(f"exact NN done in {nn_seconds:.1f}s", flush=True)

    # ---- exactness audit ----------------------------------------------------
    audit_corpus = corpus if args.pilot else corpus[:200_000]
    if args.pilot:
        audit = audit_exactness(q, corpus, res, frac=args.audit_frac)
        audit["scope"] = "full corpus used in this run"
    else:
        # Audit against the full corpus would need a per-query full popcount pass;
        # do exactly that on the pre-declared 5% subsample.
        audit = audit_exactness(q, corpus, res, frac=args.audit_frac)
        audit["scope"] = "full ChEMBL corpus, 5% random subsample of queries"
    print("audit:", json.dumps(audit), flush=True)
    assert audit["max_abs_deviation_vs_popcount"] == 0.0, "NN not exact vs popcount — FAIL"
    assert audit["max_abs_deviation_vs_rdkit"] == 0.0, "NN not exact vs RDKit — FAIL"

    # ---- own-train-split reference NN ---------------------------------------
    own: dict[str, np.ndarray] = {}
    for b in bench_names:
        te = fp_data[(b, "test")]["fps"]
        tr = fp_data[(b, "train_val")]["fps"]
        own[b] = exact_nn_tanimoto(te, tr, chunk=8192).max_sim

    # ---- per-molecule rows ---------------------------------------------------
    chembl_smiles_set = set(chembl_can)
    ev = PAPER / "evidence"
    ev.mkdir(exist_ok=True)
    per_mol_path = ev / f"{args.out_prefix}_nn_per_molecule.jsonl"

    counters: dict[tuple[str, str], list] = {k: [] for k in fp_data}
    own_pos = {b: 0 for b in bench_names}
    with open(per_mol_path, "w") as fh:
        for row_i, (b, s, i) in enumerate(order):
            d = fp_data[(b, s)]
            sp = d["split"]
            src_i = d["keep"][i]
            can = d["cansmi"][i]
            sim = float(res.max_sim[row_i])
            j = int(res.argmax[row_i])
            rec = {
                "benchmark": b,
                "split": s,
                "drug_id": sp.drug_id[src_i],
                "smiles": sp.smiles[src_i],
                "cansmi": can,
                "y": sp.y[src_i],
                "nn_chembl_tanimoto": sim,
                "nn_chembl_molregno": int(chembl_regno[j]) if j >= 0 else None,
                "nn_chembl_cansmi": chembl_can[j] if j >= 0 else None,
                "exact_smiles_dup_in_chembl": can in chembl_smiles_set,
                "scaffold": bemis_murcko(can),
            }
            if s == "test":
                rec["nn_own_trainval_tanimoto"] = float(own[b][own_pos[b]])
                own_pos[b] += 1
            counters[(b, s)].append(rec)
            fh.write(json.dumps(rec) + "\n")

    # ---- aggregates ----------------------------------------------------------
    agg = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "pilot": bool(args.pilot),
        "chembl_corpus": {
            "n_fingerprints": int(corpus.shape[0]),
            "build": build_meta,
        },
        "nn_seconds": round(nn_seconds, 1),
        "exactness_audit": audit,
        "benchmarks": {},
    }
    for b in bench_names:
        rows_te = counters[(b, "test")]
        sims_te = np.array([r["nn_chembl_tanimoto"] for r in rows_te])
        own_te = np.array([r["nn_own_trainval_tanimoto"] for r in rows_te])
        n_exact = sum(1 for r in rows_te if r["exact_smiles_dup_in_chembl"])
        lo, hi = wilson_ci(n_exact, len(rows_te))
        scaf_leaked = {r["scaffold"] for r in rows_te if r["nn_chembl_tanimoto"] >= 0.9}
        scaf_all = {r["scaffold"] for r in rows_te}
        # scaffold-clustered rate: fraction of distinct test scaffolds with >=1
        # member at NN >= 0.9 (guards against one over-represented scaffold
        # carrying the whole leakage rate)
        agg["benchmarks"][b] = {
            "display": BENCHMARKS[b]["display"],
            "task": BENCHMARKS[b]["task"],
            "n_test": len(rows_te),
            "n_train_val": len(counters[(b, "train_val")]),
            "test_parse_failures": fp_data[(b, "test")]["n_parse_failures"],
            "train_val_parse_failures": fp_data[(b, "train_val")]["n_parse_failures"],
            "test_vs_chembl": rates(sims_te),
            "test_vs_own_trainval": rates(own_te),
            "trainval_vs_chembl": rates(
                np.array([r["nn_chembl_tanimoto"] for r in counters[(b, "train_val")]])
            ),
            "exact_smiles_duplicate": {
                "k": n_exact,
                "rate": n_exact / len(rows_te),
                "ci_low": lo,
                "ci_high": hi,
                "definition": "identical RDKit canonical SMILES present in ChEMBL 36 compound_structures",
            },
            "scaffold_clustering": {
                "n_distinct_test_scaffolds": len(scaf_all),
                "n_distinct_scaffolds_with_leak_ge_0p9": len(scaf_leaked),
                "scaffold_level_leak_rate": len(scaf_leaked) / max(len(scaf_all), 1),
                "largest_scaffold_share_of_test": max(
                    Counter(r["scaffold"] for r in rows_te).values()
                ) / len(rows_te),
            },
        }

    out = ev / f"{args.out_prefix}_census_aggregates.json"
    out.write_text(json.dumps(agg, indent=2) + "\n")
    (ev / f"{args.out_prefix}_exactness_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(f"\nwrote {per_mol_path}\nwrote {out}\ntotal {time.time()-t_start:.1f}s")
    for b in bench_names:
        a = agg["benchmarks"][b]
        print(f"  {a['display']:<22} n={a['n_test']:<5} "
              f">=0.9 {a['test_vs_chembl']['ge_0p9']['rate']*100:5.1f}%  "
              f">=0.95 {a['test_vs_chembl']['ge_0p95']['rate']*100:5.1f}%  "
              f"exact-dup {a['exact_smiles_duplicate']['rate']*100:5.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
