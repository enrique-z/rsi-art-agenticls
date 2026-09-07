#!/usr/bin/env python3
"""build_chembl_fps.py — canonicalise + fingerprint the ChEMBL 36 compound corpus.

Reads `compound_structures.canonical_smiles` from the real ChEMBL 36 SQLite
database (resolved through `dataset_paths.resolve`, never hardcoded), parses every
SMILES with RDKit, and writes:

    <work>/chembl_fps.npy        (N, 256) uint8   packed Morgan r=2 / 2048-bit
    <work>/chembl_molregno.npy   (N,)     int64   ChEMBL molregno per row
    <work>/chembl_cansmi.txt     N lines          RDKit canonical SMILES per row
    <work>/chembl_build.json                      provenance + counts

Usage:
    python3 build_chembl_fps.py --limit 100000 --out <dir>     # pilot
    python3 build_chembl_fps.py --out <dir>                    # full corpus
"""

from __future__ import annotations

import argparse
import json
import os
import sqlite3
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from data_access import is_trap, resolve  # noqa: E402
from leakage_lib import FP_BITS, FP_BYTES, FP_RADIUS, morgan_batch  # noqa: E402

os.environ.setdefault("PYTHONHASHSEED", "0")


def rows(db_path: Path, limit: int | None):
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    q = ("SELECT molregno, canonical_smiles FROM compound_structures "
         "WHERE canonical_smiles IS NOT NULL ORDER BY molregno")
    if limit:
        q += f" LIMIT {limit}"
    cur = conn.execute(q)
    while True:
        batch = cur.fetchmany(5000)
        if not batch:
            break
        yield batch
    conn.close()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--procs", type=int, default=10)
    args = ap.parse_args()

    db_path = resolve("chembl_36")
    trap = is_trap(db_path)
    assert not trap, f"TRAP dataset: {trap}"
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    n_total = conn.execute("SELECT count(*) FROM compound_structures").fetchone()[0]
    n_nonnull = conn.execute(
        "SELECT count(*) FROM compound_structures WHERE canonical_smiles IS NOT NULL"
    ).fetchone()[0]
    conn.close()
    assert n_total > 0 and n_nonnull > 0, "empty ChEMBL query — FAIL"

    fps: list[bytes] = []
    regnos: list[int] = []
    cans: list[str] = []
    n_read = 0
    n_failed = 0
    t0 = time.time()

    with Pool(args.procs) as pool:
        pending = []
        for batch in rows(db_path, args.limit):
            n_read += len(batch)
            pending.append(batch)
            if len(pending) >= args.procs:
                results = pool.map(morgan_batch, [[r[1] for r in b] for b in pending])
                for b, res in zip(pending, results):
                    n_failed += len(b) - len(res)
                    for i, can, fp in res:
                        regnos.append(b[i][0])
                        cans.append(can)
                        fps.append(fp)
                pending = []
                el = time.time() - t0
                print(f"  {n_read:,} read  {len(fps):,} fp  {el:.0f}s "
                      f"({n_read/max(el,1e-9):,.0f}/s)", flush=True)
        if pending:
            results = pool.map(morgan_batch, [[r[1] for r in b] for b in pending])
            for b, res in zip(pending, results):
                n_failed += len(b) - len(res)
                for i, can, fp in res:
                    regnos.append(b[i][0])
                    cans.append(can)
                    fps.append(fp)

    assert len(fps) > 0, "no fingerprints produced — FAIL"
    arr = np.frombuffer(b"".join(fps), dtype=np.uint8).reshape(len(fps), FP_BYTES)
    np.save(out / "chembl_fps.npy", arr)
    np.save(out / "chembl_molregno.npy", np.asarray(regnos, dtype=np.int64))
    (out / "chembl_cansmi.txt").write_text("\n".join(cans) + "\n")

    meta = {
        "db_path": str(db_path),
        "table": "compound_structures",
        "n_rows_table": int(n_total),
        "n_rows_smiles_nonnull": int(n_nonnull),
        "limit": args.limit,
        "n_read": int(n_read),
        "n_fingerprinted": int(len(fps)),
        "n_rdkit_parse_failures": int(n_failed),
        "fp_radius": FP_RADIUS,
        "fp_bits": FP_BITS,
        "fp_generator": "rdkit.Chem.rdFingerprintGenerator.GetMorganGenerator",
        "elapsed_s": round(time.time() - t0, 1),
    }
    (out / "chembl_build.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
