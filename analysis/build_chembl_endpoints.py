#!/usr/bin/env python3
"""build_chembl_endpoints.py — extract REAL ChEMBL 36 measurements for the two
TDC ADMET endpoints that ChEMBL actually carries.

This is the training pool a ChEMBL-drawing LLM curator would have available.

    lipophilicity  <- activities.standard_type = 'LogD', relation '=', assay
                      description mentioning pH 7.4, value in [-5, 10] (log units)
    solubility     <- activities.standard_type = 'Solubility', relation '=',
                      units nM or ug.mL-1, converted to log10(mol/L) using
                      compound_properties.full_mwt; value in [-12, 2]

BBB penetration is deliberately ABSENT: ChEMBL 36 carries only ~164 activity rows
of logBB / Log BB type, far too few to train on.  That is recorded as a
pre-declaration deviation in DESIGN_NOTES.md, not silently worked around.

Writes <work>/chembl_endpoint_<name>.jsonl and <work>/chembl_endpoints.json.
"""

from __future__ import annotations

import argparse
import json
import math
import sqlite3
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from data_access import is_trap, resolve  # noqa: E402

LOGD_SQL = """
SELECT act.molregno, act.standard_value, cs.canonical_smiles
FROM activities act
JOIN assays a  ON act.assay_id = a.assay_id
JOIN compound_structures cs ON cs.molregno = act.molregno
WHERE act.standard_type = 'LogD'
  AND act.standard_relation = '='
  AND act.standard_value IS NOT NULL
  AND act.standard_value BETWEEN -5 AND 10
  AND a.description LIKE '%7.4%'
  AND cs.canonical_smiles IS NOT NULL
"""

SOL_SQL = """
SELECT act.molregno, act.standard_value, act.standard_units,
       cp.full_mwt, cs.canonical_smiles
FROM activities act
JOIN compound_structures cs ON cs.molregno = act.molregno
JOIN compound_properties cp ON cp.molregno = act.molregno
WHERE act.standard_type = 'Solubility'
  AND act.standard_relation = '='
  AND act.standard_value IS NOT NULL
  AND act.standard_value > 0
  AND act.standard_units IN ('nM', 'ug.mL-1')
  AND cp.full_mwt IS NOT NULL
  AND cp.full_mwt > 0
  AND cs.canonical_smiles IS NOT NULL
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="<REPO_ROOT>/.work/agenticls/full")
    args = ap.parse_args()
    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)

    db = resolve("chembl_36")
    trap = is_trap(db)
    assert not trap, f"TRAP dataset: {trap}"
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)

    meta: dict = {"db_path": str(db), "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}

    # ---- lipophilicity (logD 7.4) -----------------------------------------
    vals: dict[int, list[float]] = defaultdict(list)
    smi: dict[int, str] = {}
    rows = conn.execute(LOGD_SQL).fetchall()
    assert len(rows) > 0, "empty LogD query — FAIL"
    for regno, v, s in rows:
        vals[regno].append(float(v))
        smi[regno] = s
    _write(work, "lipophilicity", vals, smi)
    meta["lipophilicity"] = {
        "sql": " ".join(LOGD_SQL.split()),
        "n_activity_rows": len(rows),
        "n_compounds": len(vals),
        "aggregation": "median of standard_value per molregno",
        "units": "log D at pH 7.4 (dimensionless)",
    }
    print(f"lipophilicity: {len(rows):,} rows -> {len(vals):,} compounds", flush=True)

    # ---- aqueous solubility (log S, mol/L) --------------------------------
    vals = defaultdict(list)
    smi = {}
    n_raw = 0
    n_dropped_range = 0
    for regno, v, units, mwt, s in conn.execute(SOL_SQL):
        n_raw += 1
        v = float(v)
        mwt = float(mwt)
        if units == "nM":
            mol_per_l = v * 1e-9
        else:  # ug.mL-1  ->  mg/L -> g/L = v*1e-3 ; mol/L = (v*1e-3)/MW
            mol_per_l = (v * 1e-3) / mwt
        if mol_per_l <= 0:
            n_dropped_range += 1
            continue
        log_s = math.log10(mol_per_l)
        if not (-12.0 <= log_s <= 2.0):
            n_dropped_range += 1
            continue
        vals[regno].append(log_s)
        smi[regno] = s
    assert n_raw > 0, "empty Solubility query — FAIL"
    _write(work, "solubility", vals, smi)
    meta["solubility"] = {
        "sql": " ".join(SOL_SQL.split()),
        "n_activity_rows_raw": n_raw,
        "n_dropped_out_of_range": n_dropped_range,
        "n_compounds": len(vals),
        "aggregation": "median of log10(mol/L) per molregno",
        "unit_conversion": "nM -> mol/L = v*1e-9 ; ug.mL-1 -> mol/L = (v*1e-3)/full_mwt",
        "range_filter": "log10(mol/L) in [-12, 2]",
        "units": "log S (mol/L)",
    }
    print(f"solubility: {n_raw:,} rows ({n_dropped_range:,} out of range) -> {len(vals):,} compounds",
          flush=True)

    # ---- BBB: record the ABSENCE, with the count that justifies it --------
    n_bbb = conn.execute(
        "SELECT count(*) FROM activities WHERE standard_type IN "
        "('logBB','Log BB','pBBB','Brain/Blood concentration','Kbb')"
    ).fetchone()[0]
    meta["bbb_martins"] = {
        "available": False,
        "n_candidate_activity_rows": int(n_bbb),
        "reason": "ChEMBL 36 carries only a few hundred blood-brain-barrier "
                  "distribution activity rows across incompatible standard_types; "
                  "no ChEMBL-derived BBB training pool of usable size exists. "
                  "S1b arm A is therefore NOT RUN for BBB-Martins; arm B "
                  "(train on the benchmark's own train_val split) still runs.",
    }
    print(f"bbb: {n_bbb} candidate rows -> arm A not runnable", flush=True)

    (work / "chembl_endpoints.json").write_text(json.dumps(meta, indent=2) + "\n")
    conn.close()
    return 0


def _write(work: Path, name: str, vals: dict, smi: dict) -> None:
    import statistics

    path = work / f"chembl_endpoint_{name}.jsonl"
    with open(path, "w") as fh:
        for regno in sorted(vals):
            fh.write(json.dumps({
                "molregno": int(regno),
                "y": float(statistics.median(vals[regno])),
                "n_measurements": len(vals[regno]),
                "smiles": smi[regno],
            }) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
