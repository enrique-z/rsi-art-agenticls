"""bench_data.py — loaders for the three TDC ADMET benchmarks (REAL files only).

Two views of the same three datasets are used and must not be confused:

* the flat `*.tab` release files (`tdc_bbb_martins`, `tdc_lipophilicity`,
  `tdc_solubility`) — the FULL dataset, quoted fields, no split; and
* the TDC **ADMET Benchmark Group** scaffold splits (`tdc_admet_group`) —
  `train_val.csv` / `test.csv` per task, the split the public leaderboard uses.

S1's "held-out test molecules" are the benchmark-group `test.csv` rows.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from data_access import is_trap, resolve  # noqa: E402

BENCHMARKS = {
    "bbb_martins": {
        "tab_key": "tdc_bbb_martins",
        "task": "classification",
        "metric": "roc_auc",
        "label": "BBB penetration (binary)",
        "display": "BBB-Martins",
    },
    "lipophilicity_astrazeneca": {
        "tab_key": "tdc_lipophilicity",
        "task": "regression",
        "metric": "mae",
        "label": "logD7.4",
        "display": "Lipophilicity-AZ",
    },
    "solubility_aqsoldb": {
        "tab_key": "tdc_solubility",
        "task": "regression",
        "metric": "mae",
        "label": "log S (mol/L)",
        "display": "Solubility-AqSolDB",
    },
}


@dataclass
class Split:
    benchmark: str
    split: str
    drug_id: list[str]
    smiles: list[str]
    y: list[float]

    def __len__(self) -> int:
        return len(self.smiles)


def load_group_split(benchmark: str, split: str) -> Split:
    """Load the canonical TDC ADMET Benchmark Group split (`train_val` | `test`)."""
    assert benchmark in BENCHMARKS, benchmark
    assert split in ("train_val", "test"), split
    root = resolve("tdc_admet_group")
    trap = is_trap(root)
    assert not trap, f"TRAP dataset: {trap}"
    path = Path(root) / benchmark / f"{split}.csv"
    assert path.exists(), f"missing split file {path}"
    ids, smi, ys = [], [], []
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh)
        assert rd.fieldnames == ["Drug_ID", "Drug", "Y"], f"unexpected header {rd.fieldnames}"
        for row in rd:
            ids.append(row["Drug_ID"])
            smi.append(row["Drug"])
            ys.append(float(row["Y"]))
    assert len(smi) > 0, f"empty split {path} — FAIL"
    return Split(benchmark, split, ids, smi, ys)


def load_tab(benchmark: str) -> Split:
    """Load the flat `*.tab` full release file (quoted fields — quotes stripped)."""
    key = BENCHMARKS[benchmark]["tab_key"]
    path = resolve(key)
    trap = is_trap(path)
    assert not trap, f"TRAP dataset: {trap}"
    ids, smi, ys = [], [], []
    with open(path, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t", quotechar='"')
        header = next(rd)
        assert header == ["Drug_ID", "Drug", "Y"], f"unexpected header {header}"
        for row in rd:
            if not row:
                continue
            ids.append(row[0])
            smi.append(row[1])
            ys.append(float(row[2]))
    assert len(smi) > 0, f"empty file {path} — FAIL"
    return Split(benchmark, "full_tab", ids, smi, ys)
