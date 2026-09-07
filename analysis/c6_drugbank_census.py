#!/usr/bin/env python3
"""c6_drugbank_census.py — C6, the second-corpus census.

The manuscript's own limitation says ChEMBL is one of several corpora a curator
could query, that DrugBank would need its own census, and that our exposure
rates are therefore a lower bound.  This turns that hedge into a measurement.

What it does
    1. Stream-parses DrugBank 5.1 `full database.xml` (1.5 GB) with
       `iterparse`, pulling the calculated SMILES for every drug record.
    2. Canonicalises and fingerprints with the SAME toolkit, radius and bit
       length as the ChEMBL census (`leakage_lib.morgan_packed`), so the two
       corpora are measured on one scale.
    3. Runs the SAME exact nearest-neighbour Tanimoto search
       (`leakage_lib.exact_nn_tanimoto`) from each of the three TDC ADMET test
       splits against DrugBank -- no index, no sampling.
    4. Reports, per benchmark: the DrugBank-only exposure rate at >= 0.9 with a
       Wilson interval, and the UNION rate against ChEMBL 36 OR DrugBank, taken
       as the per-molecule maximum over both corpora.
    5. Reports the quantity that actually matters for the paper's method: how
       much the leakage-clean subset SHRINKS when one more public corpus is
       added.  If a leakage-clean evaluation subset is materially smaller once a
       second corpus is counted, then "leakage-clean" is corpus-relative and the
       paper's own rates are demonstrably a lower bound rather than arguably so.

Writes evidence/c6_drugbank_census.json.  Modifies nothing else.
"""

from __future__ import annotations

import json
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bench_data import BENCHMARKS, load_group_split  # noqa: E402
from leakage_lib import (  # noqa: E402
    FP_BYTES,
    exact_nn_tanimoto,
    morgan_packed,
    sha256_file,
    wilson_ci,
)

PAPER = Path(__file__).resolve().parents[1]
EV = PAPER / "evidence"
DRUGBANK = Path("<DATA_ROOT>/drugbank/5.1/full database.xml")
LEAK_T = 0.9
NS = "{http://www.drugbank.ca}"


def extract_drugbank_smiles(path: Path) -> tuple[list[str], dict]:
    """Stream the XML and pull the calculated SMILES of every <drug> record."""
    smiles, n_drugs, n_with_smiles = [], 0, 0
    ctx = ET.iterparse(str(path), events=("end",))
    for _, el in ctx:
        if el.tag != f"{NS}drug":
            continue
        n_drugs += 1
        s = None
        cp = el.find(f"{NS}calculated-properties")
        if cp is not None:
            for prop in cp.findall(f"{NS}property"):
                kind = prop.findtext(f"{NS}kind")
                if kind == "SMILES":
                    s = prop.findtext(f"{NS}value")
                    break
        if s:
            n_with_smiles += 1
            smiles.append(s.strip())
        el.clear()               # keep memory flat over a 1.5 GB file
        if n_drugs % 2000 == 0:
            print(f"  parsed {n_drugs} drug records, {n_with_smiles} with SMILES",
                  flush=True)
    return smiles, {"n_drug_records": n_drugs, "n_with_calculated_smiles": n_with_smiles}


def fp_matrix(smiles_list):
    packed, cans, n_fail = [], [], 0
    for s in smiles_list:
        r = morgan_packed(s)
        if r is None:
            n_fail += 1
            continue
        cans.append(r[0])
        packed.append(r[1])
    arr = np.frombuffer(b"".join(packed), dtype=np.uint8).reshape(len(packed), FP_BYTES)
    return arr, cans, n_fail


def main() -> int:
    assert DRUGBANK.exists(), f"missing {DRUGBANK}"
    t0 = time.time()
    print(f"hashing {DRUGBANK.name} …", flush=True)
    db_sha = sha256_file(DRUGBANK)
    print(f"  sha256 {db_sha}", flush=True)

    print("parsing DrugBank XML …", flush=True)
    raw, parse_stats = extract_drugbank_smiles(DRUGBANK)
    print(f"  {len(raw)} SMILES in {time.time()-t0:.0f}s", flush=True)

    db_fp, db_can, db_fail = fp_matrix(raw)
    uniq = sorted(set(db_can))
    print(f"  {len(db_can)} RDKit-parseable ({db_fail} failures), "
          f"{len(uniq)} distinct canonical structures", flush=True)
    # deduplicate to distinct canonical structures, as the ChEMBL census does not
    # (ChEMBL rows are already distinct structures) -- report both counts
    keep = {}
    for i, c in enumerate(db_can):
        keep.setdefault(c, i)
    idx = np.array(sorted(keep.values()))
    db_fp_u = db_fp[idx]

    out = {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "corpus": "DrugBank 5.1",
        "corpus_path": str(DRUGBANK),
        "corpus_bytes": DRUGBANK.stat().st_size,
        "corpus_sha256": db_sha,
        "parse": {**parse_stats, "n_rdkit_parse_failures": db_fail,
                  "n_distinct_canonical": len(uniq)},
        "leak_threshold": LEAK_T,
        "fingerprint": "Morgan r=2, 2048 bits (identical to the ChEMBL census)",
        "benchmarks": {},
    }

    nn_rows = [json.loads(l) for l in (EV / "s1_nn_per_molecule.jsonl").read_text().splitlines()]
    for b, cfg in BENCHMARKS.items():
        test = load_group_split(b, "test")
        te_fp, te_can, _ = fp_matrix(test.smiles)
        by_can = {r["cansmi"]: r for r in nn_rows
                  if r["benchmark"] == b and r["split"] == "test"}
        chembl_sim = np.array([by_can[c]["nn_chembl_tanimoto"] for c in te_can], float)

        r = exact_nn_tanimoto(te_fp, db_fp_u, chunk=8192)
        db_sim = r.max_sim
        n = len(te_can)

        db_leak = db_sim >= LEAK_T
        ch_leak = chembl_sim >= LEAK_T
        union_leak = np.maximum(db_sim, chembl_sim) >= LEAK_T

        k_db, k_ch, k_un = int(db_leak.sum()), int(ch_leak.sum()), int(union_leak.sum())
        lo_db, hi_db = wilson_ci(k_db, n)
        lo_un, hi_un = wilson_ci(k_un, n)

        clean_ch = int((~ch_leak).sum())
        clean_un = int((~union_leak).sum())
        newly_lost = int((~ch_leak & db_leak).sum())

        out["benchmarks"][b] = {
            "display": cfg["display"], "n_test": n,
            "drugbank_rate_0p9": round(100.0 * k_db / n, 2),
            "drugbank_ci": [round(100.0 * lo_db, 2), round(100.0 * hi_db, 2)],
            "n_leaked_drugbank": k_db,
            "chembl_rate_0p9": round(100.0 * k_ch / n, 2),
            "union_rate_0p9": round(100.0 * k_un / n, 2),
            "union_ci": [round(100.0 * lo_un, 2), round(100.0 * hi_un, 2)],
            "n_leaked_union": k_un,
            "clean_subset_chembl_only": clean_ch,
            "clean_subset_union": clean_un,
            "clean_molecules_lost_to_second_corpus": newly_lost,
            "clean_shrinkage_percent": (round(100.0 * newly_lost / clean_ch, 2)
                                        if clean_ch else None),
            "median_nn_drugbank": round(float(np.median(db_sim)), 4),
        }
        e = out["benchmarks"][b]
        print(f"[C6] {b}: DrugBank {e['drugbank_rate_0p9']}% | ChEMBL "
              f"{e['chembl_rate_0p9']}% | union {e['union_rate_0p9']}% | clean "
              f"{clean_ch} -> {clean_un} (lost {newly_lost}, "
              f"{e['clean_shrinkage_percent']}%)", flush=True)

    p = EV / "c6_drugbank_census.json"
    p.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {p} in {time.time()-t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
