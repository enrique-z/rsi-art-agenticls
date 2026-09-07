#!/usr/bin/env python3
"""build_chembl_scaffolds.py — Bemis-Murcko scaffold frequency over ChEMBL 36.

Provides the `scaffold frequency` covariate the pre-declaration names for the
subset-composition rival mechanism: how common is a test molecule's Bemis-Murcko
scaffold in the public corpus a curator draws from?

Reads <work>/chembl_cansmi.txt (produced by build_chembl_fps.py) and writes
<work>/chembl_scaffold_counts.json : {scaffold_smiles: count}.
"""

from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from multiprocessing import Pool
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from leakage_lib import bemis_murcko  # noqa: E402


def batch(smis: list[str]) -> list[str]:
    return [bemis_murcko(s) or "" for s in smis]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="<REPO_ROOT>/.work/agenticls/full")
    ap.add_argument("--procs", type=int, default=10)
    args = ap.parse_args()
    work = Path(args.work)
    smis = (work / "chembl_cansmi.txt").read_text().splitlines()
    assert smis, "empty ChEMBL canonical SMILES file — FAIL"
    t0 = time.time()
    counts: Counter[str] = Counter()
    chunks = [smis[i:i + 5000] for i in range(0, len(smis), 5000)]
    with Pool(args.procs) as pool:
        for k, res in enumerate(pool.imap(batch, chunks, chunksize=1)):
            counts.update(res)
            if k % 50 == 0:
                print(f"  {min((k+1)*5000, len(smis)):,}/{len(smis):,} "
                      f"{time.time()-t0:.0f}s", flush=True)
    n_unparsed = counts.pop("", 0)
    out = work / "chembl_scaffold_counts.json"
    out.write_text(json.dumps({
        "n_molecules": len(smis),
        "n_scaffold_failures": int(n_unparsed),
        "n_distinct_scaffolds": len(counts),
        "counts": dict(counts),
    }))
    print(f"wrote {out}: {len(counts):,} distinct scaffolds, {time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
