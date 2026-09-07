"""leakage_lib.py — shared primitives for the ChEMBL-vs-TDC leakage census.

Everything here operates on REAL data only.  No synthetic molecules, no random
fingerprints, no fixtures.  The only randomness in this module is the bootstrap
resampler, which resamples *observed* rows.

Design decisions that are load-bearing for the paper (see DESIGN_NOTES.md):

* Fingerprint  : RDKit Morgan, radius 2, 2048 bits, `rdFingerprintGenerator`
                 (the modern generator API; `GetMorganFingerprintAsBitVect` is
                 deprecated in RDKit 2026.03).  Counts are folded to bits.
* Canonical    : `Chem.MolToSmiles(mol)` of the parsed molecule — RDKit canonical
                 SMILES.  Molecules that fail to parse are dropped and COUNTED.
* NN search    : EXACT, not approximate.  Similarity is computed as a dense
                 float32 GEMM over the 0/1 fingerprint matrices.  Every entry of
                 the product is an integer intersection count in [0, 2048]; all
                 such integers are exactly representable in float32 (< 2**24), so
                 the GEMM is exact integer arithmetic, not an approximation.
                 `audit_exactness()` re-verifies a random subsample against
                 RDKit's own `BulkTanimotoSimilarity`.
"""

from __future__ import annotations

import hashlib
import zlib
from dataclasses import dataclass

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

FP_RADIUS = 2
FP_BITS = 2048
FP_BYTES = FP_BITS // 8

_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=FP_RADIUS, fpSize=FP_BITS)


# --------------------------------------------------------------------------- #
# fingerprints
# --------------------------------------------------------------------------- #
def morgan_packed(smiles: str) -> tuple[str, bytes] | None:
    """(canonical_smiles, 256-byte packed Morgan fp) or None if unparseable."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    can = Chem.MolToSmiles(mol)
    arr = _GEN.GetFingerprintAsNumPy(mol).astype(np.uint8)
    return can, np.packbits(arr).tobytes()


def morgan_batch(smiles_list: list[str]) -> list[tuple[int, str, bytes]]:
    """Worker entry point: [(index_in_batch, canonical_smiles, packed_fp), ...]."""
    out = []
    for i, smi in enumerate(smiles_list):
        r = morgan_packed(smi)
        if r is not None:
            out.append((i, r[0], r[1]))
    return out


def unpack_f32(packed: np.ndarray) -> np.ndarray:
    """(N, 256) uint8 -> (N, 2048) float32 0/1 matrix."""
    return np.unpackbits(packed, axis=1).astype(np.float32)


def popcounts(packed: np.ndarray) -> np.ndarray:
    """(N, 256) uint8 -> (N,) int32 bit counts."""
    return np.bitwise_count(packed).sum(axis=1).astype(np.int32)


def bemis_murcko(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception:
        return None
    return Chem.MolToSmiles(core)


# --------------------------------------------------------------------------- #
# exact nearest-neighbour search
# --------------------------------------------------------------------------- #
@dataclass
class NNResult:
    max_sim: np.ndarray      # (Q,) float64 — exact max Tanimoto
    argmax: np.ndarray       # (Q,) int64  — row index in the corpus
    n_corpus: int


def exact_nn_tanimoto(
    q_packed: np.ndarray,
    c_packed: np.ndarray,
    chunk: int = 8192,
    progress=None,
) -> NNResult:
    """Exact nearest-neighbour Tanimoto of every query against the whole corpus.

    Tanimoto(a,b) = |a & b| / (|a| + |b| - |a & b|).  The intersection matrix is
    obtained by a float32 GEMM on the 0/1 expansions; entries are integers <= 2048
    and therefore exact in float32 (all integers below 2**24 are).  The single
    division is then done in float64, so the result is bit-identical to the
    integer-popcount and RDKit reference paths (`audit_exactness` asserts this).
    """
    Q = q_packed.shape[0]
    N = c_packed.shape[0]
    qf = unpack_f32(q_packed)                       # (Q, 2048)
    qa = popcounts(q_packed).astype(np.float64)     # (Q,)

    best = np.zeros(Q, dtype=np.float64)
    best_idx = np.full(Q, -1, dtype=np.int64)
    ar = np.arange(Q)

    for start in range(0, N, chunk):
        stop = min(start + chunk, N)
        blk = c_packed[start:stop]
        cf = unpack_f32(blk)                        # (c, 2048)
        cb = popcounts(blk).astype(np.float64)      # (c,)
        inter = (qf @ cf.T).astype(np.float64)      # (Q, c) exact integers
        denom = qa[:, None] + cb[None, :] - inter
        np.maximum(denom, 1.0, out=denom)           # both fps empty -> sim 0
        sim = np.divide(inter, denom, out=inter)
        loc = np.argmax(sim, axis=1)
        val = sim[ar, loc]
        upd = val > best
        best[upd] = val[upd]
        best_idx[upd] = start + loc[upd]
        if progress is not None:
            progress(stop, N)

    return NNResult(max_sim=best, argmax=best_idx, n_corpus=N)


def bv_from_packed(row: np.ndarray):
    """Packed 256-byte fingerprint -> RDKit ExplicitBitVect (for cross-checks)."""
    bits = np.unpackbits(np.asarray(row, dtype=np.uint8))
    bv = DataStructs.ExplicitBitVect(FP_BITS)
    bv.SetBitsFromList([int(i) for i in np.nonzero(bits)[0]])
    return bv


def popcount_nn_max(q_row: np.ndarray, c_u64: np.ndarray, c_pop: np.ndarray) -> tuple[float, int]:
    """Full-corpus max Tanimoto for ONE query via integer popcount (no GEMM).

    Independent of `exact_nn_tanimoto`: uses `np.bitwise_count` on uint64 words
    instead of a float32 matrix product.  Used by `audit_exactness`.
    """
    q_u64 = np.frombuffer(np.asarray(q_row, dtype=np.uint8).tobytes(), dtype=np.uint64)
    qa = int(np.bitwise_count(np.asarray(q_row, dtype=np.uint8)).sum())
    inter = np.bitwise_count(c_u64 & q_u64).sum(axis=1).astype(np.int64)
    denom = qa + c_pop.astype(np.int64) - inter
    np.maximum(denom, 1, out=denom)
    sims = inter / denom
    j = int(np.argmax(sims))
    return float(sims[j]), j


def audit_exactness(
    q_packed: np.ndarray,
    c_packed: np.ndarray,
    res: NNResult,
    frac: float = 0.05,
    seed: int = 20260907,
) -> dict:
    """Audit the GEMM nearest-neighbour result on a random `frac` of queries.

    Two independent re-implementations are used on the FULL corpus / the reported
    neighbour:

    1. `popcount_nn_max` — integer popcount over uint64 words, full corpus, no
       floating-point matrix product.  Compares the max Tanimoto.
    2. RDKit `DataStructs.TanimotoSimilarity` on the reported nearest-neighbour
       pair — a third-party C++ implementation.

    A GEMM result that is genuinely exact yields max_abs_deviation == 0.0 on both.
    """
    rng = np.random.default_rng(seed)
    Q = q_packed.shape[0]
    k = max(1, int(round(frac * Q)))
    idx = rng.choice(Q, size=k, replace=False)

    c_u64 = np.frombuffer(np.ascontiguousarray(c_packed).tobytes(), dtype=np.uint64)
    c_u64 = c_u64.reshape(c_packed.shape[0], FP_BYTES // 8)
    c_pop = popcounts(c_packed)

    dev_popcount = 0.0
    dev_rdkit = 0.0
    argmax_mismatch_nontie = 0
    for qi in idx:
        m_pc, j_pc = popcount_nn_max(q_packed[qi], c_u64, c_pop)
        dev_popcount = max(dev_popcount, abs(m_pc - float(res.max_sim[qi])))
        j = int(res.argmax[qi])
        if j >= 0:
            s_rd = DataStructs.TanimotoSimilarity(
                bv_from_packed(q_packed[qi]), bv_from_packed(c_packed[j])
            )
            dev_rdkit = max(dev_rdkit, abs(s_rd - float(res.max_sim[qi])))
            if abs(s_rd - m_pc) > 1e-12:
                argmax_mismatch_nontie += 1
    return {
        "audited_queries": int(k),
        "audit_fraction": float(k) / Q,
        "corpus_size": int(c_packed.shape[0]),
        "max_abs_deviation_vs_popcount": float(dev_popcount),
        "max_abs_deviation_vs_rdkit": float(dev_rdkit),
        "argmax_mismatch_nontie": int(argmax_mismatch_nontie),
        "reference_implementations": [
            "numpy.bitwise_count uint64 popcount (full corpus)",
            "rdkit.DataStructs.TanimotoSimilarity (reported NN pair)",
        ],
        "seed": seed,
    }


# --------------------------------------------------------------------------- #
# statistics
# --------------------------------------------------------------------------- #
def wilson_ci(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = (z / d) * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, centre - half), min(1.0, centre + half))


def bootstrap_ci(
    values: np.ndarray,
    stat=np.mean,
    n_boot: int = 10000,
    seed_key: str = "boot",
    alpha: float = 0.05,
) -> tuple[float, float]:
    """Percentile bootstrap CI.  Seed is derived deterministically via crc32."""
    rng = np.random.default_rng(zlib.crc32(seed_key.encode()))
    values = np.asarray(values, dtype=np.float64)
    n = len(values)
    draws = rng.integers(0, n, size=(n_boot, n))
    boots = stat(values[draws], axis=1)
    return (
        float(np.quantile(boots, alpha / 2)),
        float(np.quantile(boots, 1 - alpha / 2)),
    )


def sha256_file(path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()
