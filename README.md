# Reproducibility Artifact — agenticls-rsi-leakage

**Paper:** *An LLM-Curated Self-Improvement Loop for ADMET Models, and the Leakage
Census Its Evaluation Requires* — Enrique Zueco, David Scott Lewis, Haley Yi
(AIXC Research, AIXC Bio).

**Venue:** AgenticLS — Agentic AI for Biological Discovery, NeurIPS 2026 (Workshop),
Track I (Building Agentic Systems for Life Science). Accepted as a **Spotlight**.
The workshop is non-archival.

**Camera-ready release.** This repository holds the reproducibility artifacts for the
camera-ready paper; the state the paper cites is tag `v1.0-camera-ready`. Absolute
paths are placeholders (`<REPO_ROOT>`, `<DATA_ROOT>`). The manuscript itself
(`main.tex` / `main.pdf`) is not part of this repository.

**Contact:** reports@aixcresearch.com

---

## What the paper claims

**Stage 1 (leakage census):** For each of three TDC ADMET benchmarks (BBB-Martins,
Lipophilicity-AstraZeneca, Solubility-AqSolDB), more than 20% of held-out test molecules
have a near-duplicate (Tanimoto ≥ 0.9, Morgan radius 2 / 2048 bits) in ChEMBL 36 — the
public corpus an LLM curator can draw from when generating training data.

**Stage 1b (does leakage matter):** A standard ADMET predictor trained on ChEMBL-derived
data scores better on the near-duplicated test subset than on the leakage-clean subset.
As reported in the paper, the primary solubility estimate is MAE: +0.4093 raw, +0.1595
after covariate matching, and a +0.0811 [+0.0104, +0.1567] ChEMBL-specific residual
against a corpus-free control — i.e. most of the raw gap is subset composition, not
leakage. (The pre-declaration wrote this as an "AUC/RMSE gap"; the realised endpoint is
MAE on Solubility-AqSolDB.)

**Stage 2 (RSI gains):** In a 3-round LLM-curator RSI loop, the apparent improvement on
the full TDC test set overstates the improvement on the leakage-clean subset by a margin
whose 95% CI excludes zero (+0.0216 MAE [+0.0085, +0.0353], against +0.0021 for a
random-curator control). Note the loop did **not** improve: every arm, control included,
ends three rounds worse on the pre-committed scaffold hold-out.

These claims are pre-registered in `A2_PREDECLARATION_2026-09-07.md` before any data were
collected. The observed effect is much larger than the threshold: 93.8% / 99.5% / 50.8%
of test molecules have a Tanimoto ≥ 0.9 near-duplicate in ChEMBL 36, versus 0.25% / 2.1%
/ 1.5% in the benchmark's own training split.

---

## How to reproduce

### Prerequisites

- **Python 3.14** with **RDKit** (`pip install rdkit`)
- **ChEMBL 36** SQLite database (~28 GB) — must be downloaded separately (see below)
- **TDC** package: `pip install PyTDC` (fetches the ADMET files)
- For Stage 2 (LLM arms): **Ollama** with `qwen2.5-coder:14b` and `qwen3:8b` served
  locally (thinking disabled — see `evidence/DESIGN_NOTES.md` A8)
- For the C6 second-corpus census and the C8 remediation ladder: **DrugBank 5.1**
  `full database.xml` (registered access; not redistributable — point `DRUGBANK_XML` at
  your own copy)
- `pip install scikit-learn numpy pandas scipy`

### Data that must be downloaded separately

| Dataset | Size | URL | SHA-256 |
|---------|------|-----|---------|
| ChEMBL 36 SQLite | ~28 GB | https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_36/chembl_36_sqlite.tar.gz | `cf8547cad7da721e3cfaa79b60d38f616e21e05df18059efa34428c794b16d2b` (of `.db` after extraction) |
| TDC ADMET (auto-downloaded by PyTDC) | < 10 MB per file | via `PyTDC` package | SHA-256 recorded in `evidence/DESIGN_NOTES.md` §1.2 |

### Reproduce in order

```bash
# 1. Build ChEMBL fingerprint index (run once; output ~900 MB)
cd analysis
python3 build_chembl_fps.py
python3 build_chembl_scaffolds.py
python3 build_chembl_endpoints.py

# 2. Stage 1 — leakage census
python3 s1_census.py

# 3. Stage 1b — predictive gap (leakage matters)
python3 s1b_gap.py
python3 s1b_rivals.py     # rival mechanism check

# 4. Stage 2 — leakage-conditioned RSI gains (requires Ollama)
python3 s2_rsi_loop.py
python3 s2_shuffled_feedback.py   # control arm

# 5. Supporting analyses
python3 matching_ladder.py
python3 c6_drugbank_census.py     # DrugBank overlap audit (C6)
python3 c8_remediation.py         # C8 remediation ladder (requires DrugBank 5.1 XML)
python3 ci_degeneracy_audit.py    # interval-degeneracy audit
python3 make_ledger.py            # assemble EGCL ledger

# 6. Figures and guard
python3 make_figure.py
python3 guard.py --strict     # verify every reported number against the ledger
```

---

## File inventory

| Path | Description |
|------|-------------|
| `A2_PREDECLARATION_2026-09-07.md` | Pre-registration (binding design, before any data) |
| `C7_PRECOMMITMENT_2026-09-07.md` | C7 pre-commitment (feedback-shuffled control arm) |
| `C8_PRECOMMITMENT_2026-09-07.md` | C8 pre-commitment (does the R6 remedy work?) |
| `CLAIMS.md` | Machine-readable claim list |
| `DOMAIN_AUDIT.md` | Domain validity audit |
| `EXPERIMENT_REPORT.md` | Post-hoc experiment report (headline results) |
| `analysis/s1_census.py` | Stage 1 leakage census |
| `analysis/s1b_gap.py` | Stage 1b predictive gap |
| `analysis/s1b_rivals.py` | Rival mechanism analysis |
| `analysis/s2_rsi_loop.py` | Stage 2 RSI loop |
| `analysis/s2_shuffled_feedback.py` | Stage 2 control (shuffled feedback) |
| `analysis/matching_ladder.py` | Tanimoto threshold sensitivity |
| `analysis/build_chembl_fps.py` | ChEMBL fingerprint index builder |
| `analysis/build_chembl_scaffolds.py` | ChEMBL scaffold index builder |
| `analysis/build_chembl_endpoints.py` | ChEMBL endpoint table builder |
| `analysis/c6_drugbank_census.py` | DrugBank overlap audit (C6) |
| `analysis/c8_remediation.py` | C8 remediation ladder — re-scores the loop on progressively stricter decontaminations |
| `analysis/ci_degeneracy_audit.py` | Audits every reported interval for small-n percentile degeneracy |
| `analysis/leakage_lib.py` | Shared leakage computation utilities |
| `analysis/make_ledger.py` | EGCL ledger assembler |
| `analysis/make_figure.py` | Figure generation |
| `analysis/guard.py` | Anti-fabrication guard (E1–E7) |
| `analysis/bench_data.py` | TDC benchmark data loader |
| `analysis/review_round.py` | Second-round supporting analyses |
| `evidence/DESIGN_NOTES.md` | Full provenance record, file hashes, amendments |
| `evidence/agenticls-rsi-leakage_ledger.jsonl` | EGCL ledger |
| `evidence/agenticls-rsi-leakage_stamp_manifest.json` | Stamp manifest |
| `evidence/s1_census_aggregates.json` | Stage 1 per-benchmark aggregates |
| `evidence/s1_exactness_audit.json` | Stage 1 exactness audit |
| `evidence/s1_nn_per_molecule.jsonl` | Per-molecule nearest-neighbour records (6.6 MB) |
| `evidence/s1b_gap_results.json` | Stage 1b predictive gap results |
| `evidence/s1b_gap_results_predictions.jsonl` | Per-molecule predictions (6.7 MB) |
| `evidence/s2_rsi_results.json` | Stage 2 RSI results |
| `evidence/s2_model_calls.jsonl` | Stage 2 LLM call log |
| `evidence/matching_ladder.json` | Tanimoto threshold table |
| `evidence/c6_drugbank_census.json` | DrugBank census results (C6) |
| `evidence/c8_remediation.json` | C8 remediation ladder, per-unit gains and the pre-committed reproduction check |
| `evidence/ci_degeneracy_audit.json` | Per-interval degeneracy audit |
| `evidence/figure_claims.json` | Figure data for guard |
| `refs.bib` | Bibliography (BibTeX) |
| `build.guard.log` | Anti-fabrication guard output for the submitted manuscript: every stamped cell resolved against the manifest |
| `lit/` | Literature notes, novelty check, citation verification log |

---

## Software versions (from `evidence/DESIGN_NOTES.md`)

| Component | Version |
|-----------|---------|
| Python | 3.14.7 |
| RDKit | (installed; version in `evidence/DESIGN_NOTES.md`) |
| ChEMBL | 36 (SQLite) |
| TDC ADMET | PyTDC package release as of 2026-09-07 |
| Ollama models | qwen2.5-coder:14b and qwen3:8b (Stage 2 curator arms) |
| scikit-learn | (version in `evidence/DESIGN_NOTES.md`) |
| Platform | macOS, arm64 |

---

## What is NOT included and why

| Excluded | Reason |
|----------|--------|
| `main.tex` / `main.pdf` | The manuscript is not part of the artifact repository |
| ChEMBL 36 SQLite | ~28 GB — far above this repository's size budget. Download it from the URL above, verify the SHA-256, and point `CHEMBL_DB` at the extracted `.db`; `analysis/build_chembl_fps.py` regenerates every derived index from it. |
| TDC ADMET source files | Auto-downloaded by PyTDC; SHA-256 in `evidence/DESIGN_NOTES.md` |
| DrugBank 5.1 `full database.xml` | Registered-access download, not redistributable. Required only by `c6_drugbank_census.py` and `c8_remediation.py`; the derived counts are shipped in `evidence/c6_drugbank_census.json` and `evidence/c8_remediation.json`. |
| `data/` | No downloaded corpora copied into paper directory |
| `__pycache__/` | Python bytecode |
| `.git/` | Repository history |
