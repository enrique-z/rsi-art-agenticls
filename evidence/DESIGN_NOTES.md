# DESIGN_NOTES — agenticls-rsi-leakage evidence provenance

Written 2026-09-07. Every number in this paper traces
back through this file to a real file on disk with the SHA256 recorded below.
Deviations from `A2_PREDECLARATION_2026-09-07.md` are listed as explicit
**amendments** at the end; nothing is silently changed.

---

## 1. Data — real, on disk, hashed

### 1.1 ChEMBL 36 (corpus the curator draws from)

| field | value |
|---|---|
| resolve key | `chembl_36` (via the shared dataset-path registry) |
| path | `<DATA_ROOT>/chembl/chembl_36/chembl_36_sqlite/chembl_36.db` |
| bytes | 29,739,835,392 |
| SHA256 | `cf8547cad7da721e3cfaa79b60d38f616e21e05df18059efa34428c794b16d2b` |
| `compound_structures` rows | 2,854,815 |
| rows with non-null `canonical_smiles` | 2,854,815 |
| RDKit-parseable | 2,854,800 (15 parse failures, all recorded) |
| licence | CC BY-SA 3.0, EMBL-EBI |
| `is_trap()` | returns `None` (not a known-junk file) |

The 15 unparseable structures are dropped from the corpus and counted in
`.work/agenticls/full/chembl_build.json`; they are not silently discarded.

### 1.2 TDC ADMET — flat release files (full datasets, no split)

| file | rows (excl. header) | SHA256 |
|---|---|---|
| `datasets/tdc/bbb_martins.tab` | 2,039 | `eb32676c7bf6e2bcd98b7625186bbb2ba4e4644c1d1468dfa18b6be665bd5e98` |
| `datasets/tdc/lipophilicity_astrazeneca.tab` | 4,200 | `14add1edb356e62f78256ff964c0cbba03c007be3ff48d6d4c18679f307f1d4c` |
| `datasets/tdc/solubility_aqsoldb.tab` | 9,982 | `0b593f7b55abbe477ba4284aefd7b5450551f4ed21ac4aed80a097931141c9a4` |

These match the row counts the ideation agent recorded (2039 / 4200 / 9982).

### 1.3 TDC ADMET Benchmark Group — the canonical scaffold splits (what S1 uses)

The flat `.tab` files carry **no split**. "Held-out test molecules" in the
pre-declaration means the TDC **ADMET Benchmark Group** test partition — the
split the public leaderboard scores on. It was re-downloaded live from the TDC
release rather than reconstructed:

```
pip install PyTDC==1.1.15          # python 3.11 venv
from tdc.benchmark_group import admet_group
g = admet_group(path=...)          # downloads admet_group.zip from the TDC dataverse
```

* PyTDC version: **1.1.15**, downloaded 2026-09-07.
* `admet_group.zip` SHA256 `bd0005246cb6f1672333c28ad202112a8bf071ab65d847848715030308faf1f1`
* `tdc.metadata.admet_splits` confirms all three tasks use `scaffold` splitting.
* Installed centrally at `datasets/tdc/admet_group/` and registered as resolve
  key `tdc_admet_group` in the shared dataset-path registry (no dataset is
  copied into the paper folder).

| task | train_val | test | SHA256 (test.csv) |
|---|---|---|---|
| `bbb_martins` | 1,624 | 406 | `82d4fa96f671a8382e9d8e6de3dde8385b9f337854c83e6fde2e69a7be93754f` |
| `lipophilicity_astrazeneca` | 3,360 | 840 | `69f3870c9b6f2cec79a992a31d59d2e940962f7e2c4b27dabaed1fe1e767e2d8` |
| `solubility_aqsoldb` | 7,985 | 1,997 | `18cd02deec4e06c611a6e7e38fb44a1731bc0c0a96bc5d4ded098ff29da645fe` |

train_val.csv SHA256: `5035139072311f4cb0a1a54cffe11bc727004886ace4629b98f6db20c431e4af`
(bbb), `a54636f85cf31eb4c914ea15e2c3a15364b0b8b51986262ea229203637e08f31` (lipo),
`adbee092423fe0cb0f4219520f2f93730b6ba5bb1d0ef8dd7764c2d32b3393b9` (sol).

**Reconciliation with the flat files.** lipo 3,360+840 = 4,200 = the `.tab` row
count; sol 7,985+1,997 = 9,982 = the `.tab` row count. BBB is 1,624+406 = 2,030
against 2,039 in the `.tab` file — the benchmark group drops 9 rows. This
discrepancy is recorded rather than papered over; S1 reports on the 406-molecule
benchmark-group test split, which is what the leaderboard uses.

### 1.4 ChEMBL-derived endpoint pools (S1b training data)

Extracted by `analysis/build_chembl_endpoints.py` straight from `activities`:

| endpoint | filter | activity rows | compounds |
|---|---|---|---|
| logD 7.4 | `standard_type='LogD'`, `standard_relation='='`, value ∈ [-5,10], assay description contains `7.4` | 19,007 | 17,895 |
| log S (mol/L) | `standard_type='Solubility'`, `standard_relation='='`, units ∈ {nM, ug.mL-1}, converted with `compound_properties.full_mwt`, log10(mol/L) ∈ [-12,2] | 42,819 (9 dropped out of range) | 34,658 |
| BBB penetration | — | **202 candidate rows total** | **not runnable** |

Per-compound label = **median** of that compound's measurements.

Unit conversion for solubility: `nM -> mol/L = v*1e-9`;
`ug.mL-1 -> mol/L = (v*1e-3)/full_mwt`. AqSolDB's `Y` is already log10(mol/L),
so the two are on the same scale.

---

## 2. Method

### 2.1 Fingerprints and canonicalisation

* RDKit **2026.03.1**, `rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)`.
  (The older `GetMorganFingerprintAsBitVect` is deprecated in this RDKit.)
* Canonical form = `Chem.MolToSmiles(Chem.MolFromSmiles(s))`.
* Fingerprints stored packed, 256 bytes per molecule (`np.packbits`).
* Molecules RDKit cannot parse are dropped **and counted** at every stage.

### 2.2 Exact nearest-neighbour Tanimoto — and why it is exact

Tanimoto(a,b) = |a∧b| / (|a|+|b|−|a∧b|).

The intersection matrix |a∧b| is computed as a **float32 GEMM** on the 0/1
expansions of the fingerprints. Every entry is an integer in [0, 2048]; all
integers below 2²⁴ are exactly representable in float32, so the GEMM is exact
integer arithmetic, not an approximation. The single division is then performed
in float64, so the result is bit-identical to an integer implementation.

**This is a full exact search, not an approximate index**: every one of the
16,210 benchmark molecules — the test *and* train_val splits of all three tasks,
3,241 + 12,969 after canonicalisation — is compared against all 2,854,800 ChEMBL
fingerprints, i.e. **46,276,308,000 pairs**.

> *Corrected 2026-09-07 by the review round.* This paragraph previously read
> "15,212 benchmark molecules ... (4.34 × 10¹⁰ pairs)". That was written against
> an early draft of `s1_census.py` that censused the test splits only; the
> shipped run censuses the training splits as well, because the own-train-split
> rate is the reference scale the paper's central contrast depends on. The
> ledger row `s1_corpus::n_pairwise_comparisons` has always carried the correct
> value (46,276,308,000 = 16,210 × 2,854,800) and the manuscript quotes it from
> there through a stamped cell, so no number in the PDF was ever affected — only
> this provenance note was stale. No ledger row changed.

**Exactness audit** (`evidence/s1_exactness_audit.json`). Because the
pre-declaration required an audit if any approximation were used, the audit was
run anyway on a random **5%** of queries, against two *independent*
reimplementations over the **full** corpus:

1. `numpy.bitwise_count` integer popcount over uint64 words (no matrix product);
2. RDKit's own C++ `DataStructs.TanimotoSimilarity` on the reported NN pair.

Acceptance criterion, asserted in code (`s1_census.py` aborts otherwise):
`max_abs_deviation == 0.0` against both. Results are in the JSON.

### 2.3 Operational definitions

| term | definition |
|---|---|
| **leaked** test molecule | max Tanimoto to **any** ChEMBL 36 compound ≥ 0.9 |
| **non-leaked** | the same quantity < 0.9 |
| **exact duplicate** | identical RDKit canonical SMILES present in ChEMBL 36 `compound_structures` |
| **Tanimoto = 1.0** | fingerprint-identical (weaker than exact duplicate — folded 2048-bit fingerprints can collide); both are reported separately |
| **scaffold** | Bemis-Murcko, `MurckoScaffold.GetScaffoldForMol` |
| **scaffold-level leak rate** | fraction of *distinct* test scaffolds with ≥1 member at NN ≥ 0.9 (guards against one over-represented scaffold carrying the rate) |
| **gap** (regression) | MAE(non-leaked) − MAE(leaked); positive ⇒ leaked subset looks better |
| **gap** (classification) | AUC(leaked) − AUC(non-leaked); positive ⇒ leaked subset looks better |

### 2.4 Statistics

* Binomial rates: Wilson score 95% intervals.
* S1b gap: two CIs, both reported, never one silently chosen —
  * **seed-level** (the pre-registered unit of independence): percentile
    bootstrap over the per-seed gaps;
  * **molecule-level**: bootstrap resampling test molecules inside each subset,
    on seed-averaged predictions.
* All bootstrap seeds derived with `zlib.crc32(key)` — never Python `hash()`,
  which is per-process randomised.

### 2.5 Software

python 3.14.7 · numpy 2.3.4 (OpenBLAS 0.3.30) · scipy 1.17.1 · scikit-learn 1.8.0 ·
rdkit 2026.03.1 · PyTDC 1.1.15 (python 3.11 venv, splits only) ·
macOS 15.6.1 arm64 (Apple M3, 16 cores, 64 GB).

---

## 3. Pilots (pre-declaration §Feasibility bound)

| stage | pilot | wall time | outcome |
|---|---|---|---|
| ChEMBL fingerprinting | 100,000 compounds | 88 s | PASS — 0 parse failures |
| S1 census | BBB × 100 test molecules × 100,000 ChEMBL | 25 s | PASS — exactness deviation 0.0 vs both reference implementations |
| S1b | Lipophilicity, 2 seeds, small forests | — | PASS on the first arm (pipeline verified end to end), then killed to free CPU for the full S1 census; re-run at full size |

The first S1 pilot **failed its own exactness assertion** at 2.9e-8 (float32
division rounding). Rather than relax the assertion, the division was moved to
float64; the re-run gave 0.0. This is recorded because it is exactly the kind of
gate that must never be loosened to make a run finish.

---

## 4. Amendments to the pre-declaration

**A1 — "held-out test" is the TDC ADMET Benchmark Group scaffold test split.**
The pre-declaration cites the flat `.tab` files (2039/4200/9982 rows) but those
carry no split. S1 uses the canonical benchmark-group test partitions
(406/840/1997). This is the stricter and more meaningful reading — it is the
partition the public leaderboard reports on. Both are hashed above.

**A2 — S1b arm A is not runnable for BBB-Martins.** The pre-declaration assumed
a ChEMBL-derived training pool for "the same endpoint" for all three benchmarks.
ChEMBL 36 contains only **202** blood-brain-barrier distribution activity rows
across incompatible `standard_type`s (`logBB`, `Log BB`, `pBBB`, ...), which
cannot train a predictor. BBB therefore runs **arm B only** (train on the
benchmark's own train_val split). This is reported as a limitation, not hidden.

**A3 — S1b runs three arms, not one.** The pre-declaration specified one
ChEMBL-trained arm. Three are run so the pre-declared rival mechanism
"chemistry-overlap-is-normal" is *measured* rather than argued:
* **A2 (primary)** ChEMBL pool with exact TDC-test canonical-SMILES duplicates
  removed — the literal reading of "never on TDC test";
* **A1 (supplementary)** the same pool *without* that removal — what a naive LLM
  curator actually trains on;
* **B (control)** the benchmark's own train_val split — if the leaked/non-leaked
  gap appears here too, the gap is not ChEMBL-training leakage.

**A4 — 10 seeds, not 5.** The pre-declaration set a floor of 5.

**A5 — the second model family is not gradient boosting.** The pre-declaration
said "random forest on ECFP4, plus a second model family". `HistGradientBoosting`
was the first choice and was **dropped after measurement**: one 200-iteration fit
on the 1,462-molecule BBB training pool took **1,259 s** on this machine
(sklearn 1.8.0, 2048 binary columns), putting 10 seeds × 7 arm/model combinations
at roughly 24 h. It is replaced by two faster families that are, if anything,
more informative:
* `linear_ecfp` — L2-regularised linear model on the same bits (Ridge for
  regression, logistic regression for classification): a standard cheminformatics
  baseline and genuinely not a tree ensemble;
* `knn_tanimoto` — k=5 similarity-weighted nearest neighbours under the *same*
  Tanimoto metric the census uses. This is the mechanism probe: it is the model
  that can exploit a near-duplicate most directly, so it upper-bounds what a
  near-duplicate is worth. It is recorded with `must_print=False` and is never
  the headline.

Timing evidence for the replacement, all measured on this machine on 2026-09-07
against the real ChEMBL logD pool (17,895 × 2048, 200 trees): random forest dense
13.7 s, the same fit on a CSC sparse matrix 65.4 s, and with sklearn's default
`max_features=1.0` 510.2 s. The code uses dense input and `max_features='sqrt'`,
and says so in a comment beside the choice.

**A6 — the unit of independence for the S1b gap is the MOLECULE, not the seed.**
The pre-declaration named the seed. That is wrong for this quantity: the test
molecules are identical across seeds, so a seed-level bootstrap resamples only
model-fitting noise and yields intervals roughly an order of magnitude too
narrow — on Lipophilicity-AZ it "excludes zero" on a four-molecule non-leaked
subset. **Both** intervals are computed and both are in the ledger
(`gap_ci_seed_bootstrap`, `gap_ci_molecule_bootstrap`); every verdict in
`EXPERIMENT_REPORT.md` and every interval drawn in the figure uses the
molecule-level one. The seed-level interval is retained, labelled, and explicitly
described as understating uncertainty, rather than quietly deleted.

`MIN_SUBSET_N = 20` was set in code *before* the S1b run as the size below which a
subset metric is flagged `underpowered`; it has not been changed since, even
though BBB-Martins' 25 non-leaked molecules clear it and are still too few for a
stable AUC. The molecule-level CI is what carries that judgement.

**A7 — S2 runs on Solubility-AqSolDB, not Lipophilicity.** S2 was designed against
the lipophilicity endpoint. S1 then measured **4** non-leaked molecules in the
Lipophilicity-AZ test split and 25 in BBB-Martins, which makes "the gain on the
non-leaked subset" unmeasurable there. Solubility-AqSolDB (1,013 leaked / 982
non-leaked) is the only benchmark on which the pre-registered S2 comparison can
be made at all, so S2 runs there. The switch was forced by an S1 measurement, not
by an S2 result.

**A8 — S2 curator models.** The pre-declaration allowed "local Ollama models
and/or an API model — recorded before the run". Probed on 2026-09-07 with a
12-cluster prompt: `qwen2.5-coder:14b` 22.7 s and a well-formed selection;
`qwen3:8b` 1.0 s with thinking disabled (86.9 s with thinking on, and the whole
token budget consumed inside `<think>`); `gpt-oss:20b` 277.6 s with thinking
disabled. `gpt-oss:20b` was dropped as too slow for 27 calls. Curators used:
**`qwen2.5-coder:14b` and `qwen3:8b`**, plus a random-curator control arm with an
identical budget and an identical number of curation steps. Reasoning models are
called with thinking off and `num_predict=800`.

---

## 6. Pilot failures, and what was done about them

Recorded because "the gate fired" is the system working, and because each was
fixed at the source rather than by loosening the check.

| pilot | failure | fix |
|---|---|---|
| S1 census | exactness assertion tripped at 2.93e-8 | division moved from float32 to float64; re-run gave exactly 0.0. The assertion was not relaxed. |
| S1b | `AssertionError: prediction rows missing for bbb_martins/…` | the predictions join was keyed on canonical SMILES, and the BBB-Martins test split contains **12 internally duplicated structures** (406 rows → 394 distinct), so rows were silently dropped. Join is now positional with an asserted SMILES-sequence match. The duplicate count is itself a ledger row. |
| S1b | one HGB fit took 1,259 s | second model family replaced — amendment A5, with timings. |
| S2 | `qwen3:8b` timed out at 300 s inside `<think>`; the arm silently fell back to "select every cluster" | curators re-probed and re-chosen (A8); thinking disabled; token cap added; fallback use is now counted in the results (final run: 18/18 calls succeeded, 0 fallbacks). |

---

## 5. Files

| artefact | what |
|---|---|
| `evidence/s1_nn_per_molecule.jsonl` | one row per benchmark molecule: NN Tanimoto to ChEMBL, NN identity, exact-duplicate flag, scaffold, own-train-split NN |
| `evidence/s1_census_aggregates.json` | per-benchmark rates, Wilson CIs, full similarity histograms |
| `evidence/s1_exactness_audit.json` | the NN exactness audit |
| `evidence/s1b_gap_results.json` | per-seed and aggregated leaked/non-leaked scores and gaps |
| `evidence/s1b_gap_results_predictions.jsonl` | seed-averaged per-molecule predictions |
| `evidence/s1b_rival_mechanisms.json` | label agreement, covariate distributions, matched gap |
| `evidence/agenticls_rsi_leakage_ledger.jsonl` | EGCL ledger rows |
| `evidence/agenticls-rsi-leakage_stamp_manifest.json` | stamped-cell manifest |
| `.work/agenticls/full/` | derived fingerprint matrices (NOT committed; regenerate with `analysis/build_chembl_fps.py`) |
