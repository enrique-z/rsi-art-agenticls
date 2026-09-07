# EXPERIMENT REPORT — agenticls-rsi-leakage

Written 2026-09-07. Provenance, hashes, operational
definitions and the four pre-declaration amendments are in
`evidence/DESIGN_NOTES.md`.

**Headline.** The near-duplicate leakage the pre-declaration hypothesised is not
a marginal 20% effect. On the canonical TDC ADMET Benchmark Group scaffold test
splits, **93.8% / 99.5% / 50.8%** of test molecules have a Tanimoto ≥ 0.9
near-duplicate somewhere in ChEMBL 36 — the corpus an LLM curator draws from —
while the same molecules have a ≥ 0.9 neighbour in the benchmark's *own* training
split only **0.25% / 2.1% / 1.5%** of the time. The scaffold split protects the
benchmark from itself and does nothing at all against the public corpus.

---

## Stage 1 — leakage census (S1)

### Pilot
60-second foreground pilot: BBB-Martins × 100 test molecules × 100,000 ChEMBL
compounds. **First attempt FAILED** its own exactness assertion at max deviation
2.93e-8 (float32 division rounding). The assertion was not relaxed — the division
was moved to float64 and the re-run gave **exactly 0.0** against both reference
implementations, in 25 s. PASS.

### Method
Morgan r=2, 2048-bit fingerprints (RDKit 2026.03.1
`rdFingerprintGenerator.GetMorganGenerator`) over RDKit-canonicalised SMILES.
**Exact** nearest-neighbour Tanimoto of every benchmark molecule against the
whole ChEMBL 36 compound corpus — no approximate index, no sampling:
15,212 benchmark molecules × 2,854,800 ChEMBL fingerprints
= **4.34 × 10¹⁰ pairwise comparisons**, 2,413 s of dense float32 GEMM.

The intersection matrix is computed as a float32 matrix product of the 0/1
fingerprint expansions; every entry is an integer ≤ 2048 and therefore exactly
representable, and the single division is done in float64. The result is
bit-identical to an integer implementation.

**Exactness audit** (`evidence/s1_exactness_audit.json`), run even though the
search is exact, on a random **5%** (810) of queries against the **full** corpus
using two independent reimplementations:

| reference | max abs deviation |
|---|---|
| `numpy.bitwise_count` uint64 popcount, full corpus | **0.0** |
| RDKit `DataStructs.TanimotoSimilarity` on the reported NN pair | **0.0** |

Non-tie argmax mismatches: **0**.

### Results

Test split = TDC ADMET Benchmark Group scaffold test partition.

| benchmark | n test | NN≥0.9 vs ChEMBL (Wilson 95%) | NN≥0.95 | NN=1.0 | exact canonical-SMILES duplicate | **NN≥0.9 vs OWN train split** |
|---|---|---|---|---|---|---|
| BBB-Martins | 406 | **93.84%** [91.07, 95.79] | 92.86% | 91.63% | 85.47% | **0.25%** |
| Lipophilicity-AZ | 840 | **99.52%** [98.78, 99.81] | 99.52% | 99.52% | 99.29% | **2.14%** |
| Solubility-AqSolDB | 1,995 | **50.78%** [48.58, 52.97] | 49.67% | 48.97% | 42.26% | **1.50%** |

Median nearest-neighbour Tanimoto to ChEMBL: 1.000 (BBB), 1.000 (lipo),
0.936 (solubility). The 5th percentile is 0.729 / 1.000 / 0.420 — for
Lipophilicity-AZ, *even the least-similar 5% of test molecules are exact
fingerprint matches*.

Scaffold-clustered (guards against one over-represented scaffold carrying the
rate): 93.72% / 99.39% / 58.55% of **distinct Bemis-Murcko test scaffolds** have
at least one member at NN ≥ 0.9. The largest single test scaffold is 24.4% of
BBB, 2.4% of lipo, 7.6% of solubility, so the rates are not scaffold artefacts.

The benchmarks' own *train_val* splits leak at the same rate (91.56% / 99.70% /
53.36% at ≥0.9), so this is a property of the datasets, not of the test partition.

### Verdict — S1: **HELD**, by a wide margin

The pre-declaration set the bar at ">20% of held-out test molecules have a
Tanimoto ≥ 0.9 near-duplicate in ChEMBL 36". Measured: 93.8%, 99.5%, 50.8% —
the lower Wilson bound of the *smallest* is 48.6%, still 2.4× the pre-registered
threshold. All three benchmarks clear it.

One consequence is worth stating plainly, because it is the paper's sharpest
finding and it also constrains S1b: for Lipophilicity-AstraZeneca the
benchmark-group test split contains **4 molecules** that are not near-duplicates
of a ChEMBL compound. 838 of its 840 test molecules are present in ChEMBL 36 *by
their own ChEMBL identifier*, and 149 of them carry a logD7.4 measurement in
ChEMBL 36 — the label itself is one SQL query away.

---

## Stage 1b — does the leakage inflate scores? (S1b)

### Pilot
Full-coverage pilot (all three benchmarks, both tasks, all arms, all three model
families, 2 seeds, small models): 575 s, every code path exercised. PASS.
An earlier pilot exposed two real defects that were fixed rather than worked
around: a canonical-SMILES join that silently dropped rows (BBB-Martins' test
split contains 12 internally duplicated structures — 406 rows, 394 distinct), and
a `HistGradientBoosting` fit that took 1,259 s (amendment A5).

### Method
Three training arms, three model families, 10 seeds (pre-declared floor 5), 90%
subsample of the training pool per seed. Leaked = NN-to-ChEMBL ≥ 0.9 (the S1
definition). Gap = MAE(non-leaked) − MAE(leaked) for regression,
AUC(leaked) − AUC(non-leaked) for classification; positive means the
near-duplicate subset looks better.

| arm | training pool |
|---|---|
| **A2 (primary)** | ChEMBL endpoint pool with exact TDC-test canonical-SMILES duplicates removed |
| A1 | the same pool without that removal (naive curator) |
| B (control) | the benchmark's own train_val split — never touches ChEMBL |

Families: `random_forest` (200 trees, `max_features='sqrt'`), `linear_ecfp`
(Ridge / logistic regression on the same bits), `knn_tanimoto` (k=5,
similarity-weighted, the mechanism probe).

### Result — the binary split is degenerate on two of three benchmarks

| benchmark | n leaked | n non-leaked |
|---|---|---|
| BBB-Martins | 381 | **25** |
| Lipophilicity-AZ | 836 | **4** |
| Solubility-AqSolDB | 1,013 | 982 |

That is a finding, not an inconvenience: on Lipophilicity-AZ a "leakage-free"
evaluation subset does not exist. Only Solubility-AqSolDB supports the
pre-registered comparison.

### Solubility-AqSolDB — the estimable case

MAE, mean of 10 seeds; CI is the **molecule-level** bootstrap.

| arm | model | leaked | non-leaked | gap | 95% CI |
|---|---|---|---|---|---|
| **A2 ChEMBL dedup** | RF | 1.543 | 1.952 | **+0.409** | [+0.295, +0.526] |
| A2 ChEMBL dedup | linear | 1.412 | 1.802 | +0.389 | [+0.281, +0.504] |
| A2 ChEMBL dedup | kNN-Tc | 1.660 | 2.260 | +0.599 | [+0.464, +0.739] |
| A1 ChEMBL naive | RF | 1.516 | 1.945 | +0.429 | [+0.317, +0.545] |
| **B TDC own train** | RF | 1.198 | 1.447 | **+0.250** | [+0.164, +0.338] |
| B TDC own train | linear | 1.518 | 1.446 | −0.072 | [−0.175, +0.032] |
| B TDC own train | kNN-Tc | 1.254 | 1.493 | +0.239 | [+0.142, +0.333] |

All three families agree on the ChEMBL-trained arms and every CI excludes zero.
Removing exact test duplicates from the pool (A1 → A2) barely moves the gap
(+0.429 → +0.409): **the effect is near-duplicates, not verbatim copies.**

### Verdict — S1b: **MIXED**

* **Solubility-AqSolDB — HELD.** Gap +0.409 MAE, molecule-level 95% CI
  [+0.295, +0.526], reproduced by three model families, and surviving covariate
  matching (below).
* **BBB-Martins — NOT RESOLVED.** 25 non-leaked molecules. The sign is not even
  stable across families (RF −0.050, kNN −0.087, linear +0.054) and every
  molecule-level CI spans zero. Reported as unresolved, not as a null.
* **Lipophilicity-AZ — NOT ESTIMABLE.** 4 non-leaked molecules; the CI
  [−0.49, +1.43] is uninformative by construction.

**A methodological correction we are obliged to state (amendment A6).** The
pre-declaration named the *seed* as the unit of independence. That is wrong for
this quantity: the test molecules are identical across seeds, so a seed-level
bootstrap captures only model-fitting noise and produces intervals roughly ten
times too narrow — on Lipophilicity-AZ it "excludes zero" on a four-molecule
subset. Both intervals are in the ledger; **the molecule-level interval is the
one to trust**, and it is what the figure and every verdict above use.

### The continuous alternative (reported because the binary split fails)

Spearman ρ between NN-to-ChEMBL similarity and per-molecule absolute error,
over all test molecules (seed-averaged predictions):

| benchmark | A2 RF | B RF |
|---|---|---|
| Solubility-AqSolDB | **−0.219** (p = 4e−23) | −0.105 (p = 2e−6) |
| Lipophilicity-AZ | −0.014 (p = 0.68) | −0.050 (p = 0.15) |
| BBB-Martins | — | −0.087 (p = 0.08) |

For solubility the ChEMBL-trained model's error falls monotonically with
similarity to ChEMBL, twice as steeply as the control arm's.

### Rival mechanisms (all three pre-declared, all three measured)

**R1 — "chemistry overlap is normal."** Partly true, and quantified rather than
argued. On solubility the control arm B — trained only on the benchmark's own
train split, never on ChEMBL — still shows a +0.250 gap, i.e. **61% of the RF gap
is generic** ("molecules that resemble ChEMBL are easier for any model"), leaving
a ChEMBL-training-specific excess of ≈ +0.16 MAE. Neither the naive "it's all
leakage" nor the dismissive "it's all chemistry" reading survives.

**R2 — endpoint-value mismatch (activity cliffs). The falsifier fails to fire.**
Mean |Δ label| for near-duplicate pairs vs random pairs, paired bootstrap:

| pairing | near-duplicate | random | difference (95% CI) |
|---|---|---|---|
| Lipophilicity test ↔ ChEMBL logD pool (157 pairs) | **0.081** | 1.618 | −1.537 [−1.748, −1.334] |
| Solubility test ↔ ChEMBL solubility pool (124 pairs) | **0.674** | 1.567 | −0.893 [−1.139, −0.648] |
| Solubility test ↔ own train split (30 pairs) | 1.671 | 2.708 | −1.037 [−1.939, −0.114] |

A structural near-duplicate in ChEMBL carries a logD value **0.08 log units**
from the benchmark's own label — inside experimental error. The pre-declared
discriminator ("if near-duplicate pairs agree no better than random pairs the
mechanism is falsified") is decisively not met: the duplicates carry the answer.

**R3 — subset composition.** Leaked and non-leaked solubility molecules do differ
(MW 282 vs 339, SMD −0.25; TPSA 64 vs 78, SMD −0.19; label mean −3.10 vs −3.72,
SMD +0.27), though their ChEMBL scaffold frequencies are nearly identical
(SMD −0.04). Coarsened-exact matching on MW quintile × label quintile
(796 molecules per side) shrinks the gap but does not remove it:

| arm · model | raw gap | matched gap (95% CI) |
|---|---|---|
| A2 ChEMBL dedup · RF | +0.409 | **+0.160** [+0.133, +0.189] |
| A2 ChEMBL dedup · kNN-Tc | +0.599 | +0.293 [+0.251, +0.334] |
| A2 ChEMBL dedup · linear | +0.389 | +0.138 [+0.110, +0.175] |
| B TDC own train · RF | +0.250 | +0.069 [+0.039, +0.094] |

So ≈ 61% of the raw gap is composition, and of the matched remainder roughly half
again is generic — leaving a small but robustly non-zero ChEMBL-training-specific
component. For BBB (25 non-leaked) every matched CI spans zero, and for
Lipophilicity-AZ matching is impossible (4 molecules per side).

---

## Stage 2 — leakage-conditioned RSI gains (S2)

### Pilot
Two pilots. The **first failed**: `qwen3:8b` spent its whole token budget inside
`<think>` and the call timed out at 300 s, so the arm silently fell back to
"select everything". Rather than accept that, the curators were re-probed
(`qwen2.5-coder:14b` 22.7 s and a well-formed selection; `qwen3:8b` 1.0 s with
thinking disabled; `gpt-oss:20b` 277.6 s — too slow for 27 calls, dropped) and
the harness now disables thinking and caps `num_predict`. The second pilot
returned a real curated selection in 13.7 s. PASS.

### Method
Ran on **Solubility-AqSolDB** (amendment A7) — the only benchmark with a
non-leaked subset large enough to measure a gain on at all.

* Curation universe: the 34,658-compound ChEMBL solubility pool, minus a
  **pre-committed scaffold hold-out** (deterministic rule
  `sha256("agenticls-rsi-leakage|" + scaffold)[:8] % 5 == 0` → 3,609 of 17,612
  scaffolds, 7,338 compounds, hold-out SHA256 recorded in
  `evidence/s2_scaffold_holdout_commit.json` before the loop ran), minus the 90
  compounds that are exact TDC-test duplicates. Universe: 27,230 compounds,
  MiniBatchKMeans into 40 clusters.
* Round 0 = fixed-seed random draw of 6,000 compounds. Rounds 1–3: the curator
  sees each cluster's size, label statistics, an example SMILES, and **the current
  model's MAE on the benchmark molecules nearest that cluster** — real benchmark
  feedback — and returns the clusters to draw the next 6,000 from.
* Arms: `qwen2.5-coder:14b`, `qwen3:8b`, and a **random-curator control** with an
  identical budget and an identical number of curation steps. 3 seeds each.
* **18/18 LLM calls succeeded**, zero parse fallbacks, median latency 46 s. Every
  prompt hash, raw response, latency and selection is in
  `evidence/s2_model_calls.jsonl`.

### Results — round 3 minus round 0 (positive = MAE improved)

| arm | full test | non-leaked subset | scaffold hold-out | full − non-leaked |
|---|---|---|---|---|
| random control | −0.0104 | −0.0125 | −0.0235 | +0.0021 |
| qwen2.5-coder:14b | −0.0059 | −0.0241 | −0.0482 | +0.0182 |
| qwen3:8b | **+0.0094** | −0.0156 | −0.0204 | +0.0250 |

Pooled over the 6 LLM curator × seed units: overstatement
**+0.0216**, bootstrap 95% CI **[+0.0085, +0.0353]**, excludes zero.

### Verdict — S2: **MIXED, and the caveat is larger than the result**

The pre-registered directional claim is met: the improvement measured on the full
TDC test set exceeds the improvement measured on the non-leaked subset, with a
paired CI excluding zero, and the random-curator control shows almost none of it
(+0.0021 vs +0.0216) — so the bias is specific to benchmark-feedback-driven
curation, not to reshuffling the pool.

**But there is no gain to inflate.** Every arm *loses* accuracy over three rounds
on the pre-committed scaffold hold-out (−0.020 to −0.048 MAE), and only
`qwen3:8b` improves at all on the full test set (+0.0094). The honest statement
is therefore about *measurement*, not about performance: benchmark-feedback
curation degrades the model, and the leaked benchmark hides part of that
degradation. Reporting this as "RSI gains are leakage-driven" would be an
overclaim, and the paper must not make it.

Two further limits we state rather than bury: n = 6 curator × seed units is a
small basis for a percentile bootstrap, and two of three `qwen2.5-coder:14b` seeds
converged to an identical 3-cluster selection by round 1 and never moved again —
this RSI loop stalls at temperature 0.

---

## Summary of verdicts

| stage | claim | verdict |
|---|---|---|
| **S1** | >20% of test molecules have a ≥0.9 near-duplicate in ChEMBL 36 | **HELD** — 93.8% / 99.5% / 50.8% |
| **S1b** | ChEMBL-trained predictors score better on the leaked subset | **MIXED** — held on Solubility-AqSolDB (+0.409 MAE, CI [+0.295, +0.526], ≈39% surviving covariate matching and ≈39% ChEMBL-specific); unresolved on BBB (n=25); not estimable on Lipophilicity-AZ (n=4) |
| **R2 falsifier** | near-duplicates carry no exploitable label signal | **NOT TRIGGERED** — near-duplicate logD differs by 0.081 vs 1.618 for random pairs |
| **S2** | apparent RSI gain overstates the non-leaked gain | **MIXED** — direction held (+0.0216, CI [+0.0085, +0.0353]) but all arms lose accuracy on the pre-committed hold-out, so the finding is measurement bias, not inflated gains |

---

## Reproducing

```bash
cd analysis
python3 build_chembl_fps.py       --out <work>          # ~5 min, 12 procs
python3 build_chembl_endpoints.py --work <work>         # ~30 s
python3 build_chembl_scaffolds.py --work <work>         # ~8 min
python3 s1_census.py --pilot --work <pilot_work>        # 25 s pilot, must pass first
python3 s1_census.py        --work <work>               # ~50 min
python3 s1b_gap.py   --pilot                            # pipeline check
python3 s1b_gap.py   --seeds 10                         # ~2 h
python3 s1b_rivals.py
python3 s2_rsi_loop.py --pilot                          # then without --pilot
python3 make_ledger.py && python3 make_figure.py
```

`<work>` defaults to `<REPO_ROOT>/.work/agenticls/full` and
holds only derived arrays (fingerprints, scaffold counts) — no dataset is copied
into the paper folder, and every dataset read goes through
`dataset_paths.resolve()`.
