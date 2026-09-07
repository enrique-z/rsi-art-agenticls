# Claim–Evidence–Risk Registry — agenticls-rsi-leakage

Every number in `main.tex` is a `\StampedCell` resolving against
`evidence/agenticls-rsi-leakage_stamp_manifest.json` (205 rows, all
`verdict=SUPPORTED`, 45 flagged `must_print`; all 45 appear in the PDF or in
`evidence/figure_claims.json`). Interval bounds that live in the `ci_low` /
`ci_high` **fields** of a ledger row rather than as their own manifest keys are
written as `\SafeNum` and quoted from that same row; they are listed in the
"Unexported interval bounds" section at the end.

## Registry

| # | Claim text (verbatim or close) | Type | Evidence (ledger rows) | Statistical treatment | Reviewer-objection risk | Limitation sentence location |
|---|-------------------------------|------|------------------------|-----------------------|-------------------------|------------------------------|
| 1 | "93.84% / 99.52% / 50.78% of scaffold-split test molecules have a Tanimoto ≥ 0.9 near-duplicate in ChEMBL 36" | positive (pre-registered, HELD) | ledger rows `s1_bbb_martins::leak_rate_0p9`, `s1_lipophilicity_astrazeneca::leak_rate_0p9`, `s1_solubility_aqsoldb::leak_rate_0p9` | Wilson 95% intervals; unit of independence = molecule; exact (not approximate) NN search, audited at max deviation 0.0 on a random 5% of queries | (1) "similarity threshold is arbitrary" — whole survival curve reported in Fig. 1A, plus rates at 0.95 and 1.0; (2) "fingerprint folding collides" — exact canonical-SMILES duplicate rate reported separately; (3) "one over-represented scaffold carries it" — scaffold-level rates reported | §9 Limitations (one corpus, one fingerprint, one threshold); §1 "What we claim, and what we do not" |
| 2 | "against 0.25% / 2.14% / 1.50% in each benchmark's own training split" | diagnostic (the load-bearing contrast) | rows `s1_*::own_trainval_rate_0p9` | Same Wilson treatment, same molecules, different reference corpus | "the comparison is unfair, ChEMBL is far larger" — conceded and reframed: the claim is about the corpus the curator can reach, not about corpus size | §9 Limitations; §7 R1 |
| 3 | "ChEMBL-trained RF: MAE 1.5432 leaked vs 1.9524 non-leaked, gap +0.4093 [+0.2949, +0.5262]" | positive (pre-registered S1b, HELD on one benchmark) | rows `s1b_solubility_aqsoldb_A2_chembl_dedup_random_forest::{score_leaked,score_nonleaked,gap,gap_mol_ci_low,gap_mol_ci_high}` | Molecule-level bootstrap 95% CI (amendment A6: the pre-declared seed-level unit was wrong and produces intervals ~10× too narrow); 10 seeds, 3 model families; CI excludes zero | (1) "seeds are not independent units" — corrected in A6, both intervals in the ledger, molecule-level used everywhere; (2) "cherry-picked benchmark" — the other two are reported as unresolved/not estimable, not dropped | §9 Limitations; §1 "What we claim, and what we do not" |
| 4 | "the effect is near-duplicates, not verbatim copies (A1 +0.4293 vs A2 +0.4093)" | diagnostic | rows `s1b_solubility_aqsoldb_A1_chembl_naive_random_forest::gap`, `s1b_solubility_aqsoldb_A2_chembl_dedup_random_forest::gap` | Same bootstrap; arms differ only by removal of exact TDC-test SMILES from the training pool | "exact-duplicate removal is the standard fix" — that is the point: it does not fix this | §9 Limitations; §1 "What we claim, and what we do not" |
| 5 | "BBB-Martins is unresolved (n non-leaked = 25), Lipophilicity-AZ not estimable (n = 4)" | negative / honest non-result | rows `s1b_bbb_martins::n_nonleaked`, `s1b_lipophilicity_astrazeneca::n_nonleaked`, `s1b_bbb_martins_B_tdc_trainval_*::gap*` | Molecule-level CIs all span zero; sign unstable across three families | "why report an unpowered arm at all" — because dropping it would be selection; reported as unresolved, explicitly not as a null | §9 Limitations; §4 "BBB-Martins: unresolved"; §1 claim boundary |
| 6 | "matched gap +0.1595 [+0.1327, +0.189]; control arm matched gap +0.0689 [+0.0393, +0.0941]" | diagnostic (argues AGAINST the strong form of our own claim) | rows `r3_solubility_aqsoldb_A2_chembl_dedup_random_forest::matched_gap`, `r3_solubility_aqsoldb_B_tdc_trainval_random_forest::matched_gap` (+ `ci_low`/`ci_high` fields) | Coarsened-exact matching, MW quintile × label quintile, 796 molecules per side; bootstrap 95% CIs; the two matched intervals do not overlap (weaker than a formal difference test — stated as such in the paper) | "composition confound" — this row IS the composition analysis; "non-overlapping CIs is not a test" — the paper says so explicitly | §9 Limitations (last clause); §7 R3 |
| 7 | "near-duplicate pairs differ by 0.0808 log units vs 1.6181 for random pairs" | diagnostic (pre-declared falsifier, NOT triggered) | rows `r2_lipophilicity_astrazeneca_vs_chembl_endpoint_pool::{near_pair_stat,random_pair_stat,pair_diff}` | Paired bootstrap 95% CI on the difference, 157 near-duplicate pairs vs the same number of random pairs | "activity cliffs make near-duplicates useless" — that is exactly the pre-declared falsifier, and it did not fire | §5 R2 (BBB case explicitly excluded, n=1); §7 R4 |
| 8 | "pooled full-minus-non-leaked overstatement +0.0216 [+0.0085, +0.0353]; random control ≈ a tenth of it" | positive (pre-registered S2 direction, MIXED) | row `s2_pooled::overstatement` (+ `ci_low`/`ci_high` fields); arm rows `s2_*::gain_full`, `s2_*::gain_nonleaked` | Percentile bootstrap over 6 curator × seed units — a small basis, stated as such; random-curator control with identical budget and identical number of curation steps | (1) "n=6 is too small for a bootstrap" — conceded in the body and in §Limitations; (2) "3 seeds" — pre-declared floor met, stated | §9 Limitations (loop size) |
| 9 | "every arm loses accuracy on the pre-committed scaffold hold-out (−0.0204 to −0.0482); there were no gains to inflate" | negative (the caveat that bounds claim 8) | rows `s2_random_curator_control::gain_holdout`, `s2_qwen2_5-coder_14b::gain_holdout`, `s2_qwen3_8b::gain_holdout` | Hold-out defined by a deterministic scaffold hash rule, SHA-256 written to `evidence/s2_scaffold_holdout_commit.json` before the first model call | "then your headline is about a broken loop" — yes; the paper claims measurement bias and explicitly refuses the performance claim | §9 Limitations; §6 "There is no gain to inflate"; §1 claim boundary |

## Abstract claims

| # | Claim text (verbatim from abstract) | claim_id | Evidence locator |
|---|-------------------------------------|----------|-------------------|
| 1 | "Our loop runs 3 rounds × 3 seeds × (two local curator models and a random-curator control) over ChEMBL 36 for a solubility predictor" | `s2_pooled::overstatement` | ledger and stamp manifest in `evidence/` |
| 2 | "an exact — not approximate — nearest-neighbour census of 46,276,308,000 Tanimoto comparisons" | `s1_corpus::n_pairwise_comparisons` | ledger and stamp manifest in `evidence/` |
| 3 | "and all 2,854,800 ChEMBL 36 structures" | `s1_corpus::n_chembl_fingerprints` | ledger and stamp manifest in `evidence/` |
| 4 | "audited on a random 5% of queries against two independent implementations at a maximum deviation of 0.0" | `s1_corpus::audit_max_abs_deviation` | ledger and stamp manifest in `evidence/` |
| 5 | "93.84%, 99.52% and 50.78% of the scaffold-split test molecules ... have a near-duplicate (Tanimoto at least 0.9) in the curator's own corpus" (BBB term) | `s1_bbb_martins::leak_rate_0p9` | ledger and stamp manifest in `evidence/` |
| 6 | "93.84%, 99.52% and 50.78% of the scaffold-split test molecules ... have a near-duplicate (Tanimoto at least 0.9) in the curator's own corpus" (lipophilicity term) | `s1_lipophilicity_astrazeneca::leak_rate_0p9` | ledger and stamp manifest in `evidence/` |
| 7 | "93.84%, 99.52% and 50.78% of the scaffold-split test molecules ... have a near-duplicate (Tanimoto at least 0.9) in the curator's own corpus" (solubility term) | `s1_solubility_aqsoldb::leak_rate_0p9` | ledger and stamp manifest in `evidence/` |
| 8 | "against 0.25%, 2.14% and 1.50% in each benchmark's own training split" (BBB term) | `s1_bbb_martins::own_trainval_rate_0p9` | ledger and stamp manifest in `evidence/` |
| 9 | "against 0.25%, 2.14% and 1.50% in each benchmark's own training split" (lipophilicity term) | `s1_lipophilicity_astrazeneca::own_trainval_rate_0p9` | ledger and stamp manifest in `evidence/` |
| 10 | "against 0.25%, 2.14% and 1.50% in each benchmark's own training split" (solubility term) | `s1_solubility_aqsoldb::own_trainval_rate_0p9` | ledger and stamp manifest in `evidence/` |
| 11 | "a ChEMBL-trained forest scores 1.5432 MAE on near-duplicated test molecules" | `s1b_solubility_aqsoldb_A2_chembl_dedup_random_forest::score_leaked` | ledger and stamp manifest in `evidence/` |
| 12 | "against 1.9524 on the rest" | `s1b_solubility_aqsoldb_A2_chembl_dedup_random_forest::score_nonleaked` | ledger and stamp manifest in `evidence/` |
| 13 | "(+0.4093, molecule-level bootstrap 95% CI [+0.2949, +0.5262], three model families)" — point estimate | `s1b_solubility_aqsoldb_A2_chembl_dedup_random_forest::gap` | ledger and stamp manifest in `evidence/` |
| 14 | "(+0.4093, molecule-level bootstrap 95% CI [+0.2949, +0.5262], three model families)" — lower bound | `s1b_solubility_aqsoldb_A2_chembl_dedup_random_forest::gap_mol_ci_low` | ledger and stamp manifest in `evidence/` |
| 15 | "(+0.4093, molecule-level bootstrap 95% CI [+0.2949, +0.5262], three model families)" — upper bound | `s1b_solubility_aqsoldb_A2_chembl_dedup_random_forest::gap_mol_ci_high` | ledger and stamp manifest in `evidence/` |
| 16 | "a near-duplicate carries a logD label 0.0808 log units from the benchmark's own" | `r2_lipophilicity_astrazeneca_vs_chembl_endpoint_pool::near_pair_stat` | ledger and stamp manifest in `evidence/` |
| 17 | "against 1.6181 for random pairs" | `r2_lipophilicity_astrazeneca_vs_chembl_endpoint_pool::random_pair_stat` | ledger and stamp manifest in `evidence/` |
| 18 | "coarsened-exact matching cuts the gap to +0.1595" | `r3_solubility_aqsoldb_A2_chembl_dedup_random_forest::matched_gap` | ledger and stamp manifest in `evidence/` |
| 19 | "while a control that never sees ChEMBL retains +0.0689" | `r3_solubility_aqsoldb_B_tdc_trainval_random_forest::matched_gap` | ledger and stamp manifest in `evidence/` |
| 20 | "the loop's reported improvement exceeds its improvement on the leakage-clean subset by +0.0216 MAE" | `s2_pooled::overstatement` | ledger and stamp manifest in `evidence/` |
| 21 | "yet every arm loses accuracy on the pre-committed hold-out" (control arm term) | `s2_random_curator_control::gain_holdout` | ledger and stamp manifest in `evidence/` |
| 22 | "yet every arm loses accuracy on the pre-committed hold-out" (qwen2.5-coder term) | `s2_qwen2_5-coder_14b::gain_holdout` | ledger and stamp manifest in `evidence/` |
| 23 | "yet every arm loses accuracy on the pre-committed hold-out" (qwen3 term) | `s2_qwen3_8b::gain_holdout` | ledger and stamp manifest in `evidence/` |

Every abstract sentence that carries a number names a manifest key above. The non-numeric abstract sentences ("We implement the loop this workshop's
Track I describes…", "Under a staged pre-registration fixed before any data was
collected…", "so this is a measurement bias in how such loops are scored and not
inflated gains", "We close with five design rules…") are, respectively, a
description of the apparatus in §2, a statement about
`A2_PREDECLARATION_2026-09-07.md` (commit `5121253`), the paper's explicit
refusal of a stronger claim, and a forward reference to §7; none asserts an
unevidenced quantity.

## Rival mechanisms (gate 13b)

The pre-declaration named three rival mechanisms and one discriminator
falsifier. All were measured; none was argued away.

- **R1 — "chemistry overlap is normal" (a leaked molecule is just an easier
  molecule)** → **partly true, quantified, and conceded in the body.** Control
  arm B, trained only on the benchmark's own training split and never on ChEMBL,
  still shows a solubility gap of `+0.2496`
  (`s1b_solubility_aqsoldb_B_tdc_trainval_random_forest::gap`, molecule-level CI
  [+0.1641, +0.338]) against the ChEMBL-trained `+0.4093`. Roughly three fifths
  of the raw gap is generic. Addressed in §5 R1, §7 R3, and §9 Limitations.
- **R2 — endpoint-value mismatch / activity cliffs (near-duplicates carry no
  usable label)** → **the pre-declared falsifier DID NOT FIRE.** Mean |Δ label|
  for near-duplicate pairs vs matched random pairs:
  `r2_lipophilicity_astrazeneca_vs_chembl_endpoint_pool::near_pair_stat` =
  `0.0808` against `random_pair_stat` = `1.6181`, paired difference
  `pair_diff` = `-1.5372` [−1.7484, −1.3339], 157 pairs; and for solubility
  `0.6739` vs `1.5673`, difference `-0.8934` [−1.1391, −0.6478], 124 pairs.
  The duplicates carry the answer. The BBB-Martins row
  (`r2_bbb_martins_vs_own_trainval::*`, near `1.0` vs random `0.0`) rests on
  **one** pair and is explicitly reported as a count, not an estimate.
  Addressed in §5 R2 and §7 R4.
- **R3 — subset-composition confound (leaked and non-leaked subsets differ in
  covariates)** → **confirmed, and it explains most of the effect.** SMDs before
  matching: `r3_solubility_aqsoldb_mw::smd` = `-0.2529`,
  `r3_solubility_aqsoldb_chembl_scaffold_count::smd` = `-0.0441`. Coarsened-exact
  matching on MW quintile × label quintile (796 molecules per side) shrinks the
  primary gap from `+0.4093` to
  `r3_solubility_aqsoldb_A2_chembl_dedup_random_forest::matched_gap` = `+0.1595`
  [+0.1327, +0.189], with the control arm's matched gap at
  `r3_solubility_aqsoldb_B_tdc_trainval_random_forest::matched_gap` = `+0.0689`
  [+0.0393, +0.0941]. A small ChEMBL-training-specific component survives.
  Addressed in §5 R3 and §7 R3, conceded in §9 Limitations.
- **Discriminator falsifier (pre-declared): "if near-duplicate pairs show label
  agreement no better than random pairs at matched property distributions, the
  leakage mechanism is falsified"** → **not met.** See R2 above. Reported in the
  abstract, in §5 R2, and in §7 R4.
- **Refutation criteria for the three stages** → S1 **held** (bar was 20%,
  smallest lower Wilson bound 48.58%); S1b **mixed** (held on
  Solubility-AqSolDB, unresolved on BBB-Martins with n=25, not estimable on
  Lipophilicity-AZ with n=4 — the last two reported as such and never as nulls);
  S2 **mixed** (direction held with a CI excluding zero, but no arm improved on
  the pre-committed hold-out, so the paper claims measurement bias and
  explicitly not inflated gains).

## Per-citation aptness (load-bearing citations only)

| # | Citekey | Claim clause (verbatim) | Aptness line |
|---|--------|-------------------------|--------------|
| 1 | `huang2021tdc` | "built TDC on scaffold splits precisely to suppress train--test structural similarity" | supports the clause because Huang et al. introduce the Therapeutics Data Commons and specify scaffold splitting as the ADMET Benchmark Group's default split strategy, which is the design our own-split rates (0.25–2.14%) confirm is working. |
| 2 | `schuh2026auditing` | "audit widely used biomolecular benchmarks and report systematic annotation-level inconsistencies --- duplicated and conflicting labels across source datasets" | supports the clause because Schuh, Daniluk and Sieber's audit reports systematic data inconsistencies at the level of records and their annotations across widely used biomolecular benchmarks; this is the paper we must distinguish ourselves from by name, and the distinguishing feature is their unit of analysis (a record's label) versus ours (a structural neighbour in an external corpus). |
| 3 | `vantilborg2022activity` | "expose molecular ML with activity cliffs --- structurally similar, experimentally divergent pairs within a dataset" | supports the clause because van Tilborg et al. define activity cliffs as structurally similar molecule pairs with large activity differences and show models fail on them; our R2 measures the same pair-level quantity across corpora and finds the opposite regime (mean |Δ logD| = 0.0808) dominant. |
| 4 | `faraji2026leakage` | "run a leakage-aware benchmark study for deep-eutectic-solvent property prediction" | supports the clause because Faraji et al.'s ACS Omega study explicitly re-benchmarks ML models for DES property prediction under leakage-aware splits, establishing the methodology while operating in a different chemical domain with no external-corpus census. |
| 5 | `golbraikh2002beware` | "warned that internal validation statistics overstate QSAR performance when the test set sits inside the training set's applicability domain" | supports the clause because Golbraikh and Tropsha's paper is the canonical demonstration that leave-one-out q² does not predict external predictive power, i.e. that internal statistics overstate performance without chemical-space control. |
| 6 | `palavalli2024taxonomy` | "give a taxonomy of data contamination for text benchmarks; we instantiate one of its modes in a numeric molecular benchmark" | supports the clause because Palavalli, Bertsch and Gormley enumerate contamination modes for LLM evaluation (including derived/analogue forms beyond exact string overlap), which is the taxonomy slot our structural-analogue leakage occupies. |
| 7 | `shumailov2024collapse` | "show that models trained on recursively generated data collapse" | supports the clause because Shumailov et al. demonstrate degeneration across generations when training data is model-generated; our loop's degradation occurs with entirely real ChEMBL measurements, which is what makes it a distinct mechanism. |
| 8 | `gao2022scaling` | "the same structure as reward-model overoptimisation" | supports the clause because Gao, Schulman and Hilton show true reward declining while proxy reward rises as optimisation pressure increases, the same proxy/true divergence our full-test versus hold-out trajectories show. |
| 9 | `manheim2018goodhart` | "Framed as an optimisation problem, the loop is a Goodhart instance" | supports the clause because Manheim and Garrabrant's taxonomy of Goodhart variants includes regressional and extremal cases in which optimising a proxy selects for the proxy's error term, which is what benchmark-feedback curation does here. |
| 10 | `bajusz2015tanimoto` | "the similarity coefficient and its use as a clustering threshold" | supports the clause because Bajusz, Rácz and Héberger compare similarity coefficients on binary fingerprints and conclude Tanimoto is the appropriate default, which is our justification for using it as the leakage metric. |
| 11 | `rogers2010ecfp` | "Morgan fingerprints, radius 2, 2048 bits" | supports the clause because Rogers and Hahn define extended-connectivity fingerprints and the radius/bit-length parameterisation we use (ECFP4 = Morgan radius 2). |
| 12 | `bemis1996molecular` | "Scaffolds are Bemis--Murcko" | supports the clause because Bemis and Murcko define the molecular framework decomposition that both TDC's scaffold split and our scaffold-level leak rate and Stage 2 hold-out rule use. |
| 13 | `delong1988comparing` | "the correlated-ROC machinery of DeLong et al. does not apply and we bootstrap instead" | supports the clause because DeLong, DeLong and Clarke-Pearson's test is derived for ROC curves computed on the *same* sample, whereas our leaked and non-leaked subsets are disjoint molecule sets, so their variance estimator's correlation structure is absent. |
| 14 | `gaulton2017chembl2017` | "ChEMBL is the corpus we census" | supports the clause because Gaulton et al. describe the ChEMBL database release line and its bioactivity content, which is the corpus (release 36) our census enumerates. |
| 15 | `wishart2017drugbank` | "one of several a curator could use" | supports the clause because DrugBank 5.0 is an independently curated public compound/bioactivity resource of comparable accessibility, so a curator restricted to neither ChEMBL nor DrugBank alone would face at least the exposure we measure. |
| 16 | `ott2022mapping` | "document benchmark saturation across AI at the level of publication records without proposing a mechanism" | supports the clause because Ott et al.'s analysis is bibliometric — mapping benchmark creation and saturation dynamics from publication data — and does not identify corpus overlap or any other data-level mechanism. |
| 17 | `wu2018moleculenet` | "the tasks themselves descend from MoleculeNet" | supports the clause because MoleculeNet assembled the original molecular ML benchmark suite (including the solubility and BBB task lineage) that TDC later restated. |
| 18 | `boiko2023autonomous` | "give language models tools and let them plan and execute chemistry" | supports the clause because Boiko et al. demonstrate an LLM-driven system that plans and executes real chemical syntheses using search and instrument tools. |
| 19 | `bai2022constitutional` | "feedback from the model's own criticism" | supports the clause because Constitutional AI trains on model-generated critiques and revisions rather than human labels, which is the self-evaluation ingredient of the RSI template we instrument. |
| 20 | `zheng2023judging` | "model-generated evaluation" | supports the clause because Zheng et al. establish LLM-as-a-judge as an evaluation mechanism and characterise its biases, which is the automated-evaluation ingredient of the loop template. |
| 21 | `wang2023selfinstruct` | "self-generated training data" | supports the clause because Self-Instruct generates the model's own instruction-tuning data from the model itself, the data-generation ingredient of the loop template. |

## Unexported interval bounds (written as `\SafeNum`, quoted from ledger rows)

These are real measured bounds stored in the `ci_low` / `ci_high` fields of the
ledger rows named beside them; `analysis/make_ledger.py` does not export them as
separate manifest keys, so `\StampedCell` cannot address them. **If the
the analysis regenerates the manifest with `<metric>_ci_low` / `_ci_high` keys
for these rows, they should be converted to `\StampedCell`.**

- Wilson bounds for `s1_*::leak_rate_0p9` (91.07/95.79, 98.78/99.81, 48.58/52.97)
- Bootstrap bounds for `r3_*::matched_gap` (all seven rows in Table 3, plus the BBB row)
- Paired-bootstrap bounds for `r2_*::pair_diff`
- Bootstrap bounds for `s2_pooled::overstatement` (0.0085 / 0.0353)
- Per-arm S2 `full − non-leaked` differences (0.0021 / 0.0182 / 0.025), recorded
  per arm as `mean_overstatement_full_minus_nonleaked` in
  `evidence/s2_rsi_results.json`; the paper states in the Table 4 caption that
  this column is the arithmetic difference of the two stamped columns beside it.
- `0.4079` (§5 R3), the molecular-weight-only matched gap for
  `A2_chembl_dedup|random_forest`, quoted from
  `evidence/s1b_rival_mechanisms.json` →
  `benchmarks.solubility_aqsoldb.R3_matched_gap["A2_chembl_dedup|random_forest"].matched_on_mw_only.matched_gap_mean`.
  Added 2026-09-07 by the review round to disclose that the matched shrinkage is
  driven by the label quintile and not by molecular weight (DOMAIN_AUDIT row 25).
- `157` and `124` (§5 R2), the near-duplicate pair counts against the ChEMBL
  endpoint pools, from the same file's
  `R2_label_agreement.vs_chembl_endpoint_pool.n_near_pairs`. Their denominators
  (`s1b_*::n_leaked` = 836 and 1,013) are stamped cells and are now printed
  beside them (DOMAIN_AUDIT row 22).

## Fourth analysis round, 2026-09-07 (C8 remediation ladder + interval-degeneracy audit)

**This registry was written before the rows below were added; `evidence/agenticls-rsi-leakage_ledger.jsonl`
is the authority.** The section that precedes this one reports 309 ledger rows and 60 `must_print`;
both counts are **superseded** — the shipped ledger and stamp manifest carry **322 rows, 65
`must_print`**, and the guard resolves 276 stamped cells against them (`build.guard.log`). The 13
rows added after the third round are transcribed below directly from the ledger; no verdict is
restated or reinterpreted here.

### Rows ADDED (values quoted from the ledger)

- **C8, does the R6 remedy work? — pre-registered, and it did not reproduce.** The design,
  primary quantity and the reporting rule for every outcome were fixed in
  `C8_PRECOMMITMENT_2026-09-07.md` (SHA-256 `79d502ef…`, recorded in
  `evidence/c8_remediation.json`) before `analysis/c8_remediation.py` was written. The loop was to
  be re-scored on progressively stricter decontaminations (`clean@0.9/ChEMBL`,
  `clean@0.9/union`, `clean@0.8/union`, `clean@0.7/union`) from the recorded curator selections,
  with no new model calls. The pre-committed reproduction check failed:
  `c8::holdout_repro_deviation` = **0.0062** against a tolerance of **0.0005**, so
  **outcome D of the pre-commitment fired and no re-scored number is published**. The full
  per-unit table is shipped in `evidence/c8_remediation.json` for audit.
- **Interval-degeneracy audit.** `ci_audit::n_intervals` = **111** intervals reaching the PDF or
  the ledger were resolved to their resampling unit and unit count by
  `analysis/ci_degeneracy_audit.py`; `ci_audit::n_degenerate` = **3** have percentile bounds that
  coincide with the [min, max] of their units or have zero width
  (`r2_bbb_martins_vs_own_trainval::near_pair_stat`, `::pair_diff`, `::random_pair_stat`).
- **Paired-t intervals beside the percentile bootstraps.** A percentile bootstrap is degenerate at
  very small n, so a paired-t interval is reported alongside the pooled bootstrap for the same
  units. `s2_pooled::overstatement_t_ci_low/high` = **[+0.0021, +0.0411]** over 6 curator x seed
  units; the five C7 paired quantities get the same treatment —
  `c7_paired::gain_scaffold_holdout_t` **[−0.0236, +0.0370]**,
  `gain_full_test_t` **[−0.0155, +0.0196]**,
  `gain_nonleaked_test_t` **[−0.0426, +0.0308]**,
  `overstatement_t` **[−0.0183, +0.0342]**.

## Third analysis round, 2026-09-07 (experiments C7 and C6 + design rule R6)

Ledger re-diffed against the previous revision: **205 pre-existing rows unchanged, 0 removed,
104 added** (309 rows, 60 `must_print`). Two new scripts, both reusing artefacts on
disk: `analysis/s2_shuffled_feedback.py` (C7) and `analysis/c6_drugbank_census.py` (C6).

### Claims ADDED

- **C7, feedback-shuffled control arm — a pre-registered NULL.** Design, primary
  quantities and the reporting rule for every outcome were fixed in
  `C7_PRECOMMITMENT_2026-09-07.md` and hashed (`f885a47c…`) before the implementation
  was written. Paired real-minus-shuffled: hold-out
  `c7_paired::gain_scaffold_holdout` = +0.0067 [−0.0132, +0.0297]; full test +0.0020;
  clean subset −0.0059; overstatement +0.0080 — **every interval spans zero**.
  Selection agreement `c7_selection::jaccard` = 0.2494 [0.1593, 0.3460]. 18 calls,
  0 fallbacks. **Claim: benchmark feedback changed what this curator selected but not
  where the loop ended.** In the abstract, per the pre-committed rule.
- **C6, DrugBank second-corpus census.** `c6_drugbank::corpus_size` = 12,293;
  `c6_*::union_rate_0p9` = 94.58 / 99.76 / 51.23%; clean subsets 982→973, 25→22, 4→2.
  **Claim: our rates are a lower bound, measured rather than asserted, and
  "leakage-clean" is corpus-relative — but the increments are small, so the headline
  is not an artefact of one corpus choice.**
- **Design rule R6 (Decision 10).** "Fix the split, not the curator." Every clause
  restates a stamped row; no extrapolation. "Five design rules" → six in the abstract,
  the contribution list and the conclusion.

### Pre-commitment deviation, disclosed

The C7 pre-commitment's outcome-A row was a conjunction — interval includes zero *and*
selection agreement is high. **The performance half held; the agreement half did not**
(Jaccard 0.2494 is low). The paper therefore states the sharper claim the data supports
rather than the weaker one row A anticipated, and the failed conjunct is recorded in §7
of the pre-commitment file with both the pre-run and post-run hashes published.

### Claims unchanged

Both claims weakened in round one **stay weakened**: no arm-versus-control separation
(p = 0.131), and no non-overlap inference. The raw/matched **bracket stays a bracket** —
+0.1595 lower bound, +0.4093 upper bound, neither claimed as the estimate.

### Structural change

The §7 closing cost paragraph and the §3 C6 paragraph were removed to pay for R6, the
C7 paragraph and the abstract revision within the 9-page cap; the C6 result moved into
§9 Limitations, where the hedge it corrects already lived. All 60 `must_print` rows are
printed (guard verifies on every build).

## Second analysis round, 2026-09-07 (experiments C1, C2, C5 + matching ladder)

Three authorised analyses and one required reframing were run. **No pre-existing
ledger row changed and none was removed** — verified by diffing the regenerated
ledger against `git HEAD`: 205 rows in, 205 unchanged, **69 added**, 0 removed
(`must_print` 45 → 54). New numbers come from two scripts that reuse artefacts
already on disk:

- `analysis/review_round.py` — C1 endpoint-pool census, C2 nested bootstrap and
  paired A2−B difference, C5 cluster-level leak readout. Writes
  `evidence/review_round_2026_09_07.json`.
- `analysis/matching_ladder.py` — the missing label-only rung. It calls
  `s1b_rivals.cem_matched_gap` itself with the same seed-key convention, and
  **asserts that both published rungs reproduce exactly** before emitting the new
  one, so a drift would fail loudly rather than reshape a table.

### Claims ADDED

- **Endpoint-pool exposure (C1).** `s1_*::endpoint_pool_rate_0p9` = 18.69% (Lipo,
  Wilson [16.2, 21.47]) and 6.22% (Sol, [5.24, 7.36]) against structural rates of
  99.52% and 50.78%. Table 1 now carries both. This is the measured answer to what
  the first round could only concede in prose.
- **Matching ladder.** `r3_*::matched_gap_mw_only` and
  `::matched_gap_label_only` for all nine cells. Primary arm: raw +0.4093 →
  MW-only +0.4079 → label-only +0.1440 → MW×label +0.1595. Molecular weight does
  essentially nothing; the label quintile does all of it.
- **Nested bootstrap (C2).** `r3_*::matched_gap_nested` with `_ci_low`/`_ci_high`,
  1000 outer draws over molecules with re-stratification and re-matching inside
  each. Replaces the redraw-only interval in the table.
- **Paired ChEMBL-specific residual (C2).** `r3_solubility_aqsoldb_a2_minus_b_*::
  paired_difference`: +0.0811 [+0.0104, +0.1567] (forest), +0.2749 (linear),
  +0.1778 (kNN), all excluding zero.
- **Cluster leak readout (C5).** `s2_cluster_readout::*` — 16 near-duplicates in a
  27,230-compound universe; LLM selection-weighted reproduction rate 0.0006 against
  a universe mean of 0.0006; Spearman −0.3201.

### Claims RESTORED

- **A ChEMBL-training-specific residual bounded away from zero.** Withdrawn in the
  first round because non-overlapping marginal intervals are not a difference test.
  Now computed directly and paired over shared molecules, and restored with its
  lower bound stated as a fortieth of the raw gap.

### Claims WITHDRAWN or WEAKENED

- **The S2 steering mechanism — REFUTED by our own measurement.** §6 previously
  asserted that a curator following benchmark feedback drifts toward regions where
  the test set is reproduced. C5 measured it and it does not hold: the curation
  universe contains almost no reproducible test chemistry, the curators showed no
  selection preference, and the feedback signal correlates *negatively* with
  reproduction. The paragraph is now a reported negative, and the surviving reading
  is that the leaked and clean halves of the split are not exchangeable so any pool
  change separates the two numbers. The abstract's claim ("a measurement bias in how
  such loops are scored") was already the claim this supports and did not change in
  substance. Recorded in §9 Limitations as well.
- **Control arm B's matched gap.** Under the nested bootstrap its interval
  (+0.0711 [−0.0049, +0.1428]) no longer excludes zero, where the redraw-only
  interval had. Reported.
- The two claims weakened in the first round — arm-versus-control separation
  (p = 0.131) and the raw/matched bracket — **stay weakened**.

### Structural change

`Table 2` (per-cell MAE leaked/clean) was **removed**; its gap column is reproduced
in the ladder table's new raw-gap column, and the primary arm's MAE pair is in §4
prose. This paid for the ladder within the 9-page cap. Every `must_print` row it
carried is still printed (guard verifies).

## Review round, 2026-09-07 (post-RSI-reframe)

An adversarial review plus an independent verification pass against `analysis/`
produced ten fixes, audited in `DOMAIN_AUDIT.md` rows 22–37.
**No ledger row changed and no number was recomputed**; every fix is a prose
correction, a disclosure of a value that was already in the evidence tree, or a
withdrawal of an inference the evidence did not support. Claims materially
weakened by this round, recorded so the change is not silent:

- **Arm-versus-control separation (§6, §1) — WITHDRAWN.** The pooled LLM
  interval still excludes zero, but the paper no longer asserts that the bias
  belongs to LLM curation rather than to reshuffling: an exact permutation test
  over the 6 LLM and 3 control units gives p = 0.131.
- **"ChEMBL-training-specific component bounded away from zero" (§5 R3) —
  WITHDRAWN.** Non-overlap of two marginal intervals is not a difference test,
  and no interval for `G_A2,matched − G_B,matched` exists in the ledger.
- **Matched-interval semantics (Table 3, §5 R3) — RELABELLED.** The bounds are
  quantiles over redraws of the within-stratum down-sample
  (`cem_matched_gap`, `n_rep=200`), i.e. matching variability, not a molecule
  bootstrap, and are no longer presented as comparable to Table 2.
- **kNN "upper bound" (§8) and hold-out "told the truth" (§7 R5) —
  DOWNGRADED** to what the design actually supports.
- **Architecture invariance (§2, §8) — NARROWED** from the gap magnitude to the
  exposure, since the gap varies ~1.5× across our own three families.

## Checklist

- [x] Every row's `claim_id` resolves in the ledger with `verdict=SUPPORTED`
- [x] No bare prose percentage or ratio in the PDF lacks a row here (guard gate 3, strict, E7 clean)
- [x] CI includes-zero status noted for every CI (BBB and lipophilicity rows: spans zero, reported as unresolved / not estimable)
- [x] Reviewer-objection risk column filled
- [x] Limitation sentence location cited for every row
- [x] Per-citation aptness filled for every load-bearing citation
