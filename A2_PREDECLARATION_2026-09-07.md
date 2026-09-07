# A2 Pre-declaration — agenticls-rsi-leakage

**Date:** 2026-09-07 (UTC) · **Venue:** AgenticLS — Agentic AI for Biological Discovery, NeurIPS 2026
(`NeurIPS.cc/2026/Workshop/AgenticLS`; site deadline 2026-09-18 AoE, OpenReview probe 2026-09-19
19:00 AoE — the earlier date is treated as binding) · **Blind:** double.
**Track:** Track I (Building Agentic Systems for Life Science) — the CFP lists verbatim
"Recursive self-improvement for bio-specialized models: using general-purpose LLMs to iteratively
improve specialized models through data generation/curation, automated evaluation, and
feedback-driven refinement" — with the evaluation contribution also serving Track II.
**Target length:** staged (see below): Stage 1 alone supports the 4-page short track; Stage 1+2
support the 9-page full track, which is the target.

**Provenance:** a new hypothesis, fixed 2026-09-07, extending a leakage-conditioned RSI
measurement; not a reframe of any prior paper. No row depends on any pre-existing local
archive: all evidence comes from live databases and new code written for this study.

## Pre-registered claims

**Stage 1 — leakage census (S1).** For each of ≥3 TDC ADMET benchmarks (BBB-Martins,
Lipophilicity-AstraZeneca, Solubility-AqSolDB), **more than 20%** of held-out *test* molecules have
a near-duplicate (Morgan radius 2, 2048 bits, Tanimoto ≥ 0.9) somewhere in ChEMBL 36 — the public
corpus an LLM curator can freely draw from when generating training data for a bio-specialized
model.

**Stage 1b — does leakage matter (S1b).** A standard ADMET predictor (random forest on ECFP4,
plus a second model family) trained on ChEMBL-derived data scores **higher on the leaked test
subset than on the non-leaked subset**, with the AUC/RMSE gap's 95% CI excluding zero.

**Stage 2 — leakage-conditioned RSI gains (S2).** In a 3-round LLM-curator recursive
self-improvement loop (a general-purpose LLM iteratively selects/curates ChEMBL training data for
the specialized ADMET model, guided by benchmark feedback), the **apparent** round-3 improvement
measured on the full TDC test set overstates the improvement measured on the **non-leaked subset**
by a margin whose paired 95% CI excludes zero.

## Design (fixed before any run)

- **Data (real, live):** ChEMBL 36 local SQLite (`chembl_36.db`, 27.7 GB — verified on disk:
  `compound_structures`, `molecule_dictionary`, `activities`); TDC ADMET benchmark files (verified
  on disk: 2039 / 4200 / 9982 rows) re-verified against the TDC package/release with SHA256 logged.
  RDKit for fingerprints and canonicalisation. Every file's path, row count and SHA256 recorded in
  `evidence/DESIGN_NOTES.md` before analysis.
- **S1 method:** canonicalise all SMILES; Morgan r=2, 2048-bit fingerprints; for each test molecule,
  nearest-neighbour Tanimoto against the ChEMBL corpus (blocked/indexed search — exact NN, not
  approximate, or approximate with an exactness audit on a random 5% subsample). Report per
  benchmark: exact-duplicate rate, ≥0.9 rate, ≥0.95 rate, and the full similarity distribution.
  Also report the same rates against the benchmark's *own training* split as a reference scale.
- **S1b method:** train predictors on ChEMBL-derived data for the same endpoint (never on TDC test),
  evaluate separately on leaked (NN ≥ 0.9) and non-leaked (NN < 0.9) test subsets; ≥5 seeds;
  paired bootstrap CIs on the gap; report subset sizes and endpoint-distribution differences.
- **S2 method:** RSI loop, 3 rounds, ≥3 seeds, 2 curator LLMs (local Ollama models and/or an API
  model — recorded before the run). Each round the curator proposes which ChEMBL records to add or
  drop from the training pool, given the previous round's benchmark feedback; the specialized model
  is retrained; the loop reports its gain. We measure gains on (i) the full test set, (ii) the
  non-leaked subset, (iii) a pre-committed scaffold-split held-out set. Identical curation budget
  per round across arms.
- **Unit of independence:** molecule for S1 (with scaffold-level clustering reported), seed for
  S1b/S2. All numbers enter the paper only via ledger rows → stamped cells.

## Refutation criteria / honest-negative fallbacks

- S1 refuted if the ≥0.9 rate is < 20% for a benchmark → report the actual leakage curve per
  benchmark; a measured low-leakage result is still the first such census for the RSI-curation
  setting and is published as a defence of the benchmarks.
- S1b refuted if the leaked/non-leaked performance gap CI includes zero → report that near-duplicate
  presence does not inflate this endpoint's scores (an important negative for the census's
  interpretation, and reported as such rather than buried).
- S2 refuted if the apparent-vs-non-leaked gain gap CI includes zero → report that LLM-curated RSI
  gains on these benchmarks are not leakage-driven.
- If S2 cannot complete before the deadline, the paper ships as the 4-page short version with S1+S1b
  only, and says so explicitly; S2 is never reported as complete when it is not.

## Rival mechanisms (gate 13b)

1. **Chemistry-overlap-is-normal** — any test set inherits similar chemistry; near-duplicates need
   not inflate scores. Directly tested by S1b (the leaked/non-leaked gap) rather than assumed.
2. **Endpoint-value mismatch** — a structural near-duplicate may carry a different measured value
   (activity cliffs), so leakage need not be informative; we report the *label agreement* of
   near-duplicate pairs alongside the structural rate.
3. **Subset-composition confound** — leaked and non-leaked subsets may differ in molecular weight,
   scaffold frequency, or label balance; we report covariate distributions for both subsets and a
   covariate-matched re-estimate of the gap.

**Discriminator/falsifier:** if near-duplicate pairs show label agreement no better than random
pairs at matched property distributions, the leakage mechanism is falsified — the duplicates carry
no exploitable signal — and the paper reports that.

## Feasibility bound

S1 is the guaranteed deliverable (fingerprinting + NN search over a 27.7 GB SQLite corpus; 1–2 days).
A 60-second foreground pilot (one benchmark × 100 test molecules × a ChEMBL subset) must succeed
before the full census. S1b and S2 proceed only if S1 completes with ≥2 days of runway left.
