# C8 pre-commitment — does the R6 remedy actually work?

**Written and hashed BEFORE the analysis was implemented or run**, in the same discipline
as `A2_PREDECLARATION_2026-09-07.md` and `C7_PRECOMMITMENT_2026-09-07.md`. It fixes the
design, the primary quantity and the reporting rule **for every outcome, including the one
that undercuts our own design rule R6**.

Fixed 2026-09-07, following the C6 and C7 results.

---

## 1. Why this analysis exists

Rules R1, R2 and R6 tell a practitioner to fix the *evaluation split*: census the corpus
the loop can reach, construct a leakage-clean subset, and score the loop on that. The
paper currently shows this is *necessary* — the full-benchmark number and the clean-subset
number diverge by +0.0216 MAE. It does not show the remedy is *sufficient*, or even that
it is stable. Two things could be wrong with it:

1. **The clean subset may not be clean.** It is defined against ChEMBL only. C6 measured a
   second corpus and found the solubility clean subset drops from 982 to 973 once DrugBank
   is counted. If those residual molecules were carrying the effect, the divergence we
   report is an artefact of an under-specified "clean".
2. **The threshold may be doing the work.** "Clean" means nearest-neighbour Tanimoto
   below 0.9. A more conservative cut would remove far more molecules. If the divergence
   evaporates under a stricter cut, the remedy is threshold-dependent and R6 is weaker
   than stated.

C8 tests both by re-scoring the loop on progressively stricter decontaminations. **No new
model calls are needed:** every curator selection is already logged in
`evidence/s2_model_calls.jsonl` and `evidence/s2_rsi_results.json`, so the arms are
re-fitted and re-scored from recorded decisions, not re-run.

## 2. Design

* Universe, clustering, hold-out, budget, trees, seeds: reproduced exactly as in S2, with
  the pre-committed hold-out SHA-256 re-asserted at run time (abort on mismatch).
* Each arm × seed × round is refitted on **the pool implied by the recorded cluster
  selection**, using the same `zlib.crc32` seed key, so the pools are byte-identical to the
  original run. Round 0 is the original fixed-seed random draw.
* Each fitted model is scored on **five** sets:
  1. the full TDC test split;
  2. `clean@0.9 / ChEMBL` — the paper's current leakage-clean subset;
  3. `clean@0.9 / ChEMBL ∪ DrugBank` — decontaminated against both accessible corpora;
  4. `clean@0.8 / union` and 5. `clean@0.7 / union` — progressively more conservative.
* DrugBank nearest neighbours for the solubility test split are recomputed with the same
  fingerprint and the same exact search used everywhere else.

## 3. Primary quantity, fixed in advance

**Overstatement at each decontamination level** = (round-3 minus round-0 gain on the full
test set) − (the same gain on that clean subset), pooled over the **6 LLM curator × seed
units**, percentile bootstrap, reported as directional evidence rather than as a calibrated
interval (6 units cannot support nominal coverage — as everywhere else in this paper).

Secondary, and reported regardless: **each clean subset's size**, and the **hold-out
trajectory**, which must remain unchanged because it does not depend on the test split at
all. If the hold-out numbers move, the reproduction is broken and the run is void.

## 4. Reporting rule — fixed before the run, for every outcome

| Outcome | What we will report, and where |
|---|---|
| **A — the remedy is stable.** The overstatement still excludes zero at the union-clean level and is close to the ChEMBL-clean value. | §6 and the abstract: scoring on a subset decontaminated against *both* accessible corpora gives the same corrected picture, so the divergence is not an artefact of an under-specified "clean". R6 stands, and is now a **tested** intervention rather than a diagnosis. |
| **B — the effect disappears.** The overstatement interval covers zero once the union-clean subset is used. | **This undercuts R6 and will be reported exactly that plainly**, in §6, in the abstract, and in the Limitations, stated as: our reported divergence was an artefact of defining "clean" against one corpus, and the design rule that follows from it is correspondingly weaker. The rule would be rewritten or withdrawn in the same revision, not defended. |
| **C — the effect grows** materially at stricter thresholds. | Reported in §6 as-is, with the note that a conservative threshold makes the bias *more* visible, and with the shrinking subset size stated beside it so the reader can see the power cost. |
| **D — the reproduction fails**: hold-out gains do not match the original S2 run to printing precision. | The analysis is void and is reported as void. No re-scored number is published. Nothing is "fixed" to make it agree. |

**Threshold ladder.** Whatever the outcome, all five subset sizes and all four
overstatement values are published together as a ladder, so a reader can see where the
estimate becomes unstable rather than being shown only the level that suits us. Where a
subset falls below `MIN_SUBSET_N = 20` — the threshold fixed in code before the S1b run —
the estimate is reported as underpowered and no verdict is drawn from it.

**Prediction, recorded so it can be scored against.** My last two acceptance-effect
predictions were both too high, in the same direction. I expect outcome A, because C6
measured only 9 molecules of difference at the 0.9 level. That expectation is written here
so that if B or C occurs it is visibly a surprise and not retrofitted. **The outcome will
be reported, not the prediction.**

## 5. What this analysis cannot show

It re-scores the same six units on stricter splits; it does not add seeds, curators,
rounds or benchmarks, and it cannot make a null about this loop into a general claim. It
also cannot show that a *decontaminated* loop would improve — no arm improved on the
hold-out, and nothing here changes that.

## 6. Provenance

* Implementation: `analysis/c8_remediation.py`, written after this file.
* Output: `evidence/c8_remediation.json`.
* Ledger rows via `analysis/make_ledger.py`; no existing row is modified.
* This file's SHA-256 is recorded in the run output.
