# C7 pre-commitment — feedback-shuffled control arm

**Written and hashed BEFORE the arm was implemented or run.** This file fixes the design,
the primary quantities, and — critically — **the reporting rule for every outcome**, so that
what the paper says about C7 cannot be chosen after seeing the result. It follows the same
discipline as `A2_PREDECLARATION_2026-09-07.md`, which fixed the S1/S1b/S2 design before any
data was collected.

Fixed 2026-09-07, following the C5 result.

---

## 1. Why this arm exists

C5 measured the mechanism the paper had asserted for its S2 overstatement and **refuted it**:
only 16 of 27,230 curation-universe compounds are Tanimoto ≥ 0.9 near-duplicates of a test
molecule; the clusters the LLM curators selected carry a selection-weighted test-reproduction
rate of 0.0006 against a universe mean of 0.0006 (no preference); and the rank correlation
between a cluster's reproduction rate and the per-cluster benchmark error the curator was
shown is **negative** (ρ = −0.3201, consistent across three round-0 seeds).

That leaves an obvious question, and a reviewer will ask it: **does the benchmark feedback do
anything at all?** If the curator's selections and the loop's trajectory are unchanged when
the feedback column is replaced by noise, then "feedback-driven refinement" — the third of the
three things Track I names — is not operating in this loop, and that is a finding about
LLM-curated RSI loops rather than a footnote about our implementation.

## 2. Design

One new arm per existing LLM curator, identical to the real-feedback arm in every respect
except the feedback column.

* **Curators:** `qwen2.5-coder:14b` and `qwen3:8b`, temperature 0.0, thinking off,
  `num_predict=800` — unchanged from S2.
* **Seeds:** 3 per curator, using the same seed keys as S2, so round 0 is the *identical*
  fixed-seed random draw and the two arms start from the same model and the same pool.
* **Rounds:** 3. **Budget:** 6,000 compounds. **Clusters:** 40. **Specialized model:**
  random forest, 200 trees, `max_features='sqrt'`, retrained from scratch each round.
* **Universe:** the same ChEMBL solubility pool minus the pre-committed scaffold hold-out
  and minus exact TDC-test duplicates; the same MiniBatchKMeans clustering
  (`random_state=0`). The hold-out SHA-256 is unchanged and is re-asserted at run time.
* **The single manipulation:** the per-cluster `benchmark MAE on the benchmark molecules
  closest to this cluster` column is **randomly permuted across clusters** before the prompt
  is built. Every other field the curator sees — cluster id, size, label mean, label spread,
  example SMILES — is byte-identical to the real-feedback arm, and the permutation is drawn
  from a `zlib.crc32`-derived seed so it is reproducible. Cluster ids are **not** permuted,
  so the curator's selection still addresses real clusters and the drawn pool is real
  chemistry; only the *signal* about which clusters need work is destroyed.
* Every call is logged (prompt, prompt SHA-256, raw response, latency, parse outcome) into
  the same `s2_model_calls.jsonl` schema. Fallbacks to "select everything" are counted.

## 3. Primary quantities, fixed in advance

For each curator × seed unit, round 3 minus round 0 (positive = MAE improved), scored on the
three sets S2 already uses:

1. `gain_scaffold_holdout` — **the primary outcome.** It is the only scored set that is not
   the benchmark, and S2 already established it is the one that moved.
2. `gain_full_test` and `gain_nonleaked_test`, and their difference (the overstatement).
3. **Selection agreement:** mean Jaccard similarity between the cluster set chosen under real
   feedback and the set chosen under shuffled feedback, matched on curator, seed and round.
   This is the direct behavioural read: a curator that is using the feedback column should
   select differently when it is noise.

Contrasts are **paired** on (curator, seed), because the two arms share round 0 exactly.
Interval: percentile bootstrap over the 6 paired units, reported — as everywhere else in this
paper — as directional evidence and not as a calibrated interval, since 6 units cannot
support nominal coverage.

## 4. Reporting rule — fixed before the run, for every outcome

The rule is fixed in advance and adopted verbatim: *if shuffled
feedback performs indistinguishably from real feedback, that is a finding about LLM-curated
RSI loops that belongs in the abstract, not a footnote.*

| Outcome | What we will report, and where |
|---|---|
| **A — indistinguishable.** The paired real-minus-shuffled interval on hold-out gain includes zero, and selection agreement is high. | **Abstract, plus a design rule.** Stated as: benchmark feedback did not measurably change what the curator selected or where the loop ended, so the "feedback-driven refinement" component of this loop was inert. Reported as a property of *this* loop at *this* scale, not as a general claim about all curator LLMs. |
| **B — real feedback is better** on the hold-out (interval excludes zero, favouring real). | Abstract and §6: the feedback carries usable signal, which strengthens the loop's status as a genuine Track I instance and makes the S2 overstatement a more serious finding. |
| **C — real feedback is worse** on the hold-out (interval excludes zero, favouring shuffled). | Abstract and §6, reported plainly: following benchmark feedback actively harmed the model relative to noise. This would be the most striking outcome and it will not be softened. |
| **D — parse failures or fallbacks** in the shuffled arm materially exceed the real arm's zero. | The arm is reported as compromised and the comparison is withheld; the failure count is published either way. |

In **all** outcomes the arm is reported. It will not be dropped, re-run with different
settings, or relegated to an appendix if the result is inconvenient. The selection-agreement
statistic is reported in all outcomes regardless of the hold-out result, because it is the
cleanest behavioural measurement and does not depend on the noisy downstream gain.

**Bar for the abstract.** Outcome A enters the abstract as a stated null. Outcomes B and C
enter the abstract as stated effects. No outcome is reported as a trend or a suggestion.

## 5. What this arm cannot show

Six paired units, two small local curators, one benchmark, one specialized-model family. A
null here is a null about *this* loop at *this* scale — a frontier curator, a larger budget,
or more rounds could use the feedback where these do not. That scoping sentence ships with
whatever the result is.

## 6. Provenance

* Implementation: `analysis/s2_shuffled_feedback.py`, written after this file.
* Output: `evidence/s2_shuffled_feedback.json`, `evidence/s2_shuffled_model_calls.jsonl`.
* Ledger rows are added by `analysis/make_ledger.py`; no existing row is modified.
* This file is hashed into the run output, so the reporting rule above is verifiable against
  the version that existed before the first call.

---

## 7. Post-run note — how the observed outcome maps onto §4 (added after the run)

Recorded because the outcome did not land cleanly inside one row of the table above,
and the honest thing is to say so rather than pick the nearest row.

**Observed.** All four paired real-minus-shuffled intervals include zero:
hold-out gain **+0.0067** [−0.0132, +0.0297]; full-test gain +0.0020 [−0.0098, +0.0144];
leakage-clean gain −0.0059 [−0.0268, +0.0230]; overstatement +0.0080 [−0.0123, +0.0243].
18 calls, **0 fallbacks**, so outcome D does not apply.

**The deviation.** Row **A** was written as a conjunction: the interval includes zero
*"and selection agreement is high"*. The performance half held; **the selection-agreement
half did not.** Mean Jaccard between the real-feedback and shuffled-feedback cluster
selections is **0.2494** [0.1593, 0.3460] over 18 matched calls — that is *low*
agreement. The curator plainly reads the benchmark-error column and acts on it: destroy
the signal and it picks substantially different clusters.

**How this is reported.** The performance criterion is the one §4 names — *"if shuffled feedback performs indistinguishably from real feedback,
that is a finding that belongs in the abstract"* — and it is met on every scored set, so
the finding goes in the abstract as outcome A prescribes. But it is stated in the sharper
form the data actually supports, and *not* as the weaker "the curator ignores the
feedback" that a high-agreement result would have licensed:

> Benchmark feedback changed what the curator selected — only about a quarter of clusters
> overlap with the real-feedback arm — but it did not change where the loop ended on any
> of the three scored sets.

Writing it the other way round would have been a claim the selection-agreement statistic
refutes. The conjunct that failed is published here and the statistic is reported in the
paper, exactly as §4 required of every outcome ("the selection-agreement statistic is
reported in all outcomes regardless of the hold-out result").

**Scope, as pre-committed in §5.** Six paired units, two small local curators, one
benchmark, one specialized-model family. This is a null about *this* loop at *this*
scale, not about curator LLMs in general.
