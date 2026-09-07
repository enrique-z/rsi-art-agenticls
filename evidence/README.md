# Evidence ledger — conventions

Every numeric literal reported in the manuscript (except layout/architecture
constants) resolves to a row in one of this folder's ledgers. The
anti-fabrication guard (`../analysis/guard.py`) fails if it cannot resolve one.

## Ledger row schema

One JSONL file per ledger; one row per claim:

```json
{
  "claim_id": "cohort_size_n",
  "metric": "count",
  "value": 500,
  "source": "evidence/p1_cohort/D1.jsonl",
  "verdict": "SUPPORTED",
  "verdict_reason": "D1.jsonl line count = 500 (target 500)",
  "ci_low": null,
  "ci_high": null,
  "must_print": true,
  "timestamp": "2026-08-24T17:34:16Z"
}
```

Verdicts:

- `SUPPORTED` — ships. The value traces to a real measurement, a real
  database, or a cited source.
- `UNSUPPORTED` — blocks. The corresponding literal must be withdrawn.
- `NOT_APPLICABLE` — informational only; not a release path.
- `CONTESTED` — needs human review; blocks release until resolved.

## The `must_print` field

`must_print` controls which ledger rows the guard requires to appear in the
manuscript.

```
must_print: true   -- the check FAILS if this row is omitted from the manuscript.
                      Used for every pre-registered outcome and headline claim.
must_print: false  -- the row is recorded for audit but need not be typeset.
                      Used for decomposition/robustness rows that back a printed
                      RANGE: the range they support is what gets typeset, not
                      the decomposition row itself.
```

A manifest with **no** `must_print` flags at all falls back to the older
behaviour in which every ledger row must be typeset. Adding even one
`must_print: true` row activates the selective mode for the whole manifest, so
the field is set explicitly on every row from the first run.

Typical `must_print=False` rows:

- per-cell breakdown rows behind a summary range claim
- robustness checks (the same model on a held-out slice)
- convergence diagnostics (`n_fits`, `n_converged`, `final_grad_norm`)
- any row whose value is not directly typeset in prose or tables

## The stamped-cell manifest

For per-cell value checking, the paper also emits a
`<paper>_stamp_manifest.json` mapping `\StampedCell{k1}{k2}{v}` keys to the
ledger row they must equal:

```json
{
  "ledger_rows": {
    "main_result::accuracy": {"value": 0.811, "source": "ledger.jsonl:42"},
    "main_result::ci_low":   {"value": 0.718, "source": "ledger.jsonl:43"}
  }
}
```

The guard walks the manuscript for every `\StampedCell{k1}{k2}{v}` and asserts
`manifest.ledger_rows["k1::k2"].value == v` within numeric tolerance.

## The `quantity_descriptor` field (optional, `QD_ENFORCE=1`)

An additive, optional sub-object on any ledger row. It records the extraction
provenance behind a stamped number so the guard can check that the surrounding
prose describes the *same* slice of data:

```json
{
  "metric_class": "RATIO",
  "base_artifact": "cohort_v3.jsonl",
  "numerator_predicates": ["status == active"],
  "denominator_predicates": ["all_rows"],
  "silent_drops": ["left_censored"],
  "population_n_before": 500,
  "population_n_after": 483,
  "extracted_by": "aggregates.py::compute_active_rate",
  "extractor_version": "0.1"
}
```

`metric_class` ∈ {COUNT, RATIO, MEAN, VARIANCE}. The predicate lists give the
filters actually applied in code, in a shape a linter can grep (`==`, `>=`,
`!=`, …). `silent_drops` names every row-eliminating filter the code applied —
an empty list is a positive claim that no drops happened.

When `QD_ENFORCE=1` is set, the guard checks, for every stamped cell whose
manifest row carries a descriptor:

- **(a) predicate presence** — every `numerator_predicate` appears in the
  sentence containing the macro, either literally or via the synonym table at
  the top of `guard.py` (`==` ↔ "called"/"labeled"/"is", `>=` ↔ "at least"/"≥").
  String-set logic only; no model call.
- **(b) undeclared drop** — `population_n_before - population_n_after > 0`
  with `silent_drops == []` is a hard failure: rows were dropped and the ledger
  denies it.
- **(c) declared-but-hidden drop** — every `silent_drops` tag appears in the
  sentence or in a same-paragraph footnote, so a declared exclusion is not
  hidden from the prose.

## GRIM lint (warning level, always on)

Independently of the descriptor checks, the guard runs a GRIM
(Granularity-Related Inconsistency of Means) check on every stamped cell whose
manifest row records an integer sample size — keys `n`, `n_units`, `sample_n`
or `population_n_after`. It verifies the reported mean is reachable as `k/n` at
the printed decimal precision. Rows without an integer `n`, and rows whose
`metric_class` is not `MEAN`, are skipped. The warning quotes the mean, `n`,
and the nearest reachable values above and below.

## Reading the ledgers safely

- **Claim-row status is effective, not ADD-time.** A claim's status is the
  `new_state` of its last transaction; reading the first row fabricates metrics.
- **Refutation transactions resolve on `claim_id`**, not on a template id.
- **Ledger keys can drift semantically.** Recompute subset counts from the raw
  cells rather than trusting a key name, and check whether the same quantity is
  stated twice.
- **A `SUPPORTED` verdict is a label.** Read the function that produced it
  before relying on the row.
- **Search the whole evidence tree** before calling a number unsupported; a
  partial search produces false fabrication alarms.
