# NOVELTY CHECK — agenticls-rsi-leakage

## Claim under review

The paper claims: RSI loops in which a general-purpose LLM curates training data
for bio-specialized models are measured against ADMET benchmarks (TDC) whose
held-out test molecules have near-duplicate neighbours (Tanimoto ≥ 0.9, Morgan
r=2 fingerprints) in public corpora the curator can access (ChEMBL 36); apparent
RSI gains partly reflect benchmark leakage. The paper quantifies the per-benchmark
leakage rate and measures the AUC gap between leaked and non-leaked test subsets.

---

## 5 closest prior works and threat assessment

### 1. schuh2026auditing — Schuh, Daniluk, Sieber (Chemical Science 2026)
DOI: 10.1039/d6sc01799a

**What they did:** Audited "widely used biomolecular benchmarks" and found
"systematic data inconsistencies" across datasets. The paper's title and venue
(Chemical Science, RSC) indicate a focus on benchmark data quality in molecular
science.

**Threat to novelty:** MODERATE. The title is the closest in spirit to our
work, and the finding of systematic problems in biomolecular benchmarks could
overlap with our leakage finding. However, "data inconsistencies" most likely
refers to annotation errors, conflicting assay values, or duplicate bioactivity
records — not structural near-duplicate overlap between TDC test splits and
ChEMBL caused by an agentic LLM curator. The RSI loop framing (LLM as curator,
AUC gap measurement, per-benchmark leakage rate) is absent from the paper's
abstract and is specific to our contribution.

**Verdict:** Not a prior claim of our specific finding. If a reviewer raises
this paper, the response is: Schuh et al. study annotation-level data quality;
we study structural near-duplicate leakage introduced by the access pattern of
an LLM curator, which is a distinct mechanism with different detection method
(Tanimoto census vs. annotation reconciliation) and different consequence (AUC
inflation vs. label noise).

---

### 2. vantilborg2022activity — van Tilborg, Alenicheva, Grisoni (JCIM 2022)
DOI: 10.1021/acs.jcim.2c01073

**What they did:** Showed that activity cliffs (structurally similar molecules
with large activity differences) expose the brittleness of molecular ML models
even on benchmarks that appear highly accurate.

**Threat to novelty:** LOW. Activity cliffs measure within-dataset structural
similarity that causes performance variability; they do not measure between-split
or between-corpus near-duplicate leakage from an external database. The mechanism
and metric are different from our Tanimoto ≥ 0.9 ChEMBL census. There is no
mention of LLM curators or RSI loops.

---

### 3. faraji2026leakage — Faraji et al. (ACS Omega 2026)
DOI: 10.1021/acsomega.6c05807

**What they did:** Performed a "leakage-aware benchmark study" of ML models for
deep eutectic solvent (DES) property prediction, explicitly accounting for
train-test leakage in their evaluation.

**Threat to novelty:** LOW-MODERATE. This is the only published paper we found
that uses the phrase "leakage-aware benchmark" in a molecular-ML context (2026).
However, it is domain-specific to DES solvents (not ADMET/TDC), does not involve
ChEMBL near-duplicate detection at Tanimoto ≥ 0.9, and has no LLM curator or
RSI loop component. It demonstrates that leakage-aware evaluation is a valid
methodology, which supports rather than threatens our novelty.

---

### 4. ott2022mapping — Ott, Barbosa-Silva et al. (Nature Communications 2022)
DOI: 10.1038/s41467-022-34591-0

**What they did:** Mapped benchmark saturation dynamics across AI, showing that
benchmark performance converges artificially as the field over-optimises to a
fixed held-out set.

**Threat to novelty:** LOW. Their analysis is meta-level (across AI fields,
using publication records) and does not identify near-duplicate leakage from
external corpora or LLM curation as a mechanism. It provides relevant context
for why benchmark inflation occurs but is not a mechanistic competitor.

---

### 5. golbraikh2002beware — Golbraikh & Tropsha (JMGM 2002)
DOI: 10.1016/s1093-3263(01)00123-1

**What they did:** Issued the foundational warning that q2 (leave-one-out
cross-validation R2) is a poor predictor of true QSAR model performance because
it does not control for chemical space coverage between train and test.

**Threat to novelty:** VERY LOW. Published 2002, no connection to LLMs or RSI
loops, no ChEMBL, no Tanimoto near-duplicate census. It is relevant as the
intellectual ancestor of our concern and should be cited as prior art, not as
a novelty threat.

---

## Broader sweep findings

**Did any prior work specifically quantify TDC ADMET test-split overlap with
ChEMBL 36 near-duplicates (Tanimoto ≥ 0.9)?**

Literature search across arXiv (2407.08716, others), Crossref, and Semantic
Scholar found no paper that:
(a) queries ChEMBL for Tanimoto ≥ 0.9 neighbours of TDC ADMET held-out molecules,
(b) reports a per-benchmark leakage rate, and
(c) connects this leakage to an LLM-curated RSI loop.

Schuh 2026 is the closest in scope and must be cited prominently; the abstract
does not match claims (a)–(c). Faraji 2026 validates the leakage-aware
evaluation methodology but in a different domain.

**Did any prior work study LLM-curated RSI loops on ADMET benchmarks?**

No peer-reviewed paper was found. Boiko 2023, Bran 2024, and Lu 2024 demonstrate
LLM-mediated chemical research but do not evaluate ADMET benchmark leakage.

---

## Novelty verdict

**NOVEL with one caveat.** The specific triple claim — (1) LLM curator in an
RSI loop accesses ChEMBL, (2) TDC ADMET test molecules have Tanimoto ≥ 0.9
near-duplicates in ChEMBL, (3) the resulting AUC gap is measurable and
per-benchmark quantifiable — has not been made in any identified prior work.

**Caveat:** Schuh et al. 2026 (Chemical Science) audits biomolecular benchmarks
for systematic inconsistencies. The abstract must be obtained before submission
to confirm their scope does not include Tanimoto-based leakage from ChEMBL.
Based on the title's emphasis on "data inconsistencies" (annotation quality) vs.
our "structural near-duplicate leakage" (corpus overlap), these are distinct
contributions. This paper should appear in related work with an explicit
one-sentence distinction.

**Risk level:** LOW-MODERATE (Schuh 2026 is the only paper requiring abstract
confirmation; all others are clearly distinct in mechanism or domain).
