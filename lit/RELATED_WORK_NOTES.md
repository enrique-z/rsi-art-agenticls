# RELATED WORK NOTES — agenticls-rsi-leakage

One comparative sentence per entry: what they did / what we do differently.
Format: **bibkey** — comparison sentence.

---

## Core comparisons (~8 most important)

**huang2021tdc** (TDC, NeurIPS Datasets & Benchmarks 2021) — Huang et al. introduced the Therapeutics Data Commons benchmark suite with scaffold-based splits intended to mitigate train-test leakage; we demonstrate that even with scaffold splits, a non-trivial fraction of TDC ADMET held-out test molecules retain Tanimoto ≥ 0.9 near-duplicate neighbours in ChEMBL 36, enabling an LLM curator operating on ChEMBL to effectively observe test-analogous structure-activity data before evaluation.

**schuh2026auditing** (Schuh et al., Chemical Science 2026) — Schuh et al. audit widely used biomolecular benchmarks and find systematic data annotation inconsistencies (duplicate or conflicting labels across data sources); our contribution is orthogonal: we measure structural near-duplicate overlap between TDC test molecules and the ChEMBL corpus that an LLM curator can access, and quantify the resulting AUC gap in an RSI loop rather than studying annotation quality.

**shumailov2024collapse** (Shumailov et al., Nature 2024) — Shumailov et al. show that models trained on recursively self-generated synthetic data collapse in quality over generations; we study the complementary failure mode where a general-purpose LLM curates real (ChEMBL) data for bio-specialized ADMET models, inflating apparent benchmark gains without synthetic-data collapse.

**vantilborg2022activity** (van Tilborg et al., JCIM 2022) — van Tilborg et al. reveal that activity cliffs (pairs of structurally similar but activity-divergent molecules) expose molecular ML limitations even when benchmarks appear saturated; our work identifies a distinct leakage mechanism (Tanimoto ≥ 0.9 near-duplicates spanning the train/test split boundary via ChEMBL), which inflates AUC rather than exposes it as too variable.

**palavalli2024taxonomy** (Palavalli et al., arXiv 2024) — Palavalli et al. provide a taxonomy of data contamination modes for large language models evaluated on text benchmarks; we operationalize one contamination mode (structural analogue contamination via an external corpus) in the molecular property-prediction domain under an agentic RSI loop, making the taxonomy concrete for ADMET-style numeric benchmarks.

**boiko2023autonomous** (Boiko et al., Nature 2023) — Boiko et al. demonstrated that an LLM can autonomously perform chemical experiments and retrieve literature to guide synthesis decisions; our paper uses the same agentic-loop architecture but focuses on the data curation sub-task, showing that when the LLM selects training examples from ChEMBL for a downstream ADMET model, its access to the full ChEMBL corpus introduces benchmark leakage.

**ott2022mapping** (Ott et al., Nature Communications 2022) — Ott et al. map macro-level benchmark saturation dynamics across AI fields, finding that performance on many benchmarks converges irrespective of true capability gains; we provide a mechanistic account for one source of artificial saturation in ADMET benchmarks: near-duplicate leakage from public chemical databases amplified by an LLM curator.

**gao2022scaling** (Gao et al., arXiv 2022) — Gao et al. show that reward models over-optimized by RL agents degrade in true value even as proxy scores rise (Goodhart's Law in RLHF); we document the analogous effect in RSI loops for bio-ML: an LLM data curator that optimises for benchmark AUC (its implicit reward signal) achieves gains partly by surfacing test-analogous ChEMBL compounds, making reported progress an unreliable indicator of generalisation.

---

## Supporting entries

**wu2018moleculenet** — MoleculeNet established the original ADMET benchmarks that TDC extended; we show the leakage phenomenon is present in the TDC restatement with its updated ChEMBL-aligned chemical space.

**yang2019chemprop** — Chemprop (message-passing neural network) is the downstream model in our RSI loop; its state-of-the-art ADMET AUC makes leaked gains harder to attribute to architecture rather than data quality.

**ross2022molformer** — MolFormer is the alternative bio-specialized model we evaluate under RSI; AUC gaps between leaked and non-leaked subsets are consistent across MolFormer and Chemprop, ruling out model-specific artefacts.

**golbraikh2002beware** — Golbraikh and Tropsha's foundational warning about over-optimistic cross-validation metrics in QSAR is the intellectual ancestor of our analysis; we extend their concern to a new failure mode (LLM-mediated corpus leakage) not present in the classical QSAR setting.

**rogers2010ecfp** and **bajusz2015tanimoto** — ECFP4 Morgan fingerprints and Tanimoto similarity are the operational definitions of "near-duplicate" in our leakage census; we justify the Tanimoto ≥ 0.9 threshold with Bajusz et al.'s analysis of fingerprint-based similarity limits.

**bemis1996molecular** — Bemis-Murcko scaffold splits are the TDC default split strategy; we verify that scaffold splits do not fully resolve ChEMBL near-duplicate leakage when the leakage operates at the individual-molecule level rather than at the scaffold level.
