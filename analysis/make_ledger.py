#!/usr/bin/env python3
"""make_ledger.py — turn the S1 / S1b / S2 result JSONs into EGCL ledger rows
and the stamped-cell manifest.

Idempotent: the ledger file is rewritten from scratch on every run, so the same
real inputs always produce the same rows.  Nothing is computed here that is not
already in a result file produced by a real run — this script only reshapes.

Outputs
    evidence/agenticls-rsi-leakage_ledger.jsonl
    evidence/agenticls-rsi-leakage_stamp_manifest.json
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

PAPER_SLUG = "agenticls-rsi-leakage"
HERE = Path(__file__).resolve().parent
EV = HERE.parent / "evidence"
LEDGER = EV / f"{PAPER_SLUG}_ledger.jsonl"
MANIFEST = EV / f"{PAPER_SLUG}_stamp_manifest.json"

ROWS: list[dict] = []


def add(claim_id, value, source, verdict_reason, metric="value",
        verdict="SUPPORTED", ci_low=None, ci_high=None, must_print=True, **extra):
    ROWS.append({
        "claim_id": claim_id, "metric": metric, "value": value,
        "source": source, "verdict": verdict, "verdict_reason": verdict_reason,
        "ci_low": ci_low, "ci_high": ci_high, "must_print": must_print,
        "timestamp": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        **extra,
    })


def pct(x):
    return round(100.0 * x, 2)


def main() -> int:
    # ---------------- S1 ----------------------------------------------------
    s1_path = EV / "s1_census_aggregates.json"
    assert s1_path.exists(), f"missing {s1_path} — run s1_census.py first"
    s1 = json.loads(s1_path.read_text())
    src1 = str(s1_path)
    assert not s1["pilot"], "S1 aggregates are from a PILOT run — refusing to ledger"

    corpus_n = s1["chembl_corpus"]["n_fingerprints"]
    add("s1_corpus::n_chembl_fingerprints", corpus_n, src1,
        f"ChEMBL 36 compound_structures parsed by RDKit "
        f"({s1['chembl_corpus']['build']['n_rdkit_parse_failures']} parse failures "
        f"of {s1['chembl_corpus']['build']['n_rows_table']} rows)", metric="count")
    a = s1["exactness_audit"]
    add("s1_corpus::audit_max_abs_deviation", a["max_abs_deviation_vs_popcount"], src1,
        f"exactness audit on {a['audited_queries']} queries "
        f"({a['audit_fraction']*100:.1f}%) vs integer popcount and RDKit; "
        f"{a['argmax_mismatch_nontie']} non-tie argmax mismatches",
        metric="max_abs_deviation", must_print=False)

    n_q = 0
    for b, e in s1["benchmarks"].items():
        n_q += e["n_test"] + e["n_train_val"]
        k = f"s1_{b}"
        t = e["test_vs_chembl"]
        add(f"{k}::n_test", e["n_test"], src1,
            f"{e['display']} TDC ADMET Benchmark Group test split", metric="count")
        for th, lab in (("ge_0p9", "0p9"), ("ge_0p95", "0p95"), ("ge_1p0", "1p0")):
            r = t[th]
            add(f"{k}::leak_rate_{lab}", pct(r["rate"]), src1,
                f"{r['k']}/{e['n_test']} test molecules with max Tanimoto to ChEMBL 36 "
                f">= {lab.replace('p', '.')}; Wilson 95% CI",
                metric="percent", ci_low=pct(r["ci_low"]), ci_high=pct(r["ci_high"]),
                must_print=(lab in ("0p9", "0p95")))
        d = e["exact_smiles_duplicate"]
        add(f"{k}::exact_dup_rate", pct(d["rate"]), src1,
            f"{d['k']}/{e['n_test']} test molecules whose RDKit canonical SMILES is "
            f"present verbatim in ChEMBL 36", metric="percent",
            ci_low=pct(d["ci_low"]), ci_high=pct(d["ci_high"]))
        add(f"{k}::median_nn_sim", round(t["median"], 4), src1,
            "median nearest-neighbour Tanimoto of test molecules against ChEMBL 36",
            metric="median", must_print=False)
        own = e["test_vs_own_trainval"]["ge_0p9"]
        add(f"{k}::own_trainval_rate_0p9", pct(own["rate"]), src1,
            f"reference scale: {own['k']}/{e['n_test']} test molecules with a "
            f">=0.9 near-duplicate in the benchmark's OWN train_val split",
            metric="percent", ci_low=pct(own["ci_low"]), ci_high=pct(own["ci_high"]))
        sc = e["scaffold_clustering"]
        add(f"{k}::scaffold_leak_rate", pct(sc["scaffold_level_leak_rate"]), src1,
            f"{sc['n_distinct_scaffolds_with_leak_ge_0p9']}/"
            f"{sc['n_distinct_test_scaffolds']} distinct Bemis-Murcko test scaffolds "
            f"have >=1 member at NN>=0.9", metric="percent", must_print=False)
    # Internal duplication inside the test splits themselves — recomputed here
    # from the raw per-molecule rows, not read from an aggregate.
    from collections import Counter
    per_mol = EV / "s1_nn_per_molecule.jsonl"
    if per_mol.exists():
        by_bench: dict[str, list[str]] = {}
        for line in per_mol.read_text().splitlines():
            r = json.loads(line)
            if r["split"] == "test":
                by_bench.setdefault(r["benchmark"], []).append(r["cansmi"])
        for b, cans in by_bench.items():
            c = Counter(cans)
            n_dup_rows = sum(v - 1 for v in c.values() if v > 1)
            add(f"s1_{b}::internal_duplicate_rows", n_dup_rows, str(per_mol),
                f"{len(cans)} test rows collapse to {len(c)} distinct RDKit "
                f"canonical structures: {n_dup_rows} rows are internal duplicates "
                f"of another test row", metric="count", must_print=False)

    add("s1_corpus::n_pairwise_comparisons", n_q * corpus_n, src1,
        f"{n_q} benchmark molecules x {corpus_n} ChEMBL fingerprints, exact search",
        metric="count", must_print=False)

    # ---------------- S1b ---------------------------------------------------
    s1b_path = EV / "s1b_gap_results.json"
    if s1b_path.exists():
        s1b = json.loads(s1b_path.read_text())
        assert not s1b["pilot"], "S1b results are from a PILOT run — refusing to ledger"
        src2 = str(s1b_path)
        for b, e in s1b["benchmarks"].items():
            # BBB has no ChEMBL-derived arm (see DESIGN_NOTES A2), so its primary
            # arm is the control arm B.
            primary_arm = ("B_tdc_trainval" if b == "bbb_martins"
                           else "A2_chembl_dedup")
            for arm, av in e["arms"].items():
                for model, mv in av["models"].items():
                    key = f"s1b_{b}_{arm}_{model}"
                    under = bool(mv.get("underpowered"))
                    primary = (arm == primary_arm and model == "random_forest"
                               and not under)
                    lo, hi = mv["gap_ci_seed_bootstrap"]
                    add(f"{key}::gap", round(mv["mean_gap"], 4), src2,
                        f"{e['display']} {e['gap_definition']}, mean over "
                        f"{mv['n_seeds']} seeds; arm={av['description']}",
                        metric="gap", ci_low=round(lo, 4), ci_high=round(hi, 4),
                        must_print=primary, n_seeds=mv["n_seeds"])
                    mlo, mhi = mv["gap_ci_molecule_bootstrap"]
                    add(f"{key}::gap_mol_ci_low", round(mlo, 4), src2,
                        "molecule-level bootstrap CI lower bound on the same gap",
                        metric="ci_bound", must_print=False)
                    add(f"{key}::gap_mol_ci_high", round(mhi, 4), src2,
                        "molecule-level bootstrap CI upper bound on the same gap",
                        metric="ci_bound", must_print=False)
                    add(f"{key}::score_leaked", round(mv["mean_score_leaked"], 4), src2,
                        f"{e['metric']} on the leaked subset (n={e['n_leaked']})",
                        metric=e["metric"], must_print=primary)
                    add(f"{key}::score_nonleaked", round(mv["mean_score_nonleaked"], 4), src2,
                        f"{e['metric']} on the non-leaked subset (n={e['n_nonleaked']})",
                        metric=e["metric"], must_print=primary)
                    sr = mv["similarity_response"]
                    add(f"{key}::spearman_rho_sim_error",
                        round(sr["spearman_rho"], 4), src2,
                        f"Spearman rho between NN-to-ChEMBL Tanimoto and "
                        f"per-molecule absolute error over all {sr['n_molecules']} "
                        f"test molecules (p={sr['spearman_p']:.2e}); the binary "
                        f"leaked/non-leaked split is degenerate on two of the three "
                        f"benchmarks, this statistic is not",
                        metric="spearman_rho", must_print=False)
            add(f"s1b_{b}::n_leaked", e["n_leaked"], src2,
                f"{e['display']} test molecules at NN>=0.9", metric="count",
                must_print=False)
            add(f"s1b_{b}::n_nonleaked", e["n_nonleaked"], src2,
                f"{e['display']} test molecules at NN<0.9", metric="count",
                must_print=(e["n_nonleaked"] < 20))

    # ---------------- rival mechanisms -------------------------------------
    riv_path = EV / "s1b_rival_mechanisms.json"
    if riv_path.exists():
        riv = json.loads(riv_path.read_text())
        src3 = str(riv_path)
        for b, e in riv["benchmarks"].items():
            for pairing, v in e["R2_label_agreement"].items():
                if not v.get("n_near_pairs"):
                    continue
                key = f"r2_{b}_{pairing}"
                dv = v["paired_difference_near_minus_random"]
                add(f"{key}::near_pair_stat", round(v["near_duplicate_pairs"]["mean"], 4),
                    src3, f"{v['statistic']} for {v['n_near_pairs']} near-duplicate "
                          f"pairs (Tanimoto>=0.9)", metric=v["statistic"],
                    ci_low=round(v["near_duplicate_pairs"]["ci"][0], 4),
                    ci_high=round(v["near_duplicate_pairs"]["ci"][1], 4))
                add(f"{key}::random_pair_stat", round(v["random_pairs"]["mean"], 4),
                    src3, f"{v['statistic']} for the same number of random pairs",
                    metric=v["statistic"],
                    ci_low=round(v["random_pairs"]["ci"][0], 4),
                    ci_high=round(v["random_pairs"]["ci"][1], 4))
                add(f"{key}::pair_diff", round(dv["mean"], 4), src3,
                    "near-duplicate minus random, paired bootstrap 95% CI",
                    metric="difference", ci_low=round(dv["ci"][0], 4),
                    ci_high=round(dv["ci"][1], 4))
            for combo, v in e["R3_matched_gap"].items():
                m = v.get("matched_on_mw_and_label", {})
                if "matched_gap_mean" not in m:
                    continue
                add(f"r3_{b}_{combo.replace('|', '_')}::matched_gap",
                    round(m["matched_gap_mean"], 4), src3,
                    f"coarsened-exact-matched gap ({v['strata_definition']}), "
                    f"{m['n_matched_per_side']} molecules per side",
                    metric="gap", ci_low=round(m["matched_gap_ci"][0], 4),
                    ci_high=round(m["matched_gap_ci"][1], 4), must_print=False)
            for cov in ("mw", "chembl_scaffold_count"):
                c = e["R3_covariates"][cov]
                add(f"r3_{b}_{cov}::smd", round(c["standardised_mean_difference"], 4),
                    src3, f"standardised mean difference in {cov}, leaked minus "
                          f"non-leaked test subsets", metric="smd", must_print=False)

    # ---------------- S2 ----------------------------------------------------
    s2_path = EV / "s2_rsi_results.json"
    if s2_path.exists():
        s2 = json.loads(s2_path.read_text())
        src4 = str(s2_path)
        if s2.get("pilot"):
            print("NOTE: S2 results are a PILOT — not ledgered")
        else:
            for arm, v in s2["arms"].items():
                key = f"s2_{arm.replace(':', '_').replace('.', '_')}"
                add(f"{key}::gain_full", round(v["mean_gain_full_test"], 4), src4,
                    f"round-3 minus round-0 MAE improvement on the full TDC test set, "
                    f"curator={v['curator']}, mean over {len(v['per_seed'])} seeds",
                    metric="gain")
                add(f"{key}::gain_nonleaked", round(v["mean_gain_nonleaked_test"], 4),
                    src4, "same improvement measured on the non-leaked test subset",
                    metric="gain")
                add(f"{key}::gain_holdout", round(v["mean_gain_scaffold_holdout"], 4),
                    src4, "same improvement on the pre-committed ChEMBL scaffold "
                          "hold-out", metric="gain", must_print=False)
            p = s2.get("pooled_llm_overstatement")
            if p:
                add("s2_pooled::overstatement", round(p["mean"], 4), src4,
                    f"full-test gain minus non-leaked gain, pooled over "
                    f"{p['n_units']} {p['unit']} units, bootstrap 95% CI",
                    metric="difference", ci_low=round(p["ci"][0], 4),
                    ci_high=round(p["ci"][1], 4))

    # ---------------- 2026-09-07 review round -------------------------------
    # C1 endpoint-pool census, matching ladder, C2 nested bootstrap,
    # C5 cluster leak readout.  All produced by analysis/review_round.py and
    # analysis/matching_ladder.py from artefacts that already existed; no row
    # above is modified.
    rr_path = EV / "review_round_2026_09_07.json"
    ml_path = EV / "matching_ladder.json"
    if rr_path.exists():
        rr = json.loads(rr_path.read_text())
        src5 = str(rr_path)

        # --- C1: exposure against the endpoint-labelled pool only -----------
        for b, e in rr["A_endpoint_pool_census"].items():
            if not e.get("available"):
                continue
            add(f"s1_{b}::endpoint_pool_rate_0p9", e["rate_0p9_vs_endpoint_pool"], src5,
                f"{e['n_leaked_vs_endpoint_pool']} of {e['n_test']} test molecules have a "
                f"Tanimoto >= 0.9 neighbour in the {e['pool_size']}-compound ChEMBL "
                f"endpoint-labelled pool for this task (Wilson 95% CI); this is the "
                f"label-exploitable exposure, against the structural rate measured over "
                f"all of ChEMBL 36", metric="percent",
                ci_low=e["ci_low"], ci_high=e["ci_high"])
            add(f"s1_{b}::endpoint_pool_rate_ci_low", e["ci_low"], src5,
                "Wilson lower bound for the endpoint-pool exposure rate",
                metric="ci_bound", must_print=False)
            add(f"s1_{b}::endpoint_pool_rate_ci_high", e["ci_high"], src5,
                "Wilson upper bound for the endpoint-pool exposure rate",
                metric="ci_bound", must_print=False)
            add(f"s1_{b}::endpoint_pool_size", e["pool_size"], src5,
                "compounds in the ChEMBL endpoint-labelled pool for this task",
                metric="count", must_print=False)
            add(f"s1_{b}::n_leaked_endpoint_pool", e["n_leaked_vs_endpoint_pool"], src5,
                "test molecules with a >= 0.9 neighbour in that pool",
                metric="count", must_print=False)

        # --- C2: nested bootstrap and the paired A2-minus-B difference ------
        for key, v in rr["BC_matching"]["nested_bootstrap"].items():
            arm, model = key.split("|")
            cid = f"r3_solubility_aqsoldb_{arm}_{model}"
            add(f"{cid}::matched_gap_nested", v["nested_bootstrap_mean"], src5,
                f"matched gap under an OUTER bootstrap over test molecules "
                f"({v['n_outer_draws']} draws, re-stratified and re-matched inside each), "
                f"so the interval covers sampling uncertainty and not only matching "
                f"variability", metric="gap", must_print=False,
                ci_low=v["nested_bootstrap_ci"][0], ci_high=v["nested_bootstrap_ci"][1])
            add(f"{cid}::matched_gap_nested_ci_low", v["nested_bootstrap_ci"][0], src5,
                "nested-bootstrap lower bound", metric="ci_bound", must_print=False)
            add(f"{cid}::matched_gap_nested_ci_high", v["nested_bootstrap_ci"][1], src5,
                "nested-bootstrap upper bound", metric="ci_bound", must_print=False)
        for model, v in rr["BC_matching"]["paired_a2_minus_b"].items():
            cid = f"r3_solubility_aqsoldb_a2_minus_b_{model}"
            add(f"{cid}::paired_difference", v["paired_difference_mean"], src5,
                f"matched gap of arm A2 minus matched gap of control arm B, PAIRED inside "
                f"each of {v['n_outer_draws']} outer bootstrap draws over the same "
                f"resampled test molecules — the ChEMBL-training-specific residual, "
                f"computed directly rather than inferred from two marginal intervals",
                metric="difference",
                ci_low=v["paired_difference_ci"][0], ci_high=v["paired_difference_ci"][1])
            add(f"{cid}::paired_difference_ci_low", v["paired_difference_ci"][0], src5,
                "paired-difference lower bound", metric="ci_bound", must_print=False)
            add(f"{cid}::paired_difference_ci_high", v["paired_difference_ci"][1], src5,
                "paired-difference upper bound", metric="ci_bound", must_print=False)

        # --- C5: cluster-level leak readout (a NEGATIVE) --------------------
        D = rr["D_cluster_leak_readout"]
        repro = D["per_cluster_test_reproduction_rate"]
        size = D["per_cluster_size"]
        n_dup = int(round(sum(r * s for r, s in zip(repro, size))))
        add("s2_cluster_readout::n_universe_near_duplicates", n_dup, src5,
            f"compounds in the {D['n_universe_compounds']}-compound S2 curation universe "
            f"that are Tanimoto >= 0.9 near-duplicates of a TDC test molecule; the "
            f"universe already excludes exact TDC-test duplicates and the pre-committed "
            f"hold-out, so the loop's action space contains almost no reproducible test "
            f"chemistry at all", metric="count")
        add("s2_cluster_readout::n_universe_compounds", D["n_universe_compounds"], src5,
            "curation-universe size after removing the hold-out and exact test duplicates",
            metric="count", must_print=False)
        sw = D["selection_weighted_reproduction_rate"]
        add("s2_cluster_readout::llm_selection_weighted_repro", sw["llm_curators"], src5,
            "selection-weighted mean test-reproduction rate of the clusters the LLM "
            "curators chose across all 18 calls", metric="rate")
        add("s2_cluster_readout::universe_mean_repro", sw["universe_mean"], src5,
            "test-reproduction rate of the curation universe as a whole; the random "
            "control reproduces this value by construction", metric="rate")
        sp = D["spearman_repro_vs_cluster_benchmark_mae"]
        add("s2_cluster_readout::spearman_repro_vs_bench_mae", sp["mean"], src5,
            f"Spearman rank correlation between a cluster's test-reproduction rate and "
            f"the per-cluster benchmark error the curator was shown, mean over "
            f"{len(sp['per_seed'])} round-0 seeds; NEGATIVE, i.e. the feedback signal "
            f"points a curator AWAY from reproduced neighbourhoods, not toward them",
            metric="spearman_rho")

    if ml_path.exists():
        ml = json.loads(ml_path.read_text())
        src6 = str(ml_path)
        for key, v in ml["rungs"].items():
            arm, model = key.split("|")
            cid = f"r3_solubility_aqsoldb_{arm}_{model}"
            add(f"{cid}::matched_gap_mw_only", v["mw_only"]["matched_gap"], src6,
                "gap after coarsened-exact matching on the MOLECULAR-WEIGHT quintile "
                "only", metric="gap", must_print=False)
            add(f"{cid}::matched_gap_label_only", v["label_only"]["matched_gap"], src6,
                "gap after coarsened-exact matching on the LABEL quintile only; compare "
                "with the MW-only rung to see which covariate carries the composition "
                "difference", metric="gap", must_print=False)

    # ---------------- C6: second-corpus census (DrugBank 5.1) ---------------
    c6_path = EV / "c6_drugbank_census.json"
    if c6_path.exists():
        c6 = json.loads(c6_path.read_text())
        src7 = str(c6_path)
        add("c6_drugbank::corpus_size", c6["parse"]["n_distinct_canonical"], src7,
            f"distinct RDKit-canonical structures extracted from DrugBank 5.1 "
            f"({c6['parse']['n_drug_records']} drug records, "
            f"{c6['parse']['n_rdkit_parse_failures']} parse failures), fingerprinted "
            f"with the same Morgan r=2 / 2048-bit settings as the ChEMBL census",
            metric="count")
        for b, e in c6["benchmarks"].items():
            add(f"c6_{b}::drugbank_rate_0p9", e["drugbank_rate_0p9"], src7,
                f"{e['n_leaked_drugbank']} of {e['n_test']} test molecules have a "
                f"Tanimoto >= 0.9 neighbour in DrugBank 5.1 alone (Wilson 95% CI)",
                metric="percent", ci_low=e["drugbank_ci"][0], ci_high=e["drugbank_ci"][1],
                must_print=False)
            add(f"c6_{b}::union_rate_0p9", e["union_rate_0p9"], src7,
                f"exposure against ChEMBL 36 OR DrugBank 5.1, taken as the per-molecule "
                f"maximum over both corpora; this is the measured version of the "
                f"manuscript's 'our rates are a lower bound' claim",
                metric="percent", ci_low=e["union_ci"][0], ci_high=e["union_ci"][1])
            add(f"c6_{b}::clean_subset_union", e["clean_subset_union"], src7,
                f"leakage-clean subset size once a second public corpus is counted, "
                f"against {e['clean_subset_chembl_only']} for ChEMBL alone",
                metric="count", must_print=False)
            add(f"c6_{b}::clean_lost_to_second_corpus", e["clean_molecules_lost_to_second_corpus"],
                src7, "leakage-clean molecules that stop being clean when DrugBank is "
                      "added — leakage-clean is corpus-relative, and this measures by how "
                      "much", metric="count", must_print=False)

    # ---------------- C7: feedback-shuffled control arm ---------------------
    c7_path = EV / "s2_shuffled_feedback.json"
    if c7_path.exists():
        c7 = json.loads(c7_path.read_text())
        src8 = str(c7_path)
        for arm, v in c7["arms"].items():
            key = f"c7_{arm.replace(':', '_').replace('.', '_')}"
            add(f"{key}::gain_full", round(v["mean_gain_full_test"], 4), src8,
                f"shuffled-feedback arm, round-3 minus round-0 on the full TDC test set, "
                f"curator={arm}", metric="gain", must_print=False)
            add(f"{key}::gain_nonleaked", round(v["mean_gain_nonleaked_test"], 4), src8,
                "shuffled-feedback arm, same on the leakage-clean subset",
                metric="gain", must_print=False)
            add(f"{key}::gain_holdout", round(v["mean_gain_scaffold_holdout"], 4), src8,
                "shuffled-feedback arm, same on the pre-committed scaffold hold-out",
                metric="gain", must_print=False)
        for k, v in c7["paired_real_minus_shuffled"].items():
            add(f"c7_paired::{k}", round(v["mean"], 4), src8,
                f"real-feedback minus shuffled-feedback, PAIRED on curator x seed "
                f"(both arms share round 0 exactly), {v['n_units']} units, percentile "
                f"bootstrap; interval "
                f"{'EXCLUDES' if v['ci_excludes_zero'] else 'INCLUDES'} zero",
                metric="difference",
                ci_low=round(v["ci"][0], 4), ci_high=round(v["ci"][1], 4),
                must_print=(k == "gain_scaffold_holdout"))
            add(f"c7_paired::{k}_ci_low", round(v["ci"][0], 4), src8,
                "paired bootstrap lower bound", metric="ci_bound", must_print=False)
            add(f"c7_paired::{k}_ci_high", round(v["ci"][1], 4), src8,
                "paired bootstrap upper bound", metric="ci_bound", must_print=False)
        sa = c7["selection_agreement"]
        add("c7_selection::jaccard", round(sa["mean"], 4), src8,
            f"mean Jaccard between the clusters chosen under real feedback and under "
            f"shuffled feedback, over {sa['n_calls']} calls matched on curator, seed and "
            f"round; LOW agreement means the curator does read the feedback column and "
            f"acts on it", metric="jaccard",
            ci_low=round(sa["ci"][0], 4), ci_high=round(sa["ci"][1], 4))
        add("c7_selection::jaccard_ci_low", round(sa["ci"][0], 4), src8,
            "bootstrap lower bound", metric="ci_bound", must_print=False)
        add("c7_selection::jaccard_ci_high", round(sa["ci"][1], 4), src8,
            "bootstrap upper bound", metric="ci_bound", must_print=False)
        add("c7::n_fallbacks", c7["n_fallbacks"], src8,
            f"parse fallbacks in the shuffled-feedback arm across {c7['n_calls']} calls",
            metric="count", must_print=False)

    # ---------------- small-n interval robustness (2026-09-07) --------------
    # A percentile bootstrap at n=3 is exactly [min, max] of the three values, so
    # "excludes zero" would be forced by the sign alone.  None of our intervals is
    # over n=3 -- the pooled ones are n=6 and n=18 -- but we report a paired-t
    # interval beside the pooled bootstrap because the t interval has no such
    # degeneracy at any n.  Computed here from the per-seed values already in
    # s2_rsi_results.json; no new measurement.
    s2p = EV / "s2_rsi_results.json"
    if s2p.exists():
        import math
        s2b = json.loads(s2p.read_text())
        per = []
        for arm, v in s2b["arms"].items():
            if v["curator"] == "random":
                continue
            per += [s["gain_full_test"] - s["gain_nonleaked_test"] for s in v["per_seed"]]
        n = len(per)
        mean = sum(per) / n
        sd = math.sqrt(sum((x - mean) ** 2 for x in per) / (n - 1))
        # two-sided 97.5th percentile of Student t at n-1 = 5 df
        tcrit = 2.570581835636197
        half = tcrit * sd / math.sqrt(n)
        add("s2_pooled::overstatement_t_ci_low", round(mean - half, 4), str(s2p),
            f"paired-t 95% lower bound for the pooled overstatement over {n} "
            f"curator x seed units; reported beside the percentile bootstrap because "
            f"a t interval is not degenerate at small n", metric="ci_bound")
        add("s2_pooled::overstatement_t_ci_high", round(mean + half, 4), str(s2p),
            "paired-t 95% upper bound for the pooled overstatement",
            metric="ci_bound", must_print=False)

    # The same paired-t treatment for the shuffled-feedback arm's four paired
    # quantities, which sit on the same 6 curator x seed units.  Reported beside
    # every c7 bootstrap so that a null does not rest on a small-n percentile.
    c7p = EV / "s2_shuffled_feedback.json"
    if c7p.exists():
        import math
        c7b = json.loads(c7p.read_text())
        for k, v in c7b["paired_real_minus_shuffled"].items():
            per = list(v["per_unit"])
            n = len(per)
            mean = sum(per) / n
            sd = math.sqrt(sum((x - mean) ** 2 for x in per) / (n - 1))
            tcrit = 2.570581835636197  # 97.5th pct of Student t, 5 df
            half = tcrit * sd / math.sqrt(n)
            add(f"c7_paired::{k}_t_ci_low", round(mean - half, 4), str(c7p),
                f"paired-t 95% lower bound over the same {n} curator x seed units "
                f"as the percentile bootstrap for {k}", metric="ci_bound",
                must_print=(k == "gain_scaffold_holdout"))
            add(f"c7_paired::{k}_t_ci_high", round(mean + half, 4), str(c7p),
                f"paired-t 95% upper bound for {k}", metric="ci_bound",
                must_print=False)

    # ---------------- interval degeneracy audit (2026-09-07) ----------------
    # analysis/ci_degeneracy_audit.py walks every interval that reaches the PDF
    # or the ledger, resolves the resampling unit behind it, and flags any whose
    # percentile bounds coincide with the [min, max] of its units or have zero
    # width.  Its counts are stamped here so the paper can quote them.
    audit_p = EV / "ci_degeneracy_audit.json"
    if audit_p.exists():
        aud = json.loads(audit_p.read_text())
        add("ci_audit::n_intervals", aud["n_intervals_audited"], str(audit_p),
            "intervals reaching the PDF or the ledger, each resolved to its "
            "resampling unit and unit count by analysis/ci_degeneracy_audit.py",
            metric="count")
        add("ci_audit::n_degenerate", aud["n_degenerate"], str(audit_p),
            "intervals whose percentile bounds coincide with the [min, max] of "
            "their units or have zero width: " +
            (", ".join(aud["degenerate_intervals"]) or "none"),
            metric="count")

    # ---------------- C8, and why nothing from it is published --------------
    c8p = EV / "c8_remediation.json"
    if c8p.exists():
        c8 = json.loads(c8p.read_text())
        rep = c8["reproduction_check"]
        add("c8::holdout_repro_deviation", round(rep["max_holdout_gain_deviation"], 4),
            str(c8p),
            "maximum absolute deviation between the C8 re-fit hold-out gains and "
            "the published S2 hold-out gains, against a pre-committed tolerance of "
            f"{rep['tolerance']}; the hold-out does not depend on the test split, so "
            "a mismatch voids the re-scoring under outcome D of the C8 "
            "pre-commitment and no re-scored number is published",
            metric="max_abs_deviation",
            must_print=True)

    # ---------------- write -------------------------------------------------
    with open(LEDGER, "w") as fh:
        for r in ROWS:
            fh.write(json.dumps(r) + "\n")
    manifest = {"ledger_rows": {
        r["claim_id"]: {"value": r["value"],
                        "source": f"{LEDGER.name}:{i+1}",
                        "must_print": r["must_print"]}
        for i, r in enumerate(ROWS)}}
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")

    from collections import Counter
    c = Counter(r["verdict"] for r in ROWS)
    print(f"wrote {LEDGER} ({len(ROWS)} rows)")
    print(f"wrote {MANIFEST}")
    print("verdicts:", dict(c))
    print("must_print rows:", sum(1 for r in ROWS if r["must_print"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
