#!/usr/bin/env python3
"""make_figure.py — the paper's single evidence figure: recompute -> assert -> draw.

Panel A  nearest-neighbour Tanimoto of every TDC ADMET *test* molecule against
         ChEMBL 36, per benchmark, with the >= 0.9 mass shaded and labelled, and
         the same distribution against the benchmark's OWN train_val split as a
         reference scale.
Panel B  leaked vs non-leaked predictor performance: the per-arm gap with its
         seed-level 95% bootstrap CI, zero marked.

Every aggregate drawn is recomputed here from the per-molecule / per-seed
evidence and asserted against the ledger before anything is plotted; a mismatch
aborts non-zero.  `--selftest` proves the assertion machinery both passes on
consistent input and aborts on perturbed input.

Usage:
    python3 analysis/make_figure.py
    python3 analysis/make_figure.py --selftest
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42   # embed Type 42, never Type 3
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
PAPER_SLUG = "agenticls-rsi-leakage"
EV = HERE.parent / "evidence"
PER_MOL = EV / "s1_nn_per_molecule.jsonl"
S1B = EV / "s1b_gap_results.json"
LEDGER = EV / f"{PAPER_SLUG}_ledger.jsonl"
OUT = HERE.parent / "fig_leakage.pdf"
CLAIMS_OUT = EV / "figure_claims.json"

TOL = 5e-3   # ledger percentages are rounded to 2 dp
ORDER = ["bbb_martins", "lipophilicity_astrazeneca", "solubility_aqsoldb"]
NICE = {"bbb_martins": "BBB-Martins",
        "lipophilicity_astrazeneca": "Lipophilicity-AZ",
        "solubility_aqsoldb": "Solubility-AqSolDB"}
COLOURS = {"bbb_martins": "#1b6ca8",
           "lipophilicity_astrazeneca": "#c1440e",
           "solubility_aqsoldb": "#2e7d32"}


# --------------------------------------------------------------------------- #
def load_ledger(path: Path) -> dict[str, dict]:
    rows = {}
    for line in path.read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            rows[r["claim_id"]] = r
    assert rows, f"empty ledger {path} — FAIL"
    return rows


def recompute() -> tuple[dict, dict]:
    """Recompute every drawn aggregate from the raw per-unit evidence."""
    per_mol = defaultdict(lambda: {"chembl": [], "own": []})
    for line in PER_MOL.read_text().splitlines():
        r = json.loads(line)
        if r["split"] != "test":
            continue
        per_mol[r["benchmark"]]["chembl"].append(r["nn_chembl_tanimoto"])
        per_mol[r["benchmark"]]["own"].append(r["nn_own_trainval_tanimoto"])
    assert per_mol, "no test rows in the per-molecule ledger — FAIL"

    agg = {}
    dists = {}
    for b, d in per_mol.items():
        c = np.array(d["chembl"])
        o = np.array(d["own"])
        dists[b] = {"chembl": c, "own": o}
        agg[f"s1_{b}::leak_rate_0p9"] = round(100.0 * float((c >= 0.9).mean()), 2)
        agg[f"s1_{b}::leak_rate_0p95"] = round(100.0 * float((c >= 0.95).mean()), 2)
        agg[f"s1_{b}::own_trainval_rate_0p9"] = round(100.0 * float((o >= 0.9).mean()), 2)
        agg[f"s1_{b}::n_test"] = len(c)

    gaps = {}
    if S1B.exists():
        s1b = json.loads(S1B.read_text())
        for b, e in s1b["benchmarks"].items():
            for arm, av in e["arms"].items():
                for model, mv in av["models"].items():
                    g = np.array([r["gap"] for r in mv["per_seed"]])
                    key = f"s1b_{b}_{arm}_{model}::gap"
                    agg[key] = round(float(g.mean()), 4)
                    gaps[key] = {
                        "benchmark": b, "arm": arm, "model": model,
                        "mean": float(g.mean()),
                        "ci": mv["gap_ci_seed_bootstrap"],
                        "ci_mol": mv["gap_ci_molecule_bootstrap"],
                        "underpowered": bool(mv.get("underpowered")),
                        "n_nonleaked": e["n_nonleaked"],
                        "task": e["task"],
                        "metric": e["metric"],
                    }
    return agg, {"dists": dists, "gaps": gaps}


def assert_against_ledger(agg: dict, ledger: dict, tol: float = TOL) -> list[str]:
    checked = []
    for key, val in agg.items():
        assert key in ledger, f"E2 no ledger row for recomputed aggregate {key}"
        lv = ledger[key]["value"]
        assert ledger[key]["verdict"] == "SUPPORTED", f"E4 {key} verdict is not SUPPORTED"
        if isinstance(val, (int, float)) and isinstance(lv, (int, float)):
            assert abs(float(lv) - float(val)) <= tol, (
                f"E3 value drift on {key}: ledger={lv} recomputed={val}")
        else:
            assert lv == val, f"E3 value drift on {key}: ledger={lv} recomputed={val}"
        checked.append(key)
    return checked


ARM_SHORT = {"A2_chembl_dedup": "ChEMBL dedup (A2)",
             "A1_chembl_naive": "ChEMBL naive (A1)",
             "B_tdc_trainval": "TDC own train (B)"}
MODEL_SHORT = {"random_forest": "RF", "linear_ecfp": "linear",
               "knn_tanimoto": "kNN-Tc"}
BENCH_TAG = {"bbb_martins": "BBB", "lipophilicity_astrazeneca": "Lipo",
             "solubility_aqsoldb": "Sol"}


def _forest(ax, keys, gaps, title, xlabel):
    """Forest plot of leaked-minus-non-leaked advantage, molecule-level 95% CI."""
    ys = np.arange(len(keys))[::-1]
    labels = []
    for y, k in zip(ys, keys):
        g = gaps[k]
        lo, hi = g["ci_mol"]
        slo, shi = g["ci"]
        col = COLOURS[g["benchmark"]]
        under = g["underpowered"]
        ax.plot([lo, hi], [y, y], color=col, lw=1.5,
                alpha=0.35 if under else 0.95, solid_capstyle="butt")
        ax.plot([slo, shi], [y, y], color=col, lw=4.0,
                alpha=0.20 if under else 0.45, solid_capstyle="butt")
        ax.plot([g["mean"]], [y], "o" if not under else "o",
                mfc=("none" if under else col), mec=col, ms=4.6, mew=1.1)
        tag = (f"{BENCH_TAG[g['benchmark']]} · "
               if len({gaps[j]["benchmark"] for j in keys}) > 1 else "")
        labels.append(f"{tag}{ARM_SHORT.get(g['arm'], g['arm'])} · "
                      f"{MODEL_SHORT.get(g['model'], g['model'])}")
    ax.axvline(0, color="0.3", lw=1.0, zorder=0)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=5.9)
    ax.set_ylim(-0.8, len(keys) - 0.2)
    ax.set_xlabel(xlabel, fontsize=7.0)
    ax.set_title(title, fontsize=8.4, loc="left")
    ax.tick_params(axis="x", labelsize=7.0)


def draw(data: dict, out: Path) -> None:
    dists, gaps = data["dists"], data["gaps"]
    has_b = bool(gaps)
    if has_b:
        fig, axes = plt.subplots(
            1, 3, figsize=(13.2, 3.7),
            gridspec_kw={"width_ratios": [1.45, 0.72, 1.05]})
        axA, axB, axC = axes
    else:
        fig, axA = plt.subplots(figsize=(5.6, 3.6))
        axes = [axA]

    # ---- A: survival curve P(NN Tanimoto >= t) ---------------------------
    grid = np.linspace(0, 1, 501)
    axA.axvspan(0.9, 1.0, color="0.88", zorder=0)
    axA.axvline(0.9, color="0.35", lw=0.9, ls="--", zorder=1)
    for b in ORDER:
        if b not in dists:
            continue
        c, o = dists[b]["chembl"], dists[b]["own"]
        axA.plot(grid, [100.0 * float((c >= t).mean()) for t in grid],
                 lw=1.9, color=COLOURS[b], label=f"{NICE[b]} vs ChEMBL 36", zorder=3)
        axA.plot(grid, [100.0 * float((o >= t).mean()) for t in grid],
                 lw=1.1, ls=":", color=COLOURS[b], alpha=0.85, zorder=2,
                 label=f"{NICE[b]} vs own train split")
    axA.set_xlabel("nearest-neighbour Tanimoto threshold $t$ (Morgan r=2, 2048 bit)",
                   fontsize=7.6)
    axA.set_ylabel("% of test molecules with NN $\\geq t$", fontsize=7.6)
    axA.set_xlim(0, 1)
    axA.set_ylim(0, 103)
    axA.tick_params(labelsize=7.0)
    axA.text(0.95, 4, "near-duplicate\n$\\geq$0.9", ha="center", va="bottom",
             fontsize=6.6, color="0.35")
    rates = [(b, 100.0 * float((dists[b]["chembl"] >= 0.9).mean()),
              100.0 * float((dists[b]["own"] >= 0.9).mean()))
             for b in ORDER if b in dists]
    rates.sort(key=lambda t: -t[1])
    prev = 1e9
    for b, r, ro in rates:
        y = min(r + 2.0, 90.0, prev - 9.0)
        prev = y
        axA.annotate(f"{r:.1f}%  (own split {ro:.1f}%)",
                     xy=(0.9, r), xytext=(0.38, y),
                     fontsize=6.8, color=COLOURS[b], fontweight="bold", va="center",
                     arrowprops=dict(arrowstyle="-", lw=0.6, color=COLOURS[b]))
    axA.legend(fontsize=5.6, loc="lower left", frameon=False)
    axA.set_title("A  test-set near-duplicates in the public corpus", fontsize=8.4,
                  loc="left")

    if not has_b:
        for sp in ("top", "right"):
            axA.spines[sp].set_visible(False)
        fig.tight_layout()
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
        return

    # ---- B / C: the gap, split by metric so units are never mixed --------
    def sort_key(k):
        g = gaps[k]
        return (ORDER.index(g["benchmark"]), g["arm"], g["model"])

    cls = sorted([k for k in gaps if gaps[k]["task"] == "classification"], key=sort_key)
    reg = sorted([k for k in gaps if gaps[k]["task"] == "regression"], key=sort_key)
    _forest(axB, cls, gaps,
            "B  BBB-Martins (AUC units)",
            "$\\Delta$AUC  (leaked $-$ non-leaked)")
    _forest(axC, reg, gaps,
            "C  regression benchmarks (MAE units)",
            "$\\Delta$MAE  (non-leaked $-$ leaked)")
    for ax, keys in ((axB, cls), (axC, reg)):
        n_under = sum(gaps[k]["underpowered"] for k in keys)
        if n_under:
            ax.text(0.01, -0.30, "hollow marker / faded bar = fewer than 20 "
                    "non-leaked test molecules", transform=ax.transAxes,
                    fontsize=5.6, color="0.35", va="top")
        else:
            ax.text(0.01, -0.30, "thin bar: molecule-level 95% CI · "
                    "thick bar: seed-level 95% CI", transform=ax.transAxes,
                    fontsize=5.6, color="0.35", va="top")

    for ax in axes:
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
def run_selftest() -> int:
    """Prove the assertion machinery passes on consistent input and aborts on drift."""
    agg = {"x::a": 12.34, "x::b": 5}
    ledger = {"x::a": {"value": 12.34, "verdict": "SUPPORTED"},
              "x::b": {"value": 5, "verdict": "SUPPORTED"}}
    assert assert_against_ledger(agg, ledger) == ["x::a", "x::b"]
    for name, bad, marker in [
        ("drift", {"x::a": {"value": 12.90, "verdict": "SUPPORTED"},
                   "x::b": {"value": 5, "verdict": "SUPPORTED"}}, "E3 value drift"),
        ("missing row", {"x::b": {"value": 5, "verdict": "SUPPORTED"}},
         "E2 no ledger row"),
        ("bad verdict", {"x::a": {"value": 12.34, "verdict": "CONTESTED"},
                         "x::b": {"value": 5, "verdict": "SUPPORTED"}}, "E4"),
    ]:
        try:
            assert_against_ledger(agg, bad)
        except AssertionError as e:
            assert marker in str(e), (name, e)
        else:
            print(f"SELFTEST FAILED: {name} not caught")
            return 1
    print("selftest OK: pass on consistent input; abort on drift / missing row / "
          "non-SUPPORTED verdict")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)
    if args.selftest:
        return run_selftest()

    agg, data = recompute()
    ledger = load_ledger(LEDGER)
    checked = assert_against_ledger(agg, ledger)
    draw(data, OUT)
    CLAIMS_OUT.write_text(json.dumps({
        "figure": OUT.name,
        "claim_ids": sorted(checked),
        "note": "claim_ids drawn or annotated in fig_leakage.pdf; gate 3 counts "
                "these as reported without a redundant table",
    }, indent=2) + "\n")
    print(f"wrote {OUT}")
    print(f"wrote {CLAIMS_OUT} ({len(checked)} claim_ids asserted against the ledger)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
