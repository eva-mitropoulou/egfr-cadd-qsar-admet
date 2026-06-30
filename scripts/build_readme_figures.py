#!/usr/bin/env python3
"""Build compact README figures from committed EGFR report metrics."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
METRICS = ROOT / "reports" / "metrics"
OUT = ROOT / "docs" / "assets"


COLORS = {
    "blue": "#2563eb",
    "navy": "#1e3a8a",
    "teal": "#0f766e",
    "green": "#16a34a",
    "orange": "#ea580c",
    "red": "#dc2626",
    "slate": "#475569",
    "gray": "#94a3b8",
}


def load_json(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def setup() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 220,
            "font.size": 11,
            "axes.titlesize": 14,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.22,
            "grid.linewidth": 0.7,
        }
    )


def annotate_bars(ax, bars, fmt="{:.2f}", dy=0.02, fontsize=9) -> None:
    ylim = ax.get_ylim()
    span = ylim[1] - ylim[0]
    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height + dy * span,
            fmt.format(height),
            ha="center",
            va="bottom",
            fontsize=fontsize,
        )


def save(fig, filename: str) -> None:
    fig.tight_layout(pad=2.2)
    fig.savefig(OUT / filename, bbox_inches="tight")
    plt.close(fig)


def build_model_benchmark() -> None:
    status = load_json(METRICS / "final_egfr_project_status.json")
    assay = load_json(METRICS / "egfr_assay_aware_validation_metrics.json")
    gnn = load_json(METRICS / "egfr_gnn_benchmark_metrics.json")

    rf_rows = [
        ("Random split", status["best_random_split_model"]),
        ("Scaffold split", status["best_scaffold_split_model"]),
    ]
    for row in assay["validation_rows"]:
        label = "Assay group" if row["split"] == "assay_group_split" else "Document group"
        rf_rows.append((label, row))

    labels = [label for label, _ in rf_rows]
    rmse = [row["RMSE"] for _, row in rf_rows]
    r2 = [row["R2"] for _, row in rf_rows]
    mae = [row["MAE"] for _, row in rf_rows]

    fig = plt.figure(figsize=(15.5, 7.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.92], hspace=0.38, wspace=0.28)
    ax_rmse = fig.add_subplot(gs[0, 0])
    ax_r2 = fig.add_subplot(gs[0, 1])
    ax_gnn = fig.add_subplot(gs[1, 0])
    ax_note = fig.add_subplot(gs[1, 1])

    colors = [COLORS["green"], COLORS["blue"], COLORS["orange"], COLORS["red"]]
    bars = ax_rmse.bar(labels, rmse, color=colors)
    ax_rmse.set_title("Morgan Random Forest validation error")
    ax_rmse.set_ylabel("RMSE pIC50, lower is better")
    ax_rmse.set_ylim(0, max(rmse) * 1.25)
    ax_rmse.tick_params(axis="x", rotation=18)
    annotate_bars(ax_rmse, bars)

    bars = ax_r2.bar(labels, r2, color=colors)
    ax_r2.set_title("Validation strength drops under harder splits")
    ax_r2.set_ylabel("R2, higher is better")
    ax_r2.set_ylim(0, max(r2) * 1.25)
    ax_r2.tick_params(axis="x", rotation=18)
    annotate_bars(ax_r2, bars)

    x = range(2)
    rf_gnn_labels = ["Random", "Scaffold"]
    rf_vals = [gnn["rf_reference"]["random_split"]["RMSE"], gnn["rf_reference"]["scaffold_split"]["RMSE"]]
    gnn_vals = [gnn["random_split"]["RMSE"], gnn["scaffold_split"]["RMSE"]]
    width = 0.34
    ax_gnn.bar([i - width / 2 for i in x], rf_vals, width=width, color=COLORS["blue"], label="Morgan RF")
    ax_gnn.bar([i + width / 2 for i in x], gnn_vals, width=width, color=COLORS["slate"], label="Dense GCN")
    ax_gnn.set_xticks(list(x), rf_gnn_labels)
    ax_gnn.set_title("Exploratory GCN was retained as benchmark evidence")
    ax_gnn.set_ylabel("RMSE pIC50")
    ax_gnn.set_ylim(0, max(gnn_vals) * 1.25)
    ax_gnn.legend(frameon=False, ncols=2, loc="upper left")

    ax_note.axis("off")
    note = (
        "Selected scorer: Morgan-fingerprint Random Forest\n\n"
        f"Random split: MAE {mae[0]:.3f}, RMSE {rmse[0]:.3f}, R2 {r2[0]:.3f}\n"
        f"Scaffold split: MAE {mae[1]:.3f}, RMSE {rmse[1]:.3f}, R2 {r2[1]:.3f}\n\n"
        "Interpretation:\n"
        "random split is optimistic interpolation;\n"
        "scaffold, assay, and document splits are the\n"
        "chemistry/public-data robustness checks."
    )
    ax_note.text(
        0.02,
        0.96,
        note,
        va="top",
        ha="left",
        fontsize=12,
        linespacing=1.45,
        bbox=dict(boxstyle="round,pad=0.6", facecolor="#f8fafc", edgecolor="#cbd5e1"),
    )

    fig.suptitle("EGFR QSAR model validation summary", fontsize=17, fontweight="bold", y=1.02)
    save(fig, "egfr_model_benchmark.png")


def build_applicability_uncertainty() -> None:
    ad = load_json(METRICS / "applicability_domain_metrics.json")
    conformal = load_json(METRICS / "egfr_conformal_uncertainty_metrics.json")
    pred = pd.read_csv(ROOT / "reports" / "egfr_applicability_domain_predictions.csv")

    order = ["low", "medium", "high"]
    ad_rows = {row["similarity_bin"]: row for row in ad["summary_by_similarity_bin"]}
    mae = [ad_rows[key]["MAE"] for key in order]
    counts = [ad_rows[key]["count"] for key in order]

    fig = plt.figure(figsize=(15.5, 7.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.25, 1.2], wspace=0.34)
    ax_mae = fig.add_subplot(gs[0, 0])
    ax_scatter = fig.add_subplot(gs[0, 1])
    ax_cov = fig.add_subplot(gs[0, 2])

    bars = ax_mae.bar(order, mae, color=[COLORS["red"], COLORS["orange"], COLORS["green"]])
    ax_mae.set_title("Applicability domain matters")
    ax_mae.set_ylabel("MAE pIC50")
    ax_mae.set_ylim(0, max(mae) * 1.28)
    annotate_bars(ax_mae, bars)
    for bar, count in zip(bars, counts):
        ax_mae.text(
            bar.get_x() + bar.get_width() / 2,
            0.03,
            f"n={count:,}",
            ha="center",
            va="bottom",
            color="white",
            fontsize=9,
            fontweight="bold",
        )

    sample = pred.sample(min(3500, len(pred)), random_state=11)
    ax_scatter.scatter(
        sample["max_tanimoto_to_train"],
        sample["absolute_error"],
        s=8,
        alpha=0.22,
        color=COLORS["blue"],
        edgecolors="none",
    )
    ax_scatter.set_title("Prediction error versus nearest-neighbor support")
    ax_scatter.set_xlabel("Max Tanimoto similarity to training set")
    ax_scatter.set_ylabel("Absolute error pIC50")
    ax_scatter.set_ylim(0, min(4.5, max(2.5, pred["absolute_error"].quantile(0.995))))
    ax_scatter.axvline(0.3, color=COLORS["red"], linestyle="--", linewidth=1)
    ax_scatter.axvline(0.7, color=COLORS["green"], linestyle="--", linewidth=1)
    ax_scatter.text(0.305, ax_scatter.get_ylim()[1] * 0.93, "low/medium", color=COLORS["red"], fontsize=9)
    ax_scatter.text(0.705, ax_scatter.get_ylim()[1] * 0.93, "medium/high", color=COLORS["green"], fontsize=9)

    coverage_rows = []
    for split in conformal["splits"]:
        split_label = "Random" if split["split"].startswith("random") else "Scaffold"
        coverage_rows.append((split_label, "overall", split["empirical_coverage"]))
        for row in split["coverage_by_similarity_bin"]:
            coverage_rows.append((split_label, row["similarity_bin"], row["empirical_coverage"]))
    cov = pd.DataFrame(coverage_rows, columns=["split", "bin", "coverage"])
    xlabels = ["overall", "low", "medium", "high"]
    x = range(len(xlabels))
    width = 0.36
    for idx, split in enumerate(["Random", "Scaffold"]):
        values = [float(cov[(cov["split"] == split) & (cov["bin"] == label)]["coverage"].iloc[0]) for label in xlabels]
        offset = -width / 2 if idx == 0 else width / 2
        ax_cov.bar(
            [i + offset for i in x],
            values,
            width=width,
            label=split,
            color=COLORS["teal"] if idx == 0 else COLORS["slate"],
        )
    ax_cov.axhline(0.90, color=COLORS["red"], linestyle="--", linewidth=1.2, label="Target 0.90")
    ax_cov.set_xticks(list(x), xlabels)
    ax_cov.set_ylim(0, 1.08)
    ax_cov.set_ylabel("Empirical interval coverage")
    ax_cov.set_title("Conformal-style uncertainty check")
    ax_cov.legend(frameon=False, fontsize=9)

    fig.suptitle("Applicability-domain and uncertainty checks", fontsize=17, fontweight="bold", y=1.02)
    save(fig, "egfr_applicability_uncertainty.png")


def build_review_structure() -> None:
    triage = load_json(METRICS / "egfr_candidate_triage_metrics.json")
    alerts = load_json(METRICS / "egfr_medchem_alerts_metrics.json")
    redock = load_json(METRICS / "egfr_redocking_metrics.json")
    top5 = pd.read_csv(ROOT / "reports" / "egfr_top5_docking_scores.csv")

    fig = plt.figure(figsize=(15.5, 7.4))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.32)
    ax_triage = fig.add_subplot(gs[0, 0])
    ax_alerts = fig.add_subplot(gs[0, 1])
    ax_dock = fig.add_subplot(gs[1, 0])
    ax_note = fig.add_subplot(gs[1, 1])

    triage_counts = triage["triage_counts"]
    labels = ["prioritize", "review", "review high risk"]
    values = [triage_counts["prioritize"], triage_counts["review"], triage_counts["review_high_risk"]]
    bars = ax_triage.bar(labels, values, color=[COLORS["green"], COLORS["orange"], COLORS["slate"]])
    ax_triage.set_title("Existing-molecule triage categories")
    ax_triage.set_ylabel("Molecules")
    ax_triage.tick_params(axis="x", rotation=12)
    for bar in bars:
        ax_triage.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(values) * 0.025,
            f"{int(bar.get_height()):,}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    alert_labels = ["PAINS", "Brenk", "NIH", "Any alert", "No alert"]
    alert_values = [
        alerts["pains_flagged_fraction"],
        alerts["brenk_flagged_fraction"],
        alerts["nih_flagged_fraction"],
        alerts["combined_medchem_alert_fraction"],
        alerts["no_alert_count"] / alerts["input_row_count"],
    ]
    bars = ax_alerts.barh(alert_labels, alert_values, color=[COLORS["red"], COLORS["orange"], COLORS["slate"], COLORS["blue"], COLORS["green"]])
    ax_alerts.set_xlim(0, 1)
    ax_alerts.set_xlabel("Fraction of model-ready molecules")
    ax_alerts.set_title("Medicinal-chemistry alert review")
    for bar in bars:
        ax_alerts.text(
            bar.get_width() + 0.015,
            bar.get_y() + bar.get_height() / 2,
            f"{bar.get_width() * 100:.1f}%",
            va="center",
            fontsize=9,
        )

    top5 = top5.sort_values("rank_before_docking")
    ids = top5["molecule_id"].str.replace("CHEMBL", "C", regex=False).tolist()
    scores = top5["vina_score_kcal_mol"].astype(float).tolist()
    bars = ax_dock.bar(ids, scores, color=COLORS["teal"])
    ax_dock.set_title("Top-5 docking sanity check")
    ax_dock.set_ylabel("Vina score, kcal/mol")
    ax_dock.axhline(0, color="#0f172a", linewidth=0.8)
    ax_dock.tick_params(axis="x", rotation=20)
    ymin = min(scores) - 0.8
    ax_dock.set_ylim(ymin, 0)
    for bar, score in zip(bars, scores):
        ax_dock.text(
            bar.get_x() + bar.get_width() / 2,
            score - 0.18,
            f"{score:.2f}",
            ha="center",
            va="top",
            color="white",
            fontsize=9,
            fontweight="bold",
        )

    ax_note.axis("off")
    note = (
        "Structure-based sanity check\n\n"
        f"Redocking case: {redock['pdb_id']} / {redock['ligand_id']}\n"
        f"Pose-recovery RMSD: {redock['pose_recovery_rmsd_angstrom']:.3f} A\n"
        f"Top-5 docked: {len(top5)}/5\n"
        f"Score range: {min(scores):.3f} to {max(scores):.3f} kcal/mol\n\n"
        "Docking is used only as retrospective\n"
        "sanity context. It does not confirm\n"
        "binding, potency, or discovery status."
    )
    ax_note.text(
        0.02,
        0.96,
        note,
        va="top",
        ha="left",
        fontsize=12,
        linespacing=1.45,
        bbox=dict(boxstyle="round,pad=0.6", facecolor="#f8fafc", edgecolor="#cbd5e1"),
    )

    fig.suptitle("Ranking, med-chem review, and structure sanity checks", fontsize=17, fontweight="bold", y=1.02)
    save(fig, "egfr_review_and_structure.png")


def main() -> None:
    setup()
    build_model_benchmark()
    build_applicability_uncertainty()
    build_review_structure()


if __name__ == "__main__":
    main()
