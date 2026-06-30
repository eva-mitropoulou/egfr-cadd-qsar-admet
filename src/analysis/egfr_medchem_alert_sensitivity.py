"""Sensitivity analyses for medicinal-chemistry alert annotations.

The main EGFR QSAR benchmark remains the full model-ready dataset. This script
asks how PAINS/Brenk/unwanted-substructure alert filters affect dataset
composition and matched Morgan RF baseline performance.
"""

from __future__ import annotations

import sys
from pathlib import Path
import os

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import (  # noqa: E402
    FIGURES_DIR,
    METRICS_DIR,
    RANDOM_STATE,
    REPORTS_DIR,
    TARGET_COLUMN,
    add_scaffold_hashes,
    calculate_morgan_matrix,
    markdown_table,
    regression_metrics,
    save_figure,
    save_json,
    setup_matplotlib,
    write_text,
)


setup_matplotlib()
import matplotlib.pyplot as plt  # noqa: E402


ANNOTATED_PATH = PROJECT_ROOT / "data" / "processed" / "egfr_model_ready_with_medchem_alerts.csv"
FEATURE_MATRIX_PATH = PROJECT_ROOT / "data" / "processed" / "features_morgan_fingerprints.npz"
FEATURE_INDEX_PATH = PROJECT_ROOT / "data" / "processed" / "features_morgan_index.csv"
RANKED_PATH = REPORTS_DIR / "egfr_ranked_existing_molecules.csv"
REPORT_PATH = REPORTS_DIR / "egfr_medchem_alert_sensitivity_report.md"
METRICS_PATH = METRICS_DIR / "egfr_medchem_alert_sensitivity_metrics.json"
SENSITIVITY_RF_TREES = int(os.environ.get("EGFR_MEDCHEM_SENSITIVITY_TREES", "30"))


def subset_definitions(df: pd.DataFrame) -> dict[str, pd.Series]:
    """Return masks for full and alert-filtered sensitivity subsets."""
    return {
        "main_all_model_ready": pd.Series(True, index=df.index),
        "pains_excluded": ~df["pains_flag"].fillna(False).astype(bool),
        "pains_brenk_excluded": ~df["pains_flag"].fillna(False).astype(bool)
        & ~df["brenk_flag"].fillna(False).astype(bool),
        "strict_medchem_clean": ~df["medchem_alert_flag"].fillna(False).astype(bool),
    }


def assign_scaffold_split(df: pd.DataFrame, test_fraction: float = 0.2) -> tuple[np.ndarray, np.ndarray]:
    """Create a deterministic scaffold-disjoint split by molecule count."""
    if "scaffold_hash" not in df.columns:
        df = add_scaffold_hashes(df)
    group_sizes = df.groupby("scaffold_hash", sort=False).size().sort_values(ascending=False)
    test_scaffolds: set[str] = set()
    test_count = 0
    target_count = int(round(len(df) * test_fraction))
    for scaffold, size in group_sizes.items():
        if test_count >= target_count and test_scaffolds:
            break
        test_scaffolds.add(str(scaffold))
        test_count += int(size)
    test_mask = df["scaffold_hash"].astype(str).isin(test_scaffolds).to_numpy()
    indices = np.arange(len(df))
    return indices[~test_mask], indices[test_mask]


def load_or_build_features(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, int, str]:
    """Load aligned Morgan fingerprints, falling back to local calculation."""
    if FEATURE_MATRIX_PATH.exists() and FEATURE_INDEX_PATH.exists():
        index = pd.read_csv(FEATURE_INDEX_PATH)
        if "molecule_chembl_id" in index.columns and list(index["molecule_chembl_id"]) == list(df["molecule_chembl_id"]):
            matrix = sparse.load_npz(FEATURE_MATRIX_PATH)
            return matrix.toarray().astype(np.uint8, copy=False), df[TARGET_COLUMN].astype(float).to_numpy(), 0, "loaded_existing_feature_matrix"

    smiles_column = "standardized_smiles" if "standardized_smiles" in df.columns else "canonical_smiles"
    matrix, _, invalid_count = calculate_morgan_matrix(df[smiles_column])
    return matrix.toarray().astype(np.uint8, copy=False), df[TARGET_COLUMN].astype(float).to_numpy(), int(invalid_count), "calculated_from_smiles"


def fit_evaluate(df: pd.DataFrame, x: np.ndarray, y: np.ndarray, split: str, feature_invalid_count: int) -> dict[str, object]:
    """Fit/evaluate the standard Morgan RF baseline for one split."""
    if split == "random_split":
        train_idx, test_idx = train_test_split(
            np.arange(len(df)),
            test_size=0.2,
            random_state=RANDOM_STATE,
        )
        train_scaffold_count = None
        test_scaffold_count = None
        overlap_count = None
    elif split == "scaffold_split":
        train_idx, test_idx = assign_scaffold_split(df)
        train_scaffolds = set(df.iloc[train_idx]["scaffold_hash"].astype(str))
        test_scaffolds = set(df.iloc[test_idx]["scaffold_hash"].astype(str))
        train_scaffold_count = len(train_scaffolds)
        test_scaffold_count = len(test_scaffolds)
        overlap_count = len(train_scaffolds & test_scaffolds)
    else:
        raise ValueError(f"Unsupported split: {split}")

    model = RandomForestRegressor(n_estimators=SENSITIVITY_RF_TREES, random_state=RANDOM_STATE, n_jobs=-1)
    model.fit(x[train_idx], y[train_idx])
    pred = model.predict(x[test_idx])
    metrics = regression_metrics(y[test_idx], pred)
    return {
        "split": split,
        "model": "Random Forest Morgan",
        "train_size": int(len(train_idx)),
        "test_size": int(len(test_idx)),
        "train_scaffold_count": train_scaffold_count,
        "test_scaffold_count": test_scaffold_count,
        "scaffold_overlap_count": overlap_count,
        "invalid_molecule_count": int(feature_invalid_count),
        "n_estimators": int(SENSITIVITY_RF_TREES),
        **metrics,
    }


def subset_summary(name: str, subset: pd.DataFrame, total_count: int, full_mean: float) -> dict[str, object]:
    """Summarize one alert-filtered subset."""
    values = subset[TARGET_COLUMN].astype(float)
    return {
        "subset": name,
        "row_count": int(len(subset)),
        "fraction_retained": float(len(subset) / total_count) if total_count else 0.0,
        "molecule_count": int(subset["molecule_chembl_id"].nunique()) if "molecule_chembl_id" in subset.columns else int(len(subset)),
        "scaffold_count": int(subset["scaffold_hash"].nunique()) if "scaffold_hash" in subset.columns else None,
        "pIC50_mean": float(values.mean()),
        "pIC50_median": float(values.median()),
        "pIC50_std": float(values.std()),
        "pIC50_mean_shift_vs_full": float(values.mean() - full_mean),
    }


def top20_composition(ranked: pd.DataFrame | None) -> dict[str, object]:
    """Summarize medchem-alert composition for the top-ranked molecules."""
    if ranked is None or ranked.empty:
        return {"status": "ranked_table_unavailable"}
    top20 = ranked.head(20).copy()
    diverse = ranked.drop_duplicates("scaffold_hash").head(20).copy() if "scaffold_hash" in ranked.columns else top20
    return {
        "status": "available",
        "top20_count": int(len(top20)),
        "top20_clean_count": int((~top20["medchem_alert_flag"].fillna(False).astype(bool)).sum()),
        "top20_alert_flagged_count": int(top20["medchem_alert_flag"].fillna(False).astype(bool).sum()),
        "top20_pains_flagged_count": int(top20["pains_flag"].fillna(False).astype(bool).sum()),
        "top20_brenk_flagged_count": int(top20["brenk_flag"].fillna(False).astype(bool).sum()),
        "top20_unwanted_flagged_count": int(top20["unwanted_substructure_flag"].fillna(False).astype(bool).sum()),
        "diverse_top20_clean_count": int((~diverse["medchem_alert_flag"].fillna(False).astype(bool)).sum()),
        "diverse_top20_alert_flagged_count": int(diverse["medchem_alert_flag"].fillna(False).astype(bool).sum()),
    }


def save_plots(df: pd.DataFrame, summaries: list[dict[str, object]], metrics_rows: list[dict[str, object]], composition: dict[str, object]) -> None:
    """Save sensitivity figures."""
    subsets = subset_definitions(df)

    plt.figure(figsize=(8, 5))
    for name, mask in subsets.items():
        subset = df.loc[mask, TARGET_COLUMN].astype(float)
        if subset.empty:
            continue
        plt.hist(subset, bins=30, alpha=0.35, label=name)
    plt.xlabel("median pIC50")
    plt.ylabel("Molecule count")
    plt.title("Activity Distribution By Medchem-Alert Sensitivity Subset")
    plt.legend(fontsize=8)
    save_figure(FIGURES_DIR / "medchem_alert_activity_distribution.png")

    metrics_df = pd.DataFrame(metrics_rows)
    scaffold = metrics_df[metrics_df["split"] == "scaffold_split"].copy()
    if not scaffold.empty:
        plt.figure(figsize=(8, 5))
        plt.bar(scaffold["subset"], scaffold["RMSE"], color="#4C78A8")
        plt.ylabel("Scaffold-split RMSE")
        plt.title("Morgan RF Sensitivity To Alert Exclusions")
        plt.xticks(rotation=25, ha="right")
        save_figure(FIGURES_DIR / "medchem_alert_subset_model_performance.png")

    if composition.get("status") == "available":
        plt.figure(figsize=(6, 4))
        labels = ["top20 clean", "top20 alert", "diverse clean", "diverse alert"]
        values = [
            composition["top20_clean_count"],
            composition["top20_alert_flagged_count"],
            composition["diverse_top20_clean_count"],
            composition["diverse_top20_alert_flagged_count"],
        ]
        plt.bar(labels, values, color=["#59A14F", "#E15759", "#59A14F", "#E15759"])
        plt.ylabel("Molecule count")
        plt.title("Top-Ranked Alert Composition")
        plt.xticks(rotation=20, ha="right")
        save_figure(FIGURES_DIR / "top_ranked_medchem_alert_composition.png")


def main() -> None:
    """Run medchem-alert sensitivity analyses."""
    if not ANNOTATED_PATH.exists():
        raise FileNotFoundError(f"Missing annotated alert table: {ANNOTATED_PATH}")

    df = pd.read_csv(ANNOTATED_PATH)
    if "scaffold_hash" not in df.columns:
        smiles_column = "standardized_smiles" if "standardized_smiles" in df.columns else "canonical_smiles"
        df = add_scaffold_hashes(df, smiles_column=smiles_column)

    total_count = len(df)
    full_mean = float(df[TARGET_COLUMN].astype(float).mean())
    full_x, full_y, feature_invalid_count, feature_source = load_or_build_features(df)
    subset_masks = subset_definitions(df)
    summaries: list[dict[str, object]] = []
    metrics_rows: list[dict[str, object]] = []
    availability: dict[str, str] = {}

    for name, mask in subset_masks.items():
        subset = df.loc[mask].reset_index(drop=True)
        summaries.append(subset_summary(name, subset, total_count, full_mean))
        if len(subset) < 500 or subset["scaffold_hash"].nunique() < 20:
            availability[name] = "skipped_too_small_for_stable_random_and_scaffold_evaluation"
            continue
        availability[name] = "evaluated"
        subset_indices = np.flatnonzero(mask.to_numpy())
        subset_x = full_x[subset_indices]
        subset_y = full_y[subset_indices]
        for split in ["random_split", "scaffold_split"]:
            row = fit_evaluate(subset, subset_x, subset_y, split, feature_invalid_count)
            row["subset"] = name
            row["row_count"] = int(len(subset))
            row["fraction_retained"] = float(len(subset) / total_count)
            metrics_rows.append(row)

    ranked = pd.read_csv(RANKED_PATH) if RANKED_PATH.exists() else None
    composition = top20_composition(ranked)
    save_plots(df, summaries, metrics_rows, composition)

    metrics = {
        "input_path": str(ANNOTATED_PATH.relative_to(PROJECT_ROOT)),
        "input_row_count": int(total_count),
        "feature_source": feature_source,
        "sensitivity_model_note": (
            "Sensitivity analysis uses the same Morgan fingerprint features, split logic, and random state as the primary RF setup, "
            f"with {SENSITIVITY_RF_TREES} trees for runtime. The official 100-tree primary benchmark is preserved unchanged."
        ),
        "subset_summaries": summaries,
        "performance_rows": metrics_rows,
        "subset_availability": availability,
        "top20_candidate_composition": composition,
        "primary_benchmark_preserved": True,
        "interpretation": "Sensitivity analysis only; alert-containing molecules are not removed from the primary EGFR QSAR benchmark by default.",
    }
    save_json(METRICS_PATH, metrics)

    summary_df = pd.DataFrame(summaries)
    perf_df = pd.DataFrame(metrics_rows)
    if not perf_df.empty:
        perf_display = perf_df[
            [
                "subset",
                "split",
                "train_size",
                "test_size",
                "MAE",
                "RMSE",
                "R2",
                "Pearson",
                "Spearman",
            ]
        ].copy()
    else:
        perf_display = pd.DataFrame()

    full_scaffold = perf_df[(perf_df["subset"] == "main_all_model_ready") & (perf_df["split"] == "scaffold_split")]
    strict_scaffold = perf_df[(perf_df["subset"] == "strict_medchem_clean") & (perf_df["split"] == "scaffold_split")]
    if not full_scaffold.empty and not strict_scaffold.empty:
        rmse_delta = float(strict_scaffold["RMSE"].iloc[0] - full_scaffold["RMSE"].iloc[0])
        conclusion = (
            "The strict medchem-clean subset changed scaffold-split RMSE by "
            f"{rmse_delta:.3f} pIC50 units relative to the full sensitivity baseline."
        )
    else:
        conclusion = "Strict medchem-clean scaffold-split performance was unavailable or skipped."

    report = [
        "# EGFR Medchem-Alert Sensitivity Report",
        "",
        "PAINS, Brenk, and external unwanted-substructure SMARTS alerts were used as medicinal-chemistry risk annotations and sensitivity-analysis filters, not automatic exclusions from the primary EGFR QSAR benchmark.",
        "",
        f"This sensitivity run uses Morgan fingerprints and the same split logic/random state as the primary Random Forest baseline, with {SENSITIVITY_RF_TREES} trees for runtime. The official 100-tree primary benchmark files are preserved unchanged.",
        "",
        "## Subset Composition",
        "",
        markdown_table(summary_df),
        "",
        "## Morgan RF Sensitivity Metrics",
        "",
        markdown_table(perf_display) if not perf_display.empty else "Sensitivity model metrics were unavailable.",
        "",
        "## Top-Ranked Alert Composition",
        "",
        f"- Top-20 clean count: {composition.get('top20_clean_count', 'unavailable')}",
        f"- Top-20 alert-flagged count: {composition.get('top20_alert_flagged_count', 'unavailable')}",
        f"- Diverse top-20 clean count: {composition.get('diverse_top20_clean_count', 'unavailable')}",
        f"- Diverse top-20 alert-flagged count: {composition.get('diverse_top20_alert_flagged_count', 'unavailable')}",
        "",
        "## Interpretation",
        "",
        conclusion,
        "This is a sensitivity analysis, not a replacement for the primary full-dataset benchmark.",
        "",
    ]
    write_text(REPORT_PATH, "\n".join(report))

    print(f"Sensitivity report: {REPORT_PATH}")
    print(f"Input rows: {total_count}")
    print(f"Subsets attempted: {len(subset_masks)}")
    print(f"Performance rows: {len(metrics_rows)}")
    print(f"Metrics: {METRICS_PATH}")


if __name__ == "__main__":
    main()
