"""Create final portfolio summary assets for the EGFR QSAR/CADD project."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


MODEL_READY_PATH = Path("data/processed/egfr_model_ready.csv")
DESCRIPTOR_METRICS_PATH = Path("results/descriptor_baseline_metrics.csv")
FINGERPRINT_METRICS_PATH = Path("results/fingerprint_baseline_metrics.csv")
COMBINED_METRICS_PATH = Path("results/combined_baseline_metrics.csv")
SCAFFOLD_METRICS_PATH = Path("results/scaffold_fingerprint_metrics.csv")
CROSS_VALIDATION_PATH = Path("results/cross_validation_metrics.csv")
APPLICABILITY_DOMAIN_PATH = Path("results/applicability_domain_summary.csv")
TOP_20_PATH = Path("results/top_20_candidates.csv")
TOP_20_DIVERSE_PATH = Path("results/top_20_diverse_candidates.csv")

SUMMARY_REPORT_PATH = Path("reports/project_results_summary.md")
TOP_20_DIVERSE_REPORT_PATH = Path("reports/top_20_diverse_candidates.md")
MODEL_COMPARISON_PLOT_PATH = Path("figures/model_performance_comparison.png")
VALIDATION_COMPARISON_PLOT_PATH = Path("figures/random_vs_scaffold_validation.png")
APPLICABILITY_DOMAIN_PLOT_PATH = Path("figures/applicability_domain_bins.png")

RF_MODEL_NAME = "Random forest"


def format_float(value: float, digits: int = 3) -> str:
    """Format numeric values for compact markdown tables."""
    return f"{value:.{digits}f}"


def markdown_table(df: pd.DataFrame, float_digits: int = 3) -> str:
    """Convert a small DataFrame to a GitHub-flavored Markdown table."""
    display = df.copy()

    for column in display.columns:
        if pd.api.types.is_float_dtype(display[column]):
            display[column] = display[column].map(lambda value: format_float(value, float_digits))

    columns = list(display.columns)
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join(["---"] * len(columns)) + " |"
    rows = [
        "| " + " | ".join(str(row[column]) for column in columns) + " |"
        for _, row in display.iterrows()
    ]

    return "\n".join([header, separator, *rows])


def read_random_forest_row(path: Path) -> pd.Series:
    """Read one metrics file and return the Random Forest result row."""
    metrics = pd.read_csv(path)
    random_forest_rows = metrics[metrics["model"] == RF_MODEL_NAME]

    if random_forest_rows.empty:
        raise ValueError(f"Could not find {RF_MODEL_NAME!r} in {path}")

    return random_forest_rows.iloc[0]


def build_feature_comparison() -> pd.DataFrame:
    """Collect descriptor, fingerprint, and combined Random Forest metrics."""
    rows = []
    for feature_set, path in [
        ("RDKit descriptors", DESCRIPTOR_METRICS_PATH),
        ("Morgan fingerprints", FINGERPRINT_METRICS_PATH),
        ("Fingerprints + descriptors", COMBINED_METRICS_PATH),
    ]:
        row = read_random_forest_row(path)
        rows.append(
            {
                "feature_set": feature_set,
                "model": RF_MODEL_NAME,
                "MAE": row["MAE"],
                "RMSE": row["RMSE"],
                "R2": row["R2"],
            }
        )

    return pd.DataFrame(rows)


def build_random_vs_scaffold_comparison() -> pd.DataFrame:
    """Collect random-split and scaffold-split Random Forest metrics."""
    random_row = read_random_forest_row(FINGERPRINT_METRICS_PATH)
    scaffold_row = read_random_forest_row(SCAFFOLD_METRICS_PATH)

    return pd.DataFrame(
        [
            {
                "validation": "Random split",
                "feature_set": "Morgan fingerprints",
                "MAE": random_row["MAE"],
                "RMSE": random_row["RMSE"],
                "R2": random_row["R2"],
            },
            {
                "validation": "Scaffold split",
                "feature_set": "Morgan fingerprints",
                "MAE": scaffold_row["MAE"],
                "RMSE": scaffold_row["RMSE"],
                "R2": scaffold_row["R2"],
            },
        ]
    )


def build_cross_validation_summary() -> pd.DataFrame:
    """Summarize Random Forest cross-validation as mean plus standard deviation."""
    cv = pd.read_csv(CROSS_VALIDATION_PATH)
    rf_cv = cv[cv["model"] == RF_MODEL_NAME].copy()

    summary = (
        rf_cv.groupby("validation_scheme")
        .agg(
            MAE_mean=("MAE", "mean"),
            MAE_std=("MAE", "std"),
            RMSE_mean=("RMSE", "mean"),
            RMSE_std=("RMSE", "std"),
            R2_mean=("R2", "mean"),
            R2_std=("R2", "std"),
        )
        .reset_index()
    )

    labels = {
        "random_kfold": "Random KFold",
        "scaffold_groupkfold": "Scaffold GroupKFold",
    }
    summary["validation_scheme"] = summary["validation_scheme"].replace(labels)

    return summary


def build_candidate_summary(
    top_20: pd.DataFrame,
    top_20_diverse: pd.DataFrame,
    ranked_count: int,
) -> dict[str, object]:
    """Calculate compact summary values for the final ranked candidates."""
    return {
        "ranked_count": ranked_count,
        "top20_pred_min": top_20["predicted_pIC50"].min(),
        "top20_pred_max": top_20["predicted_pIC50"].max(),
        "diverse_top20_pred_min": top_20_diverse["predicted_pIC50"].min(),
        "diverse_top20_pred_max": top_20_diverse["predicted_pIC50"].max(),
        "unique_scaffolds": top_20_diverse["scaffold"].nunique(),
        "risk_counts": top_20_diverse["model_risk_category"].value_counts().to_dict(),
        "lipinski_clean_count": int((top_20_diverse["lipinski_violations"] == 0).sum()),
    }


def save_metric_bar_plot(df: pd.DataFrame, label_column: str, title: str, output_path: Path) -> None:
    """Save a simple three-panel bar chart for MAE, RMSE, and R2."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    metrics = ["MAE", "RMSE", "R2"]

    for axis, metric in zip(axes, metrics):
        axis.bar(df[label_column], df[metric], color="#4C78A8")
        axis.set_title(metric)
        axis.tick_params(axis="x", labelrotation=25)
        axis.grid(axis="y", alpha=0.25)

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def save_applicability_domain_plot(summary: pd.DataFrame) -> None:
    """Save a bar plot of error by max-Tanimoto similarity bin."""
    APPLICABILITY_DOMAIN_PLOT_PATH.parent.mkdir(parents=True, exist_ok=True)

    fig, axis = plt.subplots(figsize=(7, 4.5))
    axis.bar(summary["similarity_bin"], summary["MAE"], color="#59A14F")
    axis.set_xlabel("Max Tanimoto similarity to training set")
    axis.set_ylabel("MAE in pIC50")
    axis.set_title("Applicability Domain: Error by Similarity Bin")
    axis.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(APPLICABILITY_DOMAIN_PLOT_PATH, dpi=200)
    plt.close(fig)


def create_top_20_report(top_20_diverse: pd.DataFrame) -> str:
    """Create a markdown table for the diverse top 20 candidates."""
    columns = [
        "molecule_chembl_id",
        "predicted_pIC50",
        "max_tanimoto_to_train",
        "model_risk_category",
        "QED",
        "MolWt",
        "MolLogP",
        "TPSA",
        "lipinski_violations",
        "final_score",
    ]

    report_lines = [
        "# Top 20 Diverse EGFR Candidates",
        "",
        "These are portfolio-style ranked compounds selected with one molecule per scaffold.",
        "The ranking uses predicted potency, drug-likeness proxies, and model-risk penalties.",
        "",
        markdown_table(top_20_diverse[columns]),
        "",
        "Note: this is retrospective ChEMBL triage, not prospective candidate nomination.",
        "",
    ]

    return "\n".join(report_lines)


def create_project_summary_report(
    model_ready: pd.DataFrame,
    feature_comparison: pd.DataFrame,
    validation_comparison: pd.DataFrame,
    cv_summary: pd.DataFrame,
    applicability_summary: pd.DataFrame,
    top_20_diverse: pd.DataFrame,
    candidate_summary: dict[str, object],
) -> str:
    """Create the final project summary markdown report."""
    cv_display = cv_summary.copy()
    cv_display["MAE"] = cv_display.apply(
        lambda row: f"{row['MAE_mean']:.3f} +/- {row['MAE_std']:.3f}",
        axis=1,
    )
    cv_display["RMSE"] = cv_display.apply(
        lambda row: f"{row['RMSE_mean']:.3f} +/- {row['RMSE_std']:.3f}",
        axis=1,
    )
    cv_display["R2"] = cv_display.apply(
        lambda row: f"{row['R2_mean']:.3f} +/- {row['R2_std']:.3f}",
        axis=1,
    )
    cv_display = cv_display[["validation_scheme", "MAE", "RMSE", "R2"]]

    applicability_display = applicability_summary.copy()
    applicability_display["count"] = applicability_display["count"].astype(int)

    risk_counts = candidate_summary["risk_counts"]
    low_risk = risk_counts.get("low", 0)
    medium_risk = risk_counts.get("medium", 0)

    report_lines = [
        "# EGFR QSAR/CADD Project Results Summary",
        "",
        "## Dataset",
        "",
        f"- Model-ready molecules: {len(model_ready):,}",
        "- Activity endpoint: ChEMBL EGFR IC50 converted to median pIC50 per molecule.",
        "- Broad sanity filters removed extreme activity/property artifacts while retaining most molecules.",
        "",
        "## Feature-Set Comparison",
        "",
        markdown_table(feature_comparison),
        "",
        "Morgan fingerprints captured EGFR-relevant substructure information better than broad descriptors alone.",
        "",
        "## Random vs Scaffold Split",
        "",
        markdown_table(validation_comparison),
        "",
        "Random split was optimistic; scaffold split gave a harder estimate of generalization to new chemotypes.",
        "",
        "## Cross-Validation",
        "",
        markdown_table(cv_display),
        "",
        "Model performance was stable under random CV but dropped under scaffold-aware CV.",
        "",
        "## Applicability Domain",
        "",
        markdown_table(applicability_display),
        "",
        "Prediction error decreased as maximum Tanimoto similarity to the training chemistry increased.",
        "",
        "## Candidate Ranking",
        "",
        f"- Ranked molecules: {candidate_summary['ranked_count']:,}",
        (
            "- Top 20 predicted pIC50 range: "
            f"{candidate_summary['top20_pred_min']:.3f}-{candidate_summary['top20_pred_max']:.3f}"
        ),
        (
            "- Diverse top 20 predicted pIC50 range: "
            f"{candidate_summary['diverse_top20_pred_min']:.3f}-"
            f"{candidate_summary['diverse_top20_pred_max']:.3f}"
        ),
        f"- Diverse top 20 unique scaffolds: {candidate_summary['unique_scaffolds']}",
        f"- Diverse top 20 model risk: {low_risk} low-risk, {medium_risk} medium-risk",
        f"- Lipinski-clean diverse top 20: {candidate_summary['lipinski_clean_count']}/20",
        "",
        "Top diverse candidates are saved in `reports/top_20_diverse_candidates.md`.",
        "",
        "## Limitations",
        "",
        "- This is a retrospective ChEMBL portfolio project, not prospective drug discovery.",
        "- IC50 values come from heterogeneous assays and publications.",
        "- The ADMET layer is proxy/drug-likeness triage, not true ADMET prediction.",
        "- Docking was not included.",
        "- Predictions are less reliable outside the model applicability domain.",
        "",
        "## Generated Figures",
        "",
        "- `figures/model_performance_comparison.png`",
        "- `figures/random_vs_scaffold_validation.png`",
        "- `figures/applicability_domain_bins.png`",
        "",
    ]

    return "\n".join(report_lines)


def main() -> None:
    """Create final markdown reports and summary figures."""
    model_ready = pd.read_csv(MODEL_READY_PATH)
    top_20 = pd.read_csv(TOP_20_PATH)
    top_20_diverse = pd.read_csv(TOP_20_DIVERSE_PATH)
    applicability_summary = pd.read_csv(APPLICABILITY_DOMAIN_PATH)

    feature_comparison = build_feature_comparison()
    validation_comparison = build_random_vs_scaffold_comparison()
    cv_summary = build_cross_validation_summary()
    candidate_summary = build_candidate_summary(top_20, top_20_diverse, ranked_count=len(model_ready))

    save_metric_bar_plot(
        feature_comparison,
        label_column="feature_set",
        title="Random-Split Random Forest Performance",
        output_path=MODEL_COMPARISON_PLOT_PATH,
    )
    save_metric_bar_plot(
        validation_comparison,
        label_column="validation",
        title="Random vs Scaffold Validation",
        output_path=VALIDATION_COMPARISON_PLOT_PATH,
    )
    save_applicability_domain_plot(applicability_summary)

    SUMMARY_REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_REPORT_PATH.write_text(
        create_project_summary_report(
            model_ready,
            feature_comparison,
            validation_comparison,
            cv_summary,
            applicability_summary,
            top_20_diverse,
            candidate_summary,
        ),
        encoding="utf-8",
    )
    TOP_20_DIVERSE_REPORT_PATH.write_text(create_top_20_report(top_20_diverse), encoding="utf-8")

    print(f"Saved {SUMMARY_REPORT_PATH}")
    print(f"Saved {TOP_20_DIVERSE_REPORT_PATH}")
    print(f"Saved {MODEL_COMPARISON_PLOT_PATH}")
    print(f"Saved {VALIDATION_COMPARISON_PLOT_PATH}")
    print(f"Saved {APPLICABILITY_DOMAIN_PLOT_PATH}")


if __name__ == "__main__":
    main()
