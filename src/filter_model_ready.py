"""Create a model-ready EGFR QSAR dataset with broad sanity filters."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


INPUT_PATH = Path("data/processed/egfr_descriptors.csv")
OUTPUT_PATH = Path("data/processed/egfr_model_ready.csv")
REPORT_PATH = Path("reports/model_ready_summary.txt")
FIGURES_DIR = Path("figures")

SUMMARY_COLUMNS = ["median_pIC50", "MolWt", "MolLogP", "TPSA", "HeavyAtomCount", "QED"]


def apply_filter(df: pd.DataFrame, mask: pd.Series, description: str, count_log: list[tuple[str, int]]) -> pd.DataFrame:
    """Apply one filter and print the remaining row count."""
    filtered = df[mask].copy()
    print(f"After {description}: {len(filtered)}")
    count_log.append((description, len(filtered)))
    return filtered


def build_report(df: pd.DataFrame, count_log: list[tuple[str, int]]) -> str:
    """Build a plain-text report documenting the model-ready filter."""
    lines = []

    lines.append("EGFR Model-Ready Dataset Summary")
    lines.append("=" * 32)
    lines.append("")

    lines.append("Filter Row Counts")
    lines.append("-" * 17)
    for description, count in count_log:
        lines.append(f"{description}: {count}")
    lines.append("")

    lines.append("Final Summary Statistics")
    lines.append("-" * 24)
    lines.append(df[SUMMARY_COLUMNS].describe().to_string())
    lines.append("")

    bins = [-float("inf"), 5, 6, 7, 8, float("inf")]
    labels = ["weak (<5)", "low/moderate (5-6)", "moderate (6-7)", "strong (7-8)", "very strong (>8)"]
    activity_bins = pd.cut(df["median_pIC50"], bins=bins, labels=labels, right=False)

    lines.append("Model-Ready pIC50 Activity Bins")
    lines.append("-" * 32)
    lines.append(activity_bins.value_counts(sort=False).to_string())
    lines.append("")

    lines.append("Filtering Rationale")
    lines.append("-" * 19)
    lines.append("median_pIC50 between 3 and 11 removes extreme activity artifacts.")
    lines.append("MolWt between 150 and 900 keeps a broad small-molecule-like range.")
    lines.append("MolLogP between -2 and 8 removes extreme polarity/lipophilicity outliers.")
    lines.append("TPSA <= 250 removes very high-polar-surface-area compounds.")
    lines.append("HeavyAtomCount <= 70 removes unusually large non-small-molecule records.")
    lines.append("")

    return "\n".join(lines)


def save_histogram(df: pd.DataFrame, column: str, xlabel: str, output_path: Path) -> None:
    """Save a simple histogram for one model-ready numeric column."""
    plt.figure(figsize=(7, 5))
    plt.hist(df[column].dropna(), bins=40, edgecolor="black")
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    plt.title(f"Model-Ready {xlabel} Distribution")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def save_scatter(df: pd.DataFrame, x_column: str, y_column: str, xlabel: str, ylabel: str, output_path: Path) -> None:
    """Save a simple scatter plot for model-ready numeric columns."""
    plt.figure(figsize=(7, 5))
    plt.scatter(df[x_column], df[y_column], s=12, alpha=0.5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f"Model-Ready {ylabel} vs {xlabel}")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def save_figures(df: pd.DataFrame) -> None:
    """Save model-ready EDA figures without overwriting full-dataset EDA plots."""
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    save_histogram(df, "median_pIC50", "median pIC50", FIGURES_DIR / "model_ready_pic50_distribution.png")
    save_histogram(df, "MolWt", "Molecular Weight", FIGURES_DIR / "model_ready_molwt_distribution.png")
    save_histogram(df, "MolLogP", "MolLogP", FIGURES_DIR / "model_ready_logp_distribution.png")
    save_histogram(df, "QED", "QED", FIGURES_DIR / "model_ready_qed_distribution.png")
    save_scatter(df, "QED", "median_pIC50", "QED", "median pIC50", FIGURES_DIR / "model_ready_pic50_vs_qed.png")
    save_scatter(
        df,
        "MolLogP",
        "median_pIC50",
        "MolLogP",
        "median pIC50",
        FIGURES_DIR / "model_ready_pic50_vs_logp.png",
    )


def main() -> None:
    """Apply broad activity and chemistry sanity filters before modeling."""
    df = pd.read_csv(INPUT_PATH)
    count_log = [("starting rows", len(df))]
    print(f"Starting rows: {len(df)}")

    df = apply_filter(
        df,
        df["median_pIC50"].between(3, 11, inclusive="both"),
        "keeping 3 <= median_pIC50 <= 11",
        count_log,
    )
    df = apply_filter(
        df,
        df["MolWt"].between(150, 900, inclusive="both"),
        "keeping 150 <= MolWt <= 900",
        count_log,
    )
    df = apply_filter(
        df,
        df["MolLogP"].between(-2, 8, inclusive="both"),
        "keeping -2 <= MolLogP <= 8",
        count_log,
    )
    df = apply_filter(
        df,
        df["TPSA"] <= 250,
        "keeping TPSA <= 250",
        count_log,
    )
    df = apply_filter(
        df,
        df["HeavyAtomCount"] <= 70,
        "keeping HeavyAtomCount <= 70",
        count_log,
    )

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_PATH, index=False)

    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    report_text = build_report(df, count_log)
    REPORT_PATH.write_text(report_text + "\n", encoding="utf-8")
    save_figures(df)

    print(f"Saved model-ready dataset to {OUTPUT_PATH}")
    print(f"Saved model-ready report to {REPORT_PATH}")
    print(f"Saved model-ready figures to {FIGURES_DIR}")
    print("")
    print("Final summary statistics:")
    print(df[SUMMARY_COLUMNS].describe().to_string())


if __name__ == "__main__":
    main()
