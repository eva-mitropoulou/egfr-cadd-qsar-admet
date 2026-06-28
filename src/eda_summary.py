"""Generate EDA sanity-check summaries and plots for EGFR descriptors."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


INPUT_PATH = Path("data/processed/egfr_descriptors.csv")
REPORT_PATH = Path("reports/eda_summary.txt")
FIGURES_DIR = Path("figures")

NUMERIC_COLUMNS = [
    "median_pIC50",
    "median_IC50_nM",
    "n_measurements",
    "MolWt",
    "MolLogP",
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
    "RingCount",
    "HeavyAtomCount",
    "QED",
]


def build_report(df: pd.DataFrame) -> str:
    """Build a plain-text EDA summary for the descriptor dataset."""
    lines = []

    lines.append("EGFR Descriptor EDA Summary")
    lines.append("=" * 28)
    lines.append("")
    lines.append(f"Rows: {df.shape[0]}")
    lines.append(f"Columns: {df.shape[1]}")
    lines.append("")

    lines.append("Missing Values Per Column")
    lines.append("-" * 25)
    lines.append(df.isna().sum().to_string())
    lines.append("")

    lines.append("Summary Statistics")
    lines.append("-" * 18)
    lines.append(df[NUMERIC_COLUMNS].describe().to_string())
    lines.append("")

    bins = [-float("inf"), 5, 6, 7, 8, float("inf")]
    labels = ["weak (<5)", "low/moderate (5-6)", "moderate (6-7)", "strong (7-8)", "very strong (>8)"]
    activity_bins = pd.cut(df["median_pIC50"], bins=bins, labels=labels, right=False)

    lines.append("pIC50 Activity Bins")
    lines.append("-" * 19)
    lines.append(activity_bins.value_counts(sort=False).to_string())
    lines.append("")

    display_columns = [
        "molecule_chembl_id",
        "canonical_smiles",
        "median_pIC50",
        "median_IC50_nM",
        "n_measurements",
        "MolWt",
        "MolLogP",
        "QED",
    ]

    lines.append("Top 10 Strongest Compounds")
    lines.append("-" * 26)
    strongest = df.sort_values("median_pIC50", ascending=False).head(10)
    lines.append(strongest[display_columns].to_string(index=False))
    lines.append("")

    lines.append("Top 10 Weakest Compounds")
    lines.append("-" * 24)
    weakest = df.sort_values("median_pIC50", ascending=True).head(10)
    lines.append(weakest[display_columns].to_string(index=False))
    lines.append("")

    return "\n".join(lines)


def save_histogram(df: pd.DataFrame, column: str, xlabel: str, output_path: Path) -> None:
    """Save a simple histogram for one numeric column."""
    plt.figure(figsize=(7, 5))
    plt.hist(df[column].dropna(), bins=40, edgecolor="black")
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    plt.title(f"{xlabel} Distribution")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def save_scatter(df: pd.DataFrame, x_column: str, y_column: str, xlabel: str, ylabel: str, output_path: Path) -> None:
    """Save a simple scatter plot for two numeric columns."""
    plt.figure(figsize=(7, 5))
    plt.scatter(df[x_column], df[y_column], s=12, alpha=0.5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f"{ylabel} vs {xlabel}")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def main() -> None:
    """Run the EDA sanity check and save report/figures."""
    df = pd.read_csv(INPUT_PATH)

    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    report_text = build_report(df)
    REPORT_PATH.write_text(report_text + "\n", encoding="utf-8")

    save_histogram(df, "median_pIC50", "median pIC50", FIGURES_DIR / "pic50_distribution.png")
    save_histogram(df, "MolWt", "Molecular Weight", FIGURES_DIR / "molwt_distribution.png")
    save_histogram(df, "MolLogP", "MolLogP", FIGURES_DIR / "logp_distribution.png")
    save_histogram(df, "QED", "QED", FIGURES_DIR / "qed_distribution.png")

    save_scatter(df, "QED", "median_pIC50", "QED", "median pIC50", FIGURES_DIR / "pic50_vs_qed.png")
    save_scatter(
        df,
        "MolLogP",
        "median_pIC50",
        "MolLogP",
        "median pIC50",
        FIGURES_DIR / "pic50_vs_logp.png",
    )

    print(report_text)
    print(f"Saved EDA report to {REPORT_PATH}")
    print(f"Saved figures to {FIGURES_DIR}")


if __name__ == "__main__":
    main()
