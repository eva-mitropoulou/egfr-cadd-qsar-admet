"""Rank existing EGFR molecules with model-risk-aware ADMET-style triage."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import (  # noqa: E402
    DESCRIPTOR_COLUMNS,
    FIGURES_DIR,
    METRICS_DIR,
    PROCESSED_DIR,
    REPORTS_DIR,
    load_standardized_or_model_ready,
    markdown_table,
    save_figure,
    save_json,
    setup_matplotlib,
    write_text,
)
from triage.admet_risk_scoring import (  # noqa: E402
    final_triage_category,
    lipinski_violations,
    model_risk_from_similarity,
    property_penalty,
    risk_penalty,
    triage_reason,
    triage_risk_bin,
    uncertainty_penalty,
)


setup_matplotlib()
import matplotlib.pyplot as plt  # noqa: E402


AD_PREDICTIONS_PATH = REPORTS_DIR / "egfr_applicability_domain_predictions.csv"
UNCERTAINTY_PATH = REPORTS_DIR / "egfr_uncertainty_predictions.csv"
MEDCHEM_ALERTS_PATH = PROCESSED_DIR / "egfr_model_ready_with_medchem_alerts.csv"
RANKED_PATH = REPORTS_DIR / "egfr_ranked_existing_molecules.csv"
REPORT_PATH = REPORTS_DIR / "egfr_candidate_triage_report.md"
METRICS_PATH = METRICS_DIR / "egfr_candidate_triage_metrics.json"

ALERT_COLUMNS = [
    "pains_flag",
    "pains_alert_count",
    "brenk_flag",
    "brenk_alert_count",
    "nih_alert_flag",
    "nih_alert_count",
    "unwanted_substructure_flag",
    "unwanted_substructure_count",
    "medchem_alert_flag",
    "medchem_alert_count",
    "medchem_alert_summary",
]


def load_molecule_table_with_alerts() -> tuple[pd.DataFrame, str]:
    """Load molecule properties plus alert annotations when available."""
    if MEDCHEM_ALERTS_PATH.exists():
        return pd.read_csv(MEDCHEM_ALERTS_PATH), "medchem_alerts_available"
    df = load_standardized_or_model_ready()
    for column in ALERT_COLUMNS:
        if column.endswith("_flag"):
            df[column] = False
        elif column.endswith("_summary"):
            df[column] = "none"
        else:
            df[column] = 0
    return df, "medchem_alerts_missing_default_zero"


def main() -> None:
    """Create model-risk-aware ranked table for existing EGFR molecules."""
    df, alert_catalog_status = load_molecule_table_with_alerts()
    ad = pd.read_csv(AD_PREDICTIONS_PATH)
    uncertainty = pd.read_csv(UNCERTAINTY_PATH)

    join_columns = ["molecule_chembl_id"]
    base_columns = join_columns + ["molecule_hash", "scaffold_hash", *DESCRIPTOR_COLUMNS, *ALERT_COLUMNS]
    available_base_columns = [column for column in base_columns if column in df.columns]
    base = df[available_base_columns].copy()
    for column in ALERT_COLUMNS:
        if column not in base.columns:
            if column.endswith("_flag"):
                base[column] = False
            elif column.endswith("_summary"):
                base[column] = "none"
            else:
                base[column] = 0

    ranked = ad[
        join_columns
        + [
            "true_pIC50",
            "predicted_pIC50",
            "absolute_error",
            "max_tanimoto_to_train",
            "similarity_bin",
            "out_of_domain_flag",
        ]
    ].merge(base, on=join_columns, how="inner", validate="one_to_one")

    uncertainty_cols = join_columns + ["rf_prediction_std", "interval_lower_90", "interval_upper_90"]
    ranked = ranked.merge(uncertainty[uncertainty_cols], on=join_columns, how="left", validate="one_to_one")

    ranked["synthetic_accessibility_status"] = "not_available"
    ranked["synthetic_accessibility_score"] = pd.NA
    ranked["lipinski_violations"] = lipinski_violations(ranked)
    ranked["model_risk_category"] = ranked["max_tanimoto_to_train"].map(model_risk_from_similarity)
    ranked["model_risk_penalty"] = ranked["model_risk_category"].map(risk_penalty)
    ranked["uncertainty_penalty"] = ranked["rf_prediction_std"].map(uncertainty_penalty)
    ranked["property_penalty"] = property_penalty(ranked)
    ranked["final_score"] = (
        ranked["predicted_pIC50"]
        + ranked["QED"]
        - ranked["model_risk_penalty"]
        - ranked["uncertainty_penalty"]
        - ranked["property_penalty"]
    )
    ranked["triage_risk_bin"] = ranked.apply(triage_risk_bin, axis=1)
    ranked["triage_reason"] = ranked.apply(triage_reason, axis=1)
    ranked["final_triage_category"] = ranked.apply(final_triage_category, axis=1)

    output_columns = [
        "molecule_chembl_id",
        "molecule_hash",
        "scaffold_hash",
        "predicted_pIC50",
        "max_tanimoto_to_train",
        "similarity_bin",
        "out_of_domain_flag",
        "rf_prediction_std",
        "interval_lower_90",
        "interval_upper_90",
        "model_risk_category",
        "QED",
        "MolWt",
        "MolLogP",
        "TPSA",
        "NumHDonors",
        "NumHAcceptors",
        "NumRotatableBonds",
        "lipinski_violations",
        "pains_flag",
        "pains_alert_count",
        "brenk_flag",
        "brenk_alert_count",
        "unwanted_substructure_flag",
        "unwanted_substructure_count",
        "medchem_alert_flag",
        "medchem_alert_count",
        "medchem_alert_summary",
        "synthetic_accessibility_status",
        "synthetic_accessibility_score",
        "property_penalty",
        "model_risk_penalty",
        "uncertainty_penalty",
        "final_score",
        "triage_risk_bin",
        "triage_reason",
        "final_triage_category",
        "true_pIC50",
        "absolute_error",
    ]
    ranked = ranked[output_columns].sort_values("final_score", ascending=False).reset_index(drop=True)
    ranked.to_csv(RANKED_PATH, index=False)

    diverse_top = ranked.drop_duplicates("scaffold_hash").head(20).copy()
    risk_counts = ranked["model_risk_category"].value_counts().to_dict()
    triage_counts = ranked["final_triage_category"].value_counts().to_dict()
    triage_risk_counts = ranked["triage_risk_bin"].value_counts().to_dict()
    lipinski_counts = ranked["lipinski_violations"].value_counts().sort_index().to_dict()
    top20 = ranked.head(20).copy()

    metrics = {
        "ranked_molecule_count": int(len(ranked)),
        "ranked_table_preserved_all_molecules": int(len(ranked)) == int(len(df)),
        "diverse_top20_unique_scaffolds": int(diverse_top["scaffold_hash"].nunique()),
        "diverse_top20_low_or_medium_risk_count": int(diverse_top["model_risk_category"].isin(["low", "medium"]).sum()),
        "diverse_top20_lipinski_clean_count": int((diverse_top["lipinski_violations"] == 0).sum()),
        "top20_clean_medchem_alert_count": int((~top20["medchem_alert_flag"]).sum()),
        "top20_alert_flagged_count": int(top20["medchem_alert_flag"].sum()),
        "diverse_top20_clean_medchem_alert_count": int((~diverse_top["medchem_alert_flag"]).sum()),
        "diverse_top20_alert_flagged_count": int(diverse_top["medchem_alert_flag"].sum()),
        "risk_counts": {str(key): int(value) for key, value in risk_counts.items()},
        "triage_risk_bin_counts": {str(key): int(value) for key, value in triage_risk_counts.items()},
        "triage_counts": {str(key): int(value) for key, value in triage_counts.items()},
        "lipinski_violation_counts": {str(key): int(value) for key, value in lipinski_counts.items()},
        "pains_flagged_count": int(ranked["pains_flag"].sum()),
        "brenk_flagged_count": int(ranked["brenk_flag"].sum()),
        "unwanted_substructure_flagged_count": int(ranked["unwanted_substructure_flag"].sum()),
        "combined_medchem_alert_count": int(ranked["medchem_alert_flag"].sum()),
        "alert_catalog_status": alert_catalog_status,
        "retrospective_limitation": "Existing-molecule ChEMBL triage with drug-likeness and model-risk proxy rules.",
    }
    save_json(METRICS_PATH, metrics)

    plt.figure(figsize=(7, 4.5))
    funnel_labels = ["ranked", "prioritize", "diverse_top20"]
    funnel_values = [
        len(ranked),
        int((ranked["final_triage_category"] == "prioritize").sum()),
        len(diverse_top),
    ]
    plt.bar(funnel_labels, funnel_values, color="#4C78A8")
    plt.ylabel("Molecule count")
    plt.title("EGFR Triage Funnel")
    save_figure(FIGURES_DIR / "egfr_triage_funnel.png")

    plt.figure(figsize=(7, 4.5))
    categories = ["low", "medium", "high"]
    colors = {"low": "#59A14F", "medium": "#F28E2B", "high": "#E15759"}
    for category in categories:
        subset = ranked[ranked["model_risk_category"] == category]
        plt.scatter(subset["predicted_pIC50"], subset["final_score"], s=10, alpha=0.35, label=category, color=colors[category])
    plt.xlabel("Predicted pIC50")
    plt.ylabel("Final triage score")
    plt.title("Predicted Activity vs Model Risk")
    plt.legend(title="Model risk")
    save_figure(FIGURES_DIR / "predicted_activity_vs_risk.png")

    plt.figure(figsize=(7, 4.5))
    top_ranked = ranked.head(500).copy()
    plt.plot(range(1, len(top_ranked) + 1), top_ranked["max_tanimoto_to_train"])
    plt.xlabel("Rank")
    plt.ylabel("Max Tanimoto to training")
    plt.title("Applicability Domain by Rank")
    save_figure(FIGURES_DIR / "applicability_domain_by_rank.png")

    top_display = diverse_top[
        [
            "molecule_chembl_id",
            "predicted_pIC50",
            "model_risk_category",
            "triage_risk_bin",
            "QED",
            "lipinski_violations",
            "medchem_alert_flag",
            "medchem_alert_count",
            "final_score",
            "final_triage_category",
        ]
    ].copy()
    report = [
        "# EGFR Candidate Triage Report",
        "",
        "This is ADMET-style, drug-likeness, uncertainty, and model-risk triage over existing ChEMBL/project molecules.",
        "It is drug-likeness and model-risk triage over existing molecules.",
        "",
        f"- Ranked molecules: {len(ranked):,}",
        f"- Diverse top-20 unique scaffolds: {metrics['diverse_top20_unique_scaffolds']}",
        f"- Diverse top-20 low/medium model risk: {metrics['diverse_top20_low_or_medium_risk_count']}/20",
        f"- Diverse top-20 Lipinski-clean: {metrics['diverse_top20_lipinski_clean_count']}/20",
        f"- Top-20 with no medicinal-chemistry alert: {metrics['top20_clean_medchem_alert_count']}/20",
        f"- Diverse top-20 with no medicinal-chemistry alert: {metrics['diverse_top20_clean_medchem_alert_count']}/20",
        f"- PAINS-flagged ranked molecules: {metrics['pains_flagged_count']:,}",
        f"- Brenk-flagged ranked molecules: {metrics['brenk_flagged_count']:,}",
        f"- External unwanted-substructure flagged ranked molecules: {metrics['unwanted_substructure_flagged_count']:,}",
        f"- Medicinal-chemistry alert status: {alert_catalog_status}",
        "",
        "PAINS, Brenk, and external unwanted-substructure SMARTS alerts are risk annotations and sensitivity-analysis filters, not automatic exclusions from the primary EGFR QSAR benchmark.",
        "",
        "## Diverse Top 20 Summary",
        "",
        markdown_table(top_display),
        "",
        "Detailed ranked output is saved to `reports/egfr_ranked_existing_molecules.csv`.",
        "",
    ]
    write_text(REPORT_PATH, "\n".join(report))

    print(f"Ranked molecules: {len(ranked)}")
    print(f"Diverse top-20 scaffolds: {metrics['diverse_top20_unique_scaffolds']}")
    print(f"Top-20 clean medchem-alert count: {metrics['top20_clean_medchem_alert_count']}")
    print(f"Metrics: {METRICS_PATH}")


if __name__ == "__main__":
    main()
