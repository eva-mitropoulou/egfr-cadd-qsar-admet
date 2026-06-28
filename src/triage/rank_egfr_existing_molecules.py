"""Rank existing EGFR molecules with model-risk-aware ADMET-style triage."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import (  # noqa: E402
    DESCRIPTOR_COLUMNS,
    FIGURES_DIR,
    METRICS_DIR,
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
    uncertainty_penalty,
)


setup_matplotlib()
import matplotlib.pyplot as plt  # noqa: E402


AD_PREDICTIONS_PATH = REPORTS_DIR / "egfr_applicability_domain_predictions.csv"
UNCERTAINTY_PATH = REPORTS_DIR / "egfr_uncertainty_predictions.csv"
RANKED_PATH = REPORTS_DIR / "egfr_ranked_existing_molecules.csv"
REPORT_PATH = REPORTS_DIR / "egfr_candidate_triage_report.md"
METRICS_PATH = METRICS_DIR / "egfr_candidate_triage_metrics.json"


def alert_flags(smiles: pd.Series) -> tuple[list[bool], list[bool], str]:
    """Calculate PAINS/Brenk alerts when RDKit filter catalogs are available."""
    try:
        from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

        pains_params = FilterCatalogParams()
        pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
        pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
        pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
        pains_catalog = FilterCatalog(pains_params)

        brenk_params = FilterCatalogParams()
        brenk_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        brenk_catalog = FilterCatalog(brenk_params)
    except Exception:
        return [False] * len(smiles), [False] * len(smiles), "degraded_filter_catalog_unavailable"

    pains: list[bool] = []
    brenk: list[bool] = []
    for value in smiles:
        mol = Chem.MolFromSmiles(str(value))
        if mol is None:
            pains.append(False)
            brenk.append(False)
            continue
        pains.append(bool(pains_catalog.HasMatch(mol)))
        brenk.append(bool(brenk_catalog.HasMatch(mol)))
    return pains, brenk, "available"


def main() -> None:
    """Create model-risk-aware ranked table for existing EGFR molecules."""
    df = load_standardized_or_model_ready()
    ad = pd.read_csv(AD_PREDICTIONS_PATH)
    uncertainty = pd.read_csv(UNCERTAINTY_PATH)

    join_columns = ["molecule_chembl_id"]
    base_columns = join_columns + ["molecule_hash", "scaffold_hash", *DESCRIPTOR_COLUMNS]
    smiles_column = "standardized_smiles" if "standardized_smiles" in df.columns else "canonical_smiles"
    base = df[base_columns + [smiles_column]].copy()

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

    pains, brenk, catalog_status = alert_flags(ranked[smiles_column])
    ranked["pains_alert"] = pains
    ranked["brenk_alert"] = brenk
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
        "pains_alert",
        "brenk_alert",
        "synthetic_accessibility_status",
        "synthetic_accessibility_score",
        "property_penalty",
        "model_risk_penalty",
        "uncertainty_penalty",
        "final_score",
        "final_triage_category",
        "true_pIC50",
        "absolute_error",
    ]
    ranked = ranked[output_columns].sort_values("final_score", ascending=False).reset_index(drop=True)
    ranked.to_csv(RANKED_PATH, index=False)

    diverse_top = ranked.drop_duplicates("scaffold_hash").head(20).copy()
    risk_counts = ranked["model_risk_category"].value_counts().to_dict()
    triage_counts = ranked["final_triage_category"].value_counts().to_dict()
    lipinski_counts = ranked["lipinski_violations"].value_counts().sort_index().to_dict()

    metrics = {
        "ranked_molecule_count": int(len(ranked)),
        "diverse_top20_unique_scaffolds": int(diverse_top["scaffold_hash"].nunique()),
        "diverse_top20_low_or_medium_risk_count": int(diverse_top["model_risk_category"].isin(["low", "medium"]).sum()),
        "diverse_top20_lipinski_clean_count": int((diverse_top["lipinski_violations"] == 0).sum()),
        "risk_counts": {str(key): int(value) for key, value in risk_counts.items()},
        "triage_counts": {str(key): int(value) for key, value in triage_counts.items()},
        "lipinski_violation_counts": {str(key): int(value) for key, value in lipinski_counts.items()},
        "pains_alert_count": int(ranked["pains_alert"].sum()),
        "brenk_alert_count": int(ranked["brenk_alert"].sum()),
        "alert_catalog_status": catalog_status,
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
            "QED",
            "lipinski_violations",
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
        f"- PAINS/Brenk catalog status: {catalog_status}",
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
    print(f"Metrics: {METRICS_PATH}")


if __name__ == "__main__":
    main()
