"""Rank EGFR compounds with potency, drug-likeness proxies, and model risk."""

from pathlib import Path

import pandas as pd


MODEL_READY_PATH = Path("data/processed/egfr_model_ready.csv")
PREDICTIONS_PATH = Path("results/applicability_domain_predictions.csv")
VALIDATION_RANKED_PATH = Path("results/ranked_candidates_with_validation.csv")
PORTFOLIO_RANKED_PATH = Path("results/ranked_candidates_portfolio.csv")
TOP_20_PATH = Path("results/top_20_candidates.csv")
TOP_20_DIVERSE_PATH = Path("results/top_20_diverse_candidates.csv")

JOIN_COLUMNS = ["molecule_chembl_id", "canonical_smiles"]

VALIDATION_OUTPUT_COLUMNS = [
    "molecule_chembl_id",
    "canonical_smiles",
    "scaffold",
    "predicted_pIC50",
    "true_pIC50",
    "absolute_error",
    "max_tanimoto_to_train",
    "model_risk_category",
    "model_risk_penalty",
    "QED",
    "MolWt",
    "MolLogP",
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
    "lipinski_violations",
    "flag_tpsa_gt_140",
    "flag_rotatable_bonds_gt_10",
    "flag_qed_lt_030",
    "property_penalty",
    "final_score",
]

PORTFOLIO_OUTPUT_COLUMNS = [
    column
    for column in VALIDATION_OUTPUT_COLUMNS
    if column not in {"true_pIC50", "absolute_error"}
]


def model_risk_category(max_tanimoto: float) -> str:
    """Categorize model risk from nearest training-set Tanimoto similarity."""
    if max_tanimoto > 0.7:
        return "low"
    if max_tanimoto >= 0.5:
        return "medium"
    if max_tanimoto >= 0.3:
        return "high"
    return "very_high"


def model_risk_penalty(category: str) -> float:
    """Map model-risk category to ranking penalty."""
    penalties = {
        "low": 0.0,
        "medium": 0.25,
        "high": 0.50,
        "very_high": 1.00,
    }
    return penalties[category]


def add_lipinski_and_property_flags(df: pd.DataFrame) -> pd.DataFrame:
    """Add simple drug-likeness proxy flags and penalties."""
    df = df.copy()

    df["lipinski_molwt_violation"] = df["MolWt"] > 500
    df["lipinski_logp_violation"] = df["MolLogP"] > 5
    df["lipinski_hbd_violation"] = df["NumHDonors"] > 5
    df["lipinski_hba_violation"] = df["NumHAcceptors"] > 10

    lipinski_columns = [
        "lipinski_molwt_violation",
        "lipinski_logp_violation",
        "lipinski_hbd_violation",
        "lipinski_hba_violation",
    ]
    df["lipinski_violations"] = df[lipinski_columns].sum(axis=1)

    df["flag_tpsa_gt_140"] = df["TPSA"] > 140
    df["flag_rotatable_bonds_gt_10"] = df["NumRotatableBonds"] > 10
    df["flag_qed_lt_030"] = df["QED"] < 0.30

    df["property_penalty"] = (
        0.25 * df["lipinski_violations"]
        + 0.25 * df["flag_tpsa_gt_140"].astype(int)
        + 0.25 * df["flag_rotatable_bonds_gt_10"].astype(int)
        + 0.50 * df["flag_qed_lt_030"].astype(int)
    )

    return df


def build_ranked_table(model_ready: pd.DataFrame, predictions: pd.DataFrame) -> pd.DataFrame:
    """Join descriptor data with out-of-fold predictions and calculate ranking score."""
    columns_from_model_ready = JOIN_COLUMNS + [
        "MolWt",
        "MolLogP",
        "TPSA",
        "NumHDonors",
        "NumHAcceptors",
        "NumRotatableBonds",
        "QED",
    ]
    columns_from_predictions = JOIN_COLUMNS + [
        "scaffold",
        "true_pIC50",
        "predicted_pIC50",
        "absolute_error",
        "max_tanimoto_to_train",
    ]

    ranked = predictions[columns_from_predictions].merge(
        model_ready[columns_from_model_ready],
        on=JOIN_COLUMNS,
        how="inner",
        validate="one_to_one",
    )

    ranked["model_risk_category"] = ranked["max_tanimoto_to_train"].apply(model_risk_category)
    ranked["model_risk_penalty"] = ranked["model_risk_category"].apply(model_risk_penalty)
    ranked = add_lipinski_and_property_flags(ranked)

    ranked["final_score"] = (
        ranked["predicted_pIC50"]
        + ranked["QED"]
        - ranked["model_risk_penalty"]
        - ranked["property_penalty"]
    )

    ranked = ranked[VALIDATION_OUTPUT_COLUMNS].sort_values("final_score", ascending=False).reset_index(drop=True)

    return ranked


def make_portfolio_table(validation_ranked: pd.DataFrame) -> pd.DataFrame:
    """Remove retrospective validation columns for prospective-style ranking."""
    return validation_ranked[PORTFOLIO_OUTPUT_COLUMNS].copy()


def select_diverse_top_candidates(portfolio_ranked: pd.DataFrame, n_candidates: int = 20) -> pd.DataFrame:
    """Select top candidates while keeping only one molecule per scaffold."""
    diverse = portfolio_ranked.drop_duplicates(subset="scaffold", keep="first")
    return diverse.head(n_candidates).copy()


def print_summary(ranked: pd.DataFrame, top_20: pd.DataFrame) -> None:
    """Print a short ranking summary."""
    print(f"Total ranked molecules: {len(ranked)}")
    print(
        "Top 20 predicted_pIC50 range: "
        f"{top_20['predicted_pIC50'].min():.3f} to {top_20['predicted_pIC50'].max():.3f}"
    )
    print("")
    print("Model risk category counts:")
    print(ranked["model_risk_category"].value_counts().to_string())
    print("")
    print("Lipinski violation counts:")
    print(ranked["lipinski_violations"].value_counts().sort_index().to_string())


def main() -> None:
    """Create ranked EGFR compound tables for triage."""
    model_ready = pd.read_csv(MODEL_READY_PATH)
    predictions = pd.read_csv(PREDICTIONS_PATH)

    print(f"Model-ready rows: {len(model_ready)}")
    print(f"Applicability-domain prediction rows: {len(predictions)}")

    validation_ranked = build_ranked_table(model_ready, predictions)
    portfolio_ranked = make_portfolio_table(validation_ranked)
    top_20 = portfolio_ranked.head(20).copy()
    top_20_diverse = select_diverse_top_candidates(portfolio_ranked, n_candidates=20)

    VALIDATION_RANKED_PATH.parent.mkdir(parents=True, exist_ok=True)
    validation_ranked.to_csv(VALIDATION_RANKED_PATH, index=False)
    portfolio_ranked.to_csv(PORTFOLIO_RANKED_PATH, index=False)
    top_20.to_csv(TOP_20_PATH, index=False)
    top_20_diverse.to_csv(TOP_20_DIVERSE_PATH, index=False)

    print_summary(portfolio_ranked, top_20)
    print("")
    print("Saved ADMET-style/model-risk triage tables:")
    print(f"- {VALIDATION_RANKED_PATH}")
    print(f"- {PORTFOLIO_RANKED_PATH}")
    print(f"- {TOP_20_PATH}")
    print(f"- {TOP_20_DIVERSE_PATH}")
    print("")
    print("Note: this is drug-likeness triage, not true ADMET prediction or candidate nomination.")


if __name__ == "__main__":
    main()
