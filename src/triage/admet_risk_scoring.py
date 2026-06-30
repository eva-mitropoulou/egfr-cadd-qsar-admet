"""ADMET-style and model-risk scoring helpers for existing EGFR molecules."""

from __future__ import annotations

import pandas as pd


def _bool_column(df: pd.DataFrame, column: str) -> pd.Series:
    """Return a boolean column with a False default."""
    if column not in df.columns:
        return pd.Series(False, index=df.index)
    return df[column].fillna(False).astype(bool)


def lipinski_violations(df: pd.DataFrame) -> pd.Series:
    """Count simple Lipinski rule-of-five violations."""
    return (
        (df["MolWt"] > 500).astype(int)
        + (df["MolLogP"] > 5).astype(int)
        + (df["NumHDonors"] > 5).astype(int)
        + (df["NumHAcceptors"] > 10).astype(int)
    )


def model_risk_from_similarity(value: float) -> str:
    """Convert nearest-neighbor similarity to model-risk category."""
    if value > 0.7:
        return "low"
    if value >= 0.3:
        return "medium"
    return "high"


def risk_penalty(category: str) -> float:
    """Return model-risk penalty."""
    return {"low": 0.0, "medium": 0.35, "high": 0.9}.get(category, 0.9)


def uncertainty_penalty(std_value: float) -> float:
    """Convert RF prediction standard deviation to a smooth penalty."""
    if pd.isna(std_value):
        return 0.25
    if std_value < 0.35:
        return 0.0
    if std_value < 0.60:
        return 0.20
    return 0.45


def property_penalty(df: pd.DataFrame) -> pd.Series:
    """Calculate simple drug-likeness and alert-risk property penalty."""
    return (
        0.25 * df["lipinski_violations"]
        + 0.25 * (df["TPSA"] > 140).astype(int)
        + 0.25 * (df["NumRotatableBonds"] > 10).astype(int)
        + 0.50 * (df["QED"] < 0.30).astype(int)
        + 0.75 * _bool_column(df, "pains_flag").astype(int)
        + 0.50 * _bool_column(df, "brenk_flag").astype(int)
        + 0.50 * _bool_column(df, "unwanted_substructure_flag").astype(int)
    )


def triage_reasons(row: pd.Series) -> list[str]:
    """Return compact triage reasons without implying causality."""
    reasons: list[str] = []
    if bool(row.get("pains_flag", False)):
        reasons.append("PAINS alert annotation")
    if bool(row.get("brenk_flag", False)):
        reasons.append("Brenk alert annotation")
    if bool(row.get("unwanted_substructure_flag", False)):
        reasons.append("external unwanted-substructure annotation")
    if bool(row.get("out_of_domain_flag", False)) or row.get("model_risk_category") == "high":
        reasons.append("outside/low applicability-domain similarity")
    if row.get("uncertainty_penalty", 0) >= 0.45:
        reasons.append("higher uncertainty proxy")
    if row.get("lipinski_violations", 0) > 1:
        reasons.append("multiple Lipinski violations")
    if row.get("TPSA", 0) > 140:
        reasons.append("high TPSA")
    if row.get("NumRotatableBonds", 0) > 10:
        reasons.append("high rotatable-bond count")
    if row.get("QED", 1) < 0.30:
        reasons.append("low QED")
    if not reasons:
        reasons.append("no major triage flags")
    return reasons


def triage_risk_bin(row: pd.Series) -> str:
    """Assign a transparent low/medium/high triage risk bin."""
    alert_flag = any(
        bool(row.get(column, False))
        for column in ["pains_flag", "brenk_flag", "unwanted_substructure_flag"]
    )
    high_model_risk = bool(row.get("out_of_domain_flag", False)) or row.get("model_risk_category") == "high"
    high_uncertainty = row.get("uncertainty_penalty", 0) >= 0.45
    liability_count = sum(
        [
            alert_flag,
            high_model_risk,
            high_uncertainty,
            row.get("lipinski_violations", 0) > 1,
            row.get("TPSA", 0) > 140,
            row.get("NumRotatableBonds", 0) > 10,
            row.get("QED", 1) < 0.30,
        ]
    )
    if alert_flag or high_model_risk or high_uncertainty or liability_count >= 2:
        return "high"
    if liability_count == 1 or row.get("lipinski_violations", 0) == 1:
        return "medium"
    return "low"


def triage_reason(row: pd.Series) -> str:
    """Return semicolon-separated risk reasons."""
    return "; ".join(triage_reasons(row))


def final_triage_category(row: pd.Series) -> str:
    """Assign a compact final triage category."""
    if row.get("triage_risk_bin") == "high":
        return "review_high_risk"
    if row["model_risk_category"] == "high" or row["uncertainty_penalty"] >= 0.45:
        return "deprioritize_high_model_risk"
    if row["property_penalty"] >= 1.0:
        return "deprioritize_property_risk"
    if row["predicted_pIC50"] >= 8.0 and row["QED"] >= 0.3 and row["lipinski_violations"] <= 1:
        return "prioritize"
    return "review"
