"""ADMET-style and model-risk scoring helpers for existing EGFR molecules."""

from __future__ import annotations

import pandas as pd


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
    """Calculate simple drug-likeness property penalty."""
    return (
        0.25 * df["lipinski_violations"]
        + 0.25 * (df["TPSA"] > 140).astype(int)
        + 0.25 * (df["NumRotatableBonds"] > 10).astype(int)
        + 0.50 * (df["QED"] < 0.30).astype(int)
        + 0.50 * df.get("pains_alert", False).astype(int)
        + 0.25 * df.get("brenk_alert", False).astype(int)
    )


def final_triage_category(row: pd.Series) -> str:
    """Assign a compact final triage category."""
    if row["model_risk_category"] == "high" or row["uncertainty_penalty"] >= 0.45:
        return "deprioritize_high_model_risk"
    if row["property_penalty"] >= 1.0:
        return "deprioritize_property_risk"
    if row["predicted_pIC50"] >= 8.0 and row["QED"] >= 0.3 and row["lipinski_violations"] <= 1:
        return "prioritize"
    return "review"

