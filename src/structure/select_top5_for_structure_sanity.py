"""Select top clean, diverse ranked EGFR molecules for structure sanity docking."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import METRICS_DIR, REPORTS_DIR, save_json, write_text  # noqa: E402


RANKED_PATH = REPORTS_DIR / "egfr_ranked_existing_molecules.csv"
MOLECULE_TABLE_PATH = PROJECT_ROOT / "data" / "processed" / "egfr_model_ready_with_medchem_alerts.csv"
SELECTION_PATH = REPORTS_DIR / "egfr_top5_structure_selection.csv"
METRICS_PATH = METRICS_DIR / "egfr_top5_structure_selection_metrics.json"


def bool_series(df: pd.DataFrame, column: str, default: bool = False) -> pd.Series:
    """Return a robust boolean series."""
    if column not in df.columns:
        return pd.Series(default, index=df.index)
    return df[column].fillna(default).astype(bool)


def interval_width(df: pd.DataFrame) -> pd.Series:
    """Return conformal interval width when available."""
    lower = pd.to_numeric(df.get("interval_lower_90"), errors="coerce")
    upper = pd.to_numeric(df.get("interval_upper_90"), errors="coerce")
    return upper - lower


def select_unique_scaffolds(candidates: pd.DataFrame, needed: int) -> list[int]:
    """Select ranked rows while enforcing unique scaffold hashes when possible."""
    selected: list[int] = []
    seen_scaffolds: set[str] = set()
    for idx, row in candidates.iterrows():
        scaffold = str(row.get("scaffold_hash", "missing"))
        if scaffold in seen_scaffolds:
            continue
        selected.append(idx)
        seen_scaffolds.add(scaffold)
        if len(selected) == needed:
            break
    return selected


def selection_rounds(df: pd.DataFrame) -> list[tuple[str, pd.Series]]:
    """Return strict-to-relaxed selection masks."""
    clean_alerts = (
        ~bool_series(df, "pains_flag")
        & ~bool_series(df, "brenk_flag")
        & ~bool_series(df, "unwanted_substructure_flag")
        & ~bool_series(df, "medchem_alert_flag")
    )
    in_domain = ~bool_series(df, "out_of_domain_flag")
    low_model_risk = df.get("model_risk_category", "").astype(str).isin(["low"])
    low_or_medium_model_risk = df.get("model_risk_category", "").astype(str).isin(["low", "medium"])
    acceptable_lipinski = pd.to_numeric(df.get("lipinski_violations"), errors="coerce").fillna(99) <= 1
    acceptable_qed = pd.to_numeric(df.get("QED"), errors="coerce").fillna(0) >= 0.30
    low_uncertainty = pd.to_numeric(df.get("rf_prediction_std"), errors="coerce").fillna(999) <= pd.to_numeric(
        df.get("rf_prediction_std"), errors="coerce"
    ).median()

    return [
        (
            "strict_clean_unique_scaffold",
            clean_alerts & in_domain & low_model_risk & low_uncertainty & acceptable_lipinski & acceptable_qed,
        ),
        (
            "relaxed_medium_model_risk_no_alerts",
            clean_alerts & in_domain & low_or_medium_model_risk & acceptable_lipinski & acceptable_qed,
        ),
        (
            "relaxed_same_scaffold_if_needed_no_alerts",
            clean_alerts & in_domain & acceptable_lipinski & acceptable_qed,
        ),
        (
            "relaxed_uncertainty_no_alerts",
            clean_alerts & acceptable_lipinski & acceptable_qed,
        ),
    ]


def main() -> None:
    """Select five existing ranked molecules for top-5 structure sanity docking."""
    if not RANKED_PATH.exists():
        raise FileNotFoundError(f"Missing ranked table: {RANKED_PATH}")
    if not MOLECULE_TABLE_PATH.exists():
        raise FileNotFoundError(f"Missing molecule table: {MOLECULE_TABLE_PATH}")

    ranked = pd.read_csv(RANKED_PATH)
    molecule_table = pd.read_csv(MOLECULE_TABLE_PATH, usecols=["molecule_chembl_id"])
    ranked = ranked[ranked["molecule_chembl_id"].isin(set(molecule_table["molecule_chembl_id"]))].copy()
    ranked["rank_before_docking"] = range(1, len(ranked) + 1)
    ranked["conformal_interval_width"] = interval_width(ranked)

    sort_columns = [column for column in ["final_score", "predicted_pIC50"] if column in ranked.columns]
    ranked = ranked.sort_values(sort_columns, ascending=False).reset_index(drop=True)
    ranked["rank_before_docking"] = range(1, len(ranked) + 1)

    selected_indices: list[int] = []
    selected_reasons: dict[int, str] = {}
    relaxation_used: list[str] = []
    for reason, mask in selection_rounds(ranked):
        candidates = ranked.loc[mask].copy()
        candidates = candidates.sort_values(sort_columns, ascending=False)
        unique_indices = select_unique_scaffolds(candidates, 5 - len(selected_indices))
        for idx in unique_indices:
            if idx not in selected_indices:
                selected_indices.append(idx)
                selected_reasons[idx] = reason
        if len(selected_indices) >= 5:
            relaxation_used.append(reason)
            break
        relaxation_used.append(reason)

    if len(selected_indices) < 5:
        clean_alerts = (
            ~bool_series(ranked, "pains_flag")
            & ~bool_series(ranked, "brenk_flag")
            & ~bool_series(ranked, "unwanted_substructure_flag")
            & ~bool_series(ranked, "medchem_alert_flag")
        )
        fallback = ranked.loc[clean_alerts].sort_values(sort_columns, ascending=False)
        for idx in fallback.index:
            if idx in selected_indices:
                continue
            selected_indices.append(idx)
            selected_reasons[idx] = "fallback_ranked_clean_same_scaffold_allowed"
            if len(selected_indices) >= 5:
                break

    selected = ranked.loc[selected_indices[:5]].copy()
    selected["selection_reason"] = [selected_reasons.get(idx, "selected") for idx in selected.index]
    selected["molecule_id"] = selected["molecule_chembl_id"].astype(str)
    selected["scaffold_id"] = selected["scaffold_hash"].astype(str)
    selected["applicability_domain_bin"] = selected.get("similarity_bin", "unavailable")

    output_columns = [
        "molecule_id",
        "rank_before_docking",
        "scaffold_id",
        "predicted_pIC50",
        "conformal_interval_width",
        "applicability_domain_bin",
        "pains_flag",
        "brenk_flag",
        "unwanted_substructure_flag",
        "medchem_alert_flag",
        "selection_reason",
    ]
    SELECTION_PATH.parent.mkdir(parents=True, exist_ok=True)
    selected[output_columns].to_csv(SELECTION_PATH, index=False)

    metrics = {
        "selection_status": "completed" if len(selected) == 5 else "degraded_less_than_5_selected",
        "ranked_table": str(RANKED_PATH.relative_to(PROJECT_ROOT)),
        "molecule_table": str(MOLECULE_TABLE_PATH.relative_to(PROJECT_ROOT)),
        "ranked_row_count": int(len(ranked)),
        "selected_count": int(len(selected)),
        "unique_scaffold_count": int(selected["scaffold_id"].nunique()) if not selected.empty else 0,
        "relaxation_rounds_attempted": relaxation_used,
        "medchem_alert_clean_selected_count": int((~selected["medchem_alert_flag"].fillna(False).astype(bool)).sum())
        if not selected.empty
        else 0,
        "selection_output": str(SELECTION_PATH.relative_to(PROJECT_ROOT)),
        "raw_smiles_in_public_selection_table": False,
    }
    save_json(METRICS_PATH, metrics)
    write_text(
        REPORTS_DIR / "egfr_top5_structure_selection_note.md",
        "\n".join(
            [
                "# EGFR Top-5 Structure Selection Note",
                "",
                "Selected molecules are existing ranked EGFR records. Public selection output uses molecule IDs and scaffold IDs only.",
                "",
                f"- Selected molecules: {metrics['selected_count']}",
                f"- Unique scaffolds: {metrics['unique_scaffold_count']}",
                f"- Medchem-alert-clean selected molecules: {metrics['medchem_alert_clean_selected_count']}",
                "",
            ]
        ),
    )

    print(f"Selection output: {SELECTION_PATH}")
    print(f"Selected molecule count: {metrics['selected_count']}")
    print(f"Unique scaffold count: {metrics['unique_scaffold_count']}")


if __name__ == "__main__":
    main()
