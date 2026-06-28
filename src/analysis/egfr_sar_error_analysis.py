"""SAR-support and interpretable error analysis for EGFR QSAR."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.ensemble import RandomForestRegressor


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import (  # noqa: E402
    DESCRIPTOR_COLUMNS,
    METRICS_DIR,
    MODEL_READY_PATH,
    PROCESSED_DIR,
    RANDOM_STATE,
    REPORTS_DIR,
    stable_hash,
    save_json,
    setup_matplotlib,
    write_text,
)


setup_matplotlib()
import matplotlib.pyplot as plt  # noqa: E402


REPORT_PATH = REPORTS_DIR / "egfr_sar_interpretability_report.md"
METRICS_PATH = METRICS_DIR / "egfr_sar_interpretability_metrics.json"
CLIFFS_PATH = REPORTS_DIR / "egfr_activity_cliffs.csv"
SCAFFOLD_ERROR_PATH = REPORTS_DIR / "egfr_scaffold_error_table.csv"
FIGURES_DIR = REPORTS_DIR / "figures"


def max_tanimoto_pairs(x_matrix, y: np.ndarray, index: pd.DataFrame, batch_size: int = 128) -> pd.DataFrame:
    """Find nearest-neighbor activity cliff candidates without writing structures."""
    counts = np.asarray(x_matrix.sum(axis=1)).ravel().astype(float)
    rows: list[dict] = []
    n_rows = x_matrix.shape[0]
    for start in range(0, n_rows, batch_size):
        stop = min(start + batch_size, n_rows)
        batch = x_matrix[start:stop]
        intersections = (batch @ x_matrix.T).astype(float).toarray()
        batch_counts = counts[start:stop]
        denominators = batch_counts[:, None] + counts[None, :] - intersections
        with np.errstate(divide="ignore", invalid="ignore"):
            similarities = np.divide(intersections, denominators, out=np.zeros_like(intersections), where=denominators > 0)
        for local_i, row_i in enumerate(range(start, stop)):
            similarities[local_i, row_i] = -1.0
            row = similarities[local_i]
            candidate_order = np.argpartition(row, -5)[-5:]
            for row_j in candidate_order:
                sim = float(row[row_j])
                delta = float(abs(y[row_i] - y[row_j]))
                if sim >= 0.85 and delta >= 1.0:
                    first, second = sorted([int(row_i), int(row_j)])
                    rows.append(
                        {
                            "pair_key": stable_hash(f"{first}_{second}"),
                            "molecule_a_id": str(index.iloc[first]["molecule_chembl_id"]),
                            "molecule_b_id": str(index.iloc[second]["molecule_chembl_id"]),
                            "molecule_a_hash": str(index.iloc[first]["molecule_hash"]),
                            "molecule_b_hash": str(index.iloc[second]["molecule_hash"]),
                            "scaffold_a_hash": str(index.iloc[first].get("scaffold_hash", "not_available")),
                            "scaffold_b_hash": str(index.iloc[second].get("scaffold_hash", "not_available")),
                            "tanimoto_similarity": round(sim, 4),
                            "delta_pIC50": round(delta, 4),
                        }
                    )
    if not rows:
        return pd.DataFrame(
            columns=[
                "pair_key",
                "molecule_a_id",
                "molecule_b_id",
                "molecule_a_hash",
                "molecule_b_hash",
                "scaffold_a_hash",
                "scaffold_b_hash",
                "tanimoto_similarity",
                "delta_pIC50",
            ]
        )
    cliffs = pd.DataFrame(rows).drop_duplicates("pair_key")
    return cliffs.sort_values(["delta_pIC50", "tanimoto_similarity"], ascending=[False, False]).head(1000)


def nearest_neighbors_for_ranked(x_matrix, index: pd.DataFrame) -> list[dict]:
    """Summarize nearest existing analogue for top-ranked molecule IDs."""
    ranked_path = REPORTS_DIR / "egfr_ranked_existing_molecules.csv"
    if not ranked_path.exists():
        return []
    ranked = pd.read_csv(ranked_path, usecols=["molecule_chembl_id", "molecule_hash", "model_risk_category"]).head(20)
    id_to_row = {str(mid): row for row, mid in enumerate(index["molecule_chembl_id"].astype(str))}
    counts = np.asarray(x_matrix.sum(axis=1)).ravel().astype(float)
    summaries: list[dict] = []
    for _, ranked_row in ranked.iterrows():
        molecule_id = str(ranked_row["molecule_chembl_id"])
        if molecule_id not in id_to_row:
            continue
        row_i = id_to_row[molecule_id]
        intersections = (x_matrix[row_i] @ x_matrix.T).astype(float).toarray().ravel()
        denominators = counts[row_i] + counts - intersections
        with np.errstate(divide="ignore", invalid="ignore"):
            sims = np.divide(intersections, denominators, out=np.zeros_like(intersections), where=denominators > 0)
        sims[row_i] = -1.0
        neighbor = int(np.argmax(sims))
        summaries.append(
            {
                "molecule_id": molecule_id,
                "molecule_hash": str(ranked_row["molecule_hash"]),
                "nearest_neighbor_id": str(index.iloc[neighbor]["molecule_chembl_id"]),
                "nearest_neighbor_hash": str(index.iloc[neighbor]["molecule_hash"]),
                "nearest_neighbor_similarity": round(float(sims[neighbor]), 4),
                "model_risk_category": str(ranked_row["model_risk_category"]),
            }
        )
    return summaries


def main() -> None:
    """Run SAR-support and error analysis."""
    model_ready = pd.read_csv(MODEL_READY_PATH)
    index = pd.read_csv(PROCESSED_DIR / "features_morgan_index.csv")
    x_fp = sparse.load_npz(PROCESSED_DIR / "features_morgan_fingerprints.npz").astype(np.float32)
    y = index["median_pIC50"].astype(float).to_numpy()

    descriptor_rf = RandomForestRegressor(n_estimators=100, random_state=RANDOM_STATE, n_jobs=-1)
    descriptor_rf.fit(model_ready[DESCRIPTOR_COLUMNS], model_ready["median_pIC50"].astype(float))
    descriptor_importances = pd.DataFrame(
        {
            "descriptor": DESCRIPTOR_COLUMNS,
            "importance": descriptor_rf.feature_importances_,
        }
    ).sort_values("importance", ascending=False)

    morgan_rf = RandomForestRegressor(n_estimators=100, random_state=RANDOM_STATE, n_jobs=-1)
    morgan_rf.fit(x_fp, y)
    morgan_importances = pd.DataFrame(
        {
            "bit_id": np.arange(x_fp.shape[1]),
            "importance": morgan_rf.feature_importances_,
        }
    ).sort_values("importance", ascending=False).head(20)

    cliffs = max_tanimoto_pairs(x_fp, y, index)
    CLIFFS_PATH.parent.mkdir(parents=True, exist_ok=True)
    cliffs.to_csv(CLIFFS_PATH, index=False)

    predictions_path = REPORTS_DIR / "egfr_applicability_domain_predictions.csv"
    pred = pd.read_csv(predictions_path)
    scaffold_errors = (
        pred.groupby("scaffold_hash")
        .agg(
            molecule_count=("absolute_error", "size"),
            MAE=("absolute_error", "mean"),
            RMSE=("absolute_error", lambda values: float(np.sqrt(np.mean(np.square(values))))),
            mean_similarity=("max_tanimoto_to_train", "mean"),
        )
        .reset_index()
    )
    scaffold_errors = scaffold_errors[scaffold_errors["molecule_count"] >= 5].sort_values("MAE", ascending=False)
    scaffold_errors.head(200).to_csv(SCAFFOLD_ERROR_PATH, index=False)

    nearest_neighbor_summary = nearest_neighbors_for_ranked(x_fp, index)

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(7, 4.5))
    top_desc = descriptor_importances.head(9).sort_values("importance")
    plt.barh(top_desc["descriptor"], top_desc["importance"], color="#4C78A8")
    plt.xlabel("Random Forest importance")
    plt.title("Top RDKit Descriptor Importances")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "top_descriptor_importances.png", dpi=200)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    top_scaffolds = scaffold_errors.head(15).sort_values("MAE")
    plt.barh(top_scaffolds["scaffold_hash"], top_scaffolds["MAE"], color="#E15759")
    plt.xlabel("MAE")
    plt.title("Count-Filtered Scaffold-Level Error")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "scaffold_level_error.png", dpi=200)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    if not cliffs.empty:
        plt.scatter(cliffs["tanimoto_similarity"], cliffs["delta_pIC50"], s=18, alpha=0.45)
    plt.xlabel("Tanimoto similarity")
    plt.ylabel("Absolute pIC50 difference")
    plt.title("Activity Cliff Candidates")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "activity_cliff_similarity_vs_delta.png", dpi=200)
    plt.close()

    metrics = {
        "status": "completed",
        "descriptor_importance_status": "completed",
        "morgan_importance_status": "completed_bit_level_only",
        "activity_cliff_count": int(len(cliffs)),
        "activity_cliff_thresholds": {"tanimoto_similarity_min": 0.85, "delta_pIC50_min": 1.0},
        "scaffold_error_rows": int(len(scaffold_errors)),
        "scaffold_error_min_count": 5,
        "top_descriptors": descriptor_importances.head(9).to_dict(orient="records"),
        "top_morgan_bits": morgan_importances.to_dict(orient="records"),
        "nearest_neighbor_summary_count": int(len(nearest_neighbor_summary)),
        "nearest_neighbor_summary": nearest_neighbor_summary,
        "interpretation_scope": "SAR-support analysis using model behavior and chemical-neighborhood checks.",
    }
    save_json(METRICS_PATH, metrics)

    lines = [
        "# EGFR SAR-Support and Interpretable Error Analysis",
        "",
        "This report adds model interpretation and error-analysis evidence for the retrospective EGFR QSAR workflow.",
        "The analysis is SAR-supportive and descriptive for retrospective error analysis.",
        "",
        "## Descriptor Importance",
        "",
    ]
    for _, row in descriptor_importances.head(9).iterrows():
        lines.append(f"- {row['descriptor']}: {row['importance']:.4f}")
    lines.extend(
        [
            "",
            "## Morgan Fingerprint Importance",
            "",
            "Top Morgan fingerprint bits are reported as bit IDs only because reliable substructure reconstruction was not required for this evidence layer.",
        ]
    )
    for _, row in morgan_importances.head(10).iterrows():
        lines.append(f"- bit {int(row['bit_id'])}: {row['importance']:.5f}")
    lines.extend(
        [
            "",
            "## Activity Cliff Candidates",
            "",
            f"- Similarity threshold: >= 0.85",
            f"- Activity-difference threshold: >= 1.0 pIC50",
            f"- Activity cliff candidate pairs saved: {len(cliffs)}",
            f"- Detailed table: `{CLIFFS_PATH.relative_to(PROJECT_ROOT)}`",
            "",
            "## Scaffold-Level Error",
            "",
            f"- Count-filtered scaffold rows: {len(scaffold_errors)}",
            f"- Detailed table: `{SCAFFOLD_ERROR_PATH.relative_to(PROJECT_ROOT)}`",
            "",
            "## Nearest-Neighbor Evidence for Top-Ranked Existing Molecules",
            "",
            f"- Top-ranked molecules summarized: {len(nearest_neighbor_summary)}",
            "- Reports use molecule IDs/hashes and scaffold hashes only.",
            "",
        ]
    )
    write_text(REPORT_PATH, "\n".join(lines))

    print("SAR analysis status: completed")
    print(f"Activity cliff count: {len(cliffs)}")
    print(f"Scaffold error rows: {len(scaffold_errors)}")


if __name__ == "__main__":
    main()
