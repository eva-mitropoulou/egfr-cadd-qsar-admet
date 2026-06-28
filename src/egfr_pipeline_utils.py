"""Shared utilities for the EGFR CADD and QSAR decision workflow."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
PROCESSED_DIR = DATA_DIR / "processed"
RAW_DIR = DATA_DIR / "raw"
REPORTS_DIR = PROJECT_ROOT / "reports"
METRICS_DIR = REPORTS_DIR / "metrics"
FIGURES_DIR = REPORTS_DIR / "figures"
RUN_LOGS_DIR = REPORTS_DIR / "run_logs"
MODELS_DIR = PROJECT_ROOT / "models"

MODEL_READY_PATH = PROCESSED_DIR / "egfr_model_ready.csv"
RAW_ACTIVITY_PATH = RAW_DIR / "egfr_chembl_ic50_raw.csv"
CLEAN_ACTIVITY_PATH = PROCESSED_DIR / "egfr_ic50_clean.csv"
STANDARDIZED_PATH = PROCESSED_DIR / "egfr_standardized_molecules.csv"

DESCRIPTOR_COLUMNS = [
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

TARGET_COLUMN = "median_pIC50"
RANDOM_STATE = 42
FINGERPRINT_RADIUS = 2
FINGERPRINT_BITS = 2048


def ensure_project_dirs() -> None:
    """Create standard project output directories."""
    for path in [
        PROCESSED_DIR,
        REPORTS_DIR,
        METRICS_DIR,
        FIGURES_DIR,
        RUN_LOGS_DIR,
        MODELS_DIR,
        PROJECT_ROOT / "scripts",
        PROJECT_ROOT / "portfolio_assets",
    ]:
        path.mkdir(parents=True, exist_ok=True)


def save_json(path: Path, payload: dict) -> None:
    """Save pretty JSON with stable formatting."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def read_json(path: Path, default: dict | None = None) -> dict:
    """Read JSON if present, otherwise return a default dictionary."""
    if not path.exists():
        return {} if default is None else default
    return json.loads(path.read_text(encoding="utf-8"))


def write_text(path: Path, content: str) -> None:
    """Write UTF-8 text to a path."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def stable_hash(value: object, length: int = 12) -> str:
    """Return a short stable hash for molecule identifiers or SMILES strings."""
    text = "" if value is None else str(value)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:length]


def molecule_key(row: pd.Series) -> str:
    """Return a safe public-facing molecule key."""
    if "molecule_chembl_id" in row and pd.notna(row["molecule_chembl_id"]):
        return str(row["molecule_chembl_id"])
    if "standardized_smiles" in row and pd.notna(row["standardized_smiles"]):
        return stable_hash(row["standardized_smiles"])
    if "canonical_smiles" in row and pd.notna(row["canonical_smiles"]):
        return stable_hash(row["canonical_smiles"])
    return stable_hash(row.name)


def load_model_ready() -> pd.DataFrame:
    """Load the model-ready EGFR table."""
    if not MODEL_READY_PATH.exists():
        raise FileNotFoundError(f"Missing model-ready dataset: {MODEL_READY_PATH}")
    return pd.read_csv(MODEL_READY_PATH)


def load_standardized_or_model_ready() -> pd.DataFrame:
    """Prefer standardized molecules, falling back to the existing model-ready table."""
    if STANDARDIZED_PATH.exists():
        return pd.read_csv(STANDARDIZED_PATH)
    return load_model_ready()


def markdown_table(df: pd.DataFrame, float_digits: int = 3) -> str:
    """Convert a compact DataFrame to a Markdown table without extra dependencies."""
    display = df.copy()
    for column in display.columns:
        if pd.api.types.is_float_dtype(display[column]):
            display[column] = display[column].map(lambda value: f"{value:.{float_digits}f}")
    columns = list(display.columns)
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join(["---"] * len(columns)) + " |"
    rows = [
        "| " + " | ".join(str(row[column]) for column in columns) + " |"
        for _, row in display.iterrows()
    ]
    return "\n".join([header, separator, *rows])


def value_counts_dict(series: pd.Series, max_items: int = 20) -> dict[str, int]:
    """Return a compact value-count dictionary, including missing values."""
    counts = series.value_counts(dropna=False).head(max_items)
    result: dict[str, int] = {}
    for key, value in counts.items():
        label = "missing" if pd.isna(key) else str(key)
        result[label] = int(value)
    return result


def numeric_summary(df: pd.DataFrame, columns: Iterable[str]) -> dict[str, dict[str, float | int | None]]:
    """Return safe numeric summary statistics for selected columns."""
    summary: dict[str, dict[str, float | int | None]] = {}
    for column in columns:
        if column not in df.columns:
            continue
        values = pd.to_numeric(df[column], errors="coerce")
        summary[column] = {
            "count": int(values.notna().sum()),
            "missing": int(values.isna().sum()),
            "min": _clean_float(values.min()),
            "median": _clean_float(values.median()),
            "mean": _clean_float(values.mean()),
            "max": _clean_float(values.max()),
        }
    return summary


def _clean_float(value: object) -> float | None:
    """Return JSON-safe floats."""
    if value is None or pd.isna(value):
        return None
    as_float = float(value)
    if math.isfinite(as_float):
        return as_float
    return None


def setup_matplotlib() -> None:
    """Configure matplotlib for non-interactive report generation."""
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-egfr")
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)


def scaffold_from_smiles(smiles: str, row_index: int = 0) -> str:
    """Generate a Bemis-Murcko scaffold label from SMILES, with acyclic fallback."""
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return f"invalid_{row_index}"
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    if not scaffold:
        return f"acyclic_{row_index}"
    return scaffold


def add_scaffold_hashes(df: pd.DataFrame, smiles_column: str = "standardized_smiles") -> pd.DataFrame:
    """Add scaffold and scaffold hash columns without exposing scaffolds in reports."""
    if smiles_column not in df.columns:
        smiles_column = "canonical_smiles"
    result = df.copy()
    result["scaffold"] = [
        scaffold_from_smiles(smiles, row_index)
        for row_index, smiles in enumerate(result[smiles_column])
    ]
    result["scaffold_hash"] = result["scaffold"].map(stable_hash)
    return result


def regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    """Calculate standard regression metrics for pIC50 predictions."""
    from scipy.stats import pearsonr, spearmanr
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    pearson = pearsonr(y_true, y_pred).statistic if len(y_true) > 1 else np.nan
    spearman = spearmanr(y_true, y_pred).statistic if len(y_true) > 1 else np.nan

    return {
        "MAE": float(mean_absolute_error(y_true, y_pred)),
        "RMSE": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "R2": float(r2_score(y_true, y_pred)),
        "Pearson": _nan_to_none(pearson),
        "Spearman": _nan_to_none(spearman),
    }


def _nan_to_none(value: float) -> float | None:
    """Convert NaN correlation outputs to None for JSON safety."""
    if value is None:
        return None
    value = float(value)
    if math.isfinite(value):
        return value
    return None


def calculate_morgan_matrix(smiles: Iterable[str], n_bits: int = FINGERPRINT_BITS):
    """Calculate a sparse Morgan fingerprint matrix and RDKit bit vectors."""
    from scipy import sparse
    from rdkit import Chem
    from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

    generator = GetMorganGenerator(radius=FINGERPRINT_RADIUS, fpSize=n_bits)
    rows: list[int] = []
    cols: list[int] = []
    bit_vectors = []
    invalid_count = 0

    for row_idx, smiles_value in enumerate(smiles):
        mol = Chem.MolFromSmiles(str(smiles_value))
        if mol is None:
            invalid_count += 1
            bit_vectors.append(None)
            continue
        bit_vector = generator.GetFingerprint(mol)
        bit_vectors.append(bit_vector)
        on_bits = list(bit_vector.GetOnBits())
        rows.extend([row_idx] * len(on_bits))
        cols.extend(on_bits)

    matrix = sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.uint8), (rows, cols)),
        shape=(len(bit_vectors), n_bits),
        dtype=np.uint8,
    )
    return matrix, bit_vectors, invalid_count


def max_tanimoto_to_train(test_fp, train_fps: list) -> float:
    """Calculate maximum Tanimoto similarity from one fingerprint to training fingerprints."""
    from rdkit import DataStructs

    if test_fp is None or not train_fps:
        return float("nan")
    return float(max(DataStructs.BulkTanimotoSimilarity(test_fp, train_fps)))


def save_figure(path: Path) -> None:
    """Save current matplotlib figure consistently."""
    import matplotlib.pyplot as plt

    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()

