"""Annotate EGFR molecules with medicinal-chemistry alert flags.

Alerts are annotations for triage and sensitivity analysis. They are not used as
automatic deletions from the primary EGFR QSAR benchmark.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterable

import pandas as pd
from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import (  # noqa: E402
    METRICS_DIR,
    MODEL_READY_PATH,
    PROCESSED_DIR,
    REPORTS_DIR,
    STANDARDIZED_PATH,
    save_json,
    stable_hash,
    write_text,
)


OUTPUT_PATH = PROCESSED_DIR / "egfr_model_ready_with_medchem_alerts.csv"
REPORT_PATH = REPORTS_DIR / "egfr_medchem_alerts_report.md"
METRICS_PATH = METRICS_DIR / "egfr_medchem_alerts_metrics.json"
SCHEMA_PATH = PROJECT_ROOT / "data" / "reference" / "unwanted_substructures_schema.md"

NAME_COLUMNS = ["name", "alert", "description", "rule", "substructure"]
SMARTS_COLUMNS = ["smarts", "SMARTS", "pattern"]
EXTERNAL_CANDIDATES = [
    PROJECT_ROOT / "data" / "reference" / "unwanted_substructures.csv",
    PROJECT_ROOT / "data" / "raw" / "unwanted_substructures.csv",
    PROJECT_ROOT / "data" / "processed" / "unwanted_substructures.csv",
    PROJECT_ROOT / "unwanted_substructures.csv",
]


def load_input_table() -> tuple[pd.DataFrame, Path]:
    """Load standardized molecules when present, otherwise model-ready molecules."""
    if STANDARDIZED_PATH.exists():
        return pd.read_csv(STANDARDIZED_PATH), STANDARDIZED_PATH
    if MODEL_READY_PATH.exists():
        return pd.read_csv(MODEL_READY_PATH), MODEL_READY_PATH
    raise FileNotFoundError("No model-ready EGFR molecule table was found.")


def detect_id_column(df: pd.DataFrame) -> str:
    """Detect or create a stable molecule identifier column."""
    if "molecule_chembl_id" in df.columns:
        return "molecule_chembl_id"
    for column in df.columns:
        lowered = column.lower()
        if "molecule" in lowered and "id" in lowered:
            return column
    if "molecule_hash" in df.columns:
        return "molecule_hash"
    df["molecule_hash"] = [
        stable_hash(value)
        for value in df.get("standardized_smiles", df.get("canonical_smiles", pd.Series(df.index)))
    ]
    return "molecule_hash"


def detect_smiles_column(df: pd.DataFrame) -> str:
    """Detect the canonical/standardized molecule representation column."""
    for column in ["standardized_smiles", "canonical_smiles", "smiles", "SMILES"]:
        if column in df.columns:
            return column
    for column in df.columns:
        if "smiles" in column.lower():
            return column
    raise ValueError("No SMILES-like column was found in the molecule table.")


def build_filter_catalog(catalog_names: Iterable[str]) -> FilterCatalog | None:
    """Build an RDKit filter catalog from stable catalog names."""
    params = FilterCatalogParams()
    added = False
    catalogs = FilterCatalogParams.FilterCatalogs
    for name in catalog_names:
        if not hasattr(catalogs, name):
            continue
        params.AddCatalog(getattr(catalogs, name))
        added = True
    if not added:
        return None
    return FilterCatalog(params)


def catalog_alert_names(mol: Chem.Mol | None, catalog: FilterCatalog | None) -> list[str]:
    """Return unique RDKit FilterCatalog descriptions for one molecule."""
    if mol is None or catalog is None:
        return []
    if not catalog.HasMatch(mol):
        return []
    names: set[str] = set()
    for match in catalog.GetMatches(mol):
        try:
            names.add(str(match.GetDescription()))
        except Exception:
            names.add("unnamed_filter_catalog_alert")
    return sorted(names)


def find_external_unwanted_csv() -> Path | None:
    """Find an optional external unwanted-substructure SMARTS catalog."""
    candidates = [path for path in EXTERNAL_CANDIDATES if path.exists()]
    candidates.extend(sorted((PROJECT_ROOT / "data").glob("**/unwanted_substructures.csv")))
    unique: list[Path] = []
    seen: set[Path] = set()
    for path in candidates:
        resolved = path.resolve()
        if resolved not in seen:
            unique.append(path)
            seen.add(resolved)
    return unique[0] if unique else None


def create_external_schema_note() -> None:
    """Write expected external SMARTS CSV schema when the catalog is absent."""
    write_text(
        SCHEMA_PATH,
        "\n".join(
            [
                "# External Unwanted-Substructure SMARTS Schema",
                "",
                "Optional file name: `unwanted_substructures.csv`",
                "",
                "Expected columns:",
                "",
                "- `name`: human-readable alert name",
                "- `smarts`: SMARTS pattern",
                "- `description`: optional note",
                "",
                "Accepted aliases include `alert`, `rule`, or `substructure` for the name field and `SMARTS` or `pattern` for the SMARTS field.",
                "",
                "These SMARTS are used as medicinal-chemistry risk annotations and sensitivity-analysis filters, not automatic exclusions from the primary EGFR QSAR benchmark.",
                "",
            ]
        ),
    )


def load_external_patterns(path: Path | None) -> tuple[list[tuple[str, Chem.Mol]], dict[str, object]]:
    """Load and validate optional unwanted-substructure SMARTS patterns."""
    if path is None:
        create_external_schema_note()
        return [], {
            "external_csv_found": False,
            "external_csv_path": None,
            "external_smarts_total_count": 0,
            "external_smarts_valid_count": 0,
            "external_smarts_invalid_count": 0,
            "external_schema_path": str(SCHEMA_PATH.relative_to(PROJECT_ROOT)),
        }

    catalog = pd.read_csv(path)
    name_column = next((column for column in NAME_COLUMNS if column in catalog.columns), None)
    smarts_column = next((column for column in SMARTS_COLUMNS if column in catalog.columns), None)
    if smarts_column is None:
        return [], {
            "external_csv_found": True,
            "external_csv_path": str(path.relative_to(PROJECT_ROOT)),
            "external_smarts_total_count": int(len(catalog)),
            "external_smarts_valid_count": 0,
            "external_smarts_invalid_count": int(len(catalog)),
            "external_csv_error": "No SMARTS-like column found.",
        }

    patterns: list[tuple[str, Chem.Mol]] = []
    invalid_count = 0
    for idx, row in catalog.iterrows():
        smarts = row.get(smarts_column)
        pattern = Chem.MolFromSmarts("" if pd.isna(smarts) else str(smarts))
        if pattern is None:
            invalid_count += 1
            continue
        if name_column is not None and pd.notna(row.get(name_column)):
            name = str(row.get(name_column))
        else:
            name = f"external_unwanted_{idx + 1}"
        patterns.append((name, pattern))

    return patterns, {
        "external_csv_found": True,
        "external_csv_path": str(path.relative_to(PROJECT_ROOT)),
        "external_smarts_total_count": int(len(catalog)),
        "external_smarts_valid_count": int(len(patterns)),
        "external_smarts_invalid_count": int(invalid_count),
    }


def external_alert_names(mol: Chem.Mol | None, patterns: list[tuple[str, Chem.Mol]]) -> list[str]:
    """Return external unwanted-substructure alert names for one molecule."""
    if mol is None:
        return []
    return sorted(name for name, pattern in patterns if mol.HasSubstructMatch(pattern))


def join_names(names: list[str]) -> str:
    """Join alert names for CSV storage."""
    return ";".join(names)


def family_overlap_metrics(df: pd.DataFrame) -> dict[str, int]:
    """Summarize overlap between PAINS, Brenk, and external unwanted alerts."""
    pains = df["pains_flag"]
    brenk = df["brenk_flag"]
    unwanted = df["unwanted_substructure_flag"]
    return {
        "pains_only_count": int((pains & ~brenk & ~unwanted).sum()),
        "brenk_only_count": int((~pains & brenk & ~unwanted).sum()),
        "unwanted_csv_only_count": int((~pains & ~brenk & unwanted).sum()),
        "pains_brenk_count": int((pains & brenk & ~unwanted).sum()),
        "brenk_unwanted_csv_count": int((~pains & brenk & unwanted).sum()),
        "pains_unwanted_csv_count": int((pains & ~brenk & unwanted).sum()),
        "all_alert_families_count": int((pains & brenk & unwanted).sum()),
        "no_alert_count": int((~pains & ~brenk & ~unwanted).sum()),
    }


def fraction(count: int, total: int) -> float:
    """Return count fraction with zero-safe behavior."""
    return float(count / total) if total else 0.0


def main() -> None:
    """Annotate molecules with PAINS, Brenk, NIH, and optional SMARTS alerts."""
    df, input_path = load_input_table()
    id_column = detect_id_column(df)
    smiles_column = detect_smiles_column(df)

    pains_catalog = build_filter_catalog(["PAINS_A", "PAINS_B", "PAINS_C"])
    brenk_catalog = build_filter_catalog(["BRENK"])
    nih_catalog = build_filter_catalog(["NIH"])
    external_path = find_external_unwanted_csv()
    external_patterns, external_metrics = load_external_patterns(external_path)

    molecules = [Chem.MolFromSmiles(str(value)) for value in df[smiles_column]]
    invalid_molecule_count = sum(mol is None for mol in molecules)

    pains_names = [catalog_alert_names(mol, pains_catalog) for mol in molecules]
    brenk_names = [catalog_alert_names(mol, brenk_catalog) for mol in molecules]
    nih_names = [catalog_alert_names(mol, nih_catalog) for mol in molecules]
    unwanted_names = [external_alert_names(mol, external_patterns) for mol in molecules]

    annotated = df.copy()
    if id_column != "molecule_chembl_id" and "molecule_chembl_id" not in annotated.columns:
        annotated["molecule_id"] = annotated[id_column]

    annotated["pains_alert_count"] = [len(names) for names in pains_names]
    annotated["pains_flag"] = annotated["pains_alert_count"] > 0
    annotated["pains_alert_names"] = [join_names(names) for names in pains_names]
    annotated["brenk_alert_count"] = [len(names) for names in brenk_names]
    annotated["brenk_flag"] = annotated["brenk_alert_count"] > 0
    annotated["brenk_alert_names"] = [join_names(names) for names in brenk_names]
    annotated["nih_alert_count"] = [len(names) for names in nih_names]
    annotated["nih_alert_flag"] = annotated["nih_alert_count"] > 0
    annotated["unwanted_substructure_count"] = [len(names) for names in unwanted_names]
    annotated["unwanted_substructure_flag"] = annotated["unwanted_substructure_count"] > 0
    annotated["unwanted_substructure_names"] = [join_names(names) for names in unwanted_names]
    annotated["medchem_alert_count"] = (
        annotated["pains_alert_count"]
        + annotated["brenk_alert_count"]
        + annotated["nih_alert_count"]
        + annotated["unwanted_substructure_count"]
    )
    annotated["medchem_alert_flag"] = annotated["medchem_alert_count"] > 0

    def summarize_alerts(row: pd.Series) -> str:
        parts: list[str] = []
        for family, column in [
            ("PAINS", "pains_alert_count"),
            ("Brenk", "brenk_alert_count"),
            ("NIH", "nih_alert_count"),
            ("unwanted", "unwanted_substructure_count"),
        ]:
            count = int(row[column])
            if count:
                parts.append(f"{family}:{count}")
        return ";".join(parts) if parts else "none"

    annotated["medchem_alert_summary"] = annotated.apply(summarize_alerts, axis=1)
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(OUTPUT_PATH, index=False)

    total = len(annotated)
    counts = {
        "input_row_count": int(len(df)),
        "output_row_count": int(total),
        "invalid_molecule_count": int(invalid_molecule_count),
        "pains_flagged_count": int(annotated["pains_flag"].sum()),
        "brenk_flagged_count": int(annotated["brenk_flag"].sum()),
        "nih_flagged_count": int(annotated["nih_alert_flag"].sum()),
        "unwanted_substructure_flagged_count": int(annotated["unwanted_substructure_flag"].sum()),
        "combined_medchem_alert_count": int(annotated["medchem_alert_flag"].sum()),
    }
    metrics = {
        **counts,
        "input_path": str(input_path.relative_to(PROJECT_ROOT)),
        "output_path": str(OUTPUT_PATH.relative_to(PROJECT_ROOT)),
        "id_column": id_column,
        "smiles_column": smiles_column,
        "pains_flagged_fraction": fraction(counts["pains_flagged_count"], total),
        "brenk_flagged_fraction": fraction(counts["brenk_flagged_count"], total),
        "nih_flagged_fraction": fraction(counts["nih_flagged_count"], total),
        "unwanted_substructure_flagged_fraction": fraction(counts["unwanted_substructure_flagged_count"], total),
        "combined_medchem_alert_fraction": fraction(counts["combined_medchem_alert_count"], total),
        "rdkit_filter_catalogs": {
            "PAINS": pains_catalog is not None,
            "BRENK": brenk_catalog is not None,
            "NIH": nih_catalog is not None,
        },
        **external_metrics,
        **family_overlap_metrics(annotated),
        "interpretation": "PAINS, Brenk, NIH, and external SMARTS alerts are medicinal-chemistry risk annotations and sensitivity-analysis filters, not automatic exclusions from the primary EGFR QSAR benchmark.",
    }
    save_json(METRICS_PATH, metrics)

    report = [
        "# EGFR Medicinal-Chemistry Alert Report",
        "",
        "PAINS, Brenk, NIH, and external unwanted-substructure SMARTS alerts were used as medicinal-chemistry risk annotations and sensitivity-analysis filters, not automatic exclusions from the primary EGFR QSAR benchmark.",
        "",
        "## Input And Output",
        "",
        f"- Input table: `{metrics['input_path']}`",
        f"- Output table: `{metrics['output_path']}`",
        f"- Input rows: {counts['input_row_count']:,}",
        f"- Output rows: {counts['output_row_count']:,}",
        f"- Invalid molecules during alert parsing: {counts['invalid_molecule_count']:,}",
        "",
        "## Alert Counts",
        "",
        f"- PAINS-flagged molecules: {counts['pains_flagged_count']:,} ({metrics['pains_flagged_fraction']:.1%})",
        f"- Brenk-flagged molecules: {counts['brenk_flagged_count']:,} ({metrics['brenk_flagged_fraction']:.1%})",
        f"- NIH-flagged molecules: {counts['nih_flagged_count']:,} ({metrics['nih_flagged_fraction']:.1%})",
        f"- External unwanted-substructure flagged molecules: {counts['unwanted_substructure_flagged_count']:,} ({metrics['unwanted_substructure_flagged_fraction']:.1%})",
        f"- Any medicinal-chemistry alert: {counts['combined_medchem_alert_count']:,} ({metrics['combined_medchem_alert_fraction']:.1%})",
        "",
        "## External SMARTS Catalog",
        "",
        f"- External CSV found: {metrics['external_csv_found']}",
        f"- External CSV path: {metrics.get('external_csv_path') or 'not found'}",
        f"- SMARTS total: {metrics['external_smarts_total_count']:,}",
        f"- SMARTS valid: {metrics['external_smarts_valid_count']:,}",
        f"- SMARTS invalid: {metrics['external_smarts_invalid_count']:,}",
        "",
        "## Family Overlap",
        "",
        f"- PAINS only: {metrics['pains_only_count']:,}",
        f"- Brenk only: {metrics['brenk_only_count']:,}",
        f"- External unwanted CSV only: {metrics['unwanted_csv_only_count']:,}",
        f"- PAINS + Brenk: {metrics['pains_brenk_count']:,}",
        f"- Brenk + external unwanted CSV: {metrics['brenk_unwanted_csv_count']:,}",
        f"- PAINS + external unwanted CSV: {metrics['pains_unwanted_csv_count']:,}",
        f"- All alert families: {metrics['all_alert_families_count']:,}",
        f"- No PAINS/Brenk/external unwanted alert: {metrics['no_alert_count']:,}",
        "",
        "## Interpretation",
        "",
        "The alert layer is used for risk labeling, triage penalties, and sensitivity checks. It does not prove that a molecule is inactive, toxic, or an assay artifact.",
        "",
    ]
    write_text(REPORT_PATH, "\n".join(report))

    print(f"Annotated table: {OUTPUT_PATH}")
    print(f"Input rows: {counts['input_row_count']}")
    print(f"Output rows: {counts['output_row_count']}")
    print(f"PAINS flagged: {counts['pains_flagged_count']} ({metrics['pains_flagged_fraction']:.4f})")
    print(f"Brenk flagged: {counts['brenk_flagged_count']} ({metrics['brenk_flagged_fraction']:.4f})")
    print(
        "Unwanted-substructure flagged: "
        f"{counts['unwanted_substructure_flagged_count']} ({metrics['unwanted_substructure_flagged_fraction']:.4f})"
    )
    print(f"Combined medchem alerts: {counts['combined_medchem_alert_count']} ({metrics['combined_medchem_alert_fraction']:.4f})")
    print(f"External CSV found: {metrics['external_csv_found']}")
    print(f"Metrics: {METRICS_PATH}")


if __name__ == "__main__":
    main()
