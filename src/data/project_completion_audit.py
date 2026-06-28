"""Audit existing EGFR project artifacts without printing raw records."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import METRICS_DIR, REPORTS_DIR, save_json, write_text  # noqa: E402


def list_files(relative_dir: str) -> list[dict[str, object]]:
    """List files with size metadata for an artifact directory."""
    path = PROJECT_ROOT / relative_dir
    if not path.exists():
        return []
    return [
        {
            "path": str(file.relative_to(PROJECT_ROOT)),
            "size_bytes": file.stat().st_size,
        }
        for file in sorted(path.rglob("*"))
        if file.is_file() and "__pycache__" not in file.parts
    ]


def csv_shape(path: Path) -> dict[str, object]:
    """Return safe CSV shape and column metadata."""
    try:
        df = pd.read_csv(path, nrows=0)
        full_rows = sum(1 for _ in path.open("r", encoding="utf-8", errors="ignore")) - 1
        return {"rows": max(full_rows, 0), "columns": list(df.columns), "column_count": len(df.columns)}
    except Exception as exc:
        return {"error": type(exc).__name__}


def main() -> None:
    """Create project completion audit artifacts."""
    artifact_groups = {
        "raw_data": list_files("data/raw"),
        "processed_data": list_files("data/processed"),
        "results": list_files("results"),
        "reports": list_files("reports"),
        "figures": list_files("figures"),
        "report_figures": list_files("reports/figures"),
        "models": list_files("models"),
        "notebooks": list_files("notebooks"),
        "scripts": list_files("scripts"),
        "src": list_files("src"),
    }

    critical_expected = [
        "data/raw/egfr_chembl_ic50_raw.csv",
        "data/processed/egfr_model_ready.csv",
        "results/fingerprint_baseline_metrics.csv",
        "results/scaffold_fingerprint_metrics.csv",
        "results/applicability_domain_summary.csv",
        "results/top_20_diverse_candidates.csv",
        "reports/project_results_summary.md",
    ]
    missing_critical = [item for item in critical_expected if not (PROJECT_ROOT / item).exists()]

    csv_metadata = {}
    for group_name in ["raw_data", "processed_data", "results"]:
        for item in artifact_groups[group_name]:
            path = PROJECT_ROOT / str(item["path"])
            if path.suffix.lower() == ".csv":
                csv_metadata[str(item["path"])] = csv_shape(path)

    metrics = {
        "artifact_counts": {group: len(items) for group, items in artifact_groups.items()},
        "artifacts": artifact_groups,
        "csv_metadata": csv_metadata,
        "missing_critical_artifacts": missing_critical,
        "reproducibility_status": "usable" if not missing_critical else "degraded",
        "unsupported_claims_added": False,
    }
    save_json(METRICS_DIR / "project_completion_audit.json", metrics)

    lines = [
        "# Project Completion Audit",
        "",
        "This audit lists available artifacts and missing critical outputs without exposing raw molecule records.",
        "",
        "## Artifact Counts",
        "",
    ]
    for group, count in metrics["artifact_counts"].items():
        lines.append(f"- {group}: {count}")
    lines.extend(["", "## Missing Critical Artifacts", ""])
    if missing_critical:
        lines.extend([f"- `{item}`" for item in missing_critical])
    else:
        lines.append("- None")
    lines.extend(["", "## CSV Shape Metadata", ""])
    for path, meta in csv_metadata.items():
        if "error" in meta:
            lines.append(f"- `{path}`: could not inspect ({meta['error']})")
        else:
            lines.append(f"- `{path}`: {meta['rows']} rows, {meta['column_count']} columns")
    lines.extend(["", "## Reproducibility Status", "", f"`{metrics['reproducibility_status']}`", ""])
    write_text(REPORTS_DIR / "project_completion_audit.md", "\n".join(lines))

    print(f"Audit report: {REPORTS_DIR / 'project_completion_audit.md'}")
    print(f"Audit metrics: {METRICS_DIR / 'project_completion_audit.json'}")


if __name__ == "__main__":
    main()
