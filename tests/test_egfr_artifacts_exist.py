from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_required_public_artifacts_exist():
    required = [
        "reports/egfr_assay_aware_validation_report.md",
        "reports/egfr_conformal_uncertainty_report.md",
        "reports/egfr_sar_interpretability_report.md",
        "reports/egfr_redocking_audit_report.md",
        "reports/final_egfr_cadd_qsar_report.md",
        "reports/final_egfr_cv_bullets.md",
        "docs/project_card.md",
    ]
    missing = [path for path in required if not (ROOT / path).exists()]
    assert not missing


def test_redocking_artifacts_exist_and_nonempty():
    required = [
        "data/structure_prepared/5UG9_receptor.pdbqt",
        "data/structure_prepared/5UG9_8AM_ligand.pdbqt",
        "data/structure_prepared/5UG9_8AM_redocked_out.pdbqt",
    ]
    for path in required:
        artifact = ROOT / path
        assert artifact.exists()
        assert artifact.stat().st_size > 0
