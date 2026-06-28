"""Small fallback smoke-test runner when pytest is unavailable."""

from __future__ import annotations

import importlib.util
import inspect
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TEST_DIR = ROOT / "tests"
METRICS_PATH = ROOT / "reports" / "metrics" / "egfr_smoke_test_metrics.json"
REPORT_PATH = ROOT / "reports" / "egfr_smoke_test_report.md"


def load_module(path: Path):
    """Load a test module from a path."""
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main() -> None:
    """Run test_ functions directly."""
    failures: list[str] = []
    passed = 0
    for path in sorted(TEST_DIR.glob("test_*.py")):
        module = load_module(path)
        for name, function in inspect.getmembers(module, inspect.isfunction):
            if not name.startswith("test_"):
                continue
            try:
                function()
                passed += 1
            except Exception as exc:
                failures.append(f"{path.name}:{name}:{exc.__class__.__name__}")
    payload = {"status": "passed" if not failures else "failed", "passed": passed, "failures": failures}
    METRICS_PATH.parent.mkdir(parents=True, exist_ok=True)
    METRICS_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    report = [
        "# EGFR Smoke Test Report",
        "",
        f"- Direct smoke-test status: {payload['status']}",
        f"- Direct smoke-test functions passed: {passed}",
        f"- Direct smoke-test failures: {len(failures)}",
        "- Pytest status: unavailable in this environment; direct test runner was used as the documented fallback.",
        "",
    ]
    if failures:
        report.append("## Failures")
        report.append("")
        report.extend(f"- {failure}" for failure in failures)
    REPORT_PATH.write_text("\n".join(report), encoding="utf-8")
    if failures:
        raise SystemExit(1)
    print("Smoke tests status: passed")


if __name__ == "__main__":
    main()
