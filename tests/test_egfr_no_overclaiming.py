from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PUBLIC_PATHS = list((ROOT / "reports").glob("*.md")) + list((ROOT / "portfolio_assets").glob("*.md")) + [ROOT / "README.md"]
BANNED_CLAIM_PHRASES = [
    "clinical candidate",
    "prospective discovery",
    "binding free energy",
]


def test_public_reports_do_not_make_banned_claims():
    offenders = []
    for path in PUBLIC_PATHS:
        if not path.exists():
            continue
        for number, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1):
            lowered = line.lower().replace("-", " ")
            for phrase in BANNED_CLAIM_PHRASES:
                if phrase in lowered:
                    offenders.append(f"{path.relative_to(ROOT)}:{number}:{phrase}")
    assert not offenders
