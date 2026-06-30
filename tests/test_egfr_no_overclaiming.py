from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PUBLIC_PATHS = list((ROOT / "reports").glob("*.md")) + list((ROOT / "portfolio_assets").glob("*.md")) + [ROOT / "README.md"]
BANNED_CLAIM_PHRASES = [
    "clinical candidate",
    "prospective discovery",
    "binding free energy",
]
ALLOWED_NEGATED_LIMITATIONS = [
    "Docking of top-ranked molecules was used as a structure-aware sanity check, not as proof of binding affinity, therapeutic efficacy, or prospective discovery.",
]


def normalized_public_text(text: str) -> str:
    lowered = text.lower().replace("-", " ")
    for allowed in ALLOWED_NEGATED_LIMITATIONS:
        lowered = lowered.replace(allowed.lower().replace("-", " "), "")
    return lowered


def test_public_reports_do_not_make_banned_claims():
    offenders = []
    for path in PUBLIC_PATHS:
        if not path.exists():
            continue
        for number, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), start=1):
            lowered = normalized_public_text(line)
            for phrase in BANNED_CLAIM_PHRASES:
                if phrase in lowered:
                    offenders.append(f"{path.relative_to(ROOT)}:{number}:{phrase}")
    assert not offenders
