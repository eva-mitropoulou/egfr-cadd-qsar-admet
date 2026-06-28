import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PUBLIC_PATHS = list((ROOT / "reports").glob("*.md")) + list((ROOT / "portfolio_assets").glob("*.md")) + [ROOT / "README.md"]
SMILES_LIKE = re.compile(r"(?<![A-Za-z0-9_])[BCNOFPSIclbr@+\-\[\]\(\)=#$\\/0-9]{30,}(?![A-Za-z0-9_])")


def test_public_markdown_does_not_expose_long_raw_smiles():
    offenders = []
    for path in PUBLIC_PATHS:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        if SMILES_LIKE.search(text):
            offenders.append(str(path.relative_to(ROOT)))
    assert not offenders
