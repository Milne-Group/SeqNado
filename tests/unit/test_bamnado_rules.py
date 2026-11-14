"""Static checks on bamnado invocations in the Snakemake rules."""

import re
from pathlib import Path

import pytest

RULES_DIR = Path(__file__).resolve().parents[2] / "seqnado" / "workflow" / "rules"


@pytest.mark.unit
def test_bamnado_scale_factor_uses_equals_form():
    """bamnado's CLI parser treats a separate negative value (`--scale-factor -1.5`)
    as an unknown flag, and minus-strand spike-in factors are negative. Every
    templated `--scale-factor` must use the `--scale-factor=VALUE` form."""
    offenders = []
    for smk in RULES_DIR.rglob("*.smk"):
        for lineno, line in enumerate(smk.read_text().splitlines(), start=1):
            if re.search(r"--scale-factor\s+\{", line):
                offenders.append(f"{smk.relative_to(RULES_DIR)}:{lineno}")
    assert not offenders, "Use --scale-factor={...}: " + ", ".join(offenders)
