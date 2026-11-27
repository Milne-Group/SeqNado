"""Guards for Snakemake ``script:`` files in seqnado/workflow/scripts/."""

import ast
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "seqnado" / "workflow" / "scripts"
SCRIPTS = sorted(SCRIPTS_DIR.rglob("*.py"))


def _seqnado_imports(path: Path) -> list[str]:
    imports = []
    for node in ast.walk(ast.parse(path.read_text(), filename=str(path))):
        if isinstance(node, ast.Import):
            names = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            names = [node.module]
        else:
            continue
        imports += [n for n in names if n == "seqnado" or n.startswith("seqnado.")]
    return imports


@pytest.mark.unit
@pytest.mark.parametrize("path", SCRIPTS, ids=lambda p: str(p.relative_to(SCRIPTS_DIR)))
def test_workflow_scripts_do_not_import_seqnado(path: Path):
    """Scripts run inside rule containers, where seqnado is not importable for
    editable installs (only site-packages is bound, not the source tree). Pass
    anything needed from seqnado through the rule's ``params`` instead."""
    imports = _seqnado_imports(path)
    if path.name == "create_hub_with_tracknado.py":
        # This rule explicitly runs on the host, where the installed package is available.
        rule = SCRIPTS_DIR.parent / "rules" / "visualise" / "hub.smk"
        assert "containerized: False" in rule.read_text()
        return
    assert not imports, (
        f"{path.name} imports {imports}; pass the values via rule params instead"
    )
