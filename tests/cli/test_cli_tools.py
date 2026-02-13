import re
import subprocess

ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")


def _strip_ansi(text: str) -> str:
    return ANSI_RE.sub("", text)


def test_seqnado_tools_categories_order():
    res = subprocess.run(
        ["seqnado", "tools", "--list"],
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stderr
    output = _strip_ansi(res.stdout)

    expected = [
        "1. Download",
        "2. Quality Control",
        "3. Preprocessing",
        "4. Alignment",
        "5. Analysis",
        "6. Visualization",
        "7. Reporting",
        "8. Quantification",
        "9. Utilities",
    ]

    positions = []
    for item in expected:
        assert item in output
        positions.append(output.index(item))

    assert positions == sorted(positions)


def test_seqnado_tools_category_number_filters():
    res = subprocess.run(
        ["seqnado", "tools", "--category", "1"],
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stderr
    output = _strip_ansi(res.stdout)
    assert "Download" in output


def test_seqnado_tools_category_prompt_requires_tty():
    res = subprocess.run(
        ["seqnado", "tools", "--category"],
        capture_output=True,
        text=True,
    )
    assert res.returncode != 0
    combined = _strip_ansi(res.stdout + res.stderr)
    assert "requires a TTY" in combined
