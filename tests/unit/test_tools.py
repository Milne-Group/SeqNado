"""Unit tests for tools module."""

import pytest

import seqnado.tools.tools as tools


@pytest.fixture
def sample_tools(monkeypatch):
    fake_tools = {
        "alpha": {
            "description": "Alpha tool",
            "command": "alpha",
            "category": "Download",
            "container": "docker://fake/container:latest",
        },
        "beta": {
            "description": "Beta tool",
            "command": "beta",
            "category": "Quality Control",
            "container": "docker://fake/container:latest",
        },
        "gamma": {
            "description": "Gamma tool",
            "command": "gamma",
            "category": "Custom",
            "subcommands": ["run", "check"],
            "container": "docker://fake/container:latest",
        },
    }
    monkeypatch.setattr(tools, "AVAILABLE_TOOLS", fake_tools, raising=False)
    return fake_tools


def test_get_categories_order(sample_tools):
    categories = tools.get_categories()
    assert categories[:2] == ["Download", "Quality Control"]
    assert "Custom" in categories


def test_get_tool_subcommands(sample_tools):
    assert tools.get_tool_subcommands("gamma") == ["run", "check"]
    assert tools.get_tool_subcommands("alpha") == []


def test_resolve_tool_command_with_subcommand(sample_tools):
    tool_info = tools.get_tool_info("gamma")
    assert tool_info is not None
    assert tools._resolve_tool_command(tool_info, "run") == "run"
    assert tools._resolve_tool_command(tool_info, "unknown") == "gamma"


def test_first_non_info_line():
    from seqnado.tools.version_parsing import first_non_info_line

    output = "INFO: cached\nTool v1.2.3\n"
    assert first_non_info_line(output) == "Tool v1.2.3"


def test_get_tool_version_container_filters_info(monkeypatch, sample_tools):
    def fake_run_command_in_container(*args, **kwargs):
        return 0, "INFO: cached\nBeta 2.0", ""

    monkeypatch.setattr(
        tools, "run_command_in_container", fake_run_command_in_container
    )
    version = tools.get_tool_version("beta", use_container=True)
    # Version extraction now returns just the version number, not the tool name
    assert version == "2.0"
    assert "INFO:" not in version


def test_get_tool_help_container_filters_info(monkeypatch, sample_tools):
    def fake_run_command_in_container(*args, **kwargs):
        return 0, "INFO: cached\nUsage: gamma [options]", ""

    monkeypatch.setattr(
        tools, "run_command_in_container", fake_run_command_in_container
    )
    help_text = tools.get_tool_help("gamma", use_container=True)
    assert "INFO:" not in help_text
    assert "Usage: gamma" in help_text


def test_get_tool_citation_returns_bibtex(sample_tools, monkeypatch):
    monkeypatch.setitem(sample_tools["alpha"], "citation", "samtools")
    citation = tools.get_tool_citation("alpha")
    assert citation is not None
    assert "@article{samtools," in citation


def test_get_tool_citation_returns_none_when_no_citation(sample_tools):
    citation = tools.get_tool_citation("alpha")
    assert citation is None


def test_get_tool_citation_returns_none_for_unknown_tool(sample_tools):
    citation = tools.get_tool_citation("nonexistent")
    assert citation is None


def test_run_tool_help_in_container_filters_info(monkeypatch, sample_tools):
    def fake_run_command_in_container(*args, **kwargs):
        return 0, "INFO: cached\nUsage: alpha", ""

    monkeypatch.setattr(
        tools, "run_command_in_container", fake_run_command_in_container
    )
    help_text = tools.run_tool_help_in_container("alpha")
    assert "INFO:" not in help_text
    assert "Usage: alpha" in help_text
