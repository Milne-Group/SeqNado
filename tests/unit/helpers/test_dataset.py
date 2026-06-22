"""Tests for QuantNado dataset workflow helpers."""

import importlib.util
from pathlib import Path
from unittest.mock import Mock

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
HELPER_PATH = REPO_ROOT / "seqnado" / "workflow" / "helpers" / "dataset.py"
RULE_PATH = REPO_ROOT / "seqnado" / "workflow" / "rules" / "dataset" / "dataset.smk"


def _load_dataset_helper():
    spec = importlib.util.spec_from_file_location("dataset_helper", HELPER_PATH)
    assert spec is not None
    assert spec.loader is not None

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


get_quantnado_pairedness_flag = _load_dataset_helper().get_quantnado_pairedness_flag


@pytest.mark.unit
def test_get_quantnado_pairedness_flag_omits_paired_default():
    input_files = Mock()
    input_files.is_paired_end.return_value = True

    assert get_quantnado_pairedness_flag(input_files, "sample1") == ""


@pytest.mark.unit
def test_get_quantnado_pairedness_flag_marks_single_end():
    input_files = Mock()
    input_files.is_paired_end.return_value = False

    assert get_quantnado_pairedness_flag(input_files, "sample1") == "--single-end"


@pytest.mark.unit
def test_get_quantnado_pairedness_flag_defaults_to_paired_when_unknown():
    input_files = object()

    assert get_quantnado_pairedness_flag(input_files, "sample1") == ""


@pytest.mark.unit
def test_get_quantnado_pairedness_flag_defaults_to_paired_on_missing_sample():
    input_files = Mock()
    input_files.is_paired_end.side_effect = KeyError("sample1")

    assert get_quantnado_pairedness_flag(input_files, "sample1") == ""


@pytest.mark.unit
def test_dataset_rule_uses_pairedness_flag_not_boolean_option():
    rule = RULE_PATH.read_text()

    assert 'paired_end="' not in rule
    assert "--paired {params.paired_end}" not in rule
    assert "{params.pairedness_flag}" in rule
