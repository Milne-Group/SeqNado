"""Tests for seqnado.config.third_party_tools module."""

import pytest

from seqnado.config.third_party_tools import (
    CommandLineArguments,
    ThirdPartyToolsConfig,
    ToolConfig,
    none_str_to_none,
)
from seqnado import Assay


class TestNoneStrToNone:
    """Tests for the none_str_to_none converter function."""

    def test_none_string_converts_to_none(self):
        """Test that string 'None' converts to None."""
        assert none_str_to_none("None") is None

    def test_none_string_lowercase_converts(self):
        """Test that 'none' converts to None."""
        assert none_str_to_none("none") is None

    def test_none_string_with_whitespace_converts(self):
        """Test that ' None ' converts to None."""
        assert none_str_to_none(" None ") is None

    def test_non_none_string_unchanged(self):
        """Test that non-None strings are unchanged."""
        assert none_str_to_none("value") == "value"

    def test_non_string_unchanged(self):
        """Test that non-string values are unchanged."""
        assert none_str_to_none(42) == 42
        assert none_str_to_none([1, 2]) == [1, 2]


class TestCommandLineArguments:
    """Tests for CommandLineArguments validation and filtering."""

    def test_empty_value(self):
        """Test CommandLineArguments with empty string."""
        args = CommandLineArguments(value="")
        assert args.value == ""
        assert args.option_string_filtered == ""

    def test_leading_space_added(self):
        """Test that leading space is added if missing."""
        args = CommandLineArguments(value="--opt value")
        assert args.value.startswith(" ")

    def test_leading_space_preserved(self):
        """Test that existing leading space is preserved."""
        args = CommandLineArguments(value=" --opt value")
        assert args.value == " --opt value"

    def test_invalid_syntax_raises_error(self):
        """Test that invalid CLI syntax raises ValueError."""
        with pytest.raises(ValueError, match="Invalid CLI option string"):
            CommandLineArguments(value="--opt 'unclosed quote")

    def test_unsafe_characters_rejected(self):
        """Test that newlines and null bytes are rejected."""
        with pytest.raises(ValueError, match="unsafe characters"):
            CommandLineArguments(value="--opt\nvalue")
        
        with pytest.raises(ValueError, match="unsafe characters"):
            CommandLineArguments(value="--opt\x00value")

    def test_non_string_value_rejected(self):
        """Test that non-string values are rejected."""
        with pytest.raises(ValueError, match="must be a string"):
            CommandLineArguments(value=123)

    def test_exclude_single_option(self):
        """Test excluding a single option."""
        args = CommandLineArguments(
            value="--opt1 value1 --opt2 value2",
            exclude={"--opt1"}
        )
        filtered = args.option_string_filtered
        assert "--opt1" not in filtered
        assert "value1" not in filtered
        assert "--opt2" in filtered
        assert "value2" in filtered

    def test_exclude_option_with_equals(self):
        """Test excluding option with equals syntax."""
        args = CommandLineArguments(
            value="--opt1=value1 --opt2=value2",
            exclude={"--opt1"}
        )
        filtered = args.option_string_filtered
        assert "--opt1" not in filtered
        assert "--opt2=value2" in filtered

    def test_include_ensures_pattern_present(self):
        """Test that include patterns are added if missing."""
        args = CommandLineArguments(
            value="--existing value",
            include={"--required"}
        )
        filtered = args.option_string_filtered
        assert "--required" in filtered
        assert "--existing" in filtered

    def test_include_does_not_duplicate(self):
        """Test that include patterns already present are not duplicated."""
        args = CommandLineArguments(
            value="--opt value",
            include={"--opt"}
        )
        filtered = args.option_string_filtered
        # Should only appear once
        assert filtered.count("--opt") == 1

    def test_exclude_removes_then_include_adds(self):
        """Test that exclude removes, then include can re-add pattern."""
        args = CommandLineArguments(
            value="--opt value",
            include={"--opt"},
            exclude={"--opt"}
        )
        filtered = args.option_string_filtered
        # Pattern is excluded, then re-added by include (without the value)
        assert "--opt" in filtered
        # But the original value is excluded
        tokens = filtered.split()
        # Include adds "--opt" back as a standalone pattern
        assert tokens == ["--opt"]

    def test_add_include_method(self):
        """Test add_include method."""
        args = CommandLineArguments(value="")
        args.add_include("--opt1", "--opt2")
        assert "--opt1" in args.include
        assert "--opt2" in args.include

    def test_remove_include_method(self):
        """Test remove_include method."""
        args = CommandLineArguments(value="", include={"--opt1", "--opt2"})
        args.remove_include("--opt1")
        assert "--opt1" not in args.include
        assert "--opt2" in args.include

    def test_clear_include_method(self):
        """Test clear_include method."""
        args = CommandLineArguments(value="", include={"--opt1", "--opt2"})
        args.clear_include()
        assert len(args.include) == 0

    def test_set_include_method(self):
        """Test set_include method."""
        args = CommandLineArguments(value="", include={"--opt1"})
        args.set_include({"--opt2", "--opt3"})
        assert args.include == {"--opt2", "--opt3"}

    def test_add_exclude_method(self):
        """Test add_exclude method."""
        args = CommandLineArguments(value="")
        args.add_exclude("--opt1", "--opt2")
        assert "--opt1" in args.exclude
        assert "--opt2" in args.exclude

    def test_remove_exclude_method(self):
        """Test remove_exclude method."""
        args = CommandLineArguments(value="", exclude={"--opt1", "--opt2"})
        args.remove_exclude("--opt1")
        assert "--opt1" not in args.exclude
        assert "--opt2" in args.exclude

    def test_clear_exclude_method(self):
        """Test clear_exclude method."""
        args = CommandLineArguments(value="", exclude={"--opt1", "--opt2"})
        args.clear_exclude()
        assert len(args.exclude) == 0

    def test_set_exclude_method(self):
        """Test set_exclude method."""
        args = CommandLineArguments(value="", exclude={"--opt1"})
        args.set_exclude({"--opt2", "--opt3"})
        assert args.exclude == {"--opt2", "--opt3"}

    def test_method_chaining(self):
        """Test that methods return self for chaining."""
        args = CommandLineArguments(value="")
        result = args.add_include("--opt1").add_exclude("--opt2").clear_include()
        assert result is args

    def test_option_string_raw_property(self):
        """Test option_string_raw returns unfiltered value."""
        args = CommandLineArguments(
            value="--opt1 value1 --opt2 value2",
            exclude={"--opt1"}
        )
        raw = args.option_string_raw
        assert "--opt1" in raw
        assert "--opt2" in raw

    def test_str_method_returns_filtered_stripped(self):
        """Test __str__ returns filtered and stripped string."""
        args = CommandLineArguments(
            value="--opt1 value1 --opt2 value2",
            exclude={"--opt1"}
        )
        result = str(args)
        assert "--opt1" not in result
        assert "--opt2" in result
        assert not result.startswith(" ")

    def test_complex_filtering_scenario(self):
        """Test complex scenario with multiple includes and excludes."""
        args = CommandLineArguments(
            value="--keep1 v1 --remove1 v2 --keep2=v3",
            include={"--keep1", "--keep2", "--add-new"},
            exclude={"--remove1"}
        )
        filtered = args.option_string_filtered
        assert "--keep1" in filtered
        assert "--remove1" not in filtered
        assert "--keep2" in filtered  # Kept because it's in the original value
        assert "--add-new" in filtered

    def test_model_post_init_initializes_sets(self):
        """Test that None exclude/include are initialized to empty sets."""
        args = CommandLineArguments(value="--opt")
        assert isinstance(args.exclude, set)
        assert isinstance(args.include, set)


class TestToolConfig:
    """Tests for ToolConfig base class."""

    def test_default_values(self):
        """Test ToolConfig default values."""
        config = ToolConfig()
        assert config.threads == 1
        assert isinstance(config.command_line_arguments, CommandLineArguments)

    def test_threads_must_be_positive(self):
        """Test that threads must be >= 1."""
        with pytest.raises(ValueError):
            ToolConfig(threads=0)
        
        with pytest.raises(ValueError):
            ToolConfig(threads=-1)

    def test_threads_set_correctly(self):
        """Test setting threads value."""
        config = ToolConfig(threads=8)
        assert config.threads == 8

    def test_command_line_arguments_from_string(self):
        """Test creating ToolConfig with CLI args as string."""
        config = ToolConfig(command_line_arguments="--opt value")
        assert isinstance(config.command_line_arguments, CommandLineArguments)
        assert "--opt" in config.command_line_arguments.value

    def test_command_line_arguments_from_object(self):
        """Test creating ToolConfig with CLI args as object."""
        args = CommandLineArguments(value="--opt value", exclude={"--opt"})
        config = ToolConfig(command_line_arguments=args)
        assert config.command_line_arguments is args
        assert "--opt" not in config.command_line_arguments.option_string_filtered

    def test_serialization_uses_filtered_string(self):
        """Test that serialization uses filtered option string."""
        args = CommandLineArguments(value="--opt1 v1 --opt2 v2", exclude={"--opt1"})
        config = ToolConfig(command_line_arguments=args)
        
        # Serialize to dict
        data = config.model_dump()
        
        # The serialized command_line_arguments should be the filtered string
        assert "--opt1" not in data["command_line_arguments"]
        assert "--opt2" in data["command_line_arguments"]


class TestThirdPartyToolsConfig:
    """Tests for assay-specific third-party tool defaults."""

    def test_atac_bowtie2_defaults_include_max_insert_size(self):
        """ATAC defaults should include Bowtie2 max insert size constraint."""
        config = ThirdPartyToolsConfig.for_assay(Assay.ATAC)

        assert config.bowtie2 is not None
        options = str(config.bowtie2.align.command_line_arguments)
        assert "--maxins 2000" in options

    def test_atac_bowtie2_max_insert_size_is_not_duplicated(self):
        """ATAC defaults should not duplicate max insert size when already present."""
        config = ThirdPartyToolsConfig.for_assay(
            Assay.ATAC,
            bowtie2={
                "align": {
                    "command_line_arguments": "--very-sensitive --maxins 2000"
                }
            },
        )

        assert config.bowtie2 is not None
        options = str(config.bowtie2.align.command_line_arguments)
        assert options.count("--maxins 2000") == 1

    def test_non_atac_bowtie2_defaults_do_not_force_max_insert_size(self):
        """Non-ATAC defaults should not inject Bowtie2 max insert size by default."""
        config = ThirdPartyToolsConfig.for_assay(Assay.CHIP)

        assert config.bowtie2 is not None
        options = str(config.bowtie2.align.command_line_arguments)
        assert "-X" not in options
        assert "--maxins" not in options
