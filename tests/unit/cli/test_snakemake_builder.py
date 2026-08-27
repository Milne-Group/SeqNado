"""Unit tests for Snakemake command builder."""

import os
from pathlib import Path
from unittest.mock import patch

import pytest

from seqnado.cli.snakemake_builder import SnakemakeCommandBuilder


class TestSnakemakeCommandBuilder:
    """Tests for SnakemakeCommandBuilder class."""

    def test_builder_basic_command(self, tmp_path):
        """Test building a basic Snakemake command."""
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        builder = SnakemakeCommandBuilder(snakefile, cores=4)
        cmd = builder.build()

        assert "snakemake" in cmd
        assert "--snakefile" in cmd
        assert str(snakefile) in cmd
        assert "-c" in cmd
        assert "4" in cmd

    def test_builder_add_configfile(self, tmp_path):
        """Test adding configfile to command."""
        snakefile = tmp_path / "Snakefile"
        config = tmp_path / "config.yaml"
        snakefile.touch()
        config.touch()

        builder = SnakemakeCommandBuilder(snakefile)
        cmd = builder.add_configfile(config).build()

        assert "--configfile" in cmd
        assert str(config) in cmd

    def test_builder_add_config(self, tmp_path):
        """Test adding config key-value pairs."""
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        builder = SnakemakeCommandBuilder(snakefile)
        cmd = builder.add_config(genome="hg38", output_dir="/tmp").build()

        assert "--config" in cmd
        assert "genome=hg38" in cmd
        assert "output_dir=/tmp" in cmd

    def test_builder_add_dry_run(self, tmp_path):
        """Test adding --dry-run flag."""
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        builder = SnakemakeCommandBuilder(snakefile)
        cmd = builder.add_dry_run().build()

        assert "--dry-run" in cmd

    def test_builder_add_unlock(self, tmp_path):
        """Test adding --unlock flag."""
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        builder = SnakemakeCommandBuilder(snakefile)
        cmd = builder.add_unlock().build()

        assert "--unlock" in cmd

    def test_builder_chaining(self, tmp_path):
        """Test fluent API chaining."""
        snakefile = tmp_path / "Snakefile"
        config = tmp_path / "config.yaml"
        snakefile.touch()
        config.touch()

        cmd = (
            SnakemakeCommandBuilder(snakefile, cores=8)
            .add_configfile(config)
            .add_dry_run()
            .add_config(test="value")
            .build()
        )

        assert len(cmd) > 0
        assert "snakemake" in cmd
        assert "--dry-run" in cmd

    def test_builder_add_use_conda(self, tmp_path):
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        cmd = SnakemakeCommandBuilder(snakefile).add_use_conda().build()
        assert "--use-conda" in cmd

    def test_builder_add_use_envmodules(self, tmp_path):
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        cmd = SnakemakeCommandBuilder(snakefile).add_use_envmodules().build()
        assert "--use-envmodules" in cmd

    def test_builder_add_forceall(self, tmp_path):
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        cmd = SnakemakeCommandBuilder(snakefile).add_forceall().build()
        assert "--forceall" in cmd

    def test_builder_add_touch(self, tmp_path):
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        cmd = SnakemakeCommandBuilder(snakefile).add_touch().build()
        assert "--touch" in cmd

    def test_builder_add_notemp(self, tmp_path):
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        cmd = SnakemakeCommandBuilder(snakefile).add_notemp().build()
        assert "--notemp" in cmd

    def test_builder_add_report(self, tmp_path):
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()
        report_path = tmp_path / "report.html"

        cmd = SnakemakeCommandBuilder(snakefile).add_report(report_path).build()
        assert "--report" in cmd
        assert str(report_path) in cmd
        assert "--report-stylesheet" not in cmd

    def test_builder_add_queue_slurm_preset(self, tmp_path):
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        cmd = SnakemakeCommandBuilder(snakefile).add_queue("short", "ss").build()
        assert "--default-resources" in cmd
        assert "slurm_partition=short" in cmd

    def test_builder_add_queue_non_slurm_preset_warns(self, tmp_path, caplog):
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        with patch("seqnado.cli.snakemake_builder.logger") as mock_logger:
            cmd = SnakemakeCommandBuilder(snakefile).add_queue("short", "le").build()

        assert "--default-resources" not in cmd
        assert "slurm_partition=short" not in cmd
        mock_logger.warning.assert_called_once()

    def test_builder_add_queue_no_queue_noop(self, tmp_path):
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        cmd = SnakemakeCommandBuilder(snakefile).add_queue("", "ss").build()
        assert "--default-resources" not in cmd

    def test_builder_add_report_with_stylesheet(self, tmp_path):
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()
        report_path = tmp_path / "report.html"
        stylesheet = tmp_path / "style.css"

        cmd = (
            SnakemakeCommandBuilder(snakefile)
            .add_report(report_path, stylesheet)
            .build()
        )
        assert "--report" in cmd
        assert str(report_path) in cmd
        assert "--report-stylesheet" in cmd
        assert str(stylesheet) in cmd


class TestSnakemakeCommandBuilderRun:
    """Tests for SnakemakeCommandBuilder.run(), mocking Snakemake's own entrypoint."""

    def test_run_success_returns_zero(self, tmp_path):
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        with (
            patch("seqnado.cli.snakemake_builder.parse_args") as mock_parse_args,
            patch("seqnado.cli.snakemake_builder.args_to_api") as mock_args_to_api,
        ):
            mock_parse_args.return_value = ("parser", "args")
            mock_args_to_api.return_value = True

            exit_code = SnakemakeCommandBuilder(snakefile).run(str(tmp_path))

        assert exit_code == 0

    def test_run_failure_returns_one(self, tmp_path):
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        with (
            patch("seqnado.cli.snakemake_builder.parse_args") as mock_parse_args,
            patch("seqnado.cli.snakemake_builder.args_to_api") as mock_args_to_api,
        ):
            mock_parse_args.return_value = ("parser", "args")
            mock_args_to_api.return_value = False

            exit_code = SnakemakeCommandBuilder(snakefile).run(str(tmp_path))

        assert exit_code == 1

    def test_run_resolves_cwd_and_sets_pwd(self, tmp_path):
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()
        original_cwd = os.getcwd()

        try:
            with (
                patch("seqnado.cli.snakemake_builder.parse_args") as mock_parse_args,
                patch("seqnado.cli.snakemake_builder.args_to_api") as mock_args_to_api,
            ):
                mock_parse_args.return_value = ("parser", "args")
                mock_args_to_api.return_value = True

                SnakemakeCommandBuilder(snakefile).run(str(tmp_path))

            resolved = str(tmp_path.resolve())
            assert os.getcwd() == resolved
            assert os.environ["PWD"] == resolved
        finally:
            os.chdir(original_cwd)

    def test_run_injects_directory_flag_when_absent(self, tmp_path):
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        with (
            patch("seqnado.cli.snakemake_builder.parse_args") as mock_parse_args,
            patch("seqnado.cli.snakemake_builder.args_to_api") as mock_args_to_api,
        ):
            mock_parse_args.return_value = ("parser", "args")
            mock_args_to_api.return_value = True

            SnakemakeCommandBuilder(snakefile).run(str(tmp_path))

            argv = mock_parse_args.call_args[0][0]
            resolved = str(tmp_path.resolve())
            assert argv.count("--directory") == 1
            idx = argv.index("--directory")
            assert argv[idx + 1] == resolved

    def test_run_does_not_duplicate_directory_flag(self, tmp_path):
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        with (
            patch("seqnado.cli.snakemake_builder.parse_args") as mock_parse_args,
            patch("seqnado.cli.snakemake_builder.args_to_api") as mock_args_to_api,
        ):
            mock_parse_args.return_value = ("parser", "args")
            mock_args_to_api.return_value = True

            builder = SnakemakeCommandBuilder(snakefile).add_directory(".")
            builder.run(str(tmp_path))

            argv = mock_parse_args.call_args[0][0]
            assert argv.count("--directory") == 1

    def test_run_catches_systemexit_from_parse_args(self, tmp_path):
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        with patch("seqnado.cli.snakemake_builder.parse_args") as mock_parse_args:
            mock_parse_args.side_effect = SystemExit(2)

            exit_code = SnakemakeCommandBuilder(snakefile).run(str(tmp_path))

        assert exit_code == 2

    def test_run_catches_systemexit_with_non_int_code(self, tmp_path):
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        with patch("seqnado.cli.snakemake_builder.parse_args") as mock_parse_args:
            mock_parse_args.side_effect = SystemExit("some error message")

            exit_code = SnakemakeCommandBuilder(snakefile).run(str(tmp_path))

        assert exit_code == 1


@pytest.mark.snakemake
class TestSnakemakeCommandBuilderRunIntegration:
    """Real (non-mocked) parse_args/args_to_api round-trips through bundled profile YAMLs."""

    FIXTURE_SNAKEFILE = """
rule all:
    input: "done.txt"

rule touch_done:
    output: "done.txt"
    shell: "touch {output}"
"""

    def _bundled_profile_dir(self, dirname: str) -> Path:
        import seqnado

        profile_dir = (
            Path(seqnado.__file__).parent
            / "workflow"
            / "envs"
            / "profiles"
            / dirname
        )
        assert profile_dir.exists(), f"Expected bundled profile at {profile_dir}"
        return profile_dir

    def test_dry_run_with_test_profile(self, tmp_path):
        snakefile = tmp_path / "Snakefile"
        snakefile.write_text(self.FIXTURE_SNAKEFILE)

        profile_dir = self._bundled_profile_dir("profile_test")

        exit_code = (
            SnakemakeCommandBuilder(snakefile)
            .add_profile_from_path(profile_dir)
            .add_dry_run()
            .run(str(tmp_path))
        )

        assert exit_code == 0
        assert not (tmp_path / "done.txt").exists()

    def test_dry_run_with_slurm_yte_profile(self, tmp_path):
        pytest.importorskip(
            "snakemake_executor_plugin_slurm",
            reason="slurm executor plugin (seqnado[slurm] extra) not installed",
        )
        snakefile = tmp_path / "Snakefile"
        snakefile.write_text(self.FIXTURE_SNAKEFILE)

        profile_dir = self._bundled_profile_dir("profile_slurm_singularity")

        exit_code = (
            SnakemakeCommandBuilder(snakefile)
            .add_profile_from_path(profile_dir)
            .add_dry_run()
            .run(str(tmp_path))
        )

        assert exit_code == 0
        assert not (tmp_path / "done.txt").exists()

