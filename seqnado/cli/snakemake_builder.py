"""Centralized Snakemake command builder for consistent CLI command construction."""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import List, Optional

from loguru import logger
from snakemake.cli import args_to_api, parse_args

from seqnado.utils import get_preset_profiles, resolve_profile_path


class SnakemakeCommandBuilder:
    """
    Builder class for constructing Snakemake commands with common options.
    
    Provides a fluent API for composing Snakemake commands with:
    - Base options (snakefile, cores, profiles)
    - Configuration files
    - Container support
    - Pass-through arguments
    - Resource scaling
    """

    def __init__(self, snakefile_path: Path, cores: int = 1):
        """
        Initialize the builder with required options.
        
        Args:
            snakefile_path: Path to the Snakefile
            cores: Number of parallel cores (default: 1)
        """
        self.cmd: List[str] = [
            "snakemake",
            "--snakefile",
            str(snakefile_path),
            "--show-failed-logs",
        ]
        self.cores = cores
        self._add_cores()
        self.targets = []

    def _add_cores(self) -> None:
        """Add cores option to command."""
        self.cmd.extend(["-c", str(self.cores)])

    def add_configfile(self, config_file: Path | str) -> SnakemakeCommandBuilder:
        """Add configfile to command."""
        self.cmd.extend(["--configfile", str(config_file)])
        return self

    def add_config(self, **kwargs) -> SnakemakeCommandBuilder:
        """
        Add key=value config pairs to command.
        
        Example:
            builder.add_config(genome="hg38", output_dir="/tmp")
        """
        config_parts = [f"{k}={v}" for k, v in kwargs.items()]
        self.cmd.extend(["--config"] + config_parts)
        return self

    def add_profile_from_path(self, profile_path: Path | str) -> SnakemakeCommandBuilder:
        """
        Add a pre-resolved Snakemake profile path to command.
        
        Args:
            profile_path: Resolved path to profile directory
        
        Returns:
            self for method chaining
        """
        if profile_path:
            self.cmd.extend(["--profile", str(profile_path)])
        return self

    def add_profile_with_logging(
        self,
        profile_path: Optional[Path],
        preset: str,
        is_custom: bool,
    ) -> SnakemakeCommandBuilder:
        """
        Add profile path to command and log appropriate message.
        
        Args:
            profile_path: Resolved profile path (may be None)
            preset: Profile preset name (for logging)
            is_custom: Whether this is a custom profile path
        
        Returns:
            self for method chaining
        """
        if profile_path:
            self.add_profile_from_path(profile_path)
            if is_custom:
                logger.info(f"Using custom Snakemake profile: {profile_path}")
            else:
                logger.info(
                    f"Using Snakemake profile preset '{preset}' -> {profile_path}"
                )
        return self

    def add_unlock(self) -> SnakemakeCommandBuilder:
        """Add --unlock flag."""
        self.cmd.append("--unlock")
        return self

    def add_dry_run(self) -> SnakemakeCommandBuilder:
        """Add --dry-run flag."""
        self.cmd.append("--dry-run")
        return self

    def add_show_failed_logs(self) -> SnakemakeCommandBuilder:
        """Add --show-failed-logs flag."""
        if "--show-failed-logs" not in self.cmd:
            self.cmd.append("--show-failed-logs")
        return self

    def add_container_support(self) -> SnakemakeCommandBuilder:
        """
        Add container support (apptainer or singularity).
        """
        if shutil.which("apptainer"):
            self.cmd.append("--use-apptainer")
            logger.info("Using Apptainer for container support")
        elif shutil.which("singularity"):
            self.cmd.append("--use-singularity")
            logger.info("Using Singularity for container support")
        else:
            logger.warning(
                "No container runtime (apptainer/singularity) found. "
                "Some workflows may fail without containerization."
            )
        return self

    def add_pass_through_args(self, args: List[str], filter_fn=None) -> SnakemakeCommandBuilder:
        """
        Add pass-through arguments to Snakemake command.
        
        Args:
            args: List of arguments to pass through
            filter_fn: Optional callable to filter which args to include
        
        Returns:
            self for method chaining
        """
        if not args:
            return self
        
        if filter_fn:
            filtered = [arg for arg in args if filter_fn(arg)]
        else:
            filtered = args
        
        if filtered:
            self.cmd.extend(filtered)
        return self

    def add_target(self, target: str) -> SnakemakeCommandBuilder:
        """
        Add a specific target/rule to execute.
        
        Args:
            target: Target rule name (e.g., "geo_download_all")
        
        Returns:
            self for method chaining
        """

        # Target rules should be added at the end of the command, so we store them separately and append at the end
        self.targets.append(target)
        return self

    def add_directory(self, directory: str = ".") -> SnakemakeCommandBuilder:
        """
        Add --directory flag to run in a specific directory.
        
        Args:
            directory: Target directory (default: ".")
        
        Returns:
            self for method chaining
        """
        self.cmd.extend(["--directory", directory])
        return self

    def add_default_resources(self, **kwargs) -> SnakemakeCommandBuilder:
        """
        Add --default-resources flag with key=value pairs.
        
        Args:
            **kwargs: Resource key=value pairs (e.g., slurm_partition="short")
        
        Returns:
            self for method chaining
        """
        if kwargs:
            resource_parts = [f"{k}={v}" for k, v in kwargs.items()]
            self.cmd.extend(["--default-resources"] + resource_parts)
        return self

    def add_queue(self, queue: str, preset: str) -> SnakemakeCommandBuilder:
        """
        Add queue/partition to default resources (for Slurm presets).
        
        Args:
            queue: Queue name (e.g., "short", "long")
            preset: Profile preset name to check if it's a Slurm preset
        
        Returns:
            self for method chaining
        """
        if not queue:
            return self
        if preset and preset.startswith("s"):
            self.cmd.extend(["--default-resources", f"slurm_partition={queue}"])
        else:
            logger.warning(
                f"--queue '{queue}' ignored: only applies to Slurm profiles (preset 'ss'), got '{preset}'. "
                "Use --default-resources to set resource values for non-Slurm presets."
            )
        return self

    def add_workflow_args(self, workflow_args: List[str]) -> SnakemakeCommandBuilder:
        """
        Add workflow_args config for nested Snakemake runs (multiomics mode).
        
        Args:
            workflow_args: List of arguments to pass to nested Snakemake runs
        
        Returns:
            self for method chaining
        """
        if workflow_args:
            workflow_args_str = " ".join(workflow_args)
            self.cmd.extend(["--config", f"workflow_args={workflow_args_str}"])
        return self

    def add_use_conda(self) -> SnakemakeCommandBuilder:
        """Add --use-conda flag."""
        self.cmd.append("--use-conda")
        return self

    def add_use_envmodules(self) -> SnakemakeCommandBuilder:
        """Add --use-envmodules flag."""
        self.cmd.append("--use-envmodules")
        return self

    def add_forceall(self) -> SnakemakeCommandBuilder:
        """Add --forceall flag."""
        self.cmd.append("--forceall")
        return self

    def add_touch(self) -> SnakemakeCommandBuilder:
        """Add --touch flag."""
        self.cmd.append("--touch")
        return self

    def add_notemp(self) -> SnakemakeCommandBuilder:
        """Add --notemp flag."""
        self.cmd.append("--notemp")
        return self

    def add_report(
        self, path: Path | str, stylesheet: Path | str | None = None
    ) -> SnakemakeCommandBuilder:
        """Add --report flag (and optional --report-stylesheet)."""
        self.cmd.extend(["--report", str(path)])
        if stylesheet:
            self.cmd.extend(["--report-stylesheet", str(stylesheet)])
        return self

    def build(self) -> List[str]:
        """Return the complete command as a list of strings."""
        return self.cmd + self.targets

    def build_str(self) -> str:
        """Return the complete command as a shell command string."""
        return " ".join(self.cmd + self.targets)

    def run(self, cwd: str, print_cmd: bool = False) -> int:
        """
        Resolve cwd, then invoke Snakemake in-process via its own CLI entrypoint
        (parse_args + args_to_api — the same functions snakemake.cli.main() calls).

        This reuses Snakemake's own profile-YAML/YTE parsing and full flag
        validation instead of reimplementing any of it.
        """
        resolved_cwd = str(Path(cwd).resolve())
        os.chdir(resolved_cwd)
        os.environ["PWD"] = resolved_cwd

        argv = self.build()[1:]  # drop the leading "snakemake" sentinel
        if "--directory" not in argv:
            argv = argv + ["--directory", resolved_cwd]

        if print_cmd:
            logger.info("Snakemake args:\n$ snakemake " + " ".join(map(str, argv)))

        try:
            parser, args = parse_args(argv)
            success = args_to_api(args, parser)
        except SystemExit as e:
            return e.code if isinstance(e.code, int) else 1
        return 0 if success else 1
