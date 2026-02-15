"""Centralized Snakemake command builder for consistent CLI command construction."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import List, Optional

from loguru import logger

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

    def _add_cores(self) -> None:
        """Add cores option to command."""
        self.cmd.extend(["-c", str(self.cores)])

    def add_configfile(self, config_file: Path | str) -> SnakemakeCommandBuilder:
        """Add configfile to command. Returns self for chaining."""
        self.cmd.extend(["--configfile", str(config_file)])
        return self

    def add_config(self, **kwargs) -> SnakemakeCommandBuilder:
        """
        Add key=value config pairs to command. Returns self for chaining.
        
        Example:
            builder.add_config(genome="hg38", output_dir="/tmp")
        """
        config_parts = [f"{k}={v}" for k, v in kwargs.items()]
        self.cmd.extend(["--config"] + config_parts)
        return self

    def add_profile(self, preset: str, pkg_root_trav=None) -> SnakemakeCommandBuilder:
        """
        Add Snakemake profile to command. Returns self for chaining.
        
        Args:
            preset: Profile preset name (e.g., "le", "ss")
            pkg_root_trav: Optional importlib.resources Traversable for package root
        
        Returns:
            self for method chaining
        """
        if not preset:
            return self

        # If pkg_root_trav provided, use resolve_profile_path to locate profile
        if pkg_root_trav:
            from importlib import resources

            profile_trav = resolve_profile_path(preset, pkg_root_trav)
            if profile_trav:
                try:
                    with resources.as_file(profile_trav) as profile_path:
                        if Path(profile_path).exists():
                            self.cmd.extend(["--profile", str(profile_path)])
                            logger.info(f"Using Snakemake profile preset '{preset}' -> {profile_path}")
                        else:
                            logger.warning(f"Profile path does not exist: {profile_path}")
                except Exception as e:
                    logger.warning(f"Could not resolve profile {preset}: {e}")
            else:
                profiles = get_preset_profiles()
                logger.warning(
                    f"Unknown preset '{preset}'. Available: {list(profiles.keys())}"
                )
        return self

    def add_unlock(self) -> SnakemakeCommandBuilder:
        """Add --unlock flag. Returns self for chaining."""
        self.cmd.append("--unlock")
        return self

    def add_dry_run(self) -> SnakemakeCommandBuilder:
        """Add --dry-run flag. Returns self for chaining."""
        self.cmd.append("--dry-run")
        return self

    def add_show_failed_logs(self) -> SnakemakeCommandBuilder:
        """Add --show-failed-logs flag. Returns self for chaining."""
        if "--show-failed-logs" not in self.cmd:
            self.cmd.append("--show-failed-logs")
        return self

    def add_container_support(self) -> SnakemakeCommandBuilder:
        """
        Add container support (apptainer or singularity). Returns self for chaining.
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
        Add pass-through arguments to Snakemake command. Returns self for chaining.
        
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
        Add a specific target/rule to execute. Returns self for chaining.
        
        Args:
            target: Target rule name (e.g., "geo_download_all")
        
        Returns:
            self for method chaining
        """
        self.cmd.append(target)
        return self

    def build(self) -> List[str]:
        """Return the complete command as a list of strings."""
        return self.cmd

    def build_str(self) -> str:
        """Return the complete command as a shell command string."""
        return " ".join(self.cmd)
