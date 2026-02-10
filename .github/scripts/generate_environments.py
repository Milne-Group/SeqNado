#!/usr/bin/env python3
"""
Auto-generate conda environment files from pyproject.toml.

This script parses pyproject.toml and generates consistent environment files
to prevent version drift across envs/environment.yml, envs/environment_minimal.yml,
envs/testing.yml, etc.

Usage:
    python .github/scripts/generate_environments.py

Generated files:
    - envs/environment.yml (full pipeline with bioinformatics tools)
    - envs/environment_minimal.yml (minimal dev environment)
    - envs/testing.yml (testing environment)
    - containers/pipeline/environment_pipeline.yml (updated with core deps)
"""

import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple


def parse_pyproject():
    """Parse pyproject.toml and extract dependencies using regex."""
    toml_path = Path(__file__).parent.parent.parent / "pyproject.toml"
    content = toml_path.read_text()

    # Extract Python requirement
    py_match = re.search(r'requires-python\s*=\s*"([^"]+)"', content)
    py_req = py_match.group(1) if py_match else ">=3.10"

    # Extract dependencies list
    deps_match = re.search(r"dependencies\s*=\s*\[([\s\S]*?)\]", content, re.MULTILINE)

    core_deps = []
    if deps_match:
        deps_text = deps_match.group(1)
        # Extract quoted strings (dependency specifications)
        dep_lines = re.findall(r'"([^"]+)"', deps_text)
        core_deps = [d.strip() for d in dep_lines if d.strip()]

    # Extract optional dependencies (slurm)
    slurm_match = re.search(r"slurm\s*=\s*\[([\s\S]*?)\]", content, re.MULTILINE)

    slurm_deps = []
    if slurm_match:
        slurm_text = slurm_match.group(1)
        slurm_deps = re.findall(r'"([^"]+)"', slurm_text)

    return core_deps, slurm_deps, py_req


def normalize_dep(dep_str: str) -> Tuple[str, str]:
    """
    Normalize dependency string to conda format.

    Returns:
        (package_name, conda_spec)
    """
    # Find the position of the first operator
    operators = ["<=", ">=", "==", "!=", "<", ">"]
    min_pos = len(dep_str)
    
    for op in operators:
        pos = dep_str.find(op)
        if pos != -1 and pos < min_pos:
            min_pos = pos
    
    if min_pos < len(dep_str):
        name = dep_str[:min_pos].strip()
        return name, dep_str.strip()
    
    return dep_str.strip(), dep_str.strip()


def separate_conda_and_pip(deps: List[str]) -> Tuple[List[str], List[str]]:
    """
    Separate dependencies into conda and pip packages.
    Most packages are in conda, let pip be explicit.
    """
    # Packages that are typically pip-only or better via pip
    pip_only = {"tracknado", "mccnado"}

    conda_deps = []
    pip_deps = []

    for dep in deps:
        pkg_name = normalize_dep(dep)[0].replace("-", "_").lower()
        if pkg_name in pip_only:
            pip_deps.append(dep)
        else:
            conda_deps.append(dep)

    return conda_deps, pip_deps


def generate_environment_yml(core_deps: List[str], py_req: str, output_path: Path):
    """Generate envs/environment.yml with full pipeline dependencies."""

    conda_deps, pip_deps = separate_conda_and_pip(core_deps)

    env_data = {
        "name": "seqnado",
        "channels": ["conda-forge", "bioconda", "defaults"],
        "dependencies": [
            "# Core Python dependencies (synced from pyproject.toml)",
            f"python{py_req}",
        ]
        + conda_deps
        + [
            "",
            "# Bioinformatics tools",
            "bedtools",
            "bcftools",
            "bowtie2",
            "deeptools",
            "fastqc",
            "fastqsplitter",
            "homer",
            "macs2",
            "multiqc",
            "picard-slim",
            "pybedtools",
            "pybigwig",
            "pysam",
            "samtools>1.7",
            "star",
            "subread",
            "trackhub",
            "trim-galore",
            "ucsc-bedgraphtobigwig",
            "ucsc-bigbedtobed",
            "ucsc-fetchchromsizes",
        ],
    }

    if pip_deps:
        env_data["dependencies"].append({"pip": pip_deps})

    _write_yaml(output_path, env_data)
    print(f"✓ Generated {output_path}")


def generate_environment_minimal_yml(
    core_deps: List[str], py_req: str, output_path: Path
):
    """Generate envs/environment_minimal.yml for development."""

    conda_deps, pip_deps = separate_conda_and_pip(core_deps)

    # For minimal, only include essential conda packages
    essential_conda = [
        d
        for d in conda_deps
        if any(x in d.lower() for x in ["python", "snakemake", "numpy", "pandas"])
    ]

    env_data = {
        "name": "seqnado",
        "channels": ["conda-forge", "bioconda", "defaults"],
        "dependencies": [
            "# Minimal dev environment (synced from pyproject.toml)",
            f"python{py_req}",
        ]
        + essential_conda
        + [
            "pip",
        ],
    }

    if pip_deps or [d for d in conda_deps if d not in essential_conda]:
        all_pip = pip_deps + [d for d in conda_deps if d not in essential_conda]
        env_data["dependencies"].append({"pip": all_pip})

    _write_yaml(output_path, env_data)
    print(f"✓ Generated {output_path}")


def generate_testing_yml(core_deps: List[str], py_req: str, output_path: Path):
    """Generate envs/testing.yml for CI/testing."""

    conda_deps, pip_deps = separate_conda_and_pip(core_deps)

    env_data = {
        "name": "test",
        "channels": ["conda-forge", "bioconda", "defaults"],
        "dependencies": [
            "# Testing environment (synced from pyproject.toml)",
            f"python{py_req}",
        ]
        + conda_deps
        + [
            "pip",
            "pytest",
            "pytest-cov",
            "pytest-xdist",
        ],
    }

    if pip_deps:
        env_data["dependencies"].append({"pip": pip_deps})

    _write_yaml(output_path, env_data)
    print(f"✓ Generated {output_path}")


def update_pipeline_environment_yml(
    core_deps: List[str], py_req: str, output_path: Path
):
    """
    Update containers/pipeline/environment_pipeline.yml with core deps.
    Preserves existing bioinformatics and R/Bioconductor packages.
    """

    conda_deps, pip_deps = separate_conda_and_pip(core_deps)

    # Restrict Python range for containers
    py_spec = f"python{py_req},<3.13"

    env_data = {
        "name": "base",
        "channels": ["conda-forge", "bioconda", "defaults"],
        "dependencies": [
            "# Python and core dependencies (synced from pyproject.toml)",
            py_spec,
        ]
        + conda_deps
        + [
            "pip",
            "",
            "# Bioinformatics tools",
            "bcftools",
            "bedtools",
            "bowtie2",
            "cooler",
            "deeptools",
            "fastq-screen",
            "fastqc",
            "fastqsplitter",
            "flash",
            "homer",
            "macs2",
            "macs3",
            "minimap2",
            "multiqc",
            "picard",
            "pyranges",
            "qualimap",
            "salmon",
            "samtools>1.7",
            "seacr",
            "star=2.7.11b",
            "subread",
            "trackhub",
            "trim-galore",
            "ucsc-bedgraphtobigwig",
            "ucsc-bedtobigbed",
            "ucsc-bigbedtobed",
            "ucsc-fetchchromsizes",
            "methyldackel=0.6.1",
            "",
            "# Bioconductor packages",
            "bioconductor-annotationdbi",
            "bioconductor-chipseeker",
            "bioconductor-clusterprofiler",
            "bioconductor-complexheatmap",
            "bioconductor-deseq2",
            "bioconductor-edaseq",
            "bioconductor-edger",
            "bioconductor-fgsea",
            "bioconductor-genomeinfodb",
            "bioconductor-genomicalignments",
            "bioconductor-genomicfeatures",
            "bioconductor-genomicranges",
            "bioconductor-goseq",
            "bioconductor-org.hs.eg.db",
            "bioconductor-org.mm.eg.db",
            "bioconductor-pcaexplorer",
            "bioconductor-qvalue",
            "bioconductor-rgreat",
            "bioconductor-rtracklayer",
            "bioconductor-ruvseq",
            "",
            "# R packages",
            "r-base",
            "r-circlize",
            "r-dt",
            "r-ggrepel",
            "r-knitr",
            "r-msigdbr",
            "r-rcolorbrewer",
            "r-rjson",
            "r-tidyverse",
            "r-yaml",
            "",
            "# Additional utilities",
            "curl",
            "git",
            "wget",
        ],
    }

    if pip_deps:
        env_data["dependencies"].append({"pip": pip_deps})

    _write_yaml(output_path, env_data)
    print(f"✓ Generated {output_path}")


def _write_yaml(path: Path, data: dict):
    """Write YAML file with custom formatting for conda compatibility."""
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as f:
        # Custom YAML dump for better formatting
        f.write("# Environment file auto-generated from pyproject.toml\n")
        f.write(
            "# Do not edit manually - run: python .github/scripts/generate_environments.py\n\n"
        )

        if "name" in data:
            f.write(f"name: {data['name']}\n")

        if "channels" in data:
            f.write(f"channels:\n")
            for channel in data["channels"]:
                f.write(f"  - {channel}\n")

        if "dependencies" in data:
            f.write(f"dependencies:\n")
            for dep in data["dependencies"]:
                if isinstance(dep, str):
                    if dep.startswith("#"):
                        f.write(f"  {dep}\n")
                    elif dep == "":
                        f.write(f"\n")
                    else:
                        f.write(f"  - {dep}\n")
                elif isinstance(dep, dict) and "pip" in dep:
                    f.write(f"  - pip:\n")
                    for pip_pkg in dep["pip"]:
                        f.write(f"      - {pip_pkg}\n")


def main():
    """Generate all environment files."""

    print("Parsing pyproject.toml...")
    core_deps, slurm_deps, py_req = parse_pyproject()

    print(f"Found {len(core_deps)} core dependencies")
    print(f"Python requirement: {py_req}")

    repo_root = Path(__file__).parent.parent.parent

    print("\nGenerating environment files...")
    generate_environment_yml(core_deps, py_req, repo_root / "envs/environment.yml")
    generate_environment_minimal_yml(
        core_deps, py_req, repo_root / "envs/environment_minimal.yml"
    )
    generate_testing_yml(core_deps, py_req, repo_root / "envs/testing.yml")
    update_pipeline_environment_yml(
        core_deps, py_req, repo_root / "containers/pipeline/environment_pipeline.yml"
    )

    print("\n✓ All environment files generated successfully!")
    print("\nNext steps:")
    print("  1. Review the generated files")
    print(
        "  2. Commit: git add envs/*.yml containers/pipeline/environment_pipeline.yml"
    )
    print("  3. Update only pyproject.toml in the future - re-run this script")


if __name__ == "__main__":
    main()
