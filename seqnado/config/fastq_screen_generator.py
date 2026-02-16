"""Helper functions for generating FastqScreen configuration from genome configs."""

import json
import logging
from pathlib import Path
from typing import Dict

from seqnado.config.configs import GenomeConfig

logger = logging.getLogger(__name__)


def generate_fastq_screen_config(
    genome_configs: Dict[str, GenomeConfig],
    output_path: Path,
    threads: int = 8,
    include_contaminants: bool = True,
    contaminant_base_path: Path | str | None = None,
) -> None:
    """
    Generate a fastq_screen.conf file from genome configurations.

    Args:
        genome_configs: Dictionary of genome name to GenomeConfig
        output_path: Path where the config file should be written
        threads: Number of threads for Bowtie2 (default: 8)
        include_contaminants: Whether to include common contaminant databases
        contaminant_base_path: Base path for contaminant reference files (required if include_contaminants=True)
    """

    # Convert to Path if string
    if contaminant_base_path is not None:
        contaminant_base_path = Path(contaminant_base_path)

    lines = []

    # Header
    lines.append("# FastqScreen configuration file")
    lines.append("# Auto-generated from SeqNado genome configurations")
    lines.append("")
    lines.append("#" + "-" * 31)
    lines.append("# Bowtie2")
    lines.append("#" + "-" * 31)
    lines.append("# Bowtie2 should be in your PATH")
    lines.append("# BOWTIE2 bowtie2")
    lines.append("")
    lines.append("#" + "-" * 31)
    lines.append("# Threads")
    lines.append("#" + "-" * 31)
    lines.append(f"THREADS\t\t{threads}")
    lines.append("")
    lines.append("#" + "-" * 31)
    lines.append("# Databases")
    lines.append("#" + "-" * 31)
    lines.append("")

    # Add genome databases from config
    organism_map = {
        "hg": "Human",
        "mm": "Mouse",
        "dm": "Drosophila",
        "rn": "Rat",
        "ce": "Worm",
        "sc": "Yeast",
        "at": "Arabidopsis",
    }

    # Group genomes by organism prefix to get latest version only
    genome_groups = {}
    for genome_name, genome_config in genome_configs.items():
        # Skip combined genomes (they'll be covered by individual ones)
        if "_" in genome_name and any(
            x in genome_name for x in ["dm6", "mm10", "mm39", "ecoli"]
        ):
            continue

        # Determine organism prefix
        organism_prefix = None
        for prefix in organism_map.keys():
            if genome_name.startswith(prefix):
                organism_prefix = prefix
                break

        if organism_prefix:
            # Keep the latest genome for this organism (higher version number)
            if organism_prefix not in genome_groups:
                genome_groups[organism_prefix] = (genome_name, genome_config)
            else:
                # Compare versions - keep the one with higher number/later in alphabet
                existing_name = genome_groups[organism_prefix][0]
                if genome_name > existing_name:
                    genome_groups[organism_prefix] = (genome_name, genome_config)
        else:
            # Unknown organism - keep it (might be custom genome)
            genome_groups[genome_name] = (genome_name, genome_config)

    # Now add the selected genomes
    for key, (genome_name, genome_config) in sorted(genome_groups.items()):
        # Get display name
        display_name = None
        for prefix, name in organism_map.items():
            if genome_name.startswith(prefix):
                display_name = name
                break
        if display_name is None:
            display_name = genome_name.replace("_", " ").title()

        # Get bowtie2 index path
        index_path = None
        if hasattr(genome_config.index, "prefix"):
            index_path = genome_config.index.prefix
        elif isinstance(genome_config.index, dict) and "prefix" in genome_config.index:
            index_path = genome_config.index["prefix"]

        if index_path:
            # Validate that bowtie2 index files exist
            index_path_obj = Path(index_path)
            parent = index_path_obj.parent
            prefix = index_path_obj.name

            if parent.exists():
                # Check for bowtie2 index files (.bt2 or .bt2l)
                index_files = list(parent.glob(f"{prefix}*.bt2")) + list(
                    parent.glob(f"{prefix}*.bt2l")
                )
                if index_files:
                    lines.append("#" + "-" * 31)
                    lines.append(f"# {display_name}")
                    lines.append(f"DATABASE\t{display_name}\t{index_path}")
                    lines.append("")
                else:
                    logger.warning(
                        f"No bowtie2 index files found for {genome_name} at {index_path}"
                    )
            else:
                logger.warning(
                    f"Index path does not exist for {genome_name}: {parent}"
                )

    # Add common contaminants if requested
    if include_contaminants:
        if contaminant_base_path is None:
            logger.warning(
                "include_contaminants=True but no contaminant_base_path provided. "
                "Skipping contaminant databases. Use --contaminant-path to specify location."
            )
        elif contaminant_base_path.exists():
            contaminants = [
                ("Worm", "Worm/Caenorhabditis_elegans.WBcel235"),
                ("Yeast", "Yeast/Saccharomyces_cerevisiae.R64-1-1"),
                ("Arabidopsis", "Arabidopsis/Arabidopsis_thaliana.TAIR10"),
                ("Ecoli", "E_coli/Ecoli"),
                ("Sars_cov_2", "Sars_cov_2/Sars_cov_2"),
                ("rRNA", "rRNA/GRCm38_rRNA"),
                ("MT", "Mitochondria/mitochondria"),
                ("PhiX", "PhiX/phi_plus_SNPs"),
                ("Lambda", "Lambda/Lambda"),
                ("Vectors", "Vectors/Vectors"),
                ("Adapters", "Adapters/Contaminants"),
                ("Spike_in_RNA", "Spike_in_RNA/Spike_in_RNA"),
            ]

            for name, rel_path in contaminants:
                full_path = contaminant_base_path / rel_path
                # Check if any files with this prefix exist
                parent = full_path.parent
                prefix = full_path.name
                if parent.exists():
                    # Check for bowtie2 index files
                    index_files = list(parent.glob(f"{prefix}*.bt2")) + list(
                        parent.glob(f"{prefix}*.bt2l")
                    )
                    if index_files:
                        lines.append("#" + "-" * 31)
                        lines.append(f"# {name}")
                        lines.append(f"DATABASE\t{name}\t{full_path}")
                        lines.append("")
        else:
            logger.warning(
                f"Contaminant path does not exist: {contaminant_base_path}"
            )

    # Write to file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write("\n".join(lines))


def load_genome_configs_for_fastqscreen() -> Dict[str, GenomeConfig]:
    """Load all genome configs from the user's genome_config.json."""
    import os

    from seqnado.config.configs import BowtieIndex

    config_path = (
        Path(os.getenv("SEQNADO_CONFIG", Path.home()))
        / ".config/seqnado/genome_config.json"
    )

    if not config_path.exists():
        raise FileNotFoundError(
            f"Genome config not found at {config_path}. Run 'seqnado init' first."
        )

    with open(config_path, "r") as f:
        all_configs = json.load(f)

    genome_configs = {}
    for name, config_data in all_configs.items():
        config_data["name"] = name
        # Use bowtie2 index for fastqscreen
        try:
            index = BowtieIndex(prefix=config_data.get("bt2_index"))
        except Exception as e:
            logger.warning(f"Skipping genome '{name}': {e}")
            continue
        config_data["index"] = index
        genome_configs[name] = GenomeConfig(**config_data)

    return genome_configs
