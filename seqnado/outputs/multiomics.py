from pathlib import Path
from typing import Annotated, Any

from pydantic import BaseModel, BeforeValidator, Field
from seqnado import Assay


def none_str_to_none(v):
    """Convert string 'none' to None."""
    if isinstance(v, str) and v.strip().lower() == "none":
        return None
    return v


def get_assay_bigwigs(wildcards, rules, ASSAYS) -> list[str]:
    """Get all bigwigs from assay-specific 'all' rules."""
    bigwigs = []
    for assay in ASSAYS:
        rule_name = f"{assay}_all"
        inputs = getattr(rules, rule_name).input

        # Handle both single files and collections of files
        if isinstance(inputs, str):
            if inputs.endswith(".bigWig"):
                bigwigs.append(inputs)
        else:
            # InputFiles is iterable
            for file in inputs:
                if isinstance(file, str) and file.endswith(".bigWig"):
                    bigwigs.append(file)
    return bigwigs



def find_assay_config_paths(directory: Path) -> dict[Assay, Path]:
    """
    Find assay config files in the given directory.

    Note: Skips config_multiomics.yaml as it's not an assay-specific config.

    Args:
        directory (Path): Directory to search for config files.
    Returns:
        dict[Assay, Path]: Mapping of Assay to its config file path.
    
    """

    paths = dict()
    for config_file in Path(directory).glob("config_*.yaml"):
        assay_name = config_file.stem.replace("config_", "")
        assay = Assay.from_clean_name(assay_name)

        # Skip multiomics config - it's not an assay config
        if assay == Assay.MULTIOMICS:
            continue

        paths[assay] = Path(config_file)
    return paths


def find_metadata_paths(directory: Path) -> dict[Assay, Path]:
    """
    Find assay metadata files in the given directory.

    Note: Skips metadata_multiomics.csv as it's not an assay-specific metadata.

    Args:
        directory (Path): Directory to search for metadata files.
    Returns:
        dict[Assay, Path]: Mapping of Assay to its metadata file path.
    
    """

    paths = dict()
    for metadata_file in Path(directory).glob("metadata_*.csv"):
        assay_name = metadata_file.stem.replace("metadata_", "")
        assay = Assay.from_clean_name(assay_name)

        # Skip multiomics metadata - it's not an assay metadata
        if assay == Assay.MULTIOMICS:
            continue

        paths[assay] = Path(metadata_file)
    return paths


def validate_config_and_metadata(
    config_paths: dict[Assay, Path],
    metadata_paths: dict[Assay, Path],
) -> None:
    """
    Validate the existence of the config and metadata files for each assay.

    Arguments:
        config_paths (dict[Assay, Path]): Mapping of Assay to config file path.
        metadata_paths (dict[Assay, Path]): Mapping of Assay to metadata file path.
    Raises:
        FileNotFoundError: If any config or metadata file is missing.
    """
    
    assays = set([*list(config_paths.keys()), *list(metadata_paths.keys())])
    for assay in assays:
        if assay not in config_paths:
            raise FileNotFoundError(
                f"Missing config file for assay {assay.name}.\n"
                f"Please create it using 'seqnado design {assay.name}'"
            )
        if assay not in metadata_paths:
            raise FileNotFoundError(
                f"Missing metadata file for assay {assay.name}.\n"
                f"Please create it using 'seqnado config {assay.name}'"
            )


def take_browser_plotting(raw_configs: dict[Assay, dict]) -> dict | None:
    """Move per-assay plotting settings to the combined multiomics browser plot.

    In multiomics mode the genome-browser plots combine every assay's tracks in
    one figure per region, so the per-assay plots are switched off: ``plotting``
    is removed and ``plot_with_plotnado`` set False in each raw assay config
    (modified in place). Returns the first ``plotting`` section found, in
    ``raw_configs`` order, or None if no assay defines one.
    """
    browser_plotting = None
    for raw in raw_configs.values():
        assay_config = raw.get("assay_config") or {}
        plotting = assay_config.pop("plotting", None)
        if "plot_with_plotnado" in assay_config:
            assay_config["plot_with_plotnado"] = False
        if browser_plotting is None and plotting:
            browser_plotting = plotting
    return browser_plotting
