"""
User input module for SeqNado configuration using the new Pydantic models.
"""

import json
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional

import jinja2
from loguru import logger
from pydantic import ValidationError

if TYPE_CHECKING:
    from seqnado.config.multiomics import MultiomicsConfig

from seqnado import (
    Assay,
    MethylationMethod,
    MotifMethod,
    PeakCallingMethod,
    PileupMethod,
    QuantificationMethod,
    SNPCallingMethod,
    SpikeInMethod,
)
from seqnado.config import (
    AssaySpecificConfig,
    ATACAssayConfig,
    BigwigConfig,
    BowtieIndex,
    CATAssayConfig,
    ChIPAssayConfig,
    CRISPRAssayConfig,
    GenomeConfig,
    MCCAssayConfig,
    MCCConfig,
    MethylationAssayConfig,
    MethylationConfig,
    MLDatasetConfig,
    PeakCallingConfig,
    PlottingConfig,
    ProjectConfig,
    QCConfig,
    RNAAssayConfig,
    RNAQuantificationConfig,
    SeqnadoConfig,
    SNPAssayConfig,
    SNPCallingConfig,
    SpikeInConfig,
    STARIndex,
    UCSCHubConfig,
)


def get_user_input(
    prompt: str,
    default: Optional[str] = None,
    is_boolean: bool = False,
    choices: Optional[List[str]] = None,
    is_path: bool = False,
    required: bool = True,
    multi_select: bool = False,
) -> str | List[str]:
    """
    Prompt the user for input with validation for choices, boolean values, or path existence.
    Re-prompts until valid input is provided.

    Args:
        prompt: The prompt message to display
        default: Default value if user provides no input
        is_boolean: If True, expect yes/no input
        choices: List of valid choices
        is_path: If True, validate that path exists
        required: If False, allow empty input
        multi_select: If True, allow comma-separated multiple selections from choices

    Returns:
        str or List[str]: Single value or list of values (when multi_select=True)
    """
    while True:
        # Construct the prompt suffix based on choices or default
        if choices:
            if multi_select:
                prompt_suffix = f"({', '.join(choices)})"
            else:
                prompt_suffix = f"({'/'.join(choices)})"
        elif default is not None:
            prompt_suffix = f"(default: {default})"
        else:
            prompt_suffix = ""

        user_input = input(f"{prompt} {prompt_suffix}: ").strip()

        # Handle empty input and apply default if available
        if not user_input:
            if default is not None:
                user_input = default
            elif not required:
                return None
            else:
                print("Input cannot be empty. Please try again.")
                continue

        # Validate boolean input
        if is_boolean:
            if user_input.lower() in {"yes", "y", "true", "1"}:
                return True
            elif user_input.lower() in {"no", "n", "false", "0"}:
                return False
            else:
                print(
                    "Invalid boolean value. Please enter yes/no, y/n, true/false, or 1/0."
                )
                continue

        # Handle multi-select from choices
        if multi_select and choices:
            # Split by comma and validate each selection
            selections = [s.strip() for s in user_input.split(",")]
            invalid_selections = [s for s in selections if s not in choices]

            if invalid_selections:
                print(f"Invalid choices: {', '.join(invalid_selections)}")
                print(f"Please choose from: {', '.join(choices)}")
                continue

            if not selections:
                print("No valid selections made. Please try again.")
                continue

            return selections

        # Validate against allowed choices (single select)
        if choices and user_input not in choices:
            print(f"Invalid choice. Please choose from: {', '.join(choices)}")
            continue

        # Validate path existence if required
        if is_path and user_input and not os.path.exists(user_input):
            print(f"The path '{user_input}' does not exist. Please try again.")
            continue

        return user_input


def load_genome_configs(assay: Assay) -> Dict[str, GenomeConfig]:
    """Load genome configurations from the config file."""
    config_path = (
        Path(os.getenv("SEQNADO_CONFIG", Path.home()))
        / ".config/seqnado/genome_config.json"
    )
    if not config_path.exists():
        logger.error("Genome config not found. Run 'seqnado-init' first.")
        sys.exit(1)

    with open(config_path, "r") as f:
        all_genome_configs = json.load(f)

    genome_configs: dict[str, GenomeConfig] = dict()
    for genome_name, config_data in all_genome_configs.items():
        try:
            # Ensure required fields are present
            config_data["name"] = genome_name

            # Select appropriate index based on assay requirements
            if assay in [Assay.RNA]:
                index = STARIndex(prefix=config_data.get("star_index"))
            else:
                index = BowtieIndex(prefix=config_data.get("bt2_index"))

            config_data["index"] = index
            genome_configs[genome_name] = GenomeConfig(**config_data)
        except Exception as e:
            # Skip invalid genome configs (e.g., from user's personal config with placeholder paths)
            logger.debug(f"Skipping invalid genome config '{genome_name}': {e}")
            continue

    return genome_configs


def get_project_config() -> ProjectConfig:
    """Get project configuration from user input."""
    username = os.getenv("USER", "unknown_user")
    today = datetime.strftime(datetime.today(), "%Y-%m-%d")

    project_name = get_user_input(
        "Project name?", default=f"{username}_project"
    ).replace(" ", "_")

    project_dir = get_user_input(
        "Project directory?", default=os.getcwd(), is_path=False
    )

    return ProjectConfig(name=project_name, date=today, directory=Path(project_dir))


def get_qc_config() -> QCConfig:
    """Get QC configuration from user input."""
    run_fastq_screen = get_user_input(
        "Perform FastQScreen?", default="no", is_boolean=True
    )
    calculate_library_complexity = get_user_input(
        "Calculate library complexity?", default="no", is_boolean=True
    )
    calculate_fraction_of_reads_in_peaks = get_user_input(
        "Calculate Fraction of Reads in Peaks (FRiP)?", default="no", is_boolean=True
    )
    return QCConfig(
        run_fastq_screen=run_fastq_screen,
        calculate_library_complexity=calculate_library_complexity,
        calculate_fraction_of_reads_in_peaks=calculate_fraction_of_reads_in_peaks,
    )


def select_genome_config(
    genome_configs: Dict[str, GenomeConfig], assay: Assay = None
) -> GenomeConfig:
    """Allow user to select a genome configuration.

    If assay is provided and an assay-specific genome config exists (e.g., 'meth'),
    it will be used as the default choice.
    """
    available_genomes = list(genome_configs.keys())

    # If assay is provided, check if there's an assay-specific config
    default_genome = available_genomes[0] if available_genomes else ""
    if assay:
        assay_name = assay.value.lower()
        if assay_name in genome_configs:
            default_genome = assay_name

    genome_name = get_user_input(
        f"Genome? (Available: {', '.join(available_genomes)})",
        default=default_genome,
    )

    while genome_name not in genome_configs:
        print(
            f"Genome '{genome_name}' is not configured. Please choose from: {', '.join(available_genomes)}"
        )
        genome_name = get_user_input("Genome?", default=default_genome)

    return genome_configs[genome_name]


def get_bigwig_config(assay: Assay) -> Optional[BigwigConfig]:
    """Get bigwig configuration if user wants to create bigwigs."""
    make_bigwigs = get_user_input("Make Bigwigs?", default="yes", is_boolean=True)

    if not make_bigwigs:
        return None

    # Multi-select for pileup methods
    pileup_methods = get_user_input(
        "Bigwig method(s) (comma-separated for multiple):",
        choices=[m.value for m in PileupMethod],
        default="deeptools" if not assay == Assay.MCC else PileupMethod.BAMNADO.value,
        multi_select=True,
    )

    # Convert single string to list if default was used
    if isinstance(pileup_methods, str):
        pileup_methods = [pileup_methods]

    binsize = get_user_input("Binsize for bigwigs:", default="10", required=False)

    return BigwigConfig(
        pileup_method=[PileupMethod(m) for m in pileup_methods], binsize=binsize
    )


def get_plotting_config() -> Optional[PlottingConfig]:
    """Get plotting configuration if user wants plotting."""
    perform_plotting = get_user_input(
        "Perform plotting?", default="no", is_boolean=True
    )

    if not perform_plotting:
        return None

    coordinates = get_user_input(
        "Path to bed file with coordinates for plotting?", required=False, is_path=False
    )

    genes = get_user_input("Path to bed file with genes?", required=False, is_path=False)

    return PlottingConfig(coordinates=coordinates, genes=genes)


def get_peak_calling_config(assay: Assay) -> Optional[PeakCallingConfig]:
    """Get peak calling configuration for assays that support it."""
    if assay not in [Assay.CHIP, Assay.ATAC, Assay.CAT, Assay.MCC]:
        return None

    call_peaks = get_user_input("Call peaks?", default="yes", is_boolean=True)

    if not call_peaks:
        return None

    # Set default based on assay
    default_method = {
        Assay.CHIP: "lanceotron",
        Assay.ATAC: "lanceotron",
        Assay.CAT: "seacr",
    }.get(assay, "lanceotron")

    # Multi-select for peak calling methods
    peak_calling_methods = get_user_input(
        "Peak calling method(s) (comma-separated for multiple):",
        choices=[m.value for m in PeakCallingMethod],
        default=default_method,
        multi_select=True,
    )

    # Convert single string to list if default was used
    if isinstance(peak_calling_methods, str):
        peak_calling_methods = [peak_calling_methods]

    consensus_counts = get_user_input(
        "Generate consensus counts from Design merge column?",
        default="no",
        is_boolean=True,
    )

    # Motif analysis questions
    run_motif_analysis = get_user_input(
        "Run motif analysis on called peaks?",
        default="no",
        is_boolean=True,
    )

    motif_method = None

    if run_motif_analysis:
        motif_methods = get_user_input(
            "Motif analysis method(s) (comma-separated for multiple):",
            choices=[m.value for m in MotifMethod],
            default="homer",
            multi_select=True,
        )

        # Convert single string to list if default was used
        if isinstance(motif_methods, str):
            motif_methods = [motif_methods]

        motif_method = [MotifMethod(m) for m in motif_methods]

    return PeakCallingConfig(
        method=[PeakCallingMethod(m) for m in peak_calling_methods],
        consensus_counts=consensus_counts,
        run_motif_analysis=run_motif_analysis,
        motif_method=motif_method,
    )


def get_spikein_config(assay: Assay) -> Optional[SpikeInConfig]:
    """Get spike-in configuration for assays that support it."""

    spikein = get_user_input("Do you have spikein?", default="no", is_boolean=True)

    if not spikein:
        return None

    normalisation_method = get_user_input(
        "Normalisation method(s) (comma-separated for multiple)?",
        choices=[m.value for m in SpikeInMethod],
        default="deseq2" if assay == Assay.RNA else "orlando",
        multi_select=True,
    )

    # Convert single string to list if default was used
    if isinstance(normalisation_method, str):
        normalisation_method = [normalisation_method]

    reference_genome = get_user_input("Reference genome:", default="hg38")
    spikein_genome = get_user_input(
        "Spikein genome:", default="dm6" if assay != Assay.RNA else "spikein_rna"
    )

    # Ask for control genes if using DESeq2 or edgeR methods
    control_genes = None
    if any(m in ["deseq2", "edger"] for m in normalisation_method):
        control_genes_input = get_user_input(
            "Spike-in control gene names (comma-separated):",
            default="AmpR,Cas9_3p,Cas9_5p",
        )
        if control_genes_input:
            control_genes = [g.strip() for g in control_genes_input.split(",")]

    return SpikeInConfig(
        method=[SpikeInMethod(m) for m in normalisation_method],
        endogenous_genome=reference_genome,
        exogenous_genome=spikein_genome,
        control_genes=control_genes,
    )


def get_ucsc_hub_config() -> Optional[UCSCHubConfig]:
    """Get UCSC hub configuration if user wants it."""
    make_ucsc_hub = get_user_input("Make UCSC hub?", default="no", is_boolean=True)

    if not make_ucsc_hub:
        return None

    username = os.getenv("USER", "unknown_user")

    directory = get_user_input("UCSC hub directory:", default="seqnado_output/hub/")

    email = get_user_input(
        "What is your email address?", default=f"{username}@example.com"
    )
    genome_name = get_user_input("Genome name for UCSC hub?", default="hg38")
    color_by_input = get_user_input("Color by (for UCSC hub):", default="samplename")
    # Convert to list if it's a string
    color_by = [color_by_input] if isinstance(color_by_input, str) else color_by_input

    return UCSCHubConfig(
        directory=directory, email=email, color_by=color_by, genome_name=genome_name
    )


def get_ml_dataset_config(assay: Assay) -> Optional[MLDatasetConfig]:
    """Get ML dataset configuration for supported assays."""
    if assay not in [Assay.ATAC, Assay.CHIP, Assay.CAT]:
        return None

    make_dataset = get_user_input("Make dataset for ML?", default="no", is_boolean=True)

    if not make_dataset:
        return None

    use_regions = get_user_input(
        "Use regions BED file?", default="yes", is_boolean=True
    )

    if use_regions:
        regions_bed = get_user_input(
            "Path to regions BED file:", default="path/to/regions.bed", is_path=False
        )
        # Skip path validation for placeholder paths during interactive config
        is_placeholder = regions_bed.startswith("path/to/")
        return MLDatasetConfig.model_validate(
            {"regions_bed": Path(regions_bed)},
            context={"skip_path_validation": is_placeholder}
        )
    else:
        binsize = int(get_user_input("Binsize for dataset:", default="1000"))
        return MLDatasetConfig(binsize=binsize)


def get_rna_quantification_config() -> Optional[RNAQuantificationConfig]:
    """Get RNA quantification configuration."""
    method = get_user_input(
        "Quantification method:",
        choices=[m.value for m in QuantificationMethod],
        default="feature_counts",
    )

    salmon_index = None
    if method == "salmon":
        salmon_index = get_user_input(
            "Salmon index path:", default="path/to/salmon_index"
        )

    run_deseq2 = get_user_input("Run DESeq2?", default="no", is_boolean=True)

    # Skip path validation for placeholder paths during interactive config
    is_placeholder = salmon_index and salmon_index.startswith("path/to/")
    return RNAQuantificationConfig.model_validate(
        {
            "method": QuantificationMethod(method),
            "salmon_index": salmon_index,
            "run_deseq2": run_deseq2,
        },
        context={"skip_path_validation": is_placeholder}
    )


def get_snp_calling_config() -> Optional[SNPCallingConfig]:
    """Get SNP calling configuration."""
    call_snps = get_user_input("Call SNPs?", default="yes", is_boolean=True)

    if not call_snps:
        return None

    method = get_user_input(
        "SNP calling method:",
        choices=[m.value for m in SNPCallingMethod],
        default="bcftools",
    )

    annotate_snps = get_user_input("Annotate SNPs?", default="no", is_boolean=True)

    snp_database = None
    if annotate_snps:
        snp_database = get_user_input(
            "Path to SNP database:", default="path/to/snp_database", is_path=False
        )

    # Skip path validation for placeholder paths during interactive config
    is_placeholder = snp_database and snp_database.startswith("path/to/")
    return SNPCallingConfig.model_validate(
        {
            "method": SNPCallingMethod(method),
            "annotate_snps": annotate_snps,
            "snp_database": snp_database,
        },
        context={"skip_path_validation": is_placeholder}
    )


def get_mcc_config() -> Optional[MCCConfig]:
    """Get MCC (Capture-C) configuration."""
    viewpoints = get_user_input(
        "Path to viewpoints file:", default="path/to/viewpoints.bed", is_path=False
    )

    resolutions_str = get_user_input(
        "Resolutions for MCC cooler files (comma-separated):", default="100,1000"
    )

    resolutions = [int(r.strip()) for r in resolutions_str.split(",")]

    # Skip path validation for placeholder paths during interactive config
    is_placeholder = viewpoints.startswith("path/to/")
    return MCCConfig.model_validate(
        {"viewpoints": Path(viewpoints), "resolutions": resolutions},
        context={"skip_path_validation": is_placeholder}
    )


def get_methylation_config() -> Optional[MethylationConfig]:
    """Get methylation calling configuration."""
    call_methylation = get_user_input(
        "Call methylation?", default="no", is_boolean=True
    )
    spikein_genomes: list[str] = []
    if call_methylation:
        spikein_genomes_input = get_user_input(
            "Spike-in genomes (comma-separated):",
            default="Lambda,250bp-v1,2kb-unmod",
            required=False,
        )
        if spikein_genomes_input:
            spikein_genomes = [
                genome.strip()
                for genome in spikein_genomes_input.split(",")
                if genome.strip()
            ]
    # Always ensure spikein_genomes is a list, never None
    if not call_methylation:
        return None

    methylation_assay = get_user_input(
        "Methylation assay:",
        choices=[m.value for m in MethylationMethod],
        default="taps",
    )

    return MethylationConfig(
        method=MethylationMethod(methylation_assay),
        spikein_genomes=spikein_genomes or [],
    )


def build_assay_config(
    assay: Assay, genome_config: GenomeConfig
) -> Optional[AssaySpecificConfig]:
    """Build assay-specific configuration based on the assay type."""

    # Collect standard configuration for assays that need it
    if assay not in [Assay.CRISPR, Assay.SNP]:
        bigwigs = get_bigwig_config(assay=assay)
        plotting = get_plotting_config()
        ucsc_hub = get_ucsc_hub_config()
        create_heatmaps = get_user_input(
            "Make heatmaps?", default="no", is_boolean=True
        )
        geo_files = get_user_input(
            "Generate GEO submission files?", default="no", is_boolean=True
        )

        base_config = {
            "genome": genome_config,
            "bigwigs": bigwigs,
            "plotting": plotting,
            "ucsc_hub": ucsc_hub,
            "create_heatmaps": create_heatmaps,
            "create_geo_submission_files": geo_files,
        }

    match assay:
        case Assay.ATAC:
            tn5_shift = get_user_input(
                "Shift ATAC reads?", default="yes", is_boolean=True
            )
            peak_calling = get_peak_calling_config(assay)
            dataset_for_ml = get_ml_dataset_config(assay)

            return ATACAssayConfig(
                **base_config,
                tn5_shift=tn5_shift,
                peak_calling=peak_calling,
                dataset_for_ml=dataset_for_ml,
            )

        case Assay.CHIP:
            spikein = get_spikein_config(assay)
            peak_calling = get_peak_calling_config(assay)
            dataset_for_ml = get_ml_dataset_config(assay)

            return ChIPAssayConfig(
                **base_config,
                spikein=spikein,
                peak_calling=peak_calling,
                dataset_for_ml=dataset_for_ml,
            )

        case Assay.CAT:
            tn5_shift = get_user_input(
                "Shift CAT reads?", default="no", is_boolean=True
            )
            spikein = get_spikein_config(assay)
            peak_calling = get_peak_calling_config(assay)
            dataset_for_ml = get_ml_dataset_config(assay)

            return CATAssayConfig(
                **base_config,
                tn5_shift=tn5_shift,
                spikein=spikein,
                peak_calling=peak_calling,
                dataset_for_ml=dataset_for_ml,
            )

        case Assay.RNA:
            spikein = get_spikein_config(assay)
            rna_quantification = get_rna_quantification_config()

            return RNAAssayConfig(
                **base_config,
                spikein=spikein,
                rna_quantification=rna_quantification,
            )

        case Assay.MCC:
            mcc = get_mcc_config()
            peak_calling = get_peak_calling_config(assay)
            return MCCAssayConfig(**base_config, mcc=mcc, peak_calling=peak_calling)

        case Assay.METH:
            methylation = get_methylation_config()
            return MethylationAssayConfig(
                genome=genome_config,
                bigwigs=bigwigs,
                plotting=plotting,
                ucsc_hub=None,
                create_heatmaps=create_heatmaps,
                create_geo_submission_files=geo_files,
                methylation=methylation,
            )

        case Assay.SNP:
            # SNP assays don't generate bigwigs, so UCSC hub is not supported
            geo_files = get_user_input(
                "Generate GEO submission files?", default="no", is_boolean=True
            )
            snp_calling = get_snp_calling_config()
            return SNPAssayConfig(
                genome=genome_config,
                bigwigs=None,
                plotting=None,
                ucsc_hub=None,
                create_heatmaps=False,
                create_geo_submission_files=geo_files,
                snp_calling=snp_calling,
            )

        case Assay.CRISPR:
            geo_files = get_user_input(
                "Generate GEO submission files?", default="no", is_boolean=True
            )
            use_mageck = get_user_input(
                "Use MAGeCK for guide RNA analysis?", default="no", is_boolean=True
            )
            return CRISPRAssayConfig(
                genome=genome_config,
                bigwigs=None,
                plotting=None,
                ucsc_hub=None,
                create_heatmaps=False,
                create_geo_submission_files=geo_files,
                use_mageck=use_mageck,
            )

        case _:
            raise ValueError(f"Unsupported assay type: {assay}")


def build_default_assay_config(
    assay: Assay, genome_config: GenomeConfig
) -> Optional[AssaySpecificConfig]:
    """Build a default assay-specific configuration for non-interactive mode."""
    # Set common defaults
    bigwigs = BigwigConfig(pileup_method=[PileupMethod.DEEPTOOLS], binsize=10)

    # Set default plotting coordinates to None
    # Tests will set this explicitly to avoid referencing package test data
    default_coordinates = None

    plotting = PlottingConfig(coordinates=default_coordinates)
    ucsc_hub = UCSCHubConfig(
        directory="seqnado_output/hub/",
        genome=genome_config.name,
        email="user@example.com",
        genome_name=genome_config.name,
    )
    create_heatmaps = False
    geo_files = False

    base_config = {
        "genome": genome_config,
        "bigwigs": bigwigs,
        "plotting": plotting,
        "ucsc_hub": ucsc_hub,
        "create_heatmaps": create_heatmaps,
        "create_geo_submission_files": geo_files,
    }

    match assay:
        case Assay.ATAC:
            tn5_shift = True
            peak_calling = PeakCallingConfig(
                method=[PeakCallingMethod.LANCEOTRON], consensus_counts=False
            )
            dataset_for_ml = MLDatasetConfig(binsize=1000)

            return ATACAssayConfig(
                **base_config,
                tn5_shift=tn5_shift,
                peak_calling=peak_calling,
                dataset_for_ml=dataset_for_ml,
            )

        case Assay.CHIP:
            spikein = None
            peak_calling = PeakCallingConfig(
                method=[PeakCallingMethod.LANCEOTRON], consensus_counts=False
            )
            dataset_for_ml = MLDatasetConfig(binsize=1000)

            return ChIPAssayConfig(
                **base_config,
                spikein=spikein,
                peak_calling=peak_calling,
                dataset_for_ml=dataset_for_ml,
            )
        
        case Assay.CAT:
            tn5_shift = False
            spikein = None
            peak_calling = PeakCallingConfig(
                method=[PeakCallingMethod.SEACR], consensus_counts=False
            )
            dataset_for_ml = MLDatasetConfig(binsize=1000)
            return CATAssayConfig(
                **base_config,
                tn5_shift=tn5_shift,
                spikein=spikein,
                peak_calling=peak_calling,
                dataset_for_ml=dataset_for_ml,
            )
        
        case Assay.RNA:
            rna_quantification = RNAQuantificationConfig(
                method=QuantificationMethod.FEATURE_COUNTS,
                salmon_index=None,
                run_deseq2=False,
            )
            return RNAAssayConfig(**base_config, rna_quantification=rna_quantification)
        
        case Assay.SNP:
            # SNP assays don't use bigwigs, plotting, heatmaps, or UCSC hub
            base_config_snp = {
                k: v for k, v in base_config.items() if k not in ("bigwigs", "plotting", "ucsc_hub")
            }
            base_config_snp["bigwigs"] = None
            base_config_snp["plotting"] = None
            base_config_snp["ucsc_hub"] = None
            base_config_snp["create_heatmaps"] = False
            snp_calling = SNPCallingConfig(
                method=SNPCallingMethod.BCFTOOLS,
                annotate_snps=False,
                snp_database=None,
            )
            return SNPAssayConfig(**base_config_snp, snp_calling=snp_calling)
        
        case Assay.MCC:
            # MCC assays don't use UCSC hub
            base_config_mcc = {k: v for k, v in base_config.items() if k != "ucsc_hub"}
            base_config_mcc["ucsc_hub"] = None
            # Allow override for test environment
            import os

            viewpoints_path = os.environ.get(
                "SEQNADO_MCC_VIEWPOINTS", "path/to/viewpoints.bed"
            )
            # Skip path validation only for placeholder paths
            is_placeholder = viewpoints_path.startswith("path/to/")
            mcc = MCCConfig.model_validate(
                {
                    "viewpoints": Path(viewpoints_path),
                    "resolutions": [100, 1000],
                },
                context={"skip_path_validation": is_placeholder}
            )
            return MCCAssayConfig(**base_config_mcc, mcc=mcc)
        
        case Assay.METH:
            # Methylation assays don't use UCSC hub
            base_config_meth = {k: v for k, v in base_config.items() if k != "ucsc_hub"}
            base_config_meth["ucsc_hub"] = None
            methylation = MethylationConfig(method=MethylationMethod.TAPS)
            return MethylationAssayConfig(**base_config_meth, methylation=methylation)
        
        case Assay.CRISPR:
            # CRISPR assays don't use UCSC hub or bigwigs
            base_config_crispr = {
                k: v for k, v in base_config.items() if k not in ("ucsc_hub", "bigwigs")
            }
            base_config_crispr["ucsc_hub"] = None
            base_config_crispr["bigwigs"] = None
            return CRISPRAssayConfig(**base_config_crispr)
        
        case _:
            raise ValueError(f"Unsupported assay type: {assay}")


def build_workflow_config(assay: Assay, seqnado_version: str) -> SeqnadoConfig:
    """Build complete workflow configuration from user input.
    This function orchestrates the collection of all necessary configuration
    details from the user and constructs a SeqnadoConfig object.

    Args:
        assay (Assay): The assay type for the workflow.
        seqnado_version (str): The current version of SeqNado.
    Returns:
        SeqnadoConfig: The complete workflow configuration.
    Raises:
        SystemExit: If there are any configuration errors.
    """

    logger.info(
        f"Building workflow configuration for {assay.value} assay with SeqNado version {seqnado_version}"
    )
    logger.info(f"Current Working Directory: {os.getcwd()}")
    logger.debug(f"Assay: {assay}, SeqNado version: {seqnado_version}")

    # Load available genome configurations
    genome_configs = load_genome_configs(assay)

    # Get project configuration
    project = get_project_config()

    # Select genome configuration
    genome = select_genome_config(genome_configs, assay=assay)

    # Get metadata path
    metadata_path = get_user_input(
        "Path to metadata file:",
        default=f"metadata_{assay.clean_name}.csv",
        is_path=False,
    )

    # Build assay-specific configuration
    assay_config = build_assay_config(assay, genome)

    try:
        workflow_config = SeqnadoConfig(
            assay=assay,
            project=project,
            genome=genome,
            metadata=Path(metadata_path),
            assay_config=assay_config,
        )
        return workflow_config

    except ValidationError as e:
        logger.error(f"Configuration validation error: {e}")
        sys.exit(1)


def build_default_workflow_config(assay: Assay) -> SeqnadoConfig:
    """Build a default workflow configuration for non-interactive mode."""
    logger.info(f"Building default workflow configuration for {assay.value} assay")

    # Load available genome configurations
    genome_configs = load_genome_configs(assay)

    if not genome_configs:
        logger.error("No genome configurations available. Cannot build default config.")
        sys.exit(1)

    # Use the first available genome configuration as default
    default_genome = next(iter(genome_configs.values()))

    # Create a default project configuration
    default_project = ProjectConfig(
        name="default_project",
        date=datetime.strftime(datetime.today(), "%Y-%m-%d"),
        directory=Path(os.getcwd()),
    )

    # Build assay-specific configuration with defaults
    assay_config = build_default_assay_config(assay, default_genome)

    try:
        workflow_config = SeqnadoConfig(
            assay=assay,
            project=default_project,
            genome=default_genome,
            metadata=f"metadata_{assay.clean_name}.csv",
            assay_config=assay_config,
        )
        return workflow_config

    except ValidationError as e:
        logger.error(f"Configuration validation error: {e}")
        sys.exit(1)


def render_config(
    template: Path,
    workflow_config: SeqnadoConfig,
    outfile: Path,
    all_options: bool = False,
) -> None:
    """Render the workflow configuration to a file."""

    env = jinja2.Environment(loader=jinja2.FileSystemLoader(template.parent))
    template = env.get_template(template.name)

    # Convert the Pydantic model to a dictionary for rendering
    # Always include fields with None values to ensure required fields like ucsc_hub are present
    config_dict = workflow_config.model_dump(mode="json", exclude_none=False)

    try:
        rendered_content = template.render(**config_dict)
    except jinja2.TemplateError as e:
        logger.error(f"Template rendering error: {e}")
        sys.exit(1)

    with open(outfile, "w") as f:
        f.write(rendered_content)


def build_multiomics_config(
    seqnado_version: str, interactive: bool = True, assays: list[str] | None = None
) -> tuple["MultiomicsConfig", dict[str, SeqnadoConfig]]:
    """Build multiomics configuration with multiple assays.

    Args:
        seqnado_version: The current version of SeqNado
        interactive: Whether to prompt for user input
        assays: List of assay names to include (only used in non-interactive mode)

    Returns:
        tuple: (MultiomicsConfig, dict of assay_name -> SeqnadoConfig)
    """
    # Import here to avoid circular import
    from seqnado.config.multiomics import MultiomicsConfig

    logger.info("Building multiomics configuration")

    # Get list of assays to include (exclude multiomics itself)
    available_assays = [assay.clean_name for assay in Assay.non_multiomics_assays()]

    if interactive:
        print("\nAvailable assays:")
        for i, assay in enumerate(available_assays, 1):
            print(f"  {i}. {assay}")

        assay_selection = get_user_input(
            "Enter assay numbers separated by commas (e.g., 1,3,5) or assay names (e.g., atac,chip,rna)",
            required=True,
        )

        # Parse selection
        selected_assays = []
        for item in assay_selection.split(","):
            item = item.strip()
            # Try as number first
            try:
                idx = int(item) - 1
                if 0 <= idx < len(available_assays):
                    selected_assays.append(available_assays[idx])
                else:
                    logger.warning(f"Skipping invalid selection: {item}")
            except ValueError:
                # Try as assay name
                if item in available_assays:
                    selected_assays.append(item)
                else:
                    logger.warning(f"Skipping unknown assay: {item}")

        if not selected_assays:
            logger.error("No valid assays selected")
            sys.exit(1)

        logger.info(f"Selected assays: {', '.join(selected_assays)}")

        # Get multiomics-specific settings
        output_dir = get_user_input(
            "Base output directory for all assays", default="seqnado_output/"
        )

        create_heatmaps = get_user_input(
            "Generate multiomics heatmaps?", default="yes", is_boolean=True
        )

        create_dataset = get_user_input(
            "Generate ML-ready dataset?", default="yes", is_boolean=True
        )

        create_summary = get_user_input(
            "Generate summary report?", default="yes", is_boolean=True
        )

        # Only ask for regions_bed and binsize if creating dataset
        regions_bed = None
        binsize = None
        if create_dataset:
            regions_bed = get_user_input(
                "BED file with regions of interest (optional, press Enter to skip)",
                required=False,
                is_path=False,
            )

            binsize_input = get_user_input(
                "Binsize for genome-wide analysis (optional, press Enter to skip)",
                required=False,
            )

            binsize = int(binsize_input) if binsize_input and binsize_input.isdigit() else None

    else:
        # Non-interactive mode: still prompt for assay selection, but use defaults for everything else
        logger.info("Non-interactive mode: prompting for assays, using defaults for other settings")

        print("\nAvailable assays:")
        for i, assay in enumerate(available_assays, 1):
            print(f"  {i}. {assay}")

        assay_selection = get_user_input(
            "Enter assay numbers separated by commas (e.g., 1,3,5) or assay names (e.g., atac,chip,rna)",
            required=True,
        )

        # Parse the selection (numbers or names)
        selected_assays = []
        for item in assay_selection.split(","):
            item = item.strip()
            if item.isdigit():
                # User entered a number
                idx = int(item) - 1
                if 0 <= idx < len(available_assays):
                    selected_assays.append(available_assays[idx])
                else:
                    logger.warning(f"Invalid assay number: {item}")
            else:
                # User entered an assay name
                if item in available_assays:
                    selected_assays.append(item)
                else:
                    logger.warning(f"Invalid assay name: {item}")

        if not selected_assays:
            logger.error("No valid assays selected")
            sys.exit(1)

        logger.info(f"Selected assays: {', '.join(selected_assays)}")

        # Use defaults for all other settings
        output_dir = "seqnado_output/"
        create_heatmaps = True
        create_dataset = True
        create_summary = True
        regions_bed = None
        binsize = None

    # Build individual assay configs
    assay_configs = {}
    for assay_name in selected_assays:
        assay = Assay.from_clean_name(assay_name)
        logger.info(f"\n=== Configuring {assay.value} assay ===")

        if interactive:
            # Ask if user wants to configure this assay now or use defaults
            configure_now = get_user_input(
                f"Configure {assay_name} now?", default="yes", is_boolean=True
            )

            if configure_now:
                config = build_workflow_config(assay, seqnado_version)
            else:
                logger.info(f"Using default configuration for {assay_name}")
                config = build_default_workflow_config(assay)
        else:
            config = build_default_workflow_config(assay)

        assay_configs[assay_name] = config

    # Create MultiomicsConfig
    multiomics_config = MultiomicsConfig(
        assays=selected_assays,
        output_dir=output_dir,
        create_heatmaps=create_heatmaps,
        create_dataset=create_dataset,
        create_summary=create_summary,
        regions_bed=Path(regions_bed) if regions_bed else None,
        binsize=binsize,
    )

    return multiomics_config, assay_configs


def render_multiomics_configs(
    multiomics_config: "MultiomicsConfig",
    assay_configs: dict[str, SeqnadoConfig],
    template: Path,
    output_dir: Path,
) -> List[Path]:
    """Render all config files for multiomics analysis.

    Args:
        multiomics_config: The multiomics configuration
        assay_configs: Dictionary of assay name -> SeqnadoConfig
        template: Path to the config template file
        output_dir: Directory to write config files to

    Returns:
        List of paths to generated config files
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    generated_files = []

    # Render individual assay configs
    for assay_name, config in assay_configs.items():
        config_file = output_dir / f"config_{assay_name}.yaml"
        render_config(template, config, config_file, all_options=False)
        generated_files.append(config_file)
        logger.success(f"Generated {config_file}")

    # Render multiomics config
    multiomics_file = output_dir / "config_multiomics.yaml"

    # Create a simple YAML representation of multiomics config
    multiomics_dict = multiomics_config.model_dump(mode="json", exclude_none=True)

    import yaml

    with open(multiomics_file, "w") as f:
        yaml.dump(multiomics_dict, f, default_flow_style=False, sort_keys=False)

    generated_files.append(multiomics_file)
    logger.success(f"Generated {multiomics_file}")

    return generated_files
