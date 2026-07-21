import sys
import inspect
from pathlib import Path
import tracknado as tn
from seqnado.workflow.helpers.hub import (
    format_hub_labels,
    load_design_metadata,
    make_seqnado_extractor,
)
from loguru import logger


def create_hub_with_tracknado():
    # Configure logger to write to snakemake log file if available
    if "snakemake" in globals():
        logger.remove()
        logger.add(
            snakemake.log[0],
            format="{time} {level} {message}",
            level="INFO",
        )
        logger.add(sys.stderr, format="{time} {level} {message}", level="ERROR")

    try:
        logger.info("🌪️  TrackNado: Generating UCSC Genome Browser hub...")
        logger.info(
            f"Assay: {snakemake.params.assay}, Files: {len(snakemake.input.data)}"
        )

        # 1. Initialize the HubBuilder with the input files
        # We use the fluent API to configure everything in a single chain
        supergroup_by = snakemake.params.supergroup_by
        subgroup_by = snakemake.params.subgroup_by
        overlay_by = snakemake.params.overlay_by
        color_by = snakemake.params.color_by

        # The assay is known here (not guessable from the path), so bind the
        # assay-aware extractor to it - see seqnado.workflow.helpers.hub.
        assay = snakemake.params.assay
        assay_name = assay.clean_name if hasattr(assay, "clean_name") else str(assay)
        design_metadata = load_design_metadata(snakemake.input.metadata)

        builder = (
            tn.HubBuilder()
            .add_tracks(snakemake.input.data)
            # 2. Extract grouping metadata (method/norm/antibody/strand/viewpoint)
            # from SeqNado's real output layout.
            .with_metadata_extractor(
                make_seqnado_extractor(assay_name, design_metadata)
            )
        )

        if supergroup_by:
            builder = builder.group_by(*supergroup_by, as_supertrack=True)

        if subgroup_by:
            builder = builder.group_by(*subgroup_by)

        if overlay_by:
            builder = builder.overlay_by(*overlay_by)

        if color_by:
            builder = builder.color_by(color_by)
        else:
            builder = builder.color_by(
                "samplename"
            )  # Default coloring by method if none provided

        # If we have a custom genome, set it here
        if (
            hasattr(snakemake.params, "custom_genome")
            and snakemake.params.custom_genome
        ):
            # custom_genome=lambda wc: True if CONFIG.assay_config.ucsc_hub.two_bit else False,
            # genome_twobit=CONFIG.assay_config.ucsc_hub.two_bit,
            # genome_organism=CONFIG.assay_config.ucsc_hub.organism,
            # genome_default_position=CONFIG.assay_config.ucsc_hub.default_position,
            builder = builder.with_custom_genome(
                name=snakemake.params.genome,
                twobit_file=snakemake.params.genome_twobit,
                organism=snakemake.params.genome_organism,
                default_position=snakemake.params.genome_default_position,
            )

        # 4. Build and stage the hub
        # This automatically generates hub.txt, genomes.txt, trackDb.txt
        # and saves a 'tracknado_config.json' sidecar for easy merging later.
        outdir = Path(str(snakemake.output.hub)).parent
        hub = builder.build(
            name=snakemake.params.hub_name,
            genome=snakemake.params.genome,
            outdir=outdir,
            hub_email=snakemake.params.hub_email,
            description_html=Path(snakemake.input.report)
            if hasattr(snakemake.input, "report")
            else None,
        )
        format_hub_labels(hub, supergroup_by)

        # remove_existing so reruns don't leave stale tracks behind (e.g. a
        # sample removed from the design) - the hub dir's bigWig/bigBed
        # symlinks aren't Snakemake-tracked outputs, only the *.txt files are.
        if "remove_existing" in inspect.signature(hub.stage_hub).parameters:
            hub.stage_hub(remove_existing=True)
        else:
            logger.warning(
                "Installed TrackNado does not support removing stale staged tracks; "
                "staging the current hub contents."
            )
            hub.stage_hub()

        logger.info(f"✓ Hub successfully generated in {outdir}")
        logger.info("💡 You can merge multiple hubs later using 'tracknado merge'")

    except Exception as e:
        logger.error("=" * 80)
        logger.error(f"❌ ERROR: Failed to generate UCSC hub: {e}")
        import traceback

        logger.error(f"Traceback:\n{traceback.format_exc()}")
        logger.error("=" * 80)
        raise


if __name__ == "__main__":
    create_hub_with_tracknado()
