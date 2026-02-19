import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotnado.api as pn
import pyranges as pr
import seaborn as sns
from loguru import logger


def snakemake_setup():
    if "snakemake" not in globals():
        raise RuntimeError("This script must be run via Snakemake.")
    log_file = snakemake.log[0]
    logger.remove()
    logger.add(log_file, level="DEBUG")
    logger.add(sys.stderr, level="ERROR")
    assay = snakemake.params.assay
    input_data = snakemake.input.data
    output_plots = snakemake.output.plots
    regions = snakemake.params.regions
    outdir = Path(snakemake.output.template).parent
    plotting_format = snakemake.params.plotting_format
    genes = snakemake.params.genes if hasattr(snakemake.params, "genes") else None
    template = snakemake.output.template
    return (
        assay,
        input_data,
        output_plots,
        regions,
        outdir,
        plotting_format,
        genes,
        template,
    )


def configure_matplotlib():
    """Set matplotlib style parameters."""
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["svg.fonttype"] = "none"



def parse_bigwig_metadata(path: Path) -> tuple[str, str]:
    """
    Parse pileup method and normalisation from a bigwig path.

    Handles all path structures produced by SeqNado:
    - bigwigs/{method}/{scale}/{name}.bigWig
    - bigwigs/{method}/merged/{scale}/{name}.bigWig
    - bigwigs/{method}/spikein/{spikein_method}/{name}.bigWig
    - bigwigs/{method}/merged/spikein/{spikein_method}/{name}.bigWig

    Returns:
        Tuple of (method, normalisation) e.g. ("deeptools", "merged/spikein/orlando")
    """
    parts = path.parts
    if "bigwigs" not in parts:
        # Fallback for non-standard paths
        return parts[-3], parts[-2]
    bw_idx = parts.index("bigwigs")
    # Everything between bigwigs/ and the filename
    inner = parts[bw_idx + 1 : -1]
    method = inner[0]
    normalisation = "/".join(inner[1:]) if len(inner) > 1 else ""
    return method, normalisation


def load_tracks(input_paths: list, assay: str) -> pd.DataFrame:
    """
    Load track files into a DataFrame with metadata.

    Args:
        input_paths: List of file paths to track files
        assay: Assay type (e.g., 'ChIP')

    Returns:
        DataFrame with track metadata
    """
    logger.info("Loading input tracks...")

    df = pd.DataFrame([Path(p) for p in input_paths], columns=["path"])
    df["name"] = df["path"].apply(lambda x: x.stem)
    df["type"] = df["path"].apply(lambda x: x.suffix)
    df["type"] = pd.Categorical(
        df["type"], categories=[".bigWig", ".bed"], ordered=True
    )

    # Extract metadata from path structure — anchor on the known "bigwigs" dir
    # so that merged/spikein/merged+spikein paths are all handled correctly.
    def _method(p):
        return parse_bigwig_metadata(p)[0]

    def _normalisation(p):
        return parse_bigwig_metadata(p)[1]

    df["method"] = np.where(
        df["type"] != ".bed",
        df["path"].apply(_method),
        df["path"].apply(lambda x: x.parts[-2]),
    )
    df["normalisation"] = np.where(
        df["type"] != ".bed", df["path"].apply(_normalisation), ""
    )

    # Sort and create track names
    df = df.sort_values(by=["name", "type", "method", "normalisation"])
    df["track_name"] = (
        df["name"]
        + "-"
        + df["method"]
        + "-"
        + df["normalisation"]
        + df["type"].astype(str)
    )
    df["track_name"] = df["track_name"].str.replace("-.", ".")

    # Add antibody info for ChIP assays
    if assay == "ChIP":
        df["antibody"] = df["name"].str.split("_").str[-1]

    return df


def get_track_config(track, assay: str) -> tuple:
    """
    Determine track configuration based on file type and assay.

    Args:
        track: DataFrame row containing track info
        assay: Assay type

    Returns:
        Tuple of (track_type, style, autoscaling_group)
    """
    if track.type == ".bed":
        return "bed_simple", None, None

    # BigWig configuration
    style = "stairsfilled"
    if assay == "ChIP":
        autoscaling_group = f"{track.antibody}-{track.method}-{track.normalisation}"
    else:
        autoscaling_group = f"{track.method}-{track.normalisation}"

    return "bigwig", style, autoscaling_group


def create_color_palette(names: list) -> dict:
    """Create a color dictionary for unique sample names."""
    return dict(zip(names, sns.color_palette("tab20", n_colors=len(names))))


def build_figure(df: pd.DataFrame, assay: str, genes_file: str = None) -> pn.Figure:
    """
    Build a Plotnado figure with all tracks.

    Args:
        df: DataFrame with track metadata
        assay: Assay type
        genes_file: Optional path to genes file

    Returns:
        Configured Plotnado Figure object
    """
    logger.info("Creating Plotnado figure...")

    fig = pn.Figure(autospacing=True)
    fig.add_track("scale")

    # Add genes track if provided
    if genes_file:
        fig.add_track(
            "genes",
            file=genes_file,
            gene_style="normal",
            min_gene_length=int(1e3),
            label_y_offset=-75,
            label_loc="right",
            arrow_color="black",
            fontsize=6,
        )

    # Add data tracks
    colors_dict = create_color_palette(df["name"].unique())

    for track in df.itertuples():
        track_type, style, autoscaling_group = get_track_config(track, assay)

        t = pn.TrackWrapper(
            track_type,
            str(track.path),
            name=track.track_name,
            title=track.track_name,
            color=colors_dict[track.name],
            style=style,
            data_range_style="text",
            data_range_location="right",
            label_on_track=True,
            label_loc="left",
            autoscaling_group=autoscaling_group,
        )
        fig.add_track(t)

    return fig


def generate_region_name(region) -> str:
    """Generate a filename-friendly name for a genomic region."""
    if hasattr(region, "Name") and region.Name:
        return region.Name
    return f"{region.Chromosome}-{region.Start}-{region.End}"


def save_plots(fig: pn.Figure, coords: pr.PyRanges, outdir: Path, plotting_format: str):
    """
    Save plots for all regions.

    Args:
        fig: Configured Plotnado figure
        coords: PyRanges object with genomic coordinates
        outdir: Output directory path
        plotting_format: Output format (e.g., 'png', 'svg')
    """
    logger.info(f"Output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    for region in coords.df.itertuples():
        fig_name = generate_region_name(region)
        region_coords = f"{region.Chromosome}:{region.Start}-{region.End}"
        output_file = outdir / f"{fig_name}.{plotting_format}"

        logger.info(f"Saving plot for region {fig_name}: {output_file}")
        fig.save(output=output_file, gr=region_coords)


def main():
    (
        assay,
        input_data,
        output_plots,
        regions,
        outdir,
        plotting_format,
        genes,
        template,
    ) = snakemake_setup()
    configure_matplotlib()

    logger.info("Starting Plotnado visualization")
    logger.info(f"Processing {assay} assay")
    logger.debug(f"Input files: {input_data}")
    logger.debug(f"Output plots: {output_plots}")
    logger.debug(f"Plotting regions: {regions}")
    logger.debug(f"Output directory: {outdir}")

    # Load and process data
    df = load_tracks(input_data, assay)

    logger.info("Loading plotting regions...")
    coords = pr.read_bed(regions)
    logger.info(f"Found {len(coords)} regions to plot")

    logger.info(f"Output format: {plotting_format}")

    # Build figure and save plots
    fig = build_figure(df, assay, genes)
    save_plots(fig, coords, Path(outdir), plotting_format)

    # Save template
    logger.info("Saving Plotnado template...")
    fig.to_toml(template)
    logger.info("Plotnado visualization complete!")


if __name__ == "__main__":
    main()
