import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotnado.api as pn
import pybedtools
import pyranges as pr
import seaborn as sns
from loguru import logger


def _patch_pybedtools_tabix_detection():
    """
    Monkey-patch pybedtools so that tabix_intervals auto-detects an existing
    .tbi index when a BedTool is created from a .gz path (without going through
    .tabix()).  This is needed because plotnado creates fresh BedTool objects
    from file paths, losing the _tabix attribute set by our prepare_bed_for_tabix
    helper.
    """
    import pysam

    _original = pybedtools.BedTool.tabix_intervals

    def _patched(self, interval, **kwargs):
        # Safely check if _tabix exists and is None
        if getattr(self, "_tabix", None) is None and getattr(self, "fn", None) and self.fn.endswith(".gz"):
            tbi = self.fn + ".tbi"
            if Path(tbi).exists():
                self._tabix = pysam.TabixFile(self.fn)
        return _original(self, interval, **kwargs)

    pybedtools.BedTool.tabix_intervals = _patched


_patch_pybedtools_tabix_detection()


def snakemake_setup():
    if "snakemake" not in globals():
        raise RuntimeError("This script must be run via Snakemake.")
    log_file = snakemake.log[0]
    logger.remove()
    logger.add(log_file, level="DEBUG")
    logger.add(sys.stderr, level="ERROR")
    assay = snakemake.params.assay
    input_data = snakemake.input.data
    peak_files = snakemake.params.peak_files if hasattr(snakemake.params, "peak_files") else []
    output_plots = snakemake.output.plots
    regions = snakemake.params.regions
    outdir = Path(snakemake.output.template).parent
    plotting_format = snakemake.params.plotting_format
    genes = snakemake.params.genes if hasattr(snakemake.params, "genes") else None
    template = snakemake.output.template
    return (
        assay,
        input_data,
        peak_files,
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
    
    # Handle .bed.gz files: stem should remove one .gz, keeping the .bed part
    df["name"] = df["path"].apply(lambda x: str(x).replace(".bed.gz", "").split("/")[-1] if str(x).endswith(".bed.gz") else x.stem)
    
    # Detect file type, treating .bed.gz as .bed
    def get_type(path_str):
        if str(path_str).endswith(".bed.gz"):
            return ".bed"
        return Path(path_str).suffix
    
    df["type"] = df["path"].apply(lambda x: get_type(x))
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


def prepare_bed_for_tabix(bed_path: Path) -> str:
    """Sort, bgzip, and tabix-index a BED file. Returns path to the .gz file."""
    bt = pybedtools.BedTool(str(bed_path))
    bt_tabix = bt.sort().tabix(force=True)
    return bt_tabix.fn


def build_figure(df: pd.DataFrame, assay: str, genes_file: str = None) -> pn.Figure:
    """
    Build a Plotnado figure with all tracks.

    Args:
        df: DataFrame with track metadata
        assay: Assay type
        genes_file: Optional path to genes file (can be .bed or .bed.gz)

    Returns:
        Configured Plotnado Figure object
    """
    logger.info("Creating Plotnado figure...")

    fig = pn.Figure(autospacing=True)
    fig.add_track("scale")

    # Add genes track if provided
    if genes_file and Path(genes_file).exists() and Path(genes_file).stat().st_size > 0:
        # Use genes file directly if it's already indexed (.bed.gz with .tbi)
        # Otherwise prepare it for tabix
        if genes_file.endswith(".bed.gz"):
            # Assume it's pre-indexed if it ends with .gz
            genes_path = genes_file
        else:
            # Index on-the-fly if it's a raw .bed file
            genes_path = prepare_bed_for_tabix(Path(genes_file))
        
        fig.add_track(
            "genes",
            file=genes_path,
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

        # Use pre-indexed .bed.gz files directly, or index on-the-fly if raw .bed
        if track.type == ".bed":
            if track.path.endswith(".bed.gz"):
                # Already indexed by index_individual_peaks rule
                file_path = str(track.path)
            else:
                # Fallback for raw .bed files
                file_path = prepare_bed_for_tabix(track.path)
        else:
            file_path = str(track.path)

        t = pn.TrackWrapper(
            track_type,
            file_path,
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
        peak_files,
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
    logger.debug(f"Peak files: {peak_files}")
    logger.debug(f"Output plots: {output_plots}")
    logger.debug(f"Plotting regions: {regions}")
    logger.debug(f"Output directory: {outdir}")

    # Load and process data (combine bigwigs and pre-indexed peaks)
    all_data = list(input_data) + (peak_files if isinstance(peak_files, list) else [peak_files])
    df = load_tracks(all_data, assay)

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
    try:
        main()
    except Exception as e:
        logger.exception(f"Fatal error during Plotnado execution: {e}")
        raise
