"""Render genome-browser track plots with PlotNado (>= 0.3).

Run by Snakemake inside the PlotNado container (``CONTAINERS.PLOTNADO``), so
only dependencies shipped in that image may be imported here (no pyranges or
seaborn).
"""

from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotnado as pn
from loguru import logger

BED_COLUMNS = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
TRACK_TYPES = [".bigWig", ".bed"]


def snakemake_setup():
    if "snakemake" not in globals():
        raise RuntimeError("This script must be run via Snakemake.")
    logger.remove()
    logger.add(snakemake.log[0], level="DEBUG")
    params = snakemake.params
    input_data = list(snakemake.input.data)
    # Multiomics passes one assay per input file; single-assay rules pass one
    # assay for everything.
    track_assays = getattr(params, "track_assays", None) or [params.assay] * len(input_data)
    return (
        track_assays,
        input_data,
        snakemake.output.plots,
        params.regions,
        params.outdir,
        params.plotting_format,
        getattr(params, "genes", None),
        snakemake.output.template,
    )


def configure_matplotlib():
    """Set matplotlib style parameters."""
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["svg.fonttype"] = "none"


def _track_type(path: Path) -> str:
    """Normalise file type so ``.bed.gz`` peaks are treated as BED."""
    return ".bed" if path.name.endswith((".bed", ".bed.gz")) else path.suffix


def _track_stem(path: Path) -> str:
    name = path.name
    for suffix in (".bed.gz", ".bed", ".bigWig", ".bigwig", ".bw"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def load_tracks(input_paths: list, assays: list[str]) -> pd.DataFrame:
    """Load track files into a DataFrame with metadata parsed from their paths.

    ``assays`` gives the assay of each input path, in the same order.
    """
    logger.info("Loading input tracks...")

    if len(assays) != len(input_paths):
        raise ValueError(
            f"Got {len(assays)} track assays for {len(input_paths)} input files"
        )
    df = pd.DataFrame({"path": [Path(p) for p in input_paths], "assay": assays})
    df["name"] = df["path"].apply(_track_stem)
    df["type"] = pd.Categorical(
        df["path"].apply(_track_type), categories=TRACK_TYPES, ordered=True
    )

    is_bed = df["type"] == ".bed"
    df["normalisation"] = np.where(is_bed, "", df["path"].apply(lambda x: x.parts[-2]))
    df["method"] = np.where(
        is_bed,
        df["path"].apply(lambda x: x.parts[-2]),
        df["path"].apply(lambda x: x.parts[-3]),
    )

    df = df.sort_values(by=["assay", "name", "type", "method", "normalisation"])
    df["track_name"] = (
        df["name"] + "-" + df["method"] + "-" + df["normalisation"] + df["type"].astype(str)
    ).str.replace("-.", ".", regex=False)

    df["antibody"] = np.where(
        df["assay"] == "ChIP", df["name"].str.split("_").str[-1], ""
    )

    return df


def create_color_palette(names) -> dict:
    """Assign a tab20 colour to each unique sample name."""
    cmap = plt.get_cmap("tab20")
    return {name: mcolors.to_hex(cmap(i % cmap.N)) for i, name in enumerate(names)}


def track_kwargs(track, color: str) -> tuple[str, dict]:
    """Return the PlotNado track alias and options for one input file."""
    common = dict(
        data=str(track.path),
        title=track.track_name,
        color=color,
        label_on_track=True,
        title_location="left",
    )
    if track.type == ".bed":
        return "bed", common

    # Tracks share a y-axis only within an assay (and, for ChIP, an antibody).
    group = "-".join(
        part
        for part in (track.assay, track.antibody, track.method, track.normalisation)
        if part
    )

    return "bigwig", dict(
        common,
        style="scatter" if track.assay == "METH" else "fill",
        autoscale_group=group,
        data_range_style="text",
        scale_location="right",
    )


def _usable(path: str | None) -> bool:
    """The genes rule touches an empty file when no annotation is configured."""
    return bool(path) and Path(path).exists() and Path(path).stat().st_size > 0


def build_figure(df: pd.DataFrame, genes_file: str | None = None) -> pn.GenomicFigure:
    """Build a PlotNado figure containing a scale bar, genes and all data tracks."""
    logger.info("Creating PlotNado figure...")

    fig = pn.GenomicFigure()
    fig.add_track("scalebar")

    if _usable(genes_file):
        fig.add_track("genes", data=str(genes_file), minimum_gene_length=1000)
    elif genes_file:
        logger.warning(f"Genes file {genes_file} is missing or empty; skipping genes track")

    colors = create_color_palette(df["name"].unique())
    for track in df.itertuples():
        alias, kwargs = track_kwargs(track, colors[track.name])
        logger.info("ASSAY: {} TRACK: {} TYPE: {} OPTIONS: {}", track.assay, track.track_name, alias, kwargs)
        fig.add_track(alias, **kwargs)

    return fig


def read_regions(regions: str) -> pd.DataFrame:
    """Read the plotting BED exactly as ``seqnado.outputs`` does, so names match."""
    coords = pd.read_csv(str(regions), sep="\t", header=None, comment="#")
    coords.columns = BED_COLUMNS[: len(coords.columns)]
    return coords


def generate_region_name(region) -> str:
    """Generate a filename-friendly name for a genomic region."""
    if hasattr(region, "Name") and region.Name:
        return region.Name
    return f"{region.Chromosome}-{region.Start}-{region.End}"


def save_plots(fig: pn.GenomicFigure, coords: pd.DataFrame, outdir: Path, plotting_format: str):
    """Render and save one plot per region."""
    logger.info(f"Output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    for region in coords.itertuples():
        fig_name = generate_region_name(region)
        region_coords = f"{region.Chromosome}:{region.Start}-{region.End}"
        output_file = outdir / f"{fig_name}.{plotting_format}"
        logger.info(f"Saving plot for region {fig_name}: {output_file}")
        fig.save(output_file, region_coords)


def main():
    (
        track_assays,
        input_data,
        output_plots,
        regions,
        outdir,
        plotting_format,
        genes,
        template,
    ) = snakemake_setup()
    configure_matplotlib()

    logger.info(
        f"Starting PlotNado {pn.__version__} visualisation for {sorted(set(track_assays))}"
    )
    logger.debug(f"Input files: {input_data}")
    logger.debug(f"Output plots: {output_plots}")
    logger.debug(f"Plotting regions: {regions}")

    df = load_tracks(input_data, track_assays)
    coords = read_regions(regions)
    logger.info(f"Found {len(coords)} regions to plot")

    fig = build_figure(df, genes)
    save_plots(fig, coords, Path(outdir), plotting_format)

    logger.info("Saving PlotNado template...")
    fig.to_toml(template)
    logger.info("PlotNado visualisation complete!")


if __name__ == "__main__":
    main()
