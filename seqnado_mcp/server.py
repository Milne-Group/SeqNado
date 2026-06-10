"""
SeqNado MCP Server

Exposes SeqNado pipeline operations as MCP tools for use in agentic workflows.
Designed as the middle piece of: GEO MCP → SeqNado MCP → Plotting MCP.
"""

from __future__ import annotations

import asyncio
import json
import os
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any

import yaml
from fastmcp import FastMCP

from seqnado import Assay
from seqnado.config import (
    ATACAssayConfig,
    BigwigConfig,
    BowtieIndex,
    CATAssayConfig,
    ChIPAssayConfig,
    CRISPRAssayConfig,
    GenomeConfig,
    MCCAssayConfig,
    MethylationAssayConfig,
    PeakCallingConfig,
    ProjectConfig,
    QCConfig,
    RNAAssayConfig,
    RNAQuantificationConfig,
    SeqnadoConfig,
    SNPAssayConfig,
    STARIndex,
)
from seqnado.config.user_input import load_genome_configs

mcp = FastMCP("seqnado", instructions=(
    "SeqNado pipeline server. Use list_genomes first to see available reference genomes. "
    "Typical workflow: download_data → create_design → generate_config → launch_pipeline → monitor_pipeline → get_outputs."
))

GENOME_CONFIG_PATH = Path(
    os.getenv("SEQNADO_CONFIG", Path.home())
) / ".config/seqnado/genome_config.json"

ASSAY_CONFIG_CLS = {
    "atac": ATACAssayConfig,
    "chip": ChIPAssayConfig,
    "cat": CATAssayConfig,
    "rna": RNAAssayConfig,
    "snp": SNPAssayConfig,
    "mcc": MCCAssayConfig,
    "meth": MethylationAssayConfig,
    "crispr": CRISPRAssayConfig,
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_raw_genomes() -> dict[str, Any]:
    if not GENOME_CONFIG_PATH.exists():
        raise FileNotFoundError(
            f"Genome config not found at {GENOME_CONFIG_PATH}. Run 'seqnado init' first."
        )
    with open(GENOME_CONFIG_PATH) as f:
        return json.load(f)


def _genome_config(name: str, assay_str: str) -> GenomeConfig:
    assay = Assay(assay_str)
    genomes = load_genome_configs(assay)
    if name not in genomes:
        raise ValueError(f"Genome '{name}' not found. Available: {list(genomes)}")
    return genomes[name]


async def _run(cmd: list[str], cwd: str | None = None) -> tuple[int, str, str]:
    proc = await asyncio.create_subprocess_exec(
        *cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd=cwd,
    )
    stdout, stderr = await proc.communicate()
    return proc.returncode, stdout.decode(), stderr.decode()


# ---------------------------------------------------------------------------
# Tools
# ---------------------------------------------------------------------------

@mcp.tool()
def list_genomes(assay: str = "atac") -> dict:
    """List available reference genomes and their configured paths.

    Args:
        assay: Assay type (atac/chip/cat/rna/snp/mcc/meth/crispr). Affects which
               index type is loaded (STAR for rna, Bowtie2 for all others).
    """
    raw = _load_raw_genomes()
    result = {}
    for name, data in raw.items():
        result[name] = {
            "organism": data.get("organism"),
            "version": data.get("version"),
            "bt2_index": data.get("bt2_index"),
            "star_index": data.get("star_index"),
            "fasta": data.get("fasta"),
            "gtf": data.get("gtf"),
            "blacklist": data.get("blacklist"),
            "chromosome_sizes": data.get("chromosome_sizes"),
        }
    return result


@mcp.tool()
async def download_data(
    metadata_tsv: str,
    outdir: str = "fastqs",
    assay: str | None = None,
    cores: int = 4,
) -> dict:
    """Download sequencing data from GEO/SRA using a metadata TSV file.

    The TSV must have columns: run_accession, sample_title, library_name.
    Optionally: library_layout (PAIRED/SINGLE).

    Args:
        metadata_tsv: Path to TSV metadata file (from GEO MCP or manual).
        outdir: Directory for downloaded FASTQs.
        assay: If provided, also generates a design CSV after download.
        cores: Parallel download jobs.
    """
    cmd = ["seqnado", "download", metadata_tsv, "-o", outdir, "-c", str(cores)]
    if assay:
        cmd += ["-a", assay]

    rc, stdout, stderr = await _run(cmd)
    return {
        "success": rc == 0,
        "returncode": rc,
        "stdout": stdout[-3000:] if len(stdout) > 3000 else stdout,
        "stderr": stderr[-2000:] if len(stderr) > 2000 else stderr,
        "outdir": str(Path(outdir).resolve()),
    }


@mcp.tool()
async def create_design(
    assay: str,
    fastq_dir: str,
    output: str | None = None,
    ip_to_control: str | None = None,
) -> dict:
    """Generate a sample design CSV from FASTQ files.

    Parses FASTQ filenames to extract sample IDs, read pairs, and (for ChIP/CAT)
    IP-to-control pairings.

    Args:
        assay: Assay type (atac/chip/cat/rna/snp/mcc/meth/crispr).
        fastq_dir: Directory containing FASTQ files.
        output: Output CSV path (default: metadata_{assay}.csv).
        ip_to_control: For ChIP/CAT: comma-separated antibody:control pairs,
                       e.g. "H3K4me3:IgG,H3K27ac:IgG".
    """
    fastq_glob = str(Path(fastq_dir) / "*.fastq.gz")
    cmd = ["seqnado", "design", assay, fastq_glob, "--accept-all-defaults"]
    if output:
        cmd += ["-o", output]
    if ip_to_control:
        cmd += ["--ip-to-control", ip_to_control]

    rc, stdout, stderr = await _run(cmd)
    design_path = output or f"metadata_{assay}.csv"
    return {
        "success": rc == 0,
        "returncode": rc,
        "design_file": design_path,
        "stdout": stdout[-3000:] if len(stdout) > 3000 else stdout,
        "stderr": stderr[-2000:] if len(stderr) > 2000 else stderr,
    }


@mcp.tool()
def generate_config(
    assay: str,
    genome: str,
    project_name: str,
    workdir: str = ".",
    make_bigwigs: bool = True,
    call_peaks: bool = True,
    create_heatmaps: bool = False,
    output: str | None = None,
) -> dict:
    """Build a SeqNado config YAML from Pydantic models (no interactive prompts).

    Generates a config_{assay}.yaml with sensible defaults. Pass the returned
    config_path to launch_pipeline.

    Args:
        assay: Assay type (atac/chip/cat/rna/snp/mcc/meth/crispr).
        genome: Genome name (from list_genomes).
        project_name: Name for this project/run.
        workdir: Working directory for the run (default: current dir).
        make_bigwigs: Generate bigwig coverage tracks.
        call_peaks: Run peak calling (ATAC/ChIP/CAT/MCC only).
        create_heatmaps: Generate deepTools heatmaps.
        output: Config output path (default: {workdir}/config_{assay}.yaml).
    """
    assay_str = assay.lower()
    assay_enum = Assay(assay_str)
    genome_cfg = _genome_config(genome, assay_str)

    workdir_path = Path(workdir).resolve()
    today = datetime.today().strftime("%Y-%m-%d")

    project = ProjectConfig(
        name=project_name.replace(" ", "_"),
        date=today,
        directory=workdir_path,
    )
    qc = QCConfig(remove_blacklist=bool(genome_cfg.blacklist))

    bigwig_cfg = BigwigConfig() if make_bigwigs else None
    peak_cfg = PeakCallingConfig() if call_peaks else None

    base = dict(
        genome=genome_cfg,
        bigwigs=bigwig_cfg,
        create_heatmaps=create_heatmaps,
    )

    match assay_enum:
        case Assay.ATAC:
            assay_config = ATACAssayConfig(**base, tn5_shift=True, peak_calling=peak_cfg)
        case Assay.CHIP:
            assay_config = ChIPAssayConfig(**base, peak_calling=peak_cfg)
        case Assay.CAT:
            assay_config = CATAssayConfig(**base, peak_calling=peak_cfg)
        case Assay.RNA:
            assay_config = RNAAssayConfig(**base, rna_quantification=RNAQuantificationConfig())
        case Assay.SNP:
            assay_config = SNPAssayConfig(genome=genome_cfg, create_heatmaps=create_heatmaps)
        case Assay.MCC:
            assay_config = MCCAssayConfig(**base, peak_calling=peak_cfg)
        case Assay.METH:
            assay_config = MethylationAssayConfig(
                genome=genome_cfg, create_heatmaps=create_heatmaps
            )
        case Assay.CRISPR:
            assay_config = CRISPRAssayConfig(genome=genome_cfg, create_heatmaps=create_heatmaps)
        case _:
            raise ValueError(f"Unknown assay: {assay}")

    config = SeqnadoConfig(
        assay=assay_enum,
        project=project,
        genome=genome_cfg,
        metadata=f"metadata_{assay_str}.csv",
        qc=qc,
        assay_config=assay_config,
    )

    config_dict = config.model_dump(mode="json", exclude_none=False)

    out_path = Path(output) if output else workdir_path / f"config_{assay_str}.yaml"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)

    return {
        "config_path": str(out_path),
        "assay": assay_str,
        "genome": genome,
        "metadata_expected": f"metadata_{assay_str}.csv",
    }


@mcp.tool()
async def launch_pipeline(
    assay: str,
    workdir: str = ".",
    preset: str = "le",
    cores: int | None = None,
    dry_run: bool = False,
) -> dict:
    """Launch the SeqNado Snakemake pipeline.

    Run from the directory containing config_{assay}.yaml and metadata_{assay}.csv.

    Args:
        assay: Assay type (atac/chip/cat/rna/snp/mcc/meth/crispr).
        workdir: Directory with config and metadata files.
        preset: Execution preset — le (local), ls (local+Singularity),
                lc (local+Conda), ss (SLURM HPC).
        cores: Override number of cores (uses preset default if None).
        dry_run: Print planned jobs without running.
    """
    cmd = ["seqnado", "pipeline", assay, "--preset", preset]
    if cores:
        cmd += ["--cores", str(cores)]
    if dry_run:
        cmd += ["--dry-run"]

    rc, stdout, stderr = await _run(cmd, cwd=workdir)
    return {
        "success": rc == 0,
        "returncode": rc,
        "workdir": str(Path(workdir).resolve()),
        "stdout": stdout[-5000:] if len(stdout) > 5000 else stdout,
        "stderr": stderr[-3000:] if len(stderr) > 3000 else stderr,
        "log_dir": str(Path(workdir) / "seqnado_output" / "logs"),
    }


@mcp.tool()
def monitor_pipeline(workdir: str = ".", last_n_lines: int = 50) -> dict:
    """Check pipeline progress by scanning Snakemake logs.

    Counts job states (running/done/failed) and returns recent log lines.

    Args:
        workdir: Pipeline working directory.
        last_n_lines: Number of recent log lines to return.
    """
    log_dir = Path(workdir) / "seqnado_output" / "logs"
    snakemake_log_dir = Path(workdir) / ".snakemake" / "log"

    # Find most recent snakemake log
    log_file = None
    for d in [snakemake_log_dir, log_dir]:
        if d.exists():
            logs = sorted(d.glob("*.log"), key=lambda p: p.stat().st_mtime, reverse=True)
            if logs:
                log_file = logs[0]
                break

    if not log_file:
        return {"error": "No Snakemake log found. Has the pipeline been launched?"}

    text = log_file.read_text(errors="replace")
    lines = text.splitlines()

    # Count job states
    finished = sum(1 for l in lines if "Finished job" in l)
    errors = sum(1 for l in lines if "Error in rule" in l or "MissingOutputException" in l)
    running_jobs = [l for l in lines if l.strip().startswith("rule ") and "localrule" not in l]

    return {
        "log_file": str(log_file),
        "jobs_finished": finished,
        "jobs_failed": errors,
        "recent_lines": lines[-last_n_lines:],
        "status": "failed" if errors else ("running" if running_jobs else "done_or_idle"),
    }


@mcp.tool()
def get_outputs(assay: str, workdir: str = ".") -> dict:
    """Return structured paths to pipeline outputs for downstream use (e.g. plotting MCP).

    Scans seqnado_output/ for bigwigs, peaks, BAMs, and QC reports.

    Args:
        assay: Assay type used in the pipeline run.
        workdir: Pipeline working directory.
    """
    out = Path(workdir) / "seqnado_output"
    if not out.exists():
        return {"error": f"Output directory not found: {out}"}

    def glob_str(pattern: str) -> list[str]:
        return sorted(str(p) for p in out.glob(pattern))

    result: dict[str, Any] = {
        "output_dir": str(out.resolve()),
        "bams": glob_str("aligned/*.bam"),
        "bigwigs": {
            "deeptools_unscaled": glob_str("bigwigs/deeptools/unscaled/*.bw"),
            "deeptools_spikein": glob_str("bigwigs/deeptools/spikein/*.bw"),
            "merged": glob_str("bigwigs/deeptools/merged/*.bw"),
            "all": glob_str("bigwigs/**/*.bw"),
        },
        "peaks": {
            "macs2": glob_str("peaks/macs2/*.bed") + glob_str("peaks/macs2/*.narrowPeak"),
            "macs3": glob_str("peaks/macs3/*.bed") + glob_str("peaks/macs3/*.narrowPeak"),
            "seacr": glob_str("peaks/seacr/*.bed"),
            "homer": glob_str("peaks/homer/*.bed"),
            "all": glob_str("peaks/**/*.bed"),
        },
        "qc": {
            "fastqc": glob_str("qc/fastqc_raw/*.html"),
            "report": glob_str(f"seqnado_report_{assay}.html"),
        },
    }

    if assay == "rna":
        result["quantification"] = {
            "featurecounts": glob_str("quantification/featurecounts/*.txt"),
            "salmon": glob_str("quantification/salmon/*/quant.sf"),
            "deseq2": glob_str("quantification/deseq2/*.csv"),
        }

    if assay == "meth":
        result["methylation"] = glob_str("methylation/methyldackel/*.bedGraph")

    if assay == "snp":
        result["vcf"] = glob_str("vcf/*.vcf.gz")

    return result


if __name__ == "__main__":
    mcp.run()
