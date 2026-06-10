# SeqNado MCP Server

MCP server wrapping SeqNado pipeline operations. Middle piece of:

```
GEO MCP → SeqNado MCP → Plotting MCP
```

## Tools

| Tool | Purpose |
|------|---------|
| `list_genomes` | Show available reference genomes from `~/.config/seqnado/genome_config.json` |
| `download_data` | Download FASTQs from GEO/SRA via metadata TSV |
| `create_design` | Generate sample design CSV from FASTQ filenames |
| `generate_config` | Build `config_{assay}.yaml` from Pydantic models (no interactive prompts) |
| `launch_pipeline` | Run `seqnado pipeline <assay> --preset <preset>` |
| `monitor_pipeline` | Tail Snakemake logs, count finished/failed jobs |
| `get_outputs` | Return structured output paths (bigwigs/peaks/BAMs) for plotting MCP |

## Setup

```bash
# Install fastmcp into the seqnado conda env
conda activate seqnado
pip install fastmcp
```

## Add to Claude Code

Add to `.claude/settings.json` (or `~/.claude/settings.json`):

```json
{
  "mcpServers": {
    "seqnado": {
      "command": "conda",
      "args": ["run", "-n", "seqnado", "--no-capture-output",
               "python", "-m", "seqnado_mcp.server"],
      "cwd": "/Users/asmith/Documents/2025-oxford-postdoc/software/SeqNado"
    }
  }
}
```

## Typical workflow

```
list_genomes()
  → pick genome name

download_data(metadata_tsv="geo_metadata.tsv", assay="chip")
  → FASTQs in fastqs/

create_design(assay="chip", fastq_dir="fastqs/", ip_to_control="H3K4me3:IgG")
  → metadata_chip.csv

generate_config(assay="chip", genome="hg38", project_name="my_project")
  → config_chip.yaml

launch_pipeline(assay="chip", preset="le")
  → pipeline running

monitor_pipeline()
  → progress

get_outputs(assay="chip")
  → bigwigs/peaks/BAMs for plotting MCP
```
