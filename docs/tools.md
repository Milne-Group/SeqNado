# Third-Party Tools Reference

SeqNado integrates multiple best-in-class bioinformatics tools to provide comprehensive genomics analysis pipelines. This page documents the tools used, their purposes, and key references.

## Tool Versions & Updates

All tools are version-locked in SeqNado containers to ensure reproducibility. To check versions and get tool information, use the `seqnado tools` command:

```bash
# List all tools with versions
seqnado tools

# List tools in a specific category (e.g., Download, Alignment, Analysis)
seqnado tools --category

# View detailed information about a specific tool
seqnado tools macs2

# Show tool help/options from the container
seqnado tools macs2 --options
```

See the [CLI Reference](cli.md#cli-seqnado-tools) for complete documentation of the `tools` command.

## Citation Guidelines

When publishing results from SeqNado, please cite:

1. **SeqNado** itself (citation pending - check [GitHub releases](https://github.com/Milne-Group/SeqNado/releases) for latest version)
2. **Key tools** used in your specific analysis (see [References](#references) section below)
3. **Snakemake** workflow manager: Mölder F, Jablonski KP, Letcher B, et al. Sustainable data analysis with Snakemake. *F1000Research*. 2021;10:33.
4. **Reference genomes** used (e.g., hg38, mm10)

### Example Acknowledgment

*"Data analysis was performed using SeqNado v1.0.2 (Chahrour, C and Smith, AL, 2024), which incorporates Bowtie2 (Langmead & Salzberg, 2012) for alignment, MACS2 (Zhang et al., 2008) for peak calling, and deepTools (Ramírez et al., 2016) for coverage track generation. Workflows were managed with Snakemake (Mölder et al., 2021) within Apptainer containers."*

## Tools

Tools are organized by category, matching the structure of the `seqnado tools` CLI command. For more information about any tool, use `seqnado tools <toolname>`.

<!-- AUTO-GENERATED TOOL SECTIONS - DO NOT EDIT MANUALLY -->
<!-- This section is automatically updated by docs/scripts/generate_tool_citations.py -->
<!-- To update, run: python docs/scripts/generate_tool_citations.py --update -->

