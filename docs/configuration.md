[← Back to main page](index.md)

# Configuration

Build your analysis configuration after genome setup (see [Genomes](genomes.md) for genome configuration).

## What Configuration Does

The `seqnado config` command generates a YAML configuration file that defines all workflow parameters for your analysis. This file acts as a blueprint for the entire pipeline run:

- **Specifies tools and versions** to use for each step (alignment, peak calling, quantification, etc.)
- **Sets tool-specific parameters** (e.g., peak caller mode, strandedness, normalisation method)
- **Defines output locations** and naming conventions
- **Captures assay-specific choices** (which normalisation method? which peak caller? which pileup method?)
- **Enables reproducibility** by version-locking the configuration and making it version-controlled

The configuration file is generated based on your assay type. You can review, edit, and customize it before running the pipeline. Changes to the config do **not** affect samples already processed — each analysis can use different parameters by pointing to different config files.

## Information to Know Before Configuration

**⚠️ CRITICAL: Before running `seqnado config`, you must first run `seqnado init` to configure your reference genome(s). See [Genomes](genomes.md).**

Once you have configured genomes, you'll be asked about your configuration choices during `seqnado config`. SeqNado automatically detects some parameters from your FASTQ files (sequencing type, read length, adapter content), but you should know these details upfront:

| Information | Why It Matters | Example |
|-------------|----------------|---------|
| **Reference genome & build** | REQUIRED — Configure via `seqnado init` before running config | hg38, mm39, mm10, dm6, or custom build |
| **Strandedness** (RNA-seq only) | SeqNado will ask during config; check with sequencing facility | Unstranded (0), forward (1), or reverse (2) — usually dUTP-based kits are reverse-stranded |
| **Spike-in species** (if applicable) | You'll specify this during config; confirm what you used | *Drosophila* dm6, lambda phage, methylation controls (Lambda, 250bp-v1, 2kb-unmod), etc. |
| **Control/input samples** (ChIP-seq only) | Not asked during config, but required in [Design](design.md) file if using with-input normalisation | Do you have paired input controls? Name them as `samplename_input` in design file |
| **GTF annotation** (RNA-seq only) | You'll need this for quantification; verify it matches your genome build | GRCh38 GTF for hg38; mm10 GTF for mm10 |
| **Peak type expectations** (ChIP/ATAC/CAT) | Helps you choose peak caller during config | Narrow (TFs, H3K4me3) vs broad (H3K27me3, H3K36me3) |

---

## Questions Asked During Configuration

The `seqnado config` command guides you through a series of interactive questions tailored to your assay type. Below are the complete questions for each assay, with defaults shown in parentheses.

### Global Questions (All Assays Except SNP & CRISPR)

These questions appear for ATAC, ChIP, CAT, RNA, MCC, and Methylation:

```
Make Bigwigs? (default: yes)
  ├─ Bigwig method(s) — comma-separated (default: deeptools; options: deeptools, bamnado)
  ├─ Binsize for bigwigs (default: 10)
  └─ Bigwig scaling method(s) — comma-separated (default: unscaled; options: unscaled, csaw)

Perform plotting? (default: no)
  ├─ Path to coordinates BED file (optional)
  └─ Path to genes BED file (optional)

Make heatmaps? (default: no)

Generate GEO submission files? (default: no)

Make UCSC hub? (default: no)
  ├─ UCSC hub directory (default: seqnado_output/hub/)
  ├─ Email address (default: $USER@example.com)
  ├─ Genome name for hub (default: genome name from config)
  └─ Color by field (default: samplename)
```

### ATAC-seq Specific

```
Shift ATAC reads? (default: yes)

Call peaks? (default: yes)
  ├─ Peak calling method(s) — comma-separated (default: lanceotron; options: macs, seacr, lanceotron)
  ├─ Generate consensus counts from Design merge column? (default: no)
  └─ Run motif analysis on called peaks? (default: no)
      └─ Motif analysis method(s) — if yes (default: homer; options: homer, meme)

Make dataset for ML? (default: no)
  └─ Use regions BED file? (default: yes)
      └─ Path to regions BED file (default: path/to/regions.bed)
         OR Binsize for dataset (default: 1000)
```

### ChIP-seq Specific

```
Call peaks? (default: yes)
  ├─ Peak calling method(s) — comma-separated (default: lanceotron; options: macs, seacr, lanceotron)
  ├─ Generate consensus counts from Design merge column? (default: no)
  └─ Run motif analysis on called peaks? (default: no)
      └─ Motif analysis method(s) — if yes (default: homer)

Do you have spike-in? (default: no)
  ├─ Normalisation method(s) — comma-separated (default: orlando; options: orlando, with-input, deseq2, edger, csaw)
  ├─ Reference genome (default: hg38)
  ├─ Spike-in genome (default: dm6)
  └─ Spike-in control gene names — if using deseq2 or edger (default: AmpR,Cas9_3p,Cas9_5p)

Make dataset for ML? (default: no)
  └─ Use regions BED file? (default: yes)
      └─ Path to regions BED file (default: path/to/regions.bed)
         OR Binsize for dataset (default: 1000)
```

### CUT&Tag (CAT) Specific

```
Shift CAT reads? (default: no)

Call peaks? (default: yes)
  ├─ Peak calling method(s) — comma-separated (default: seacr; options: macs, seacr, lanceotron)
  ├─ Generate consensus counts from Design merge column? (default: no)
  └─ Run motif analysis on called peaks? (default: no)
      └─ Motif analysis method(s) — if yes (default: homer)

Do you have spike-in? (default: no)
  ├─ Normalisation method(s) — comma-separated (default: orlando)
  ├─ Reference genome (default: hg38)
  ├─ Spike-in genome (default: dm6)
  └─ Spike-in control gene names — if using deseq2 or edger (default: AmpR,Cas9_3p,Cas9_5p)

Make dataset for ML? (default: no)
  └─ Use regions BED file? (default: yes)
      └─ Path to regions BED file
         OR Binsize for dataset
```

### RNA-seq Specific

```
Do you have spike-in? (default: no)
  ├─ Normalisation method(s) — comma-separated (default: deseq2; options: orlando, with-input, deseq2, edger, csaw)
  ├─ Reference genome (default: hg38)
  ├─ Spike-in genome (default: spikein_rna)
  └─ Spike-in control gene names — if using deseq2 or edger (default: AmpR,Cas9_3p,Cas9_5p)

Quantification method? (default: feature_counts; options: feature_counts, salmon)
  └─ Strandedness? — if feature_counts (default: 0; options: 0=unstranded, 1=forward, 2=reverse)
  
Salmon index path? — if salmon method (default: path/to/salmon_index)

Run DESeq2? (default: no)
```

### MCC (Capture-C) Specific

```
Call peaks? (default: yes)
  ├─ Peak calling method(s) — comma-separated (default: lanceotronmcc; options: macs, seacr, lanceotron, lanceotronmcc)
  ├─ Generate consensus counts from Design merge column? (default: no)
  └─ Run motif analysis on called peaks? (default: no)
      └─ Motif analysis method(s) — if yes (default: homer)

Path to viewpoints file (required; default: path/to/viewpoints.bed)

Resolutions for MCC cooler files — comma-separated (default: 100,1000)
```

### Methylation Specific

```
Call methylation? (default: no)
  └─ Spike-in genomes — comma-separated, if yes (default: Lambda,250bp-v1,2kb-unmod)

Methylation assay? (default: taps; options: taps, bsseq)
```

### SNP Calling Specific

```
Call SNPs? (default: yes)
  ├─ SNP calling method? (default: bcftools; options: bcftools)
  └─ Annotate SNPs? (default: no)
      └─ Path to SNP database — if yes (default: path/to/snp_database)

Generate GEO submission files? (default: no)
```

### CRISPR Specific

```
Use MAGeCK for guide RNA analysis? (default: no)

Generate GEO submission files? (default: no)
```

---

### Understanding Multi-Select Questions

Some questions allow **comma-separated selections** for multiple values:
- **Peak calling method(s)**: Can use MACS, SEACR, and Lanceotron simultaneously
- **Normalisation method(s)**: Can run Orlando and DESeq2 on the same data
- **Motif analysis method(s)**: Can run both Homer and MEME on peaks
- **Bigwig method(s)**: Can generate with both deeptools and bamnado

Example:
```
Peak calling method(s) (comma-separated for multiple): macs, seacr, lanceotron
```

The pipeline will run all specified methods and generate outputs for each.

---

### Resource Constraints to Consider

Before finalizing your configuration, review your compute environment:

- **Memory per sample**: Does your cluster have enough RAM for large BAM files (deeptools, Salmon, variant calling)?
- **Runtime limits**: Do you have wall-clock time limits? Some tools (STAR, MACS2, DESeq2) can be slow on large files.
- **Storage**: Will intermediate files (unsorted BAMs, uncompressed bigwigs) fit on your filesystem?
- **Parallelization**: How many samples can you process in parallel without overwhelming resources?

See the [Troubleshooting guide](troubleshooting.md#resource-allocation) for guidance on resource requirements per tool.

---

## Running `seqnado config`

The `seqnado config` command interactively builds a YAML configuration file for your selected assay. It asks questions appropriate to your assay type and generates a configuration that you can review and customize before running the pipeline.

### Assay Types

| Assay | CLI name | When to Use |
|-------|----------|------------|
| ATAC-seq | `atac` | Open chromatin profiling (Tn5-based) |
| ChIP-seq | `chip` | Chromatin immunoprecipitation |
| CRISPR analysis | `crispr` | CRISPR screen analysis |
| CUT&Tag | `cat` | Chromatin profiling by tagmentation (low background) |
| MCC | `mcc` | Capture-C or similar 3C-derived methods |
| Methylation | `meth` | Bisulfite or TAPS methylation calling |
| RNA-seq | `rna` | Quantification from RNA-seq reads |
| SNP analysis | `snp` | Variant calling and annotation |

For all available arguments and flags, see: [seqnado config](cli.md#cli-seqnado-config).

### Example Usage

#### Build a Configuration for ChIP-seq
```bash
seqnado config chip --output chip_config.yaml
```

#### Interactive Multiomics Configuration (Multiple Assays)
```bash
seqnado config --make-dirs --interactive
```

## Third Party Tools

SeqNado integrates with specialized tools for each analysis step. Before running the pipeline, review the guidance below for tools critical to your assay. **Default parameters are sensible for most experiments**, but customization may be necessary for specific data, assays, and experimental designs.

### Supported Tools Reference

Below is a comprehensive list of tools integrated into SeqNado, organized by function. **⚠️ Key configuration warnings are included for critical tools.**

#### Alignment & Indexing

- **bowtie2** — Maps sequencing reads to a reference genome (ATAC, ChIP, CUT&Tag, SNP)
- **STAR** — Splice-aware aligner for RNA-seq data (RNA)
- **salmon** — Pseudoalignment-based RNA quantification (RNA; optional alternative to alignment)
  - ⚠️ **Requires compatible reference index**: GTF version must match your genome build (e.g., GRCh38 GTF for hg38)
  - Faster than alignment-based methods but less flexible for debugging strandedness
  - Generate a Salmon index before running: `salmon index -t transcripts.fa -i salmon_index_hg38`

#### Read Processing

- **samtools** — BAM/SAM file manipulation and sorting
- **picard** — Java tools for high-throughput sequencing data manipulation
- **cutadapt** — Removes adapter sequences from reads
- **trim-galore** — Wrapper for Cutadapt and FastQC for quality control

#### Peak Calling (ChIP-seq, ATAC-seq, CUT&Tag, MCC)

- **MACS2** — Model-based peak calling
  - ⚠️ **Peak calling mode depends on your mark**:
    - **Narrow mode** (default) for transcription factors and sharp marks (H3K4me3, H3K27ac)
    - **Broad mode** for diffuse marks (H3K27me3, H3K36me3). Without `--broad`, you'll miss 90% of broad domains
    - For weak ChIP-seq with few peaks, consider switching to SEACR or Lanceotron instead

- **SEACR** — Peak caller optimized for low-background data (CUT&Tag)
  - Configure stringency based on expected noise:
    - `stringent` mode (recommended for CUT&Tag) reduces false positives
    - `relaxed` mode if you're missing expected peaks
    - Threshold values control sensitivity vs. specificity trade-off

- **Lanceotron** — Deep learning–based peak caller for broad and narrow peaks
  - No hyper-parameter tuning needed; model is pre-trained
  - Can be slow on large files; consider testing on a sample first

- **lanceotronmcc** — Specialized Lanceotron variant for MCC interaction peaks

#### Quantification (RNA-seq)

- **featureCounts** (subread package) — Assigns aligned reads to genomic features for RNA-seq quantification
  - ⚠️ **Critical GTF and strandedness configuration**:
    - Default parameters count reads at the **exon level** (`-t exon`) and aggregate by **gene_id**
    - **Must verify GTF attribute match**: Does your GTF use `gene_id` or `gene_name`? featureCounts requires exact match or all counts will be zero
    - **Strandedness must match library prep**: Use `0` (unstranded), `1` (forward), or `2` (reverse). Using wrong strandedness results in 50–90% fewer counts
    - Paired-end mode (`-p --countReadPairs`) is auto-detected from FASTQ files
    - If all genes get zero counts, check GTF feature attribute match first

- **Salmon** — Faster pseudoalignment-based quantification (requires compatible index; see Alignment & Indexing above)

#### Bigwig Generation & Visualization

- **deeptools bamCoverage** — Generates BigWig files from BAM alignments
  - **Using deeptools' native `--normalizeUsing` ** — To use deeptools' built-in normalization methods (`RPKM`, `CPM`, `BPM`, `RPGC`), select `unscaled` for your bigwig scaling method, then manually edit `third_party_tools.deeptools.bam_coverage.command_line_arguments` to add `--normalizeUsing` (e.g., `"--binSize 10 --normalizeUsing RPKM"`). This disables SeqNado's scaling entirely. Example use case: comparing standard RPKM-normalized bigwigs across studies
  - **Bigwig scaling method** — You choose during config (options: `unscaled`, `csaw`):
    - **`unscaled`**: Raw read coverage without normalization; useful for visual inspection and as input for downstream analysis tools
    - **`csaw`** (ChIP-Seq Analysis with Windows): Applies scaling factors between samples to equalize library depth while preserving broad features; recommended when comparing ChIP-seq samples or using spike-in controls; see [Normalisation Methods](normalisation.md)
  - **Spike-in normalisation for bigwigs** — If you selected spike-in normalisation (Orlando, with-input, DESeq2, or edgeR) during config, scale factors are automatically calculated and applied to bigwig generation, producing spike-in–normalized bigwigs for downstream analysis; see [Normalisation Methods](normalisation.md)
  
!!! info
    **How SeqNado prevents conflicts between `--normalizeUsing` and `--scaleFactor`**: The default deeptools config includes `--normalizeUsing RPKM`. However, when you select spike-in or csaw normalization, SeqNado **automatically removes** `--normalizeUsing` and applies **only** `--scaleFactor`. This prevents deeptools from silently ignoring your scale factors (since `--normalizeUsing` overrides `--scaleFactor` in deeptools). Result: your spike-in/csaw scale factors apply correctly without unwanted RPKM normalization
  
  
- **bamnado** — Alternative tool for BigWig generation and BAM manipulation

- **deeptools heatmap** — Creates heatmaps from BigWig files around genomic coordinates

- **deeptools metaplot** — Generates metaplots of signal around features

#### Motif Analysis (peak regions)

- **Homer** — Motif discovery and annotation in peak regions
- **MEME** — Alternative motif discovery tool

#### Variant Calling (SNP)

- **bcftools** — SNP calling and VCF manipulation
- **SnpEff/SnpSift** — SNP annotation and variant functional impact prediction

#### Methylation & Specialized

- **methyldackel** — Extract methylation calls from bisulfite-seq data
- **MAGeCK** — Statistical analysis of CRISPR screen data
- **UCSC utilities** — Convert and manage genome tracks in UCSC format

#### Data Analysis & Normalisation

- **DESeq2** (R/Bioconductor) — Differential expression analysis (RNA-seq)

- **edgeR** (R/Bioconductor) — Differential expression with spike-in normalisation
  - ⚠️ **Spike-in normalisation setup**:
    - Requires correct spike-in genome specification (e.g., `dm6` for *Drosophila*)
    - Control gene names must match spike-in reference GTF (e.g., `ERCC`)
    - With-input normalisation requires ChIP and input samples paired in design file — see [Design guide](design.md)

- **Orlando** — Spike-in normalisation method for chromatin-based assays (ChIP, CAT)
  - Recommended for spike-in data with balanced endogenous and exogenous reads

- **with-input** — Normalisation using paired input controls (ChIP-seq only)
  - Requires input samples correctly paired in design file
  - Alternative to spike-in normalisation when spike-ins unavailable

- **CSAW** — Cyclic shift aware normalisation; produces merged bigwigs
  - Recommended for spike-in data with unbalanced read distributions
  - See [Normalisation Methods](normalisation.md) for detailed comparison

### Configuration Logic & Tool Interactions

See the [Pipeline Overview](pipeline.md#supported-assays) for guidance on which tools to use for each assay type.

For detailed normalisation methods comparison, see [Normalisation Methods](normalisation.md).

For configuration command options and usage patterns, see [seqnado config](cli.md#cli-seqnado-config).

---

## Common Pitfalls and Best Practices

### Pitfalls to Avoid

**Assuming defaults are correct without validation** 
   - **Problem**: Accepting default peak caller (lanceotron), GTF attributes, or strandedness without checking your specific data
   - **Fix**: For each new assay/GTF/strandedness combo, test peak calling or quantification on 1–2 samples first; spot-check BAM files and counts

**Mismatched GTF and genome builds**
   - **Problem**: Using GRCh37 GTF with hg38 genome, or mm10 GTF with mm39
   - **Fix**: Verify your annotation version matches your genome build at download time; document version numbers in your README

**Wrong peak calling mode for your biological mark**
   - **Problem**: Using MACS2 narrow mode for H3K27me3 (broad mark) → misses 90% of signal; using broad mode for H3K4me3 (sharp mark) → thousands of false positives
   - **Fix**: Know your mark. H3K4me3, H3K27ac, TF binding → narrow. H3K27me3, H3K36me3, H3K9me3 → broad. Test on a rep.

**Forgetting to pair ChIP and input samples in design file**
   - **Problem**: Configuring "with-input" normalisation but not ensuring control sample pairs in [Design](design.md)
   - **Fix**: Before running pipeline, review design file to confirm ChIP→input pairing syntax matched for all samples

**Ignoring resource constraints**
   - **Problem**: Setting memory-intensive parameters (STAR `--limitSjdbInsertNsj`, deeptools multi-threaded) on systems without available RAM
   - **Fix**: Test on 1 sample and monitor resource usage (`top`, cluster monitoring); adjust thread count or binsize if needed

### Best Practices

✅ **Test before large-scale runs**
- Run `seqnado snakemake ... -n` (dry-run) to visualize the workflow DAG before execution
- Test peak calling or quantification on 1–2 samples to validate your tool selections
- Spot-check intermediate outputs: BAM file headers (correct chromosome names?), peak files (reasonable count?), count matrices (non-zero?), BigWigs (visible in IGV?)

✅ **Version-control your configuration files**
- Commit config YAML files to git so analysis is fully reproducible
- Create assay-specific configs (e.g., `config_chip_h3k27ac.yaml`, `config_rna_pe150.yaml`) so you can reuse templates
- Different analyses can and should use different configs — this is your version control

✅ **Validate reference files upfront**
- Before running the full pipeline, verify:
  - GTF file: Spot-check a few known genes and their feature definitions
  - Spike-in genome: Confirm the spike-in species (dm6? ecoli?) is included in your genome build
  - Salmon index: If using Salmon, ensure the index matches your GTF and genome version
  - bedtools/peaks BED files: Verify chromosome names match your genome (e.g., `chr1` vs `1`)

✅ **Monitor the first few samples carefully**
- Check the first sample's outputs in detail before running the full cohort
- Verify BAM file alignment rate (expect >90% for most data)
- Check peak files are reasonable (not zero peaks, not millions of false positives)
- Confirm BigWig tracks look correct in IGV (visible signal? expected chromosome coverage?)
- This catches configuration errors early and saves days of compute time

✅ **Use meaningful output names**
- Use the `--output` flag to give your config file a descriptive name: `seqnado config rna --output lps_response_config.yaml`
- This helps keep track of different analyses without confusion

---

**See Also:**

- [Design Guide](design.md) - Create experimental design files
- [Tools Reference](tools.md) - Configure tool-specific options
- [Troubleshooting](troubleshooting.md#configuration-seqnado-config) - Configuration issues