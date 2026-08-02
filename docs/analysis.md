# Python Analysis API

[← Back to main page](index.md)

SeqNado ships a `SeqNadoProject` class for programmatic access to a completed pipeline output directory. It discovers files automatically, supports rich filtering by condition / antibody / sample, and returns `list[Path]` or `pd.DataFrame` — no path guessing required.

## Quick Start

```python
from seqnado.analysis import SeqNadoProject, open_project

# Single-assay run
proj = SeqNadoProject("seqnado_output/")
# or point at a specific config / design file
proj = SeqNadoProject(
    "seqnado_output/",
    config="chip_config.yml",
    design="design.tsv",
)

# Auto-detect single vs multiomics
proj = open_project("seqnado_output/")
```

`config` and `design` are auto-discovered when omitted: SeqNado scans the parent directory for `*_config.yml` and `design.tsv` / `design.csv`.

---

## Multiomics Runs

Multiomics pipelines place each assay under its own subdirectory:

```text
seqnado_output/
├── chip/
├── atac/
└── rna/
```

Use `SeqNadoMultiProject` (or `open_project`, which detects the layout automatically):

```python
from seqnado.analysis import SeqNadoMultiProject, open_project

mp = SeqNadoMultiProject("seqnado_output/")
# or
mp = open_project("seqnado_output/")   # returns SeqNadoMultiProject automatically

mp.assays          # ['atac', 'chip', 'rna']

# Per-assay access (full SeqNadoProject API)
mp.chip                      # SeqNadoProject for chip/
mp["atac"]                   # SeqNadoProject for atac/
mp.rna.load_counts()
mp.chip.bigwigs(scale="csaw")

# Iterate
for name, proj in mp.items():
    print(name, proj.samples)
```

### Cross-assay queries

`bigwigs()`, `peaks()`, `bams()`, and `load_peaks()` on the multi-project query **all assays at once** and return a `pd.DataFrame` with an `assay` column so provenance is always clear:

```python
# All DMSO bigwigs from every assay
mp.bigwigs(condition="DMSO")
#    assay  path
#    chip   seqnado_output/chip/bigwigs/deeptools/unscaled/sampleA_H3K27ac.bigWig
#    rna    seqnado_output/rna/bigwigs/deeptools/unscaled/sampleA.bigWig

# Restrict to specific assays
mp.bigwigs(condition="DMSO", assays=["chip", "cat"])

# Peaks from all assays, non-merged only
mp.peaks(merged=False)

# BAMs for a condition across assays
mp.bams(condition="DMSO")
```

### Cross-assay filtering

`filter()` returns a `_MultiFilteredProject` that pre-applies the mask to every method.  You can narrow by condition, antibody, sample, group, **and** assay:

```python
view = mp.filter(condition="DMSO")
view.assays            # all assays that have DMSO samples
view.bigwigs()         # DMSO bigwigs, assay column present
view.peaks()
view.bams()
view.load_peaks()

# Restrict to specific assays
view = mp.filter(condition="DMSO", assays=["chip", "cat"])
view.bigwigs(scale="unscaled")

# Chainable
mp.filter(condition="DMSO").filter(antibody="H3K27me3").bigwigs()
```

### Combined diagnostics

```python
mp.summary()       # File type × Present table, 'assay' column prepended
mp.what_exists()   # Per-sample availability table across all assays
```

Each per-assay project is a full `SeqNadoProject` with the complete API.  Discovery, loaders, and all single-assay methods work as normal within the assay subdirectory.

---

## Discovery

Inspect what the pipeline produced before fetching files.

```python
proj.samples          # ['sampleA_H3K27ac', 'sampleA_Input', ...]
proj.conditions       # ['control', 'treated']
proj.antibodies       # ['H3K27ac', 'Input']
proj.groups           # scaling groups from design

proj.pileup_methods   # ['bamnado', 'deeptools']
proj.scales           # ['csaw', 'spikein', 'unscaled']
proj.spikein_methods  # ['orlando']
proj.peak_methods     # ['lanceotron', 'macs2']
proj.consensus_groups # groups with merged bigWigs / peaks
proj.log_rules        # rule names with log files present
proj.benchmark_rules  # rule names with benchmark files present
proj.assay            # 'chip'  (from config or report filename)

# samples matching metadata criteria
proj.samples_for(condition="treated")
proj.samples_for(antibody="H3K27ac")
proj.samples_for(condition="treated", antibody="H3K27ac")

# metadata dict for one sample
proj.metadata_for("sampleA_H3K27ac")
# → {'condition': 'treated', 'antibody': 'H3K27ac', 'group': 'grp1'}

# availability table (one row per sample, boolean columns)
proj.what_exists()
#    sample               condition  antibody  bam    bigwig  peaks  ...
#    sampleA_H3K27ac      treated    H3K27ac   True   True    True   ...

# high-level summary
proj.summary()          # File type → ✓ / –
proj.summary(detail=True)  # BigWig rows split by method × scale with counts
```

---

## File Access

All methods return `list[Path]` (existing files only by default).  Pass `only_existing=False` to include paths that may not exist yet.

SeqNado enums can be used interchangeably with strings for `method`, `scale`, and `spikein_method` parameters:

```python
from seqnado.core import PileupMethod, DataScalingTechnique, PeakCallingMethod, SpikeInMethod

proj.bigwigs(method=PileupMethod.DEEPTOOLS, scale=DataScalingTechnique.UNSCALED)
# equivalent to
proj.bigwigs(method="deeptools", scale="unscaled")
```

### BigWigs

```python
proj.bigwigs()                          # all bigWigs
proj.bigwigs(method="deeptools")
proj.bigwigs(scale="csaw")
proj.bigwigs(scale="spikein", spikein_method="orlando")
proj.bigwigs(merged=True)               # consensus group tracks only
proj.bigwigs(merged=False)              # individual sample tracks only
proj.bigwigs(strand="plus")             # RNA-seq stranded tracks
proj.bigwigs(condition="treated", scale="unscaled", merged=False)

# Full metadata index
proj.bigwig_dataframe()  # DataFrame: path, sample, method, scale, spikein_method, merged, strand
```

### Peaks

```python
proj.peaks()
proj.peaks(method="macs2")
proj.peaks(method=PeakCallingMethod.LANCEOTRON)
proj.peaks(merged=True)                 # consensus peaks only
proj.peaks(condition="treated")

proj.peak_dataframe()   # DataFrame: path, sample, method, merged
```

### BAMs

```python
proj.bams()
proj.bams(antibody="H3K27ac")           # IP bams only
proj.bams(condition="treated")
proj.bams(merged=True)                  # merged BAMs from aligned/merged/
```

### Counts

```python
proj.counts()                           # Path to featureCounts matrix (or None)
proj.counts(method="salmon")            # Path to Salmon counts
```

### Assay-specific

```python
proj.methylation()                      # MethylDackel bedGraph files
proj.methylation(genome="hg38")
proj.methylation(inverted=True)         # TAPS-inverted files

proj.vcf()                              # variant VCFs
proj.vcf(annotated=True)

proj.contacts()                         # MCC .mcool files

proj.heatmaps()                         # DeepTools heatmap PDFs
proj.heatmaps(method="deeptools", scale="csaw")
proj.heatmaps(plot_type="metaplot")

proj.track_plots()                      # PlotNado visualisations
proj.track_plots(method="deeptools", scale="unscaled")

proj.normalisation_factors()            # spike-in / CSAW factor TSVs
proj.normalisation_factors(method="orlando")
```

### QC and logs

```python
proj.qc()                               # QC HTML reports
proj.qc(tool="fastqc_raw")
proj.qc(sample="sampleA_H3K27ac")

proj.report()                           # MultiQC HTML (or None)
proj.protocol()                         # protocol.txt (or None)

proj.logs(rule="bowtie2")
proj.logs(rule="bowtie2", sample="sampleA_H3K27ac")
```

---

## Data Loaders

Loaders parse files and return `pd.DataFrame`.

```python
# Gene / guide counts matrix (featureCounts or Salmon)
proj.load_counts()
proj.load_counts(method="salmon")

# All peak BED files concatenated
proj.load_peaks()
proj.load_peaks(method="macs2", condition="treated")
# → columns: chrom, start, end, name, score, strand, sample, peak_method, merged

# Per-sample alignment step statistics
proj.load_alignment_stats()
# → columns: sample, Lost at trim, Lost at align, ..., Retained

# FRiP scores (deeptools plotEnrichment output)
proj.load_frip()
# → columns: sample, peak_method, ...

# Library complexity (Picard MarkDuplicates or samtools markdup)
proj.load_library_complexity()

# Snakemake benchmark runtimes / memory
proj.load_benchmarks()
proj.load_benchmarks(rule="bowtie2")
# → columns: rule, sample, s, h:m:s, max_rss, max_vms, ...

# Spike-in read counts (from bam_split rule)
proj.load_spikein_stats()
# → columns: sample, reference_reads, spikein_reads

# Spike-in / CSAW normalisation factors
proj.load_normalisation_factors()
proj.load_normalisation_factors(method="orlando")

# DESeq2 differential expression results (RNA-seq)
proj.load_deseq2()                              # all contrasts concatenated
proj.load_deseq2(contrast="treated_vs_control") # single contrast
# → columns: contrast, gene, log2FoldChange, padj, ...

# MAGeCK CRISPR screen results
proj.load_mageck()                              # gene-level test results
proj.load_mageck(type="sgrna")
proj.load_mageck(type="gene", analysis="mle")
```

---

## Filtered Views

`filter()` returns a `_FilteredProject` that pre-applies a sample mask to every file-access and loader method.  Views are chainable.

```python
treated = proj.filter(condition="treated")
treated.samples         # only treated samples
treated.conditions      # ['treated']
treated.antibodies      # antibodies present in treated samples
treated.pileup_methods  # methods with data for treated samples
treated.scales
treated.peak_methods

treated.bigwigs(scale="unscaled")
treated.peaks(method="macs2")
treated.bams()
treated.load_peaks()

# Chaining
proj.filter(condition="treated").filter(antibody="H3K27ac").bigwigs()
```

The filter is additive: an explicit `sample=` kwarg on a filtered method call takes precedence.

---

## Utilities

### Annotating file lists with metadata

```python
bws = proj.bigwigs(scale="unscaled", merged=False)
df = proj.enrich(bws)
# → columns: path, sample, condition, antibody, group
```

### Raw file search

```python
proj.select_files("*.bam")             # rglob (recursive by default)
proj.select_files("logs/bowtie2/*.log", recursive=False)
proj.select_files("**/*.bigWig")       # explicit ** uses glob
```

### Per-sample file dict

```python
proj.sample_files("sampleA_H3K27ac")
# → {'bam': Path(...), 'bai': Path(...), 'bigwigs': [...], 'peaks': [...], ...}
```

### BAM index path

```python
proj.bai_for(bam_path)   # raises FileNotFoundError if absent
```

### Condition pairs (for contrasts)

```python
proj.condition_pairs()
# → [('treated', 'control'), ('control', 'treated')]
```

### Reload after pipeline re-run

`SeqNadoProject` caches file indexes on first access.  Call `reload()` to force re-discovery:

```python
proj.reload()
```

---

## Common Workflows

### Load all treated H3K27ac bigWigs with metadata

```python
bws = proj.filter(condition="treated", antibody="H3K27ac").bigwigs(scale="csaw")
df = proj.enrich(bws)
```

### Compare spike-in read fractions

```python
stats = proj.load_spikein_stats()
stats["spikein_fraction"] = stats["spikein_reads"] / (
    stats["reference_reads"] + stats["spikein_reads"]
)
```

### Load DESeq2 results for all contrasts

```python
de = proj.load_deseq2()
sig = de[de["padj"] < 0.05]
sig.groupby("contrast").size()
```

### Identify which samples are missing bigWigs

```python
we = proj.what_exists()
we[we["bam"] & ~we["bigwig"]]
```

### Benchmark profiling

```python
df = proj.load_benchmarks()
df.groupby("rule")[["s", "max_rss"]].mean().sort_values("s", ascending=False)
```
