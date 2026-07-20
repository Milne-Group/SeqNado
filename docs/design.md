[← Back to main page](index.md)

# Design Guide

The `seqnado design` command generates a design CSV file from FASTQ files for a specific assay. If no assay is provided, the tool operates in multiomics mode.

For full arguments and flags, see the CLI reference: [seqnado design](cli.md#cli-seqnado-design).

### FASTQ Files

After generating the configuration and project directory using `seqnado config`, you need to link your FASTQ files into the `fastqs` directory. This ensures that the pipeline can locate and process your input data.

#### Symlinking FASTQ Files

Use the following command to create symbolic links for your FASTQ files:

```bash
ln -s /path/to/your/fastq/files/* <project_directory>/fastqs/
```

Replace `/path/to/your/fastq/files/` with the directory containing your FASTQ files and `<project_directory>` with the path to the project directory created by `seqnado config`.

#### Example

If your FASTQ files are located in `/data/fastq/` and your project directory is `rna_project`, run:

```bash
ln -s /data/fastq/* rna_project/fastqs/
```

This will create symbolic links to all FASTQ files in the `fastqs` directory of your project. The glob pattern `*` expands to all files in the source directory.

#### Safe Naming Strategies for FASTQ Files (Critical)

Before linking your FASTQ files, ensure they follow a consistent naming convention. The `seqnado design` command parses filenames to infer sample metadata (replicates, antibodies, controls, groups), so proper naming is essential for successful pipeline execution.

Below are recommended naming strategies for each assay type:


- **ATAC-seq**:
  ```
  sample-name-rep1_R1.fastq.gz
  sample-name-rep1_R2.fastq.gz
  ```

- **ChIP-seq**:
  ```
  sample-name-rep1_Antibody_R1.fastq.gz
  sample-name-rep1_Antibody_R2.fastq.gz
  sample-name-rep2_Input_R1.fastq.gz
  sample-name-rep2_Input_R2.fastq.gz
  ```
  - `Antibody`: Name of the antibody used for ChIP.
  - `Input`: Control sample.

- **RNA-seq**:
  ```
  sample-name-rep1_R1.fastq.gz
  sample-name-rep1_R2.fastq.gz
  sample-name-rep2_R1.fastq.gz
  sample-name-rep2_R2.fastq.gz
  ```
  - `sample-name`: Unique identifier for the sample.
  - `rep1`, `rep2`: Biological or technical replicate number.
  - `R1`, `R2`: Read pair (forward and reverse).

Using these naming conventions ensures that the pipeline can correctly parse and process your data.

### Metadata Columns Reference

The design CSV can contain the following columns. Most are filled in automatically by `seqnado design`; you only need to edit them by hand for advanced cases (custom control pairing, batch-specific scaling, manual DESeq2 contrasts, etc.).

| Column | Required? | Applies to | Purpose |
| --- | --- | --- | --- |
| `sample_id` | Yes | All | Unique sample identifier |
| `assay` | Auto-filled | All | Assay type (`ATAC`, `ChIP`, `RNA`, …) |
| `r1` / `r2` | Yes / optional | All | Paths to the (paired) raw FASTQ files |
| `ip` | Optional | ChIP, CUT&Tag | Antibody/target name |
| `control` | Optional | ChIP, CUT&Tag | Label of the paired control sample |
| `r1_control` / `r2_control` | Optional | ChIP, CUT&Tag | FASTQ paths for the control sample |
| `scaling_group` | Optional (defaults to `default`) | All | Groups samples for normalisation-factor calculation |
| `consensus_group` | Optional | All | Groups samples for consensus peak calling / counting |
| `condition` | Optional | ATAC, ChIP, CUT&Tag, RNA | Groups samples for bigwig aggregation/subtraction comparisons |
| `group` / `deseq2` | Optional | RNA | Experimental group name and binary DESeq2 reference/treatment encoding |

**`sample_id`** — Unique per sample (or per sample+`ip` combination for IP-seq assays). It's embedded directly into output file and directory names, so it may only contain letters, numbers, underscores, and hyphens (`^[a-zA-Z0-9_-]+$`) — no spaces or other special characters. This is checked by the pipeline before running, so fix naming here rather than downstream.

**`assay`** — Normally set automatically to whatever assay you passed to `seqnado design <assay>`. You'd only touch this by hand when hand-assembling a multiomics design CSV that spans several assay types in one file.

**`r1` / `r2`** — Paths to the symlinked FASTQ files. Leave `r2` blank for single-end data.

**`ip`** — Only relevant for ChIP-seq/CUT&Tag: the antibody or target name (e.g. `H3K27ac`, `Menin`). Leave blank for ATAC/RNA/other assays that don't have an IP step. Setting it lets the same `sample_id` appear multiple times in the design (once per antibody), since uniqueness is enforced on `sample_id` + `ip` together rather than `sample_id` alone.

**`control` / `r1_control` / `r2_control`** — The paired input/IgG/mock control for a ChIP or CUT&Tag sample: `control` is the label used to build the control's output filename, `r1_control`/`r2_control` are the actual control FASTQ paths. Leave all three blank when there's no matched control available — some peak callers (e.g. LanceOtron) can run IP-only, but background subtraction quality is generally better with a matched control, so only skip it when you genuinely don't have one. If you hand-edit the CSV, keep `control` and `r1_control`/`r2_control` in sync — setting the FASTQ paths without also setting the `control` label isn't currently caught by validation and will produce a broken/`None` control filename downstream.

**`scaling_group`** — Defaults to `default` for every sample, meaning all samples are normalised together. Only change this when you have distinct batches whose scaling factors must be computed separately — for example, samples spiked in with different reference genomes, or samples from different sequencing runs that shouldn't be normalised against each other. Note this only affects how normalisation *factors* are calculated, not which samples get merged.

**`consensus_group`** — Blank by default, meaning each sample gets its own independent peak calls/counts, and no consensus set is generated. Set it when you want a shared consensus peak call / count set generated across the group, on top of each sample's individual outputs. This is what you want for comparing samples/conditions against a common set of regions (e.g. differential binding/accessibility analysis) — quantifying each sample against its own private peak set makes cross-sample counts hard to compare, whereas a consensus set gives every sample in the group the same regions to be counted against. Group samples that should share one consensus region set (e.g. all samples for a given IP/mark, or all samples you intend to compare downstream). For Micro-Capture-C (`mcc`), this defaults to `default` automatically since MCC analysis is typically run per merged group; override it only if you specifically need per-replicate consensus calls kept separate.

Rather than hand-editing this column, use `seqnado design`'s `--consensus-by` option to populate it automatically from an existing column or a regex. For example, `--consensus-by ip` groups ChIP-seq/CUT&Tag samples by antibody so each antibody gets its own consensus peak set:

```bash
seqnado design chip fastqs/* --consensus-by ip
```

You can also pass a regex (matched against `sample_id` — or `sample_id`+`ip` for ChIP/CUT&Tag) to derive the grouping from sample names, e.g. `--consensus-by "^([^-]+)"` to group by everything before the first hyphen.

**`condition`** — Blank by default and only acted on when `perform_comparisons: true` is set in the project config, and only once at least 2 unique values exist. It drives bigwig-level comparisons (aggregated mean tracks per condition, plus pairwise subtraction tracks) without merging or otherwise touching BAMs, peak calls, or counts — safe to add even if you're unsure you'll use it. Don't reach for this as a substitute for RNA-seq differential expression grouping — that's what `group`/`deseq2` are for, not `condition`.

Use `--condition-by` to populate it automatically the same way, e.g. extracting `control`/`treated` from sample names:

```bash
seqnado design atac fastqs/* --condition-by "-(control|treated)-"
```

**`group` / `deseq2`** — RNA-seq only, used to drive DESeq2's design formula (relevant when spike-in normalisation is configured). `group` holds the human-readable group name (e.g. `control`, `treated`); `deseq2` is a strictly binary encoding (`0` = reference/control group, `1` = treatment group) fed directly into the DESeq2 model — it is not a general-purpose numeric code, so don't try to represent three or more groups as `0`/`1`/`2`. For experiments with more than two groups, leave `deseq2` blank (the design tool will do this automatically and warn you) and set up contrasts manually — see [Multi-Group Comparisons](#rna-seq-grouping-for-deseq2) below.

### Example Usage

#### Generate a Design CSV for ATAC-seq
```bash
seqnado design atac fastqs/*
```

This command reads all FASTQ files in the `fastqs/` directory (the `*` glob pattern expands to every file) and generates a design CSV file named `metadata.csv` in your project directory.

#### Generate a Design CSV for ChIP-seq with explicit control pairing

For ChIP-seq experiments, the design CSV requires both IP FASTQ files and optionally control FASTQ files. The tool can infer these relationships based on file naming conventions.

#### Simple Case

For simple cases with a single control or when no control is needed, for example:

* SAMPLE1_H3K27ac_R1.fastq.gz
* SAMPLE1_H3K27ac_R2.fastq.gz
* SAMPLE1_Menin_R1.fastq.gz
* SAMPLE1_Menin_R2.fastq.gz
* SAMPLE1_input_R1.fastq.gz
* SAMPLE1_input_R2.fastq.gz
* SAMPLE_2_H3K27ac_R1.fastq.gz
* SAMPLE_2_H3K27ac_R2.fastq.gz

The command would be:

```bash
seqnado design chip fastqs/* 
```

The control will either be left blank if no appropriate files are in the directory or a single control sharing the same sample ID will be broadcast to all IP samples sharing that sample ID. e.g.:

| assay | sample_id | ip      | control | r1                           | r2                           | r1_control                | r2_control                | scaling_group |
|-------|-----------|---------|---------|------------------------------|------------------------------|---------------------------|---------------------------|---------------|
| ChIP  | SAMPLE1   | H3K27ac | input   | SAMPLE1_H3K27ac_R1.fastq.gz  | SAMPLE1_H3K27ac_R2.fastq.gz  | SAMPLE1_input_R1.fastq.gz | SAMPLE1_input_R2.fastq.gz | default       |
| ChIP  | SAMPLE1   | Menin   | input   | SAMPLE1_Menin_R1.fastq.gz    | SAMPLE1_Menin_R2.fastq.gz    | SAMPLE1_input_R1.fastq.gz | SAMPLE1_input_R2.fastq.gz | default       |
| ChIP  | SAMPLE_2  | H3K27ac |         | SAMPLE_2_H3K27ac_R1.fastq.gz | SAMPLE_2_H3K27ac_R2.fastq.gz |                           |                           | default       |



#### Complex Case with Multiple Controls and Ambiguity in Pairing

If there are multiple controls, specify which control corresponds to each IP using the `--ip-to-control` option. For example:

We want the single fixed control `sf-input` to be used for the `H3K27ac` IP, and the double fixed `df-input` control to be used for the `Menin` IP. The FASTQ files are as follows:

  * SAMPLE1_H3K27ac_R1.fastq.gz
  * SAMPLE1_H3K27ac_R2.fastq.gz 
  * SAMPLE1_sf-input_R1.fastq.gz
  * SAMPLE1_sf-input_R2.fastq.gz
  * SAMPLE1_Menin_R1.fastq.gz
  * SAMPLE1_Menin_R2.fastq.gz
  * SAMPLE1_df-input_R1.fastq.gz  
  * SAMPLE1_df-input_R2.fastq.gz

The command would be:

```bash
seqnado design chip fastqs/* --ip-to-control "H3K27ac:sf-input,Menin:df-input"
```

This will generate a design CSV file with the appropriate control pairings. e.g.,

| assay | sample_id | ip      | control  | r1                          | r2                          | r1_control                   | r2_control                   | scaling_group |
|-------|-----------|---------|----------|-----------------------------|-----------------------------|------------------------------|------------------------------|---------------|
| ChIP  | SAMPLE1   | H3K27ac | sf-input | SAMPLE1_H3K27ac_R1.fastq.gz | SAMPLE1_H3K27ac_R2.fastq.gz | SAMPLE1_sf-input_R1.fastq.gz | SAMPLE1_sf-input_R2.fastq.gz | default       |
| ChIP  | SAMPLE1   | Menin   | df-input | SAMPLE1_Menin_R1.fastq.gz   | SAMPLE1_Menin_R2.fastq.gz   | SAMPLE1_df-input_R1.fastq.gz | SAMPLE1_df-input_R2.fastq.gz | default       |

#### Condition-Based Bigwig Comparisons (All Assays)

For all assay types that support bigwigs (ATAC, ChIP, CUT&Tag, RNA), you can optionally add a `condition` column to your design file. When `perform_comparisons: true` is enabled in the configuration, the pipeline will automatically generate:

- **Aggregated condition bigwigs**: Mean signal tracks for each condition group
- **Subtraction bigwigs**: All pairwise condition comparisons (condition1 - condition2)

The `condition` column should contain the biological condition or treatment group name (e.g., control, treated, vehicle, drug). Each unique value creates one condition group. The pipeline requires **at least 2 unique condition values** to generate comparisons.

**Example design file with condition column:**

| assay | sample_id | r1 | r2 | scaling_group | condition |
|-------|-----------|----|----|---------------|-----------|
| ATAC | sample-ctrl-rep1 | ... | ... | default | control |
| ATAC | sample-ctrl-rep2 | ... | ... | default | control |
| ATAC | sample-treat-rep1 | ... | ... | default | treated |
| ATAC | sample-treat-rep2 | ... | ... | default | treated |

With `perform_comparisons: true`, this will generate:
- `bigwigs/{method}/aggregated/control.bigWig`
- `bigwigs/{method}/aggregated/treated.bigWig`
- `bigwigs/{method}/subtraction/control_vs_treated.bigWig`
- `bigwigs/{method}/subtraction/treated_vs_control.bigWig`

!!! note
    The `condition` column is independent of consensus groups. You can use both simultaneously:
    - Use `consensus_group` to merge samples for consensus peak calling (consensus)
    - Use `condition` to compare across conditions (comparison)
    - For RNA-seq, use `group`/`deseq2` for differential expression analysis in addition to condition-based bigwig comparisons

#### RNA-seq grouping for DESeq2

For RNA-seq experiments using spike-in normalization with DESeq2, the design command automatically detects experimental groups from sample names. Two columns are generated:

- **`group`**: The experimental group name (e.g., control, treated, WT, KO, vehicle, drug)
- **`deseq2`**: Binary encoding where 0 = control/reference group, 1 = treatment/comparison group

The tool detects groups using several strategies:

1. **Keyword detection**: Recognizes common keywords like control, treated, WT, KO, vehicle, DMSO
2. **Pattern extraction**: Extracts group information from sample naming patterns (e.g., `sample-GROUP-rep1`)
3. **Custom patterns**: Use `--deseq2-pattern` to specify a custom regex pattern for group extraction

**Example:**

For samples named:
- `rna-spikein-control-rep1_R1.fastq.gz`
- `rna-spikein-treated-rep1_R1.fastq.gz`

The generated design will include:

| assay | sample_id | r1 | r2 | scaling_group | group | deseq2 |
|-------|-----------|----|----|---------------|-------|--------|
| RNA | rna-spikein-control-rep1 | ... | ... | default | control | 0 |
| RNA | rna-spikein-treated-rep1 | ... | ... | default | treated | 1 |

The control/reference group is automatically identified and assigned `deseq2=0`, while treatment groups receive `deseq2=1`.

**Best Practices for Sample Naming:**

To ensure reliable automatic group detection, follow these naming conventions:

1. **Include group identifier before replicate number**:
   - Good: `sample-control-rep1`, `sample-treated-rep2`
   - Good: `batch1-WT-rep1`, `batch1-KO-rep2`
   - Avoid: `sample-rep1-control` (group after replicate)

2. **Use hyphens or underscores as separators**:
   - Good: `experiment-drug-day0-rep1` or `experiment_vehicle_day0_rep1`
   - Avoid: `experimentdrugday0rep1` (no separators)

3. **Use recognized keywords for control groups**:
   - Recognized: `control`, `ctrl`, `untreated`, `vehicle`, `mock`, `dmso`, `wt`, `wildtype`, upper or lower case.
   - Example: `sample-vehicle-rep1` will be automatically identified as the reference group

4. **Avoid ambiguous covariate naming**:
   - Good: `drug-day0-rep1`, `drug-day7-rep1` (group before timepoint)

5. **Be consistent across replicates**:
   - Good: `exp-control-rep1`, `exp-control-rep2`, `exp-treated-rep1`, `exp-treated-rep2`
   - Avoid: Mixing naming schemes between replicates

```bash
seqnado design rna fastqs/* --deseq2-pattern "-(WT|MUT)-"
```

This extracts groups from patterns like `sample-WT-day0_R1.fastq.gz` and `sample-MUT-day0_R1.fastq.gz`.

**Multi-Group Comparisons:**

The automatic binary encoding (`deseq2` column with 0/1) only works for 2-group comparisons (e.g., control vs treated). If your experiment has 3 or more groups (e.g., `DMSO-00hr`, `dTAG-00hr`, `dTAG-24hr`), the tool will:

1. Populate the `group` column with all detected groups
2. Leave the `deseq2` column empty
3. Display a warning message

For multi-group comparisons, you must manually edit the design CSV file to specify contrasts. Common approaches:

- **Reference-level coding**: Assign `0` to your reference group (e.g., control), and `1` to all other groups. This requires running separate contrasts pairwise.
- **Treatment contrasts**: Map each group to a numeric code (e.g., `0` = DMSO, `1` = dTAG-00hr, `2` = dTAG-24hr) — but note that standard DESeq2 design formulas expect binary columns, so this requires advanced configuration in the `config.yaml`.
- **Design matrix**: Use the `~0 + group` formula in your DESeq2 configuration to compare all groups simultaneously.

Consult the [Tools Reference](tools.md#deseq2) and the [Troubleshooting guide](troubleshooting.md) if you need help configuring multi-group contrasts.

#### Multiomics Mode
```bash
seqnado design 
```

For examples of additional options (auto-discovery, grouping, patterns), consult [seqnado design](cli.md#cli-seqnado-design).

---

**See Also:**

- [Pipeline Overview](pipeline.md) - Run your analysis
- [CLI Reference](cli.md#cli-seqnado-design) - Complete design command options
- [Troubleshooting](troubleshooting.md#design-files-seqnado-design) - Design file issues