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