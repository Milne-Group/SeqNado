[← Back to main page](index.md)

# Design Guide

The `seqnado design` command generates a design CSV file from FASTQ files for a specific assay. If no assay is provided, the tool operates in multiomics mode.

For full arguments and flags, see the CLI reference: [seqnado design](cli.md#cli-seqnado-design).

### Example Usage

#### Generate a Design CSV for ATAC-seq
```bash
seqnado design atac fastqs/* 
```

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



#### Complex Case with Multiple Controls and ambiguity in pairing

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
seqnado design chip fastqs/* --ip-to-control "H3K4me3:sf-input,Menin:df-input"
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
| RNA | rna-spikein-c
ontrol-rep1 | ... | ... | default | control | 0 |
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

For multi-group comparisons, you'll need to manually configure the `deseq2` column in your metadata CSV based on your specific contrasts.

#### Multiomics Mode
```bash
seqnado design 
```

For examples of additional options (auto-discovery, grouping, patterns), consult [seqnado design](cli.md#cli-seqnado-design).

## Next Steps

Once you have designed your experiment, proceed to set up the pipeline:

[Pipeline](pipeline.md)