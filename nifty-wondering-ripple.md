# Plan: Generic Condition-Based Bigwig Averaging and Comparisons

## Context

MCC already implements condition-based bigwig averaging (`compare.smk`) and subtraction comparisons, but:
1. This is MCC-specific with no equivalent for other assays (ATAC, ChIP, CAT, RNA, etc.)
2. Comparisons in MCC auto-enable when `len(condition_groups) >= 2` — no user opt-in/opt-out

The goal is to generalise this pattern to **all assays with bigwig output** using bamnado as the aggregation/comparison tool. The `condition` column already exists in metadata and is parsed into `condition_groups` in the Snakefile. Only bamnado is used for averaging/subtraction (consistent with MCC).

Rules:
- Only generate comparisons when `perform_comparisons: true` AND `len(condition_groups) >= 2`
- If spike-in is configured, use spikein-normalised bigwigs as input; otherwise use unscaled
- Output paths: `bigwigs/{pileup_method}/aggregated/{condition}.bigWig` and `bigwigs/{pileup_method}/subtraction/{cond1}_vs_{cond2}.bigWig`
- For spikein input, follow existing path convention: `bigwigs/{pileup_method}/spikein/{spikein_method}/aggregated/{condition}.bigWig`

---

## Changes Required

### 1. Config: `seqnado/config/configs.py`

Add `perform_comparisons: bool = False` to **two** existing models:

**`BigwigConfig`** (for non-MCC assays):
```python
class BigwigConfig(BaseModel):
    pileup_method: list[PileupMethod] | None = None
    binsize: int | None = None
    scale_methods: list[str] | None = None
    perform_comparisons: bool = False  # ADD THIS
```

**`MCCConfig`** (for MCC consistency):
```python
class MCCConfig(BaseModel, PathValidatorMixin):
    ...
    create_replicate_bigwigs: bool = False
    create_replicate_contact_files: bool = False
    perform_comparisons: bool = False  # ADD THIS
```

---

### 2. New Snakemake rules file: `seqnado/workflow/rules/pileup/pileup_condition_compare.smk`

Single file handling all non-MCC assays. Rules gated on `CONFIG.assay_config.bigwigs.perform_comparisons` and presence of bamnado config.

**Helper function** (at top of file):
```python
def get_condition_input_bigwigs(wildcards, pileup_method, spikein_method=None):
    """Return per-sample bigwig paths for a condition group."""
    samples = SAMPLE_GROUPINGS.get_grouping("condition").get_samples(wildcards.condition)
    if spikein_method:
        return expand(
            OUTPUT_DIR + f"/bigwigs/{pileup_method}/spikein/{spikein_method}/{{sample}}.bigWig",
            sample=samples
        )
    return expand(
        OUTPUT_DIR + f"/bigwigs/{pileup_method}/unscaled/{{sample}}.bigWig",
        sample=samples
    )
```

**Four rules** (all gated on `CONFIG.third_party_tools.bamnado is not None and CONFIG.assay_config.bigwigs and CONFIG.assay_config.bigwigs.perform_comparisons`):

**Rule: `make_bigwigs_aggregated`** (unscaled → condition mean)
- Input: `bigwigs/{pileup_method}/unscaled/{sample}.bigWig` expanded over condition group
- Output: `bigwigs/{pileup_method}/aggregated/{condition}.bigWig`
- Wildcard: `pileup_method="|".join(...)`, `condition="|".join(SAMPLE_GROUPINGS.get_grouping("condition").group_names)`
- Shell: `bamnado bigwig-aggregate --bigwigs {input.bigwigs} -o {output.bigwig} -m mean`

**Rule: `make_bigwigs_spikein_aggregated`** (spikein-normalised → condition mean; only when has_spikein)
- Input: `bigwigs/{pileup_method}/spikein/{spikein_method}/{sample}.bigWig` expanded over condition group
- Output: `bigwigs/{pileup_method}/spikein/{spikein_method}/aggregated/{condition}.bigWig`
- Additional gate: `if getattr(CONFIG.assay_config, "has_spikein", False):`

**Rule: `make_bigwigs_subtraction`** (condition average → pairwise subtraction)
- Input: `bigwigs/{pileup_method}/aggregated/{condition1}.bigWig`, `.../aggregated/{condition2}.bigWig`
- Output: `bigwigs/{pileup_method}/subtraction/{condition1}_vs_{condition2}.bigWig`
- Shell: `bamnado bigwig-compare --bw1 {input.bw1} --bw2 {input.bw2} -o {output.bigwig} -c subtraction`

**Rule: `make_bigwigs_spikein_subtraction`** (spikein-averaged → pairwise subtraction)
- Input/output follow spikein path structure: `.../spikein/{spikein_method}/aggregated/...` → `.../spikein/{spikein_method}/subtraction/...`
- Additional gate: `if getattr(CONFIG.assay_config, "has_spikein", False):`

All rules must include `log:`, `benchmark:`, `message:`, `resources:` (using `define_time_requested`, `define_memory_requested`), and `container: "docker://ghcr.io/alsmith151/bamnado:latest"`.

---

### 3. Update `seqnado/workflow/rules/pileup/dna.smk`

Add one line:
```python
include: "pileup_condition_compare.smk"
```

Also include it in `rna.smk` (RNA also has bamnado bigwigs):
```python
include: "pileup_condition_compare.smk"
```

---

### 4. Update `seqnado/workflow/rules/mcc/pileup.smk`

In `confirm_bigwigs_generated`, change the `subtractions` gate from auto-enable (based only on condition count) to require the config flag:

```python
subtractions=[
    ...
] if (
    getattr(CONFIG.assay_config.mcc, "perform_comparisons", False)
    and len(SAMPLE_GROUPINGS.get_grouping("condition").group_names) >= 2
) else [],
```

---

### 5. Update `seqnado/outputs/core.py`

Add method `add_condition_bigwig_files()` to `SeqnadoOutputBuilder`:
```python
def add_condition_bigwig_files(self) -> None:
    condition_groups = self.sample_groupings.get_grouping("condition")
    if not condition_groups or len(condition_groups.group_names) < 2:
        return
    pileup_methods = self.config.assay_config.bigwigs.pileup_method or []
    files = []
    for method in pileup_methods:
        method_name = method.value  # e.g. "bamnado"
        for cond in condition_groups.group_names:
            files.append(f"{self.output_dir}/bigwigs/{method_name}/aggregated/{cond}.bigWig")
        for c1, c2 in itertools.permutations(condition_groups.group_names, 2):
            files.append(f"{self.output_dir}/bigwigs/{method_name}/subtraction/{c1}_vs_{c2}.bigWig")
    self.file_collections.append(BasicFileCollection(files=files))
```

Call from `create_output_builder()`:
```python
if (
    self.assay_config.create_bigwigs
    and self.assay not in (Assay.MCC, Assay.SNP, Assay.CRISPR)
    and getattr(self.assay_config.bigwigs, "perform_comparisons", False)
    and self.sample_groupings
):
    builder.add_condition_bigwig_files()
```

---

## File Paths Summary

| File | Action |
|------|--------|
| `seqnado/config/configs.py` | Add `perform_comparisons` to `BigwigConfig` and `MCCConfig` |
| `seqnado/workflow/rules/pileup/pileup_condition_compare.smk` | **CREATE** |
| `seqnado/workflow/rules/pileup/dna.smk` | Add include |
| `seqnado/workflow/rules/pileup/rna.smk` | Add include |
| `seqnado/workflow/rules/mcc/pileup.smk` | Gate subtraction on `perform_comparisons` |
| `seqnado/outputs/core.py` | Add `add_condition_bigwig_files()` + call it |

---

## Output Path Conventions

| Bigwig type | Path pattern |
|-------------|-------------|
| Condition mean (unscaled source) | `bigwigs/{pileup_method}/aggregated/{condition}.bigWig` |
| Condition mean (spikein source) | `bigwigs/{pileup_method}/spikein/{spikein_method}/aggregated/{condition}.bigWig` |
| Subtraction (unscaled source) | `bigwigs/{pileup_method}/subtraction/{cond1}_vs_{cond2}.bigWig` |
| Subtraction (spikein source) | `bigwigs/{pileup_method}/spikein/{spikein_method}/subtraction/{cond1}_vs_{cond2}.bigWig` |
| MCC subtraction (unchanged) | `bigwigs/mcc/subtractions/{cond1}_vs_{cond2}_{viewpoint_group}.bigWig` |

---

## Key Reused Patterns

- `bamnado bigwig-aggregate -m mean` and `bamnado bigwig-compare -c subtraction` — from `mcc/compare.smk:29-78`
- `SAMPLE_GROUPINGS.get_grouping("condition").get_samples(wildcards.condition)` — from `mcc/compare.smk:7`
- `itertools.permutations` for all-pairs (all directed pairs, same as MCC's `itertools.product` with `group1 != group2`) — from `mcc/pileup.smk:117`
- `define_memory_requested` / `define_time_requested` — `seqnado.workflow.helpers.common`
- Container: `"docker://ghcr.io/alsmith151/bamnado:latest"`

---

## Verification

1. **Config**: `python -c "from seqnado.config.configs import BigwigConfig; print(BigwigConfig(perform_comparisons=True))"`
2. **Rule directive test**: `pytest tests/unit/test_snakemake_rule_directives.py`
3. **Dry-run** with test data that has `condition` column and `perform_comparisons: true` in config + bamnado enabled
