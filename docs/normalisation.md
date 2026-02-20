# Normalisation Methods

SeqNado supports several normalisation strategies for generating coverage bigwig files. This page documents the exact formula used by each method, how per-sample factors are applied, and how factors are derived for merged consensus-group bigwigs.

## Quick Start (TL;DR)

**ChIP-seq with spike-in?** Use `orlando` for simplicity or `with_input` if you have paired input controls.  
**RNA-seq?** Use `deseq2` or `edgeR` — they correct for compositional bias.  
**No spike-in and expect similar global levels across conditions?** Use `csaw`.  
Need help choosing? See [Choosing a Method](#choosing-a-method) below.

## Summary Table

| Method | Assay | Input data | Per-sample formula | Merged formula | Signal Unit |
|---|---|---|---|---|---|
| **unscaled** | Any | BAM files | \(1\) | \(1\) |  As set in deeptools/bamnado/homer configuration |
| **orlando** | ChIP-seq, CUT&TAG, CUT&RUN | Spike-in BAM read counts | \(10^6 / S_{\text{ip}}\) | \(10^6 / \sum_i S_{\text{ip},i}\) | RRPM (spike-in reference-adjusted RPM) |
| **with_input** | ChIP-seq, CUT&TAG, CUT&RUN | Spike-in BAM read counts (IP + input) | \((S_{\text{ctrl}} \times 10^7) / (S_{\text{ip}} \times R_{\text{ctrl}})\) | \((\sum S_{\text{ctrl},i} \times 10^7) / (\sum S_{\text{ip},i} \times \sum R_{\text{ctrl},i})\) | RRPM (IP enrichment over spike-in-normalized input) |
| **deseq2** | RNA-seq | Spike-in gene count matrix | DESeq2 median-ratio \(\hat{s}_j\) | Arithmetic mean \(\bar{\hat{s}}\) | Library-size normalized |
| **edgeR** | RNA-seq | Spike-in gene count matrix | edgeR TMM \(\hat{f}_j^{\text{TMM}}\) | Arithmetic mean \(\bar{\hat{f}}^{\text{TMM}}\) | Library-size normalized |
| **csaw** | ChIP-seq, CUT&TAG, CUT&RUN, ATAC-seq (RNA-seq: see note) | Genomic bin read counts | \(\bar{L} / L_j\) | \(1 / \sum_i (1/s_i)\) | CPM (library-size normalized, within group) |

---

## Spike-in Normalisation

Spike-in normalisation uses a known amount of exogenous material (chromatin or RNA from a different species) added to each library before sequencing. Because the amount added is constant, the number of reads mapping to the spike-in reflects sequencing depth and can be used to normalise samples to a comparable scale. SeqNado implements three spike-in strategies — two based on total spike-in read counts (orlando, with_input) and two based on gene-level counting from a spike-in-aware count matrix (deseq2, edgeR).

### Orlando (Reads-Per-Million Spike-in)

**Reference:** Orlando DA, Chen MW, Brown VE, Solanki S, Choi YJ, Olson ER, Fritz CC, Bradner JE, Guenther MG. *Quantitative ChIP-Seq Normalization Reveals Global Modulation of the Epigenome.* Cell Reports. 2014;9(3):1163–1170. doi:[10.1016/j.celrep.2014.10.018](https://doi.org/10.1016/j.celrep.2014.10.018)

**Designed for:** ChIP-seq, CUT&TAG, CUT&RUN (any assay that uses chromatin spike-in).

**Concept:** Scale each sample so that one million spike-in reads would have been sequenced, making signal proportional to the absolute amount of chromatin immunoprecipitated.

!!! note "In plain terms"
    You added a fixed amount of exogenous chromatin (e.g. *Drosophila*) to each ChIP sample before sequencing. Samples that were sequenced more deeply will produce more spike-in reads, but the amount of spike-in chromatin was the same. By dividing every sample's signal by its spike-in read count, you remove the effect of unequal sequencing depth and make signals directly comparable across samples — a bigger signal in the bigwig means more protein was genuinely bound, not just that more reads were generated.

**Per-sample formula:**

$$
\text{scale_factor} = \frac{10^6}{S_{\text{ip}}}
$$

where \(S_{\text{ip}}\) is the number of reads aligning to the spike-in genome for that sample. This scales each sample so that 1 million spike-in reads = 1 unit of scale factor. A sample with twice as many spike-in reads gets a factor of 0.5 (downscaled) so that its signal looks equivalent to a sample sequenced less deeply.

**Merged bigwig:** The merged BAM contains reads from all samples in the consensus group. The factor is computed from the pooled spike-in counts — equivalent to performing the normalisation on the merged BAM directly:

$$
\text{scale_factor}_{\text{merged}} = \frac{10^6}{\displaystyle\sum_i S_{\text{ip},i}}
$$

!!! note "In plain terms"
    For the merged track, we add up the spike-in reads from all samples in the group and use that total. This is equivalent to asking: "if we treated all these libraries as one big pooled experiment, what would the spike-in normalisation factor be?" Using the sum (rather than averaging the per-sample factors) gives the correct scaling factor when samples differ in sequencing depth.

**Output:** `seqnado_output/{assay}/resources/orlando/normalisation_factors.json`

---

### With-Input (Spike-in Input-Normalised)

**Designed for:** ChIP-seq, CUT&TAG, CUT&RUN (requires both a chromatin spike-in and a paired input control sample).

**Concept:** Corrects for both sequencing depth (via spike-in) and IP efficiency (via the paired input control), so the final signal represents the fraction of chromatin that is specifically immunoprecipitated relative to background.

!!! note "In plain terms"
    The Orlando method removes differences in sequencing depth, but two IP experiments can still appear different if one had more efficient immunoprecipitation (more protein pulled down from the same DNA input). The with-input method adds a second correction: it compares the spike-in-normalised IP signal to the spike-in-normalised input control signal. Since input reflects total chromatin (background), dividing IP by input gives you the enrichment: how much more DNA is bound in the IP versus the random background. The result is a direct measure of specific protein binding, robust to both sequencing depth and IP pulldown efficiency differences between samples.

**Per-sample formula:**

For each IP sample paired with its input control, define the spike-in-normalised read density for the IP (\(\rho_{\text{ip}}\)) and the input control (\(\rho_{\text{ctrl}}\)):

$$
\rho_{\text{ip}} = \frac{R_{\text{ip}}}{S_{\text{ip}}}, \qquad \rho_{\text{ctrl}} = \frac{R_{\text{ctrl}}}{S_{\text{ctrl}}}
$$

where \(R\) is the number of reads mapping to the reference genome and \(S\) is the number of reads mapping to the spike-in genome. The relative signal (IP enrichment over background) is then normalised to a reads-per-million scale:

$$
\text{scale_factor} = \frac{\rho_{\text{ip}}}{\rho_{\text{ctrl}}} \times \frac{10^7}{R_{\text{ip}}}
$$

Expanding:

$$
\text{scale_factor} = \frac{S_{\text{ctrl}} \times 10^7}{S_{\text{ip}} \times R_{\text{ctrl}}}
$$

**Important:** This formula applies only to IP samples. Input control samples are treated separately — they receive a scale factor of 1 (unscaled), as they do not undergo the IP/input correction. When input samples are visualised as tracks (if configured to do so), they appear unscaled.

**Merged bigwig:** The per-sample normalisation table (`resources/with_input/normalisation_factors.tsv`) already contains the IP and paired-input read counts for each sample. SeqNado sums these across all samples in the group and applies the same formula as for a single sample — mathematically equivalent to running the calculation on the merged BAM directly:

$$
\text{scale_factor}_{\text{merged}} = \frac{\displaystyle\sum_i S_{\text{ctrl},i} \times 10^7}{\displaystyle\sum_i S_{\text{ip},i} \times \displaystyle\sum_i R_{\text{ctrl},i}}
$$

!!! note "In plain terms"
    For the merged track, we pool the spike-in and reference read counts from all samples in the group and apply the same formula as for a single sample. This is mathematically identical to running the normalisation on the merged BAM file.

**Output:** `seqnado_output/{assay}/resources/with_input/normalisation_factors.json`

---

### DESeq2 Size Factors (Spike-in Genes)

**Reference:** Love MI, Huber W, Anders S. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biology. 2014;15(12):550. doi:[10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)

**Designed for:** RNA-seq.

**Concept:** Uses DESeq2's `estimateSizeFactors()` to derive normalisation factors from a count matrix of spike-in genome features. If spike-in gene names are provided in the configuration, those genes are used as `controlGenes`; otherwise the standard median-ratio method is applied to all genes.

!!! note "In plain terms"
    Instead of counting all spike-in reads in bulk, this method counts reads falling on individual spike-in genes (annotated features). For each gene, it calculates how that gene's count in each sample compares to a reference (the geometric mean across all samples), then takes the **median** of these per-gene ratios as the normalisation factor. Why median instead of mean? Because a few extremely highly expressed spike-in genes would skew a mean-based estimate; the median is robust to outliers. This approach—comparing each sample to a geometric mean reference—is the same logic used widely in RNA-seq differential expression analysis and handles compositional bias well.

**Per-sample formula:** DESeq2 size factors are the median of per-gene count ratios relative to the geometric mean across all samples:

$$
\hat{s}_j = \operatorname{median}_g \left( \frac{k_{gj}}{\left(\prod_{j'} k_{gj'}\right)^{1/J}} \right)
$$

where \(k_{gj}\) is the count for gene \(g\) in sample \(j\) and \(J\) is the total number of samples. The bigwig generator divides the raw signal by \(\hat{s}_j\) to equalise depth across samples.

**Merged bigwig:** Arithmetic mean of per-sample size factors (approximation; re-fitting the DESeq2 model is not applicable):

$$
\text{scale_factor}_{\text{merged}} = \frac{1}{N} \sum_i \hat{s}_i
$$

!!! note "In plain terms"
    DESeq2 fits a statistical model to estimate how each gene's count depends on the sample and other factors. This model requires replicate samples to work correctly—it fits the full data to estimate both per-sample size factors and per-gene variability. When you merge replicates into a single BAM, you lose the replicate structure, so re-fitting the model would not be statistically valid. Instead, SeqNado uses the arithmetic mean of the per-sample factors as a practical approximation. The merged track represents a pooled, higher-coverage view of the group signal, useful for visualisation, but is not itself a statistically normalised sample in the same sense that individual replicates are.

**Output:** `seqnado_output/{assay}/resources/deseq2/normalisation_factors.json`

---

### edgeR TMM (Spike-in Genes)

**Reference:** Robinson MD, Oshlack A. *A scaling normalization method for differential expression analysis of RNA-seq data.* Genome Biology. 2010;11(3):R25. doi:[10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25)

**Designed for:** RNA-seq.

**Concept:** Uses edgeR's Trimmed Mean of M-values (TMM) normalisation (`calcNormFactors()`) computed on spike-in gene counts. If spike-in gene names are provided, TMM is computed using only those genes as the reference subset; otherwise all genes are used. After factor calculation the spike-in genes are removed from the count matrix.

!!! note "In plain terms"
    TMM works similarly to DESeq2 but with a different approach. For each sample, compare its per-gene counts to a reference sample on a log scale (M-values = log fold-changes). Average these comparisons across genes to estimate library size differences. The "trimming" step removes the most extreme genes (those with the most dramatic fold-changes) before averaging, preventing a handful of highly expressed genes from dominating the estimate. This is particularly useful for RNA-seq where a few genes can dominate the read count. The result is a robust library size correction that handles compositional bias.

**Per-sample formula:** For each sample \(j\) relative to a reference sample \(r\), TMM computes a weighted mean of log fold changes after trimming extreme values:

$$
\log_2 \hat{f}_j^{\text{TMM}} = \frac{\displaystyle\sum_{g \in \mathcal{G}^*} w_{gj} \, M_{gj}^r}{\displaystyle\sum_{g \in \mathcal{G}^*} w_{gj}}
$$

where \(M_{gj}^r = \log_2\!\left(\tfrac{k_{gj}/N_j}{k_{gr}/N_r}\right)\) is the log fold change, \(w_{gj}\) is an inverse-variance weight, and \(\mathcal{G}^*\) is the set of genes remaining after trimming. The effective library size is \(N_j \cdot \hat{f}_j^{\text{TMM}}\).

**Merged bigwig:** Arithmetic mean of per-sample TMM factors (approximation; same rationale as DESeq2 — the merged BAM has no replicates for the model to work with):

$$
\text{scale_factor}_{\text{merged}} = \frac{1}{N} \sum_i \hat{f}_i^{\text{TMM}}
$$

!!! note "In plain terms"
    Like DESeq2, TMM factors come from a model that is designed to work across replicate samples. Merging replicates into a single BAM removes the replicate structure that the model depends on for statistical inference. SeqNado uses the average of the per-sample TMM factors as a practical approximation. The merged track serves as a high-coverage visualisation aid (cleaner, less noisy) rather than as a statistically normalised sample.

**Output:** `seqnado_output/{assay}/resources/edgeR/normalisation_factors.json`

---

## CSAW Library-Size Scaling

**Reference:** Lun ATL, Smyth GK. *csaw: a Bioconductor package for differential binding analysis of ChIP-seq data using sliding windows.* Nucleic Acids Research. 2016;44(5):e45. doi:[10.1093/nar/gkv1191](https://doi.org/10.1093/nar/gkv1191)

**Designed for:** ChIP-seq, CUT&TAG, CUT&RUN, ATAC-seq. Technically applicable to RNA-seq (it is equivalent to reads-per-million normalisation), but DESeq2 or edgeR are strongly preferred for RNA-seq because they additionally correct for compositional bias — see [Choosing a method](#choosing-a-method) below.

**Concept:** Equalise read depth across samples within a scaling group by comparing the total number of reads falling into randomly sampled genomic bins. This normalises for library size differences without relying on exogenous spike-in material.

!!! note "In plain terms"
    Sometimes you don't have a spike-in, but you still need to make samples comparable. CSAW solves this by counting reads in many small windows tiled across the genome and asking: "how many reads did each sample produce in total?" Samples that were sequenced more deeply get scaled down and shallower ones get scaled up, so that all samples appear to have the same total read count. This corrects for sequencing depth differences but **not** for genuine biological differences in the amount of immunoprecipitated chromatin — it is most appropriate when you expect the global level of the mark to be similar across conditions.

**Per-sample formula:**

Reads from each sample's BAM are counted in fixed-size genomic bins using featureCounts. Within a scaling group of \(N\) samples:

$$
\text{scale_factor}_j = \frac{\bar{L}}{L_j}, \qquad \bar{L} = \frac{1}{N}\sum_i L_i
$$

where \(L_j\) is the total bin read count (library size) for sample \(j\) and \(\bar{L}\) is the group mean.

**Interpretation of factors:**
- Samples with larger libraries receive factors \(<1\) (downscaled) — these were sequenced more deeply relative to the group average.
- Samples with smaller libraries receive factors \(>1\) (upscaled) — these were sequenced less deeply relative to the group average.
- Factors always scale towards a common library size equal to the group mean.

**Merged bigwig:** The merged BAM is the physical concatenation of all per-sample BAMs, so its total library size is \(\sum_i L_i\). Because each per-sample factor is \(s_i = \bar{L} / L_i\), we have \(L_i = \bar{L} / s_i\), and the correct scale factor for the merged BAM is:

$$
\text{scale_factor}_{\text{merged}} = \frac{\bar{L}}{\displaystyle\sum_i L_i} = \frac{1}{\displaystyle\sum_i \tfrac{1}{s_i}}
$$

The merged factor is the **harmonic mean** of the per-sample factors.

!!! note "In plain terms"
    When you merge replicates into a single BAM, that file contains all reads combined. A merged BAM of three replicates has roughly 3× the reads of any single replicate. If you scaled it like an individual replicate, it would look 3× taller — not useful for comparison. Instead, SeqNado uses the **harmonic mean** to ensure the merged track sits at the right scale.

    **Concrete example:** Three replicates sequenced at 10M, 20M, and 30M reads (group mean = 20M). Per-sample factors: [20/10, 20/20, 20/30] = [2.0, 1.0, 0.67]. Merged BAM has 60M total reads, so correct merged factor should be 20/60 = 0.33. Using arithmetic mean: (2.0+1.0+0.67)/3 = 1.22 ✗ **wrong**. Using harmonic mean: 1/(1/2.0 + 1/1.0 + 1/0.67) = 1/(0.5+1.0+1.5) = 1/3.0 = 0.33 ✓ **correct**. Harmonic mean works because scaling factors relate inversely to library size—when you concatenate BAMs, you need a formula that respects this inverse relationship.

Scaling factors are stored in `seqnado_output/{assay}/resources/{group}_scaling_factors.tsv`.

---

## Unscaled (No Normalisation)

When normalisation is not enabled, bigwigs are generated with no scaling factor applied. Signal represents raw read pileup, which is affected by sequencing depth. Use this only for exploratory analysis or when samples are known to have matched sequencing depth.

!!! note "In plain terms"
    No adjustment is made. A sample sequenced to twice the depth will appear twice as tall in the genome browser, regardless of whether there is a genuine biological difference. Only use this if your samples were sequenced to the same depth, or if you just want a quick look at the data without worrying about comparability.

---

## Choosing a Method

### By assay type

| Assay | Recommended methods | Notes |
|---|---|---|
| **ChIP-seq / CUT&TAG / CUT&RUN** | `orlando`, `with_input`, `csaw` | Prefer `orlando` or `with_input` when a chromatin spike-in was added — they are purpose-built and more direct than DESeq2/edgeR for this use case. Use `csaw` when no spike-in is available and you expect similar global levels across conditions. |
| **ATAC-seq** | `csaw` | Library-size normalization is most appropriate for ATAC-seq. |
| **RNA-seq** | `deseq2`, `edgeR` | These are the standard methods for RNA-seq. `orlando` and `with_input` are not applicable. `csaw` is technically valid but does not correct for compositional bias — see below. |

### Can CSAW be used for RNA-seq?

Technically yes — the formula is equivalent to reads-per-million (RPM) normalisation, which is a valid and widely understood approach. However, **DESeq2 or edgeR are strongly preferred** for RNA-seq for the following reason:

In an RNA-seq library, a small number of very highly expressed genes can account for a large fraction of all reads. This means two libraries can appear to have different "total read counts" even if the underlying biology is the same, simply because a few dominant genes are present. DESeq2 and edgeR correct for this **compositional bias** by comparing each gene's count to a reference derived from across many genes, rather than relying on the raw library total. CSAW (and RPM) do not account for this and can produce misleading normalisations in RNA-seq data.

### Can Orlando be used for RNA-seq?

It depends on your spike-in strategy:

- **Synthetic RNA spike-ins (e.g. ERCC):** Supported if the spike-in sequences are added to the reference genome used for alignment. Reads mapping to the spike-in sequences are then split out in the same way as chromatin spike-in reads, and the Orlando formula applies directly. Configure this by including the spike-in sequences in your genome reference.
- **Exogenous chromatin spike-in (e.g. *Drosophila* cells added to the RNA-seq library):** Conceptually yes — if reads are aligned to a combined reference+spike-in genome and the spike-in reads are split out, the Orlando formula applies in exactly the same way as for ChIP-seq. Some labs use this approach to detect global changes in transcriptional output. If your RNA-seq experiment was prepared this way, Orlando normalisation is valid and the SeqNado pipeline should handle it correctly provided the genome configuration includes both species.

---

## File Locations

| Output | Path |
|---|---|
| Orlando factors | `seqnado_output/{assay}/resources/orlando/normalisation_factors.json` |
| With-input factors | `seqnado_output/{assay}/resources/with_input/normalisation_factors.json` |
| DESeq2 factors | `seqnado_output/{assay}/resources/deseq2/normalisation_factors.json` |
| edgeR factors | `seqnado_output/{assay}/resources/edger/normalisation_factors.json` |
| CSAW factors | `seqnado_output/{assay}/resources/{group}_scaling_factors.tsv` |

---

## Background Terminology

This document assumes familiarity with high-throughput sequencing concepts. If you are new to ChIP-seq or RNA-seq normalisation, these terms may be useful:

- **Spike-in**: Known quantity of exogenous material (chromatin from another species, or synthetic RNA) added to each sample before sequencing. Because the amount is constant, reads mapping to spike-in reflect sequencing depth, not biology.
- **IP (immunoprecipitation)**: In ChIP-seq, the antibody pulls down target protein and its bound DNA. The "IP" sample contains what was pulled down; the "input" sample is unselected chromatin (background).
- **Reference genome**: The genome of interest (e.g., *Homo sapiens*). In spike-in experiments, reads are split by which genome they map to.
- **Replicates**: Multiple independent experiments. Normalisation models (DESeq2, edgeR) use replicates to estimate variability; merged BAMs (pooled replicates) have no replicate structure, so these models cannot be re-fit.
- **Compositional bias** (RNA-seq): When highly expressed genes make up a large fraction of all reads, raw library size estimates become misleading. DESeq2/edgeR handle this by comparing genes relative to a geometric mean reference.
- **Geometric mean**: The $N$-th root of the product of $N$ values — used in DESeq2 as a stable reference that is not affected by any single highly expressed gene.