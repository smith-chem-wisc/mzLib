# Statistical report

Analysis `adc8979a883969aa` (`DEF-DIFF v1`). Generated from `DifferentialResults.metadata.json` and `DifferentialResults.tsv`; nothing here is written by hand.

## 1. Methods

**Software.** mzLib 9.9.0; datarepo engine 0.1.0. Method version: diff-1.0.

**Inputs.** observation_table `observations.tsv` (sha256 `1111111111111111...`); design `design.tsv` (sha256 `2222222222222222...`).

**Transformation and normalization.** All quantities were analysed on the log2 scale. Stratum `all`, basis `mbr_kept`: each sample was shifted by its median log2 difference from the across-sample reference, over the 812 peptides with a value in every sample that has values (at least 100 required). Stratum `all`, basis `msms_only`: each sample was shifted by its median log2 difference from the across-sample reference, over the 790 peptides with a value in every sample that has values (at least 100 required). The per-sample shifts are in the metadata file.

**Missing values.** Missing values were not imputed: each model was fitted on the observed values only. A missing value is an empty cell in the results file, never 0.

**Models.**
- Stratum `all`, basis `mbr_kept`, `protein_group`, `moderated_t` (method `moderated`): `y ~ age`, fitted without REML; degrees of freedom by `moderated_residual`; a trended empirical-Bayes variance prior (`moments_legacy`) with 3.1 degrees of freedom, its variance depending on the average intensity (see the mean-variance figure), fitted on 240 features.
- Stratum `all`, basis `mbr_kept`, `protein_group`, `peptide_mixed_model` (method `moderated`): `y ~ peptide + age + (1|individual/sample)`, fitted by REML; degrees of freedom by `satterthwaite_moderated`; robust weighting: huber 1.345, MAD about 0, to convergence; an empirical-Bayes variance prior (`moments_legacy`) with 4.2 degrees of freedom and variance 0.08, fitted on 1,812 features.

Outlying values were down-weighted (robust weighting on).

**Tests.** Every comparison is tested two-sided. Each estimate is reported with its standard error, a 95% confidence interval (t-based, on the stated degrees of freedom), its test statistic, its degrees of freedom and an exact p-value.

**Multiple testing.** p-values were adjusted with the Benjamini–Hochberg method within each family of features sharing stratum, grain, quantity, quant basis, comparison and method: 2 families, their sizes in the metadata file.

**Thresholds.** No significance threshold was applied and no feature is labelled significant; the counts at adjusted p < 0.05 in section 3 are descriptive only (ASA statements, 2016 and 2019).

**Seed.** Random seed: 42.

## 2. Design

### Stratum `all`

Factors: none (the whole dataset). Chosen by: default_list.

| level | samples | individuals |
|---|---|---|
| age=old | 6 | 5 |
| age=young | 6 | 6 |

- Batches: b1, b2.
- Fractions: 1. Technical replicates: 1.
- Reference level: age = young.
- Samples without values: s12.
- Comparisons not run here: none.

### Comparisons

| id | comparison | numerator | denominator | covariate |
|---|---|---|---|---|
| c1 | age=old vs age=young | age=old | age=young |  |

### Design warnings

- biological replicate 3 of condition old is absent

## 3. Results

| stratum | basis | grain | quantity | comparison | method | features | fitted | not fitted | adjusted p < 0.05 |
|---|---|---|---|---|---|---|---|---|---|
| all | mbr_kept | protein_group | abundance | c1 | moderated | 22 | 20 | 2 | 1 |
| all | msms_only | protein_group | abundance | c1 | moderated | 22 | 20 | 2 | 1 |

Not fitted, by reason:

| stratum | basis | comparison | method | status | features |
|---|---|---|---|---|---|
| all | mbr_kept | c1 | moderated | absent_in_numerator | 1 |
| all | mbr_kept | c1 | moderated | below_support:min_2_per_side | 1 |
| all | msms_only | c1 | moderated | absent_in_numerator | 1 |
| all | msms_only | c1 | moderated | below_support:min_2_per_side | 1 |

Counts at adjusted p < 0.05 are descriptive only; no feature is labelled significant.

## 4. Diagnostics

### Figure 1. p-value histogram: stratum `all`, basis `mbr_kept`, `protein_group` `abundance`, comparison `c1`, method `moderated`

![Figure 1](figures/figure_01_pvalue_histogram.svg)

p-values of the 20 fitted features in 20 bins of width 0.05 (bin k holds 0.05k ≤ p < 0.05(k+1); p = 1 is in the last bin). With no change anywhere the bars are flat; a spike near 0 over a flat floor is expected when some features change; a U shape, or a slope rising toward 1, signals a problem with the model or the data. Data: [`figure_01_pvalue_histogram.tsv`](figures/figure_01_pvalue_histogram.tsv).

### Figure 2. Volcano: stratum `all`, basis `mbr_kept`, `protein_group` `abundance`, comparison `c1`, method `moderated`

![Figure 2](figures/figure_02_volcano.svg)

log2 effect against -log10(p) for the 20 fitted features; no threshold lines are drawn. Data: [`figure_02_volcano.tsv`](figures/figure_02_volcano.tsv).

### Figure 3. p-value histogram: stratum `all`, basis `msms_only`, `protein_group` `abundance`, comparison `c1`, method `moderated`

![Figure 3](figures/figure_03_pvalue_histogram.svg)

p-values of the 20 fitted features in 20 bins of width 0.05 (bin k holds 0.05k ≤ p < 0.05(k+1); p = 1 is in the last bin). With no change anywhere the bars are flat; a spike near 0 over a flat floor is expected when some features change; a U shape, or a slope rising toward 1, signals a problem with the model or the data. Data: [`figure_03_pvalue_histogram.tsv`](figures/figure_03_pvalue_histogram.tsv).

### Figure 4. Volcano: stratum `all`, basis `msms_only`, `protein_group` `abundance`, comparison `c1`, method `moderated`

![Figure 4](figures/figure_04_volcano.svg)

log2 effect against -log10(p) for the 20 fitted features; no threshold lines are drawn. Data: [`figure_04_volcano.tsv`](figures/figure_04_volcano.tsv).

### Figure 5. Mean-variance: stratum `all`, basis `mbr_kept`, `protein_group`, `moderated_t`

![Figure 5](figures/figure_05_mean_variance.svg)

Each feature's sqrt(sigma), the square root of its residual standard deviation, against its average log2 intensity, for 20 features, as limma's plotSA draws it. The line is the empirical-Bayes prior at the same scale ((s0^2)^(1/4)): a curve, because the prior is trended. Data: [`figure_05_mean_variance.tsv`](figures/figure_05_mean_variance.tsv).

### Figure 6. Mean-variance: stratum `all`, basis `mbr_kept`, `protein_group`, `peptide_mixed_model`

![Figure 6](figures/figure_06_mean_variance.svg)

Each feature's sqrt(sigma), the square root of its residual standard deviation, against its average log2 intensity, for 20 features, as limma's plotSA draws it. The line is the empirical-Bayes prior at the same scale ((s0^2)^(1/4)): flat, because the prior is constant. Data: [`figure_06_mean_variance.tsv`](figures/figure_06_mean_variance.tsv).

### Figure 7. Sample correlation: stratum `all`, basis `mbr_kept`

![Figure 7](figures/figure_07_sample_correlation.svg)

Pearson correlation of log2 values between each pair of the 4 samples, over the features both samples have; r needs at least 3 such features, and fewer leave the cell empty. The table gives r and the number of features. Data: [`figure_07_sample_correlation.tsv`](figures/figure_07_sample_correlation.tsv).

### Figure 8. Global shift per comparison, stratum and basis

![Figure 8](figures/figure_08_global_shift.svg)

Each comparison's global shift: the median, over peptides with a value in every sample on both sides, of the difference in mean log2 intensity, before normalization; 0 means the two sides sit level. The number of peptides each rests on is beside it. Data: [`figure_08_global_shift.tsv`](figures/figure_08_global_shift.tsv).

## 5. Reporting checklist

| item | asked for by | answered in | this run |
|---|---|---|---|
| Software and versions | JPR; PROTEOMICS; MIAPE-Quant 4.1 | 1. Methods | Y |
| Data transformation | JPR; PROTEOMICS; MIAPE-Quant 4.4 | 1. Methods | Y |
| Normalization | JPR; PROTEOMICS | 1. Methods | Y |
| Missing values | PROTEOMICS; MIAPE-Quant 4.4 | 1. Methods | Y |
| Statistical tests and sidedness | JPR; PROTEOMICS; Nature reporting summary | 1. Methods | Y |
| Degrees of freedom | PROTEOMICS; Nature reporting summary | 1. Methods; results file `df` | Y |
| Multiple-testing correction | Nature reporting summary; MIAPE-Quant 4.5 | 1. Methods | Y |
| Exact n per group | Nature reporting summary; MCP | 2. Design | Y |
| Biological and technical replicates | MCP; Nature reporting summary | 2. Design | Y |
| Effect sizes with confidence intervals | Nature reporting summary; PROTEOMICS; JPR; MIAPE-Quant 4.5 | results file `log2_effect`, `ci_low`, `ci_high` | Y |
| Exact p-values | Nature reporting summary; ASA 2016 | results file `p_value` | Y |
| Thresholds and their justification | MIAPE-Quant 4.5; PROTEOMICS; ASA 2019 | 1. Methods | Y |
| Priors | Nature reporting summary | 1. Methods | Y |
| Level of each test in nested designs | Nature reporting summary | 1. Methods (model formulas) | Y |
| Power analysis | JPR (where appropriate) | not computed | N |

An N is a gap in this run, shown rather than left out.
