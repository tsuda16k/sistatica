# Post hoc tests for the Friedman test (paired Wilcoxon + Holm p-values)

Performs pairwise post hoc comparisons for a one-way Friedman design.

## Usage

``` r
friedmanPost(
  df,
  id = "id",
  group = "group",
  y = "y",
  alpha = 0.05,
  pAdjust = "holm",
  ciAdjust = c("bonferroni", "none")
)
```

## Arguments

- df:

  A long-format data frame.

- id:

  Name of the subject identifier column (character scalar).

- group:

  Name of the within-subject condition column (character scalar).

- y:

  Name of the dependent variable column (character scalar).

- alpha:

  Family-wise alpha used for CI adjustment (default \`0.05\`).

- pAdjust:

  P-value adjustment method for multiple comparisons (default
  \`"holm"\`). Must be accepted by \`stats::p.adjust()\`.

- ciAdjust:

  CI adjustment method. One of \`"bonferroni"\` or \`"none"\`.

## Value

A tibble with one row per pairwise contrast and the columns below.

## Details

\- Difference estimate and CI: \`stats::wilcox.test()\` on paired
differences (Hodges-Lehmann location shift + CI). - Adjusted p-values:
\`stats::p.adjust()\` (default: Holm). - Effect size: rank-biserial r
via \`rstatix::wilcox_effsize()\` (paired).

## Output columns

- analysis_method:

  Analysis identifier string. Always \`"friedman_test"\`.

- effect_term:

  Effect label for schema alignment. Always \`"group"\`.

- comparison_kind:

  Comparison family label. Always \`"pairwise"\`.

- contrast_label:

  Human-readable label, e.g., \`"A vs B"\`.

- contrast_id:

  Machine-readable id, e.g., \`"group:A\|B"\`.

- group1:

  First condition label.

- group2:

  Second condition label.

- center_type:

  Summary location type for each condition. Always \`"median"\`.

- group1_center:

  Median of \`y\` in \`group1\`.

- group2_center:

  Median of \`y\` in \`group2\`.

- dispersion_type:

  Summary dispersion type for each condition. Always \`"iqr"\`.

- group1_dispersion:

  IQR of \`y\` in \`group1\`.

- group2_dispersion:

  IQR of \`y\` in \`group2\`.

- n_group1:

  Number of subjects with non-missing values in \`group1\`.

- n_group2:

  Number of subjects with non-missing values in \`group2\`.

- n_total:

  Number of complete paired observations used for the contrast.

- difference_estimate:

  Hodges-Lehmann location shift estimate (group1 - group2).

- difference_estimate_type:

  Difference estimator label. Always
  \`"hodges_lehmann_location_shift"\`.

- difference_direction:

  Difference direction label. Always \`"group1_minus_group2"\`.

- ci_low:

  Lower bound of the CI for \`difference_estimate\`.

- ci_high:

  Upper bound of the CI for \`difference_estimate\`.

- ci_level:

  Confidence level used for the CI (possibly Bonferroni-adjusted).

- ci_adjustment:

  CI adjustment label: \`"bonferroni_adjusted"\` or \`"unadjusted"\`.

- test_statistic:

  Wilcoxon signed-rank statistic value (V).

- test_statistic_type:

  Test statistic type. Always \`"wilcoxon_V"\`.

- df1:

  Reserved (always \`NA_real\_\` for this test family).

- df2:

  Reserved (always \`NA_real\_\` for this test family).

- p_value:

  Raw p-value from the paired Wilcoxon test.

- p_value_adjusted:

  Multiplicity-adjusted p-value (via \`pAdjust\`).

- p_value_adjustment:

  Adjustment label, e.g., \`"holm_adjusted"\`.

- effect_size:

  Effect size value (rank-biserial r).

- effect_size_type:

  Effect size type label. Always \`"rank_biserial_r"\`.

- effect_size_scale:

  Effect size scale label. Always \`"rank_based_r"\`.

- effect_size_ci_low:

  Reserved for effect size CI lower bound (currently \`NA_real\_\`).

- effect_size_ci_high:

  Reserved for effect size CI upper bound (currently \`NA_real\_\`).

## Examples

``` r
set.seed(22)
n = 20
u = rnorm(n, 0, 0.25)
w = data.frame(
  id = factor(seq_len(n)),
  A = rlnorm(n, meanlog = 3.7, sdlog = 0.25) * exp(u),
  B = rlnorm(n, meanlog = 3.9, sdlog = 0.25) * exp(u),
  C = rlnorm(n, meanlog = 4.0, sdlog = 0.25) * exp(u)
)
d = tidyr::pivot_longer(w, cols = c(A, B, C), names_to = "group", values_to = "y")
d$group = factor(d$group)
friedmanPost(d, id = "id", group = "group", y = "y")
#> Error in map(., .f, data, formula, method, ...): ℹ In index: 1.
#> Caused by error in `required_package()`:
#> ! coin package needed to be installed before using this function. Type this in R: install.packages('coin')
```
