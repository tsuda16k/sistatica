# Post hoc tests for the Kruskal-Wallis test (Dunn p-values + Wilcoxon HL CIs)

Performs pairwise post hoc comparisons after a Kruskal-Wallis test.

## Usage

``` r
kruskalPost(
  df,
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

- group:

  Name of the between-subject group column (character scalar).

- y:

  Name of the dependent variable column (character scalar).

- alpha:

  Family-wise alpha used for CI adjustment (default \`0.05\`).

- pAdjust:

  P-value adjustment method for multiple comparisons (default
  \`"holm"\`). Passed to \`rstatix::dunn_test()\` and recorded in
  output.

- ciAdjust:

  CI adjustment method. One of \`"bonferroni"\` or \`"none"\`.

## Value

A tibble with one row per pairwise contrast and the columns below.

## Details

\- Adjusted p-values: \`rstatix::dunn_test()\` (default: Holm). -
Difference estimate and CI: \`stats::wilcox.test()\` between groups
(Hodges-Lehmann location shift + CI). - Effect size: rank-biserial r via
\`rstatix::wilcox_effsize()\` (unpaired).

## Output columns

- analysis_method:

  Analysis identifier string. Always \`"kruskal_wallis_test"\`.

- effect_term:

  Effect label for schema alignment. Always \`"group"\`.

- comparison_kind:

  Comparison family label. Always \`"pairwise"\`.

- contrast_label:

  Human-readable label, e.g., \`"A vs B"\`.

- contrast_id:

  Machine-readable id, e.g., \`"group:A\|B"\`.

- group1:

  First group label.

- group2:

  Second group label.

- center_type:

  Summary location type for each group. Always \`"median"\`.

- group1_center:

  Median of \`y\` in \`group1\`.

- group2_center:

  Median of \`y\` in \`group2\`.

- dispersion_type:

  Summary dispersion type for each group. Always \`"iqr"\`.

- group1_dispersion:

  IQR of \`y\` in \`group1\`.

- group2_dispersion:

  IQR of \`y\` in \`group2\`.

- n_group1:

  Number of non-missing observations in \`group1\`.

- n_group2:

  Number of non-missing observations in \`group2\`.

- n_total:

  Total observations used for the contrast (\`n_group1 + n_group2\`).

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

  Dunn test statistic value (Z).

- test_statistic_type:

  Test statistic type. Always \`"dunn_Z"\`.

- df1:

  Reserved (always \`NA_real\_\` for this test family).

- df2:

  Reserved (always \`NA_real\_\` for this test family).

- p_value:

  Raw p-value from Dunn test.

- p_value_adjusted:

  Multiplicity-adjusted p-value from Dunn test.

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
set.seed(30)
d = data.frame(
  group = factor(rep(c("A","B","C"), each = 20)),
  y = c(rlnorm(20, 3.8, 0.3), rlnorm(20, 3.9, 0.35), rlnorm(20, 4.1, 0.4))
)
kruskalPost(d, group = "group", y = "y")
#> Error in map(., .f, data, formula, method, ...): ℹ In index: 1.
#> Caused by error in `required_package()`:
#> ! coin package needed to be installed before using this function. Type this in R: install.packages('coin')
```
