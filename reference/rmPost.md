# Post hoc tests for repeated-measures ANOVA (paired t-tests)

Performs pairwise paired t-tests for a one-way repeated-measures design
and returns a tibble in a standardized (method-agnostic) post-hoc
schema.

## Usage

``` r
rmPost(
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

\- Adjusted p-values: via \`stats::p.adjust()\` (default: Holm). -
Confidence intervals: unadjusted or Bonferroni-adjusted confidence
level. - Effect size: Cohen's d for paired samples (dz in practice) via
\`rstatix::cohens_d()\`.

## Output columns

- analysis_method:

  Analysis identifier string. Always \`"repeated_measures_anova"\`.

- effect_term:

  Effect label for schema alignment. Always \`"group"\`.

- comparison_kind:

  Comparison family label. Always \`"pairwise"\`.

- contrast_label:

  Human-readable label, e.g., \`"A vs B"\`.

- contrast_id:

  Machine-readable id, e.g., \`"group:A\|B"\`.

- group1:

  First condition label (maps to the first element of a pair).

- group2:

  Second condition label (maps to the second element of a pair).

- center_type:

  Summary location type for each condition. Always \`"mean"\`.

- group1_center:

  Mean of \`y\` in \`group1\`.

- group2_center:

  Mean of \`y\` in \`group2\`.

- dispersion_type:

  Summary dispersion type for each condition. Always \`"sd"\`.

- group1_dispersion:

  Standard deviation of \`y\` in \`group1\`.

- group2_dispersion:

  Standard deviation of \`y\` in \`group2\`.

- n_group1:

  Number of subjects with non-missing values in \`group1\`.

- n_group2:

  Number of subjects with non-missing values in \`group2\`.

- n_total:

  Number of complete paired observations used for the contrast.

- difference_estimate:

  Estimated mean difference (group1 - group2).

- difference_estimate_type:

  Difference estimator label. Always \`"mean_difference"\`.

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

  Paired t-test statistic value (t).

- test_statistic_type:

  Test statistic type. Always \`"t"\`.

- df1:

  Degrees of freedom for the t-test.

- df2:

  Reserved (always \`NA_real\_\` for this test family).

- p_value:

  Raw (unadjusted) p-value for the pairwise t-test.

- p_value_adjusted:

  Multiplicity-adjusted p-value (via \`pAdjust\`).

- p_value_adjustment:

  Adjustment label, e.g., \`"holm_adjusted"\`.

- effect_size:

  Effect size value (Cohen's dz).

- effect_size_type:

  Effect size type label. Always \`"cohens_dz"\`.

- effect_size_scale:

  Effect size scale label. Always \`"standardized_mean_difference"\`.

- effect_size_ci_low:

  Reserved for effect size CI lower bound (currently \`NA_real\_\`).

- effect_size_ci_high:

  Reserved for effect size CI upper bound (currently \`NA_real\_\`).

## Examples

``` r
set.seed(2)
n = 18
u = rnorm(n, 0, 3)
w = data.frame(
  id = factor(seq_len(n)),
  u  = u,
  A  = 10 + u + rnorm(n, 0, 2),
  B  = 11 + u + rnorm(n, 0, 2),
  C  = 13 + u + rnorm(n, 0, 2)
)
d = tidyr::pivot_longer(w, cols = c(A, B, C), names_to = "group", values_to = "y")
d$group = factor(d$group)
rmPost(d, id = "id", group = "group", y = "y")
#> # A tibble: 3 × 35
#>   analysis_method  effect_term comparison_kind contrast_label contrast_id group1
#>   <chr>            <chr>       <chr>           <chr>          <chr>       <chr> 
#> 1 repeated_measur… group       pairwise        A vs B         group:A|B   A     
#> 2 repeated_measur… group       pairwise        A vs C         group:A|C   A     
#> 3 repeated_measur… group       pairwise        B vs C         group:B|C   B     
#> # ℹ 29 more variables: group2 <chr>, center_type <chr>, group1_center <dbl>,
#> #   group2_center <dbl>, dispersion_type <chr>, group1_dispersion <dbl>,
#> #   group2_dispersion <dbl>, n_group1 <int>, n_group2 <int>, n_total <int>,
#> #   difference_estimate <dbl>, difference_estimate_type <chr>,
#> #   difference_direction <chr>, ci_low <dbl>, ci_high <dbl>, ci_level <dbl>,
#> #   ci_adjustment <chr>, test_statistic <dbl>, test_statistic_type <chr>,
#> #   df1 <dbl>, df2 <dbl>, p_value <dbl>, p_value_adjusted <dbl>, …
```
