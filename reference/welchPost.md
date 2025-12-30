# Post hoc tests for Welch one-way ANOVA (Games-Howell p-values + Welch t-test CIs)

Performs pairwise comparisons appropriate for heteroscedastic one-way
designs.

## Usage

``` r
welchPost(
  df,
  group = "group",
  y = "y",
  alpha = 0.05,
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

- ciAdjust:

  CI adjustment method. One of \`"bonferroni"\` or \`"none"\`.

## Value

A tibble with one row per pairwise contrast and the columns below.

## Details

\- Adjusted p-values: \`rstatix::games_howell_test()\` (Games-Howell). -
Difference estimate and CI: \`stats::t.test(var.equal = FALSE)\` for
each pair. - Effect size: Hedges' g via
\`rstatix::cohens_d(hedges.correction = TRUE)\`.

## Output columns

- analysis_method:

  Analysis identifier string. Always \`"welch_anova"\`.

- effect_term:

  Effect label for schema alignment. Always \`"group"\`.

- comparison_kind:

  Comparison family label. Always \`"pairwise"\`.

- contrast_label:

  Human-readable label, e.g., \`"G1 vs G2"\`.

- contrast_id:

  Machine-readable id, e.g., \`"group:G1\|G2"\`.

- group1:

  First group label (maps to the first element of a pair).

- group2:

  Second group label (maps to the second element of a pair).

- center_type:

  Summary location type for each group. Always \`"mean"\`.

- group1_center:

  Mean of \`y\` in \`group1\`.

- group2_center:

  Mean of \`y\` in \`group2\`.

- dispersion_type:

  Summary dispersion type for each group. Always \`"sd"\`.

- group1_dispersion:

  Standard deviation of \`y\` in \`group1\`.

- group2_dispersion:

  Standard deviation of \`y\` in \`group2\`.

- n_group1:

  Number of non-missing observations in \`group1\`.

- n_group2:

  Number of non-missing observations in \`group2\`.

- n_total:

  Total observations used for the contrast (\`n_group1 + n_group2\`).

- difference_estimate:

  Estimated mean difference (group1 - group2).

- difference_estimate_type:

  Difference estimator label. Always \`"mean_difference"\`.

- difference_direction:

  Difference direction label. Always \`"group1_minus_group2"\`.

- ci_low:

  Lower bound of the CI for \`difference_estimate\` (Welch t-test).

- ci_high:

  Upper bound of the CI for \`difference_estimate\` (Welch t-test).

- ci_level:

  Confidence level used for the CI (possibly Bonferroni-adjusted).

- ci_adjustment:

  CI adjustment label: \`"bonferroni_adjusted"\` or \`"unadjusted"\`.

- test_statistic:

  Welch t statistic value for the pair.

- test_statistic_type:

  Test statistic type. Always \`"t"\`.

- df1:

  Approximate df from the Welch t-test for the pair.

- df2:

  Reserved (always \`NA_real\_\` for this test family).

- p_value:

  Raw p-value from Welch t-test for the pair.

- p_value_adjusted:

  Games-Howell adjusted p-value for the pair.

- p_value_adjustment:

  Adjustment label. Always \`"games_howell_adjusted"\`.

- effect_size:

  Effect size value (Hedges' g).

- effect_size_type:

  Effect size type label. Always \`"hedges_g"\`.

- effect_size_scale:

  Effect size scale label. Always \`"standardized_mean_difference"\`.

- effect_size_ci_low:

  Reserved for effect size CI lower bound (currently \`NA_real\_\`).

- effect_size_ci_high:

  Reserved for effect size CI upper bound (currently \`NA_real\_\`).

## Examples

``` r
set.seed(11)
d = data.frame(
  group = factor(rep(c("G1","G2","G3","G4"), each = 20)),
  y = c(rnorm(20, 0, 1), rnorm(20, 0.3, 2), rnorm(20, 1.0, 1.2), rnorm(20, 0.6, 3))
)
welchPost(d, group = "group", y = "y")
#> # A tibble: 6 × 35
#>   analysis_method effect_term comparison_kind contrast_label contrast_id group1
#>   <chr>           <chr>       <chr>           <chr>          <chr>       <chr> 
#> 1 welch_anova     group       pairwise        G1 vs G2       group:G1|G2 G1    
#> 2 welch_anova     group       pairwise        G1 vs G3       group:G1|G3 G1    
#> 3 welch_anova     group       pairwise        G1 vs G4       group:G1|G4 G1    
#> 4 welch_anova     group       pairwise        G2 vs G3       group:G2|G3 G2    
#> 5 welch_anova     group       pairwise        G2 vs G4       group:G2|G4 G2    
#> 6 welch_anova     group       pairwise        G3 vs G4       group:G3|G4 G3    
#> # ℹ 29 more variables: group2 <chr>, center_type <chr>, group1_center <dbl>,
#> #   group2_center <dbl>, dispersion_type <chr>, group1_dispersion <dbl>,
#> #   group2_dispersion <dbl>, n_group1 <int>, n_group2 <int>, n_total <int>,
#> #   difference_estimate <dbl>, difference_estimate_type <chr>,
#> #   difference_direction <chr>, ci_low <dbl>, ci_high <dbl>, ci_level <dbl>,
#> #   ci_adjustment <chr>, test_statistic <dbl>, test_statistic_type <chr>,
#> #   df1 <dbl>, df2 <dbl>, p_value <dbl>, p_value_adjusted <dbl>, …
```
