# Welch-type post-hoc comparisons (pairwise Welch t) with Hedges' g and CIs

Performs pairwise Welch t-tests (`pool.sd = FALSE`) via
[`rstatix::pairwise_t_test()`](https://rpkgs.datanovia.com/rstatix/reference/t_test.html),
and adds Hedges' g (unpooled SD) with unadjusted and adjusted confidence
intervals.

## Usage

``` r
welch_posthoc(
  data,
  dv,
  group,
  conf_level = 0.95,
  p_adjust_method = "holm",
  ci_adjust_method = c("bonferroni", "none", "sidak"),
  missing = c("pairwise", "complete_case")
)
```

## Arguments

- data:

  A data.frame.

- dv:

  Dependent variable (numeric). Unquoted name or string.

- group:

  Grouping factor. Unquoted name or string.

- conf_level:

  Confidence level for unadjusted effect-size CIs (default: 0.95).

- p_adjust_method:

  P-value adjustment method passed to `rstatix` (default: `"holm"`).

- ci_adjust_method:

  CI multiplicity adjustment method (default: `"bonferroni"`).

- missing:

  Missing-data strategy: `"pairwise"` (default) or `"complete_case"`.

## Value

A tibble with one row per pairwise comparison. Columns:

- `design`: Design label (`"parametric_between_unequal_var"`).

- `method`: Post-hoc method label (`"welch_pairwise_t"`).

- `term`: Name of the grouping factor.

- `group1`, `group2`: Level names for the pairwise comparison.

- `statistic_type`: Type of test statistic (`"t"`).

- `statistic`: Test statistic value.

- `df`: Degrees of freedom (Welch t-test).

- `p`: Unadjusted p-value.

- `p_adj`: Adjusted p-value (according to `p_adjust_method`).

- `p_adj_mark`: `"*"` if `p_adj < 0.05`, otherwise `""` (including
  `NA`).

- `p_adjust_method`: The p-value adjustment method used.

- `effect_size`: Effect size estimate (Hedges' g, unpooled SD).

- `effect_size_type`: Effect size label (`"hedges_g"`).

- `ci_level`: Confidence level used for the unadjusted effect-size CI.

- `ci_low`, `ci_high`: Lower/upper bounds of the unadjusted effect-size
  CI.

- `ci_adjust_method`: Method used to adjust effect-size CIs.

- `ci_level_adj`: Adjusted CI level implied by `ci_adjust_method` and
  the number of comparisons.

- `ci_low_adj`, `ci_high_adj`: Lower/upper bounds of the adjusted
  effect-size CI.

- `n_pair`: Total sample size used for effect size/CIs (`n_x + n_y`).

- `missing`: Missing-data strategy used.

## Details

\*\*Post-hoc test\*\*: Welch t-tests (unpooled SD).  
\*\*Effect size\*\*: Hedges' g via
`effectsize::hedges_g(pooled_sd = FALSE)`.  
\*\*P-value adjustment\*\*: `p_adjust_method` (default: `"holm"`).  
\*\*CI adjustment\*\*: controlled by `ci_adjust_method` (independent of
p-value adjustment); supported: `"none"`, `"bonferroni"`, `"sidak"`.

**Important note on independence of p-adjustment and CI-adjustment**:
See
[`rm_posthoc`](https://tsuda16k.github.io/sistatica/reference/rm_posthoc.md)
for a detailed explanation. The same note applies here.

**Missing-data handling (`missing`)**

- `"pairwise"` and `"complete_case"` behave identically here: rows with
  missing DV/group are removed.

## Examples

``` r
welch_posthoc(
  datasets::PlantGrowth, dv = weight, group = group,
  p_adjust_method = "holm",
  ci_adjust_method = "bonferroni"
)
#> # A tibble: 3 × 23
#>   design   method term  group1 group2 statistic_type statistic    df     p p_adj
#>   <chr>    <chr>  <chr> <chr>  <chr>  <chr>              <dbl> <dbl> <dbl> <dbl>
#> 1 paramet… welch… group ctrl   trt1   t                   1.19  16.5 0.25  0.25 
#> 2 paramet… welch… group ctrl   trt2   t                  -2.13  16.8 0.048 0.096
#> 3 paramet… welch… group trt1   trt2   t                  -3.01  14.1 0.009 0.028
#> # ℹ 13 more variables: p_adj_mark <chr>, p_adjust_method <chr>,
#> #   effect_size <dbl>, effect_size_type <chr>, ci_level <dbl>, ci_low <dbl>,
#> #   ci_high <dbl>, ci_adjust_method <chr>, ci_level_adj <dbl>,
#> #   ci_low_adj <dbl>, ci_high_adj <dbl>, n_pair <dbl>, missing <chr>
```
