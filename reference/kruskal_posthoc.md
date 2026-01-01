# Kruskal-Wallis post-hoc comparisons (Dunn) with rank-biserial and CIs

Performs pairwise post-hoc comparisons using Dunn tests via
[`rstatix::dunn_test()`](https://rpkgs.datanovia.com/rstatix/reference/dunn_test.html),
and adds rank-biserial correlation (unpaired) with unadjusted and
adjusted confidence intervals.

## Usage

``` r
kruskal_posthoc(
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

- `design`: Design label (`"nonparametric_between"`).

- `method`: Post-hoc method label (`"dunn"`).

- `term`: Name of the grouping factor.

- `group1`, `group2`: Level names for the pairwise comparison.

- `statistic_type`: Type of test statistic (`"Z"`).

- `statistic`: Z statistic value (Dunn test).

- `df`: Always `NA` for this post-hoc method in this package.

- `p`: Unadjusted p-value.

- `p_adj`: Adjusted p-value (according to `p_adjust_method`).

- `p_adj_mark`: `"*"` if `p_adj < 0.05`, otherwise `""` (including
  `NA`).

- `p_adjust_method`: The p-value adjustment method used.

- `effect_size`: Effect size estimate (rank-biserial correlation;
  unpaired).

- `effect_size_type`: Effect size label (`"rank_biserial"`).

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

\*\*Post-hoc test\*\*: Dunn tests.  
\*\*Effect size\*\*: rank-biserial via
`effectsize::rank_biserial(paired = FALSE)`.  
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
kruskal_posthoc(
  datasets::PlantGrowth, dv = weight, group = group,
  p_adjust_method = "holm",
  ci_adjust_method = "bonferroni"
)
#> # A tibble: 3 × 23
#>   design       method term  group1 group2 statistic_type statistic    df       p
#>   <chr>        <chr>  <chr> <chr>  <chr>  <chr>              <dbl> <dbl>   <dbl>
#> 1 nonparametr… dunn   group ctrl   trt1   Z                  -1.12    NA 0.264  
#> 2 nonparametr… dunn   group ctrl   trt2   Z                   1.69    NA 0.0912 
#> 3 nonparametr… dunn   group trt1   trt2   Z                   2.81    NA 0.00500
#> # ℹ 14 more variables: p_adj <dbl>, p_adj_mark <chr>, p_adjust_method <chr>,
#> #   effect_size <dbl>, effect_size_type <chr>, ci_level <dbl>, ci_low <dbl>,
#> #   ci_high <dbl>, ci_adjust_method <chr>, ci_level_adj <dbl>,
#> #   ci_low_adj <dbl>, ci_high_adj <dbl>, n_pair <dbl>, missing <chr>
```
