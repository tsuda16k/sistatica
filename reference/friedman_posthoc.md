# Friedman post-hoc comparisons (paired Wilcoxon) with rank-biserial and CIs

Performs pairwise post-hoc comparisons for a one-factor within-subject
design using paired Wilcoxon signed-rank tests via
[`rstatix::pairwise_wilcox_test()`](https://rpkgs.datanovia.com/rstatix/reference/wilcox_test.html).
Also computes rank-biserial correlation (paired) with unadjusted and
adjusted CIs.

## Usage

``` r
friedman_posthoc(
  data,
  dv,
  within,
  id,
  conf_level = 0.95,
  p_adjust_method = "holm",
  ci_adjust_method = c("bonferroni", "none", "sidak"),
  missing = c("pairwise", "complete_case"),
  duplicate = c("error", "mean", "median", "first"),
  duplicate_na_rm = TRUE
)
```

## Arguments

- data:

  A data.frame in long format.

- dv:

  Dependent variable (numeric). Unquoted name or string.

- within:

  Within-subject factor. Unquoted name or string.

- id:

  Subject identifier. Unquoted name or string.

- conf_level:

  Confidence level for unadjusted effect-size CIs (default: 0.95).

- p_adjust_method:

  P-value adjustment method passed to `rstatix` (default: `"holm"`).

- ci_adjust_method:

  CI multiplicity adjustment method (default: `"bonferroni"`).

- missing:

  Missing-data strategy: `"pairwise"` (default) or `"complete_case"`.

- duplicate:

  How to handle duplicated ID x within cells: `"error"` (default),
  `"mean"`, `"median"`, `"first"`.

- duplicate_na_rm:

  Logical; if TRUE, removes NAs within duplicated cells before
  aggregation (default: TRUE).

## Value

A tibble with one row per pairwise comparison. Columns:

- `design`: Design label (`"nonparametric_within"`).

- `method`: Post-hoc method label (`"paired_wilcox"`).

- `term`: Name of the within-factor.

- `group1`, `group2`: Level names for the pairwise comparison.

- `statistic_type`: Type of test statistic (`"W"`).

- `statistic`: Wilcoxon test statistic (may be `NA` depending on rstatix
  output/version).

- `df`: Always `NA` for this post-hoc method in this package.

- `p`: Unadjusted p-value.

- `p_adj`: Adjusted p-value (according to `p_adjust_method`).

- `p_adj_mark`: `"*"` if `p_adj < 0.05`, otherwise `""` (including
  `NA`).

- `p_adjust_method`: The p-value adjustment method used.

- `effect_size`: Effect size estimate (rank-biserial correlation;
  paired).

- `effect_size_type`: Effect size label (`"rank_biserial"`).

- `ci_level`: Confidence level used for the unadjusted effect-size CI.

- `ci_low`, `ci_high`: Lower/upper bounds of the unadjusted effect-size
  CI.

- `ci_adjust_method`: Method used to adjust effect-size CIs.

- `ci_level_adj`: Adjusted CI level implied by `ci_adjust_method` and
  the number of comparisons.

- `ci_low_adj`, `ci_high_adj`: Lower/upper bounds of the adjusted
  effect-size CI.

- `n_pair`: Number of paired observations used for effect size/CIs
  (complete pairs for the two levels).

- `missing`: Missing-data strategy used.

## Details

\*\*Post-hoc test\*\*: paired Wilcoxon signed-rank tests.  
\*\*Effect size\*\*: rank-biserial via
`effectsize::rank_biserial(paired = TRUE)`.  
\*\*P-value adjustment\*\*: `p_adjust_method` (default: `"holm"`).  
\*\*CI adjustment\*\*: controlled by `ci_adjust_method` (independent of
p-value adjustment); supported: `"none"`, `"bonferroni"`, `"sidak"`.

**Important note on independence of p-adjustment and CI-adjustment**:
See
[`rm_posthoc`](https://tsuda16k.github.io/sistatica/reference/rm_posthoc.md)
for a detailed explanation. The same note applies here.

**Missing-data handling (`missing`)**

- `"pairwise"` (default): p-values are computed after removing rows with
  missing DV/within/id; effect sizes/CIs use only IDs complete for the
  two compared within-levels.

- `"complete_case"`: Restricts to IDs complete across all within-levels.

**Duplicated ID x within cells** If duplicated ID x within cells exist,
set `duplicate` to aggregate values within each cell. Aggregation is
reported via a single warning (not via return columns).

**Meaning of `duplicate_na_rm` for `duplicate="first"`**: If
`duplicate="first"` and `duplicate_na_rm=TRUE`, NAs within a duplicated
cell are removed before taking the first value (i.e., the first
non-missing value is used).

## Examples

``` r
dat = datasets::CO2 %>%
  dplyr::mutate(
    Plant = factor(Plant),
    conc = factor(conc)
  )

friedman_posthoc(
  dat,
  dv = uptake,
  within = conc,
  id = Plant,
  p_adjust_method = "holm",
  ci_adjust_method = "bonferroni",
  missing = "complete_case"
)
#> # A tibble: 21 × 23
#>    design      method term  group1 group2 statistic_type statistic    df       p
#>    <chr>       <chr>  <chr> <chr>  <chr>  <chr>              <dbl> <dbl>   <dbl>
#>  1 nonparamet… paire… conc  95     175    W                      0    NA 4.88e-4
#>  2 nonparamet… paire… conc  95     250    W                      0    NA 4.88e-4
#>  3 nonparamet… paire… conc  95     350    W                      0    NA 4.88e-4
#>  4 nonparamet… paire… conc  95     500    W                      0    NA 4.88e-4
#>  5 nonparamet… paire… conc  95     675    W                      0    NA 4.88e-4
#>  6 nonparamet… paire… conc  95     1000   W                      0    NA 4.88e-4
#>  7 nonparamet… paire… conc  175    250    W                      1    NA 9.77e-4
#>  8 nonparamet… paire… conc  175    350    W                      1    NA 9.77e-4
#>  9 nonparamet… paire… conc  175    500    W                      1    NA 9.77e-4
#> 10 nonparamet… paire… conc  175    675    W                      0    NA 4.88e-4
#> # ℹ 11 more rows
#> # ℹ 14 more variables: p_adj <dbl>, p_adj_mark <chr>, p_adjust_method <chr>,
#> #   effect_size <dbl>, effect_size_type <chr>, ci_level <dbl>, ci_low <dbl>,
#> #   ci_high <dbl>, ci_adjust_method <chr>, ci_level_adj <dbl>,
#> #   ci_low_adj <dbl>, ci_high_adj <dbl>, n_pair <dbl>, missing <chr>
```
