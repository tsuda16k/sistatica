# Repeated-measures post-hoc comparisons (paired t) with dz and CIs

Performs pairwise post-hoc comparisons for a one-factor
repeated-measures design using paired t-tests
([`rstatix::pairwise_t_test()`](https://rpkgs.datanovia.com/rstatix/reference/t_test.html)),
and adds Cohen's \\d_z\\ with both unadjusted and multiplicity-adjusted
confidence intervals.

## Usage

``` r
rm_posthoc(
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

- `design`: Design label (`"parametric_within"`).

- `method`: Post-hoc method label (`"paired_t"`).

- `term`: Name of the within-factor.

- `group1`, `group2`: Level names for the pairwise comparison.

- `statistic_type`: Type of test statistic (`"t"`).

- `statistic`: Test statistic value.

- `df`: Degrees of freedom (t-test).

- `p`: Unadjusted p-value for the pairwise test.

- `p_adj`: Adjusted p-value (according to `p_adjust_method`).

- `p_adj_mark`: Visual marker for significance based on `p_adj`: `"*"`
  if `p_adj < 0.05`, otherwise `""` (including when `p_adj` is `NA`).

- `p_adjust_method`: The p-value adjustment method used.

- `effect_size`: Effect size estimate (\\d_z\\).

- `effect_size_type`: Effect size label (`"dz"`).

- `ci_level`: Confidence level used for the unadjusted effect-size CI.

- `ci_low`, `ci_high`: Lower/upper bounds of the unadjusted effect-size
  CI.

- `ci_adjust_method`: Method used to adjust effect-size CIs (`"none"`,
  `"bonferroni"`, `"sidak"`).

- `ci_level_adj`: Adjusted CI level implied by `ci_adjust_method` and
  the number of comparisons.

- `ci_low_adj`, `ci_high_adj`: Lower/upper bounds of the adjusted
  effect-size CI.

- `n_pair`: Number of paired observations used for effect size/CIs
  (complete pairs for the two levels).

- `missing`: Missing-data strategy used.

## Details

\*\*Post-hoc test\*\*: paired t-tests.  
\*\*Effect size\*\*: \\d_z\\ via
`effectsize::repeated_measures_d(method = "z")`.  
\*\*P-value adjustment\*\*: `p_adjust_method` (default: `"holm"`).  
\*\*CI adjustment\*\*: controlled by `ci_adjust_method` (independent of
p-value adjustment); supported: `"none"`, `"bonferroni"`, `"sidak"`.

**Important note on independence of p-adjustment and CI-adjustment**:
P-value adjustment (`p_adjust_method`) and CI adjustment
(`ci_adjust_method`) are independent. Therefore, it is possible to
obtain results where `p_adj < 0.05` (significant after p-adjustment)
while the adjusted CI (`ci_low_adj`, `ci_high_adj`) still includes 0 (or
vice versa), especially when using different adjustment strategies
(e.g., `p_adjust_method="holm"` and `ci_adjust_method="none"`). This is
not a software bug but a consequence of using different error-control
strategies for p-values and confidence intervals.

**Missing-data handling (`missing`)**

- `"pairwise"` (default):

  - P-values are computed after removing rows with missing DV/within/id.

  - Effect sizes and CIs are computed pairwise: for each comparison,
    only IDs with non-missing values in the two relevant within-levels
    are used.

- `"complete_case"`: Restricts to IDs complete across all within-levels
  before analysis.

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
  dplyr::filter(Type == "Quebec", Treatment == "chilled") %>%
  dplyr::mutate(conc = factor(conc))

rm_posthoc(
  dat, dv = uptake, within = conc, id = Plant,
  p_adjust_method = "holm",
  ci_adjust_method = "bonferroni",
  missing = "pairwise"
)
#> # A tibble: 21 × 23
#>    design  method term  group1 group2 statistic_type statistic    df     p p_adj
#>    <chr>   <chr>  <chr> <chr>  <chr>  <chr>              <dbl> <dbl> <dbl> <dbl>
#>  1 parame… paire… conc  95     175    t                  -3.17     2 0.087 0.696
#>  2 parame… paire… conc  95     250    t                  -7.56     2 0.017 0.291
#>  3 parame… paire… conc  95     350    t                  -6.92     2 0.02  0.303
#>  4 parame… paire… conc  95     500    t                  -7.50     2 0.017 0.291
#>  5 parame… paire… conc  95     675    t                 -12.2      2 0.007 0.133
#>  6 parame… paire… conc  95     1000   t                 -10.7      2 0.009 0.164
#>  7 parame… paire… conc  175    250    t                  -3.03     2 0.094 0.696
#>  8 parame… paire… conc  175    350    t                 -16.1      2 0.004 0.081
#>  9 parame… paire… conc  175    500    t                  -4.46     2 0.047 0.608
#> 10 parame… paire… conc  175    675    t                  -5.07     2 0.037 0.515
#> # ℹ 11 more rows
#> # ℹ 13 more variables: p_adj_mark <chr>, p_adjust_method <chr>,
#> #   effect_size <dbl>, effect_size_type <chr>, ci_level <dbl>, ci_low <dbl>,
#> #   ci_high <dbl>, ci_adjust_method <chr>, ci_level_adj <dbl>,
#> #   ci_low_adj <dbl>, ci_high_adj <dbl>, n_pair <dbl>, missing <chr>
```
