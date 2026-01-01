# Welch one-way ANOVA (omnibus) with omega-squared

Performs the omnibus Welch one-way ANOVA for a between-subject factor
using
[`rstatix::welch_anova_test()`](https://rpkgs.datanovia.com/rstatix/reference/welch_anova_test.html),
and computes omega-squared from the resulting F statistic and degrees of
freedom using
[`effectsize::F_to_omega2()`](https://easystats.github.io/effectsize/reference/F_to_eta2.html).

## Usage

``` r
welch_anova_omnibus(data, dv, group, missing = c("pairwise", "complete_case"))
```

## Arguments

- data:

  A data.frame.

- dv:

  Dependent variable (numeric). Unquoted name or string.

- group:

  Grouping factor. Unquoted name or string.

- missing:

  Missing-data strategy: `"pairwise"` (default) or `"complete_case"`.

## Value

A one-row tibble (omnibus schema). Columns:

- `design`: Design label (`"parametric_between_unequal_var"`).

- `method`: Omnibus method label (`"welch_anova"`).

- `term`: Name of the grouping factor (the tested term).

- `statistic_type`: Type of omnibus test statistic (`"F"`).

- `statistic`: Omnibus test statistic value (F; Welch test).

- `df1`: Numerator degrees of freedom for F (Welch).

- `df2`: Denominator degrees of freedom for F (Welch).

- `df`: Single df column (always `NA` for this test; included for schema
  consistency).

- `p`: Omnibus p-value.

- `effect_size`: Effect size value (omega-squared; approximation from
  Welch F and df).

- `effect_size_type`: Effect size label (`"omega2"`).

- `missing`: Missing-data strategy used (`"pairwise"` or
  `"complete_case"`).

- `n_obs`: Number of rows (observations) used in the analysis after
  filtering for missingness.

- `n_id`: Always `NA` for between-subject designs in this package
  (reserved column).

## Details

\*\*Test\*\*: Welch one-way ANOVA.  
\*\*Statistic\*\*: F with Welch-approximated df1/df2.  
\*\*Effect size\*\*: Omega-squared computed from F and df via
[`effectsize::F_to_omega2()`](https://easystats.github.io/effectsize/reference/F_to_eta2.html).

**Note**: The omega-squared value is computed from the Welch F statistic
and associated degrees of freedom. This provides a convenient
effect-size summary, but should be interpreted as an approximation.

**Missing-data handling (`missing`)**

- `"pairwise"` and `"complete_case"` behave identically here: rows with
  missing DV/group are removed.

## Examples

``` r
welch_anova_omnibus(datasets::PlantGrowth, dv = weight, group = group)
#> # A tibble: 1 × 14
#>   design           method term  statistic_type statistic   df1   df2    df     p
#>   <chr>            <chr>  <chr> <chr>              <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1 parametric_betw… welch… group F                   5.18     2  17.1    NA 0.017
#> # ℹ 5 more variables: effect_size <dbl>, effect_size_type <chr>, missing <chr>,
#> #   n_obs <dbl>, n_id <dbl>
```
