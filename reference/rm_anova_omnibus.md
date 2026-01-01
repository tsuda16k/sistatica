# Repeated-measures ANOVA (omnibus) with generalized eta-squared

Performs the omnibus test for a one-factor repeated-measures ANOVA
(within-subject design) using
[`rstatix::anova_test()`](https://rpkgs.datanovia.com/rstatix/reference/anova_test.html)
and returns a tidy one-row tibble with a standardized omnibus schema and
generalized eta-squared (GES).

## Usage

``` r
rm_anova_omnibus(
  data,
  dv,
  within,
  id,
  missing = c("pairwise", "complete_case")
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

- missing:

  Missing-data strategy: `"pairwise"` (default) or `"complete_case"`.

## Value

A one-row tibble (omnibus schema). Columns:

- `design`: Design label (`"parametric_within"`).

- `method`: Omnibus method label (`"rm_anova"`).

- `term`: Name of the within-subject factor (the tested term).

- `statistic_type`: Type of omnibus test statistic (`"F"`).

- `statistic`: Omnibus test statistic value (F).

- `df1`: Numerator degrees of freedom for F.

- `df2`: Denominator degrees of freedom for F.

- `df`: Single df column (always `NA` for this test; included for schema
  consistency).

- `p`: Omnibus p-value.

- `effect_size`: Effect size value (generalized eta-squared).

- `effect_size_type`: Effect size label (`"ges"`).

- `missing`: Missing-data strategy used (`"pairwise"` or
  `"complete_case"`).

- `n_obs`: Number of rows (observations) used in the analysis after
  filtering for missingness.

- `n_id`: Number of unique IDs used in the analysis after filtering for
  missingness.

## Details

\*\*Test\*\*: One-way repeated-measures ANOVA (within-subject).  
\*\*Statistic\*\*: F with numerator/denominator df.  
\*\*Effect size\*\*: Generalized eta-squared (GES) computed by
`rstatix::anova_test(effect.size = "ges")`.

**Missing-data handling (`missing`)**

- `"pairwise"` (default): Removes rows with missing DV/within/id
  (row-wise filtering).

- `"complete_case"`: After row-wise filtering, keeps only IDs that have
  non-missing DV values across all within-levels (complete blocks).

## Examples

``` r
dat = datasets::CO2 %>%
  dplyr::filter(Type == "Quebec", Treatment == "chilled") %>%
  dplyr::mutate(conc = factor(conc))

dat$uptake[dat$Plant == dat$Plant[1] & dat$conc == dat$conc[1]] = NA

rm_anova_omnibus(dat, dv = uptake, within = conc, id = Plant, missing = "pairwise")
#> # A tibble: 1 × 14
#>   design          method term  statistic_type statistic   df1   df2    df      p
#>   <chr>           <chr>  <chr> <chr>              <dbl> <dbl> <dbl> <dbl>  <dbl>
#> 1 parametric_wit… rm_an… conc  F                   24.3     6     6    NA 5.8e-4
#> # ℹ 5 more variables: effect_size <dbl>, effect_size_type <chr>, missing <chr>,
#> #   n_obs <dbl>, n_id <dbl>
rm_anova_omnibus(dat, dv = uptake, within = conc, id = Plant, missing = "complete_case")
#> # A tibble: 1 × 14
#>   design          method term  statistic_type statistic   df1   df2    df      p
#>   <chr>           <chr>  <chr> <chr>              <dbl> <dbl> <dbl> <dbl>  <dbl>
#> 1 parametric_wit… rm_an… conc  F                   24.3     6     6    NA 5.8e-4
#> # ℹ 5 more variables: effect_size <dbl>, effect_size_type <chr>, missing <chr>,
#> #   n_obs <dbl>, n_id <dbl>
```
