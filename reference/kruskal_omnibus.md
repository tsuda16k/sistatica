# Kruskal-Wallis test (omnibus) with eta-squared (H-based)

Performs the omnibus Kruskal-Wallis test using
[`rstatix::kruskal_test()`](https://rpkgs.datanovia.com/rstatix/reference/kruskal_test.html)
and computes an H-based eta-squared using
[`rstatix::kruskal_effsize()`](https://rpkgs.datanovia.com/rstatix/reference/kruskal_effsize.html).

## Usage

``` r
kruskal_omnibus(data, dv, group, missing = c("pairwise", "complete_case"))
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

- `design`: Design label (`"nonparametric_between"`).

- `method`: Omnibus method label (`"kruskal_wallis"`).

- `term`: Name of the grouping factor (the tested term).

- `statistic_type`: Type of omnibus test statistic (`"H"`).

- `statistic`: Omnibus test statistic value (H).

- `df1`: Always `NA` for this test (reserved column).

- `df2`: Always `NA` for this test (reserved column).

- `df`: Degrees of freedom for the Kruskal-Wallis test (k - 1).

- `p`: Omnibus p-value.

- `effect_size`: Effect size value (eta-squared based on H).

- `effect_size_type`: Effect size label (`"eta2_H"`).

- `missing`: Missing-data strategy used (`"pairwise"` or
  `"complete_case"`).

- `n_obs`: Number of rows (observations) used in the analysis after
  filtering for missingness.

- `n_id`: Always `NA` for between-subject designs in this package
  (reserved column).

## Details

\*\*Test\*\*: Kruskal-Wallis.  
\*\*Statistic\*\*: H with df = k - 1.  
\*\*Effect size\*\*: H-based eta-squared from
[`rstatix::kruskal_effsize()`](https://rpkgs.datanovia.com/rstatix/reference/kruskal_effsize.html).

**Missing-data handling (`missing`)**

- `"pairwise"` and `"complete_case"` behave identically here: rows with
  missing DV/group are removed.

## Examples

``` r
kruskal_omnibus(datasets::PlantGrowth, dv = weight, group = group)
#> # A tibble: 1 × 14
#>   design          method term  statistic_type statistic   df1   df2    df      p
#>   <chr>           <chr>  <chr> <chr>              <dbl> <dbl> <dbl> <dbl>  <dbl>
#> 1 nonparametric_… krusk… group H                   7.99    NA    NA     2 0.0184
#> # ℹ 5 more variables: effect_size <dbl>, effect_size_type <chr>, missing <chr>,
#> #   n_obs <dbl>, n_id <dbl>
```
