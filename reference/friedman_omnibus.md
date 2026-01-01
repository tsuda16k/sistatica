# Friedman test (omnibus) with Kendall's W

Performs the omnibus Friedman test for a one-factor within-subject
design and returns Kendall's W.

## Usage

``` r
friedman_omnibus(
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

- `design`: Design label (`"nonparametric_within"`).

- `method`: Omnibus method label (`"friedman"`).

- `term`: Name of the within-subject factor (the tested term).

- `statistic_type`: Type of omnibus test statistic (`"chi_squared"`).

- `statistic`: Omnibus test statistic value (\\\chi^2\\).

- `df1`: Always `NA` for this test (reserved column).

- `df2`: Always `NA` for this test (reserved column).

- `df`: Degrees of freedom for the Friedman test (k - 1).

- `p`: Omnibus p-value.

- `effect_size`: Effect size value (Kendall's W).

- `effect_size_type`: Effect size label (`"kendalls_W"`).

- `missing`: Missing-data strategy used (`"pairwise"` or
  `"complete_case"`).

- `n_obs`: Number of rows (observations) used after row-wise filtering
  (before enforcing complete blocks).

- `n_id`: Number of complete IDs (blocks) used in the Friedman test
  after enforcing complete cases.

## Details

\*\*Test\*\*: Friedman test.  
\*\*Requirement\*\*: unreplicated complete block design.  
\*\*Statistic\*\*: Chi-squared with df = k - 1.  
\*\*Effect size\*\*: Kendall's W computed as \\W = \chi^2 / (n (k -
1))\\.

**Implementation note**: The function stops if duplicated ID x within
cells are detected. After that check, the wide ID x within matrix is
built using `tapply(...)` with `mean` as the aggregation function.
Because duplicates are ruled out, each cell is effectively a single
value; using `mean` here is a simple and stable way to map long data to
a wide matrix without introducing extra dependencies or complicated
reshaping logic.

**Missing-data handling (`missing`)**

- `"pairwise"` and `"complete_case"` both enforce complete blocks for
  the omnibus test by keeping only complete rows in the wide matrix.

## Examples

``` r
dat = datasets::CO2 %>%
  dplyr::filter(Type == "Quebec", Treatment == "chilled") %>%
  dplyr::mutate(conc = factor(conc))

friedman_omnibus(dat, dv = uptake, within = conc, id = Plant)
#> # A tibble: 1 × 14
#>   design          method term  statistic_type statistic   df1   df2    df      p
#>   <chr>           <chr>  <chr> <chr>              <dbl> <dbl> <dbl> <dbl>  <dbl>
#> 1 nonparametric_… fried… conc  chi_squared         16.1    NA    NA     6 0.0130
#> # ℹ 5 more variables: effect_size <dbl>, effect_size_type <chr>, missing <chr>,
#> #   n_obs <dbl>, n_id <dbl>
```
