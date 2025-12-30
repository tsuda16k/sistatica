# Omnibus test for the Friedman test (one-way repeated-measures, nonparametric)

Runs a one-way Friedman test using \`rstatix::friedman_test()\` and
returns a tibble in a standardized schema.

## Usage

``` r
friedmanOmni(df, id = "id", group = "group", y = "y")
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

## Value

A tibble with one row and the columns below.

## Output columns

- analysis_method:

  Analysis identifier string. Always \`"friedman_test"\`.

- effect_term:

  Effect label for schema alignment. Always \`"group"\`.

- test_statistic:

  Friedman omnibus test statistic value (chi-square).

- test_statistic_type:

  Test statistic type. Always \`"chi_square"\`.

- df1:

  Degrees of freedom.

- df2:

  Reserved (always \`NA_real\_\` for this test family).

- p_value:

  Omnibus p-value.

- effect_size:

  Effect size value (Kendall's W).

- effect_size_type:

  Effect size type label. Always \`"kendalls_W"\`.

- effect_size_scale:

  Effect size scale label. Always \`"agreement"\`.

- effect_size_ci_low:

  Reserved for effect size CI lower bound (currently \`NA_real\_\`).

- effect_size_ci_high:

  Reserved for effect size CI upper bound (currently \`NA_real\_\`).

## Examples

``` r
set.seed(21)
n = 24
u = rnorm(n, 0, 0.25)
w = data.frame(
  id = factor(seq_len(n)),
  A = rlnorm(n, meanlog = 3.8, sdlog = 0.25) * exp(u),
  B = rlnorm(n, meanlog = 3.9, sdlog = 0.25) * exp(u),
  C = rlnorm(n, meanlog = 4.0, sdlog = 0.25) * exp(u)
)
d = tidyr::pivot_longer(w, cols = c(A, B, C), names_to = "group", values_to = "y")
d$group = factor(d$group)
friedmanOmni(d, id = "id", group = "group", y = "y")
#> # A tibble: 1 × 12
#>   analysis_method effect_term test_statistic test_statistic_type   df1   df2
#>   <chr>           <chr>                <dbl> <chr>               <dbl> <dbl>
#> 1 friedman_test   group                 2.33 chi_square              2    NA
#> # ℹ 6 more variables: p_value <dbl>, effect_size <dbl>, effect_size_type <chr>,
#> #   effect_size_scale <chr>, effect_size_ci_low <dbl>,
#> #   effect_size_ci_high <dbl>
```
