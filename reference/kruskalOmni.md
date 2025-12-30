# Omnibus test for the Kruskal-Wallis test (one-way, nonparametric)

Runs a one-way Kruskal-Wallis test using \`rstatix::kruskal_test()\` and
returns a tibble in a standardized schema.

## Usage

``` r
kruskalOmni(df, group = "group", y = "y")
```

## Arguments

- df:

  A long-format data frame.

- group:

  Name of the between-subject group column (character scalar).

- y:

  Name of the dependent variable column (character scalar).

## Value

A tibble with one row and the columns below.

## Output columns

- analysis_method:

  Analysis identifier string. Always \`"kruskal_wallis_test"\`.

- effect_term:

  Effect label for schema alignment. Always \`"group"\`.

- test_statistic:

  Kruskal-Wallis omnibus test statistic value (chi-square).

- test_statistic_type:

  Test statistic type. Always \`"chi_square"\`.

- df1:

  Degrees of freedom.

- df2:

  Reserved (always \`NA_real\_\` for this test family).

- p_value:

  Omnibus p-value.

- effect_size:

  Effect size value (epsilon-squared).

- effect_size_type:

  Effect size type label. Always \`"epsilon_squared"\`.

- effect_size_scale:

  Effect size scale label. Always \`"rank_based_variance_explained"\`.

- effect_size_ci_low:

  Reserved for effect size CI lower bound (currently \`NA_real\_\`).

- effect_size_ci_high:

  Reserved for effect size CI upper bound (currently \`NA_real\_\`).

## Examples

``` r
data("PlantGrowth")
d = PlantGrowth
d$group = factor(d$group)
names(d) = c("y","group")
kruskalOmni(d, group = "group", y = "y")
#> # A tibble: 1 × 12
#>   analysis_method     effect_term test_statistic test_statistic_type   df1   df2
#>   <chr>               <chr>                <dbl> <chr>               <dbl> <dbl>
#> 1 kruskal_wallis_test group                 7.99 chi_square              2    NA
#> # ℹ 6 more variables: p_value <dbl>, effect_size <dbl>, effect_size_type <chr>,
#> #   effect_size_scale <chr>, effect_size_ci_low <dbl>,
#> #   effect_size_ci_high <dbl>
```
