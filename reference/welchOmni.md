# Omnibus test for Welch one-way ANOVA

Runs a one-way Welch ANOVA (heteroscedastic ANOVA) using
\`rstatix::welch_anova_test()\` and returns a tibble in a standardized
schema.

## Usage

``` r
welchOmni(df, group = "group", y = "y")
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

## Details

The omnibus effect size is computed as an approximation: \$\$\eta^2
\approx \frac{df_1 F}{df_1 F + df_2}\$\$ where \\F\\ is the Welch
statistic and \\df_1, df_2\\ are the approximate degrees of freedom
returned by the Welch test.

## Output columns

- analysis_method:

  Analysis identifier string. Always \`"welch_anova"\`.

- effect_term:

  Effect label for schema alignment. Always \`"group"\`.

- test_statistic:

  Welch omnibus test statistic value (F).

- test_statistic_type:

  Test statistic type. Always \`"F"\`.

- df1:

  Approximate numerator degrees of freedom.

- df2:

  Approximate denominator degrees of freedom.

- p_value:

  Omnibus p-value.

- effect_size:

  Approximate eta-squared computed from F and dfs.

- effect_size_type:

  Effect size type label. Always \`"eta_squared_approx"\`.

- effect_size_scale:

  Effect size scale label. Always \`"variance_explained"\`.

- effect_size_ci_low:

  Reserved for effect size CI lower bound (currently \`NA_real\_\`).

- effect_size_ci_high:

  Reserved for effect size CI upper bound (currently \`NA_real\_\`).

## Examples

``` r
set.seed(10)
d = data.frame(
  group = factor(rep(c("G1","G2","G3"), each = 25)),
  y = c(rnorm(25, 50, 5), rnorm(25, 55, 12), rnorm(25, 60, 7))
)
welchOmni(d, group = "group", y = "y")
#> # A tibble: 1 × 12
#>   analysis_method effect_term test_statistic test_statistic_type   df1   df2
#>   <chr>           <chr>                <dbl> <chr>               <dbl> <dbl>
#> 1 welch_anova     group                 23.6 F                       2  44.6
#> # ℹ 6 more variables: p_value <dbl>, effect_size <dbl>, effect_size_type <chr>,
#> #   effect_size_scale <chr>, effect_size_ci_low <dbl>,
#> #   effect_size_ci_high <dbl>
```
