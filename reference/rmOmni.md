# Omnibus test for repeated-measures ANOVA (rm-ANOVA)

Runs a one-way repeated-measures ANOVA using \`rstatix::anova_test()\`
and returns a tibble in a standardized (method-agnostic) schema.

## Usage

``` r
rmOmni(
  df,
  id = "id",
  group = "group",
  y = "y",
  correction = c("auto", "GG", "HF", "none")
)
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

- correction:

  Sphericity correction passed to \`rstatix::get_anova_table()\`. One of
  \`"auto"\`, \`"GG"\`, \`"HF"\`, or \`"none"\`. - \`"auto"\`: apply GG
  only when sphericity is violated. - \`"GG"\`: report
  Greenhouse-Geisser corrected dfs/p-values (when applicable). -
  \`"HF"\`: report Huynh-Feldt corrected dfs/p-values (when
  applicable). - \`"none"\`: report uncorrected dfs/p-values.

## Value

A tibble with one row (for the within-subject factor) and the columns
below.

## Details

This function is intentionally designed for a \*one-way\* within-subject
factor. The output fixes \`effect_term = "group"\` to keep a consistent
schema across multiple omnibus/post-hoc methods.

## Output columns

- analysis_method:

  Analysis identifier string. Always \`"repeated_measures_anova"\`.

- effect_term:

  Effect label for schema alignment. Always \`"group"\`.

- test_statistic:

  Omnibus test statistic value (F).

- test_statistic_type:

  Omnibus test statistic type. Always \`"F"\`.

- df1:

  Numerator degrees of freedom for the omnibus F-test (possibly
  corrected).

- df2:

  Denominator degrees of freedom for the omnibus F-test (possibly
  corrected).

- p_value:

  Omnibus p-value (possibly corrected).

- effect_size:

  Omnibus effect size value (generalized eta-squared; GES).

- effect_size_type:

  Effect size type label. Always \`"generalized_eta_squared"\`.

- effect_size_scale:

  Effect size scale label. Always \`"variance_explained"\`.

- effect_size_ci_low:

  Reserved for effect size CI lower bound (currently \`NA_real\_\`).

- effect_size_ci_high:

  Reserved for effect size CI upper bound (currently \`NA_real\_\`).

## Examples

``` r
set.seed(1)
n = 20
u = rnorm(n, 0, 4)
w = data.frame(
  id = factor(seq_len(n)),
  u  = u,
  A  = 50 + u + rnorm(n, 0, 5),
  B  = 55 + u + rnorm(n, 0, 5),
  C  = 60 + u + rnorm(n, 0, 5)
)
d = tidyr::pivot_longer(w, cols = c(A, B, C), names_to = "group", values_to = "y")
d$group = factor(d$group)
rmOmni(d, id = "id", group = "group", y = "y", correction = "auto")
#> # A tibble: 1 × 12
#>   analysis_method     effect_term test_statistic test_statistic_type   df1   df2
#>   <chr>               <chr>                <dbl> <chr>               <dbl> <dbl>
#> 1 repeated_measures_… group                 28.2 F                       2    38
#> # ℹ 6 more variables: p_value <dbl>, effect_size <dbl>, effect_size_type <chr>,
#> #   effect_size_scale <chr>, effect_size_ci_low <dbl>,
#> #   effect_size_ci_high <dbl>
rmOmni(d, id = "id", group = "group", y = "y", correction = "GG")
#> # A tibble: 1 × 12
#>   analysis_method     effect_term test_statistic test_statistic_type   df1   df2
#>   <chr>               <chr>                <dbl> <chr>               <dbl> <dbl>
#> 1 repeated_measures_… group                 28.2 F                    1.95  37.1
#> # ℹ 6 more variables: p_value <dbl>, effect_size <dbl>, effect_size_type <chr>,
#> #   effect_size_scale <chr>, effect_size_ci_low <dbl>,
#> #   effect_size_ci_high <dbl>
```
