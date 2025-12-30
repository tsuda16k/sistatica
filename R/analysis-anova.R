#' Omnibus test for repeated-measures ANOVA (rm-ANOVA)
#'
#' Runs a one-way repeated-measures ANOVA using `rstatix::anova_test()` and
#' returns a tibble in a standardized (method-agnostic) schema.
#'
#' @param df A long-format data frame.
#' @param id Name of the subject identifier column (character scalar).
#' @param group Name of the within-subject condition column (character scalar).
#' @param y Name of the dependent variable column (character scalar).
#' @param correction Sphericity correction passed to `rstatix::get_anova_table()`.
#'   One of `"auto"`, `"GG"`, `"HF"`, or `"none"`.
#'   - `"auto"`: apply GG only when sphericity is violated.
#'   - `"GG"`: report Greenhouse-Geisser corrected dfs/p-values (when applicable).
#'   - `"HF"`: report Huynh-Feldt corrected dfs/p-values (when applicable).
#'   - `"none"`: report uncorrected dfs/p-values.
#'
#' @return A tibble with one row (for the within-subject factor) and the columns below.
#'
#' @section Output columns:
#' \describe{
#'   \item{analysis_method}{Analysis identifier string. Always `"repeated_measures_anova"`.}
#'   \item{effect_term}{Effect label for schema alignment. Always `"group"`.}
#'   \item{test_statistic}{Omnibus test statistic value (F).}
#'   \item{test_statistic_type}{Omnibus test statistic type. Always `"F"`.}
#'   \item{df1}{Numerator degrees of freedom for the omnibus F-test (possibly corrected).}
#'   \item{df2}{Denominator degrees of freedom for the omnibus F-test (possibly corrected).}
#'   \item{p_value}{Omnibus p-value (possibly corrected).}
#'   \item{effect_size}{Omnibus effect size value (generalized eta-squared; GES).}
#'   \item{effect_size_type}{Effect size type label. Always `"generalized_eta_squared"`.}
#'   \item{effect_size_scale}{Effect size scale label. Always `"variance_explained"`.}
#'   \item{effect_size_ci_low}{Reserved for effect size CI lower bound (currently `NA_real_`).}
#'   \item{effect_size_ci_high}{Reserved for effect size CI upper bound (currently `NA_real_`).}
#' }
#'
#' @details
#' This function is intentionally designed for a *one-way* within-subject factor.
#' The output fixes `effect_term = "group"` to keep a consistent schema across
#' multiple omnibus/post-hoc methods.
#'
#' @examples
#' set.seed(1)
#' n = 20
#' u = rnorm(n, 0, 4)
#' w = data.frame(
#'   id = factor(seq_len(n)),
#'   u  = u,
#'   A  = 50 + u + rnorm(n, 0, 5),
#'   B  = 55 + u + rnorm(n, 0, 5),
#'   C  = 60 + u + rnorm(n, 0, 5)
#' )
#' d = tidyr::pivot_longer(w, cols = c(A, B, C), names_to = "group", values_to = "y")
#' d$group = factor(d$group)
#' rmOmni(d, id = "id", group = "group", y = "y", correction = "auto")
#' rmOmni(d, id = "id", group = "group", y = "y", correction = "GG")
#'
#' @importFrom dplyr rename mutate transmute
#' @importFrom dplyr all_of
#' @importFrom tibble as_tibble
#' @export
rmOmni = function(df, id = "id", group = "group", y = "y", correction = c("auto","GG","HF","none")){
  correction = match.arg(correction)
  # --- Standardize column names inside the function
  x = df %>% dplyr::rename(i = dplyr::all_of(id), g = dplyr::all_of(group), y = dplyr::all_of(y)) %>%
    dplyr::mutate(g = factor(g))
  # --- Run rm-ANOVA and normalize to a flat tibble (apply sphericity correction)
  x %>% rstatix::anova_test(dv = y, wid = i, within = g, effect.size = "ges") %>%
    rstatix::get_anova_table(correction = correction) %>%
    tibble::as_tibble() %>%
    dplyr::transmute(
      analysis_method = "repeated_measures_anova",
      effect_term = "group",
      test_statistic = as.numeric(F),
      test_statistic_type = "F",
      df1 = as.numeric(DFn),
      df2 = as.numeric(DFd),
      p_value = as.numeric(p),
      effect_size = as.numeric(ges),
      effect_size_type = "generalized_eta_squared",
      effect_size_scale = "variance_explained",
      effect_size_ci_low = NA_real_,
      effect_size_ci_high = NA_real_
    )
}


#' Post hoc tests for repeated-measures ANOVA (paired t-tests)
#'
#' Performs pairwise paired t-tests for a one-way repeated-measures design and
#' returns a tibble in a standardized (method-agnostic) post-hoc schema.
#'
#' - Adjusted p-values: via `stats::p.adjust()` (default: Holm).
#' - Confidence intervals: unadjusted or Bonferroni-adjusted confidence level.
#' - Effect size: Cohen's d for paired samples (dz in practice) via `rstatix::cohens_d()`.
#'
#' @param df A long-format data frame.
#' @param id Name of the subject identifier column (character scalar).
#' @param group Name of the within-subject condition column (character scalar).
#' @param y Name of the dependent variable column (character scalar).
#' @param alpha Family-wise alpha used for CI adjustment (default `0.05`).
#' @param pAdjust P-value adjustment method for multiple comparisons (default `"holm"`).
#'   Must be accepted by `stats::p.adjust()`.
#' @param ciAdjust CI adjustment method. One of `"bonferroni"` or `"none"`.
#'
#' @return A tibble with one row per pairwise contrast and the columns below.
#'
#' @section Output columns:
#' \describe{
#'   \item{analysis_method}{Analysis identifier string. Always `"repeated_measures_anova"`.}
#'   \item{effect_term}{Effect label for schema alignment. Always `"group"`.}
#'   \item{comparison_kind}{Comparison family label. Always `"pairwise"`.}
#'   \item{contrast_label}{Human-readable label, e.g., `"A vs B"`.}
#'   \item{contrast_id}{Machine-readable id, e.g., `"group:A|B"`.}
#'   \item{group1}{First condition label (maps to the first element of a pair).}
#'   \item{group2}{Second condition label (maps to the second element of a pair).}
#'   \item{center_type}{Summary location type for each condition. Always `"mean"`.}
#'   \item{group1_center}{Mean of `y` in `group1`.}
#'   \item{group2_center}{Mean of `y` in `group2`.}
#'   \item{dispersion_type}{Summary dispersion type for each condition. Always `"sd"`.}
#'   \item{group1_dispersion}{Standard deviation of `y` in `group1`.}
#'   \item{group2_dispersion}{Standard deviation of `y` in `group2`.}
#'   \item{n_group1}{Number of subjects with non-missing values in `group1`.}
#'   \item{n_group2}{Number of subjects with non-missing values in `group2`.}
#'   \item{n_total}{Number of complete paired observations used for the contrast.}
#'   \item{difference_estimate}{Estimated mean difference (group1 - group2).}
#'   \item{difference_estimate_type}{Difference estimator label. Always `"mean_difference"`.}
#'   \item{difference_direction}{Difference direction label. Always `"group1_minus_group2"`.}
#'   \item{ci_low}{Lower bound of the CI for `difference_estimate`.}
#'   \item{ci_high}{Upper bound of the CI for `difference_estimate`.}
#'   \item{ci_level}{Confidence level used for the CI (possibly Bonferroni-adjusted).}
#'   \item{ci_adjustment}{CI adjustment label: `"bonferroni_adjusted"` or `"unadjusted"`.}
#'   \item{test_statistic}{Paired t-test statistic value (t).}
#'   \item{test_statistic_type}{Test statistic type. Always `"t"`.}
#'   \item{df1}{Degrees of freedom for the t-test.}
#'   \item{df2}{Reserved (always `NA_real_` for this test family).}
#'   \item{p_value}{Raw (unadjusted) p-value for the pairwise t-test.}
#'   \item{p_value_adjusted}{Multiplicity-adjusted p-value (via `pAdjust`).}
#'   \item{p_value_adjustment}{Adjustment label, e.g., `"holm_adjusted"`.}
#'   \item{effect_size}{Effect size value (Cohen's dz).}
#'   \item{effect_size_type}{Effect size type label. Always `"cohens_dz"`.}
#'   \item{effect_size_scale}{Effect size scale label. Always `"standardized_mean_difference"`.}
#'   \item{effect_size_ci_low}{Reserved for effect size CI lower bound (currently `NA_real_`).}
#'   \item{effect_size_ci_high}{Reserved for effect size CI upper bound (currently `NA_real_`).}
#' }
#'
#' @examples
#' set.seed(2)
#' n = 18
#' u = rnorm(n, 0, 3)
#' w = data.frame(
#'   id = factor(seq_len(n)),
#'   u  = u,
#'   A  = 10 + u + rnorm(n, 0, 2),
#'   B  = 11 + u + rnorm(n, 0, 2),
#'   C  = 13 + u + rnorm(n, 0, 2)
#' )
#' d = tidyr::pivot_longer(w, cols = c(A, B, C), names_to = "group", values_to = "y")
#' d$group = factor(d$group)
#' rmPost(d, id = "id", group = "group", y = "y")
#'
#' @importFrom dplyr rename mutate group_by summarise left_join transmute select
#' @importFrom dplyr all_of any_of
#' @importFrom tidyr pivot_wider
#' @importFrom tibble tibble
#' @importFrom stats sd p.adjust
#' @importFrom utils combn
#' @export
rmPost = function(df, id = "id", group = "group", y = "y", alpha = 0.05, pAdjust = "holm", ciAdjust = c("bonferroni","none")){
  ciAdjust = match.arg(ciAdjust)
  # --- Private helpers (not exported)
  makePairs = function(lv) utils::combn(lv, 2, simplify = FALSE)
  bonfConfLevel = function(k, a = 0.05) 1 - a/(k*(k-1)/2)
  pairedN = function(x){
    lv = levels(x$g); pr = makePairs(lv)
    w = x %>% dplyr::select(i, g, y) %>% tidyr::pivot_wider(names_from = g, values_from = y)
    dplyr::bind_rows(lapply(pr, function(z){
      g1 = z[1]; g2 = z[2]
      tibble::tibble(group1 = g1, group2 = g2, n_total = sum(!is.na(w[[g1]]) & !is.na(w[[g2]])))
    }))
  }

  # --- Standardize column names and ensure factor levels
  x = df %>%
    dplyr::rename(i = dplyr::all_of(id), g = dplyr::all_of(group), y = dplyr::all_of(y)) %>%
    dplyr::mutate(g = factor(g))
  lv = levels(x$g); pr = makePairs(lv); k = length(lv)
  # --- CI level (Bonferroni level if requested)
  cl = if (ciAdjust == "bonferroni") bonfConfLevel(k, alpha) else 1 - alpha

  # --- Group summaries (mean/sd and per-condition subject counts)
  sm = x %>% dplyr::group_by(g) %>% dplyr::summarise(
    mean = mean(y, na.rm = TRUE),
    sd = stats::sd(y, na.rm = TRUE),
    n_subjects = dplyr::n_distinct(i[!is.na(y)]),
    .groups = "drop"
  )

  # --- Complete-pair counts per contrast
  nt = pairedN(x)

  # --- Paired t-tests (estimate + CI + t + df + raw p)
  tt = x %>% rstatix::t_test(y ~ g, paired = TRUE, comparisons = pr, detailed = TRUE, conf.level = cl) %>%
    dplyr::select(dplyr::any_of(c("group1","group2","estimate","conf.low","conf.high","statistic","df","p")))

  # --- Adjust p-values (stable across rstatix versions)
  tt = tt %>% dplyr::mutate(p_value_adjusted = stats::p.adjust(.data[["p"]], method = pAdjust))

  # --- Effect sizes (Cohen's dz for paired samples)
  es = x %>% rstatix::cohens_d(y ~ g, paired = TRUE, comparisons = pr) %>% dplyr::select(group1, group2, effsize)

  # --- Join components and emit standardized post-hoc schema
  tt %>%
    dplyr::left_join(sm, by = c("group1" = "g")) %>%
    dplyr::rename(group1_center = mean, group1_dispersion = sd, n_group1 = n_subjects) %>%
    dplyr::left_join(sm, by = c("group2" = "g")) %>%
    dplyr::rename(group2_center = mean, group2_dispersion = sd, n_group2 = n_subjects) %>%
    dplyr::left_join(nt, by = c("group1","group2")) %>%
    dplyr::left_join(es, by = c("group1","group2")) %>%
    dplyr::transmute(
      analysis_method = "repeated_measures_anova",
      effect_term = "group",
      comparison_kind = "pairwise",
      contrast_label = paste(group1, "vs", group2),
      contrast_id = paste0("group:", group1, "|", group2),
      group1, group2,
      center_type = "mean",
      group1_center, group2_center,
      dispersion_type = "sd",
      group1_dispersion, group2_dispersion,
      n_group1, n_group2, n_total,
      difference_estimate = as.numeric(estimate),
      difference_estimate_type = "mean_difference",
      difference_direction = "group1_minus_group2",
      ci_low = as.numeric(conf.low),
      ci_high = as.numeric(conf.high),
      ci_level = cl,
      ci_adjustment = if (ciAdjust == "bonferroni") "bonferroni_adjusted" else "unadjusted",
      test_statistic = as.numeric(statistic),
      test_statistic_type = "t",
      df1 = as.numeric(df),
      df2 = NA_real_,
      p_value = as.numeric(p),
      p_value_adjusted = as.numeric(p_value_adjusted),
      p_value_adjustment = paste0(pAdjust, "_adjusted"),
      effect_size = as.numeric(effsize),
      effect_size_type = "cohens_dz",
      effect_size_scale = "standardized_mean_difference",
      effect_size_ci_low = NA_real_,
      effect_size_ci_high = NA_real_
    )
}


#' Omnibus test for Welch one-way ANOVA
#'
#' Runs a one-way Welch ANOVA (heteroscedastic ANOVA) using
#' `rstatix::welch_anova_test()` and returns a tibble in a standardized schema.
#'
#' @param df A long-format data frame.
#' @param group Name of the between-subject group column (character scalar).
#' @param y Name of the dependent variable column (character scalar).
#'
#' @return A tibble with one row and the columns below.
#'
#' @section Output columns:
#' \describe{
#'   \item{analysis_method}{Analysis identifier string. Always `"welch_anova"`.}
#'   \item{effect_term}{Effect label for schema alignment. Always `"group"`.}
#'   \item{test_statistic}{Welch omnibus test statistic value (F).}
#'   \item{test_statistic_type}{Test statistic type. Always `"F"`.}
#'   \item{df1}{Approximate numerator degrees of freedom.}
#'   \item{df2}{Approximate denominator degrees of freedom.}
#'   \item{p_value}{Omnibus p-value.}
#'   \item{effect_size}{Approximate eta-squared computed from F and dfs.}
#'   \item{effect_size_type}{Effect size type label. Always `"eta_squared_approx"`.}
#'   \item{effect_size_scale}{Effect size scale label. Always `"variance_explained"`.}
#'   \item{effect_size_ci_low}{Reserved for effect size CI lower bound (currently `NA_real_`).}
#'   \item{effect_size_ci_high}{Reserved for effect size CI upper bound (currently `NA_real_`).}
#' }
#'
#' @details
#' The omnibus effect size is computed as an approximation:
#' \deqn{\eta^2 \approx \frac{df_1 F}{df_1 F + df_2}}
#' where \eqn{F} is the Welch statistic and \eqn{df_1, df_2} are the approximate
#' degrees of freedom returned by the Welch test.
#'
#' @examples
#' set.seed(10)
#' d = data.frame(
#'   group = factor(rep(c("G1","G2","G3"), each = 25)),
#'   y = c(rnorm(25, 50, 5), rnorm(25, 55, 12), rnorm(25, 60, 7))
#' )
#' welchOmni(d, group = "group", y = "y")
#'
#' @importFrom dplyr rename mutate transmute
#' @importFrom dplyr all_of
#' @importFrom tibble as_tibble
#' @export
welchOmni = function(df, group = "group", y = "y"){
  # --- Standardize column names inside the function
  x = df %>% dplyr::rename(g = dplyr::all_of(group), y = dplyr::all_of(y)) %>% dplyr::mutate(g = factor(g))
  # --- Welch ANOVA (handle rstatix output name differences)
  z = x %>% rstatix::welch_anova_test(y ~ g) %>% tibble::as_tibble()
  if ("DFn" %in% names(z)) z = dplyr::rename(z, df1 = DFn)
  if ("DFd" %in% names(z)) z = dplyr::rename(z, df2 = DFd)
  if (!("df1" %in% names(z)) && ("df" %in% names(z))) z = dplyr::rename(z, df1 = df)
  if (!("df2" %in% names(z))) z = dplyr::mutate(z, df2 = NA_real_)
  # --- Standardized omnibus schema (force effect_term = "group")
  z %>% dplyr::transmute(
    analysis_method = "welch_anova",
    effect_term = "group",
    test_statistic = as.numeric(statistic),
    test_statistic_type = "F",
    df1 = as.numeric(df1),
    df2 = as.numeric(df2),
    p_value = as.numeric(p),
    effect_size = (df1 * test_statistic) / (df1 * test_statistic + df2),
    effect_size_type = "eta_squared_approx",
    effect_size_scale = "variance_explained",
    effect_size_ci_low = NA_real_,
    effect_size_ci_high = NA_real_
  )
}


#' Post hoc tests for Welch one-way ANOVA (Games-Howell p-values + Welch t-test CIs)
#'
#' Performs pairwise comparisons appropriate for heteroscedastic one-way designs.
#'
#' - Adjusted p-values: `rstatix::games_howell_test()` (Games-Howell).
#' - Difference estimate and CI: `stats::t.test(var.equal = FALSE)` for each pair.
#' - Effect size: Hedges' g via `rstatix::cohens_d(hedges.correction = TRUE)`.
#'
#' @param df A long-format data frame.
#' @param group Name of the between-subject group column (character scalar).
#' @param y Name of the dependent variable column (character scalar).
#' @param alpha Family-wise alpha used for CI adjustment (default `0.05`).
#' @param ciAdjust CI adjustment method. One of `"bonferroni"` or `"none"`.
#'
#' @return A tibble with one row per pairwise contrast and the columns below.
#'
#' @section Output columns:
#' \describe{
#'   \item{analysis_method}{Analysis identifier string. Always `"welch_anova"`.}
#'   \item{effect_term}{Effect label for schema alignment. Always `"group"`.}
#'   \item{comparison_kind}{Comparison family label. Always `"pairwise"`.}
#'   \item{contrast_label}{Human-readable label, e.g., `"G1 vs G2"`.}
#'   \item{contrast_id}{Machine-readable id, e.g., `"group:G1|G2"`.}
#'   \item{group1}{First group label (maps to the first element of a pair).}
#'   \item{group2}{Second group label (maps to the second element of a pair).}
#'   \item{center_type}{Summary location type for each group. Always `"mean"`.}
#'   \item{group1_center}{Mean of `y` in `group1`.}
#'   \item{group2_center}{Mean of `y` in `group2`.}
#'   \item{dispersion_type}{Summary dispersion type for each group. Always `"sd"`.}
#'   \item{group1_dispersion}{Standard deviation of `y` in `group1`.}
#'   \item{group2_dispersion}{Standard deviation of `y` in `group2`.}
#'   \item{n_group1}{Number of non-missing observations in `group1`.}
#'   \item{n_group2}{Number of non-missing observations in `group2`.}
#'   \item{n_total}{Total observations used for the contrast (`n_group1 + n_group2`).}
#'   \item{difference_estimate}{Estimated mean difference (group1 - group2).}
#'   \item{difference_estimate_type}{Difference estimator label. Always `"mean_difference"`.}
#'   \item{difference_direction}{Difference direction label. Always `"group1_minus_group2"`.}
#'   \item{ci_low}{Lower bound of the CI for `difference_estimate` (Welch t-test).}
#'   \item{ci_high}{Upper bound of the CI for `difference_estimate` (Welch t-test).}
#'   \item{ci_level}{Confidence level used for the CI (possibly Bonferroni-adjusted).}
#'   \item{ci_adjustment}{CI adjustment label: `"bonferroni_adjusted"` or `"unadjusted"`.}
#'   \item{test_statistic}{Welch t statistic value for the pair.}
#'   \item{test_statistic_type}{Test statistic type. Always `"t"`.}
#'   \item{df1}{Approximate df from the Welch t-test for the pair.}
#'   \item{df2}{Reserved (always `NA_real_` for this test family).}
#'   \item{p_value}{Raw p-value from Welch t-test for the pair.}
#'   \item{p_value_adjusted}{Games-Howell adjusted p-value for the pair.}
#'   \item{p_value_adjustment}{Adjustment label. Always `"games_howell_adjusted"`.}
#'   \item{effect_size}{Effect size value (Hedges' g).}
#'   \item{effect_size_type}{Effect size type label. Always `"hedges_g"`.}
#'   \item{effect_size_scale}{Effect size scale label. Always `"standardized_mean_difference"`.}
#'   \item{effect_size_ci_low}{Reserved for effect size CI lower bound (currently `NA_real_`).}
#'   \item{effect_size_ci_high}{Reserved for effect size CI upper bound (currently `NA_real_`).}
#' }
#'
#' @examples
#' set.seed(11)
#' d = data.frame(
#'   group = factor(rep(c("G1","G2","G3","G4"), each = 20)),
#'   y = c(rnorm(20, 0, 1), rnorm(20, 0.3, 2), rnorm(20, 1.0, 1.2), rnorm(20, 0.6, 3))
#' )
#' welchPost(d, group = "group", y = "y")
#'
#' @importFrom dplyr rename mutate group_by summarise left_join transmute select filter
#' @importFrom dplyr all_of any_of
#' @importFrom tibble tibble
#' @importFrom stats sd
#' @importFrom utils combn
#' @export
welchPost = function(df, group = "group", y = "y", alpha = 0.05, ciAdjust = c("bonferroni","none")){
  ciAdjust = match.arg(ciAdjust)
  # --- Private helpers (not exported)
  makePairs = function(lv) utils::combn(lv, 2, simplify = FALSE)
  bonfConfLevel = function(k, a = 0.05) 1 - a/(k*(k-1)/2)

  # --- Standardize column names inside the function
  x = df %>% dplyr::rename(g = dplyr::all_of(group), y = dplyr::all_of(y)) %>% dplyr::mutate(g = factor(g))
  lv = levels(x$g); pr = makePairs(lv); k = length(lv)
  # --- CI level (Bonferroni level if requested)
  cl = if (ciAdjust == "bonferroni") bonfConfLevel(k, alpha) else 1 - alpha

  # --- Group summaries (mean/sd and n)
  sm = x %>% dplyr::group_by(g) %>% dplyr::summarise(
    mean = mean(y, na.rm = TRUE),
    sd = stats::sd(y, na.rm = TRUE),
    n = sum(!is.na(y)),
    .groups = "drop"
  )

  # --- Pairwise Welch t-tests to get estimate + CI + t + df + raw p (stable across rstatix versions)
  tt = dplyr::bind_rows(lapply(pr, function(z){
    g1 = z[1]; g2 = z[2]
    d2 = x %>% dplyr::filter(g %in% c(g1, g2)) %>% dplyr::mutate(g = droplevels(g))
    t0 = stats::t.test(y ~ g, data = d2, var.equal = FALSE, conf.level = cl)
    m1 = mean(d2$y[d2$g == g1], na.rm = TRUE); m2 = mean(d2$y[d2$g == g2], na.rm = TRUE)
    tibble::tibble(
      group1 = g1, group2 = g2,
      difference_estimate = m1 - m2,
      ci_low = unname(t0$conf.int[1]),
      ci_high = unname(t0$conf.int[2]),
      test_statistic = unname(t0$statistic),
      df1 = unname(t0$parameter),
      p_value = unname(t0$p.value)
    )
  }))

  # --- Games-Howell adjusted p-values (multiplicity-aware by design)
  gh = x %>% rstatix::games_howell_test(y ~ g) %>% tibble::as_tibble() %>%
    dplyr::select(dplyr::any_of(c("group1","group2","p.adj"))) %>%
    dplyr::rename(p_value_adjusted = p.adj)
  if (!("p_value_adjusted" %in% names(gh))) gh = gh %>% dplyr::mutate(p_value_adjusted = NA_real_)

  # --- Effect sizes (Hedges' g)
  es = x %>% rstatix::cohens_d(y ~ g, comparisons = pr, var.equal = FALSE, hedges.correction = TRUE) %>%
    dplyr::select(group1, group2, effsize) %>% dplyr::rename(effect_size = effsize)

  # --- Join summaries + GH p-values + effect sizes and emit standardized schema
  tt %>%
    dplyr::left_join(sm, by = c("group1" = "g")) %>%
    dplyr::rename(group1_center = mean, group1_dispersion = sd, n_group1 = n) %>%
    dplyr::left_join(sm, by = c("group2" = "g")) %>%
    dplyr::rename(group2_center = mean, group2_dispersion = sd, n_group2 = n) %>%
    dplyr::mutate(n_total = n_group1 + n_group2) %>%
    dplyr::left_join(gh, by = c("group1","group2")) %>%
    dplyr::left_join(es, by = c("group1","group2")) %>%
    dplyr::transmute(
      analysis_method = "welch_anova",
      effect_term = "group",
      comparison_kind = "pairwise",
      contrast_label = paste(group1, "vs", group2),
      contrast_id = paste0("group:", group1, "|", group2),
      group1, group2,
      center_type = "mean",
      group1_center, group2_center,
      dispersion_type = "sd",
      group1_dispersion, group2_dispersion,
      n_group1, n_group2, n_total,
      difference_estimate = as.numeric(difference_estimate),
      difference_estimate_type = "mean_difference",
      difference_direction = "group1_minus_group2",
      ci_low = as.numeric(ci_low),
      ci_high = as.numeric(ci_high),
      ci_level = cl,
      ci_adjustment = if (ciAdjust == "bonferroni") "bonferroni_adjusted" else "unadjusted",
      test_statistic = as.numeric(test_statistic),
      test_statistic_type = "t",
      df1 = as.numeric(df1),
      df2 = NA_real_,
      p_value = as.numeric(p_value),
      p_value_adjusted = as.numeric(p_value_adjusted),
      p_value_adjustment = "games_howell_adjusted",
      effect_size = as.numeric(effect_size),
      effect_size_type = "hedges_g",
      effect_size_scale = "standardized_mean_difference",
      effect_size_ci_low = NA_real_,
      effect_size_ci_high = NA_real_
    )
}


#' Omnibus test for the Friedman test (one-way repeated-measures, nonparametric)
#'
#' Runs a one-way Friedman test using `rstatix::friedman_test()` and returns
#' a tibble in a standardized schema.
#'
#' @param df A long-format data frame.
#' @param id Name of the subject identifier column (character scalar).
#' @param group Name of the within-subject condition column (character scalar).
#' @param y Name of the dependent variable column (character scalar).
#'
#' @return A tibble with one row and the columns below.
#'
#' @section Output columns:
#' \describe{
#'   \item{analysis_method}{Analysis identifier string. Always `"friedman_test"`.}
#'   \item{effect_term}{Effect label for schema alignment. Always `"group"`.}
#'   \item{test_statistic}{Friedman omnibus test statistic value (chi-square).}
#'   \item{test_statistic_type}{Test statistic type. Always `"chi_square"`.}
#'   \item{df1}{Degrees of freedom.}
#'   \item{df2}{Reserved (always `NA_real_` for this test family).}
#'   \item{p_value}{Omnibus p-value.}
#'   \item{effect_size}{Effect size value (Kendall's W).}
#'   \item{effect_size_type}{Effect size type label. Always `"kendalls_W"`.}
#'   \item{effect_size_scale}{Effect size scale label. Always `"agreement"`.}
#'   \item{effect_size_ci_low}{Reserved for effect size CI lower bound (currently `NA_real_`).}
#'   \item{effect_size_ci_high}{Reserved for effect size CI upper bound (currently `NA_real_`).}
#' }
#'
#' @examples
#' set.seed(21)
#' n = 24
#' u = rnorm(n, 0, 0.25)
#' w = data.frame(
#'   id = factor(seq_len(n)),
#'   A = rlnorm(n, meanlog = 3.8, sdlog = 0.25) * exp(u),
#'   B = rlnorm(n, meanlog = 3.9, sdlog = 0.25) * exp(u),
#'   C = rlnorm(n, meanlog = 4.0, sdlog = 0.25) * exp(u)
#' )
#' d = tidyr::pivot_longer(w, cols = c(A, B, C), names_to = "group", values_to = "y")
#' d$group = factor(d$group)
#' friedmanOmni(d, id = "id", group = "group", y = "y")
#'
#' @importFrom dplyr rename mutate transmute
#' @importFrom dplyr all_of
#' @importFrom tibble as_tibble
#' @export
friedmanOmni = function(df, id = "id", group = "group", y = "y"){
  # --- Standardize column names inside the function
  x = df %>% dplyr::rename(i = dplyr::all_of(id), g = dplyr::all_of(group), y = dplyr::all_of(y)) %>% dplyr::mutate(g = factor(g))
  # --- Run Friedman test and Kendall's W effect size
  ft = x %>% rstatix::friedman_test(y ~ g | i) %>% tibble::as_tibble()
  ew = x %>% rstatix::friedman_effsize(y ~ g | i) %>% tibble::as_tibble()
  # --- Standardized omnibus schema (force effect_term = "group")
  ft %>% dplyr::transmute(
    analysis_method = "friedman_test",
    effect_term = "group",
    test_statistic = as.numeric(statistic),
    test_statistic_type = "chi_square",
    df1 = as.numeric(df),
    df2 = NA_real_,
    p_value = as.numeric(p),
    effect_size = as.numeric(ew$effsize[1]),
    effect_size_type = "kendalls_W",
    effect_size_scale = "agreement",
    effect_size_ci_low = NA_real_,
    effect_size_ci_high = NA_real_
  )
}


#' Post hoc tests for the Friedman test (paired Wilcoxon + Holm p-values)
#'
#' Performs pairwise post hoc comparisons for a one-way Friedman design.
#'
#' - Difference estimate and CI: `stats::wilcox.test()` on paired differences
#'   (Hodges-Lehmann location shift + CI).
#' - Adjusted p-values: `stats::p.adjust()` (default: Holm).
#' - Effect size: rank-biserial r via `rstatix::wilcox_effsize()` (paired).
#'
#' @param df A long-format data frame.
#' @param id Name of the subject identifier column (character scalar).
#' @param group Name of the within-subject condition column (character scalar).
#' @param y Name of the dependent variable column (character scalar).
#' @param alpha Family-wise alpha used for CI adjustment (default `0.05`).
#' @param pAdjust P-value adjustment method for multiple comparisons (default `"holm"`).
#'   Must be accepted by `stats::p.adjust()`.
#' @param ciAdjust CI adjustment method. One of `"bonferroni"` or `"none"`.
#'
#' @return A tibble with one row per pairwise contrast and the columns below.
#'
#' @section Output columns:
#' \describe{
#'   \item{analysis_method}{Analysis identifier string. Always `"friedman_test"`.}
#'   \item{effect_term}{Effect label for schema alignment. Always `"group"`.}
#'   \item{comparison_kind}{Comparison family label. Always `"pairwise"`.}
#'   \item{contrast_label}{Human-readable label, e.g., `"A vs B"`.}
#'   \item{contrast_id}{Machine-readable id, e.g., `"group:A|B"`.}
#'   \item{group1}{First condition label.}
#'   \item{group2}{Second condition label.}
#'   \item{center_type}{Summary location type for each condition. Always `"median"`.}
#'   \item{group1_center}{Median of `y` in `group1`.}
#'   \item{group2_center}{Median of `y` in `group2`.}
#'   \item{dispersion_type}{Summary dispersion type for each condition. Always `"iqr"`.}
#'   \item{group1_dispersion}{IQR of `y` in `group1`.}
#'   \item{group2_dispersion}{IQR of `y` in `group2`.}
#'   \item{n_group1}{Number of subjects with non-missing values in `group1`.}
#'   \item{n_group2}{Number of subjects with non-missing values in `group2`.}
#'   \item{n_total}{Number of complete paired observations used for the contrast.}
#'   \item{difference_estimate}{Hodges-Lehmann location shift estimate (group1 - group2).}
#'   \item{difference_estimate_type}{Difference estimator label. Always `"hodges_lehmann_location_shift"`.}
#'   \item{difference_direction}{Difference direction label. Always `"group1_minus_group2"`.}
#'   \item{ci_low}{Lower bound of the CI for `difference_estimate`.}
#'   \item{ci_high}{Upper bound of the CI for `difference_estimate`.}
#'   \item{ci_level}{Confidence level used for the CI (possibly Bonferroni-adjusted).}
#'   \item{ci_adjustment}{CI adjustment label: `"bonferroni_adjusted"` or `"unadjusted"`.}
#'   \item{test_statistic}{Wilcoxon signed-rank statistic value (V).}
#'   \item{test_statistic_type}{Test statistic type. Always `"wilcoxon_V"`.}
#'   \item{df1}{Reserved (always `NA_real_` for this test family).}
#'   \item{df2}{Reserved (always `NA_real_` for this test family).}
#'   \item{p_value}{Raw p-value from the paired Wilcoxon test.}
#'   \item{p_value_adjusted}{Multiplicity-adjusted p-value (via `pAdjust`).}
#'   \item{p_value_adjustment}{Adjustment label, e.g., `"holm_adjusted"`.}
#'   \item{effect_size}{Effect size value (rank-biserial r).}
#'   \item{effect_size_type}{Effect size type label. Always `"rank_biserial_r"`.}
#'   \item{effect_size_scale}{Effect size scale label. Always `"rank_based_r"`.}
#'   \item{effect_size_ci_low}{Reserved for effect size CI lower bound (currently `NA_real_`).}
#'   \item{effect_size_ci_high}{Reserved for effect size CI upper bound (currently `NA_real_`).}
#' }
#'
#' @examples
#' set.seed(22)
#' n = 20
#' u = rnorm(n, 0, 0.25)
#' w = data.frame(
#'   id = factor(seq_len(n)),
#'   A = rlnorm(n, meanlog = 3.7, sdlog = 0.25) * exp(u),
#'   B = rlnorm(n, meanlog = 3.9, sdlog = 0.25) * exp(u),
#'   C = rlnorm(n, meanlog = 4.0, sdlog = 0.25) * exp(u)
#' )
#' d = tidyr::pivot_longer(w, cols = c(A, B, C), names_to = "group", values_to = "y")
#' d$group = factor(d$group)
#' friedmanPost(d, id = "id", group = "group", y = "y")
#'
#' @importFrom dplyr rename mutate group_by summarise left_join transmute select filter
#' @importFrom dplyr all_of
#' @importFrom tidyr pivot_wider
#' @importFrom tibble tibble
#' @importFrom stats median IQR p.adjust
#' @importFrom utils combn
#' @export
friedmanPost = function(df, id = "id", group = "group", y = "y", alpha = 0.05, pAdjust = "holm", ciAdjust = c("bonferroni","none")){
  ciAdjust = match.arg(ciAdjust)
  # --- Private helpers (not exported)
  makePairs = function(lv) utils::combn(lv, 2, simplify = FALSE)
  bonfConfLevel = function(k, a = 0.05) 1 - a/(k*(k-1)/2)
  pairedN = function(x){
    lv = levels(x$g); pr = makePairs(lv)
    w = x %>% dplyr::select(i, g, y) %>% tidyr::pivot_wider(names_from = g, values_from = y)
    dplyr::bind_rows(lapply(pr, function(z){
      g1 = z[1]; g2 = z[2]
      tibble::tibble(group1 = g1, group2 = g2, n_total = sum(!is.na(w[[g1]]) & !is.na(w[[g2]])))
    }))
  }

  # --- Standardize column names and ensure factor levels
  x = df %>%
    dplyr::rename(i = dplyr::all_of(id), g = dplyr::all_of(group), y = dplyr::all_of(y)) %>%
    dplyr::mutate(g = factor(g))
  lv = levels(x$g); pr = makePairs(lv); k = length(lv)
  # --- CI level (Bonferroni level if requested)
  cl = if (ciAdjust == "bonferroni") bonfConfLevel(k, alpha) else 1 - alpha

  # --- Group summaries (median/iqr and per-condition subject counts)
  sm = x %>% dplyr::group_by(g) %>% dplyr::summarise(
    median = stats::median(y, na.rm = TRUE),
    iqr = stats::IQR(y, na.rm = TRUE),
    n_subjects = dplyr::n_distinct(i[!is.na(y)]),
    .groups = "drop"
  )

  # --- Complete-pair counts per contrast
  nt = pairedN(x)

  # --- Pairwise Wilcoxon on paired differences (HL estimate + CI + V + raw p)
  wt = dplyr::bind_rows(lapply(pr, function(z){
    g1 = z[1]; g2 = z[2]
    w = x %>%
      dplyr::filter(g %in% c(g1, g2)) %>%
      dplyr::select(i, g, y) %>%
      tidyr::pivot_wider(names_from = g, values_from = y)
    d = w[[g1]] - w[[g2]]
    d = d[!is.na(d)]
    if (length(d) == 0) return(tibble::tibble(group1 = g1, group2 = g2, difference_estimate = NA_real_,
                                              ci_low = NA_real_, ci_high = NA_real_,
                                              test_statistic = NA_real_, p_value = NA_real_))
    t0 = suppressWarnings(stats::wilcox.test(d, mu = 0, conf.int = TRUE, conf.level = cl, exact = FALSE))
    tibble::tibble(
      group1 = g1,
      group2 = g2,
      difference_estimate = unname(t0$estimate),
      ci_low = unname(t0$conf.int[1]),
      ci_high = unname(t0$conf.int[2]),
      test_statistic = unname(t0$statistic),
      p_value = unname(t0$p.value)
    )
  })) %>% dplyr::mutate(p_value_adjusted = stats::p.adjust(.data[["p_value"]], method = pAdjust))

  # --- Effect sizes (rank-biserial r)
  es = x %>% rstatix::wilcox_effsize(y ~ g, paired = TRUE, comparisons = pr) %>% dplyr::select(group1, group2, effsize)

  # --- Join components and emit standardized post-hoc schema
  wt %>%
    dplyr::left_join(sm, by = c("group1" = "g")) %>%
    dplyr::rename(group1_center = median, group1_dispersion = iqr, n_group1 = n_subjects) %>%
    dplyr::left_join(sm, by = c("group2" = "g")) %>%
    dplyr::rename(group2_center = median, group2_dispersion = iqr, n_group2 = n_subjects) %>%
    dplyr::left_join(nt, by = c("group1","group2")) %>%
    dplyr::left_join(es, by = c("group1","group2")) %>%
    dplyr::transmute(
      analysis_method = "friedman_test",
      effect_term = "group",
      comparison_kind = "pairwise",
      contrast_label = paste(group1, "vs", group2),
      contrast_id = paste0("group:", group1, "|", group2),
      group1, group2,
      center_type = "median",
      group1_center, group2_center,
      dispersion_type = "iqr",
      group1_dispersion, group2_dispersion,
      n_group1, n_group2, n_total,
      difference_estimate = as.numeric(difference_estimate),
      difference_estimate_type = "hodges_lehmann_location_shift",
      difference_direction = "group1_minus_group2",
      ci_low = as.numeric(ci_low),
      ci_high = as.numeric(ci_high),
      ci_level = cl,
      ci_adjustment = if (ciAdjust == "bonferroni") "bonferroni_adjusted" else "unadjusted",
      test_statistic = as.numeric(test_statistic),
      test_statistic_type = "wilcoxon_V",
      df1 = NA_real_,
      df2 = NA_real_,
      p_value = as.numeric(p_value),
      p_value_adjusted = as.numeric(p_value_adjusted),
      p_value_adjustment = paste0(pAdjust, "_adjusted"),
      effect_size = as.numeric(effsize),
      effect_size_type = "rank_biserial_r",
      effect_size_scale = "rank_based_r",
      effect_size_ci_low = NA_real_,
      effect_size_ci_high = NA_real_
    )
}


#' Omnibus test for the Kruskal-Wallis test (one-way, nonparametric)
#'
#' Runs a one-way Kruskal-Wallis test using `rstatix::kruskal_test()` and returns
#' a tibble in a standardized schema.
#'
#' @param df A long-format data frame.
#' @param group Name of the between-subject group column (character scalar).
#' @param y Name of the dependent variable column (character scalar).
#'
#' @return A tibble with one row and the columns below.
#'
#' @section Output columns:
#' \describe{
#'   \item{analysis_method}{Analysis identifier string. Always `"kruskal_wallis_test"`.}
#'   \item{effect_term}{Effect label for schema alignment. Always `"group"`.}
#'   \item{test_statistic}{Kruskal-Wallis omnibus test statistic value (chi-square).}
#'   \item{test_statistic_type}{Test statistic type. Always `"chi_square"`.}
#'   \item{df1}{Degrees of freedom.}
#'   \item{df2}{Reserved (always `NA_real_` for this test family).}
#'   \item{p_value}{Omnibus p-value.}
#'   \item{effect_size}{Effect size value (epsilon-squared).}
#'   \item{effect_size_type}{Effect size type label. Always `"epsilon_squared"`.}
#'   \item{effect_size_scale}{Effect size scale label. Always `"rank_based_variance_explained"`.}
#'   \item{effect_size_ci_low}{Reserved for effect size CI lower bound (currently `NA_real_`).}
#'   \item{effect_size_ci_high}{Reserved for effect size CI upper bound (currently `NA_real_`).}
#' }
#'
#' @examples
#' data("PlantGrowth")
#' d = PlantGrowth
#' d$group = factor(d$group)
#' names(d) = c("y","group")
#' kruskalOmni(d, group = "group", y = "y")
#'
#' @importFrom dplyr rename mutate transmute
#' @importFrom dplyr all_of
#' @importFrom tibble as_tibble
#' @export
kruskalOmni = function(df, group = "group", y = "y"){
  # --- Standardize column names inside the function
  x = df %>% dplyr::rename(g = dplyr::all_of(group), y = dplyr::all_of(y)) %>% dplyr::mutate(g = factor(g))
  # --- Run Kruskal-Wallis test and effect size
  kt = x %>% rstatix::kruskal_test(y ~ g) %>% tibble::as_tibble()
  es = x %>% rstatix::kruskal_effsize(y ~ g) %>% tibble::as_tibble()
  # --- Standardized omnibus schema (force effect_term = "group")
  kt %>% dplyr::transmute(
    analysis_method = "kruskal_wallis_test",
    effect_term = "group",
    test_statistic = as.numeric(statistic),
    test_statistic_type = "chi_square",
    df1 = as.numeric(df),
    df2 = NA_real_,
    p_value = as.numeric(p),
    effect_size = as.numeric(es$effsize[1]),
    effect_size_type = "epsilon_squared",
    effect_size_scale = "rank_based_variance_explained",
    effect_size_ci_low = NA_real_,
    effect_size_ci_high = NA_real_
  )
}


#' Post hoc tests for the Kruskal-Wallis test (Dunn p-values + Wilcoxon HL CIs)
#'
#' Performs pairwise post hoc comparisons after a Kruskal-Wallis test.
#'
#' - Adjusted p-values: `rstatix::dunn_test()` (default: Holm).
#' - Difference estimate and CI: `stats::wilcox.test()` between groups
#'   (Hodges-Lehmann location shift + CI).
#' - Effect size: rank-biserial r via `rstatix::wilcox_effsize()` (unpaired).
#'
#' @param df A long-format data frame.
#' @param group Name of the between-subject group column (character scalar).
#' @param y Name of the dependent variable column (character scalar).
#' @param alpha Family-wise alpha used for CI adjustment (default `0.05`).
#' @param pAdjust P-value adjustment method for multiple comparisons (default `"holm"`).
#'   Passed to `rstatix::dunn_test()` and recorded in output.
#' @param ciAdjust CI adjustment method. One of `"bonferroni"` or `"none"`.
#'
#' @return A tibble with one row per pairwise contrast and the columns below.
#'
#' @section Output columns:
#' \describe{
#'   \item{analysis_method}{Analysis identifier string. Always `"kruskal_wallis_test"`.}
#'   \item{effect_term}{Effect label for schema alignment. Always `"group"`.}
#'   \item{comparison_kind}{Comparison family label. Always `"pairwise"`.}
#'   \item{contrast_label}{Human-readable label, e.g., `"A vs B"`.}
#'   \item{contrast_id}{Machine-readable id, e.g., `"group:A|B"`.}
#'   \item{group1}{First group label.}
#'   \item{group2}{Second group label.}
#'   \item{center_type}{Summary location type for each group. Always `"median"`.}
#'   \item{group1_center}{Median of `y` in `group1`.}
#'   \item{group2_center}{Median of `y` in `group2`.}
#'   \item{dispersion_type}{Summary dispersion type for each group. Always `"iqr"`.}
#'   \item{group1_dispersion}{IQR of `y` in `group1`.}
#'   \item{group2_dispersion}{IQR of `y` in `group2`.}
#'   \item{n_group1}{Number of non-missing observations in `group1`.}
#'   \item{n_group2}{Number of non-missing observations in `group2`.}
#'   \item{n_total}{Total observations used for the contrast (`n_group1 + n_group2`).}
#'   \item{difference_estimate}{Hodges-Lehmann location shift estimate (group1 - group2).}
#'   \item{difference_estimate_type}{Difference estimator label. Always `"hodges_lehmann_location_shift"`.}
#'   \item{difference_direction}{Difference direction label. Always `"group1_minus_group2"`.}
#'   \item{ci_low}{Lower bound of the CI for `difference_estimate`.}
#'   \item{ci_high}{Upper bound of the CI for `difference_estimate`.}
#'   \item{ci_level}{Confidence level used for the CI (possibly Bonferroni-adjusted).}
#'   \item{ci_adjustment}{CI adjustment label: `"bonferroni_adjusted"` or `"unadjusted"`.}
#'   \item{test_statistic}{Dunn test statistic value (Z).}
#'   \item{test_statistic_type}{Test statistic type. Always `"dunn_Z"`.}
#'   \item{df1}{Reserved (always `NA_real_` for this test family).}
#'   \item{df2}{Reserved (always `NA_real_` for this test family).}
#'   \item{p_value}{Raw p-value from Dunn test.}
#'   \item{p_value_adjusted}{Multiplicity-adjusted p-value from Dunn test.}
#'   \item{p_value_adjustment}{Adjustment label, e.g., `"holm_adjusted"`.}
#'   \item{effect_size}{Effect size value (rank-biserial r).}
#'   \item{effect_size_type}{Effect size type label. Always `"rank_biserial_r"`.}
#'   \item{effect_size_scale}{Effect size scale label. Always `"rank_based_r"`.}
#'   \item{effect_size_ci_low}{Reserved for effect size CI lower bound (currently `NA_real_`).}
#'   \item{effect_size_ci_high}{Reserved for effect size CI upper bound (currently `NA_real_`).}
#' }
#'
#' @examples
#' set.seed(30)
#' d = data.frame(
#'   group = factor(rep(c("A","B","C"), each = 20)),
#'   y = c(rlnorm(20, 3.8, 0.3), rlnorm(20, 3.9, 0.35), rlnorm(20, 4.1, 0.4))
#' )
#' kruskalPost(d, group = "group", y = "y")
#'
#' @importFrom dplyr rename mutate group_by summarise left_join transmute select filter
#' @importFrom dplyr all_of any_of
#' @importFrom tibble tibble
#' @importFrom stats median IQR
#' @importFrom utils combn
#' @export
kruskalPost = function(df, group = "group", y = "y", alpha = 0.05, pAdjust = "holm", ciAdjust = c("bonferroni","none")){
  ciAdjust = match.arg(ciAdjust)
  # --- Private helpers (not exported)
  makePairs = function(lv) utils::combn(lv, 2, simplify = FALSE)
  bonfConfLevel = function(k, a = 0.05) 1 - a/(k*(k-1)/2)

  # --- Standardize column names and ensure factor levels
  x = df %>% dplyr::rename(g = dplyr::all_of(group), y = dplyr::all_of(y)) %>% dplyr::mutate(g = factor(g))
  lv = levels(x$g); pr = makePairs(lv); k = length(lv)
  # --- CI level (Bonferroni level if requested)
  cl = if (ciAdjust == "bonferroni") bonfConfLevel(k, alpha) else 1 - alpha

  # --- Group summaries (median/iqr and n per group)
  sm = x %>% dplyr::group_by(g) %>% dplyr::summarise(
    median = stats::median(y, na.rm = TRUE),
    iqr = stats::IQR(y, na.rm = TRUE),
    n_observations = sum(!is.na(y)),
    .groups = "drop"
  )

  # --- Dunn test for p-values and Z statistic
  dn = x %>% rstatix::dunn_test(y ~ g, p.adjust.method = pAdjust) %>% tibble::as_tibble() %>%
    dplyr::select(dplyr::any_of(c("group1","group2","statistic","p","p.adj")))
  if ("p.adj" %in% names(dn)) {
    dn = dn %>% dplyr::mutate(p_value_adjusted = .data[["p.adj"]])
  } else {
    dn = dn %>% dplyr::mutate(p_value_adjusted = stats::p.adjust(.data[["p"]], method = pAdjust))
  }

  # --- Wilcoxon (unpaired) to obtain HL estimate + CI per contrast (robust to rstatix version differences)
  wl = dplyr::bind_rows(lapply(pr, function(z){
    g1 = z[1]; g2 = z[2]
    d2 = x %>% dplyr::filter(g %in% c(g1, g2)) %>% dplyr::mutate(g = droplevels(g))
    y1 = d2$y[d2$g == g1]; y2 = d2$y[d2$g == g2]
    if (sum(!is.na(y1)) == 0 || sum(!is.na(y2)) == 0) {
      return(tibble::tibble(group1 = g1, group2 = g2, difference_estimate = NA_real_, ci_low = NA_real_,
                            ci_high = NA_real_, p_value_wilcox = NA_real_))
    }
    t0 = suppressWarnings(stats::wilcox.test(y1, y2, conf.int = TRUE, conf.level = cl, exact = FALSE))
    tibble::tibble(
      group1 = g1,
      group2 = g2,
      difference_estimate = unname(t0$estimate),
      ci_low = unname(t0$conf.int[1]),
      ci_high = unname(t0$conf.int[2]),
      p_value_wilcox = unname(t0$p.value)
    )
  }))

  # --- Effect sizes (rank-biserial r; unpaired)
  es = x %>% rstatix::wilcox_effsize(y ~ g, paired = FALSE, comparisons = pr) %>%
    dplyr::select(group1, group2, effsize)

  # --- Join summaries + HL estimate/CI + Dunn results + effect sizes; emit standardized schema
  wl %>%
    dplyr::left_join(sm, by = c("group1" = "g")) %>%
    dplyr::rename(group1_center = median, group1_dispersion = iqr, n_group1 = n_observations) %>%
    dplyr::left_join(sm, by = c("group2" = "g")) %>%
    dplyr::rename(group2_center = median, group2_dispersion = iqr, n_group2 = n_observations) %>%
    dplyr::mutate(n_total = n_group1 + n_group2) %>%
    dplyr::left_join(dn, by = c("group1","group2")) %>%
    dplyr::left_join(es, by = c("group1","group2")) %>%
    dplyr::transmute(
      analysis_method = "kruskal_wallis_test",
      effect_term = "group",
      comparison_kind = "pairwise",
      contrast_label = paste(group1, "vs", group2),
      contrast_id = paste0("group:", group1, "|", group2),
      group1, group2,
      center_type = "median",
      group1_center, group2_center,
      dispersion_type = "iqr",
      group1_dispersion, group2_dispersion,
      n_group1, n_group2, n_total,
      difference_estimate = as.numeric(difference_estimate),
      difference_estimate_type = "hodges_lehmann_location_shift",
      difference_direction = "group1_minus_group2",
      ci_low = as.numeric(ci_low),
      ci_high = as.numeric(ci_high),
      ci_level = cl,
      ci_adjustment = if (ciAdjust == "bonferroni") "bonferroni_adjusted" else "unadjusted",
      test_statistic = as.numeric(statistic),
      test_statistic_type = "dunn_Z",
      df1 = NA_real_,
      df2 = NA_real_,
      p_value = as.numeric(p),
      p_value_adjusted = as.numeric(p_value_adjusted),
      p_value_adjustment = paste0(pAdjust, "_adjusted"),
      effect_size = as.numeric(effsize),
      effect_size_type = "rank_biserial_r",
      effect_size_scale = "rank_based_r",
      effect_size_ci_low = NA_real_,
      effect_size_ci_high = NA_real_
    )
}

