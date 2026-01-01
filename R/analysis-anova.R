###############################################################################
# Internal helpers (unexported)
###############################################################################

.first_present = function(candidates, available) {
  out = candidates[candidates %in% available][1]
  if (length(out) == 0 || is.na(out)) return(NA_character_)
  out
}

.get_num1 = function(df, candidates) {
  nm = .first_present(candidates, names(df))
  if (is.na(nm)) return(NA_real_)
  as.numeric(df[[nm]][1])
}

.get_num_vec = function(df, candidates, n) {
  nm = .first_present(candidates, names(df))
  if (is.na(nm)) return(rep(NA_real_, n))
  as.numeric(df[[nm]])
}

.coerce_factor_cols = function(data, cols) {
  for (nm in cols) {
    if (!is.null(nm) && !is.na(nm) && nm %in% names(data)) {
      data[[nm]] = base::droplevels(as.factor(data[[nm]]))
    }
  }
  data
}

.complete_ids_within = function(data, id_nm, within_nm, dv_nm) {
  dat_cc =
    data %>%
    dplyr::filter(!is.na(.data[[id_nm]]), !is.na(.data[[within_nm]]), !is.na(.data[[dv_nm]]))

  lv = sort(unique(dat_cc[[within_nm]]))

  complete_ids =
    dat_cc %>%
    dplyr::distinct(.data[[id_nm]], .data[[within_nm]]) %>%
    dplyr::count(.data[[id_nm]], name = "k") %>%
    dplyr::filter(.data$k == length(lv)) %>%
    dplyr::select(.data[[id_nm]])

  complete_ids
}

.omega2_from_F = function(F_val, df1_val, df2_val) {
  es_df = base::try(as.data.frame(effectsize::F_to_omega2(F_val, df1_val, df2_val)), silent = TRUE)
  if (inherits(es_df, "try-error")) return(NA_real_)

  if ("Omega2" %in% names(es_df)) return(as.numeric(es_df$Omega2[1]))

  cand = grep("omega|Omega", names(es_df), ignore.case = TRUE, value = TRUE)
  if (length(cand) >= 1) return(as.numeric(es_df[[cand[1]]][1]))

  NA_real_
}

# CI adjustment is independent from p-value adjustment.
# Only methods that map cleanly to a single adjusted CI level are supported.
.ci_adjust_info = function(conf_level, ci_adjust_method, m) {
  # Validate conf_level
  if (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level)) {
    stop("conf_level must be a single finite numeric value.", call. = FALSE)
  }
  if (!(conf_level > 0 && conf_level < 1)) {
    stop("conf_level must be between 0 and 1 (exclusive).", call. = FALSE)
  }

  ci_adjust_method = tolower(ci_adjust_method)
  if (!ci_adjust_method %in% c("none", "bonferroni", "sidak")) {
    stop("ci_adjust_method must be one of: 'none', 'bonferroni', 'sidak'.", call. = FALSE)
  }

  alpha = 1 - conf_level
  m = max(as.integer(m), 1L)

  if (ci_adjust_method == "none") {
    list(ci_level_adj = conf_level, ci_adjust_method = "none")
  } else if (ci_adjust_method == "bonferroni") {
    list(ci_level_adj = 1 - (alpha / m), ci_adjust_method = "bonferroni")
  } else {
    # sidak
    alpha_adj = 1 - (1 - alpha)^(1 / m)
    list(ci_level_adj = 1 - alpha_adj, ci_adjust_method = "sidak")
  }
}

# Extract effect size + CI robustly:
# - Caller provides ordered regex patterns (preference order).
# - We select the first pattern that matches at least one numeric column.
# - If multiple columns match that same pattern, we take the first by column-name order and warn.
# - If nothing matches, warn and return NA.
.extract_es_ci = function(es_df, es_regex, warn_ctx = "") {
  num_cols = names(es_df)[vapply(es_df, is.numeric, logical(1))]

  low_nm = grep("CI_low|CI\\.low|Lower", names(es_df), ignore.case = TRUE, value = TRUE)[1]
  high_nm = grep("CI_high|CI\\.high|Upper", names(es_df), ignore.case = TRUE, value = TRUE)[1]

  es_nm = NA_character_
  ambiguous = FALSE

  for (rx in es_regex) {
    hits = grep(rx, num_cols, perl = TRUE, value = TRUE)
    if (length(hits) >= 1) {
      es_nm = hits[1]
      ambiguous = length(hits) > 1
      break
    }
  }

  if (is.na(es_nm)) {
    msg = "Could not identify effect-size column from effectsize output; returning NA."
    if (nzchar(warn_ctx)) msg = paste0(msg, " Context: ", warn_ctx)
    warning(msg, call. = FALSE)
  } else if (ambiguous) {
    msg = paste0(
      "Multiple effect-size columns matched; using the first match '", es_nm, "'."
    )
    if (nzchar(warn_ctx)) msg = paste0(msg, " Context: ", warn_ctx)
    warning(msg, call. = FALSE)
  }

  list(
    es = if (!is.na(es_nm)) as.numeric(es_df[[es_nm]][1]) else NA_real_,
    low = if (!is.na(low_nm)) as.numeric(es_df[[low_nm]][1]) else NA_real_,
    high = if (!is.na(high_nm)) as.numeric(es_df[[high_nm]][1]) else NA_real_
  )
}

.warn_aggregate_once = function(messages, header) {
  messages = unique(messages[!is.na(messages) & nzchar(messages)])
  if (length(messages) == 0) return(invisible(NULL))
  warning(paste0(header, " ", paste(messages, collapse = " | ")), call. = FALSE)
  invisible(NULL)
}

# For within-subject posthoc: build wide matrix for a subset (two levels or all levels),
# optionally aggregating duplicated ID-within cells.
.build_wide_mat = function(dat, id_nm, within_nm, dv_nm,
                           duplicate = c("mean", "first", "median", "error"),
                           duplicate_na_rm = TRUE) {
  duplicate = match.arg(duplicate)

  tab = base::table(dat[[id_nm]], dat[[within_nm]])
  has_dup = any(tab > 1L)

  agg_warning = character(0)

  if (has_dup && duplicate == "error") {
    stop("Duplicated ID-within cells detected. Set duplicate = 'mean'/'median'/'first' to aggregate.", call. = FALSE)
  }

  if (duplicate == "first") {
    FUN = function(z) {
      if (duplicate_na_rm) z = z[!is.na(z)]
      if (length(z) == 0) return(NA_real_)
      z[1]
    }
  } else if (duplicate == "mean") {
    FUN = function(z) base::mean(z, na.rm = duplicate_na_rm)
  } else if (duplicate == "median") {
    FUN = function(z) stats::median(z, na.rm = duplicate_na_rm)
  } else {
    # error already handled
    FUN = function(z) base::mean(z, na.rm = duplicate_na_rm)
  }

  if (has_dup && duplicate != "error") {
    agg_warning = paste0("Duplicated ID x within cells were aggregated using duplicate='", duplicate, "'.")
  }

  mat = base::tapply(dat[[dv_nm]], list(dat[[id_nm]], dat[[within_nm]]), FUN)
  list(mat = mat, agg_warning = agg_warning)
}

.empty_posthoc_tibble = function() {
  tibble::tibble(
    design = character(),
    method = character(),
    term = character(),
    group1 = character(),
    group2 = character(),
    statistic_type = character(),
    statistic = numeric(),
    df = numeric(),
    p = numeric(),
    p_adj = numeric(),
    p_adj_mark = character(),
    p_adjust_method = character(),
    effect_size = numeric(),
    effect_size_type = character(),
    ci_level = numeric(),
    ci_low = numeric(),
    ci_high = numeric(),
    ci_adjust_method = character(),
    ci_level_adj = numeric(),
    ci_low_adj = numeric(),
    ci_high_adj = numeric(),
    n_pair = numeric(),
    missing = character()
  )
}


###############################################################################
# Omnibus tests (4 functions)
###############################################################################

#' Repeated-measures ANOVA (omnibus) with generalized eta-squared
#'
#' Performs the omnibus test for a one-factor repeated-measures ANOVA (within-subject design)
#' using \code{rstatix::anova_test()} and returns a tidy one-row tibble with a standardized
#' omnibus schema and generalized eta-squared (GES).
#'
#' @details
#' **Test**: One-way repeated-measures ANOVA (within-subject). \cr
#' **Statistic**: F with numerator/denominator df. \cr
#' **Effect size**: Generalized eta-squared (GES) computed by \code{rstatix::anova_test(effect.size = "ges")}.
#'
#' \strong{Missing-data handling (\code{missing})}
#' \itemize{
#'   \item \code{"pairwise"} (default): Removes rows with missing DV/within/id (row-wise filtering).
#'   \item \code{"complete_case"}: After row-wise filtering, keeps only IDs that have non-missing DV values
#'   across all within-levels (complete blocks).
#' }
#'
#' @param data A data.frame in long format.
#' @param dv Dependent variable (numeric). Unquoted name or string.
#' @param within Within-subject factor. Unquoted name or string.
#' @param id Subject identifier. Unquoted name or string.
#' @param missing Missing-data strategy: \code{"pairwise"} (default) or \code{"complete_case"}.
#'
#' @return
#' A one-row tibble (omnibus schema).
#' Columns:
#' \itemize{
#'   \item \code{design}: Design label (\code{"parametric_within"}).
#'   \item \code{method}: Omnibus method label (\code{"rm_anova"}).
#'   \item \code{term}: Name of the within-subject factor (the tested term).
#'   \item \code{statistic_type}: Type of omnibus test statistic (\code{"F"}).
#'   \item \code{statistic}: Omnibus test statistic value (F).
#'   \item \code{df1}: Numerator degrees of freedom for F.
#'   \item \code{df2}: Denominator degrees of freedom for F.
#'   \item \code{df}: Single df column (always \code{NA} for this test; included for schema consistency).
#'   \item \code{p}: Omnibus p-value.
#'   \item \code{effect_size}: Effect size value (generalized eta-squared).
#'   \item \code{effect_size_type}: Effect size label (\code{"ges"}).
#'   \item \code{missing}: Missing-data strategy used (\code{"pairwise"} or \code{"complete_case"}).
#'   \item \code{n_obs}: Number of rows (observations) used in the analysis after filtering for missingness.
#'   \item \code{n_id}: Number of unique IDs used in the analysis after filtering for missingness.
#' }
#'
#' @examples
#' dat = datasets::CO2 %>%
#'   dplyr::filter(Type == "Quebec", Treatment == "chilled") %>%
#'   dplyr::mutate(conc = factor(conc))
#'
#' dat$uptake[dat$Plant == dat$Plant[1] & dat$conc == dat$conc[1]] = NA
#'
#' rm_anova_omnibus(dat, dv = uptake, within = conc, id = Plant, missing = "pairwise")
#' rm_anova_omnibus(dat, dv = uptake, within = conc, id = Plant, missing = "complete_case")
#'
#' @export
#' @importFrom magrittr %>%
rm_anova_omnibus = function(data, dv, within, id, missing = c("pairwise", "complete_case")) {
  missing = match.arg(missing)

  dv_sym = rlang::ensym(dv)
  within_sym = rlang::ensym(within)
  id_sym = rlang::ensym(id)

  dv_nm = rlang::as_name(dv_sym)
  within_nm = rlang::as_name(within_sym)
  id_nm = rlang::as_name(id_sym)

  dat0 =
    data %>%
    dplyr::filter(!is.na(.data[[dv_nm]]), !is.na(.data[[within_nm]]), !is.na(.data[[id_nm]]))

  dat0 = .coerce_factor_cols(dat0, c(within_nm, id_nm))

  if (missing == "complete_case") {
    complete_ids = .complete_ids_within(dat0, id_nm = id_nm, within_nm = within_nm, dv_nm = dv_nm)
    dat0 = dat0 %>% dplyr::semi_join(complete_ids, by = stats::setNames(id_nm, id_nm))
  }

  aov_res =
    rstatix::anova_test(
      data = dat0,
      dv = !!dv_sym,
      wid = !!id_sym,
      within = !!within_sym,
      effect.size = "ges"
    )

  aov_tbl =
    rstatix::get_anova_table(aov_res) %>%
    tibble::as_tibble()

  effect_col = .first_present(c("Effect", "effect", "term"), names(aov_tbl))
  if (is.na(effect_col)) stop("Unexpected anova table format: could not find term/effect column.")

  aov_tbl[[effect_col]] = as.character(aov_tbl[[effect_col]])
  row = aov_tbl[aov_tbl[[effect_col]] == within_nm, , drop = FALSE]
  if (nrow(row) != 1) stop("Could not uniquely identify the within-term row in the ANOVA table.")

  tibble::tibble(
    design = "parametric_within",
    method = "rm_anova",
    term = within_nm,
    statistic_type = "F",
    statistic = .get_num1(row, c("F", "statistic")),
    df1 = .get_num1(row, c("DFn", "df1")),
    df2 = .get_num1(row, c("DFd", "df2")),
    df = NA_real_,
    p = .get_num1(row, c("p", "p.value")),
    effect_size = .get_num1(row, c("ges", "GES")),
    effect_size_type = "ges",
    missing = missing,
    n_obs = as.numeric(nrow(dat0)),
    n_id = as.numeric(dplyr::n_distinct(dat0[[id_nm]]))
  )
}


#' Welch one-way ANOVA (omnibus) with omega-squared
#'
#' Performs the omnibus Welch one-way ANOVA for a between-subject factor using
#' \code{rstatix::welch_anova_test()}, and computes omega-squared from the resulting
#' F statistic and degrees of freedom using \code{effectsize::F_to_omega2()}.
#'
#' @details
#' **Test**: Welch one-way ANOVA. \cr
#' **Statistic**: F with Welch-approximated df1/df2. \cr
#' **Effect size**: Omega-squared computed from F and df via \code{effectsize::F_to_omega2()}.
#'
#' \strong{Note}:
#' The omega-squared value is computed from the Welch F statistic and associated degrees of freedom.
#' This provides a convenient effect-size summary, but should be interpreted as an approximation.
#'
#' \strong{Missing-data handling (\code{missing})}
#' \itemize{
#'   \item \code{"pairwise"} and \code{"complete_case"} behave identically here:
#'   rows with missing DV/group are removed.
#' }
#'
#' @param data A data.frame.
#' @param dv Dependent variable (numeric). Unquoted name or string.
#' @param group Grouping factor. Unquoted name or string.
#' @param missing Missing-data strategy: \code{"pairwise"} (default) or \code{"complete_case"}.
#'
#' @return
#' A one-row tibble (omnibus schema).
#' Columns:
#' \itemize{
#'   \item \code{design}: Design label (\code{"parametric_between_unequal_var"}).
#'   \item \code{method}: Omnibus method label (\code{"welch_anova"}).
#'   \item \code{term}: Name of the grouping factor (the tested term).
#'   \item \code{statistic_type}: Type of omnibus test statistic (\code{"F"}).
#'   \item \code{statistic}: Omnibus test statistic value (F; Welch test).
#'   \item \code{df1}: Numerator degrees of freedom for F (Welch).
#'   \item \code{df2}: Denominator degrees of freedom for F (Welch).
#'   \item \code{df}: Single df column (always \code{NA} for this test; included for schema consistency).
#'   \item \code{p}: Omnibus p-value.
#'   \item \code{effect_size}: Effect size value (omega-squared; approximation from Welch F and df).
#'   \item \code{effect_size_type}: Effect size label (\code{"omega2"}).
#'   \item \code{missing}: Missing-data strategy used (\code{"pairwise"} or \code{"complete_case"}).
#'   \item \code{n_obs}: Number of rows (observations) used in the analysis after filtering for missingness.
#'   \item \code{n_id}: Always \code{NA} for between-subject designs in this package (reserved column).
#' }
#'
#' @examples
#' welch_anova_omnibus(datasets::PlantGrowth, dv = weight, group = group)
#'
#' @export
#' @importFrom magrittr %>%
welch_anova_omnibus = function(data, dv, group, missing = c("pairwise", "complete_case")) {
  missing = match.arg(missing)

  dv_sym = rlang::ensym(dv)
  group_sym = rlang::ensym(group)

  dv_nm = rlang::as_name(dv_sym)
  group_nm = rlang::as_name(group_sym)

  dat0 =
    data %>%
    dplyr::filter(!is.na(.data[[dv_nm]]), !is.na(.data[[group_nm]]))

  dat0 = .coerce_factor_cols(dat0, c(group_nm))

  w =
    rstatix::welch_anova_test(
      data = dat0,
      formula = stats::as.formula(paste(dv_nm, "~", group_nm))
    ) %>%
    tibble::as_tibble()

  F_val = .get_num1(w, c("statistic", "F"))
  df1_val = .get_num1(w, c("DFn", "df1"))
  df2_val = .get_num1(w, c("DFd", "df2"))
  p_val = .get_num1(w, c("p", "p.value"))

  omega2 = .omega2_from_F(F_val, df1_val, df2_val)

  tibble::tibble(
    design = "parametric_between_unequal_var",
    method = "welch_anova",
    term = group_nm,
    statistic_type = "F",
    statistic = F_val,
    df1 = df1_val,
    df2 = df2_val,
    df = NA_real_,
    p = p_val,
    effect_size = omega2,
    effect_size_type = "omega2",
    missing = missing,
    n_obs = as.numeric(nrow(dat0)),
    n_id = NA_real_
  )
}


#' Friedman test (omnibus) with Kendall's W
#'
#' Performs the omnibus Friedman test for a one-factor within-subject design and returns
#' Kendall's W.
#'
#' @details
#' **Test**: Friedman test. \cr
#' **Requirement**: unreplicated complete block design. \cr
#' **Statistic**: Chi-squared with df = k - 1. \cr
#' **Effect size**: Kendall's W computed as \eqn{W = \chi^2 / (n (k - 1))}.
#'
#' \strong{Implementation note}:
#' The function stops if duplicated ID x within cells are detected. After that check, the wide
#' ID x within matrix is built using \code{tapply(...)} with \code{mean} as the aggregation function.
#' Because duplicates are ruled out, each cell is effectively a single value; using \code{mean}
#' here is a simple and stable way to map long data to a wide matrix without introducing extra
#' dependencies or complicated reshaping logic.
#'
#' \strong{Missing-data handling (\code{missing})}
#' \itemize{
#'   \item \code{"pairwise"} and \code{"complete_case"} both enforce complete blocks for the omnibus
#'   test by keeping only complete rows in the wide matrix.
#' }
#'
#' @param data A data.frame in long format.
#' @param dv Dependent variable (numeric). Unquoted name or string.
#' @param within Within-subject factor. Unquoted name or string.
#' @param id Subject identifier. Unquoted name or string.
#' @param missing Missing-data strategy: \code{"pairwise"} (default) or \code{"complete_case"}.
#'
#' @return
#' A one-row tibble (omnibus schema).
#' Columns:
#' \itemize{
#'   \item \code{design}: Design label (\code{"nonparametric_within"}).
#'   \item \code{method}: Omnibus method label (\code{"friedman"}).
#'   \item \code{term}: Name of the within-subject factor (the tested term).
#'   \item \code{statistic_type}: Type of omnibus test statistic (\code{"chi_squared"}).
#'   \item \code{statistic}: Omnibus test statistic value (\eqn{\chi^2}).
#'   \item \code{df1}: Always \code{NA} for this test (reserved column).
#'   \item \code{df2}: Always \code{NA} for this test (reserved column).
#'   \item \code{df}: Degrees of freedom for the Friedman test (k - 1).
#'   \item \code{p}: Omnibus p-value.
#'   \item \code{effect_size}: Effect size value (Kendall's W).
#'   \item \code{effect_size_type}: Effect size label (\code{"kendalls_W"}).
#'   \item \code{missing}: Missing-data strategy used (\code{"pairwise"} or \code{"complete_case"}).
#'   \item \code{n_obs}: Number of rows (observations) used after row-wise filtering (before enforcing complete blocks).
#'   \item \code{n_id}: Number of complete IDs (blocks) used in the Friedman test after enforcing complete cases.
#' }
#'
#' @examples
#' dat = datasets::CO2 %>%
#'   dplyr::filter(Type == "Quebec", Treatment == "chilled") %>%
#'   dplyr::mutate(conc = factor(conc))
#'
#' friedman_omnibus(dat, dv = uptake, within = conc, id = Plant)
#'
#' @export
#' @importFrom magrittr %>%
friedman_omnibus = function(data, dv, within, id, missing = c("pairwise", "complete_case")) {
  missing = match.arg(missing)

  dv_sym = rlang::ensym(dv)
  within_sym = rlang::ensym(within)
  id_sym = rlang::ensym(id)

  dv_nm = rlang::as_name(dv_sym)
  within_nm = rlang::as_name(within_sym)
  id_nm = rlang::as_name(id_sym)

  dat0 =
    data %>%
    dplyr::filter(!is.na(.data[[dv_nm]]), !is.na(.data[[within_nm]]), !is.na(.data[[id_nm]]))

  dat0 = .coerce_factor_cols(dat0, c(within_nm, id_nm))

  tab = base::table(dat0[[id_nm]], dat0[[within_nm]])
  if (any(tab > 1L)) {
    stop("Friedman omnibus requires an unreplicated complete block design: found duplicated ID-within cells.", call. = FALSE)
  }

  mat = base::tapply(dat0[[dv_nm]], list(dat0[[id_nm]], dat0[[within_nm]]), mean)
  mat = mat[stats::complete.cases(mat), , drop = FALSE]

  n_id = nrow(mat)
  k = ncol(mat)
  if (is.null(n_id) || n_id < 2L || k < 2L) {
    stop("Not enough complete blocks to run Friedman test after enforcing complete cases.", call. = FALSE)
  }

  ft = stats::friedman.test(mat)

  chi2 = unname(ft$statistic[[1]])
  df_val = unname(ft$parameter[[1]])
  p_val = unname(ft$p.value)

  W = chi2 / (n_id * (k - 1))

  tibble::tibble(
    design = "nonparametric_within",
    method = "friedman",
    term = within_nm,
    statistic_type = "chi_squared",
    statistic = as.numeric(chi2),
    df1 = NA_real_,
    df2 = NA_real_,
    df = as.numeric(df_val),
    p = as.numeric(p_val),
    effect_size = as.numeric(W),
    effect_size_type = "kendalls_W",
    missing = missing,
    n_obs = as.numeric(nrow(dat0)),
    n_id = as.numeric(n_id)
  )
}


#' Kruskal-Wallis test (omnibus) with eta-squared (H-based)
#'
#' Performs the omnibus Kruskal-Wallis test using \code{rstatix::kruskal_test()} and computes
#' an H-based eta-squared using \code{rstatix::kruskal_effsize()}.
#'
#' @details
#' **Test**: Kruskal-Wallis. \cr
#' **Statistic**: H with df = k - 1. \cr
#' **Effect size**: H-based eta-squared from \code{rstatix::kruskal_effsize()}.
#'
#' \strong{Missing-data handling (\code{missing})}
#' \itemize{
#'   \item \code{"pairwise"} and \code{"complete_case"} behave identically here:
#'   rows with missing DV/group are removed.
#' }
#'
#' @param data A data.frame.
#' @param dv Dependent variable (numeric). Unquoted name or string.
#' @param group Grouping factor. Unquoted name or string.
#' @param missing Missing-data strategy: \code{"pairwise"} (default) or \code{"complete_case"}.
#'
#' @return
#' A one-row tibble (omnibus schema).
#' Columns:
#' \itemize{
#'   \item \code{design}: Design label (\code{"nonparametric_between"}).
#'   \item \code{method}: Omnibus method label (\code{"kruskal_wallis"}).
#'   \item \code{term}: Name of the grouping factor (the tested term).
#'   \item \code{statistic_type}: Type of omnibus test statistic (\code{"H"}).
#'   \item \code{statistic}: Omnibus test statistic value (H).
#'   \item \code{df1}: Always \code{NA} for this test (reserved column).
#'   \item \code{df2}: Always \code{NA} for this test (reserved column).
#'   \item \code{df}: Degrees of freedom for the Kruskal-Wallis test (k - 1).
#'   \item \code{p}: Omnibus p-value.
#'   \item \code{effect_size}: Effect size value (eta-squared based on H).
#'   \item \code{effect_size_type}: Effect size label (\code{"eta2_H"}).
#'   \item \code{missing}: Missing-data strategy used (\code{"pairwise"} or \code{"complete_case"}).
#'   \item \code{n_obs}: Number of rows (observations) used in the analysis after filtering for missingness.
#'   \item \code{n_id}: Always \code{NA} for between-subject designs in this package (reserved column).
#' }
#'
#' @examples
#' kruskal_omnibus(datasets::PlantGrowth, dv = weight, group = group)
#'
#' @export
#' @importFrom magrittr %>%
kruskal_omnibus = function(data, dv, group, missing = c("pairwise", "complete_case")) {
  missing = match.arg(missing)

  dv_sym = rlang::ensym(dv)
  group_sym = rlang::ensym(group)

  dv_nm = rlang::as_name(dv_sym)
  group_nm = rlang::as_name(group_sym)

  dat0 =
    data %>%
    dplyr::filter(!is.na(.data[[dv_nm]]), !is.na(.data[[group_nm]]))

  dat0 = .coerce_factor_cols(dat0, c(group_nm))

  kt =
    rstatix::kruskal_test(
      data = dat0,
      formula = stats::as.formula(paste(dv_nm, "~", group_nm))
    ) %>%
    tibble::as_tibble()

  es =
    rstatix::kruskal_effsize(
      data = dat0,
      formula = stats::as.formula(paste(dv_nm, "~", group_nm))
    ) %>%
    tibble::as_tibble()

  H_val = .get_num1(kt, c("statistic", "H"))
  df_val = .get_num1(kt, c("df"))
  p_val = .get_num1(kt, c("p", "p.value"))

  eff = .get_num1(es, c("effsize", "eta2", "eta2_h", "eta2H"))
  if (is.na(eff)) {
    num_cols = names(es)[vapply(es, is.numeric, logical(1))]
    if (length(num_cols) >= 1) eff = as.numeric(es[[num_cols[1]]][1])
  }

  tibble::tibble(
    design = "nonparametric_between",
    method = "kruskal_wallis",
    term = group_nm,
    statistic_type = "H",
    statistic = H_val,
    df1 = NA_real_,
    df2 = NA_real_,
    df = df_val,
    p = p_val,
    effect_size = eff,
    effect_size_type = "eta2_H",
    missing = missing,
    n_obs = as.numeric(nrow(dat0)),
    n_id = NA_real_
  )
}


###############################################################################
# Post-hoc tests (4 functions)
###############################################################################

#' Repeated-measures post-hoc comparisons (paired t) with dz and CIs
#'
#' Performs pairwise post-hoc comparisons for a one-factor repeated-measures design
#' using paired t-tests (\code{rstatix::pairwise_t_test()}), and adds Cohen's \eqn{d_z}
#' with both unadjusted and multiplicity-adjusted confidence intervals.
#'
#' @details
#' **Post-hoc test**: paired t-tests. \cr
#' **Effect size**: \eqn{d_z} via \code{effectsize::repeated_measures_d(method = "z")}. \cr
#' **P-value adjustment**: \code{p_adjust_method} (default: \code{"holm"}). \cr
#' **CI adjustment**: controlled by \code{ci_adjust_method} (independent of p-value adjustment);
#' supported: \code{"none"}, \code{"bonferroni"}, \code{"sidak"}.
#'
#' \strong{Important note on independence of p-adjustment and CI-adjustment}:
#' P-value adjustment (\code{p_adjust_method}) and CI adjustment (\code{ci_adjust_method}) are independent.
#' Therefore, it is possible to obtain results where \code{p_adj < 0.05} (significant after p-adjustment)
#' while the adjusted CI (\code{ci_low_adj}, \code{ci_high_adj}) still includes 0 (or vice versa),
#' especially when using different adjustment strategies (e.g., \code{p_adjust_method="holm"} and
#' \code{ci_adjust_method="none"}). This is not a software bug but a consequence of using different
#' error-control strategies for p-values and confidence intervals.
#'
#' \strong{Missing-data handling (\code{missing})}
#' \itemize{
#'   \item \code{"pairwise"} (default):
#'   \itemize{
#'     \item P-values are computed after removing rows with missing DV/within/id.
#'     \item Effect sizes and CIs are computed pairwise: for each comparison, only IDs with non-missing values
#'     in the two relevant within-levels are used.
#'   }
#'   \item \code{"complete_case"}: Restricts to IDs complete across all within-levels before analysis.
#' }
#'
#' \strong{Duplicated ID x within cells}
#' If duplicated ID x within cells exist, set \code{duplicate} to aggregate values within each cell.
#' Aggregation is reported via a single warning (not via return columns).
#'
#' \strong{Meaning of \code{duplicate_na_rm} for \code{duplicate="first"}}:
#' If \code{duplicate="first"} and \code{duplicate_na_rm=TRUE}, NAs within a duplicated cell are removed
#' before taking the first value (i.e., the first non-missing value is used).
#'
#' @param data A data.frame in long format.
#' @param dv Dependent variable (numeric). Unquoted name or string.
#' @param within Within-subject factor. Unquoted name or string.
#' @param id Subject identifier. Unquoted name or string.
#' @param conf_level Confidence level for unadjusted effect-size CIs (default: 0.95).
#' @param p_adjust_method P-value adjustment method passed to \code{rstatix} (default: \code{"holm"}).
#' @param ci_adjust_method CI multiplicity adjustment method (default: \code{"bonferroni"}).
#' @param missing Missing-data strategy: \code{"pairwise"} (default) or \code{"complete_case"}.
#' @param duplicate How to handle duplicated ID x within cells: \code{"mean"} (default), \code{"median"}, \code{"first"}, \code{"error"}.
#' @param duplicate_na_rm Logical; if TRUE, removes NAs within duplicated cells before aggregation (default: TRUE).
#'
#' @return
#' A tibble with one row per pairwise comparison.
#' Columns:
#' \itemize{
#'   \item \code{design}: Design label (\code{"parametric_within"}).
#'   \item \code{method}: Post-hoc method label (\code{"paired_t"}).
#'   \item \code{term}: Name of the within-factor.
#'   \item \code{group1}, \code{group2}: Level names for the pairwise comparison.
#'   \item \code{statistic_type}: Type of test statistic (\code{"t"}).
#'   \item \code{statistic}: Test statistic value.
#'   \item \code{df}: Degrees of freedom (t-test).
#'   \item \code{p}: Unadjusted p-value for the pairwise test.
#'   \item \code{p_adj}: Adjusted p-value (according to \code{p_adjust_method}).
#'   \item \code{p_adj_mark}: Visual marker for significance based on \code{p_adj}: \code{"*"} if \code{p_adj < 0.05},
#'   otherwise \code{""} (including when \code{p_adj} is \code{NA}).
#'   \item \code{p_adjust_method}: The p-value adjustment method used.
#'   \item \code{effect_size}: Effect size estimate (\eqn{d_z}).
#'   \item \code{effect_size_type}: Effect size label (\code{"dz"}).
#'   \item \code{ci_level}: Confidence level used for the unadjusted effect-size CI.
#'   \item \code{ci_low}, \code{ci_high}: Lower/upper bounds of the unadjusted effect-size CI.
#'   \item \code{ci_adjust_method}: Method used to adjust effect-size CIs (\code{"none"}, \code{"bonferroni"}, \code{"sidak"}).
#'   \item \code{ci_level_adj}: Adjusted CI level implied by \code{ci_adjust_method} and the number of comparisons.
#'   \item \code{ci_low_adj}, \code{ci_high_adj}: Lower/upper bounds of the adjusted effect-size CI.
#'   \item \code{n_pair}: Number of paired observations used for effect size/CIs (complete pairs for the two levels).
#'   \item \code{missing}: Missing-data strategy used.
#' }
#'
#' @examples
#' dat = datasets::CO2 %>%
#'   dplyr::filter(Type == "Quebec", Treatment == "chilled") %>%
#'   dplyr::mutate(conc = factor(conc))
#'
#' rm_posthoc(
#'   dat, dv = uptake, within = conc, id = Plant,
#'   p_adjust_method = "holm",
#'   ci_adjust_method = "bonferroni",
#'   missing = "pairwise"
#' )
#'
#' @export
#' @importFrom magrittr %>%
rm_posthoc = function(
    data, dv, within, id,
    conf_level = 0.95,
    p_adjust_method = "holm",
    ci_adjust_method = c("bonferroni", "none", "sidak"),
    missing = c("pairwise", "complete_case"),
    duplicate = c("mean", "median", "first", "error"),
    duplicate_na_rm = TRUE
) {
  missing = match.arg(missing)
  ci_adjust_method = match.arg(ci_adjust_method)
  duplicate = match.arg(duplicate)

  dv_sym = rlang::ensym(dv)
  within_sym = rlang::ensym(within)
  id_sym = rlang::ensym(id)

  dv_nm = rlang::as_name(dv_sym)
  within_nm = rlang::as_name(within_sym)
  id_nm = rlang::as_name(id_sym)

  dat0 =
    data %>%
    dplyr::filter(!is.na(.data[[dv_nm]]), !is.na(.data[[within_nm]]), !is.na(.data[[id_nm]]))

  dat0 = .coerce_factor_cols(dat0, c(within_nm, id_nm))

  if (missing == "complete_case") {
    complete_ids = .complete_ids_within(dat0, id_nm = id_nm, within_nm = within_nm, dv_nm = dv_nm)
    dat0 = dat0 %>% dplyr::semi_join(complete_ids, by = stats::setNames(id_nm, id_nm))
  }

  ttab =
    rstatix::pairwise_t_test(
      data = dat0,
      formula = stats::as.formula(paste(dv_nm, "~", within_nm)),
      paired = TRUE,
      detailed = TRUE,
      p.adjust.method = p_adjust_method
    ) %>%
    tibble::as_tibble()

  m = nrow(ttab)
  if (m == 0) return(.empty_posthoc_tibble())

  adj = .ci_adjust_info(conf_level, ci_adjust_method, m)
  conf_level_adj = adj$ci_level_adj

  stat = .get_num_vec(ttab, c("statistic"), m)
  df_val = .get_num_vec(ttab, c("df"), m)
  p_raw = .get_num_vec(ttab, c("p"), m)
  p_adj = .get_num_vec(ttab, c("p.adj", "p_adj", "p.adjusted", "padj"), m)

  es = rep(NA_real_, m)
  ci_low = rep(NA_real_, m)
  ci_high = rep(NA_real_, m)
  ci_low_adj = rep(NA_real_, m)
  ci_high_adj = rep(NA_real_, m)
  n_pair = rep(NA_real_, m)

  warn_dup = character(0)
  warn_skip = character(0)

  # Specific -> general (avoid overly broad partial matches)
  es_regex = c("^d_z$", "^dz$", "^Cohens_d$", "Cohens_d", "^d$")

  for (i in seq_len(m)) {
    g1 = as.character(ttab$group1[i])
    g2 = as.character(ttab$group2[i])

    sub =
      dat0 %>%
      dplyr::filter(.data[[within_nm]] %in% c(g1, g2)) %>%
      dplyr::filter(!is.na(.data[[dv_nm]]))

    wide = .build_wide_mat(
      sub, id_nm = id_nm, within_nm = within_nm, dv_nm = dv_nm,
      duplicate = duplicate, duplicate_na_rm = duplicate_na_rm
    )
    if (length(wide$agg_warning) > 0) warn_dup = c(warn_dup, wide$agg_warning)

    mat = wide$mat
    if (!all(c(g1, g2) %in% colnames(mat))) {
      warn_skip = c(warn_skip, paste0(g1, " vs ", g2, ": missing level after wide conversion"))
      next
    }

    x = mat[, g1]
    y = mat[, g2]
    ok = stats::complete.cases(x, y)
    x = x[ok]
    y = y[ok]
    n_pair[i] = length(x)

    if (length(x) < 2L) {
      warn_skip = c(warn_skip, paste0(g1, " vs ", g2, ": insufficient paired observations"))
      next
    }

    d1 = as.data.frame(effectsize::repeated_measures_d(x, y, method = "z", ci = conf_level))
    d2 = as.data.frame(effectsize::repeated_measures_d(x, y, method = "z", ci = conf_level_adj))

    e1 = .extract_es_ci(d1, es_regex = es_regex, warn_ctx = paste0("rm_posthoc: ", g1, " vs ", g2))
    e2 = .extract_es_ci(d2, es_regex = es_regex, warn_ctx = paste0("rm_posthoc(adj): ", g1, " vs ", g2))

    es[i] = e1$es
    ci_low[i] = e1$low
    ci_high[i] = e1$high
    ci_low_adj[i] = e2$low
    ci_high_adj[i] = e2$high
  }

  .warn_aggregate_once(warn_dup, header = "[rm_posthoc]")

  if (length(unique(warn_skip)) > 0) {
    .warn_aggregate_once(unique(warn_skip), header = "[rm_posthoc] Some comparisons could not compute effect sizes/CI:")
  }

  tibble::tibble(
    design = "parametric_within",
    method = "paired_t",
    term = within_nm,
    group1 = as.character(ttab$group1),
    group2 = as.character(ttab$group2),
    statistic_type = "t",
    statistic = stat,
    df = df_val,
    p = p_raw,
    p_adj = p_adj,
    p_adj_mark = ifelse(is.na(p_adj), "", ifelse(p_adj < 0.05, "*", "")),
    p_adjust_method = rep(p_adjust_method, m),
    effect_size = es,
    effect_size_type = "dz",
    ci_level = conf_level,
    ci_low = ci_low,
    ci_high = ci_high,
    ci_adjust_method = ci_adjust_method,
    ci_level_adj = conf_level_adj,
    ci_low_adj = ci_low_adj,
    ci_high_adj = ci_high_adj,
    n_pair = n_pair,
    missing = missing
  )
}


#' Welch-type post-hoc comparisons (pairwise Welch t) with Hedges' g and CIs
#'
#' Performs pairwise Welch t-tests (\code{pool.sd = FALSE}) via \code{rstatix::pairwise_t_test()},
#' and adds Hedges' g (unpooled SD) with unadjusted and adjusted confidence intervals.
#'
#' @details
#' **Post-hoc test**: Welch t-tests (unpooled SD). \cr
#' **Effect size**: Hedges' g via \code{effectsize::hedges_g(pooled_sd = FALSE)}. \cr
#' **P-value adjustment**: \code{p_adjust_method} (default: \code{"holm"}). \cr
#' **CI adjustment**: controlled by \code{ci_adjust_method} (independent of p-value adjustment);
#' supported: \code{"none"}, \code{"bonferroni"}, \code{"sidak"}.
#'
#' \strong{Important note on independence of p-adjustment and CI-adjustment}:
#' See \code{\link{rm_posthoc}} for a detailed explanation. The same note applies here.
#'
#' \strong{Missing-data handling (\code{missing})}
#' \itemize{
#'   \item \code{"pairwise"} and \code{"complete_case"} behave identically here:
#'   rows with missing DV/group are removed.
#' }
#'
#' @param data A data.frame.
#' @param dv Dependent variable (numeric). Unquoted name or string.
#' @param group Grouping factor. Unquoted name or string.
#' @param conf_level Confidence level for unadjusted effect-size CIs (default: 0.95).
#' @param p_adjust_method P-value adjustment method passed to \code{rstatix} (default: \code{"holm"}).
#' @param ci_adjust_method CI multiplicity adjustment method (default: \code{"bonferroni"}).
#' @param missing Missing-data strategy: \code{"pairwise"} (default) or \code{"complete_case"}.
#'
#' @return
#' A tibble with one row per pairwise comparison.
#' Columns:
#' \itemize{
#'   \item \code{design}: Design label (\code{"parametric_between_unequal_var"}).
#'   \item \code{method}: Post-hoc method label (\code{"welch_pairwise_t"}).
#'   \item \code{term}: Name of the grouping factor.
#'   \item \code{group1}, \code{group2}: Level names for the pairwise comparison.
#'   \item \code{statistic_type}: Type of test statistic (\code{"t"}).
#'   \item \code{statistic}: Test statistic value.
#'   \item \code{df}: Degrees of freedom (Welch t-test).
#'   \item \code{p}: Unadjusted p-value.
#'   \item \code{p_adj}: Adjusted p-value (according to \code{p_adjust_method}).
#'   \item \code{p_adj_mark}: \code{"*"} if \code{p_adj < 0.05}, otherwise \code{""} (including \code{NA}).
#'   \item \code{p_adjust_method}: The p-value adjustment method used.
#'   \item \code{effect_size}: Effect size estimate (Hedges' g, unpooled SD).
#'   \item \code{effect_size_type}: Effect size label (\code{"hedges_g"}).
#'   \item \code{ci_level}: Confidence level used for the unadjusted effect-size CI.
#'   \item \code{ci_low}, \code{ci_high}: Lower/upper bounds of the unadjusted effect-size CI.
#'   \item \code{ci_adjust_method}: Method used to adjust effect-size CIs.
#'   \item \code{ci_level_adj}: Adjusted CI level implied by \code{ci_adjust_method} and the number of comparisons.
#'   \item \code{ci_low_adj}, \code{ci_high_adj}: Lower/upper bounds of the adjusted effect-size CI.
#'   \item \code{n_pair}: Total sample size used for effect size/CIs (\code{n_x + n_y}).
#'   \item \code{missing}: Missing-data strategy used.
#' }
#'
#' @examples
#' welch_posthoc(
#'   datasets::PlantGrowth, dv = weight, group = group,
#'   p_adjust_method = "holm",
#'   ci_adjust_method = "bonferroni"
#' )
#'
#' @export
#' @importFrom magrittr %>%
welch_posthoc = function(
    data, dv, group,
    conf_level = 0.95,
    p_adjust_method = "holm",
    ci_adjust_method = c("bonferroni", "none", "sidak"),
    missing = c("pairwise", "complete_case")
) {
  missing = match.arg(missing)
  ci_adjust_method = match.arg(ci_adjust_method)

  dv_sym = rlang::ensym(dv)
  group_sym = rlang::ensym(group)

  dv_nm = rlang::as_name(dv_sym)
  group_nm = rlang::as_name(group_sym)

  dat0 =
    data %>%
    dplyr::filter(!is.na(.data[[dv_nm]]), !is.na(.data[[group_nm]]))

  dat0 = .coerce_factor_cols(dat0, c(group_nm))

  ttab =
    rstatix::pairwise_t_test(
      data = dat0,
      formula = stats::as.formula(paste(dv_nm, "~", group_nm)),
      paired = FALSE,
      pool.sd = FALSE,
      detailed = TRUE,
      p.adjust.method = p_adjust_method
    ) %>%
    tibble::as_tibble()

  m = nrow(ttab)
  if (m == 0) return(.empty_posthoc_tibble())

  adj = .ci_adjust_info(conf_level, ci_adjust_method, m)
  conf_level_adj = adj$ci_level_adj

  stat = .get_num_vec(ttab, c("statistic"), m)
  df_val = .get_num_vec(ttab, c("df"), m)
  p_raw = .get_num_vec(ttab, c("p"), m)
  p_adj = .get_num_vec(ttab, c("p.adj", "p_adj", "p.adjusted", "padj"), m)

  es = rep(NA_real_, m)
  ci_low = rep(NA_real_, m)
  ci_high = rep(NA_real_, m)
  ci_low_adj = rep(NA_real_, m)
  ci_high_adj = rep(NA_real_, m)
  n_pair = rep(NA_real_, m)

  warn_skip = character(0)

  es_regex = c("^Hedges_g$", "Hedges[_\\.]*g", "^g$")

  for (i in seq_len(m)) {
    g1 = as.character(ttab$group1[i])
    g2 = as.character(ttab$group2[i])

    x = dat0[[dv_nm]][dat0[[group_nm]] == g1]
    y = dat0[[dv_nm]][dat0[[group_nm]] == g2]
    x = x[!is.na(x)]
    y = y[!is.na(y)]
    n_pair[i] = length(x) + length(y)

    if (length(x) < 2L || length(y) < 2L) {
      warn_skip = c(warn_skip, paste0(g1, " vs ", g2, ": insufficient observations"))
      next
    }

    d1 = as.data.frame(effectsize::hedges_g(x, y, pooled_sd = FALSE, ci = conf_level))
    d2 = as.data.frame(effectsize::hedges_g(x, y, pooled_sd = FALSE, ci = conf_level_adj))

    e1 = .extract_es_ci(d1, es_regex = es_regex, warn_ctx = paste0("welch_posthoc: ", g1, " vs ", g2))
    e2 = .extract_es_ci(d2, es_regex = es_regex, warn_ctx = paste0("welch_posthoc(adj): ", g1, " vs ", g2))

    es[i] = e1$es
    ci_low[i] = e1$low
    ci_high[i] = e1$high
    ci_low_adj[i] = e2$low
    ci_high_adj[i] = e2$high
  }

  if (length(unique(warn_skip)) > 0) {
    .warn_aggregate_once(unique(warn_skip), header = "[welch_posthoc] Some comparisons could not compute effect sizes/CI:")
  }

  tibble::tibble(
    design = "parametric_between_unequal_var",
    method = "welch_pairwise_t",
    term = group_nm,
    group1 = as.character(ttab$group1),
    group2 = as.character(ttab$group2),
    statistic_type = "t",
    statistic = stat,
    df = df_val,
    p = p_raw,
    p_adj = p_adj,
    p_adj_mark = ifelse(is.na(p_adj), "", ifelse(p_adj < 0.05, "*", "")),
    p_adjust_method = rep(p_adjust_method, m),
    effect_size = es,
    effect_size_type = "hedges_g",
    ci_level = conf_level,
    ci_low = ci_low,
    ci_high = ci_high,
    ci_adjust_method = ci_adjust_method,
    ci_level_adj = conf_level_adj,
    ci_low_adj = ci_low_adj,
    ci_high_adj = ci_high_adj,
    n_pair = n_pair,
    missing = missing
  )
}


#' Friedman post-hoc comparisons (paired Wilcoxon) with rank-biserial and CIs
#'
#' Performs pairwise post-hoc comparisons for a one-factor within-subject design using
#' paired Wilcoxon signed-rank tests via \code{rstatix::pairwise_wilcox_test()}.
#' Also computes rank-biserial correlation (paired) with unadjusted and adjusted CIs.
#'
#' @details
#' **Post-hoc test**: paired Wilcoxon signed-rank tests. \cr
#' **Effect size**: rank-biserial via \code{effectsize::rank_biserial(paired = TRUE)}. \cr
#' **P-value adjustment**: \code{p_adjust_method} (default: \code{"holm"}). \cr
#' **CI adjustment**: controlled by \code{ci_adjust_method} (independent of p-value adjustment);
#' supported: \code{"none"}, \code{"bonferroni"}, \code{"sidak"}.
#'
#' \strong{Important note on independence of p-adjustment and CI-adjustment}:
#' See \code{\link{rm_posthoc}} for a detailed explanation. The same note applies here.
#'
#' \strong{Missing-data handling (\code{missing})}
#' \itemize{
#'   \item \code{"pairwise"} (default):
#'   p-values are computed after removing rows with missing DV/within/id; effect sizes/CIs use only IDs
#'   complete for the two compared within-levels.
#'   \item \code{"complete_case"}: Restricts to IDs complete across all within-levels.
#' }
#'
#' \strong{Duplicated ID x within cells}
#' If duplicated ID x within cells exist, set \code{duplicate} to aggregate values within each cell.
#' Aggregation is reported via a single warning (not via return columns).
#'
#' \strong{Meaning of \code{duplicate_na_rm} for \code{duplicate="first"}}:
#' If \code{duplicate="first"} and \code{duplicate_na_rm=TRUE}, NAs within a duplicated cell are removed
#' before taking the first value (i.e., the first non-missing value is used).
#'
#' @param data A data.frame in long format.
#' @param dv Dependent variable (numeric). Unquoted name or string.
#' @param within Within-subject factor. Unquoted name or string.
#' @param id Subject identifier. Unquoted name or string.
#' @param conf_level Confidence level for unadjusted effect-size CIs (default: 0.95).
#' @param p_adjust_method P-value adjustment method passed to \code{rstatix} (default: \code{"holm"}).
#' @param ci_adjust_method CI multiplicity adjustment method (default: \code{"bonferroni"}).
#' @param missing Missing-data strategy: \code{"pairwise"} (default) or \code{"complete_case"}.
#' @param duplicate How to handle duplicated ID x within cells: \code{"mean"} (default), \code{"median"}, \code{"first"}, \code{"error"}.
#' @param duplicate_na_rm Logical; if TRUE, removes NAs within duplicated cells before aggregation (default: TRUE).
#'
#' @return
#' A tibble with one row per pairwise comparison.
#' Columns:
#' \itemize{
#'   \item \code{design}: Design label (\code{"nonparametric_within"}).
#'   \item \code{method}: Post-hoc method label (\code{"paired_wilcox"}).
#'   \item \code{term}: Name of the within-factor.
#'   \item \code{group1}, \code{group2}: Level names for the pairwise comparison.
#'   \item \code{statistic_type}: Type of test statistic (\code{"W"}).
#'   \item \code{statistic}: Wilcoxon test statistic (may be \code{NA} depending on rstatix output/version).
#'   \item \code{df}: Always \code{NA} for this post-hoc method in this package.
#'   \item \code{p}: Unadjusted p-value.
#'   \item \code{p_adj}: Adjusted p-value (according to \code{p_adjust_method}).
#'   \item \code{p_adj_mark}: \code{"*"} if \code{p_adj < 0.05}, otherwise \code{""} (including \code{NA}).
#'   \item \code{p_adjust_method}: The p-value adjustment method used.
#'   \item \code{effect_size}: Effect size estimate (rank-biserial correlation; paired).
#'   \item \code{effect_size_type}: Effect size label (\code{"rank_biserial"}).
#'   \item \code{ci_level}: Confidence level used for the unadjusted effect-size CI.
#'   \item \code{ci_low}, \code{ci_high}: Lower/upper bounds of the unadjusted effect-size CI.
#'   \item \code{ci_adjust_method}: Method used to adjust effect-size CIs.
#'   \item \code{ci_level_adj}: Adjusted CI level implied by \code{ci_adjust_method} and the number of comparisons.
#'   \item \code{ci_low_adj}, \code{ci_high_adj}: Lower/upper bounds of the adjusted effect-size CI.
#'   \item \code{n_pair}: Number of paired observations used for effect size/CIs (complete pairs for the two levels).
#'   \item \code{missing}: Missing-data strategy used.
#' }
#'
#' @examples
#' dat = datasets::CO2 %>%
#'   dplyr::mutate(
#'     Plant = factor(Plant),
#'     conc = factor(conc)
#'   )
#'
#' friedman_posthoc(
#'   dat,
#'   dv = uptake,
#'   within = conc,
#'   id = Plant,
#'   p_adjust_method = "holm",
#'   ci_adjust_method = "bonferroni",
#'   missing = "complete_case"
#' )
#'
#' @export
#' @importFrom magrittr %>%
friedman_posthoc = function(
    data, dv, within, id,
    conf_level = 0.95,
    p_adjust_method = "holm",
    ci_adjust_method = c("bonferroni", "none", "sidak"),
    missing = c("pairwise", "complete_case"),
    duplicate = c("mean", "median", "first", "error"),
    duplicate_na_rm = TRUE
) {
  missing = match.arg(missing)
  ci_adjust_method = match.arg(ci_adjust_method)
  duplicate = match.arg(duplicate)

  dv_sym = rlang::ensym(dv)
  within_sym = rlang::ensym(within)
  id_sym = rlang::ensym(id)

  dv_nm = rlang::as_name(dv_sym)
  within_nm = rlang::as_name(within_sym)
  id_nm = rlang::as_name(id_sym)

  dat0 =
    data %>%
    dplyr::filter(!is.na(.data[[dv_nm]]), !is.na(.data[[within_nm]]), !is.na(.data[[id_nm]]))

  dat0 = .coerce_factor_cols(dat0, c(within_nm, id_nm))

  if (missing == "complete_case") {
    complete_ids = .complete_ids_within(dat0, id_nm = id_nm, within_nm = within_nm, dv_nm = dv_nm)
    dat0 = dat0 %>% dplyr::semi_join(complete_ids, by = stats::setNames(id_nm, id_nm))
  }

  wtab =
    rstatix::pairwise_wilcox_test(
      data = dat0,
      formula = stats::as.formula(paste(dv_nm, "~", within_nm)),
      paired = TRUE,
      detailed = TRUE,
      p.adjust.method = p_adjust_method
    ) %>%
    tibble::as_tibble()

  m = nrow(wtab)
  if (m == 0) return(.empty_posthoc_tibble())

  adj = .ci_adjust_info(conf_level, ci_adjust_method, m)
  conf_level_adj = adj$ci_level_adj

  p_raw = .get_num_vec(wtab, c("p"), m)
  p_adj = .get_num_vec(wtab, c("p.adj", "p_adj", "p.adjusted", "padj"), m)
  stat = .get_num_vec(wtab, c("statistic"), m)

  es = rep(NA_real_, m)
  ci_low = rep(NA_real_, m)
  ci_high = rep(NA_real_, m)
  ci_low_adj = rep(NA_real_, m)
  ci_high_adj = rep(NA_real_, m)
  n_pair = rep(NA_real_, m)

  warn_dup = character(0)
  warn_skip = character(0)

  es_regex = c("^rank_biserial$", "rank[_\\.]*biserial", "biserial", "^r$")

  for (i in seq_len(m)) {
    g1 = as.character(wtab$group1[i])
    g2 = as.character(wtab$group2[i])

    sub =
      dat0 %>%
      dplyr::filter(.data[[within_nm]] %in% c(g1, g2)) %>%
      dplyr::filter(!is.na(.data[[dv_nm]]))

    wide = .build_wide_mat(
      sub, id_nm = id_nm, within_nm = within_nm, dv_nm = dv_nm,
      duplicate = duplicate, duplicate_na_rm = duplicate_na_rm
    )
    if (length(wide$agg_warning) > 0) warn_dup = c(warn_dup, wide$agg_warning)

    mat = wide$mat
    if (!all(c(g1, g2) %in% colnames(mat))) {
      warn_skip = c(warn_skip, paste0(g1, " vs ", g2, ": missing level after wide conversion"))
      next
    }

    x = mat[, g1]
    y = mat[, g2]
    ok = stats::complete.cases(x, y)
    x = x[ok]
    y = y[ok]
    n_pair[i] = length(x)

    if (length(x) < 2L) {
      warn_skip = c(warn_skip, paste0(g1, " vs ", g2, ": insufficient paired observations"))
      next
    }

    r1 = as.data.frame(effectsize::rank_biserial(x, y, paired = TRUE, ci = conf_level))
    r2 = as.data.frame(effectsize::rank_biserial(x, y, paired = TRUE, ci = conf_level_adj))

    e1 = .extract_es_ci(r1, es_regex = es_regex, warn_ctx = paste0("friedman_posthoc: ", g1, " vs ", g2))
    e2 = .extract_es_ci(r2, es_regex = es_regex, warn_ctx = paste0("friedman_posthoc(adj): ", g1, " vs ", g2))

    es[i] = e1$es
    ci_low[i] = e1$low
    ci_high[i] = e1$high
    ci_low_adj[i] = e2$low
    ci_high_adj[i] = e2$high
  }

  .warn_aggregate_once(warn_dup, header = "[friedman_posthoc]")

  if (length(unique(warn_skip)) > 0) {
    .warn_aggregate_once(unique(warn_skip), header = "[friedman_posthoc] Some comparisons could not compute effect sizes/CI:")
  }

  tibble::tibble(
    design = "nonparametric_within",
    method = "paired_wilcox",
    term = within_nm,
    group1 = as.character(wtab$group1),
    group2 = as.character(wtab$group2),
    statistic_type = "W",
    statistic = stat,
    df = NA_real_,
    p = p_raw,
    p_adj = p_adj,
    p_adj_mark = ifelse(is.na(p_adj), "", ifelse(p_adj < 0.05, "*", "")),
    p_adjust_method = rep(p_adjust_method, m),
    effect_size = es,
    effect_size_type = "rank_biserial",
    ci_level = conf_level,
    ci_low = ci_low,
    ci_high = ci_high,
    ci_adjust_method = ci_adjust_method,
    ci_level_adj = conf_level_adj,
    ci_low_adj = ci_low_adj,
    ci_high_adj = ci_high_adj,
    n_pair = n_pair,
    missing = missing
  )
}


#' Kruskal-Wallis post-hoc comparisons (Dunn) with rank-biserial and CIs
#'
#' Performs pairwise post-hoc comparisons using Dunn tests via \code{rstatix::dunn_test()},
#' and adds rank-biserial correlation (unpaired) with unadjusted and adjusted confidence intervals.
#'
#' @details
#' **Post-hoc test**: Dunn tests. \cr
#' **Effect size**: rank-biserial via \code{effectsize::rank_biserial(paired = FALSE)}. \cr
#' **P-value adjustment**: \code{p_adjust_method} (default: \code{"holm"}). \cr
#' **CI adjustment**: controlled by \code{ci_adjust_method} (independent of p-value adjustment);
#' supported: \code{"none"}, \code{"bonferroni"}, \code{"sidak"}.
#'
#' \strong{Important note on independence of p-adjustment and CI-adjustment}:
#' See \code{\link{rm_posthoc}} for a detailed explanation. The same note applies here.
#'
#' \strong{Missing-data handling (\code{missing})}
#' \itemize{
#'   \item \code{"pairwise"} and \code{"complete_case"} behave identically here:
#'   rows with missing DV/group are removed.
#' }
#'
#' @param data A data.frame.
#' @param dv Dependent variable (numeric). Unquoted name or string.
#' @param group Grouping factor. Unquoted name or string.
#' @param conf_level Confidence level for unadjusted effect-size CIs (default: 0.95).
#' @param p_adjust_method P-value adjustment method passed to \code{rstatix} (default: \code{"holm"}).
#' @param ci_adjust_method CI multiplicity adjustment method (default: \code{"bonferroni"}).
#' @param missing Missing-data strategy: \code{"pairwise"} (default) or \code{"complete_case"}.
#'
#' @return
#' A tibble with one row per pairwise comparison.
#' Columns:
#' \itemize{
#'   \item \code{design}: Design label (\code{"nonparametric_between"}).
#'   \item \code{method}: Post-hoc method label (\code{"dunn"}).
#'   \item \code{term}: Name of the grouping factor.
#'   \item \code{group1}, \code{group2}: Level names for the pairwise comparison.
#'   \item \code{statistic_type}: Type of test statistic (\code{"Z"}).
#'   \item \code{statistic}: Z statistic value (Dunn test).
#'   \item \code{df}: Always \code{NA} for this post-hoc method in this package.
#'   \item \code{p}: Unadjusted p-value.
#'   \item \code{p_adj}: Adjusted p-value (according to \code{p_adjust_method}).
#'   \item \code{p_adj_mark}: \code{"*"} if \code{p_adj < 0.05}, otherwise \code{""} (including \code{NA}).
#'   \item \code{p_adjust_method}: The p-value adjustment method used.
#'   \item \code{effect_size}: Effect size estimate (rank-biserial correlation; unpaired).
#'   \item \code{effect_size_type}: Effect size label (\code{"rank_biserial"}).
#'   \item \code{ci_level}: Confidence level used for the unadjusted effect-size CI.
#'   \item \code{ci_low}, \code{ci_high}: Lower/upper bounds of the unadjusted effect-size CI.
#'   \item \code{ci_adjust_method}: Method used to adjust effect-size CIs.
#'   \item \code{ci_level_adj}: Adjusted CI level implied by \code{ci_adjust_method} and the number of comparisons.
#'   \item \code{ci_low_adj}, \code{ci_high_adj}: Lower/upper bounds of the adjusted effect-size CI.
#'   \item \code{n_pair}: Total sample size used for effect size/CIs (\code{n_x + n_y}).
#'   \item \code{missing}: Missing-data strategy used.
#' }
#'
#' @examples
#' kruskal_posthoc(
#'   datasets::PlantGrowth, dv = weight, group = group,
#'   p_adjust_method = "holm",
#'   ci_adjust_method = "bonferroni"
#' )
#'
#' @export
#' @importFrom magrittr %>%
kruskal_posthoc = function(
    data, dv, group,
    conf_level = 0.95,
    p_adjust_method = "holm",
    ci_adjust_method = c("bonferroni", "none", "sidak"),
    missing = c("pairwise", "complete_case")
) {
  missing = match.arg(missing)
  ci_adjust_method = match.arg(ci_adjust_method)

  dv_sym = rlang::ensym(dv)
  group_sym = rlang::ensym(group)

  dv_nm = rlang::as_name(dv_sym)
  group_nm = rlang::as_name(group_sym)

  dat0 =
    data %>%
    dplyr::filter(!is.na(.data[[dv_nm]]), !is.na(.data[[group_nm]]))

  dat0 = .coerce_factor_cols(dat0, c(group_nm))

  dtab =
    rstatix::dunn_test(
      data = dat0,
      formula = stats::as.formula(paste(dv_nm, "~", group_nm)),
      detailed = TRUE,
      p.adjust.method = p_adjust_method
    ) %>%
    tibble::as_tibble()

  m = nrow(dtab)
  if (m == 0) return(.empty_posthoc_tibble())

  adj = .ci_adjust_info(conf_level, ci_adjust_method, m)
  conf_level_adj = adj$ci_level_adj

  p_raw = .get_num_vec(dtab, c("p"), m)
  p_adj = .get_num_vec(dtab, c("p.adj", "p_adj", "p.adjusted", "padj"), m)
  stat = .get_num_vec(dtab, c("Z", "statistic"), m)

  es = rep(NA_real_, m)
  ci_low = rep(NA_real_, m)
  ci_high = rep(NA_real_, m)
  ci_low_adj = rep(NA_real_, m)
  ci_high_adj = rep(NA_real_, m)
  n_pair = rep(NA_real_, m)

  warn_skip = character(0)

  es_regex = c("^rank_biserial$", "rank[_\\.]*biserial", "biserial", "^r$")

  for (i in seq_len(m)) {
    g1 = as.character(dtab$group1[i])
    g2 = as.character(dtab$group2[i])

    x = dat0[[dv_nm]][dat0[[group_nm]] == g1]
    y = dat0[[dv_nm]][dat0[[group_nm]] == g2]
    x = x[!is.na(x)]
    y = y[!is.na(y)]
    n_pair[i] = length(x) + length(y)

    if (length(x) < 2L || length(y) < 2L) {
      warn_skip = c(warn_skip, paste0(g1, " vs ", g2, ": insufficient observations"))
      next
    }

    r1 = as.data.frame(effectsize::rank_biserial(x, y, paired = FALSE, ci = conf_level))
    r2 = as.data.frame(effectsize::rank_biserial(x, y, paired = FALSE, ci = conf_level_adj))

    e1 = .extract_es_ci(r1, es_regex = es_regex, warn_ctx = paste0("kruskal_posthoc: ", g1, " vs ", g2))
    e2 = .extract_es_ci(r2, es_regex = es_regex, warn_ctx = paste0("kruskal_posthoc(adj): ", g1, " vs ", g2))

    es[i] = e1$es
    ci_low[i] = e1$low
    ci_high[i] = e1$high
    ci_low_adj[i] = e2$low
    ci_high_adj[i] = e2$high
  }

  if (length(unique(warn_skip)) > 0) {
    .warn_aggregate_once(unique(warn_skip), header = "[kruskal_posthoc] Some comparisons could not compute effect sizes/CI:")
  }

  tibble::tibble(
    design = "nonparametric_between",
    method = "dunn",
    term = group_nm,
    group1 = as.character(dtab$group1),
    group2 = as.character(dtab$group2),
    statistic_type = "Z",
    statistic = stat,
    df = NA_real_,
    p = p_raw,
    p_adj = p_adj,
    p_adj_mark = ifelse(is.na(p_adj), "", ifelse(p_adj < 0.05, "*", "")),
    p_adjust_method = rep(p_adjust_method, m),
    effect_size = es,
    effect_size_type = "rank_biserial",
    ci_level = conf_level,
    ci_low = ci_low,
    ci_high = ci_high,
    ci_adjust_method = ci_adjust_method,
    ci_level_adj = conf_level_adj,
    ci_low_adj = ci_low_adj,
    ci_high_adj = ci_high_adj,
    n_pair = n_pair,
    missing = missing
  )
}

