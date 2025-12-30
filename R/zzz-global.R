utils::globalVariables(c(
  # tidy-eval pronoun
  ".data",

  # internal standardized column names used inside functions
  "i","g","y",

  # rstatix::anova_test / get_anova_table() typical columns
  "Effect","F","DFn","DFd","p","ges",

  # rstatix omnibus outputs (welch/kw/friedman) commonly seen
  "term","statistic","df","df1","df2",

  # pairwise test outputs (t_test / wilcox_test / dunn_test / games_howell_test)
  "group1","group2","estimate","conf.low","conf.high",
  "p.adj","p.adj.signif","p.adj.method",

  # effect size outputs
  "effsize",

  # summary columns that might appear in pipelines
  "mean","sd","median","iqr","q1","q3",
  "n","se","tcrit","n_subjects","n_observations","n_complete_pairs",

  # standardized output columns created by sistatica functions
  "analysis_method","effect_term","comparison_kind","contrast_label","contrast_id",
  "center_type","group1_center","group2_center","dispersion_type","group1_dispersion","group2_dispersion",
  "n_group1","n_group2","n_total",
  "difference_estimate","difference_estimate_type","difference_direction",
  "ci_low","ci_high","ci_level","ci_adjustment",
  "test_statistic","test_statistic_type","p_value","p_value_adjusted","p_value_adjustment",
  "effect_size","effect_size_type","effect_size_scale","effect_size_ci_low","effect_size_ci_high",

  # extra intermediate names used in some versions/variants
  "p_value_wilcox"
))
