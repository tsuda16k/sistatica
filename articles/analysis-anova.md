# analysis-anova

``` r
library(sistatica)
```

## このビネットについて

このビネットでは **sistatica**
の分散分析関連の関数の基本的な使い方を解説します。

------------------------------------------------------------------------

## 反復測定分散分析（rm-ANOVA）

この節では **sistatica**
を用いて、1要因の反復測定分散分析（repeated-measures
ANOVA）を実行する方法を解説します。  
基本的には次の2つの関数を使います。

- [`rmOmni()`](https://tsuda16k.github.io/sistatica/reference/rmOmni.md)：オムニバステスト（条件の主効果）
- [`rmPost()`](https://tsuda16k.github.io/sistatica/reference/rmPost.md)：事後検定（対応ありのペア比較：paired
  t-test）

### 前提：データ形式（long形式）

[`rmOmni()`](https://tsuda16k.github.io/sistatica/reference/rmOmni.md)
と
[`rmPost()`](https://tsuda16k.github.io/sistatica/reference/rmPost.md)
は **long形式**のデータフレームを想定します。必要な列は次の3つです。

- `id`：被験者ID（同一被験者を識別できる列）
- `group`：条件（within-subject factor）
- `y`：従属変数（数値）

このビネットでは、まずワイド形式で疑似データを作成し、long形式に変換して解析します。

``` r
library(sistatica)
library(dplyr)
library(tidyr)

set.seed(1)

n = 24
u = rnorm(n, 0, 4)

w = data.frame(
  id = factor(seq_len(n)),
  A = 50 + u + rnorm(n, 0, 5),
  B = 55 + u + rnorm(n, 0, 5),
  C = 60 + u + rnorm(n, 0, 5)
)

dat = pivot_longer(w, cols = c(A, B, C), names_to = "group", values_to = "y") %>%
  mutate(group = factor(group))
```

------------------------------------------------------------------------

## 1) オムニバステスト：`rmOmni()`

[`rmOmni()`](https://tsuda16k.github.io/sistatica/reference/rmOmni.md)
は、1要因の反復測定分散分析（rm-ANOVA）の
**オムニバステスト（条件の主効果）** を実行します。  
結果は「手法間で共通化された列名」を持つ **1行の tibble**
として返されます。

### 実行例

``` r
omni = rmOmni(dat, id = "id", group = "group", y = "y", correction = "auto")
```

球面性補正を常に Greenhouse–Geisser（GG）にしたい場合は
`correction = "GG"` を指定します。

``` r
omni_gg = rmOmni(dat, id = "id", group = "group", y = "y", correction = "GG")
```

### 引数の説明

- `df`  
  long形式データフレーム（`id`・`group`・`y` を含む）

- `id`（デフォルト `"id"`）  
  被験者ID列名（文字列）

- `group`（デフォルト `"group"`）  
  条件列名（文字列、within-subject factor）

- `y`（デフォルト `"y"`）  
  従属変数列名（文字列、数値）

- `correction`（デフォルト `"auto"`）  
  球面性補正（sphericity correction）の指定（`"auto"`, `"GG"`, `"HF"`,
  `"none"`）

  - `"auto"`：球面性が破れている場合のみGG補正（推奨のデフォルト運用）
  - `"GG"`：常にGreenhouse–Geisser補正（適用可能な場合）
  - `"HF"`：常にHuynh–Feldt補正（適用可能な場合）
  - `"none"`：補正なし（未補正の自由度・p値）

### 出力（`omni`）の列の意味

[`rmOmni()`](https://tsuda16k.github.io/sistatica/reference/rmOmni.md)
の出力は tibble で、主に以下の列を含みます。

- `analysis_method`：手法識別子（`"repeated_measures_anova"`）
- `effect_term`：効果の名前（スキーマ統一のため常に `"group"`）
- `test_statistic`：検定統計量（F値）
- `test_statistic_type`：統計量の種類（常に `"F"`）
- `df1`：自由度（分子、補正される場合あり）
- `df2`：自由度（分母、補正される場合あり）
- `p_value`：p値（補正される場合あり）
- `effect_size`：効果量（一般化η²：GES）
- `effect_size_type`：効果量の種類（`"generalized_eta_squared"`）
- `effect_size_scale`：効果量のスケール（`"variance_explained"`）
- `effect_size_ci_low`, `effect_size_ci_high`：効果量CI（現状は予約列で
  `NA`）

------------------------------------------------------------------------

## 2) 事後検定：`rmPost()`

[`rmPost()`](https://tsuda16k.github.io/sistatica/reference/rmPost.md)
は、条件間のペア比較（対応ありの t 検定）を実行し、  
**多重比較補正後 p 値**、**差の推定量と信頼区間**、**効果量**、さらに
**条件ごとの平均・SD・サンプルサイズ**
など、論文報告に有用な列をまとめて返します。

### 実行例（p値：Holm補正、CI：調整なし）

``` r
post = rmPost(
  dat,
  id = "id", group = "group", y = "y",
  pAdjust = "holm",
  ciAdjust = "none"
)
```

### 実行例（p値：Holm補正、CI：Bonferroni調整）

``` r
post_bonf_ci = rmPost(
  dat,
  id = "id", group = "group", y = "y",
  pAdjust = "holm",
  ciAdjust = "bonferroni"
)
```

### 引数の説明

- `df`  
  long形式データフレーム（`id`・`group`・`y` を含む）

- `id`（デフォルト `"id"`）  
  被験者ID列名（文字列）

- `group`（デフォルト `"group"`）  
  条件列名（文字列、within-subject factor）

- `y`（デフォルト `"y"`）  
  従属変数列名（文字列、数値）

- `alpha`（デフォルト `0.05`）  
  ファミリー水準の有意水準（CI調整に用いる）

- `pAdjust`（デフォルト `"holm"`）  
  p値の多重比較補正法。[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)
  が受け付ける文字列（例：`"holm"`, `"bonferroni"`, `"BH"` など）

- `ciAdjust`（`"none"` または `"bonferroni"`）  
  信頼区間の調整法。`"bonferroni"` の場合、比較数に基づいて CI
  の信頼水準を調整します。

### 出力（`post`）の列の意味（要点）

[`rmPost()`](https://tsuda16k.github.io/sistatica/reference/rmPost.md)
は **ペア比較ごとに1行**の tibble を返します。主な列は以下です。

#### 比較の識別

- `analysis_method`：手法識別子（`"repeated_measures_anova"`）
- `effect_term`：効果の名前（常に `"group"`）
- `comparison_kind`：比較の種類（常に `"pairwise"`）
- `contrast_label`：表示用ラベル（例：`"A vs B"`）
- `contrast_id`：機械可読なID（例：`"group:A|B"`）
- `group1`, `group2`：比較する条件名

#### 条件ごとの記述統計（報告向け）

- `center_type`：代表値の種類（常に `"mean"`）
- `group1_center`, `group2_center`：各条件の平均
- `dispersion_type`：ばらつきの種類（常に `"sd"`）
- `group1_dispersion`, `group2_dispersion`：各条件の標準偏差
- `n_group1`, `n_group2`：各条件の有効サンプルサイズ
- `n_total`：その比較で利用できた「ペア」の数（欠測があると減る可能性あり）

#### 差の推定・信頼区間

- `difference_estimate`：平均との差（`group1 - group2`）
- `ci_low`, `ci_high`：差の信頼区間
- `ci_level`：信頼水準（`ciAdjust = "bonferroni"` の場合は調整後）
- `ci_adjustment`：CI調整の有無（`"unadjusted"` /
  `"bonferroni_adjusted"`）

#### 検定結果（t検定）

- `test_statistic`：t値
- `df1`：自由度
- `p_value`：未補正p値
- `p_value_adjusted`：多重比較補正後p値
- `p_value_adjustment`：補正法ラベル（例：`"holm_adjusted"`）

#### 効果量

- `effect_size`：効果量（Cohen’s dz：対応ありの標準化平均差）
- `effect_size_type`：効果量の種類（`"cohens_dz"`）
- `effect_size_scale`：スケール（`"standardized_mean_difference"`）
- `effect_size_ci_low`, `effect_size_ci_high`：効果量CI（現状は予約列で
  `NA`）

------------------------------------------------------------------------

## レポート（表）作成に便利な抜き出し例

### オムニバステスト（rm-ANOVA）

``` r
omni %>%
  select(
    effect_term,
    test_statistic, df1, df2, p_value,
    effect_size, effect_size_type
  )
#> # A tibble: 1 × 7
#>   effect_term test_statistic   df1   df2  p_value effect_size effect_size_type  
#>   <chr>                <dbl> <dbl> <dbl>    <dbl>       <dbl> <chr>             
#> 1 group                 37.2     2    46 2.41e-10       0.363 generalized_eta_s…
```
