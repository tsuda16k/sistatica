# sistatica: 4種類の一要因デザインに対するオムニバス検定と事後検定

## 0. イントロダクション

このページでは、**sistatica**
パッケージを用いて1要因の分散分析を行う方法を説明します。

4種類の1要因分散分析デザインごとに、**オムニバス検定**と**事後検定（ペア比較）**を行う関数が用意されています。

### 関数一覧

| 手法                                         | オムニバス検定                                                                                   | 事後検定                                                                                   |
|----------------------------------------------|--------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------|
| パラメトリック・対応あり（反復測定）         | [`rm_anova_omnibus()`](https://tsuda16k.github.io/sistatica/reference/rm_anova_omnibus.md)       | [`rm_posthoc()`](https://tsuda16k.github.io/sistatica/reference/rm_posthoc.md)             |
| パラメトリック・対応なし（Welch）            | [`welch_anova_omnibus()`](https://tsuda16k.github.io/sistatica/reference/welch_anova_omnibus.md) | [`welch_posthoc()`](https://tsuda16k.github.io/sistatica/reference/welch_posthoc.md)       |
| 非パラメトリック・対応あり（Friedman）       | [`friedman_omnibus()`](https://tsuda16k.github.io/sistatica/reference/friedman_omnibus.md)       | [`friedman_posthoc()`](https://tsuda16k.github.io/sistatica/reference/friedman_posthoc.md) |
| 非パラメトリック・対応なし（Kruskal-Wallis） | [`kruskal_omnibus()`](https://tsuda16k.github.io/sistatica/reference/kruskal_omnibus.md)         | [`kruskal_posthoc()`](https://tsuda16k.github.io/sistatica/reference/kruskal_posthoc.md)   |

``` r
library(sistatica)
library(dplyr)
library(magrittr)
```

### 入力データ形式（ロング形式）

本パッケージは、全て **ロング形式**のデータを前提とします。

- 反復測定 / Friedman: `dv`（数値）, `within`（条件）, `id`（被験者ID）
- Welch / Kruskal-Wallis: `dv`（数値）, `group`（群）

### “有意”表示記号

事後検定の返り値には `p_adj_mark` 列が含まれます。

- `p_adj < 0.05` のとき `"*"`
- それ以外（および `p_adj` が `NA`）は `""`（空文字）

※ 有意水準 0.05 は固定仕様です。

### p値補正とCI補正の扱い

事後検定では、p値補正と効果量CI補正は **独立**して指定します。

- p値補正: `p_adjust_method`（推奨例: `"holm"`）
- CI補正: `ci_adjust_method`（推奨例: `"bonferroni"`、選択肢は `"none"`,
  `"bonferroni"`, `"sidak"`）

この独立性のため、例えば  
`p_adj < 0.05`（有意）だが `ci_low_adj`〜`ci_high_adj` が 0
を含む、のような
**不一致**が起こりえます。これはバグではなく、異なる誤差制御戦略を併用した帰結です。

本ビネットでは、推奨例として常に以下を用います。

- `p_adjust_method = "holm"`
- `ci_adjust_method = "bonferroni"`

なお、CI補正は **その関数が返すペア数（比較数）**に基づいて行われます。

### 欠損の扱い（`missing` 引数）

`missing` 引数は `"pairwise"`（デフォルト）と `"complete_case"`
から選びます。

- 反復測定 / Friedman（対応あり）
  - `"pairwise"`: 行単位で `dv/within/id` の欠損を除いた後、効果量とCIは
    **各ペアで両条件が揃うID**のみを用いて計算します。
  - `"complete_case"`:
    全条件にわたって欠損のないID（完全ブロック）のみで分析します。
- Welch / Kruskal-Wallis（対応なし）
  - 実質的にどちらも「欠損行を落とす」挙動になりやすいですが、APIとして引数は残しています。

### 重複セル（ID×within）の扱い（「対応あり」の事後検定を行う場合）

対応あり（[`rm_posthoc()`](https://tsuda16k.github.io/sistatica/reference/rm_posthoc.md),
[`friedman_posthoc()`](https://tsuda16k.github.io/sistatica/reference/friedman_posthoc.md)）では、`id × within`
セルが重複するデータに対し、次の選択肢があります。

- `duplicate = "mean"`（デフォルト）:
  各参加者（id）が各条件（within）で複数回測定されているデータは、いったん「参加者×条件」に1値へ集約（平均値を使用）してから、通常の一要因・反復測定として解析されます。
- `duplicate = median" / "first"`:
  中央値で集約する、あるいはデータ内で最初に登場する値を使う。
- `duplicate = "error"`（デフォルト）: 重複があれば停止

`duplicate = "first"` かつ `duplicate_na_rm = TRUE`
の場合は、重複セル内の NA を除いてから最初の値（first
non-NA）を取ります。

------------------------------------------------------------------------

## 1. パラメトリック・対応あり：反復測定分散分析

### 1.1 この手法での採用方法（本パッケージの設計）

- オムニバス:
  [`rm_anova_omnibus()`](https://tsuda16k.github.io/sistatica/reference/rm_anova_omnibus.md)（反復測定ANOVA）
  - 効果量: **GES（generalized eta-squared）**
- 事後:
  [`rm_posthoc()`](https://tsuda16k.github.io/sistatica/reference/rm_posthoc.md)（対応のある
  t 検定）
  - 効果量: **dz（Cohen’s d for paired / repeated measures）**
  - 効果量CI: `effectsize::repeated_measures_d(method="z")`
    による解析的CI
  - p値補正 / CI補正は独立（イントロ参照）

### 1.2 データ

ここでは [`datasets::CO2`](https://rdrr.io/r/datasets/zCO2.html)
を用います（ロング形式、`Plant` がID、`conc` が within）。

``` r
dat_rm = datasets::CO2 %>%
  dplyr::mutate(
    Plant = factor(Plant),
    conc = factor(conc)
  )

# 必要に応じて、対象列だけを残す（他の列があっても本関数は問題なく動作）
dat_rm_small = dat_rm %>%
  dplyr::select(Plant, conc, uptake)
```

### 1.3 オムニバス検定：`rm_anova_omnibus()`

``` r
rm_omni = rm_anova_omnibus(
  data = dat_rm_small,
  dv = uptake,
  within = conc,
  id = Plant,
  missing = "pairwise"
)
```

#### 出力の解釈（要点）

- `statistic_type` / `statistic` / `df1` / `df2` / `p` がオムニバス検定
- `effect_size_type = "ges"` と `effect_size` が GES
- `missing`, `n_obs`, `n_id` はデータ処理の要約

#### 結果データ内の主要情報の表示

``` r
rm_omni %>% dplyr::select(term, statistic, df1, df2, p, effect_size, effect_size_type, n_id, missing)
```

### 1.4 事後検定：`rm_posthoc()`

``` r
rm_post = rm_posthoc(
  data = dat_rm_small,
  dv = uptake,
  within = conc,
  id = Plant,
  p_adjust_method = "holm",
  ci_adjust_method = "bonferroni",
  missing = "pairwise",
  duplicate = "mean"
)
```

#### 出力の解釈（要点）

- `group1`, `group2` がペア
- `p` と `p_adj`、および `p_adj_mark`
- 効果量: `effect_size_type = "dz"` と `effect_size`
- CI:
  - 未補正: `ci_low`, `ci_high`（`ci_level`）
  - 補正: `ci_low_adj`, `ci_high_adj`（`ci_level_adj`,
    `ci_adjust_method`）
- `n_pair` はそのペアで効果量・CIに使えた対応データ数

#### 結果データ内の主要情報の表示

``` r
rm_post %>%
  dplyr::select(group1, group2, p_adj, p_adj_mark, effect_size, ci_low, ci_high, ci_low_adj, ci_high_adj, n_pair) %>%
  dplyr::arrange(p_adj)
```

------------------------------------------------------------------------

## 2. パラメトリック・対応なし：Welch分散分析

### 2.1 この手法での採用方法（本パッケージの設計）

- オムニバス:
  [`welch_anova_omnibus()`](https://tsuda16k.github.io/sistatica/reference/welch_anova_omnibus.md)（Welch一要因ANOVA）
  - 効果量: **omega²**（WelchのFと自由度から計算した近似値）
- 事後:
  [`welch_posthoc()`](https://tsuda16k.github.io/sistatica/reference/welch_posthoc.md)（Welch型のペア比較：`pool.sd = FALSE`）
  - 効果量: **Hedges’ g（pooled_sd = FALSE）**
  - 効果量CI: `effectsize::hedges_g(pooled_sd = FALSE)` の解析的CI
  - p値補正 / CI補正は独立（イントロ参照）

### 2.2 データ

ここでは
[`datasets::PlantGrowth`](https://rdrr.io/r/datasets/PlantGrowth.html)
を用います（ロング形式）。

``` r
dat_welch = datasets::PlantGrowth %>%
  dplyr::mutate(group = factor(group)) %>%
  dplyr::rename(y = weight)
```

### 2.3 オムニバス検定：`welch_anova_omnibus()`

``` r
welch_omni = welch_anova_omnibus(
  data = dat_welch,
  dv = y,
  group = group,
  missing = "pairwise"
)
```

#### 出力の解釈（要点）

- Welchの `F` と近似自由度（`df1`, `df2`）と `p`
- 効果量: `effect_size_type = "omega2"` と `effect_size`
  - Welch ANOVA に対する omega²
    は便宜的に提供される**近似**であり、厳密な等分散前提のomega²と同一視しない（パッケージの注記）。

#### 結果データ内の主要情報の表示

``` r
welch_omni %>% dplyr::select(term, statistic, df1, df2, p, effect_size, effect_size_type)
```

### 2.4 事後検定：`welch_posthoc()`

``` r
welch_post = welch_posthoc(
  data = dat_welch,
  dv = y,
  group = group,
  p_adjust_method = "holm",
  ci_adjust_method = "bonferroni",
  missing = "pairwise"
)
```

#### 出力の解釈（要点）

- `group1`, `group2` がペア
- `p_adj` と `p_adj_mark`
- 効果量: `effect_size_type = "hedges_g"` と `effect_size`
- CI（未補正／補正）と `ci_adjust_method`

#### 結果データ内の主要情報の表示

``` r
welch_post %>%
  dplyr::select(group1, group2, p_adj, p_adj_mark, effect_size, ci_low, ci_high, ci_low_adj, ci_high_adj) %>%
  dplyr::arrange(p_adj)
```

------------------------------------------------------------------------

## 3. 非パラメトリック・対応あり：フリードマン検定

### 3.1 この手法での採用方法（本パッケージの設計）

- オムニバス:
  [`friedman_omnibus()`](https://tsuda16k.github.io/sistatica/reference/friedman_omnibus.md)
  - 効果量: **Kendall’s W**（`W = chi^2 / (n*(k-1))`）
- 事後:
  [`friedman_posthoc()`](https://tsuda16k.github.io/sistatica/reference/friedman_posthoc.md)
  - 事後検定: **対応のある Wilcoxon（符号付順位）**のペア比較（rstatix）
  - 効果量: **rank-biserial correlation（paired = TRUE）**
  - 効果量CI: `effectsize::rank_biserial(paired=TRUE)` の解析的CI
  - p値補正 / CI補正は独立（イントロ参照）
- 重複セル（ID×within）:
  - オムニバスは unreplicated complete block を要求（重複は停止）
  - 事後は `duplicate` で集約可能

### 3.2 データ

反復測定と同様に [`datasets::CO2`](https://rdrr.io/r/datasets/zCO2.html)
を用い、within = `conc`, id = `Plant` とします。

``` r
dat_fr = datasets::CO2 %>%
  dplyr::mutate(
    Plant = factor(Plant),
    conc = factor(conc)
  ) %>%
  dplyr::select(Plant, conc, uptake)
```

### 3.3 オムニバス検定：`friedman_omnibus()`

``` r
fr_omni = friedman_omnibus(
  data = dat_fr,
  dv = uptake,
  within = conc,
  id = Plant,
  missing = "pairwise"
)
```

#### 出力の解釈（要点）

- `statistic_type = "chi_squared"` と `statistic`, `df`, `p`
- 効果量: `effect_size_type = "kendalls_W"` と `effect_size`

#### 結果データ内の主要情報の表示

``` r
fr_omni %>% dplyr::select(term, statistic, df, p, effect_size, effect_size_type)
```

### 3.4 事後検定：`friedman_posthoc()`

``` r
fr_post = friedman_posthoc(
  data = dat_fr,
  dv = uptake,
  within = conc,
  id = Plant,
  p_adjust_method = "holm",
  ci_adjust_method = "bonferroni",
  missing = "pairwise",
  duplicate = "mean"
)
```

#### 出力の解釈（要点）

- `p_adj` と `p_adj_mark`
- 効果量: `effect_size_type = "rank_biserial"` と `effect_size`
- CI（未補正／補正）と `n_pair`

#### 結果データ内の主要情報の表示

``` r
fr_post %>%
  dplyr::select(group1, group2, p_adj, p_adj_mark, effect_size, ci_low, ci_high, ci_low_adj, ci_high_adj, n_pair) %>%
  dplyr::arrange(p_adj)
```

------------------------------------------------------------------------

## 4. 非パラメトリック・対応なし：クラスカル・ウォリス検定

### 4.1 この手法での採用方法（本パッケージの設計）

- オムニバス:
  [`kruskal_omnibus()`](https://tsuda16k.github.io/sistatica/reference/kruskal_omnibus.md)
  - 効果量: **eta²（H-based）**（rstatix の `kruskal_effsize()` による）
- 事後:
  [`kruskal_posthoc()`](https://tsuda16k.github.io/sistatica/reference/kruskal_posthoc.md)
  - 事後検定: **Dunn検定**（rstatix）
  - 効果量: **rank-biserial correlation（paired = FALSE）**
  - 効果量CI: `effectsize::rank_biserial(paired=FALSE)` の解析的CI
  - p値補正 / CI補正は独立（イントロ参照）

### 4.2 データ

Welchと同様に
[`datasets::PlantGrowth`](https://rdrr.io/r/datasets/PlantGrowth.html)
を用います。

``` r
dat_kw = datasets::PlantGrowth %>%
  dplyr::mutate(group = factor(group)) %>%
  dplyr::rename(y = weight)
```

### 4.3 オムニバス検定：`kruskal_omnibus()`

``` r
kw_omni = kruskal_omnibus(
  data = dat_kw,
  dv = y,
  group = group,
  missing = "pairwise"
)
```

#### 出力の解釈（要点）

- `statistic_type = "H"` と `statistic`, `df`, `p`
- 効果量: `effect_size_type = "eta2_H"` と `effect_size`

#### 結果データ内の主要情報の表示

``` r
kw_omni %>% dplyr::select(term, statistic, df, p, effect_size, effect_size_type)
```

### 4.4 事後検定：`kruskal_posthoc()`

``` r
kw_post = kruskal_posthoc(
  data = dat_kw,
  dv = y,
  group = group,
  p_adjust_method = "holm",
  ci_adjust_method = "bonferroni",
  missing = "pairwise"
)
```

#### 出力の解釈（要点）

- `group1`, `group2`
- `p_adj` と `p_adj_mark`
- 効果量: `effect_size_type = "rank_biserial"` と `effect_size`
- CI（未補正／補正）と `ci_adjust_method`

#### 結果データ内の主要情報の表示

``` r
kw_post %>%
  dplyr::select(group1, group2, p_adj, p_adj_mark, effect_size, ci_low, ci_high, ci_low_adj, ci_high_adj) %>%
  dplyr::arrange(p_adj)
```
