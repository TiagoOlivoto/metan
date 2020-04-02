# metan 1.4.0.9000
## New functions
* `select_rows_na()` and `select_cols_na()` to select rows or columns with with `NA` values
* `mgidi()` to compute the multi-trait genotype-ideotype distance index

## Minor changes
* Remove dependency on dendextend
* Update package site with [pkgdown v1.5.0](https://pkgdown.r-lib.org/news/index.html).
* Update documentation in `ge_plot()`
* Allow using `fai_blup()` with `gamem()`
* Improve checking process with `inspect()`
* Improve feedback for results, indicating random and fixed effects. Thanks to [@NelsonJunior](https://scholar.google.com.br/citations?user=i2F6X04AAAAJ&hl=pt-BR) for his suggestion.

## Bug fixes
* `get_model_data()` now fills rows that don't matches across columns with `NA`
*

# metan 1.4.0
## Bug fixes
* Factor columns can now have custom names rather than `ENV`, `GEN`, and `REP` only ([#2](https://github.com/TiagoOlivoto/metan/issues/2)).

## New functions
* `gmd()` a shortcut to `get_model_data()`
* `gtb()` to generate a genotype-by-trait biplot.
* `gamem_met()` to analyze genotypes in multi-environment trials using mixed- or random-effect models allowing unbalanced data. Thanks to [@EderOliveira](https://www.embrapa.br/en/web/portal/team/-/empregado/321725/eder-jorge-de-oliveira) for his e-mail.
* `has_class()` to check if a class exists.
* `impute_missing_val()` to impute missing values in a two-way table based on Expectation-Maximization algoritms. 
* `non_collinear_vars()` to select a set of predictors with minimal multicollinearity.
* `replace_na()` to replace `NA` values quicly.
* `random_na()` to generate random `NA` values based on a desired proportion.


## Minor changes
* `gge()`, `performs_ammi()`, `waas()`, and `waasb()` now handle with unbalanced data by implementing a low-rank matrix approximation using singular value decomposition to impute missing entires. Imputation generates a warning message.
* `NA` values are checked and removed with a warning when computing stability indexes. Thanks to [@MdFarhad](https://www.researchgate.net/profile/Md_Farhad) for alerting me.
* New argument `plot_res` in `path_coeff()` to create a residual plot of the multiple regression model.
* Update the citation file to include the [published official reference](https://doi.org/10.1111/2041-210X.13384).
* Argument `verbose` deprecated in functions `anova_ind()` and `split_factors()`
* Argument `rep` deprecated in functions `Fox()`, `Huehn()`, `superiority()`, and `Thennarasu()`.
* Deprecated argument `means_by` removed in functions `can_corr()` and `clustering()`.
* Deprecated argument `verbose` removed in functions `colindiag()` and `split_factors()`.
* Deprecated argument `values` removed in functions `desc_stat()` and `find_outliers()`.
* Deprecated argument `var` removed in function `desc_wider()`.
* Remove dependency on lattice by using ggplot2 in `plot.resp_surf()`.
* An up-to-date cheat sheet was included.



# metan 1.3.0

## New functions
   - `alpha_color()` To get a semi-transparent color
   - `gafem()` To analyze genotypes using fixed-effect models.
   - `residual_plots()` A helper function to create residuals plots.
   - `stars_pval()` To generate significance stars from p-values
   - `doo()` An alternative to `dplyr::do` for doing anything
   
### utils_stats
   - `cv_by()` For computing coefficient of variation by levels of a factor.
   - `max_by()` For computing maximum values by levels of a factor.
   - `means_by()` For computing arithmetic means by levels of a factor.
   - `min_by()` For computing minimum values by levels of a factor.
   - `n_by()` For getting the length.
   - `sd_by()` For computing sample standard deviation.
   - `sem_by()` For computing standard error of the mean by levels of a factor.
   - `av_dev()` computes the average absolute deviation.
   - `ci_mean()` computes the confidence interval for the mean.
   - `cv()` computes the coefficient of variation.
   - `hm_mean()`, `gm_mean()` computes the harmonic and geometric means, respectively. The harmonic mean is the reciprocal of the arithmetic mean of the reciprocals. The geometric mean is the nth root of n products.
   - `kurt()` computes the kurtosis like used in SAS and SPSS.
   - `range_data()` Computes the range of the values.
   - `sd_amo()`, `sd_pop()` Computes sample and populational standard deviation, respectively.
   - `sem()` computes the standard error of the mean.
   - `skew()` computes the skewness like used in SAS and SPSS.
   - `sum_dev()` computes the sum of the absolute deviations.
   - `sum_sq_dev()` computes the sum of the squared deviations.
   - `var_amo()`, `var_pop()` computes sample and populational variance.
   - `valid_n()` Return the valid (not NA) length of a data.

### utils_rows_cols
   - `colnames_to_lower()`, `colnames_to_upper()`, and `colnames_to_title()` to translate column names to lower, upper and title cases quickly.

### utils_num_str
   - `all_lower_case()`, `all_upper_case()`, and `all_title_case()` to translate strings vectors or character columns of a data frame to lower, upper and title cases, respectively.
   - `tidy_strings()` Tidy up characters strings, non-numeric columns, or any selected columns in a data frame by putting all word in upper case, replacing any space, tabulation, punctuation characters by `'_'`, and putting `'_'` between lower and upper cases.
   - `find_text_in_num()` Find text fragments in columns assumed to be numeric. This is especially useful when `everything()` is used in argument `resp` to select the response variables.


## New arguments
   - `anova_ind()`, `anova_joint()`, `performs_ammi()`, `waas()` and `waasb()`,  now have the argument `block` to analyze data from trials conducted in an alpha-lattice design. Thanks to [@myaseen208](https://twitter.com/myaseen208?lang=en) for his suggestion regarding multi-environment trial analysis with alpha-lattice designs.
   - argument `repel` included in `plot_scores()` to control wheater the labels are repelled or not to avoid overlapping.

## Deprecated arguments
   Argument `means_by` was deprecated in functions `can_corr()` and `clustering()`. Use `means_by()` to pass data based on means of factor to these functions.
   
## Minor changes

   - Change "#000000FF" with "#FFFFFF00" in `transparent_color()`
   - `desc_stat()` now handles grouped data passed from `dplyr::group_by()`
   - `plot_scores()` now support objects of class `waas_mean`.
   - Include inst/CITATION to return a reference paper with `citation("metan")`.
   - Change 'PC2' with 'PC1' in y-axis of `plot_scores(type = 2)` ([#1](https://github.com/TiagoOlivoto/metan/issues/1))
   - `get_model_data()` now support models of class `anova_joint` and `gafem` and extract random effects of models fitted with `waasb()` and `gamem()`.
   - Update `plot.waasb()` and `plot.gamem()` to show distribution of random effects.
   - `inspect()`, `cv_blup()`, `cv_ammif()`, and `cv_ammi()` now generate a warning message saying that is not possible to compute cross-validation procedures in experiments with two replicates only. Thanks to [@Vlatko](https://www.researchgate.net/profile/Vlatko_Galic2) for his email.
   - `plot.wsmp()` now returns heatmaps created with ggplot2. Thus, we removed dependency on `gplots`.
   - Vignettes updated

# metan 1.2.1
* References describing the methods implemented in the package were included in description field of DESCRIPTION file as suggested by the CRAN team.

# metan 1.2.0
* Minor changes
   * `corr_plot()` now don't write a warning message to the console by default.
   * `select_numeric_cols()` now is used as a helper function in `metan`.
   * `metan` now reexports `mutate()` from `dplyr` package.
   * `get_model_data()` now set default values for each class of models.
   * Argument `by` that calls internally `split_factors()` included to facilitate the application of the functions to each level of one grouping variable.

* New functions
   * `add_cols()`, and `add_rows()` for adding columns and rows, respectively.
   * `remove_cols()`, and `remove_rows()` for removing columns and rows, respectively.
   * `select_cols()` and `select_rows()` for selecting columns and rows, respectively.
   * `select_numeric_cols()`, and `select_non_numeric_cols()` for selecting numeric and non-numeric variables quickly.
   * `round_cols()` for rounding a whole data frame to significant figures.
   * `all_lower_case()`, and `all_upper_case()` for handling with cases.
   * `extract_number()`, `extract_string()`, `remove_strings()`, `replace_number()`, and  `replace_string()`, for handling with numbers and strings.
   * `get_level_size()`, and `get_levels()` for getting size of levels and levels of a factor.
   * `means_by()` for computing means by one or more factors quickly.
   * `ge_means()` for computing genotype-environment means
   * `ge_winners()` for getting winner genotypes or ranking genotypes within environments.
   * `env_dissimilarity()` for computing dissimilarity between test environments.

# metan 1.1.2

* Reexport select_helpers `starts_with()`, `ends_with()`, `contains()`, `contains()`, `num_range()`, `one_of()`, `everything()`, and `last_col()`.
* When possible, argument `resp` (response variable(s) now support select helpers.
* New helper function `sem()` for computing standard error of mean.
* New helper functions `remove_rows_na()` and `remove_cols_na()` for removing rows or columns with `NA` values quickly.
* New select helpers `difference_var()`, `intersect_var()`, and 
`union_var()` for selecting variables that match an expression.
* New function `Schmildt()` for stability analysis.
* Plot regression slope and mean performance in objects of class `ge_reg`.
* Update `get_model_data()` to support objects of class `Schmildt`and `Annicchiarico`.


# metan 1.1.1

* Now `on.exit()` is used in S3 generic functions `print()` to ensure that the settings are reset when a function is excited.
* Computationally intensive parts in vignettes uses pre-computed results.

# metan 1.1.0
I'm very pleased to announce the release of `metan` 1.1.0, This is a minor release with bug fixes and new functions. The most important changes are described below.

* New function `corr_stab_ind()` for computing Spearman's rank correlation between stability indexes;
* New function `corr_coef()` for computing correlation coefficients and p-values;
* New S3 method `plot.corr_coef()` for creating correlation heat maps;
* New S3 method `print.corr_coef()` for printing correlation and p-values;
* New helper functions `make_lower_tri()` and `make_upper_tri()` for creating lower and upper triangular matrices, respectively.
* New helper function `reorder_cormat()` for reordering a correlation matrix according to the correlation coefficients;
* Improve usability of `get_model_data()` by supporting new classes of models. Now, `get_model_data()` can be used to get all statistics or ranks computed with the wrapper function `ge_stats()`.
* `arrange_ggplot()` now support objects of class `ggmatrix`.
* Change the default plot theme to `theme_metan()`
* Update function's documentation;
* Update vignettes.

# metan 1.0.2
* New function `arrange_ggplot()` for arranging ggplot2 graphics;
* New function `ge_effects()` for computing genotype-environment effects;
* New function `gai()` for computing the geometric adaptability index;
* New helper function `gm_mean()` for computing geometric mean;
* New helper function `hm_mean()` for computing harmonic mean;
* New helper function `Huehn()` for computing Huehn's stability statistic;
* New helper function `Thennasaru()` for computing Thennasaru's stability statistic;
* Improve usability of `get_model_data()` by supporting new classes of models;
* Update function's documentation;
* Update vignettes;

# metan 1.0.1
* New function `gamem()` for analyzing genotypes in one-way trials using mixed-effect models;
* New function `desc_wider()` to convert an output of the function `desc_stat()` to a 'wide' format;
* New function `Fox()` for Stability analysis;
* New function `Shukla()` for stability analysis;
* New function `to_factor()` to quickly convert variables to factors;
* Improve usability of `get_model_data()` function;
* Update function's documentation;
* Update vignettes;

# metan 1.0.0
The changes in this version were made based on suggestions received when metan was submitted to CRAN for the first time.

## Major changes
The documentation of the following functions was updated by including/updating the \\value section of .Rd files.

* `AMMI_indexes()`
* `Annichiarico()`
* `anova_ind()`
* `as.lpcor()`
* `as.split_factors()`
* `bind_cv()`
* `can_cor()`
* `comb_vars()`
* `corr_ci()`
* `corr_plot()`
* `covcor_design()`
* `cv_ammi()`
* `cv_ammif()`
* `cv_blup()`
* `desc_stat()`
* `ecovalence()`
* `fai_blup()`
* `ge_factanal()`
* `ge_plot()`
* `ge_reg()`
* `ge_stats()`
* `get_model_data()`
* `is.lpcorr()`
* `is.split_factors()`
* `mahala()`
* `mahala_design()`
* `make_mat()`
* `make_sym()`
* `mtsi()`
* `pairs_mantel()`
* `plot.*()` and `plot_*()` functions
* `rbind_fill()`
* `resca()`
* `resp_surf()` 
* `waas()`
* `wsmp()`
* `waasb()`
      
## Minor changes

To allow automatic testing, the examples of the following functions were unwrapped by deleting \\dontrun{}.

* `bind_cv()`
* `clustering()`
* `comb_vars()`
* `corr_ci()`
* `corr_plot()`
* `covcor_design()`
* `desc_stat()`
* `ecovalence()`
* `path_coefff()`
* `plot.fai_blup()`
* `plot.mtsi()`
* `plot.wsmp()`
* `plot_ci()`
* `wsmp()`

In the examples of the functions for cross-validation \\dontrun{} was changed with \\donttest{}
* `cv_ammi()`
* `cv_ammif()`
* `cv_blup()`
* `plot.cv_ammif()`
      
# metan 0.2.0
This is the first version that will be submitted to CRAN. In this version, deprecated functions in the last versions were defunct. Some new features were implemented.

* New functions
   * `fai_blup()` computes the FAI-BLUP index (https://onlinelibrary.wiley.com/doi/full/10.1111/gcbb.12443)
   * `gge()` computes the genotype plus genotype-vs-environment model.
   * `plot_factbars()` and `plot_factlines()` are used to create bar and line plots, respectively, considering an one- or two-way experimental design.
   * `desc_stat()` computes several descriptive statistics.
   * `can_corr()`computes canonical correlation coefficients.
   * `resp_surf()` computes response surface model using two quantitative factors.
   * `make_mat()` is used to create a two-way table using two columns (factors) and one response variable. 
   * `make_sym()` is used to create a symmetric matrix using a upper- or lower-diagonal matrix. 
   
* Minor improvements
   * New evaluation for text vectors are now used in the functions `AMMI_indexes()` and `fai_blup()` and `desc_stat()`. For example, to indicate the statistics to be computed in `desc_stat()` you must use now ` stats = c("mean, SE.mean, CV, max, min"))` instead  `stats = c("mean", "SE.mean", "CV", "max", "min"))`

# metan 0.1.5
In the latest development version, the package **METAAB** was renamed to **metan** (**m**ulti-**e**nvironment **t**rials **an**alysis). Aiming at a cleaner coding, in this version, some functions were deprecated and will be defunct in the near future. Alternative functions were implemented.

* For `WAAS.AMMI()`, use `waas()`.
* For `WAASBYratio()`, use `wsmp()`.
* For `WAASratio.AMMI()`, use `wsmp()`.
* For `autoplot.WAAS.AMMI()`, use `autoplot.waas()`.
* For `plot.WAASBYratio()`, use `plot.wsmp()`.
* For `plot.WAASratio.AMMI()`, use `plot.wsmp()`.
* For `predict.WAAS.AMMI()`, use `predict.waas()`.
* For `summary.WAAS.AMMI()`, use `summary.waas()`

Widely-known parametric and nonparametric methods were implemented, using the following functions.

* `Annicchiarico()` to compute the genotypic confidence index.
* `ecovalence()` to compute the Wricke's ecovalence.
* `ge_factanal()` to compute to compute the stability and environmental.
* `ge_reg()` to compute the joint-regression analysis.
stratification using factor analysis.
* `superiority()` to compute the nonparametric superiority index.


# METAAB 0.1.4
In the latest development version, some useful functions were included. One of the most interesting features included in this version was allowing the functions to receive data from the forward-pipe operator %>%. Bellow are the functions included in this version.

* `anova_ind()` to perform a within-environment analysis of variance easily;
* `colindiag()` to perform a collinearity diagnostic of a set of predictors;a
* `find_outliers()` to easily find possible outliers in the dataset;
* `group_factors()` to split a dataset into a list of subsets using one or more grouping factors. This function may be used befor some functions, e.g., `find_outliers()`, `colindiag()`, `path_coeff()` to compute the statistics for each level of the factor (or combination of levels of factors).
* `lpcor()` to compute linear and partial correlation coefficients.
* `pairs_mantel()` to compute a graphic pairwise Mantel's test using a set of correlation matrices;
* `path_coeff()` to compute path coefficients with minimal multicollinearity;

The following S3 Methods were also implemented:

* `is.group_factors()` and `as.group_factors()` to check or easily coerce a dataframe that has one or more factor columns to an object of `group_factors`;
* `is.lpcorr()` and `as.lpcorr()`  to check or easily coerce a list of correlation matrices to an object of `lpcorr`;


# METAAB 0.1.3
* AMMI-based stability indexes;
* Allow analyzing multiple variables at the same time;
* S3 methods such as `plot()`, `predict()`, `summary()` implemented.

# METAAB 0.1.2

* Mixed-effect model with environment random effect;
* Random-effect model.

# METAAB 0.1.1

* The first version of the package
