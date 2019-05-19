# metan 1.1.0
In the latest development version, the package **METAAB** was renammed to **metan** (**m**ulti-**e**nvironment **t**rials **an**alysis). Aiming at a cleaner coding, in this version, some functions were deprecated and will be defunct in the near future. Alternative functions were implemented.

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


# METAAB 1.0.3
In the latest development version, some useful functions were included. One of the most interesting features included in this version was allowing the functions to receive data from the forward-pipe operator %>%. Bellow are the functions included in this version.

* `anova_ind()` to perform a within-environment analysis of variance easely;
* `colindiag()` to perform a collinearity diagnostic of a set of predictors;
* `find_outliers()` to easily find possible outliers in the dataset;
* `group_factors()` to split a dataset into a list of subsets using one or more grouping factors. This function may be used befor some functions, e.g., `find_outliers()`, `colindiag()`, `path_coeff()` to compute the statistics for each level of the factor (or combination of levels of factors).
* `plcor()` to compute linear and partial correlation coefficients.
* `pairs_mantel()` to compute a graphic pairwise Mantel's test using a set of correlation matrices;
* `path_coeff()` to compute path coefficients with minimal multicollinearity;

The following S3 Methods were also implemented:

* `is.group_factors()` and `as.group_factors()` to check or easily coerce a dataframe that has one or more factor columns to an object of `group_factors`;
* `is.lpcorr()` and `as.lpcorr()`  to check or easily coerce a list of correlation matrice to an object of `lpcorr`;


# METAAB 1.0.2
* AMMI-based stability indexes;
* Allow analizing multiple variables at the same time;
* S3 methods such as `plot()`, `predict()`, `summary()` implemented.

# METAAB 1.0.1

* Mixed-effect model with environment random effect;
* Random-effect model.

# METAAB 1.0.0

* The first version of the package
