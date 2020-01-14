
<!-- README.md is generated from README.Rmd. Please edit that file -->

# metan <img src="man/figures/logo.png" align="right" height=140/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/metan)](https://CRAN.R-project.org/package=metan)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![Downloads](http://cranlogs.r-pkg.org/badges/metan)](https://CRAN.R-project.org/package=metan)
<!-- badges: end -->

The metan (**m**ulti-**e**nvironment **t**rials **an**alysis) package
provides useful functions for analyzing multi-environment trial data
using parametric and nonparametric methods, including, but not limited
to:

  - Within-environment analysis of variance;
  - Estimation using AMMI considering different numbers of interaction
    principal component axes;
  - AMMI-based stability indexes;
  - GGE biplot analysis;
  - Prediction in mixed-effect models;
  - BLUP-based stability indexes;
  - Variance components and genetic parameters in mixed-effect models;
  - Cross-validation procedures for AMMI-family and BLUP models;
  - Graphics tools for generating biplots;
  - Parametric and nonparametric stability statistics

For more details see the [complete
vignette](https://tiagoolivoto.github.io/metan/).

# Installation

Install the released version of metan from CRAN with:

``` r
install.packages("metan")
```

Or install the development version from GitHub with:

``` r
devtools::install_github("TiagoOlivoto/metan")

# To build the HTML vignette use
devtools::install_github("TiagoOlivoto/metan", build_vignettes = TRUE)
```

*Note*: If you are a Windows user, you should also first download and
install the latest version of
[Rtools](https://cran.r-project.org/bin/windows/Rtools/).

For the latest release notes on this development version, see the [NEWS
file](https://tiagoolivoto.github.io/metan/news/index.html).

# Getting started

Here, we will use the example dataset `data_ge` that contains data on
two variables assessed in 10 genotypes growing in 14 environments. For
more details see `?data_ge`.

``` r
library(metan)
```

The first step is to inspect the data with the function `inspect()`.

``` r
inspect(data_ge, plot = TRUE)
#> # A tibble: 5 x 9
#>   Variable Class   Missing Levels Valid_n   Min Median   Max Outlier
#>   <chr>    <fct>   <fct>   <fct>    <int> <dbl>  <dbl> <dbl>   <dbl>
#> 1 ENV      factor  No      14         420 NA     NA    NA         NA
#> 2 GEN      factor  No      10         420 NA     NA    NA         NA
#> 3 REP      factor  No      3          420 NA     NA    NA         NA
#> 4 GY       numeric No      -          420  0.67   2.61  5.09       0
#> 5 HM       numeric No      -          420 38     48    58          0
#> No issues detected while inspecting data.
```

<img src="man/figures/README-unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

No issues while inspecting the data. Let’s continue with the analyzes\!

## AMMI model

### Fitting the model

The AMMI model is fitted with the function `performs_ammi()`. To analyze
multiple variables at once we can pass a comma-separated vector of
unquoted variable names, or use any select helpe. Here, using
`everything()` we apply the function to all numeric variables in the
data. For more details, see the [complete
vignette](https://tiagoolivoto.github.io/metan/articles/vignettes_ammi.html).

``` r
model <- performs_ammi(data_ge,
                       env = ENV,
                       gen = GEN,
                       rep = REP,
                       resp = everything(),
                       verbose = FALSE)
```

### Predicting the response variable

The S3 method `predict()` is implemented for objects of class
`performs_ammi` and may be used to estimate the response of each
genotype in each environment considering different number of Interaction
Principal Component Axis (IPCA). As a example, to predict the variables
GY and HM we will use four and six IPCA (number of significant IPCAs,
respectively).

``` r
pred <- predict(model, naxis = c(4, 6))
pred$GY
#> # A tibble: 140 x 8
#>    ENV   GEN       Y  resOLS Ypred ResAMMI[,1] YpredAMMI[,1] AMMI0
#>    <fct> <fct> <dbl>   <dbl> <dbl>       <dbl>         <dbl> <dbl>
#>  1 E1    G1     2.37 -0.0843  2.45     0.0712           2.52  2.45
#>  2 E1    G10    1.97 -0.344   2.32    -0.354            1.96  2.32
#>  3 E1    G2     2.90  0.311   2.59     0.290            2.88  2.59
#>  4 E1    G3     2.89  0.0868  2.80    -0.0452           2.76  2.80
#>  5 E1    G4     2.59  0.100   2.49     0.0494           2.54  2.49
#>  6 E1    G5     2.19 -0.196   2.38    -0.0709           2.31  2.38
#>  7 E1    G6     2.30 -0.0797  2.38    -0.0829           2.30  2.38
#>  8 E1    G7     2.77  0.186   2.59     0.164            2.75  2.59
#>  9 E1    G8     2.90  0.0493  2.85    -0.00536          2.84  2.85
#> 10 E1    G9     2.33 -0.0307  2.36    -0.0170           2.34  2.36
#> # ... with 130 more rows
```

### Biplots

ggplot2-based graphics are easily obtained in metan package. For
example, the well-known AMMI2 biplot may be obtained as follows. Please,
note that since `waas()` function allows analyzing multiple variables at
the same time, e.g., `resp = c(v1, v2, ...)`, the output `model` is a
list, in this case with one element, GY.

``` r
a <- plot_scores(model)
b <- plot_scores(model,
                 type = 2,
                 polygon = TRUE,
                 col.env = "gray70",
                 col.segm.env = "gray70",
                 axis.expand = 1.5)
c <- plot_scores(model, type = 4)
arrange_ggplot(a, b, c, labels = letters[1:3], nrow = 1)
```

![](man/figures/README-AMMI-1.png)<!-- -->

## GGE model

The GGE model is fitted with the function `gge()`. For more details, see
the [complete
vignette](https://tiagoolivoto.github.io/metan/articles/vignettes_gge.html).

``` r
model <- gge(data_ge, ENV, GEN, GY)
model2 <- gge(data_ge, ENV, GEN, GY, svp = "genotype")
model3 <- gge(data_ge, ENV, GEN, GY, svp = "symmetrical")
d <- plot(model)
e <- plot(model2, type = 8)
f <- plot(model2,
          type = 2,
          col.gen = "black",
          col.env = "gray70",
          axis.expand = 1.5)
arrange_ggplot(d, e, f, labels = letters[4:6], nrow = 1)
```

![](man/figures/README-GGE-1.png)<!-- -->

## BLUP model

Linear-mixed effect models to predict the response variable in METs are
fitted using the function `waasb()`. Here we will obtain the predicted
means for genotypes in the variables `GY` and `HM`. For more details,
see the [complete
vignette](https://tiagoolivoto.github.io/metan/articles/vignettes_blup.html).

``` r
model2 <- waasb(data_ge,
                env = ENV,
                gen = GEN,
                rep = REP,
                resp = everything(),
                verbose = FALSE)

get_model_data(model2, what = "genpar")
#> # A tibble: 15 x 3
#>    Parameters                 GY      HM
#>    <chr>                   <dbl>   <dbl>
#>  1 GEI variance           0.0567  2.19  
#>  2 GEI (%)               31.3    39.7   
#>  3 Genotypic variance     0.0280  0.490 
#>  4 Gen (%)               15.4     8.88  
#>  5 Residual variance      0.0967  2.84  
#>  6 Res (%)               53.3    51.5   
#>  7 Phenotypic variance    0.181   5.52  
#>  8 Heritability           0.154   0.0888
#>  9 GEIr2                  0.313   0.397 
#> 10 Heribatility of means  0.815   0.686 
#> 11 Accuracy               0.903   0.828 
#> 12 rge                    0.370   0.435 
#> 13 CVg                    6.26    1.46  
#> 14 CVr                   11.6     3.50  
#> 15 CV ratio               0.538   0.415
get_model_data(model2, what = "blupg")
#> # A tibble: 10 x 3
#>    gen      GY    HM
#>    <fct> <dbl> <dbl>
#>  1 G1     2.62  47.4
#>  2 G10    2.51  48.4
#>  3 G2     2.73  47.1
#>  4 G3     2.90  47.8
#>  5 G4     2.65  48.0
#>  6 G5     2.56  48.9
#>  7 G6     2.56  48.5
#>  8 G7     2.73  48.0
#>  9 G8     2.94  48.8
#> 10 G9     2.54  48.0
```

### Plotting the BLUPs for genotypes

To produce a plot with the predicted means, use the function
`plot_blup()`.

``` r
g <- plot_blup(model2)
h <- plot_blup(model2,
               prob = 0.1,
               col.shape  =  c("gray20", "gray80")) + ggplot2::coord_flip()
arrange_ggplot(g, h, labels = letters[7:8])
```

![](man/figures/README-BLUP-1.png)<!-- -->

### BLUPS for genotype-vs-environment interaction

The object `BLUPgge` contains the blups for the genotype-vs-environment
interaction. In the following example, the values for `GY` are shown.

``` r
model2$GY$BLUPgge
#> # A tibble: 140 x 8
#>    ENV   GEN    BLUPge   BLUPg `BLUPg+ge` Predicted    LL    UL
#>    <fct> <fct>   <dbl>   <dbl>      <dbl>     <dbl> <dbl> <dbl>
#>  1 E1    G1    -0.0621 -0.0575    -0.120       2.40  2.30  2.50
#>  2 E1    G10   -0.243  -0.166     -0.409       2.11  2.01  2.22
#>  3 E1    G2     0.207   0.0570     0.264       2.78  2.68  2.89
#>  4 E1    G3     0.0885  0.229      0.318       2.84  2.73  2.94
#>  5 E1    G4     0.0601 -0.0264     0.0337      2.55  2.45  2.66
#>  6 E1    G5    -0.141  -0.112     -0.252       2.27  2.16  2.37
#>  7 E1    G6    -0.0673 -0.114     -0.182       2.34  2.24  2.44
#>  8 E1    G7     0.127   0.0543     0.181       2.70  2.60  2.81
#>  9 E1    G8     0.0702  0.269      0.339       2.86  2.76  2.96
#> 10 E1    G9    -0.0389 -0.134     -0.173       2.35  2.24  2.45
#> # ... with 130 more rows
```

When more than one variable is fitted, the predicted means for
genotype-vs-environment combination may be obtained for all variables in
the model using `get_model_data()`.

``` r
get_model_data(model2, what = "blupge")
#> # A tibble: 140 x 4
#>    ENV   GEN      GY    HM
#>    <fct> <fct> <dbl> <dbl>
#>  1 E1    G1     2.40  46.6
#>  2 E1    G10    2.11  47.2
#>  3 E1    G2     2.78  45.7
#>  4 E1    G3     2.84  46.2
#>  5 E1    G4     2.55  48.0
#>  6 E1    G5     2.27  49.4
#>  7 E1    G6     2.34  48.1
#>  8 E1    G7     2.70  47.4
#>  9 E1    G8     2.86  48.0
#> 10 E1    G9     2.35  47.6
#> # ... with 130 more rows
```

# Computing parametric and non-parametric stability indexes

``` r
stats <- ge_stats(data_ge, ENV, GEN, REP, GY)
get_model_data(stats, "stats")
#> # A tibble: 10 x 33
#>    var   gen       Y    CV   Var Shukla  Wi_g  Wi_f  Wi_u Ecoval   bij      Sij
#>    <chr> <chr> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl>    <dbl>
#>  1 GY    G1     2.60  35.2 10.9  0.0280  84.4  89.2  81.1  1.22  1.06  -0.00142
#>  2 GY    G10    2.47  42.3 14.2  0.244   59.2  64.6  54.4  7.96  1.12   0.177  
#>  3 GY    G2     2.74  34.0 11.3  0.0861  82.8  95.3  75.6  3.03  1.05   0.0497 
#>  4 GY    G3     2.96  29.9 10.1  0.0121 104.   99.7 107.   0.725 1.03  -0.0128 
#>  5 GY    G4     2.64  31.4  8.93 0.0640  85.9  79.5  91.9  2.34  0.937  0.0298 
#>  6 GY    G5     2.54  30.6  7.82 0.0480  82.7  82.2  82.4  1.84  0.887  0.00902
#>  7 GY    G6     2.53  29.7  7.34 0.0468  83.0  83.7  81.8  1.81  0.861  0.00304
#>  8 GY    G7     2.74  27.4  7.33 0.122   83.9  77.6  93.4  4.16  0.819  0.0579 
#>  9 GY    G8     3.00  30.4 10.8  0.0712  98.8  90.5 107.   2.57  1.03   0.0382 
#> 10 GY    G9     2.51  42.4 14.7  0.167   68.8  68.9  70.3  5.56  1.19   0.0938 
#> # ... with 21 more variables: R2 <dbl>, ASV <dbl>, SIPC <dbl>, EV <dbl>,
#> #   ZA <dbl>, WAAS <dbl>, HMGV <dbl>, RPGV <dbl>, HMRPGV <dbl>, Pi_a <dbl>,
#> #   Pi_f <dbl>, Pi_u <dbl>, Gai <dbl>, S1 <dbl>, S2 <dbl>, S3 <dbl>, S6 <dbl>,
#> #   N1 <dbl>, N2 <dbl>, N3 <dbl>, N4 <dbl>
get_model_data(stats, "ranks")
#> # A tibble: 10 x 32
#>    var   gen     Y_R  CV_R Var_R Shukla_R Wi_g_R Wi_f_R Wi_u_R Ecoval_R Sij_R
#>    <chr> <chr> <dbl> <dbl> <dbl>    <dbl>  <dbl>  <dbl>  <dbl>    <dbl> <dbl>
#>  1 GY    G1        6     8     7        2      4      4      7        2     1
#>  2 GY    G10      10     9     9       10     10     10     10       10    10
#>  3 GY    G2        3     7     8        7      7      2      8        7     7
#>  4 GY    G3        2     3     5        1      1      1      2        1     4
#>  5 GY    G4        5     6     4        5      3      7      4        5     5
#>  6 GY    G5        7     5     3        4      8      6      5        4     3
#>  7 GY    G6        8     2     2        3      6      5      6        3     2
#>  8 GY    G7        4     1     1        8      5      8      3        8     8
#>  9 GY    G8        1     4     6        6      2      3      1        6     6
#> 10 GY    G9        9    10    10        9      9      9      9        9     9
#> # ... with 21 more variables: R2_R <dbl>, ASV_R <dbl>, SIPC_R <dbl>,
#> #   EV_R <dbl>, ZA_R <dbl>, WAAS_R <dbl>, HMGV_R <dbl>, RPGV_R <dbl>,
#> #   HMRPGV_R <dbl>, Pi_a_R <dbl>, Pi_f_R <dbl>, Pi_u_R <dbl>, Gai_R <dbl>,
#> #   S1_R <dbl>, S2_R <dbl>, S3_R <dbl>, S6_R <dbl>, N1_R <dbl>, N2_R <dbl>,
#> #   N3_R <dbl>, N4_R <dbl>
```

# Getting help

  - If you encounter a clear bug, please file a minimal reproducible
    example on [github](https://github.com/TiagoOlivoto/metan/issues)

  - Suggestions and criticisms to improve the quality and usability of
    the package are welcome\!
