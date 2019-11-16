
<!-- README.md is generated from README.Rmd. Please edit that file -->

# metan <img src="man/figures/logo.png" align="right" height=140/>

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

The latest development version can be download from GitHub by running

``` r
# install.packages("devtools")
devtools::install_github("TiagoOlivoto/metan")
```

# Getting started

Here, we will use the example dataset `data_ge` that contains data on
two variables assessed in 10 genotypes growing in 14 environments. For
more details see `?data_ge`. The tables of this vignette were produced
with the package
[`KableExtra`](http://haozhu233.github.io/kableExtra/awesome_table_in_html.html)

``` r
library(metan)
library(ggplot2) # used to create the plots
library(kableExtra) # Used to produce HTML tables
print_table = function(table){
  kable(table, "html", digits = 3, escape = TRUE) %>%
    kable_styling(bootstrap_options = c("striped", "hover"), font_size = 12)
}
str(data_ge)
#> Classes 'tbl_df', 'tbl' and 'data.frame':    420 obs. of  5 variables:
#>  $ ENV: Factor w/ 14 levels "E1","E10","E11",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ GEN: Factor w/ 10 levels "G1","G10","G2",..: 1 1 1 3 3 3 4 4 4 5 ...
#>  $ REP: Factor w/ 3 levels "1","2","3": 1 2 3 1 2 3 1 2 3 1 ...
#>  $ GY : num  2.17 2.5 2.43 3.21 2.93 ...
#>  $ HM : num  44.9 46.9 47.8 45.2 45.3 ...
```

## AMMI model

### Fitting the model

The AMMI model is fitted with the function `performs_ammi()`. For more
details, see the [complete
vignette](https://tiagoolivoto.github.io/metan/articles/vignettes_ammi.html).

``` r
model <- performs_ammi(data_ge,
                       resp = GY,
                       gen = GEN,
                       env = ENV,
                       rep = REP,
                       verbose = FALSE)
```

### Predicting the response variable

The S3 method `predict()` is implemented for objects of class
`performs_ammi` and may be used to estimate the response of each
genotype in each environment considering different number of Interaction
Principal Component Axis (IPCA). For example, we will use four IPCA
(number of significant IPCAs) to estimate the variable GY using the
`model` object. Note that `$GY` was used because using the `predict()`
function to a model of class `waas` return a list, with one `tbl_df` for
each variable.

``` r
predict(model, naxis = 4)$GY %>% 
  head() %>% 
  print_table()
```

<table class="table table-striped table-hover" style="font-size: 12px; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

ENV

</th>

<th style="text-align:left;">

GEN

</th>

<th style="text-align:right;">

Y

</th>

<th style="text-align:right;">

resOLS

</th>

<th style="text-align:right;">

Ypred

</th>

<th style="text-align:right;">

ResAMMI

</th>

<th style="text-align:right;">

YpredAMMI

</th>

<th style="text-align:right;">

AMMI0

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

E1

</td>

<td style="text-align:left;">

G1

</td>

<td style="text-align:right;">

2.366

</td>

<td style="text-align:right;">

\-0.084

</td>

<td style="text-align:right;">

2.450

</td>

<td style="text-align:right;">

0.07115484

</td>

<td style="text-align:right;">

2.521273

</td>

<td style="text-align:right;">

2.450

</td>

</tr>

<tr>

<td style="text-align:left;">

E1

</td>

<td style="text-align:left;">

G10

</td>

<td style="text-align:right;">

1.974

</td>

<td style="text-align:right;">

\-0.344

</td>

<td style="text-align:right;">

2.318

</td>

<td style="text-align:right;">

\-0.35391141

</td>

<td style="text-align:right;">

1.963751

</td>

<td style="text-align:right;">

2.318

</td>

</tr>

<tr>

<td style="text-align:left;">

E1

</td>

<td style="text-align:left;">

G2

</td>

<td style="text-align:right;">

2.902

</td>

<td style="text-align:right;">

0.311

</td>

<td style="text-align:right;">

2.591

</td>

<td style="text-align:right;">

0.29035016

</td>

<td style="text-align:right;">

2.880939

</td>

<td style="text-align:right;">

2.591

</td>

</tr>

<tr>

<td style="text-align:left;">

E1

</td>

<td style="text-align:left;">

G3

</td>

<td style="text-align:right;">

2.889

</td>

<td style="text-align:right;">

0.087

</td>

<td style="text-align:right;">

2.802

</td>

<td style="text-align:right;">

\-0.04518795

</td>

<td style="text-align:right;">

2.756598

</td>

<td style="text-align:right;">

2.802

</td>

</tr>

<tr>

<td style="text-align:left;">

E1

</td>

<td style="text-align:left;">

G4

</td>

<td style="text-align:right;">

2.589

</td>

<td style="text-align:right;">

0.100

</td>

<td style="text-align:right;">

2.488

</td>

<td style="text-align:right;">

0.04942370

</td>

<td style="text-align:right;">

2.537781

</td>

<td style="text-align:right;">

2.488

</td>

</tr>

<tr>

<td style="text-align:left;">

E1

</td>

<td style="text-align:left;">

G5

</td>

<td style="text-align:right;">

2.188

</td>

<td style="text-align:right;">

\-0.196

</td>

<td style="text-align:right;">

2.384

</td>

<td style="text-align:right;">

\-0.07091881

</td>

<td style="text-align:right;">

2.312867

</td>

<td style="text-align:right;">

2.384

</td>

</tr>

</tbody>

</table>

### Biplots

ggplot2-based graphics are easily obtained in metan package. For
example, the well-known AMMI2 biplot may be obtained as follows. Please,
note that since `waas()` function allows analyzing multiple variables at
the same time, e.g., `resp = c(v1, v2, ...)`, the output `model` is a
list, in this case with one element, GY.

``` r
a <- plot_scores(model$GY, axis.expand = 1.5)
b <- plot_scores(model$GY,
                 type = 2,
                 polygon = TRUE,
                 col.env = "gray70",
                 col.segm.env = "gray70",
                 axis.expand = 1.5)
arrange_ggplot(a, b, labels = letters[1:2])
```

<img src="man/figures/AMMI.png"/>

## GGE model

The GGE model is fitted with the function `gge()`. For more details, see
the [complete
vignette](https://tiagoolivoto.github.io/metan/articles/vignettes_gge.html).

``` r
model <- gge(data_ge, ENV, GEN, GY)
model2 <- gge(data_ge, ENV, GEN, GY, svp = "symmetrical")
c <- plot(model)
d <- plot(model2,
          type = 2,
          col.gen = "black",
          col.env = "gray70",
          axis.expand = 1.5)
arrange_ggplot(c, d, labels = letters[3:4])
```

<img src="man/figures/GGE.png"/>

## BLUP model

Linear-mixed effect models to predict the response variable in METs are
fitted using the function `waasb()`. Here we will obtain the predicted
means for genotypes in the variables `GY` and `HM`. For more details,
see the [complete
vignette](https://tiagoolivoto.github.io/metan/articles/vignettes_blup.html).

``` r
model2 <- waasb(data_ge, ENV, GEN, REP,
                resp = c(GY, HM),
                verbose = FALSE)

get_model_data(model2, what = "genpar") %>% 
print_table()
```

<table class="table table-striped table-hover" style="font-size: 12px; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Parameters

</th>

<th style="text-align:right;">

GY

</th>

<th style="text-align:right;">

HM

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

GEI variance

</td>

<td style="text-align:right;">

0.057

</td>

<td style="text-align:right;">

2.189

</td>

</tr>

<tr>

<td style="text-align:left;">

GEI (%)

</td>

<td style="text-align:right;">

31.259

</td>

<td style="text-align:right;">

39.665

</td>

</tr>

<tr>

<td style="text-align:left;">

Genotypic variance

</td>

<td style="text-align:right;">

0.028

</td>

<td style="text-align:right;">

0.490

</td>

</tr>

<tr>

<td style="text-align:left;">

Gen (%)

</td>

<td style="text-align:right;">

15.447

</td>

<td style="text-align:right;">

8.876

</td>

</tr>

<tr>

<td style="text-align:left;">

Residual variance

</td>

<td style="text-align:right;">

0.097

</td>

<td style="text-align:right;">

2.840

</td>

</tr>

<tr>

<td style="text-align:left;">

Res (%)

</td>

<td style="text-align:right;">

53.294

</td>

<td style="text-align:right;">

51.458

</td>

</tr>

<tr>

<td style="text-align:left;">

Phenotypic variance

</td>

<td style="text-align:right;">

0.181

</td>

<td style="text-align:right;">

5.519

</td>

</tr>

<tr>

<td style="text-align:left;">

Heritability

</td>

<td style="text-align:right;">

0.154

</td>

<td style="text-align:right;">

0.089

</td>

</tr>

<tr>

<td style="text-align:left;">

GEIr2

</td>

<td style="text-align:right;">

0.313

</td>

<td style="text-align:right;">

0.397

</td>

</tr>

<tr>

<td style="text-align:left;">

Heribatility of means

</td>

<td style="text-align:right;">

0.815

</td>

<td style="text-align:right;">

0.686

</td>

</tr>

<tr>

<td style="text-align:left;">

Accuracy

</td>

<td style="text-align:right;">

0.903

</td>

<td style="text-align:right;">

0.828

</td>

</tr>

<tr>

<td style="text-align:left;">

rge

</td>

<td style="text-align:right;">

0.370

</td>

<td style="text-align:right;">

0.435

</td>

</tr>

<tr>

<td style="text-align:left;">

CVg

</td>

<td style="text-align:right;">

6.260

</td>

<td style="text-align:right;">

1.456

</td>

</tr>

<tr>

<td style="text-align:left;">

CVr

</td>

<td style="text-align:right;">

11.628

</td>

<td style="text-align:right;">

3.504

</td>

</tr>

<tr>

<td style="text-align:left;">

CV ratio

</td>

<td style="text-align:right;">

0.538

</td>

<td style="text-align:right;">

0.415

</td>

</tr>

</tbody>

</table>

``` r

get_model_data(model2, what = "blupg") %>% 
print_table()
```

<table class="table table-striped table-hover" style="font-size: 12px; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

gen

</th>

<th style="text-align:right;">

GY

</th>

<th style="text-align:right;">

HM

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

G1

</td>

<td style="text-align:right;">

2.617

</td>

<td style="text-align:right;">

47.396

</td>

</tr>

<tr>

<td style="text-align:left;">

G10

</td>

<td style="text-align:right;">

2.509

</td>

<td style="text-align:right;">

48.375

</td>

</tr>

<tr>

<td style="text-align:left;">

G2

</td>

<td style="text-align:right;">

2.731

</td>

<td style="text-align:right;">

47.108

</td>

</tr>

<tr>

<td style="text-align:left;">

G3

</td>

<td style="text-align:right;">

2.903

</td>

<td style="text-align:right;">

47.756

</td>

</tr>

<tr>

<td style="text-align:left;">

G4

</td>

<td style="text-align:right;">

2.648

</td>

<td style="text-align:right;">

48.050

</td>

</tr>

<tr>

<td style="text-align:left;">

G5

</td>

<td style="text-align:right;">

2.563

</td>

<td style="text-align:right;">

48.918

</td>

</tr>

<tr>

<td style="text-align:left;">

G6

</td>

<td style="text-align:right;">

2.560

</td>

<td style="text-align:right;">

48.530

</td>

</tr>

<tr>

<td style="text-align:left;">

G7

</td>

<td style="text-align:right;">

2.729

</td>

<td style="text-align:right;">

48.004

</td>

</tr>

<tr>

<td style="text-align:left;">

G8

</td>

<td style="text-align:right;">

2.943

</td>

<td style="text-align:right;">

48.784

</td>

</tr>

<tr>

<td style="text-align:left;">

G9

</td>

<td style="text-align:right;">

2.541

</td>

<td style="text-align:right;">

47.962

</td>

</tr>

</tbody>

</table>

### Plotting the BLUPs for genotypes

To produce a plot with the predicted means, use the function
`plot_blup()`.

``` r
e <- plot_blup(model2$GY)
f <- plot_blup(model2$GY,
               prob = 0.1,
               col.shape  =  c("gray20", "gray80")) + coord_flip()
arrange_ggplot(e, f, labels = letters[5:6])
```

<img src="man/figures/BLUP.png" height=340/>

### BLUPS for genotype-vs-environment interaction

The object `BLUPgge` contains the blups for the genotype-vs-environment
interaction. In the following example, the values for `GY` are shown.

``` r
model2$GY$BLUPgge[1:5,] %>% 
  head() %>% 
  print_table()
```

<table class="table table-striped table-hover" style="font-size: 12px; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

ENV

</th>

<th style="text-align:left;">

GEN

</th>

<th style="text-align:right;">

BLUPge

</th>

<th style="text-align:right;">

BLUPg

</th>

<th style="text-align:right;">

BLUPg+ge

</th>

<th style="text-align:right;">

Predicted

</th>

<th style="text-align:right;">

LL

</th>

<th style="text-align:right;">

UL

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

E1

</td>

<td style="text-align:left;">

G1

</td>

<td style="text-align:right;">

\-0.062

</td>

<td style="text-align:right;">

\-0.058

</td>

<td style="text-align:right;">

\-0.120

</td>

<td style="text-align:right;">

2.401

</td>

<td style="text-align:right;">

2.298

</td>

<td style="text-align:right;">

2.505

</td>

</tr>

<tr>

<td style="text-align:left;">

E1

</td>

<td style="text-align:left;">

G10

</td>

<td style="text-align:right;">

\-0.243

</td>

<td style="text-align:right;">

\-0.166

</td>

<td style="text-align:right;">

\-0.409

</td>

<td style="text-align:right;">

2.112

</td>

<td style="text-align:right;">

2.009

</td>

<td style="text-align:right;">

2.216

</td>

</tr>

<tr>

<td style="text-align:left;">

E1

</td>

<td style="text-align:left;">

G2

</td>

<td style="text-align:right;">

0.207

</td>

<td style="text-align:right;">

0.057

</td>

<td style="text-align:right;">

0.264

</td>

<td style="text-align:right;">

2.784

</td>

<td style="text-align:right;">

2.681

</td>

<td style="text-align:right;">

2.888

</td>

</tr>

<tr>

<td style="text-align:left;">

E1

</td>

<td style="text-align:left;">

G3

</td>

<td style="text-align:right;">

0.088

</td>

<td style="text-align:right;">

0.229

</td>

<td style="text-align:right;">

0.318

</td>

<td style="text-align:right;">

2.838

</td>

<td style="text-align:right;">

2.735

</td>

<td style="text-align:right;">

2.942

</td>

</tr>

<tr>

<td style="text-align:left;">

E1

</td>

<td style="text-align:left;">

G4

</td>

<td style="text-align:right;">

0.060

</td>

<td style="text-align:right;">

\-0.026

</td>

<td style="text-align:right;">

0.034

</td>

<td style="text-align:right;">

2.554

</td>

<td style="text-align:right;">

2.451

</td>

<td style="text-align:right;">

2.658

</td>

</tr>

</tbody>

</table>

When more than one variable is fitted, the predicted means for
genotype-vs-environment combination may be obtained for all variables in
the model using `get_model_data()`.

``` r
get_model_data(model2, what = "blupge") %>% 
  head() %>% 
  print_table()
```

<table class="table table-striped table-hover" style="font-size: 12px; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

ENV

</th>

<th style="text-align:left;">

GEN

</th>

<th style="text-align:right;">

GY

</th>

<th style="text-align:right;">

HM

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

E1

</td>

<td style="text-align:left;">

G1

</td>

<td style="text-align:right;">

2.401

</td>

<td style="text-align:right;">

46.587

</td>

</tr>

<tr>

<td style="text-align:left;">

E1

</td>

<td style="text-align:left;">

G10

</td>

<td style="text-align:right;">

2.112

</td>

<td style="text-align:right;">

47.152

</td>

</tr>

<tr>

<td style="text-align:left;">

E1

</td>

<td style="text-align:left;">

G2

</td>

<td style="text-align:right;">

2.784

</td>

<td style="text-align:right;">

45.664

</td>

</tr>

<tr>

<td style="text-align:left;">

E1

</td>

<td style="text-align:left;">

G3

</td>

<td style="text-align:right;">

2.838

</td>

<td style="text-align:right;">

46.249

</td>

</tr>

<tr>

<td style="text-align:left;">

E1

</td>

<td style="text-align:left;">

G4

</td>

<td style="text-align:right;">

2.554

</td>

<td style="text-align:right;">

48.017

</td>

</tr>

<tr>

<td style="text-align:left;">

E1

</td>

<td style="text-align:left;">

G5

</td>

<td style="text-align:right;">

2.268

</td>

<td style="text-align:right;">

49.390

</td>

</tr>

</tbody>

</table>

# Computing parametric and non-parametric stability indexes

``` r
stats <- ge_stats(data_ge, ENV, GEN, REP, GY)
get_model_data(stats, "stats") %>% 
  print_table()
```

<table class="table table-striped table-hover" style="font-size: 12px; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

var

</th>

<th style="text-align:left;">

gen

</th>

<th style="text-align:right;">

Y

</th>

<th style="text-align:right;">

Var

</th>

<th style="text-align:right;">

Shukla

</th>

<th style="text-align:right;">

Wi\_g

</th>

<th style="text-align:right;">

Wi\_f

</th>

<th style="text-align:right;">

Wi\_u

</th>

<th style="text-align:right;">

Ecoval

</th>

<th style="text-align:right;">

bij

</th>

<th style="text-align:right;">

Sij

</th>

<th style="text-align:right;">

R2

</th>

<th style="text-align:right;">

ASV

</th>

<th style="text-align:right;">

SIPC

</th>

<th style="text-align:right;">

EV

</th>

<th style="text-align:right;">

ZA

</th>

<th style="text-align:right;">

WAAS

</th>

<th style="text-align:right;">

HMGV

</th>

<th style="text-align:right;">

RPGV

</th>

<th style="text-align:right;">

HMRPGV

</th>

<th style="text-align:right;">

Pi\_a

</th>

<th style="text-align:right;">

Pi\_f

</th>

<th style="text-align:right;">

Pi\_u

</th>

<th style="text-align:right;">

Gai

</th>

<th style="text-align:right;">

S1

</th>

<th style="text-align:right;">

S2

</th>

<th style="text-align:right;">

S3

</th>

<th style="text-align:right;">

S6

</th>

<th style="text-align:right;">

N1

</th>

<th style="text-align:right;">

N2

</th>

<th style="text-align:right;">

N3

</th>

<th style="text-align:right;">

N4

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G1

</td>

<td style="text-align:right;">

2.604

</td>

<td style="text-align:right;">

10.897

</td>

<td style="text-align:right;">

0.028

</td>

<td style="text-align:right;">

84.350

</td>

<td style="text-align:right;">

89.172

</td>

<td style="text-align:right;">

81.062

</td>

<td style="text-align:right;">

1.219

</td>

<td style="text-align:right;">

1.063

</td>

<td style="text-align:right;">

\-0.001

</td>

<td style="text-align:right;">

0.966

</td>

<td style="text-align:right;">

0.346

</td>

<td style="text-align:right;">

0.463

</td>

<td style="text-align:right;">

0.015

</td>

<td style="text-align:right;">

0.100

</td>

<td style="text-align:right;">

0.151

</td>

<td style="text-align:right;">

2.316

</td>

<td style="text-align:right;">

0.969

</td>

<td style="text-align:right;">

0.967

</td>

<td style="text-align:right;">

0.169

</td>

<td style="text-align:right;">

0.228

</td>

<td style="text-align:right;">

0.125

</td>

<td style="text-align:right;">

2.449

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:right;">

8.571

</td>

<td style="text-align:right;">

16.581

</td>

<td style="text-align:right;">

6.064

</td>

<td style="text-align:right;">

2.571

</td>

<td style="text-align:right;">

0.367

</td>

<td style="text-align:right;">

0.429

</td>

<td style="text-align:right;">

0.000

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G10

</td>

<td style="text-align:right;">

2.471

</td>

<td style="text-align:right;">

14.237

</td>

<td style="text-align:right;">

0.244

</td>

<td style="text-align:right;">

59.241

</td>

<td style="text-align:right;">

64.605

</td>

<td style="text-align:right;">

54.394

</td>

<td style="text-align:right;">

7.955

</td>

<td style="text-align:right;">

1.122

</td>

<td style="text-align:right;">

0.177

</td>

<td style="text-align:right;">

0.823

</td>

<td style="text-align:right;">

1.225

</td>

<td style="text-align:right;">

2.069

</td>

<td style="text-align:right;">

0.210

</td>

<td style="text-align:right;">

0.437

</td>

<td style="text-align:right;">

0.652

</td>

<td style="text-align:right;">

2.112

</td>

<td style="text-align:right;">

0.913

</td>

<td style="text-align:right;">

0.896

</td>

<td style="text-align:right;">

0.344

</td>

<td style="text-align:right;">

0.475

</td>

<td style="text-align:right;">

0.245

</td>

<td style="text-align:right;">

2.247

</td>

<td style="text-align:right;">

0.022

</td>

<td style="text-align:right;">

16.879

</td>

<td style="text-align:right;">

32.276

</td>

<td style="text-align:right;">

9.655

</td>

<td style="text-align:right;">

3.714

</td>

<td style="text-align:right;">

0.531

</td>

<td style="text-align:right;">

0.577

</td>

<td style="text-align:right;">

0.003

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G2

</td>

<td style="text-align:right;">

2.744

</td>

<td style="text-align:right;">

11.345

</td>

<td style="text-align:right;">

0.086

</td>

<td style="text-align:right;">

82.797

</td>

<td style="text-align:right;">

95.294

</td>

<td style="text-align:right;">

75.616

</td>

<td style="text-align:right;">

3.032

</td>

<td style="text-align:right;">

1.054

</td>

<td style="text-align:right;">

0.050

</td>

<td style="text-align:right;">

0.913

</td>

<td style="text-align:right;">

0.249

</td>

<td style="text-align:right;">

1.544

</td>

<td style="text-align:right;">

0.179

</td>

<td style="text-align:right;">

0.216

</td>

<td style="text-align:right;">

0.283

</td>

<td style="text-align:right;">

2.466

</td>

<td style="text-align:right;">

1.025

</td>

<td style="text-align:right;">

1.020

</td>

<td style="text-align:right;">

0.126

</td>

<td style="text-align:right;">

0.149

</td>

<td style="text-align:right;">

0.108

</td>

<td style="text-align:right;">

2.594

</td>

<td style="text-align:right;">

0.066

</td>

<td style="text-align:right;">

9.604

</td>

<td style="text-align:right;">

17.000

</td>

<td style="text-align:right;">

4.773

</td>

<td style="text-align:right;">

2.429

</td>

<td style="text-align:right;">

0.540

</td>

<td style="text-align:right;">

0.634

</td>

<td style="text-align:right;">

0.014

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G3

</td>

<td style="text-align:right;">

2.955

</td>

<td style="text-align:right;">

10.131

</td>

<td style="text-align:right;">

0.012

</td>

<td style="text-align:right;">

103.522

</td>

<td style="text-align:right;">

99.678

</td>

<td style="text-align:right;">

106.863

</td>

<td style="text-align:right;">

0.725

</td>

<td style="text-align:right;">

1.031

</td>

<td style="text-align:right;">

\-0.013

</td>

<td style="text-align:right;">

0.977

</td>

<td style="text-align:right;">

0.113

</td>

<td style="text-align:right;">

0.552

</td>

<td style="text-align:right;">

0.021

</td>

<td style="text-align:right;">

0.080

</td>

<td style="text-align:right;">

0.106

</td>

<td style="text-align:right;">

2.681

</td>

<td style="text-align:right;">

1.105

</td>

<td style="text-align:right;">

1.104

</td>

<td style="text-align:right;">

0.041

</td>

<td style="text-align:right;">

0.072

</td>

<td style="text-align:right;">

0.018

</td>

<td style="text-align:right;">

2.824

</td>

<td style="text-align:right;">

0.022

</td>

<td style="text-align:right;">

5.297

</td>

<td style="text-align:right;">

1.830

</td>

<td style="text-align:right;">

1.525

</td>

<td style="text-align:right;">

2.000

</td>

<td style="text-align:right;">

0.667

</td>

<td style="text-align:right;">

0.862

</td>

<td style="text-align:right;">

0.008

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G4

</td>

<td style="text-align:right;">

2.642

</td>

<td style="text-align:right;">

8.934

</td>

<td style="text-align:right;">

0.064

</td>

<td style="text-align:right;">

85.924

</td>

<td style="text-align:right;">

79.549

</td>

<td style="text-align:right;">

91.944

</td>

<td style="text-align:right;">

2.343

</td>

<td style="text-align:right;">

0.937

</td>

<td style="text-align:right;">

0.030

</td>

<td style="text-align:right;">

0.917

</td>

<td style="text-align:right;">

0.594

</td>

<td style="text-align:right;">

1.036

</td>

<td style="text-align:right;">

0.052

</td>

<td style="text-align:right;">

0.219

</td>

<td style="text-align:right;">

0.326

</td>

<td style="text-align:right;">

2.388

</td>

<td style="text-align:right;">

0.990

</td>

<td style="text-align:right;">

0.988

</td>

<td style="text-align:right;">

0.173

</td>

<td style="text-align:right;">

0.289

</td>

<td style="text-align:right;">

0.085

</td>

<td style="text-align:right;">

2.515

</td>

<td style="text-align:right;">

0.088

</td>

<td style="text-align:right;">

6.269

</td>

<td style="text-align:right;">

11.200

</td>

<td style="text-align:right;">

4.400

</td>

<td style="text-align:right;">

1.929

</td>

<td style="text-align:right;">

0.321

</td>

<td style="text-align:right;">

0.402

</td>

<td style="text-align:right;">

0.015

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G5

</td>

<td style="text-align:right;">

2.537

</td>

<td style="text-align:right;">

7.822

</td>

<td style="text-align:right;">

0.048

</td>

<td style="text-align:right;">

82.717

</td>

<td style="text-align:right;">

82.217

</td>

<td style="text-align:right;">

82.350

</td>

<td style="text-align:right;">

1.844

</td>

<td style="text-align:right;">

0.887

</td>

<td style="text-align:right;">

0.009

</td>

<td style="text-align:right;">

0.937

</td>

<td style="text-align:right;">

0.430

</td>

<td style="text-align:right;">

0.997

</td>

<td style="text-align:right;">

0.043

</td>

<td style="text-align:right;">

0.186

</td>

<td style="text-align:right;">

0.270

</td>

<td style="text-align:right;">

2.302

</td>

<td style="text-align:right;">

0.954

</td>

<td style="text-align:right;">

0.952

</td>

<td style="text-align:right;">

0.240

</td>

<td style="text-align:right;">

0.382

</td>

<td style="text-align:right;">

0.133

</td>

<td style="text-align:right;">

2.422

</td>

<td style="text-align:right;">

0.000

</td>

<td style="text-align:right;">

8.571

</td>

<td style="text-align:right;">

16.804

</td>

<td style="text-align:right;">

6.471

</td>

<td style="text-align:right;">

2.714

</td>

<td style="text-align:right;">

0.339

</td>

<td style="text-align:right;">

0.384

</td>

<td style="text-align:right;">

0.000

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G6

</td>

<td style="text-align:right;">

2.534

</td>

<td style="text-align:right;">

7.339

</td>

<td style="text-align:right;">

0.047

</td>

<td style="text-align:right;">

82.995

</td>

<td style="text-align:right;">

83.668

</td>

<td style="text-align:right;">

81.753

</td>

<td style="text-align:right;">

1.806

</td>

<td style="text-align:right;">

0.861

</td>

<td style="text-align:right;">

0.003

</td>

<td style="text-align:right;">

0.942

</td>

<td style="text-align:right;">

0.265

</td>

<td style="text-align:right;">

1.140

</td>

<td style="text-align:right;">

0.091

</td>

<td style="text-align:right;">

0.172

</td>

<td style="text-align:right;">

0.233

</td>

<td style="text-align:right;">

2.304

</td>

<td style="text-align:right;">

0.954

</td>

<td style="text-align:right;">

0.952

</td>

<td style="text-align:right;">

0.238

</td>

<td style="text-align:right;">

0.377

</td>

<td style="text-align:right;">

0.134

</td>

<td style="text-align:right;">

2.423

</td>

<td style="text-align:right;">

0.022

</td>

<td style="text-align:right;">

8.725

</td>

<td style="text-align:right;">

13.982

</td>

<td style="text-align:right;">

5.228

</td>

<td style="text-align:right;">

2.429

</td>

<td style="text-align:right;">

0.347

</td>

<td style="text-align:right;">

0.411

</td>

<td style="text-align:right;">

0.003

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G7

</td>

<td style="text-align:right;">

2.741

</td>

<td style="text-align:right;">

7.331

</td>

<td style="text-align:right;">

0.122

</td>

<td style="text-align:right;">

83.876

</td>

<td style="text-align:right;">

77.552

</td>

<td style="text-align:right;">

93.387

</td>

<td style="text-align:right;">

4.162

</td>

<td style="text-align:right;">

0.819

</td>

<td style="text-align:right;">

0.058

</td>

<td style="text-align:right;">

0.852

</td>

<td style="text-align:right;">

0.663

</td>

<td style="text-align:right;">

1.787

</td>

<td style="text-align:right;">

0.191

</td>

<td style="text-align:right;">

0.303

</td>

<td style="text-align:right;">

0.428

</td>

<td style="text-align:right;">

2.520

</td>

<td style="text-align:right;">

1.037

</td>

<td style="text-align:right;">

1.031

</td>

<td style="text-align:right;">

0.149

</td>

<td style="text-align:right;">

0.318

</td>

<td style="text-align:right;">

0.021

</td>

<td style="text-align:right;">

2.641

</td>

<td style="text-align:right;">

0.011

</td>

<td style="text-align:right;">

8.308

</td>

<td style="text-align:right;">

15.306

</td>

<td style="text-align:right;">

4.777

</td>

<td style="text-align:right;">

2.286

</td>

<td style="text-align:right;">

0.508

</td>

<td style="text-align:right;">

0.564

</td>

<td style="text-align:right;">

0.002

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G8

</td>

<td style="text-align:right;">

3.004

</td>

<td style="text-align:right;">

10.817

</td>

<td style="text-align:right;">

0.071

</td>

<td style="text-align:right;">

98.799

</td>

<td style="text-align:right;">

90.521

</td>

<td style="text-align:right;">

107.307

</td>

<td style="text-align:right;">

2.569

</td>

<td style="text-align:right;">

1.034

</td>

<td style="text-align:right;">

0.038

</td>

<td style="text-align:right;">

0.922

</td>

<td style="text-align:right;">

0.574

</td>

<td style="text-align:right;">

1.178

</td>

<td style="text-align:right;">

0.067

</td>

<td style="text-align:right;">

0.225

</td>

<td style="text-align:right;">

0.327

</td>

<td style="text-align:right;">

2.740

</td>

<td style="text-align:right;">

1.126

</td>

<td style="text-align:right;">

1.122

</td>

<td style="text-align:right;">

0.041

</td>

<td style="text-align:right;">

0.088

</td>

<td style="text-align:right;">

0.006

</td>

<td style="text-align:right;">

2.877

</td>

<td style="text-align:right;">

0.033

</td>

<td style="text-align:right;">

6.571

</td>

<td style="text-align:right;">

5.732

</td>

<td style="text-align:right;">

1.951

</td>

<td style="text-align:right;">

2.000

</td>

<td style="text-align:right;">

1.000

</td>

<td style="text-align:right;">

1.116

</td>

<td style="text-align:right;">

0.015

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G9

</td>

<td style="text-align:right;">

2.510

</td>

<td style="text-align:right;">

14.745

</td>

<td style="text-align:right;">

0.167

</td>

<td style="text-align:right;">

68.849

</td>

<td style="text-align:right;">

68.891

</td>

<td style="text-align:right;">

70.308

</td>

<td style="text-align:right;">

5.563

</td>

<td style="text-align:right;">

1.192

</td>

<td style="text-align:right;">

0.094

</td>

<td style="text-align:right;">

0.897

</td>

<td style="text-align:right;">

0.983

</td>

<td style="text-align:right;">

1.495

</td>

<td style="text-align:right;">

0.131

</td>

<td style="text-align:right;">

0.336

</td>

<td style="text-align:right;">

0.507

</td>

<td style="text-align:right;">

2.177

</td>

<td style="text-align:right;">

0.926

</td>

<td style="text-align:right;">

0.917

</td>

<td style="text-align:right;">

0.291

</td>

<td style="text-align:right;">

0.390

</td>

<td style="text-align:right;">

0.217

</td>

<td style="text-align:right;">

2.303

</td>

<td style="text-align:right;">

0.066

</td>

<td style="text-align:right;">

9.632

</td>

<td style="text-align:right;">

18.276

</td>

<td style="text-align:right;">

6.414

</td>

<td style="text-align:right;">

2.500

</td>

<td style="text-align:right;">

0.357

</td>

<td style="text-align:right;">

0.436

</td>

<td style="text-align:right;">

0.010

</td>

</tr>

</tbody>

</table>

``` r
get_model_data(stats, "ranks") %>% 
  print_table()
```

<table class="table table-striped table-hover" style="font-size: 12px; margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

var

</th>

<th style="text-align:left;">

gen

</th>

<th style="text-align:right;">

Y\_R

</th>

<th style="text-align:right;">

Var\_R

</th>

<th style="text-align:right;">

Shukla\_R

</th>

<th style="text-align:right;">

Wi\_g\_R

</th>

<th style="text-align:right;">

Wi\_f\_R

</th>

<th style="text-align:right;">

Wi\_u\_R

</th>

<th style="text-align:right;">

Ecoval\_R

</th>

<th style="text-align:right;">

Sij\_R

</th>

<th style="text-align:right;">

R2\_R

</th>

<th style="text-align:right;">

ASV\_R

</th>

<th style="text-align:right;">

SIPC\_R

</th>

<th style="text-align:right;">

EV\_R

</th>

<th style="text-align:right;">

ZA\_R

</th>

<th style="text-align:right;">

WAAS\_R

</th>

<th style="text-align:right;">

HMGV\_R

</th>

<th style="text-align:right;">

RPGV\_R

</th>

<th style="text-align:right;">

HMRPGV\_R

</th>

<th style="text-align:right;">

Pi\_a\_R

</th>

<th style="text-align:right;">

Pi\_f\_R

</th>

<th style="text-align:right;">

Pi\_u\_R

</th>

<th style="text-align:right;">

Gai\_R

</th>

<th style="text-align:right;">

S1\_R

</th>

<th style="text-align:right;">

S2\_R

</th>

<th style="text-align:right;">

S3\_R

</th>

<th style="text-align:right;">

S6\_R

</th>

<th style="text-align:right;">

N1\_R

</th>

<th style="text-align:right;">

N2\_R

</th>

<th style="text-align:right;">

N3\_R

</th>

<th style="text-align:right;">

N4\_R

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G1

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

1.5

</td>

<td style="text-align:right;">

5.5

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

8.0

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

1.5

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

5.0

</td>

<td style="text-align:right;">

10.0

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10.0

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

4.5

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G2

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

8.5

</td>

<td style="text-align:right;">

8.0

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

5.5

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

8.0

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G3

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

5.0

</td>

<td style="text-align:right;">

1.0

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

2.5

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

6.0

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G4

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

10.0

</td>

<td style="text-align:right;">

2.0

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

1.0

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

9.0

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G5

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

1.5

</td>

<td style="text-align:right;">

5.5

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9.0

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1.5

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G6

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

5.0

</td>

<td style="text-align:right;">

7.0

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

5.5

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

4.5

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G7

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

3.0

</td>

<td style="text-align:right;">

4.0

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

4.0

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

3.0

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G8

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

7.0

</td>

<td style="text-align:right;">

3.0

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2.5

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

10.0

</td>

</tr>

<tr>

<td style="text-align:left;">

GY

</td>

<td style="text-align:left;">

G9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

8.5

</td>

<td style="text-align:right;">

9.0

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

7.0

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

7.0

</td>

</tr>

</tbody>

</table>

# Getting help

  - If you encounter a clear bug, please file a minimal reproducible
    example on [github](https://github.com/TiagoOlivoto/metan/issues)

  - Suggestions and criticisms to improve the quality and usability of
    the package are welcome\!
