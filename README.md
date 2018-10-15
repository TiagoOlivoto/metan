Extending WAASB package
================
Tiago Olivoto
WAASB 1.0.1 2018-10-15

<style type = "text/css">

body{ /* Normal  */
  color: black;
  }
body .main-container {
        max-width: 100%;
    }
body {
text-align: justify}

</style>
Introducing the `WAASB` R package
=================================

The WAASB R package was developed in R language and is distributed under the GPL 3.0 licence. Our main motivation for the package development was to provide friendly, easy and reliable functions provide de singular value decomposition of BLUP-interaction effects matrix as weel as to share codes used for traditional AMMI analysis and BLUP prediction. The package can be downloaded by clicking [here](https://github.com/TiagoOlivoto/WAASB/archive/master.zip) or running the following code in the R console:

Overview of the `WAASB` package.
--------------------------------

> Package `WAASB` available at `github`. Installation

``` r
# download the package from Github
devtools::install_github("TiagoOlivoto/WAASB")
```

Main functions of `WAASB` package.
----------------------------------

<table>
<colgroup>
<col width="20%" />
<col width="80%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Functions</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><code>performs_ammi()</code></td>
<td align="left">Performs Additive Main effects and Multiplicative Interaction</td>
</tr>
<tr class="even">
<td align="left"><code>plot.eigen()</code></td>
<td align="left">Plot the eigenvalues of singular value decomposition</td>
</tr>
<tr class="odd">
<td align="left"><code>plot.blup()</code></td>
<td align="left">Plot the estimated BLUPs of genotypes</td>
</tr>
<tr class="even">
<td align="left"><code>plot.validation.AMMIF()</code></td>
<td align="left">Plot the RMSE of all AMMI-family tested models</td>
</tr>
<tr class="odd">
<td align="left"><code>plot.WAASBY()</code></td>
<td align="left">Plot WAASBY values for genotype ranking</td>
</tr>
<tr class="even">
<td align="left"><code>plot.WAASBYratio()</code></td>
<td align="left">Plot heat maps with genotype ranking</td>
</tr>
<tr class="odd">
<td align="left"><code>predict.WAAS.AMMI()</code></td>
<td align="left">Predict the means of a <code>WAAS.AMMI</code> object</td>
</tr>
<tr class="even">
<td align="left"><code>summary.WAASB()</code></td>
<td align="left">Summary a WAASB object</td>
</tr>
<tr class="odd">
<td align="left"><code>theme_waasb()</code></td>
<td align="left">A personalized theme for the WAASB ggplot2-based graphics</td>
</tr>
<tr class="even">
<td align="left"><code>validation.blup()</code></td>
<td align="left">Cross-validation for BLUP prediction</td>
</tr>
<tr class="odd">
<td align="left"><code>validation.AMMIF()</code></td>
<td align="left">Cross-validation for the AMMI-family models</td>
</tr>
<tr class="even">
<td align="left"><code>WAASBYratio()</code></td>
<td align="left">BLUP-interaction effects matrix in different scenarios of WAASB/GY ratio</td>
</tr>
<tr class="odd">
<td align="left"><code>WAAS.AMMI()</code></td>
<td align="left">Weighted Average of Absolute Scores for AMMI analysis</td>
</tr>
<tr class="even">
<td align="left"><code>WAASratio.AMMI()</code></td>
<td align="left">Weighted Average of Absolute Scores for AMMI analysis in different scenarios of WAAS/GY ratio</td>
</tr>
<tr class="odd">
<td align="left"><code>WAASB()</code></td>
<td align="left">Weighted Average of Absolute Scores for BLUP-based PCA analysis</td>
</tr>
</tbody>
</table>

Getting started
===============

A data set with data on 10 oat genotypes conducted in 16 environments is provided to make reproducible examples. This data is available in Olivoto ([2018](#ref-Olivoto:2018)). An own data set can be used provided that the columns are in the following order: environment, genotype, block/replicate and response variable(s).

``` r
# Importing data
require(WAASB)
dataset = read.csv("https://data.mendeley.com/datasets/2sjz32k3s3/1/files/07764a07-172a-4285-85db-c31bc39ae480/WAASBdata.csv?dl=1")
```

Predictive accuracy
===================

The functions `validation.AMMIF` provides a complete cross-validation for AMMI model family using replicate-based data. Automatically the first validation is carried out considering the AMMIF model (all possible axis is used). Considering this model, the original data set is split into two data sets: modelling and validating data. The "modelling" data set has all combinations (genotype x environment) with the number of replications informed in the argument `nrepval`. The data set "validating" has one replication. The splitting of the data set into modelling and validating data depends on the design informed. For Completely Randomized Block Design (default), completely blocks are selected within environments, as suggested by Piepho ([1994](#ref-Piepho:1994)). The remained block serves as validation data. If `design = "RCD"` is informed, thus declaring that the a completely randomized design was used, single observations are randomized for each treatment (genotype-by-environment combination). This is the same procedure suggested by Gauch and Zobel ([1988](#ref-Gauch:1988)). The estimated values (for each AMMI model family) are compared with the "validating" data. the Root Means Square Prediction Difference is computed according to the following equation:

$$
   RMSPD = 
   \\sqrt{\\frac{\\sum\_{i = 1}^n(\\widehat{y}\_{ij}-y\_{ij})^2} {n}}
$$

where $\\widehat{y}\_{ij}$ is the model predicted value; and *y*<sub>*i**j*</sub> is the observed value of the validating data set. At the end of boots, a list with all estimated RMSE and the average of RMSE is returned.

The function `validation.blup` provides a cross-validation of replicate-based data using mixed models. By default, complete blocks are randomly selected each environment. The procedure for computing the RSME is identical to the above function.

> The following codes computes the cross-validation of oat data set based on 1000 re-sampling procedures. This number can be changed.

``` r
# cross-validation for AMMI model family
AMMIweat = validation.AMMIF(dataset,
                            resp = GY,
                            gen = GEN,
                            env = ENV,
                            rep = BLOCK,
                            nboot = 100,
                            nrepval = 2)

# cross-validation for BLUP model
BLUPweat = validation.blup(dataset,
                            resp = GY,
                            gen = GEN,
                            env = ENV,
                            rep = BLOCK,
                            nboot = 100,
                            nrepval = 2)
```

A progress bar is shown by default (the examples are bellow). Thus, it is possible to verify the status of the cross-validation process. If necessary, the progress bar can be disabled by informing the argument `progbar = FALSE` in the function.

![Progress bar for AMMI model family cross-validation process.](printAMMI.png)

![Progress bar for BLUP cross-validation process.](printBLUP.png)

Printting the means of RMSE estimates
-------------------------------------

``` r
options(digits = 4)
RMSEweat = rbind(AMMIweat$RMSPDmean, BLUPweat$RMSPDmean)
RMSEweat = dplyr::mutate(RMSEweat, CROP = "Wheat")
RMSEweat = RMSEweat[order(RMSEweat[,2], decreasing = F),]
#print(RMSEweat)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:left;">
MODEL
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:left;">
CROP
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
11
</td>
<td style="text-align:left;">
BLUP
</td>
<td style="text-align:right;">
0.4141
</td>
<td style="text-align:left;">
Wheat
</td>
</tr>
<tr>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
AMMI4
</td>
<td style="text-align:right;">
0.4320
</td>
<td style="text-align:left;">
Wheat
</td>
</tr>
<tr>
<td style="text-align:left;">
2
</td>
<td style="text-align:left;">
AMMI3
</td>
<td style="text-align:right;">
0.4322
</td>
<td style="text-align:left;">
Wheat
</td>
</tr>
<tr>
<td style="text-align:left;">
3
</td>
<td style="text-align:left;">
AMMI5
</td>
<td style="text-align:right;">
0.4322
</td>
<td style="text-align:left;">
Wheat
</td>
</tr>
<tr>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
AMMI6
</td>
<td style="text-align:right;">
0.4345
</td>
<td style="text-align:left;">
Wheat
</td>
</tr>
<tr>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
AMMI2
</td>
<td style="text-align:right;">
0.4351
</td>
<td style="text-align:left;">
Wheat
</td>
</tr>
<tr>
<td style="text-align:left;">
6
</td>
<td style="text-align:left;">
AMMI7
</td>
<td style="text-align:right;">
0.4369
</td>
<td style="text-align:left;">
Wheat
</td>
</tr>
<tr>
<td style="text-align:left;">
7
</td>
<td style="text-align:left;">
AMMIF
</td>
<td style="text-align:right;">
0.4379
</td>
<td style="text-align:left;">
Wheat
</td>
</tr>
<tr>
<td style="text-align:left;">
8
</td>
<td style="text-align:left;">
AMMI8
</td>
<td style="text-align:right;">
0.4385
</td>
<td style="text-align:left;">
Wheat
</td>
</tr>
<tr>
<td style="text-align:left;">
9
</td>
<td style="text-align:left;">
AMMI0
</td>
<td style="text-align:right;">
0.4399
</td>
<td style="text-align:left;">
Wheat
</td>
</tr>
<tr>
<td style="text-align:left;">
10
</td>
<td style="text-align:left;">
AMMI1
</td>
<td style="text-align:right;">
0.4447
</td>
<td style="text-align:left;">
Wheat
</td>
</tr>
</tbody>
</table>
The results shown above are the mean of 100 RMSE estimates for each tested model and are presented from the most accurate model (smallest RSME value) to the lowest accurate model (highest RMSE value).

Plotting the RMSPD values
-------------------------

> It is also possible plotting the RMSPD values in a boxplot graphic. Let's to do it.

``` r
# binding AMMI and BLUP RMSEs
RMSPDweat = list(RMSPD = rbind(AMMIweat$RMSPD, 
                           BLUPweat$RMSPD))
# Plotting the RMSPD values
plot.validation.AMMIF(RMSPDweat,
                      violin = FALSE,
                      col.boxplot = "gray75")
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

Five statistics are shown in this boxplot. The median, the lower and upper hinges correspond to the first and third quartiles (the 25th and 75th percentiles, respectively). The upper whisker extends from the hinge to the largest value no further than 1.5 × *I**Q**R* from the hinge (where IQR is the inter-quartile range). The lower whisker extends from the hinge to the smallest value at most 1.5 × *I**Q**R* of the hinge. Data beyond the end of the whiskers are considered outlying points. If the condition `violin = TRUE`, a violin plot is added joint with the boxplot. A violin plot is a compact display of a continuous distribution displayed in the same way as a boxplot.

Predicting the yield by traditional AMMI model using replicate-based data
=========================================================================

The function `predict.AMMI` (or the generic function `predict()`) can be used to predict the response variable of a two-way table (for example, the yielding of the *i*th genotype in the *j*th environment) based on the traditional AMMI model. This prediction is based on the number of multiplicative terms used. If `naxis = 0`, only the main effects (AMMI0) are used. In this case, the predicted mean will be the predicted value from OLS estimation. If `naxis = 1` the AMMI1 (with one multiplicative term) is used for predicting the response variable. If `naxis = min(gen-1;env-1)`, the AMMIF is fitted and the predicted value will be the cell mean, i.e., the mean of R replicates of the *i*th genotype in the *j*th environment. The number of axes to be used must be carefully chosen. Procedures based on postdictive success (such as Gollobs's d.f.) or predictive success (such as cross-validation) should be used to do this. This package provides both. `WAAS.AMMI` function compute traditional AMMI analysis showing the number of significant axes. On the other hand, `validation.AMMIF` function provides a cross-validation, estimating the RMSE of all AMMI model family based on re-sampling procedures.

> Predicting the yield of 10 oat cultivars in 16 environments using 5 multiplicative terms.

``` r
model =  WAAS.AMMI(dataset,
                   resp = GY,
                   gen = GEN,
                   env = ENV,
                   rep = BLOCK)
predictoat = predict(model, naxis = 5)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; float: left; margin-right: 10px;">
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
CF2010
</td>
<td style="text-align:left;">
G1
</td>
<td style="text-align:right;">
2.366
</td>
<td style="text-align:right;">
-0.0878
</td>
<td style="text-align:right;">
2.454
</td>
<td style="text-align:right;">
0.04956
</td>
<td style="text-align:right;">
2.503
</td>
<td style="text-align:right;">
2.454
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G10
</td>
<td style="text-align:right;">
1.974
</td>
<td style="text-align:right;">
-0.3620
</td>
<td style="text-align:right;">
2.336
</td>
<td style="text-align:right;">
-0.35636
</td>
<td style="text-align:right;">
1.980
</td>
<td style="text-align:right;">
2.336
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G2
</td>
<td style="text-align:right;">
2.902
</td>
<td style="text-align:right;">
0.3694
</td>
<td style="text-align:right;">
2.533
</td>
<td style="text-align:right;">
0.37700
</td>
<td style="text-align:right;">
2.910
</td>
<td style="text-align:right;">
2.533
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G3
</td>
<td style="text-align:right;">
2.888
</td>
<td style="text-align:right;">
0.1175
</td>
<td style="text-align:right;">
2.771
</td>
<td style="text-align:right;">
0.03069
</td>
<td style="text-align:right;">
2.801
</td>
<td style="text-align:right;">
2.771
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G4
</td>
<td style="text-align:right;">
2.589
</td>
<td style="text-align:right;">
0.0618
</td>
<td style="text-align:right;">
2.527
</td>
<td style="text-align:right;">
-0.05530
</td>
<td style="text-align:right;">
2.472
</td>
<td style="text-align:right;">
2.527
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G5
</td>
<td style="text-align:right;">
2.189
</td>
<td style="text-align:right;">
-0.2068
</td>
<td style="text-align:right;">
2.395
</td>
<td style="text-align:right;">
-0.09272
</td>
<td style="text-align:right;">
2.303
</td>
<td style="text-align:right;">
2.395
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G6
</td>
<td style="text-align:right;">
2.301
</td>
<td style="text-align:right;">
-0.0785
</td>
<td style="text-align:right;">
2.379
</td>
<td style="text-align:right;">
-0.06152
</td>
<td style="text-align:right;">
2.318
</td>
<td style="text-align:right;">
2.379
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G7
</td>
<td style="text-align:right;">
2.774
</td>
<td style="text-align:right;">
0.2337
</td>
<td style="text-align:right;">
2.540
</td>
<td style="text-align:right;">
0.22216
</td>
<td style="text-align:right;">
2.762
</td>
<td style="text-align:right;">
2.540
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G8
</td>
<td style="text-align:right;">
2.899
</td>
<td style="text-align:right;">
0.0302
</td>
<td style="text-align:right;">
2.869
</td>
<td style="text-align:right;">
-0.01589
</td>
<td style="text-align:right;">
2.853
</td>
<td style="text-align:right;">
2.869
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G9
</td>
<td style="text-align:right;">
2.326
</td>
<td style="text-align:right;">
-0.0776
</td>
<td style="text-align:right;">
2.404
</td>
<td style="text-align:right;">
-0.09762
</td>
<td style="text-align:right;">
2.306
</td>
<td style="text-align:right;">
2.404
</td>
</tr>
</tbody>
</table>
Only the first ten values are shown. The following values are presented: **ENV** is the environment; **GEN** is the genotype; **Y** is the response variable; **resOLS** is the residual ($\\hat{z}\_{ij}$) estimated by the Ordinary Least Square (OLS), where $\\hat{z}\_{ij} = y\_{ij} - \\bar{y}\_{i.} - \\bar{y}\_{.j} + \\bar{y}\_{ij}$; **Ypred** is the predicted value by OLS ($\\hat{y}\_{ij} = y\_{ij} -\\hat{z}\_{ij}$); **ResAMMI** is the residual estimated by the AMMI model ($\\hat{a}\_{ij}$) considering the number of multiplicative terms informed in the function (in this case 5), where $\\hat{a}\_{ij} = \\lambda\_1\\alpha\_{i1}\\tau\_{j1}+...+\\lambda\_5\\alpha\_{i5}\\tau\_{j5}$; **YpredAMMI** is the predicted value by AMMI model $\\hat{ya}\_{ij} = \\bar{y}\_{i.} + \\bar{y}\_{.j} - \\bar{y}\_{ij}+\\hat{a}\_{ij}$; and **AMMI0** is the predicted value when no multiplicative terms are used, i.e., $\\hat{y}\_{ij} = \\bar{y}\_{i.} + \\bar{y}\_{.j} - \\bar{y}\_{ij}$.

Estimating the WAAS (based on traditional AMMI model)
=====================================================

The `WAAS.AMMI` function compute the weighted average of absolute scores considering: (i) all principal component axes that are significant at a given probability error level; or (ii) declaring a specific number of axes to be used, according to the following equation:

$$
        WAAS\_i  = 
        \\sum\_{k = 1}^{s} |PCA\_{ik} \\times EP\_k|/ \\sum\_{k = 1}^{s}EP\_k
$$

where *W**A**A**S*<sub>*i*</sub> is the weighted average of absolute scores of the *i*th genotype; *P**C**A*<sub>*i**k*</sub> is the score of the *i*th genotype in the *k*th PCA; and *E**P*<sub>*k*</sub> is the explained variance of the *k*th PCA for *k* = 1, 2, ..,*s*, considering *s* the number of significant PCA, or a a declared number of PCA. The following functions can be used to do this.

Assuming a given probability error for chosing the number of axes
-----------------------------------------------------------------

In this example only PCA with *P*-value ≤ 0.05 will be used in the WAAS estimation.

``` r
library(WAASB)
# Assuming equal weights for productivity and stability
WAAS1 = WAAS.AMMI(dataset,
                  resp = GY,
                  gen = GEN,
                  env = ENV,
                  rep = BLOCK,
                  p.valuePC = 0.05,
                  weight.response = 50,
                  weight.WAAS = 50)

# Priorizing productivity for genotype ranking (WAASB/GY ratio = 30/70)
# no output for this script
WAAS11 = WAAS.AMMI(dataset,
                  resp = GY,
                  gen = GEN,
                  env = ENV,
                  rep = BLOCK,
                  p.valuePC = 0.05,
                  weight.response = 70,
                  weight.WAAS = 30)
```

``` r
options (digits = 4)
# printing the WAASB object
print(WAAS1$anova)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Df
</th>
<th style="text-align:right;">
Sum Sq
</th>
<th style="text-align:right;">
Mean Sq
</th>
<th style="text-align:right;">
F value
</th>
<th style="text-align:left;">
Pr(&gt;F)
</th>
<th style="text-align:left;">
Percent
</th>
<th style="text-align:left;">
Accumul
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
285.9097
</td>
<td style="text-align:right;">
19.0606
</td>
<td style="text-align:right;">
60.239
</td>
<td style="text-align:left;">
0.0000000000000000000000
</td>
<td style="text-align:left;">
.
</td>
<td style="text-align:left;">
.
</td>
</tr>
<tr>
<td style="text-align:left;">
REP(ENV)
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
10.1253
</td>
<td style="text-align:right;">
0.3164
</td>
<td style="text-align:right;">
3.066
</td>
<td style="text-align:left;">
0.0000003177891856598596
</td>
<td style="text-align:left;">
.
</td>
<td style="text-align:left;">
.
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
13.0837
</td>
<td style="text-align:right;">
1.4537
</td>
<td style="text-align:right;">
14.085
</td>
<td style="text-align:left;">
0.0000000000000000008381
</td>
<td style="text-align:left;">
.
</td>
<td style="text-align:left;">
.
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV:GEN
</td>
<td style="text-align:right;">
135
</td>
<td style="text-align:right;">
39.7552
</td>
<td style="text-align:right;">
0.2945
</td>
<td style="text-align:right;">
2.853
</td>
<td style="text-align:left;">
0.0000000000000563937810
</td>
<td style="text-align:left;">
.
</td>
<td style="text-align:left;">
.
</td>
</tr>
<tr>
<td style="text-align:left; padding-left: 2em;" indentlevel="1">
PC1
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
13.0618
</td>
<td style="text-align:right;">
0.5679
</td>
<td style="text-align:right;">
5.500
</td>
<td style="text-align:left;">
0.0000000000000000000000
</td>
<td style="text-align:left;">
32.9
</td>
<td style="text-align:left;">
32.9
</td>
</tr>
<tr>
<td style="text-align:left; padding-left: 2em;" indentlevel="1">
PC2
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
10.6865
</td>
<td style="text-align:right;">
0.5089
</td>
<td style="text-align:right;">
4.930
</td>
<td style="text-align:left;">
0.0000000000000000000000
</td>
<td style="text-align:left;">
26.9
</td>
<td style="text-align:left;">
59.7
</td>
</tr>
<tr>
<td style="text-align:left; padding-left: 2em;" indentlevel="1">
PC3
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
6.8031
</td>
<td style="text-align:right;">
0.3581
</td>
<td style="text-align:right;">
3.470
</td>
<td style="text-align:left;">
0.0000000000000000000000
</td>
<td style="text-align:left;">
17.1
</td>
<td style="text-align:left;">
76.8
</td>
</tr>
<tr>
<td style="text-align:left; padding-left: 2em;" indentlevel="1">
PC4
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
3.6847
</td>
<td style="text-align:right;">
0.2168
</td>
<td style="text-align:right;">
2.100
</td>
<td style="text-align:left;">
0.0071999999999999998029
</td>
<td style="text-align:left;">
9.3
</td>
<td style="text-align:left;">
86.1
</td>
</tr>
<tr>
<td style="text-align:left;font-weight: bold; padding-left: 2em;" indentlevel="1">
PC5
</td>
<td style="text-align:right;font-weight: bold;">
15
</td>
<td style="text-align:right;font-weight: bold;">
2.6701
</td>
<td style="text-align:right;font-weight: bold;">
0.1780
</td>
<td style="text-align:right;font-weight: bold;">
1.720
</td>
<td style="text-align:left;font-weight: bold;">
0.0466000000000000025313
</td>
<td style="text-align:left;font-weight: bold;">
6.7
</td>
<td style="text-align:left;font-weight: bold;">
92.8
</td>
</tr>
<tr>
<td style="text-align:left; padding-left: 2em;" indentlevel="1">
PC6
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
1.4540
</td>
<td style="text-align:right;">
0.1118
</td>
<td style="text-align:right;">
1.080
</td>
<td style="text-align:left;">
0.3760999999999999898748
</td>
<td style="text-align:left;">
3.7
</td>
<td style="text-align:left;">
96.5
</td>
</tr>
<tr>
<td style="text-align:left; padding-left: 2em;" indentlevel="1">
PC7
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.7266
</td>
<td style="text-align:right;">
0.0660
</td>
<td style="text-align:right;">
0.640
</td>
<td style="text-align:left;">
0.7939000000000000500933
</td>
<td style="text-align:left;">
1.8
</td>
<td style="text-align:left;">
98.3
</td>
</tr>
<tr>
<td style="text-align:left; padding-left: 2em;" indentlevel="1">
PC8
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.4652
</td>
<td style="text-align:right;">
0.0517
</td>
<td style="text-align:right;">
0.500
</td>
<td style="text-align:left;">
0.8739999999999999991118
</td>
<td style="text-align:left;">
1.2
</td>
<td style="text-align:left;">
99.5
</td>
</tr>
<tr>
<td style="text-align:left; padding-left: 2em;" indentlevel="1">
PC9
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.2032
</td>
<td style="text-align:right;">
0.0290
</td>
<td style="text-align:right;">
0.280
</td>
<td style="text-align:left;">
0.9615000000000000213163
</td>
<td style="text-align:left;">
0.5
</td>
<td style="text-align:left;">
100
</td>
</tr>
<tr>
<td style="text-align:left;">
Residuals
</td>
<td style="text-align:right;">
288
</td>
<td style="text-align:right;">
29.7248
</td>
<td style="text-align:right;">
0.1032
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
.
</td>
<td style="text-align:left;">
.
</td>
</tr>
<tr>
<td style="text-align:left;">
Total
</td>
<td style="text-align:right;">
479
</td>
<td style="text-align:right;">
378.5987
</td>
<td style="text-align:right;">
0.7904
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
.
</td>
<td style="text-align:left;">
.
</td>
</tr>
</tbody>
</table>
The above table is the traditional AMMI analysis. Nine principal component axis were fitted and five were significant at 5% probability error. Based on this result, the AMMI5 model would be the best model to predict the yielding of the genotypes in the studied environments. This is confirmed by the cross-validation in the [section 3](#predictive-accuracy). The AMMI model with the smallest RMSE was the AMMI5. The predicted values (first ten observations) can be seen in [section 4](#predicting-the-yield-using-traditional-ammi-model)

> printing the WAASB object

``` r
options (digits = 4)
print(WAAS1$model[, c(1:3,13:17, 21:22)])
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
type
</th>
<th style="text-align:left;">
Code
</th>
<th style="text-align:right;">
Y
</th>
<th style="text-align:right;">
WAAS
</th>
<th style="text-align:right;">
PctResp
</th>
<th style="text-align:right;">
PctWAAS
</th>
<th style="text-align:right;">
OrResp
</th>
<th style="text-align:right;">
OrWAAS
</th>
<th style="text-align:right;">
WAASY
</th>
<th style="text-align:right;">
OrWAASY
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G1
</td>
<td style="text-align:right;">
2.624
</td>
<td style="text-align:right;">
0.1778
</td>
<td style="text-align:right;">
86.32
</td>
<td style="text-align:right;">
99.00
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
92.66
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G10
</td>
<td style="text-align:right;">
2.506
</td>
<td style="text-align:right;">
0.5003
</td>
<td style="text-align:right;">
82.46
</td>
<td style="text-align:right;">
97.19
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
89.82
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G2
</td>
<td style="text-align:right;">
2.703
</td>
<td style="text-align:right;">
0.3124
</td>
<td style="text-align:right;">
88.93
</td>
<td style="text-align:right;">
98.24
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
93.59
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G3
</td>
<td style="text-align:right;">
2.941
</td>
<td style="text-align:right;">
0.2037
</td>
<td style="text-align:right;">
96.77
</td>
<td style="text-align:right;">
98.85
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
97.81
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G4
</td>
<td style="text-align:right;">
2.697
</td>
<td style="text-align:right;">
0.3682
</td>
<td style="text-align:right;">
88.74
</td>
<td style="text-align:right;">
97.93
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
93.33
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G5
</td>
<td style="text-align:right;">
2.566
</td>
<td style="text-align:right;">
0.2143
</td>
<td style="text-align:right;">
84.41
</td>
<td style="text-align:right;">
98.79
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
91.60
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G6
</td>
<td style="text-align:right;">
2.549
</td>
<td style="text-align:right;">
0.2439
</td>
<td style="text-align:right;">
83.88
</td>
<td style="text-align:right;">
98.63
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
91.25
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G7
</td>
<td style="text-align:right;">
2.710
</td>
<td style="text-align:right;">
0.3843
</td>
<td style="text-align:right;">
89.17
</td>
<td style="text-align:right;">
97.84
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
93.50
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G8
</td>
<td style="text-align:right;">
3.039
</td>
<td style="text-align:right;">
0.2465
</td>
<td style="text-align:right;">
100.00
</td>
<td style="text-align:right;">
98.61
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
99.31
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G9
</td>
<td style="text-align:right;">
2.574
</td>
<td style="text-align:right;">
0.4473
</td>
<td style="text-align:right;">
84.68
</td>
<td style="text-align:right;">
97.48
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
91.08
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:right;">
2.521
</td>
<td style="text-align:right;">
0.1787
</td>
<td style="text-align:right;">
62.02
</td>
<td style="text-align:right;">
98.47
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
80.24
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2011
</td>
<td style="text-align:right;">
3.180
</td>
<td style="text-align:right;">
0.2086
</td>
<td style="text-align:right;">
78.24
</td>
<td style="text-align:right;">
98.21
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
88.22
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2012
</td>
<td style="text-align:right;">
4.064
</td>
<td style="text-align:right;">
0.3627
</td>
<td style="text-align:right;">
100.00
</td>
<td style="text-align:right;">
96.89
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
98.44
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2013
</td>
<td style="text-align:right;">
3.675
</td>
<td style="text-align:right;">
0.2564
</td>
<td style="text-align:right;">
90.43
</td>
<td style="text-align:right;">
97.80
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
94.11
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2014
</td>
<td style="text-align:right;">
2.507
</td>
<td style="text-align:right;">
0.2429
</td>
<td style="text-align:right;">
61.69
</td>
<td style="text-align:right;">
97.91
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
79.80
</td>
<td style="text-align:right;">
11
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2015
</td>
<td style="text-align:right;">
3.107
</td>
<td style="text-align:right;">
0.3486
</td>
<td style="text-align:right;">
76.45
</td>
<td style="text-align:right;">
97.01
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
86.73
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2016
</td>
<td style="text-align:right;">
3.910
</td>
<td style="text-align:right;">
0.2115
</td>
<td style="text-align:right;">
96.22
</td>
<td style="text-align:right;">
98.18
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
97.20
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2017
</td>
<td style="text-align:right;">
2.663
</td>
<td style="text-align:right;">
0.1167
</td>
<td style="text-align:right;">
65.53
</td>
<td style="text-align:right;">
99.00
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
82.26
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2010
</td>
<td style="text-align:right;">
1.989
</td>
<td style="text-align:right;">
0.2973
</td>
<td style="text-align:right;">
48.94
</td>
<td style="text-align:right;">
97.45
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
73.19
</td>
<td style="text-align:right;">
13
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2011
</td>
<td style="text-align:right;">
2.536
</td>
<td style="text-align:right;">
0.2263
</td>
<td style="text-align:right;">
62.41
</td>
<td style="text-align:right;">
98.06
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
80.23
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2012
</td>
<td style="text-align:right;">
3.057
</td>
<td style="text-align:right;">
0.4128
</td>
<td style="text-align:right;">
75.21
</td>
<td style="text-align:right;">
96.45
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
85.83
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2013
</td>
<td style="text-align:right;">
2.175
</td>
<td style="text-align:right;">
0.1623
</td>
<td style="text-align:right;">
53.52
</td>
<td style="text-align:right;">
98.61
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
76.06
</td>
<td style="text-align:right;">
12
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2014
</td>
<td style="text-align:right;">
1.368
</td>
<td style="text-align:right;">
0.1164
</td>
<td style="text-align:right;">
33.66
</td>
<td style="text-align:right;">
99.00
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
66.33
</td>
<td style="text-align:right;">
16
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2015
</td>
<td style="text-align:right;">
1.609
</td>
<td style="text-align:right;">
0.1663
</td>
<td style="text-align:right;">
39.58
</td>
<td style="text-align:right;">
98.57
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
69.07
</td>
<td style="text-align:right;">
15
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2016
</td>
<td style="text-align:right;">
2.910
</td>
<td style="text-align:right;">
0.3235
</td>
<td style="text-align:right;">
71.59
</td>
<td style="text-align:right;">
97.22
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
84.41
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2017
</td>
<td style="text-align:right;">
1.782
</td>
<td style="text-align:right;">
0.1588
</td>
<td style="text-align:right;">
43.84
</td>
<td style="text-align:right;">
98.64
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
71.24
</td>
<td style="text-align:right;">
14
</td>
</tr>
</tbody>
</table>
In this example, the scores of the nine PCA were not shown. The output generated by the `WAAS.AMMI` function shows the following results: **type**, genotype (GEN) or environment (ENV); **Code**, the code attributed to each level of the factors; **Y**, the response variable (in this case the grain yield); **WAAS** the weighted average of the absolute scores, estimated with all PCA axes with *P*-value ≤ 0.05; **PctWAAS** and **PctResp** that are the percentage values for the WAAS and Y, respectively; **OrResp** and **OrWAAS** that are the ranks attributed to the genotype and environment regarding the Y or WAAS, respectively; **WAASY** is the weighted average of absolute scores and response variable. In this case, considering equal weights for PctResp and PctWAAS, the WAASY for G1 is estimated by: *W**A**A**S*<sub>*G*1</sub> = \[(86.32 × 50)+(98.88 × 50)\]/50 + 50 = 92.60. Then the \*\*OrWAASY\* is the rank for the WAASY value. The genotype (or environment) with the largest WAASY value has the first ranked. See [Estimating the WAASBY index](#estimating-the-waasby-index) for a detailed explanation.

> Biplot WAAS x GY

``` r
plot.scores(WAAS1,
            type = 3)
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

The biplot above shows the coordinates of the genotypes and environments regarding the grain yield and the WAAS values. It is in fact the plot of the columns Y and WAAS from the above table. A detailed discussion on the interpretation of this biplot can be seen in [section 6.2.3](#biplot-type-3-gy-x-waasb).

Declaring a specific number of axis to be used
----------------------------------------------

Here it is possible to inform a specific number of multiplicative terms to be used in the WAAS estimation. In this example, the number of terms informed is used independently of its significance. Let us, for the moment, assume that two multiplicative terms should be used.

> Declaring that five PCA must be used to compute the WAAS

``` r
WAAS2 = WAAS.AMMI(dataset,
                  resp = GY,
                  gen = GEN,
                  env = ENV,
                  rep = BLOCK,
                  weight.response = 50,
                  weight.WAAS = 50)
```

``` r
options (digits = 4)
# printing the WAASB object
print(WAAS2$model[, c(1:3,13:17, 21:22)])
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
type
</th>
<th style="text-align:left;">
Code
</th>
<th style="text-align:right;">
Y
</th>
<th style="text-align:right;">
WAAS
</th>
<th style="text-align:right;">
PctResp
</th>
<th style="text-align:right;">
PctWAAS
</th>
<th style="text-align:right;">
OrResp
</th>
<th style="text-align:right;">
OrWAAS
</th>
<th style="text-align:right;">
WAASY
</th>
<th style="text-align:right;">
OrWAASY
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G1
</td>
<td style="text-align:right;">
2.624
</td>
<td style="text-align:right;">
0.1778
</td>
<td style="text-align:right;">
86.32
</td>
<td style="text-align:right;">
99.00
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
92.66
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G10
</td>
<td style="text-align:right;">
2.506
</td>
<td style="text-align:right;">
0.5003
</td>
<td style="text-align:right;">
82.46
</td>
<td style="text-align:right;">
97.19
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
89.82
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G2
</td>
<td style="text-align:right;">
2.703
</td>
<td style="text-align:right;">
0.3124
</td>
<td style="text-align:right;">
88.93
</td>
<td style="text-align:right;">
98.24
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
93.59
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G3
</td>
<td style="text-align:right;">
2.941
</td>
<td style="text-align:right;">
0.2037
</td>
<td style="text-align:right;">
96.77
</td>
<td style="text-align:right;">
98.85
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
97.81
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G4
</td>
<td style="text-align:right;">
2.697
</td>
<td style="text-align:right;">
0.3682
</td>
<td style="text-align:right;">
88.74
</td>
<td style="text-align:right;">
97.93
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
93.33
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G5
</td>
<td style="text-align:right;">
2.566
</td>
<td style="text-align:right;">
0.2143
</td>
<td style="text-align:right;">
84.41
</td>
<td style="text-align:right;">
98.79
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
91.60
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G6
</td>
<td style="text-align:right;">
2.549
</td>
<td style="text-align:right;">
0.2439
</td>
<td style="text-align:right;">
83.88
</td>
<td style="text-align:right;">
98.63
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
91.25
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G7
</td>
<td style="text-align:right;">
2.710
</td>
<td style="text-align:right;">
0.3843
</td>
<td style="text-align:right;">
89.17
</td>
<td style="text-align:right;">
97.84
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
93.50
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G8
</td>
<td style="text-align:right;">
3.039
</td>
<td style="text-align:right;">
0.2465
</td>
<td style="text-align:right;">
100.00
</td>
<td style="text-align:right;">
98.61
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
99.31
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G9
</td>
<td style="text-align:right;">
2.574
</td>
<td style="text-align:right;">
0.4473
</td>
<td style="text-align:right;">
84.68
</td>
<td style="text-align:right;">
97.48
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
91.08
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:right;">
2.521
</td>
<td style="text-align:right;">
0.1787
</td>
<td style="text-align:right;">
62.02
</td>
<td style="text-align:right;">
98.47
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
80.24
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2011
</td>
<td style="text-align:right;">
3.180
</td>
<td style="text-align:right;">
0.2086
</td>
<td style="text-align:right;">
78.24
</td>
<td style="text-align:right;">
98.21
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
88.22
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2012
</td>
<td style="text-align:right;">
4.064
</td>
<td style="text-align:right;">
0.3627
</td>
<td style="text-align:right;">
100.00
</td>
<td style="text-align:right;">
96.89
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
98.44
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2013
</td>
<td style="text-align:right;">
3.675
</td>
<td style="text-align:right;">
0.2564
</td>
<td style="text-align:right;">
90.43
</td>
<td style="text-align:right;">
97.80
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
94.11
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2014
</td>
<td style="text-align:right;">
2.507
</td>
<td style="text-align:right;">
0.2429
</td>
<td style="text-align:right;">
61.69
</td>
<td style="text-align:right;">
97.91
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
79.80
</td>
<td style="text-align:right;">
11
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2015
</td>
<td style="text-align:right;">
3.107
</td>
<td style="text-align:right;">
0.3486
</td>
<td style="text-align:right;">
76.45
</td>
<td style="text-align:right;">
97.01
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
86.73
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2016
</td>
<td style="text-align:right;">
3.910
</td>
<td style="text-align:right;">
0.2115
</td>
<td style="text-align:right;">
96.22
</td>
<td style="text-align:right;">
98.18
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
97.20
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2017
</td>
<td style="text-align:right;">
2.663
</td>
<td style="text-align:right;">
0.1167
</td>
<td style="text-align:right;">
65.53
</td>
<td style="text-align:right;">
99.00
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
82.26
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2010
</td>
<td style="text-align:right;">
1.989
</td>
<td style="text-align:right;">
0.2973
</td>
<td style="text-align:right;">
48.94
</td>
<td style="text-align:right;">
97.45
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
73.19
</td>
<td style="text-align:right;">
13
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2011
</td>
<td style="text-align:right;">
2.536
</td>
<td style="text-align:right;">
0.2263
</td>
<td style="text-align:right;">
62.41
</td>
<td style="text-align:right;">
98.06
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
80.23
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2012
</td>
<td style="text-align:right;">
3.057
</td>
<td style="text-align:right;">
0.4128
</td>
<td style="text-align:right;">
75.21
</td>
<td style="text-align:right;">
96.45
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
85.83
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2013
</td>
<td style="text-align:right;">
2.175
</td>
<td style="text-align:right;">
0.1623
</td>
<td style="text-align:right;">
53.52
</td>
<td style="text-align:right;">
98.61
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
76.06
</td>
<td style="text-align:right;">
12
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2014
</td>
<td style="text-align:right;">
1.368
</td>
<td style="text-align:right;">
0.1164
</td>
<td style="text-align:right;">
33.66
</td>
<td style="text-align:right;">
99.00
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
66.33
</td>
<td style="text-align:right;">
16
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2015
</td>
<td style="text-align:right;">
1.609
</td>
<td style="text-align:right;">
0.1663
</td>
<td style="text-align:right;">
39.58
</td>
<td style="text-align:right;">
98.57
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
69.07
</td>
<td style="text-align:right;">
15
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2016
</td>
<td style="text-align:right;">
2.910
</td>
<td style="text-align:right;">
0.3235
</td>
<td style="text-align:right;">
71.59
</td>
<td style="text-align:right;">
97.22
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
84.41
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2017
</td>
<td style="text-align:right;">
1.782
</td>
<td style="text-align:right;">
0.1588
</td>
<td style="text-align:right;">
43.84
</td>
<td style="text-align:right;">
98.64
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
71.24
</td>
<td style="text-align:right;">
14
</td>
</tr>
</tbody>
</table>
The only difference of this output compared tho those from the [section 5.1](#assuming-a-given-probability-error-for-chosing-the-number-of-axes) is that here we declared that seven PCA axes should be used for computing the WAAS value. Thus, only the values of WAAS, OrWAAS, WAASY and OrWAASY are changed. The magnitude of changes, however are too small. This is explained because only two multiplicative terms were included here compared to the previous examples. In addition, these multiplicative terms (PCA6 and PCA7) explained only 3.7% and 1.8%, respectively of the GEI pattern. See section [section 5.1](#assuming-a-given-probability-error-for-chosing-the-number-of-axes).

> Biplot WAAS x GY

``` r
plot.scores(WAAS2,
            type = 3)
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-22-1.png" style="display: block; margin: auto;" />

In this biplot, we can see that the changes in the coordinates were also too small.

Estimating the WAASB (based on SVD of BLUP-interaction effects)
===============================================================

The `WAASB` function compute the weighted average of absolute scores considering all principal component axes from Singular Value Decomposition of BLUP-interaction effects matrix, considering the following model:

$$
        WAASB\_i  = 
        \\sum\_{k = 1}^{p} |PCA\_{ik} \\times EP\_k|/ \\sum\_{k = 1}^{p}EP\_k
$$

where *W**A**A**S**B*<sub>*i*</sub> is the weighted average of absolute scores of the *i*th genotype; *P**C**A*<sub>*i**k*</sub> is the scores of the *i*th genotype in the *k*th PCA; and *E**P*<sub>*k*</sub> is the explained variance of the *k*th PCA for *k* = 1, 2, ..,*p*, and *p* = *m**i**n*(*G* − 1; *E* − 1).

As the matrix of genotypic effects (free from error) is used, all multiplicative terms can be used without a noise-adding effect. The main advantage of `WAASB` over `WAAS` function is that, by using a mixed model, random effects can be included in the model. In addition, unbalanced data sets can also be modeled.

``` r
# Assuming equal weights for productivity and stability
WAASB = WAASB(dataset,
              resp = GY,
              gen = GEN,
              env = ENV,
              rep = BLOCK,
              prob = 0.95,
              weight.response = 50,
              weight.WAAS = 50)

# Priorizing productivity for genotype ranking (WAASB/GY ratio = 30/70)
# no output for this script
WAASB2 = WAASB(dataset,
              resp = GY,
              gen = GEN,
              env = ENV,
              rep = BLOCK,
              prob = 0.95,
              weight.response = 30,
              weight.WAAS = 70)
```

Printing the model outputs
--------------------------

> Likelihood ratio tests on linear mixed-effect models

``` r
print(WAASB$LRT)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:center;">
npar
</th>
<th style="text-align:center;">
logLik
</th>
<th style="text-align:center;">
AIC
</th>
<th style="text-align:center;">
LRT
</th>
<th style="text-align:center;">
Df
</th>
<th style="text-align:center;">
Pr(&gt;Chisq)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Complete
</td>
<td style="text-align:center;">
51
</td>
<td style="text-align:center;">
-260.38
</td>
<td style="text-align:center;">
622.77
</td>
<td style="text-align:center;">
NA
</td>
<td style="text-align:center;">
NA
</td>
<td style="text-align:center;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Genotype
</td>
<td style="text-align:center;">
50
</td>
<td style="text-align:center;">
-269.04
</td>
<td style="text-align:center;">
638.08
</td>
<td style="text-align:center;">
17.305
</td>
<td style="text-align:center;">
1
</td>
<td style="text-align:center;">
3e-05
</td>
</tr>
<tr>
<td style="text-align:left;">
Gen vs Env
</td>
<td style="text-align:center;">
50
</td>
<td style="text-align:center;">
-287.89
</td>
<td style="text-align:center;">
675.78
</td>
<td style="text-align:center;">
55.005
</td>
<td style="text-align:center;">
1
</td>
<td style="text-align:center;">
0e+00
</td>
</tr>
</tbody>
</table>
The output `LRT` contains the Likelihood Ratio Tests on the random effects of the linear mixed-effect.

> Variance components and some parameters

``` r
options (digits = 4)
print(WAASB$ESTIMATES)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
Parameters
</th>
<th style="text-align:left;">
Values
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
GEI variance
</td>
<td style="text-align:left;">
0.063757 (33.36% of fenotypic variance.)
</td>
</tr>
<tr>
<td style="text-align:left;">
Genotypic variance
</td>
<td style="text-align:left;">
0.024151 (12.64% of fenotypic variance.)
</td>
</tr>
<tr>
<td style="text-align:left;">
Residual variance
</td>
<td style="text-align:left;">
0.103211 (54% of fenotypic variance.)
</td>
</tr>
<tr>
<td style="text-align:left;">
Phenotypic variance
</td>
<td style="text-align:left;">
0.191119577648797
</td>
</tr>
<tr>
<td style="text-align:left;">
Heritability
</td>
<td style="text-align:left;">
0.126366950498645
</td>
</tr>
<tr>
<td style="text-align:left;">
GEIr2
</td>
<td style="text-align:left;">
0.3335996364357
</td>
</tr>
<tr>
<td style="text-align:left;">
Heribatility of means
</td>
<td style="text-align:left;">
0.797430714258779
</td>
</tr>
<tr>
<td style="text-align:left;">
Accuracy
</td>
<td style="text-align:left;">
0.892989761564364
</td>
</tr>
<tr>
<td style="text-align:left;">
rge
</td>
<td style="text-align:left;">
0.381853269660653
</td>
</tr>
<tr>
<td style="text-align:left;">
CVg
</td>
<td style="text-align:left;">
5.77536605468609
</td>
</tr>
<tr>
<td style="text-align:left;">
CVr
</td>
<td style="text-align:left;">
11.9391409407191
</td>
</tr>
<tr>
<td style="text-align:left;">
CV ratio
</td>
<td style="text-align:left;">
0.483733803241143
</td>
</tr>
</tbody>
</table>
In the output `ESTIMATES` some important parameters is shown. **GEIr2** is the coefficient of determination of the interaction effects; **rge** is the genotype-environment correlation; **CVg** is the the genotypic coefficient of variation; **CVr** is the residual coefficient of variation; **CV ratio** is the ratio between genotypic and residual coefficient of variation.

> A detailed infromation regarding the analyzed experiment.

``` r
options (digits = 4)
print(WAASB$Details)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
Parameters
</th>
<th style="text-align:left;">
Values
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
WgtResponse
</td>
<td style="text-align:left;">
50
</td>
</tr>
<tr>
<td style="text-align:left;">
WgtWAAS
</td>
<td style="text-align:left;">
50
</td>
</tr>
<tr>
<td style="text-align:left;">
Ngen
</td>
<td style="text-align:left;">
10
</td>
</tr>
<tr>
<td style="text-align:left;">
Nenv
</td>
<td style="text-align:left;">
16
</td>
</tr>
<tr>
<td style="text-align:left;">
OVmean
</td>
<td style="text-align:left;">
2.6909
</td>
</tr>
<tr>
<td style="text-align:left;">
Min
</td>
<td style="text-align:left;">
0.899 (Genotype G10 in SF2014 )
</td>
</tr>
<tr>
<td style="text-align:left;">
Max
</td>
<td style="text-align:left;">
4.812 (Genotype G8 in CF2016 )
</td>
</tr>
<tr>
<td style="text-align:left;">
MinENV
</td>
<td style="text-align:left;">
Environment SF2014 (1.3682)
</td>
</tr>
<tr>
<td style="text-align:left;">
MaxENV
</td>
<td style="text-align:left;">
Environment CF2012 (4.0643)
</td>
</tr>
<tr>
<td style="text-align:left;">
MinGEN
</td>
<td style="text-align:left;">
Genotype G10 (2.5061)
</td>
</tr>
<tr>
<td style="text-align:left;">
MaxGEN
</td>
<td style="text-align:left;">
Genotype G8 (3.0393)
</td>
</tr>
</tbody>
</table>
The following information are showed in `Details` output. **WgtResponse** is the weight for the response variable in estimating WAASB; **WgtWAAS** the weight for stability; **Ngen** is the number of genotypes; **Nenv** is the number of environments; **OVmean** is the overall mean; **Min** is the minimum value observed (returning the genotype and environment); **Max** is the maximum observed; **MinENV** is the environment with the lower mean; **MaxENV** is the environment with the largest mean observed; **MinGEN** is the genotype with the lower mean; **MaxGEN** is the genotype with the largest mean.

> printing the WAASB object

``` r
options (digits = 4)
print(WAASB$model[, c(1:3,13:17, 21:22)])
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
type
</th>
<th style="text-align:left;">
Code
</th>
<th style="text-align:right;">
Y
</th>
<th style="text-align:right;">
WAASB
</th>
<th style="text-align:right;">
PctResp
</th>
<th style="text-align:right;">
PctWAASB
</th>
<th style="text-align:right;">
OrResp
</th>
<th style="text-align:right;">
OrWAASB
</th>
<th style="text-align:right;">
WAASY
</th>
<th style="text-align:right;">
OrWAASY
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G1
</td>
<td style="text-align:right;">
2.624
</td>
<td style="text-align:right;">
0.1586
</td>
<td style="text-align:right;">
86.32
</td>
<td style="text-align:right;">
98.88
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
92.60
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G10
</td>
<td style="text-align:right;">
2.506
</td>
<td style="text-align:right;">
0.4496
</td>
<td style="text-align:right;">
82.46
</td>
<td style="text-align:right;">
96.82
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
89.64
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G2
</td>
<td style="text-align:right;">
2.703
</td>
<td style="text-align:right;">
0.2135
</td>
<td style="text-align:right;">
88.93
</td>
<td style="text-align:right;">
98.49
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
93.71
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G3
</td>
<td style="text-align:right;">
2.941
</td>
<td style="text-align:right;">
0.1415
</td>
<td style="text-align:right;">
96.77
</td>
<td style="text-align:right;">
99.00
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
97.88
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G4
</td>
<td style="text-align:right;">
2.697
</td>
<td style="text-align:right;">
0.2996
</td>
<td style="text-align:right;">
88.74
</td>
<td style="text-align:right;">
97.88
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
93.31
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G5
</td>
<td style="text-align:right;">
2.566
</td>
<td style="text-align:right;">
0.1837
</td>
<td style="text-align:right;">
84.41
</td>
<td style="text-align:right;">
98.70
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
91.56
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G6
</td>
<td style="text-align:right;">
2.549
</td>
<td style="text-align:right;">
0.1627
</td>
<td style="text-align:right;">
83.88
</td>
<td style="text-align:right;">
98.85
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
91.37
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G7
</td>
<td style="text-align:right;">
2.710
</td>
<td style="text-align:right;">
0.2994
</td>
<td style="text-align:right;">
89.17
</td>
<td style="text-align:right;">
97.88
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
93.53
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G8
</td>
<td style="text-align:right;">
3.039
</td>
<td style="text-align:right;">
0.2063
</td>
<td style="text-align:right;">
100.00
</td>
<td style="text-align:right;">
98.54
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
99.27
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
GEN
</td>
<td style="text-align:left;">
G9
</td>
<td style="text-align:right;">
2.574
</td>
<td style="text-align:right;">
0.4017
</td>
<td style="text-align:right;">
84.68
</td>
<td style="text-align:right;">
97.16
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
90.92
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:right;">
2.521
</td>
<td style="text-align:right;">
0.1750
</td>
<td style="text-align:right;">
62.02
</td>
<td style="text-align:right;">
98.02
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
80.02
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2011
</td>
<td style="text-align:right;">
3.180
</td>
<td style="text-align:right;">
0.1435
</td>
<td style="text-align:right;">
78.24
</td>
<td style="text-align:right;">
98.38
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
88.31
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2012
</td>
<td style="text-align:right;">
4.064
</td>
<td style="text-align:right;">
0.2494
</td>
<td style="text-align:right;">
100.00
</td>
<td style="text-align:right;">
97.18
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
98.59
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2013
</td>
<td style="text-align:right;">
3.675
</td>
<td style="text-align:right;">
0.2617
</td>
<td style="text-align:right;">
90.43
</td>
<td style="text-align:right;">
97.04
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
93.73
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2014
</td>
<td style="text-align:right;">
2.507
</td>
<td style="text-align:right;">
0.2304
</td>
<td style="text-align:right;">
61.69
</td>
<td style="text-align:right;">
97.39
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
79.54
</td>
<td style="text-align:right;">
11
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2015
</td>
<td style="text-align:right;">
3.107
</td>
<td style="text-align:right;">
0.2362
</td>
<td style="text-align:right;">
76.45
</td>
<td style="text-align:right;">
97.33
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
86.89
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2016
</td>
<td style="text-align:right;">
3.910
</td>
<td style="text-align:right;">
0.2109
</td>
<td style="text-align:right;">
96.22
</td>
<td style="text-align:right;">
97.61
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
96.91
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
CF2017
</td>
<td style="text-align:right;">
2.663
</td>
<td style="text-align:right;">
0.0884
</td>
<td style="text-align:right;">
65.53
</td>
<td style="text-align:right;">
99.00
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
82.26
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2010
</td>
<td style="text-align:right;">
1.989
</td>
<td style="text-align:right;">
0.1974
</td>
<td style="text-align:right;">
48.94
</td>
<td style="text-align:right;">
97.77
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
73.35
</td>
<td style="text-align:right;">
13
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2011
</td>
<td style="text-align:right;">
2.536
</td>
<td style="text-align:right;">
0.1660
</td>
<td style="text-align:right;">
62.41
</td>
<td style="text-align:right;">
98.12
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
80.26
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2012
</td>
<td style="text-align:right;">
3.057
</td>
<td style="text-align:right;">
0.3628
</td>
<td style="text-align:right;">
75.21
</td>
<td style="text-align:right;">
95.89
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
85.55
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2013
</td>
<td style="text-align:right;">
2.175
</td>
<td style="text-align:right;">
0.1673
</td>
<td style="text-align:right;">
53.52
</td>
<td style="text-align:right;">
98.11
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
75.81
</td>
<td style="text-align:right;">
12
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2014
</td>
<td style="text-align:right;">
1.368
</td>
<td style="text-align:right;">
0.0985
</td>
<td style="text-align:right;">
33.66
</td>
<td style="text-align:right;">
98.89
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
66.28
</td>
<td style="text-align:right;">
16
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2015
</td>
<td style="text-align:right;">
1.609
</td>
<td style="text-align:right;">
0.1776
</td>
<td style="text-align:right;">
39.58
</td>
<td style="text-align:right;">
97.99
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
68.78
</td>
<td style="text-align:right;">
15
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2016
</td>
<td style="text-align:right;">
2.910
</td>
<td style="text-align:right;">
0.3004
</td>
<td style="text-align:right;">
71.59
</td>
<td style="text-align:right;">
96.60
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
84.09
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
SF2017
</td>
<td style="text-align:right;">
1.782
</td>
<td style="text-align:right;">
0.1091
</td>
<td style="text-align:right;">
43.84
</td>
<td style="text-align:right;">
98.77
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
71.30
</td>
<td style="text-align:right;">
14
</td>
</tr>
</tbody>
</table>
This output generated by the `WAASB` function is very similar to those shown in the [sections 5.1](#assuming-a-given-probability-error-for-chosing-the-number-of-axes) and [5.2](#declaring-a-specific-number-of-axis-to-be-used). The main difference here, is that the singular value decomposition is based on the BLUP interaction effect matrix and the **WAASB** in this output is the weighted average of the absolute scores, estimated with all estimated PCA axes instead **WAAS** that is estimated considering only PCA axes with *P*-value ≤ 0.05.

> printing the estimated BLUP for genotypes

``` r
options(digits = 4)
print(WAASB$BLUPgen[1:10,])
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:right;">
Rank
</th>
<th style="text-align:left;">
GEN
</th>
<th style="text-align:right;">
BLUPg
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
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
G8
</td>
<td style="text-align:right;">
0.2778
</td>
<td style="text-align:right;">
2.969
</td>
<td style="text-align:right;">
2.868
</td>
<td style="text-align:right;">
3.070
</td>
</tr>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
G3
</td>
<td style="text-align:right;">
0.1994
</td>
<td style="text-align:right;">
2.890
</td>
<td style="text-align:right;">
2.789
</td>
<td style="text-align:right;">
2.991
</td>
</tr>
<tr>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
G7
</td>
<td style="text-align:right;">
0.0153
</td>
<td style="text-align:right;">
2.706
</td>
<td style="text-align:right;">
2.605
</td>
<td style="text-align:right;">
2.807
</td>
</tr>
<tr>
<td style="text-align:right;">
4
</td>
<td style="text-align:left;">
G2
</td>
<td style="text-align:right;">
0.0095
</td>
<td style="text-align:right;">
2.700
</td>
<td style="text-align:right;">
2.599
</td>
<td style="text-align:right;">
2.801
</td>
</tr>
<tr>
<td style="text-align:right;">
5
</td>
<td style="text-align:left;">
G4
</td>
<td style="text-align:right;">
0.0049
</td>
<td style="text-align:right;">
2.696
</td>
<td style="text-align:right;">
2.595
</td>
<td style="text-align:right;">
2.797
</td>
</tr>
<tr>
<td style="text-align:right;">
6
</td>
<td style="text-align:left;">
G1
</td>
<td style="text-align:right;">
-0.0536
</td>
<td style="text-align:right;">
2.637
</td>
<td style="text-align:right;">
2.536
</td>
<td style="text-align:right;">
2.738
</td>
</tr>
<tr>
<td style="text-align:right;">
7
</td>
<td style="text-align:left;">
G9
</td>
<td style="text-align:right;">
-0.0934
</td>
<td style="text-align:right;">
2.597
</td>
<td style="text-align:right;">
2.497
</td>
<td style="text-align:right;">
2.698
</td>
</tr>
<tr>
<td style="text-align:right;">
8
</td>
<td style="text-align:left;">
G5
</td>
<td style="text-align:right;">
-0.0999
</td>
<td style="text-align:right;">
2.591
</td>
<td style="text-align:right;">
2.490
</td>
<td style="text-align:right;">
2.692
</td>
</tr>
<tr>
<td style="text-align:right;">
9
</td>
<td style="text-align:left;">
G6
</td>
<td style="text-align:right;">
-0.1128
</td>
<td style="text-align:right;">
2.578
</td>
<td style="text-align:right;">
2.477
</td>
<td style="text-align:right;">
2.679
</td>
</tr>
<tr>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
G10
</td>
<td style="text-align:right;">
-0.1473
</td>
<td style="text-align:right;">
2.543
</td>
<td style="text-align:right;">
2.443
</td>
<td style="text-align:right;">
2.644
</td>
</tr>
</tbody>
</table>
> plotting the estimated BLUP for genotypes

``` r
# No file exported
plot.blup(WAASB)
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-34-1.png" style="display: block; margin: auto;" />

This output shows the predicted means for genotypes. **BLUPg** is the genotypic effect $(\\hat{g}\_{i})$ estimated by $\\hat{g}\_{i} = h\_g^2(\\bar{y}\_{i.}-\\bar{y}\_{..})$ where *h*<sub>*g*</sub><sup>2</sup> is the shrinkage effect for genotype. **Predicted** is the predicted mean estimated by $\\hat{g}\_{i}+\\mu$ where is the grand mean. **LL** and **UL** are the lower and upper limits, respectively, estimated by $(\\hat{g}\_{i}+\\mu)\\pm{CI}$. *C**I* is the confidence interval for BLUP prediction assuming a given probability error, where $CI = t\\times\\sqrt{((1-Ac)\\times{GV)}}$ where *t* is the Student's *t* value for a two-tailed t test at a given probability error; *A**c* is the genotypic accuracy and *G**V* is the genotypic variance.

> printing the estimated BLUP for genotypes X environment

``` r
options (digits = 4)
print(WAASB$BLUPgge[1:10,])
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
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
CF2010
</td>
<td style="text-align:left;">
G1
</td>
<td style="text-align:right;">
-0.0659
</td>
<td style="text-align:right;">
-0.0536
</td>
<td style="text-align:right;">
-0.1195
</td>
<td style="text-align:right;">
2.401
</td>
<td style="text-align:right;">
2.300
</td>
<td style="text-align:right;">
2.502
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G10
</td>
<td style="text-align:right;">
-0.2594
</td>
<td style="text-align:right;">
-0.1473
</td>
<td style="text-align:right;">
-0.4067
</td>
<td style="text-align:right;">
2.114
</td>
<td style="text-align:right;">
2.013
</td>
<td style="text-align:right;">
2.215
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G2
</td>
<td style="text-align:right;">
0.2415
</td>
<td style="text-align:right;">
0.0095
</td>
<td style="text-align:right;">
0.2510
</td>
<td style="text-align:right;">
2.772
</td>
<td style="text-align:right;">
2.671
</td>
<td style="text-align:right;">
2.873
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G3
</td>
<td style="text-align:right;">
0.1092
</td>
<td style="text-align:right;">
0.1994
</td>
<td style="text-align:right;">
0.3087
</td>
<td style="text-align:right;">
2.829
</td>
<td style="text-align:right;">
2.728
</td>
<td style="text-align:right;">
2.930
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G4
</td>
<td style="text-align:right;">
0.0410
</td>
<td style="text-align:right;">
0.0049
</td>
<td style="text-align:right;">
0.0459
</td>
<td style="text-align:right;">
2.567
</td>
<td style="text-align:right;">
2.466
</td>
<td style="text-align:right;">
2.667
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G5
</td>
<td style="text-align:right;">
-0.1508
</td>
<td style="text-align:right;">
-0.0999
</td>
<td style="text-align:right;">
-0.2507
</td>
<td style="text-align:right;">
2.270
</td>
<td style="text-align:right;">
2.169
</td>
<td style="text-align:right;">
2.371
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G6
</td>
<td style="text-align:right;">
-0.0696
</td>
<td style="text-align:right;">
-0.1128
</td>
<td style="text-align:right;">
-0.1825
</td>
<td style="text-align:right;">
2.338
</td>
<td style="text-align:right;">
2.237
</td>
<td style="text-align:right;">
2.439
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G7
</td>
<td style="text-align:right;">
0.1543
</td>
<td style="text-align:right;">
0.0153
</td>
<td style="text-align:right;">
0.1697
</td>
<td style="text-align:right;">
2.690
</td>
<td style="text-align:right;">
2.590
</td>
<td style="text-align:right;">
2.791
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G8
</td>
<td style="text-align:right;">
0.0655
</td>
<td style="text-align:right;">
0.2778
</td>
<td style="text-align:right;">
0.3433
</td>
<td style="text-align:right;">
2.864
</td>
<td style="text-align:right;">
2.763
</td>
<td style="text-align:right;">
2.965
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G9
</td>
<td style="text-align:right;">
-0.0658
</td>
<td style="text-align:right;">
-0.0934
</td>
<td style="text-align:right;">
-0.1592
</td>
<td style="text-align:right;">
2.361
</td>
<td style="text-align:right;">
2.261
</td>
<td style="text-align:right;">
2.462
</td>
</tr>
</tbody>
</table>
This output shows the predicted means for each genotype and environment combination **BLUPg** is the genotypic effect described above. **BLUPge** is the genotypic effect of the *i*th genotype in the *j*th environment $(\\hat{g}\_{ij})$ estimated by $\\hat{g}\_{ij} = h\_g^2(\\bar{y}\_{i.}-\\bar{y}\_{..})+h\_{ge}^2(y\_{ij}-\\bar{y}\_{i.}-\\bar{y}\_{.j}+\\bar{y}\_{..})$, where *h*<sub>*g**e*</sub><sup>2</sup> is the shrinkage effect for genotype-by-environment interaction; **BLUPg+ge** is *B**L**U**P*<sub>*g*</sub> + *B**L**U**P*<sub>*g**e*</sub>; **Predicted** is the predicted mean ($\\hat{y}\_{ij}$) estimated by $\\hat{y}\_{ij} = \\bar{y}\_{.j}+BLUP\_{g+ge}$.

> printing the eigenvalues

``` r
options (digits = 4)

print(WAASB$PCA)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
PC
</th>
<th style="text-align:right;">
Eigenvalue
</th>
<th style="text-align:right;">
Proportion
</th>
<th style="text-align:right;">
Accumulated
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
1
</td>
<td style="text-align:right;">
1.8409
</td>
<td style="text-align:right;">
32.4894
</td>
<td style="text-align:right;">
32.49
</td>
</tr>
<tr>
<td style="text-align:left;">
2
</td>
<td style="text-align:right;">
1.5178
</td>
<td style="text-align:right;">
26.7875
</td>
<td style="text-align:right;">
59.28
</td>
</tr>
<tr>
<td style="text-align:left;">
3
</td>
<td style="text-align:right;">
0.9574
</td>
<td style="text-align:right;">
16.8973
</td>
<td style="text-align:right;">
76.17
</td>
</tr>
<tr>
<td style="text-align:left;">
4
</td>
<td style="text-align:right;">
0.5214
</td>
<td style="text-align:right;">
9.2020
</td>
<td style="text-align:right;">
85.38
</td>
</tr>
<tr>
<td style="text-align:left;">
5
</td>
<td style="text-align:right;">
0.3839
</td>
<td style="text-align:right;">
6.7746
</td>
<td style="text-align:right;">
92.15
</td>
</tr>
<tr>
<td style="text-align:left;">
6
</td>
<td style="text-align:right;">
0.2218
</td>
<td style="text-align:right;">
3.9149
</td>
<td style="text-align:right;">
96.07
</td>
</tr>
<tr>
<td style="text-align:left;">
7
</td>
<td style="text-align:right;">
0.1060
</td>
<td style="text-align:right;">
1.8712
</td>
<td style="text-align:right;">
97.94
</td>
</tr>
<tr>
<td style="text-align:left;">
8
</td>
<td style="text-align:right;">
0.0843
</td>
<td style="text-align:right;">
1.4878
</td>
<td style="text-align:right;">
99.42
</td>
</tr>
<tr>
<td style="text-align:left;">
9
</td>
<td style="text-align:right;">
0.0326
</td>
<td style="text-align:right;">
0.5752
</td>
<td style="text-align:right;">
100.00
</td>
</tr>
</tbody>
</table>
> Plotting the eigenvalues

``` r
plot.eigen(WAASB, size.lab = 14, size.tex = 14)
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-39-1.png" style="display: block; margin: auto;" />

The above output shows the eigenvalues and the proportion of variance explained by each principal component axis of the BLUP interaction effects matrix.

> printing the phenotypic means for all genotype x environment combinations

``` r
options (digits = 4)
print(WAASB$MeansGxE[1:10,])
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
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
envPC1
</th>
<th style="text-align:right;">
genPC1
</th>
<th style="text-align:right;">
nominal
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G1
</td>
<td style="text-align:right;">
2.366
</td>
<td style="text-align:right;">
0.2326
</td>
<td style="text-align:right;">
0.1097
</td>
<td style="text-align:right;">
2.649
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G10
</td>
<td style="text-align:right;">
1.974
</td>
<td style="text-align:right;">
0.2326
</td>
<td style="text-align:right;">
-0.7340
</td>
<td style="text-align:right;">
2.335
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G2
</td>
<td style="text-align:right;">
2.902
</td>
<td style="text-align:right;">
0.2326
</td>
<td style="text-align:right;">
0.0630
</td>
<td style="text-align:right;">
2.717
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G3
</td>
<td style="text-align:right;">
2.888
</td>
<td style="text-align:right;">
0.2326
</td>
<td style="text-align:right;">
-0.0830
</td>
<td style="text-align:right;">
2.922
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G4
</td>
<td style="text-align:right;">
2.589
</td>
<td style="text-align:right;">
0.2326
</td>
<td style="text-align:right;">
0.2908
</td>
<td style="text-align:right;">
2.765
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G5
</td>
<td style="text-align:right;">
2.189
</td>
<td style="text-align:right;">
0.2326
</td>
<td style="text-align:right;">
0.0776
</td>
<td style="text-align:right;">
2.584
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G6
</td>
<td style="text-align:right;">
2.301
</td>
<td style="text-align:right;">
0.2326
</td>
<td style="text-align:right;">
0.1195
</td>
<td style="text-align:right;">
2.577
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G7
</td>
<td style="text-align:right;">
2.774
</td>
<td style="text-align:right;">
0.2326
</td>
<td style="text-align:right;">
0.6690
</td>
<td style="text-align:right;">
2.866
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G8
</td>
<td style="text-align:right;">
2.899
</td>
<td style="text-align:right;">
0.2326
</td>
<td style="text-align:right;">
-0.0203
</td>
<td style="text-align:right;">
3.034
</td>
</tr>
<tr>
<td style="text-align:left;">
CF2010
</td>
<td style="text-align:left;">
G9
</td>
<td style="text-align:right;">
2.326
</td>
<td style="text-align:right;">
0.2326
</td>
<td style="text-align:right;">
-0.4923
</td>
<td style="text-align:right;">
2.459
</td>
</tr>
</tbody>
</table>
Is the phenotypic mean for each genotype and environment combination (*y*<sub>*i**j*</sub>), estimated by *y*<sub>*i**j*</sub> = ∑<sub>*k*</sub>*y*<sub>*i**j*</sub>/*B* with *k* = 1, 2, ...*B*.

plotting the scores
-------------------

For both `WAAS.AMMI` and `WAASB` functions three types of graphics can be generated: 1 = PC1 x PC2, to make inferences related to the interaction effects; 2 = GY x PC1 to make inferences related to stability and productivity; 3 = GY x WAASB and 4 = nominal yield x environment PCA1 scores. These are ggplot2-based graphics; so, the apprareance of the graphics can be be personalized by using the function `theme()`

### biplot type 1: PC1 x PC2

``` r
# No file exported
plot.scores(WAASB,
            type = 1)
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-42-1.png" style="display: block; margin: auto;" />

### biplot type 2: GY x PC1

``` r
plot.scores(WAASB,
            type = 2)
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-43-1.png" style="display: block; margin: auto;" />

``` r
# dafault theme of ggplot2
library(ggplot2)
plot.scores(WAASB,
            type = 2,
            col.gen = "black",
            col.env = "blue",
            theme = theme_gray())
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-43-2.png" style="display: block; margin: auto;" />

### biplot type 3: GY x WAASB

The quadrants proposed in the following biplot represent the four classifications proposed here regarding the joint interpretation of productivity and stability. The genotypes or environments included in the quadrant I can be considered unstable genotypes or environments with high discrimination ability, however with productivity below the grand mean. In quadrant II are included unstable genotypes, however, with productivity above the grand mean. The environments included in this quadrant deserve special attention since, in addition to providing high magnitudes of the response variable, they present a good discrimination ability. Genotypes within the quadrant III have low productivity, but can be considered stable due to lower values of WAASB. The lower this value, the more stable the genotype can be considered. The environments included in this quadrant can be considered as poorly productive and with low discrimination ability. The genotypes within the quadrant IV the genotypes can be considered an "ideal" genotype due to the high magnitude of the response variable and high stability performance (lower values of WAASB).

``` r
plot.scores(WAASB,
            type = 3)
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-44-1.png" style="display: block; margin: auto;" />

``` r
# Save to a *.tiff file with resolution of 600 dpi.
plot.scores(WAASB,
            type = 3,
            export = TRUE,
            file.type = "tiff",
            resolution = 600)
```

### Nominal yield vs environment IPCA1 scores

Aiming at identifying possible mega-environments as well as visualizing the "which-won-where" pattern of the dataset, a graphic with the nominal yield ($\\mathop {\\hat y}\\nolimits\_{ij}^\*$) as a function of the environment IPCA1 scores is also produced. In this graphic, each genotype is depicted by a straight line with the equation $\\mathop {\\hat y}\\nolimits\_{ij}^\* = {\\mu \_i} + PCA{1\_i} \\times PCA{1\_j}$, where $\\mathop {\\hat y}\\nolimits\_{ij}^\*$ is the nominal yield for the ith genotype in the jth environment; *μ*<sub>*i*</sub>is the grand mean of the ith genotype; *P**C*1<sub>*i*</sub> is the IPCA1 score of the ith genotype and *P**C**A*1<sub>*j*</sub> is the IPC1 score of the jth environment. The winner genotype in a given environment has the highest nominal yield in that environment. This is the type 4 graphic and can be obtained by the following code.

``` r
plot.scores(WAASB,
            type = 4)
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-45-1.png" style="display: block; margin: auto;" />

Estimating the WAASBY index
===========================

The function `WAASBratio` considers both stability (weighted average of absolute scores based on SVD of BLUP-interaction effects matrix) and productivity for genotype ranking considering the following model:

$$
            \\begin{aligned}
              WAASBY\_I = 
            \\frac{\\left\[W\_{GY}\\times\\left(\\left(\\frac{GY\_i}{GY\_{max}}\\right)\\times 100\\right)\\right\] + 
              \\left\[W\_{S}\\times\\left(100-\\frac{WAASB\_i}{WAASB\_{min}}\\right)\\right\]} {W\_{GY}+W\_{S}}
           \\end{aligned}
$$

where *W**A**A**S**B**Y*<sub>*i*</sub> is the weighted average of absolute scores and productivity of the *i*th genotype; *W*<sub>*G**Y*</sub> is the weight for grain yield; *G**Y*<sub>*i*</sub> is the average grain yield of the *i*th genotype; *G**Y*<sub>*m**a**x*</sub> is the largest average grain yield among the genotypes; *W*<sub>*S*</sub> is the weight for the stability; *W**A**A**S**B*<sub>*i*</sub> is the weighted average of absolute scores of the *i*th genotype; and *W**A**A**S**B*<sub>*m**i**n*</sub> is the smallest value of WAASB.

This function provide the option of attributing weights for stability and productivity in genotype ranking. This is important depending on the goal of a selection strategy. For example, if a goal of a breeding program is to select a genotype whit high yielding (independently on the stability performance), that genotype with the first rank in an WAASB/GY = 0/100 ratio should be selected. The reciprocal is true. Aiming at selecting a high-stable genotype (independently on the productivity), that genotype with the first rank in an WAASB/GY = 100/0 ratio should be selected. By default, the increment on the WAASB/GY ratio is equal to 5. In other words, twenty one different combinations are computed. Each combination, the genotypes are ranked regarding the WAASY value.

``` r
WAASBYratio = WAASBYratio(dataset,
                          resp = GY,
                          gen = GEN,
                          env = ENV,
                          rep = BLOCK,
                          increment = 10,
                          saveWAASY = 50)
```

Printing the model outputs
--------------------------

> printing the WAASY values

``` r
print(WAASBYratio$WAASY)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
Code
</th>
<th style="text-align:right;">
PesRes
</th>
<th style="text-align:right;">
PesWAAS
</th>
<th style="text-align:right;">
WAASY
</th>
<th style="text-align:left;">
Mean
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
G10
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
89.40
</td>
<td style="text-align:left;">
below
</td>
</tr>
<tr>
<td style="text-align:left;">
G9
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
90.74
</td>
<td style="text-align:left;">
below
</td>
</tr>
<tr>
<td style="text-align:left;">
G6
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
91.35
</td>
<td style="text-align:left;">
below
</td>
</tr>
<tr>
<td style="text-align:left;">
G5
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
91.52
</td>
<td style="text-align:left;">
below
</td>
</tr>
<tr>
<td style="text-align:left;">
G1
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
92.57
</td>
<td style="text-align:left;">
below
</td>
</tr>
<tr>
<td style="text-align:left;">
G4
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
93.21
</td>
<td style="text-align:left;">
below
</td>
</tr>
<tr>
<td style="text-align:left;">
G7
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
93.40
</td>
<td style="text-align:left;">
above
</td>
</tr>
<tr>
<td style="text-align:left;">
G2
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
93.68
</td>
<td style="text-align:left;">
above
</td>
</tr>
<tr>
<td style="text-align:left;">
G3
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
97.88
</td>
<td style="text-align:left;">
above
</td>
</tr>
<tr>
<td style="text-align:left;">
G8
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
99.23
</td>
<td style="text-align:left;">
above
</td>
</tr>
</tbody>
</table>
> Printing the genotype ranking for each scenario

``` r
print(WAASBYratio$hetcomb)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
100/0
</th>
<th style="text-align:right;">
90/10
</th>
<th style="text-align:right;">
80/20
</th>
<th style="text-align:right;">
70/30
</th>
<th style="text-align:right;">
60/40
</th>
<th style="text-align:right;">
50/50
</th>
<th style="text-align:right;">
40/60
</th>
<th style="text-align:right;">
30/70
</th>
<th style="text-align:right;">
20/80
</th>
<th style="text-align:right;">
10/90
</th>
<th style="text-align:right;">
0/100
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
G1
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
5
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
6
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
</tr>
<tr>
<td style="text-align:left;">
G10
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
</tr>
<tr>
<td style="text-align:left;">
G2
</td>
<td style="text-align:right;">
6
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
3
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
4
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
G3
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
2
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
2
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
G4
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
7
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
5
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
</tr>
<tr>
<td style="text-align:left;">
G5
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
6
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
8
</td>
</tr>
<tr>
<td style="text-align:left;">
G6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
5
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
8
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
9
</td>
</tr>
<tr>
<td style="text-align:left;">
G7
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
4
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
4
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
</tr>
<tr>
<td style="text-align:left;">
G8
</td>
<td style="text-align:right;">
5
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
1
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
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
G9
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
8
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
7
</td>
</tr>
</tbody>
</table>
> Printing the genotype ranking depending on the number of multiplicative terms used to estimate the WAASB index.

``` r
print(WAASBYratio$hetdata)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
9PCA
</th>
<th style="text-align:right;">
8PCA
</th>
<th style="text-align:right;">
7PCA
</th>
<th style="text-align:right;">
6PCA
</th>
<th style="text-align:right;">
5PCA
</th>
<th style="text-align:right;">
4PCA
</th>
<th style="text-align:right;">
3PCA
</th>
<th style="text-align:right;">
2PCA
</th>
<th style="text-align:right;">
1PCA
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
G1
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
3
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
2
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
</tr>
<tr>
<td style="text-align:left;">
G10
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
</tr>
<tr>
<td style="text-align:left;">
G2
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
6
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
1
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
G3
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
1
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
4
</td>
</tr>
<tr>
<td style="text-align:left;">
G4
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
7
</td>
</tr>
<tr>
<td style="text-align:left;">
G5
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
4
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
5
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
G6
</td>
<td style="text-align:right;">
2
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
3
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
6
</td>
</tr>
<tr>
<td style="text-align:left;">
G7
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
9
</td>
</tr>
<tr>
<td style="text-align:left;">
G8
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
5
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
1
</td>
</tr>
<tr>
<td style="text-align:left;">
G9
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
8
</td>
</tr>
</tbody>
</table>
plotting the WAASBY values
--------------------------

``` r
plot.WAASBY(WAASBYratio,
            theme = theme_waasb() +
                    theme(aspect.ratio = NULL,
                          legend.position = c(0.85, 0.2)))
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-53-1.png" style="display: block; margin: auto;" />

Plotting the heat map graphics
------------------------------

The first type of heatmap shows the genotype ranking depending on the number of principal component axis used for estimating the WAASB index. A euclidean distance-based dendrogram is used for grouping the genotype ranking for both genotypes and principal component axis. The second type of heatmap shows the genotype ranking depending on the WAASB/GY ratio. The ranks obtained with a ratio of 100/0 considers exclusively the stability for genotype ranking. On the other hand, a ratio of 0/100 considers exclusively the productivity for genotype ranking. Four clusters are estimated (1) unproductive and unstable genotypes; (2) productive, but unstable genotypes; (3) stable, but unproductive genotypes; and (4), productive and stable genotypes.

### Ranks of genotypes depending on the number PCAs used to estimate the WAASB

``` r
plot.WAASBYratio(WAASBYratio,
                 type = 1)
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-54-1.png" style="display: block; margin: auto;" />

``` r
# save to a *.pdf file (default)
plot.WAASBYratio(WAASBYratio,
                 type = 1,
                 export = T)
```

### Ranks of genotypes depending on the WAASB/GY ratio

``` r
plot.WAASBYratio(WAASBYratio,
                 type = 2)
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-55-1.png" style="display: block; margin: auto;" />

``` r
#save to a *tiff file
plot.WAASBYratio(WAASBYratio,
                type = 2,
               export = T,
              file.type = "tiff",
             resolution = 600)
```

    ## RStudioGD 
    ##         2

References
==========

Gauch, H.G, and Zobel R.W. 1988. “Predictive and Postdictive Success of Statistical Analyses of Yield Trials.” *Theor. Appl. Genet.* 76 (1): 1–10. doi:[10.1007/BF00288824](https://doi.org/10.1007/BF00288824).

Olivoto, T. 2018. “WAASBdata.” *Mendeley Data* V1. doi:[10.17632/2sjz32k3s3.1](https://doi.org/10.17632/2sjz32k3s3.1).

Piepho, H.P. 1994. “Best Linear Unbiased Prediction (Blup) for Regional Yield Trials: A Comparison to Additive Main Effects and Multiplicative Interaction (Ammi) Analysis.” *Theor. Appl. Genet.* 89 (5): 647–54. doi:[10.1007/BF00222462](https://doi.org/10.1007/BF00222462).
