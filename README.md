Extending WAASB package
================
Tiago Olivoto
WAASB 1.0.1 2018-05-09

<style>
.column-left{
  float: left;
  width: 49%;
  text-align: left;
}
.column-right{
  float: right;
  width: 49%;
  text-align: left;
}

</style>
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

The WAASB R package was developed in R language and is distributed under the GPL 3.0 licence. Our main motivation for the package development was to provide friendly, easy and reliable functions to reproduce the methods proposed in the paper, as well as to share codes used for traditional AMMI analysis and BLUP prediction. The package can be downloaded by clicking [here](https://github.com/TiagoOlivoto/WAASB/archive/master.zip) or running the following code in the R console:

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
<td align="left"><code>WAASB()</code></td>
<td align="left">Weighted Average of Absolute Scores for BLUP-based PCA analysis</td>
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
<td align="left"><code>validation.blup()</code></td>
<td align="left">Cross-validation for blup prediction</td>
</tr>
<tr class="even">
<td align="left"><code>validation.AMMIF()</code></td>
<td align="left">BLUP-interaction effects matrix in different scenarios of WAASB/GY ratio</td>
</tr>
<tr class="odd">
<td align="left"><code>predict.AMMI()</code></td>
<td align="left">Predict means in replicated-based MET using AMMI analysis</td>
</tr>
<tr class="even">
<td align="left"><code>predict.AMMImean()</code></td>
<td align="left">Predict means in AMMI analysis with data from single observation</td>
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
<td align="left"><code>summary.WAASB()</code></td>
<td align="left">Summary a WAASB object</td>
</tr>
<tr class="even">
<td align="left"><code>plot.eigen()</code></td>
<td align="left">Plot the eigenvalues of singular value decomposition</td>
</tr>
</tbody>
</table>

Get started
===========

The data set from experiment 1 (oat) was provided to make reproducible examples. This data is available in Olivoto ([2018](#ref-Olivoto:2018)). An own data set can be used provided that the columns are in the following order: environment, genotype, block/replicate and response variable(s).

``` r
# Importing data
require(WAASB)
dataset = read.csv("https://data.mendeley.com/datasets/2sjz32k3s3/1/files/07764a07-172a-4285-85db-c31bc39ae480/WAASBdata.csv?dl=1")
```

Predictive accuracy
===================

The functions `validation.AMMIF` provides a complete cross-validation for AMMI model family using replicate-based data. Automatically the first validation is carried out considering the AMMIF (all possible axis used). Considering this model, the original data set is split into two data sets: modelling and validating data. The data set "modelling" has all combinations (genotype x environment) with the number of replications informed in `nrepval` argument. The data set "validating" has one replication. The splitting of the data set into modelling and validating data depends on the design informed. For Completely Randomized Block Design (default), completely blocks are selected within environments, as suggested by Piepho ([1994](#ref-Piepho:1994)). The remained block serves validation data. If `design = "RCD"` is informed, thus declaring that the a completely randomized design was used, single observations are randomized for each treatment (genotype-by-environment combination). This is the same procedure suggested by Gauch and Zobel ([1988](#ref-Gauch:1988)). The estimated values (for each AMMI model family) are compared with the "validating" data. the Root Means Square error is computed according to the following equation:

$$
   RMSE = 
   \\sqrt{\\frac{\\sum\_{i = 1}^n(\\widehat{y}\_{ij}-y\_{ij})^2} {n}}
$$

where $\\widehat{y}\_{ij}$ is the model predicted value; and *y*<sub>*i**j*</sub> is the observed value of the validating data set. At the end of boots, a list with all estimated RMSE and the average of RMSE is returned.

The function `validation.blup` provides a cross-validation of replicate-based data using mixed models. By default, complete blocks are randomly selected each environment. The procedure for computing the RSME is identical to the above function.

> The following codeS computes computes the cross-validation of oat data set based on 1000 re-sampling procedures. This number can be changed.

``` r
# cross-validation for AMMI model family
AMMIweat = validation.AMMIF(dataset,
                            resp = "GY",
                            nboot = 10,
                            nrepval = 2)

# cross-validation for BLUP model
BLUPweat = validation.blup(dataset,
                            resp = "GY",
                            nboot = 10,
                            nrepval = 2)
```

A progress bar is shown by default (the examples are bellow). Thus, it is possible to verify the status of the cross-validation process. If necessary, the progress bar can be disabled by informing the argument `progbar = FALSE` in the function.

![Progress bar for AMMI model family cross-validation process.](printAMMI.png)

![Progress bar for BLUP cross-validation process.](printBLUP.png)

Printting the means of RMSE estimates
-------------------------------------

``` r
options(digits = 4)
RMSEweat = rbind(AMMIweat$RMSEmean, BLUPweat$RMSEmean)
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
0.4021
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
AMMI3
</td>
<td style="text-align:right;">
0.4156
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
AMMI4
</td>
<td style="text-align:right;">
0.4200
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
AMMI2
</td>
<td style="text-align:right;">
0.4234
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
AMMI5
</td>
<td style="text-align:right;">
0.4274
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
AMMI0
</td>
<td style="text-align:right;">
0.4284
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
AMMI6
</td>
<td style="text-align:right;">
0.4302
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
0.4309
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
0.4309
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
AMMI7
</td>
<td style="text-align:right;">
0.4314
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
0.4353
</td>
<td style="text-align:left;">
Wheat
</td>
</tr>
</tbody>
</table>
The results shown above are the mean of the 50 RMSE estimates for each tested model and are presented from the most accurate model (smallest RSME value) to the lowest accurate model (highest RMSE value).

Plotting the RMSE values
------------------------

``` r
# binding AMMI and BLUP RMSEs
RMSEweat = list(RMSE = rbind(AMMIweat$RMSE, 
                           BLUPweat$RMSE))
# Plotting the RMSE values
plot.validation.AMMIF(RMSEweat,
                      violin = FALSE,
                      col.boxplot = "gray75")
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

Five statistics are shown in this boxplot. The median, the lower and upper hinges correspond to the first and third quartiles (the 25th and 75th percentiles, respectively). The upper whisker extends from the hinge to the largest value no further than 1.5 × *I**Q**R* from the hinge (where IQR is the inter-quartile range). The lower whisker extends from the hinge to the smallest value at most 1.5 × *I**Q**R* of the hinge. Data beyond the end of the whiskers are considered outlying points. If the condition `violin = TRUE`, a violin plot is added joint with the boxplot. A violin plot is a compact display of a continuous distribution displayed in the same way as a boxplot.

Predicting the yield using traditional AMMI model
=================================================

The function `predict.AMMI` can be used to predict the response variable of a two-way table (for examples the yielding of the *i*th genotype in the *j*th environment) based on the traditional AMMI model. This prediction is based on the number of multiplicative terms used. If `naxis = 0`, only the main effects (AMMI0) are used. In this case, the predicted mean will be the predicted value from OLS estimation. If `naxis = 1` the AMMI1 (with one multiplicative term) is used for predicting the response variable. If `naxis = min(gen-1;env-1)`, the AMMIF is fitted and the predicted value will be the cell mean, i.e., the mean of R replicates of the *i*th genotype in the *j*th environment. The number of axes to be used must be carefully chosen. Procedures based on postdictive success (such as Gollobs's d.f.) or predictive success (such as cross-validation) should be used to do this. This package provides both. `WAAS.AMMI` function compute traditional AMMI analysis showing the number of significant axes. On the other hand, `validation.AMMIF` function provides a cross-validation, estimating the RMSE of all AMMI model family based on re-sampling procedures.

> Predicting the yield of 10 oat cultivars in 16 environments using 5 multiplicative terms.

``` r
predictoat = predict.AMMI(dataset,
                          resp = "GY",
                          naxis = 5)
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

<br />

Estimating the WAAS (based on traditional AMMI model)
=====================================================

The `WAAS.AMMI` function compute the weighted average of absolute scores considering (i) all principal component axes that are significant at a given probability error level; or (ii) declaring a specific number of axes to be used, according to the following equation:

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
                 resp = "GY",
                 p.valuePC = 0.05,
                 weight.response = 50,
                 weight.WAAS = 50)

# Priorizing productivity for genotype ranking (WAASB/GY ratio = 30/70)
# no output for this script
WAAS11 = WAAS.AMMI(dataset,
                 resp = "GY",
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
0.00000000000000000000000
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
0.00000031778918565985957
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
0.00000000000000000083806
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
0.00000000000005639378097
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
0.00000000000000000000000
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
0.00000000000000000000000
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
0.00000000000000000000000
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
0.2167
</td>
<td style="text-align:right;">
2.100
</td>
<td style="text-align:left;">
0.00719999999999999980294
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
0.04660000000000000253131
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
0.37609999999999998987477
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
0.79390000000000005009326
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
0.4653
</td>
<td style="text-align:right;">
0.0517
</td>
<td style="text-align:right;">
0.500
</td>
<td style="text-align:left;">
0.87399999999999999911182
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
0.96150000000000002131628
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
0.1878
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
0.5979
</td>
<td style="text-align:right;">
82.46
</td>
<td style="text-align:right;">
96.43
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
89.45
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
0.2820
</td>
<td style="text-align:right;">
88.93
</td>
<td style="text-align:right;">
98.32
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
93.62
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
0.1676
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
0.3880
</td>
<td style="text-align:right;">
88.74
</td>
<td style="text-align:right;">
97.69
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
93.21
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
0.2213
</td>
<td style="text-align:right;">
84.41
</td>
<td style="text-align:right;">
98.68
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
91.55
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
0.2070
</td>
<td style="text-align:right;">
83.88
</td>
<td style="text-align:right;">
98.76
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
91.32
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
0.4005
</td>
<td style="text-align:right;">
89.17
</td>
<td style="text-align:right;">
97.61
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
93.39
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
0.2488
</td>
<td style="text-align:right;">
100.00
</td>
<td style="text-align:right;">
98.52
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
99.26
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
0.5208
</td>
<td style="text-align:right;">
84.68
</td>
<td style="text-align:right;">
96.89
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
90.79
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
0.2078
</td>
<td style="text-align:right;">
62.02
</td>
<td style="text-align:right;">
98.07
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
80.04
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
0.2021
</td>
<td style="text-align:right;">
78.24
</td>
<td style="text-align:right;">
98.12
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
88.18
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
0.3268
</td>
<td style="text-align:right;">
100.00
</td>
<td style="text-align:right;">
96.96
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
98.48
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
0.3300
</td>
<td style="text-align:right;">
90.43
</td>
<td style="text-align:right;">
96.93
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
93.68
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
0.2925
</td>
<td style="text-align:right;">
61.69
</td>
<td style="text-align:right;">
97.28
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
79.48
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
0.3169
</td>
<td style="text-align:right;">
76.45
</td>
<td style="text-align:right;">
97.05
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
86.75
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
0.2370
</td>
<td style="text-align:right;">
96.22
</td>
<td style="text-align:right;">
97.80
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
97.01
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
0.1075
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
0.2560
</td>
<td style="text-align:right;">
48.94
</td>
<td style="text-align:right;">
97.62
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
73.28
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
0.2173
</td>
<td style="text-align:right;">
62.41
</td>
<td style="text-align:right;">
97.98
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
80.19
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
0.4817
</td>
<td style="text-align:right;">
75.21
</td>
<td style="text-align:right;">
95.52
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
85.36
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
0.1992
</td>
<td style="text-align:right;">
53.52
</td>
<td style="text-align:right;">
98.15
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
75.83
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
0.1116
</td>
<td style="text-align:right;">
33.66
</td>
<td style="text-align:right;">
98.96
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
66.31
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
0.2011
</td>
<td style="text-align:right;">
39.58
</td>
<td style="text-align:right;">
98.13
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
68.85
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
0.3728
</td>
<td style="text-align:right;">
71.59
</td>
<td style="text-align:right;">
96.53
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
84.06
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
0.1529
</td>
<td style="text-align:right;">
43.84
</td>
<td style="text-align:right;">
98.58
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
71.21
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

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

The biplot above shows the coordinates of the genotypes and environments regarding the grain yield and the WAAS values. It is in fact the plot of the columns Y and WAAS from the above table. A detailed discussion on the interpretation of this biplot can be seen in [section 6.2.3](#biplot-type-3-gy-x-waasb).

Declaring a specific number of axis to be used
----------------------------------------------

Here it is possible to inform a specific number of multiplicative terms to be used in the WAAS estimation. In this example, the number of terms informed is used independently of its significance. Let us, for the moment, assume that two multiplicative terms should be used.

> Declaring that five PCA must be used to compute the WAAS

``` r
WAAS2 = WAAS.AMMI(dataset,
                 resp = "GY",
                 naxis = 7,
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
0.1939
</td>
<td style="text-align:right;">
86.32
</td>
<td style="text-align:right;">
98.84
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
92.58
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
0.5694
</td>
<td style="text-align:right;">
82.46
</td>
<td style="text-align:right;">
96.60
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
89.53
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
0.2696
</td>
<td style="text-align:right;">
88.93
</td>
<td style="text-align:right;">
98.39
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
93.66
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
0.1675
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
0.3767
</td>
<td style="text-align:right;">
88.74
</td>
<td style="text-align:right;">
97.75
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
93.25
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
0.2299
</td>
<td style="text-align:right;">
84.41
</td>
<td style="text-align:right;">
98.63
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
91.52
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
0.2002
</td>
<td style="text-align:right;">
83.88
</td>
<td style="text-align:right;">
98.81
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
91.34
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
0.3870
</td>
<td style="text-align:right;">
89.17
</td>
<td style="text-align:right;">
97.69
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
93.43
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
0.2556
</td>
<td style="text-align:right;">
100.00
</td>
<td style="text-align:right;">
98.47
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
99.24
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
0.5070
</td>
<td style="text-align:right;">
84.68
</td>
<td style="text-align:right;">
96.97
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
90.83
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
0.1993
</td>
<td style="text-align:right;">
62.02
</td>
<td style="text-align:right;">
98.15
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
80.08
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
0.1952
</td>
<td style="text-align:right;">
78.24
</td>
<td style="text-align:right;">
98.18
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
88.21
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
0.3200
</td>
<td style="text-align:right;">
100.00
</td>
<td style="text-align:right;">
97.02
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
98.51
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
0.3180
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
13
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
0.2833
</td>
<td style="text-align:right;">
61.69
</td>
<td style="text-align:right;">
97.36
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
79.52
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
0.3012
</td>
<td style="text-align:right;">
76.45
</td>
<td style="text-align:right;">
97.20
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
86.82
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
0.2435
</td>
<td style="text-align:right;">
96.22
</td>
<td style="text-align:right;">
97.73
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
96.97
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
0.1096
</td>
<td style="text-align:right;">
65.53
</td>
<td style="text-align:right;">
98.98
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
82.25
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
0.2507
</td>
<td style="text-align:right;">
48.94
</td>
<td style="text-align:right;">
97.67
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
73.30
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
0.2204
</td>
<td style="text-align:right;">
62.41
</td>
<td style="text-align:right;">
97.95
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
80.18
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
0.4616
</td>
<td style="text-align:right;">
75.21
</td>
<td style="text-align:right;">
95.70
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
85.45
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
0.1913
</td>
<td style="text-align:right;">
53.52
</td>
<td style="text-align:right;">
98.22
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
75.87
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
0.1075
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
0.2103
</td>
<td style="text-align:right;">
39.58
</td>
<td style="text-align:right;">
98.04
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
68.81
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
0.3590
</td>
<td style="text-align:right;">
71.59
</td>
<td style="text-align:right;">
96.66
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
84.12
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
0.1537
</td>
<td style="text-align:right;">
43.84
</td>
<td style="text-align:right;">
98.57
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
71.21
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

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

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
              resp = "GY",
              prob = 0.95,
              weight.response = 50,
              weight.WAAS = 50)

# Priorizing productivity for genotype ranking (WAASB/GY ratio = 30/70)
# no output for this script
WAASB2 = WAASB(dataset,
              resp = "GY",
              prob = 0.95,
              weight.response = 30,
              weight.WAAS = 70)

# Priorizing productivity for genotype ranking (WAASB/GY ratio = 30/70)
# no output for this script
WAASB3 = WAASB(dataset,
              resp = "GY",
              random = "all",
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
<th style="border-bottom:hidden" colspan="1">
</th>
<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4">
Genotype random effect

</th>
<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="6">
Genotype and environment random effects

</th>
</tr>
<tr>
<th style="border-bottom:hidden" colspan="1">
</th>
<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">
Genotype LRT

</th>
<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">
Interaction LRT

</th>
<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">
Genotype LRT

</th>
<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">
Environment LRT

</th>
<th style="border-bottom:hidden; padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="2">
Interaction LRT

</th>
</tr>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:center;">
reducedG<sup>\*</sup>
</th>
<th style="text-align:center;">
Complete<sup>†</sup>
</th>
<th style="text-align:center;">
ReducedGE<sup>‡</sup>
</th>
<th style="text-align:center;">
Complete
</th>
<th style="text-align:center;">
reducedG
</th>
<th style="text-align:center;">
Complete
</th>
<th style="text-align:center;">
reducedE<sup>§</sup>
</th>
<th style="text-align:center;">
Complete
</th>
<th style="text-align:center;">
ReducedGE
</th>
<th style="text-align:center;">
Complete
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Df
</td>
<td style="text-align:center;">
50.00
</td>
<td style="text-align:center;">
51.00000
</td>
<td style="text-align:center;">
50.00
</td>
<td style="text-align:center;">
51.000
</td>
<td style="text-align:center;">
51.00
</td>
<td style="text-align:center;">
52.00000
</td>
<td style="text-align:center;">
51.00
</td>
<td style="text-align:center;">
52.00
</td>
<td style="text-align:center;">
51.00
</td>
<td style="text-align:center;">
52.000
</td>
</tr>
<tr>
<td style="text-align:left;">
AIC
</td>
<td style="text-align:center;">
524.48
</td>
<td style="text-align:center;">
507.25568
</td>
<td style="text-align:center;">
566.37
</td>
<td style="text-align:center;">
507.256
</td>
<td style="text-align:center;">
526.48
</td>
<td style="text-align:center;">
509.25568
</td>
<td style="text-align:center;">
507.26
</td>
<td style="text-align:center;">
509.26
</td>
<td style="text-align:center;">
568.37
</td>
<td style="text-align:center;">
509.256
</td>
</tr>
<tr>
<td style="text-align:left;">
BIC
</td>
<td style="text-align:center;">
733.17
</td>
<td style="text-align:center;">
720.11877
</td>
<td style="text-align:center;">
775.06
</td>
<td style="text-align:center;">
720.119
</td>
<td style="text-align:center;">
739.35
</td>
<td style="text-align:center;">
726.29255
</td>
<td style="text-align:center;">
720.12
</td>
<td style="text-align:center;">
726.29
</td>
<td style="text-align:center;">
781.24
</td>
<td style="text-align:center;">
726.293
</td>
</tr>
<tr>
<td style="text-align:left;">
logLik
</td>
<td style="text-align:center;">
-212.24
</td>
<td style="text-align:center;">
-202.62784
</td>
<td style="text-align:center;">
-233.19
</td>
<td style="text-align:center;">
-202.628
</td>
<td style="text-align:center;">
-212.24
</td>
<td style="text-align:center;">
-202.62784
</td>
<td style="text-align:center;">
-202.63
</td>
<td style="text-align:center;">
-202.63
</td>
<td style="text-align:center;">
-233.19
</td>
<td style="text-align:center;">
-202.628
</td>
</tr>
<tr>
<td style="text-align:left;">
deviance
</td>
<td style="text-align:center;">
424.48
</td>
<td style="text-align:center;">
405.25568
</td>
<td style="text-align:center;">
466.37
</td>
<td style="text-align:center;">
405.256
</td>
<td style="text-align:center;">
424.48
</td>
<td style="text-align:center;">
405.25568
</td>
<td style="text-align:center;">
405.26
</td>
<td style="text-align:center;">
405.26
</td>
<td style="text-align:center;">
466.37
</td>
<td style="text-align:center;">
405.256
</td>
</tr>
<tr>
<td style="text-align:left;">
Chisq
</td>
<td style="text-align:center;">
NA
</td>
<td style="text-align:center;">
19.22808
</td>
<td style="text-align:center;">
NA
</td>
<td style="text-align:center;">
61.117
</td>
<td style="text-align:center;">
NA
</td>
<td style="text-align:center;">
19.22808
</td>
<td style="text-align:center;">
NA
</td>
<td style="text-align:center;">
0.00
</td>
<td style="text-align:center;">
NA
</td>
<td style="text-align:center;">
61.117
</td>
</tr>
<tr>
<td style="text-align:left;">
Chi Df
</td>
<td style="text-align:center;">
NA
</td>
<td style="text-align:center;">
1.00000
</td>
<td style="text-align:center;">
NA
</td>
<td style="text-align:center;">
1.000
</td>
<td style="text-align:center;">
NA
</td>
<td style="text-align:center;">
1.00000
</td>
<td style="text-align:center;">
NA
</td>
<td style="text-align:center;">
1.00
</td>
<td style="text-align:center;">
NA
</td>
<td style="text-align:center;">
1.000
</td>
</tr>
<tr>
<td style="text-align:left;">
Pr(&gt;Chisq)
</td>
<td style="text-align:center;">
NA
</td>
<td style="text-align:center;">
0.00001
</td>
<td style="text-align:center;">
NA
</td>
<td style="text-align:center;">
0.000
</td>
<td style="text-align:center;">
NA
</td>
<td style="text-align:center;">
0.00001
</td>
<td style="text-align:center;">
NA
</td>
<td style="text-align:center;">
1.00
</td>
<td style="text-align:center;">
NA
</td>
<td style="text-align:center;">
0.000
</td>
</tr>
</tbody>
<tfoot>
<tr>
<td style="padding: 0; border: 0;" colspan="100%">
<sup>\*</sup> Reduced model without genotype effect;
</td>
</tr>
<tr>
<td style="padding: 0; border: 0;" colspan="100%">
<sup>†</sup> Complete model;
</td>
</tr>
<tr>
<td style="padding: 0; border: 0;" colspan="100%">
<sup>‡</sup> Reduced model without genotype-vs-environment interaction effect;
</td>
</tr>
<tr>
<td style="padding: 0; border: 0;" colspan="100%">
<sup>§</sup> Reduced model without environment effect;
</td>
</tr>
</tfoot>
</table>
The output `LRT` contains the Likelihood Ratio Tests on linear mixed-effect models depending on the random effects argument in the function `WAASB` (`random = "all"`or `random  = "gen"`). In the case when `random = "gen"` (default) the LRT is computed for genotype and genotype-vs-environment effects. When `random = "all"`, the LRT is computed for genotype, environment and genotype-vs-environment interaction.

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
0.531083326297866
</td>
</tr>
<tr>
<td style="text-align:left;">
Accuracy
</td>
<td style="text-align:left;">
0.728754640669867
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
0.1602
</td>
<td style="text-align:right;">
86.32
</td>
<td style="text-align:right;">
98.82
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
92.57
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
0.4968
</td>
<td style="text-align:right;">
82.46
</td>
<td style="text-align:right;">
96.34
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
89.40
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
0.2125
</td>
<td style="text-align:right;">
88.93
</td>
<td style="text-align:right;">
98.43
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
93.68
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
0.1356
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
0.3151
</td>
<td style="text-align:right;">
88.74
</td>
<td style="text-align:right;">
97.68
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
93.21
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
0.1858
</td>
<td style="text-align:right;">
84.41
</td>
<td style="text-align:right;">
98.63
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
91.52
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
0.1602
</td>
<td style="text-align:right;">
83.88
</td>
<td style="text-align:right;">
98.82
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
91.35
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
0.3218
</td>
<td style="text-align:right;">
89.17
</td>
<td style="text-align:right;">
97.63
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
93.40
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
0.2077
</td>
<td style="text-align:right;">
100.00
</td>
<td style="text-align:right;">
98.47
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
99.23
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
0.4338
</td>
<td style="text-align:right;">
84.68
</td>
<td style="text-align:right;">
96.80
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
90.74
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
0.1870
</td>
<td style="text-align:right;">
62.02
</td>
<td style="text-align:right;">
97.85
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
79.94
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
0.1488
</td>
<td style="text-align:right;">
78.24
</td>
<td style="text-align:right;">
98.29
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
88.27
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
0.2546
</td>
<td style="text-align:right;">
100.00
</td>
<td style="text-align:right;">
97.08
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
98.54
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
0.2903
</td>
<td style="text-align:right;">
90.43
</td>
<td style="text-align:right;">
96.67
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
93.55
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
0.2539
</td>
<td style="text-align:right;">
61.69
</td>
<td style="text-align:right;">
97.09
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
79.39
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
0.2374
</td>
<td style="text-align:right;">
76.45
</td>
<td style="text-align:right;">
97.28
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
86.86
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
0.2201
</td>
<td style="text-align:right;">
96.22
</td>
<td style="text-align:right;">
97.47
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
96.85
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
0.0872
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
0.1967
</td>
<td style="text-align:right;">
48.94
</td>
<td style="text-align:right;">
97.74
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
73.34
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
0.1663
</td>
<td style="text-align:right;">
62.41
</td>
<td style="text-align:right;">
98.09
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
80.25
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
0.3912
</td>
<td style="text-align:right;">
75.21
</td>
<td style="text-align:right;">
95.51
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
85.36
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
0.1835
</td>
<td style="text-align:right;">
53.52
</td>
<td style="text-align:right;">
97.89
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
75.71
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
0.1031
</td>
<td style="text-align:right;">
33.66
</td>
<td style="text-align:right;">
98.82
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
66.24
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
0.1883
</td>
<td style="text-align:right;">
39.58
</td>
<td style="text-align:right;">
97.84
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
68.71
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
0.3256
</td>
<td style="text-align:right;">
71.59
</td>
<td style="text-align:right;">
96.26
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
83.93
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
0.1089
</td>
<td style="text-align:right;">
43.84
</td>
<td style="text-align:right;">
98.75
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
2.808
</td>
<td style="text-align:right;">
3.129
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
2.730
</td>
<td style="text-align:right;">
3.051
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
2.546
</td>
<td style="text-align:right;">
2.867
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
2.540
</td>
<td style="text-align:right;">
2.861
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
2.535
</td>
<td style="text-align:right;">
2.856
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
2.477
</td>
<td style="text-align:right;">
2.798
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
2.437
</td>
<td style="text-align:right;">
2.758
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
2.430
</td>
<td style="text-align:right;">
2.752
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
2.417
</td>
<td style="text-align:right;">
2.739
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
2.383
</td>
<td style="text-align:right;">
2.704
</td>
</tr>
</tbody>
</table>

> plotting the estimated BLUP for genotypes

``` r
# No file exported
plot.blup(WAASB)
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-30-1.png" style="display: block; margin: auto;" />

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
2.241
</td>
<td style="text-align:right;">
2.562
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
1.953
</td>
<td style="text-align:right;">
2.275
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
2.611
</td>
<td style="text-align:right;">
2.932
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
2.669
</td>
<td style="text-align:right;">
2.990
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
2.406
</td>
<td style="text-align:right;">
2.727
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
2.109
</td>
<td style="text-align:right;">
2.431
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
2.178
</td>
<td style="text-align:right;">
2.499
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
2.530
</td>
<td style="text-align:right;">
2.851
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
2.703
</td>
<td style="text-align:right;">
3.025
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
2.201
</td>
<td style="text-align:right;">
2.522
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

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-35-1.png" style="display: block; margin: auto;" />

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
</tr>
</tbody>
</table>
Is the phenotypic mean for each genotype and environment combination (*y*<sub>*i**j*</sub>), estimated by *y*<sub>*i**j*</sub> = ∑<sub>*k*</sub>*y*<sub>*i**j*</sub>/*B* with *k* = 1, 2, ...*B*.

plotting the scores
-------------------

For both `WAAS.AMMI` and `WAASB` functions three types of graphics can be generated: 1 = PC1 x PC2, to make inferences related to the interaction effects; 2 = GY x PC1 to make inferences related to stability and productivity; and 3 = GY x WAASB. These graphics can be plotted as follows:

### biplot type 1: PC1 x PC2

``` r
# No file exported
plot.scores(WAASB,
            type = 1)
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-38-1.png" style="display: block; margin: auto;" />

### biplot type 2: GY x PC1

``` r
plot.scores(WAASB,
            type = 2)
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-39-1.png" style="display: block; margin: auto;" />

``` r
# Save to *.pdf file (default)
plot.scores(WAASB,
            type = 2,
            export = TRUE)
```

### biplot type 3: GY x WAASB

The quadrants proposed in the following biplot represent the four classifications proposed here regarding the joint interpretation of productivity and stability. The genotypes or environments included in the quadrant I can be considered unstable genotypes or environments with high discrimination ability, however with productivity below the grand mean. In quadrant II are included unstable genotypes, however, with productivity above the grand mean. The environments included in this quadrant deserve special attention since, in addition to providing high magnitudes of the response variable, they present a good discrimination ability. Genotypes within the quadrant III have low productivity, but can be considered stable due to lower values of WAASB. The lower this value, the more stable the genotype can be considered. The environments included in this quadrant can be considered as poorly productive and with low discrimination ability. The genotypes within the quadrant IV the genotypes can be considered an "ideal" genotype due to the high magnitude of the response variable and high stability performance (lower values of WAASB).

``` r
plot.scores(WAASB,
            type = 3)
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-40-1.png" style="display: block; margin: auto;" />

``` r
# Save to a *.tiff file with resolution of 600 dpi.
plot.scores(WAASB,
            type = 3,
            export = TRUE,
            file.type = "tiff",
            resolution = 600)
```

    ## png 
    ##   2

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
                          resp = "GY",
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
            legend.pos = c(0.9, 0.2))
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-48-1.png" style="display: block; margin: auto;" />

Plotting the heat map graphics
------------------------------

The first type of heatmap shows the genotype ranking depending on the number of principal component axis used for estimating the WAASB index. A euclidean distance-based dendrogram is used for grouping the genotype ranking for both genotypes and principal component axis. The second type of heatmap shows the genotype ranking depending on the WAASB/GY ratio. The ranks obtained with a ratio of 100/0 considers exclusively the stability for genotype ranking. On the other hand, a ratio of 0/100 considers exclusively the productivity for genotype ranking. Four clusters are estimated (1) unproductive and unstable genotypes; (2) productive, but unstable genotypes; (3) stable, but unproductive genotypes; and (4), productive and stable genotypes.

### Ranks of genotypes depending on the number PCAs used to estimate the WAASB

``` r
plot.WAASBYratio(WAASBYratio,
                 type = 1)
```

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-49-1.png" style="display: block; margin: auto;" />

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

<img src="D:\Desktop\WAASB\README_files/figure-markdown_github/unnamed-chunk-50-1.png" style="display: block; margin: auto;" />

``` r
#save to a *tiff file
plot.WAASBYratio(WAASBYratio,
                type = 2,
               export = T,
              file.type = "tiff",
             resolution = 600)
```

    ## png 
    ##   2

References
==========

Gauch, H.G, and Zobel R.W. 1988. “Predictive and Postdictive Success of Statistical Analyses of Yield Trials.” *Theor. Appl. Genet.* 76 (1): 1–10. doi:[10.1007/BF00288824](https://doi.org/10.1007/BF00288824).

Olivoto, T. 2018. “WAASBdata.” *Mendeley Data* V1. doi:[10.17632/2sjz32k3s3.1](https://doi.org/10.17632/2sjz32k3s3.1).

Piepho, H.P. 1994. “Best Linear Unbiased Prediction (Blup) for Regional Yield Trials: A Comparison to Additive Main Effects and Multiplicative Interaction (Ammi) Analysis.” *Theor. Appl. Genet.* 89 (5): 647–54. doi:[10.1007/BF00222462](https://doi.org/10.1007/BF00222462).
