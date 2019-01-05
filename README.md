Extending METAAB package
================
Tiago Olivoto
METAAB 1.0.0 2019-01-05

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
Introducing the `METAAB` R package
==================================

The METAAB R package was developed in R language and is distributed under the GPL 3.0 licence. Our main motivation for the package development was to provide friendly, easy and reliable functions provide de singular value decomposition of BLUP-interaction effects matrix as weel as to share codes used for traditional AMMI analysis and BLUP prediction. The package can be downloaded by clicking [here](https://github.com/TiagoOlivoto/METAAB/archive/master.zip) or running the following code in the R console:

Overview of the `METAAB` package.
---------------------------------

> Package `METAAB` available at `github`. Installation

``` r
# download the package from Github
devtools::install_github("TiagoOlivoto/METAAB")
```

Main functions of `METAAB` package.
-----------------------------------

<table>
<colgroup>
<col width="28%" />
<col width="71%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Functions</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><code>AMMI_indexes()</code></td>
<td align="left">AMMI-based stability indexes</td>
</tr>
<tr class="even">
<td align="left"><code>resca()</code></td>
<td align="left">Helper function to rescale a vector to have specified minimum and maximum values</td>
</tr>
<tr class="odd">
<td align="left"><code>Resende_indexes()</code></td>
<td align="left">Stability indexes proposed by Resende based on a mixed-effect model</td>
</tr>
<tr class="even">
<td align="left"><code>WAASB()</code></td>
<td align="left">Weighted Average of Absolute Scores from SVD of BLUP-based GEI effects</td>
</tr>
<tr class="odd">
<td align="left"><code>WAASBYratio()</code></td>
<td align="left">Different scenarios for weighting between mean performance and stability (mixed models)</td>
</tr>
<tr class="even">
<td align="left"><code>WAAS.AMMI()</code></td>
<td align="left">Weighted Average of Absolute Scores for AMMI analysis</td>
</tr>
<tr class="odd">
<td align="left"><code>WAASratio.AMMI()</code></td>
<td align="left">Different scenarios for weighting between mean performance and stability (fixed models)</td>
</tr>
<tr class="even">
<td align="left"><code>validation.blup()</code></td>
<td align="left">Cross-validation for BLUP prediction</td>
</tr>
<tr class="odd">
<td align="left"><code>validation.AMMIF()</code></td>
<td align="left">Cross-validation for all member of AMMI family models</td>
</tr>
<tr class="even">
<td align="left"><code>plot.blup()</code></td>
<td align="left">Plot the estimated BLUPs of genotypes</td>
</tr>
<tr class="odd">
<td align="left"><code>performs_ammi()</code></td>
<td align="left">Helper function for computing AMMI analysis</td>
</tr>
<tr class="even">
<td align="left"><code>plot.eigen()</code></td>
<td align="left">Plot the eigenvalues of singular value decomposition</td>
</tr>
<tr class="odd">
<td align="left"><code>theme_waasb()</code></td>
<td align="left">The default theme for the ggplot-based graphics generated in the package</td>
</tr>
<tr class="even">
<td align="left">Methods</td>
<td align="left">Description</td>
</tr>
<tr class="odd">
<td align="left"><code>autoplot()</code></td>
<td align="left">Plot for diagnostic of residuals</td>
</tr>
<tr class="even">
<td align="left"><code>plot.validation.AMMIF()</code></td>
<td align="left">Plot the RMSPD obtained in the cross-validation process</td>
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
<td align="left">Predict the response variable of an object of class <code>WAAS.AMMI</code></td>
</tr>
<tr class="even">
<td align="left"><code>summary.WAASB()</code></td>
<td align="left">Summarize an object of class <code>WAASB</code></td>
</tr>
<tr class="odd">
<td align="left"><code>summary.WAAS.AMMI()</code></td>
<td align="left">Summarize an object of class <code>WAAS.AMMI</code></td>
</tr>
</tbody>
</table>

Getting started
===============

The data set from experiment 1 (oat) was provided to make reproducible examples. This data is available in Olivoto ([2018](#ref-Olivoto:2018)). Other data sets can be used provided that the following columns are in the dataset: environment, genotype, block/replicate and response variable(s). The order of collumns doesn't need necessarily be the same since in each function the user will need to inform the columns to be evaluated.

``` r
# Importing data
require(METAAB)
require(ggplot2)
require(cowplot) # used in this material to arrange the graphics
dataset = read.csv("https://data.mendeley.com/datasets/2sjz32k3s3/2/files/1561de7f-b8fd-4093-b4d1-bfcef299dd22/WAASBdata.csv?dl=1")
```

Predictive accuracy
===================

The first steap in the article was to evaluate the predictive accucary of both AMMI and BLUP models. We will show in details how we can do that using the functions `validation.AMMIF()` and `validation.blup()` functions. The `validation.AMMIF()` function provides a complete cross-validation procedure for all member of AMMI model family (AMMI0-AMMIF) using replicate-based data, according to the diagram below: ![Diagrama demonstrando o proceso de validação cruzada dos membros da família AMMI de modelos.](validation.png). Automatically the first validation is carried out considering the AMMIF (all possible axis used). Considering this model, the original data set is split up into two sets: training set and validation set. The training set has all combinations (genotype x environment) with the number of replications informed in `nrepval` argument. The validation set has one replication that were not included in the training set. The splitting of the data set into training and validation sets depends on the design considered. For a Randomized Complete Block Design (default option) and the procedure we used in the article, completely blocks are randomly selected within environments, as suggested by Piepho ([1994](#ref-Piepho:1994)). The remaining block serves as validation data. If `design = "CRD"` is informed, thus declaring that a completely randomized design was used, single observations are randomized for each treatment (genotype-by-environment combination). This is the same procedure suggested by Gauch and Zobel ([1988](#ref-Gauch:1988)). The estimated values for each member of the AMMI model family in each re-sampling cycle are compared with the observed values in the validation data. Then, the Root Mean Square Prediction Difference is computed as follows:

$$
   RMSPD = 
   {\\left\[ {\\left( {\\sum\\nolimits\_{i = 1}^n {{{\\left( {{{\\hat y}\_{ij}} - {y\_{ij}}} \\right)}^2}} } \\right)/n} \\right\]^{0.5}}
$$

where $\\widehat{y}\_{ij}$ is the model predicted value; and *y*<sub>*i**j*</sub> is the observed value in the validation set. The number of random selection of blocks/replicates (*n*) is defined in the argument `nboot`. At the end of the *n* cycles for all models, a list with all estimated RMSPD and the average of RMSPD is returned.

The function `validation.blup` provides a cross-validation of replicate-based data using mixed models. By default, complete blocks are randomly selected for each environment. The procedure for computing the RSME is identical to the above function.

The following code computes the cross-validation of the oat data set based on 1000 re-sampling procedures. This number can be changed.

``` r

AMMIweat = validation.AMMIF(dataset,
                            resp = GY,
                            gen = GEN,
                            env = ENV,
                            rep = BLOCK,
                            nboot = 10,
                            nrepval = 2)

# cross-validation for BLUP model
BLUPweat = validation.blup(dataset,
                            resp = GY,
                            gen = GEN,
                            env = ENV,
                            rep = BLOCK,
                            nboot = 10,
                            nrepval = 2)
```

A progress bar is shown by default (the examples are below). Thus, it is possible to verify the status of the cross-validation process. If necessary, the progress bar can be disabled by informing the argument `progbar = FALSE` in the function.

![Progress bar for AMMI model family cross-validation process.](printAMMI.png)

![Progress bar for BLUP cross-validation process.](printBLUP.png)

Printing the means of RMSPD estimates
-------------------------------------

All outputs shown in tables were laid out using the R package kableExtra. The pdf file of this package can be found [here](https://haozhu233.github.io/kableExtra/awesome_table_in_pdf.pdf).

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
0.4123
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
0.4183
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
0.4195
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
AMMI7
</td>
<td style="text-align:right;">
0.4207
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
0.4215
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
AMMI6
</td>
<td style="text-align:right;">
0.4239
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
AMMI8
</td>
<td style="text-align:right;">
0.4243
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
AMMI2
</td>
<td style="text-align:right;">
0.4301
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
AMMIF
</td>
<td style="text-align:right;">
0.4317
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
AMMI1
</td>
<td style="text-align:right;">
0.4415
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
AMMI0
</td>
<td style="text-align:right;">
0.4478
</td>
<td style="text-align:left;">
Wheat
</td>
</tr>
</tbody>
</table>
The results shown above are the means of the 1000 RMSPD estimates for each tested model and are presented from the most accurate model (smallest RSME value) to the least accurate model (highest RMSPD value).

We can print the result from `RMSPDweat` in basically two distinct ways. First, printing the results in the R environment using `print(RMSPDweat)` (this is not a good idea since we probably will need work with these results). The second option is to export the results to an editable file, such as a .csv, or .xlsx file. We have R packages that facilitate this procedure. Let's do it. For example, to export the results to a .csv file, the command to be run is: `utils::write.csv(RMSPDweat, file = "myfile.csv")`. This command will create a .csv file called "myfile" in the R directory. To export the results to a .xlsx file, the package `xlsx` is needed. After properly installed, the command will be then `xlsx::write.xls(RMSPDweat, file = "myfile2.xlsx")`.

Plotting the RMSPD values
-------------------------

The values of the RMSPD estimates obtained in the cross-validation process may be plotted using the function`plot()` may be used.

``` r
# binding AMMI and BLUP RMSPDs
RMSPDweat = list(RMSPD = rbind(AMMIweat$RMSPD, 
                               BLUPweat$RMSPD))
class(RMSPDweat) = "validation.AMMIF"

# Plotting the RMSPD values
p1 = plot(RMSPDweat)
p2 = plot(RMSPDweat, width.boxplot = 0.5, col.boxplot = "transparent")
plot_grid(p1, p2, labels = c("p1", "p2"))
```

<img src="D:\Desktop\METAAB\README_files/figure-markdown_github/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

Six statistics are shown in this boxplot. The mean (black rhombus), the median (black line), the lower and upper hinges that correspond sto the first and third quartiles (the 25th and 75th percentiles, respectively). The upper whisker extends from the hinge to the largest value no further than 1.5 × *I**Q**R* from the hinge (where IQR is the inter-quartile range). The lower whisker extends from the hinge to the smallest value at most 1.5 × *I**Q**R* of the hinge. Data beyond the end of the whiskers are considered outlying points. If the condition `violin = TRUE`, a violin plot is added along with the boxplot. A violin plot is a compact display of a continuous distribution displayed in the same way as a boxplot.

Estimating the response variable using traditional AMMI model
=============================================================

An interesting feature of `WAASB` package for traditional AMMI model estimation is the implementation of the S3 method `predict()`. The response variable of a two-way table (for example, the yield of *m* genotypes in *n* environments) may be estimated using the function `predict(model)`, where `model` is an object of class `WAAS.AMMI`. This estimation is based on the number of multiplicative terms declared in the function. If `naxis = 0` is declared, only the main effects (AMMI0) are considered. In this case, the estimated mean will be the estimate from OLS estimation. If `naxis = 1`, the AMMI1 (with one multiplicative term) is used for estimating the response variable. If `naxis = min(gen-1;env-1)`, the AMMIF is fitted. A summary of all possible AMMI models is presented below.

<table style="width:100%;">
<colgroup>
<col width="24%" />
<col width="75%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Member of AMMI family</th>
<th align="left">Espected response of the <em>i</em>-th genotype in the <em>j</em>th environment</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">AMMI0</td>
<td align="left"><span class="math inline">$\hat{y}_{ij} = \bar{y}_{i.} + \bar{y}_{.j} - \bar{y}_{..}$</span></td>
</tr>
<tr class="even">
<td align="left">AMMI1</td>
<td align="left"><span class="math inline">$\hat{y}_{ij} = \bar{y}_{i.} + \bar{y}_{.j} - \bar{y}_{..} +\lambda_1 a_{i1}t_{j1}$</span></td>
</tr>
<tr class="odd">
<td align="left">AMMI2</td>
<td align="left"><span class="math inline">$\hat{y}_{ij} = \bar{y}_{i.} + \bar{y}_{.j} - \bar{y}_{..} +\lambda_1 a_{i1}t_{j1}+\lambda_2 a_{i2}t_{j2}$</span></td>
</tr>
<tr class="even">
<td align="left">...</td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left">AMMIF</td>
<td align="left"><span class="math inline">$\hat{y}_{ij} = \bar{y}_{i.} + \bar{y}_{.j} - \bar{y}_{..} +\lambda_1 a_{i1}t_{j1}+\lambda_2 a_{i2}t_{j2}+...+\lambda_p a_{ip}t_{jp}$</span></td>
</tr>
</tbody>
</table>

The number of axes used in the estimation must be carefully chosen. Procedures based on postdictive success, such as Gollobs's test (Gollob [1968](#ref-Gollob:1968)) or predictive success, such as cross-validation procedures (Piepho [1994](#ref-Piepho:1994)) should be used. This package provides both. `WAAS.AMMI` function compute traditional AMMI analysis showing the number of significant axes according Gollobs's test. On the other hand, `validation.AMMIF` function provides a cross-validation, estimating the RMSPD for all AMMI model family based on re-sampling procedures, considering an completely randomized desing or a randomized complete block design.

-   Estimating the yield of the 10 oat cultivars in 16 environments using 5 multiplicative terms.
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
    NF2010
    </td>
    <td style="text-align:left;">
    G1
    </td>
    <td style="text-align:right;">
    1.898
    </td>
    <td style="text-align:right;">
    -0.0240
    </td>
    <td style="text-align:right;">
    1.922
    </td>
    <td style="text-align:right;">
    -0.09148
    </td>
    <td style="text-align:right;">
    1.830
    </td>
    <td style="text-align:right;">
    1.922
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    NF2010
    </td>
    <td style="text-align:left;">
    G10
    </td>
    <td style="text-align:right;">
    2.244
    </td>
    <td style="text-align:right;">
    0.4395
    </td>
    <td style="text-align:right;">
    1.804
    </td>
    <td style="text-align:right;">
    0.40405
    </td>
    <td style="text-align:right;">
    2.208
    </td>
    <td style="text-align:right;">
    1.804
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    NF2010
    </td>
    <td style="text-align:left;">
    G2
    </td>
    <td style="text-align:right;">
    1.988
    </td>
    <td style="text-align:right;">
    -0.0132
    </td>
    <td style="text-align:right;">
    2.001
    </td>
    <td style="text-align:right;">
    -0.03240
    </td>
    <td style="text-align:right;">
    1.968
    </td>
    <td style="text-align:right;">
    2.001
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    NF2010
    </td>
    <td style="text-align:left;">
    G3
    </td>
    <td style="text-align:right;">
    2.159
    </td>
    <td style="text-align:right;">
    -0.0803
    </td>
    <td style="text-align:right;">
    2.239
    </td>
    <td style="text-align:right;">
    -0.03960
    </td>
    <td style="text-align:right;">
    2.199
    </td>
    <td style="text-align:right;">
    2.239
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    NF2010
    </td>
    <td style="text-align:left;">
    G4
    </td>
    <td style="text-align:right;">
    1.981
    </td>
    <td style="text-align:right;">
    -0.0140
    </td>
    <td style="text-align:right;">
    1.995
    </td>
    <td style="text-align:right;">
    -0.07142
    </td>
    <td style="text-align:right;">
    1.924
    </td>
    <td style="text-align:right;">
    1.995
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    NF2010
    </td>
    <td style="text-align:left;">
    G5
    </td>
    <td style="text-align:right;">
    1.657
    </td>
    <td style="text-align:right;">
    -0.2070
    </td>
    <td style="text-align:right;">
    1.864
    </td>
    <td style="text-align:right;">
    -0.06381
    </td>
    <td style="text-align:right;">
    1.800
    </td>
    <td style="text-align:right;">
    1.864
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    NF2010
    </td>
    <td style="text-align:left;">
    G6
    </td>
    <td style="text-align:right;">
    1.758
    </td>
    <td style="text-align:right;">
    -0.0897
    </td>
    <td style="text-align:right;">
    1.847
    </td>
    <td style="text-align:right;">
    -0.13428
    </td>
    <td style="text-align:right;">
    1.713
    </td>
    <td style="text-align:right;">
    1.847
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    NF2010
    </td>
    <td style="text-align:left;">
    G7
    </td>
    <td style="text-align:right;">
    2.551
    </td>
    <td style="text-align:right;">
    0.5429
    </td>
    <td style="text-align:right;">
    2.008
    </td>
    <td style="text-align:right;">
    0.58685
    </td>
    <td style="text-align:right;">
    2.595
    </td>
    <td style="text-align:right;">
    2.008
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    NF2010
    </td>
    <td style="text-align:left;">
    G8
    </td>
    <td style="text-align:right;">
    2.262
    </td>
    <td style="text-align:right;">
    -0.0757
    </td>
    <td style="text-align:right;">
    2.337
    </td>
    <td style="text-align:right;">
    -0.14999
    </td>
    <td style="text-align:right;">
    2.187
    </td>
    <td style="text-align:right;">
    2.337
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">
    NF2010
    </td>
    <td style="text-align:left;">
    G9
    </td>
    <td style="text-align:right;">
    1.393
    </td>
    <td style="text-align:right;">
    -0.4784
    </td>
    <td style="text-align:right;">
    1.872
    </td>
    <td style="text-align:right;">
    -0.40792
    </td>
    <td style="text-align:right;">
    1.464
    </td>
    <td style="text-align:right;">
    1.872
    </td>
    </tr>
    </tbody>
    </table>

Only the first ten values are shown. The following values are presented: **ENV** is the environment; **GEN** is the genotype; **Y** is the response variable; **resOLS** is the residual ($\\hat{z}\_{ij}$) estimated by the Ordinary Least Square (OLS), where $\\hat{z}\_{ij} = y\_{ij} - \\bar{y}\_{i.} - \\bar{y}\_{.j} + \\bar{y}\_{ij}$; **Ypred** is the predicted value by OLS ($\\hat{y}\_{ij} = y\_{ij} -\\hat{z}\_{ij}$); **ResAMMI** is the residual estimated by the AMMI model ($\\hat{a}\_{ij}$) considering the number of multiplicative terms informed in the function (in this case 5), where $\\hat{a}\_{ij} = \\lambda\_1\\alpha\_{i1}\\tau\_{j1}+...+\\lambda\_5\\alpha\_{i5}\\tau\_{j5}$; **YpredAMMI** is the predicted value by AMMI model $\\hat{ya}\_{ij} = \\bar{y}\_{i.} + \\bar{y}\_{.j} - \\bar{y}\_{ij}+\\hat{a}\_{ij}$; and **AMMI0** is the predicted value when no multiplicative terms are used, i.e., $\\hat{y}\_{ij} = \\bar{y}\_{i.} + \\bar{y}\_{.j} - \\bar{y}\_{ij}$.

Estimating the WAAS
===================

The `WAAS.AMMI` function computes the Weighted Average of Absolute Scores considering (i) all principal component axes that were significant (*p* &lt; 0.05 by default); or (ii) declaring a specific number of axes to be used, according to the following equation:

$$
        WAAS\_i  = 
        \\sum\_{k = 1}^{p} |IPCA\_{ik} \\times EP\_k|/ \\sum\_{k = 1}^{p}EP\_k
$$

where *W**A**A**S*<sub>*i*</sub> is the weighted average of absolute scores of the *i*th genotype; *P**C**A*<sub>*i**k*</sub> is the score of the *i*th genotype in the *k*th IPCA; and *E**P*<sub>*k*</sub> is the explained variance of the *k*th IPCA for *k* = 1, 2, ..,*p*, considering *p* the number of significant PCAs, or a declared number of PCAs. The following functions may be used to do that.

Number of axes based on F-test
------------------------------

In this example only IPCAs with *P*-value &lt; 0.05 will be considered in the WAAS estimation.

``` r
# Assuming equal weights for productivity and stability (default)
WAAS1 = WAAS.AMMI(dataset,
                  resp = GY,
                  gen = GEN,
                  env = ENV,
                  rep = BLOCK)
```

``` r
# printing the WAAS object
options(digits = 3)
data = WAAS1$individual$individual
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", full_width = F, position = "left", font_size = 12)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
ENV
</th>
<th style="text-align:right;">
Mean
</th>
<th style="text-align:right;">
MSblock
</th>
<th style="text-align:right;">
MSgen
</th>
<th style="text-align:right;">
MSres
</th>
<th style="text-align:right;">
Fcal(Blo)
</th>
<th style="text-align:right;">
Pr&gt;F(Blo)
</th>
<th style="text-align:right;">
Fcal(Gen)
</th>
<th style="text-align:right;">
Pr&gt;F(Gen)
</th>
<th style="text-align:right;">
CV(%)
</th>
<th style="text-align:right;">
h2
</th>
<th style="text-align:right;">
AS
</th>
<th style="text-align:right;">
R2
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:right;">
1.99
</td>
<td style="text-align:right;">
0.381
</td>
<td style="text-align:right;">
0.337
</td>
<td style="text-align:right;">
0.091
</td>
<td style="text-align:right;">
4.189
</td>
<td style="text-align:right;">
0.032
</td>
<td style="text-align:right;">
3.71
</td>
<td style="text-align:right;">
0.009
</td>
<td style="text-align:right;">
15.16
</td>
<td style="text-align:right;">
0.730
</td>
<td style="text-align:right;">
0.854
</td>
<td style="text-align:right;">
0.787
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2011
</td>
<td style="text-align:right;">
2.54
</td>
<td style="text-align:right;">
0.817
</td>
<td style="text-align:right;">
0.215
</td>
<td style="text-align:right;">
0.028
</td>
<td style="text-align:right;">
29.370
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
7.72
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
6.58
</td>
<td style="text-align:right;">
0.870
</td>
<td style="text-align:right;">
0.933
</td>
<td style="text-align:right;">
0.885
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2012
</td>
<td style="text-align:right;">
3.06
</td>
<td style="text-align:right;">
0.583
</td>
<td style="text-align:right;">
0.679
</td>
<td style="text-align:right;">
0.111
</td>
<td style="text-align:right;">
5.253
</td>
<td style="text-align:right;">
0.016
</td>
<td style="text-align:right;">
6.12
</td>
<td style="text-align:right;">
0.001
</td>
<td style="text-align:right;">
10.90
</td>
<td style="text-align:right;">
0.837
</td>
<td style="text-align:right;">
0.915
</td>
<td style="text-align:right;">
0.860
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2013
</td>
<td style="text-align:right;">
2.17
</td>
<td style="text-align:right;">
0.654
</td>
<td style="text-align:right;">
0.296
</td>
<td style="text-align:right;">
0.027
</td>
<td style="text-align:right;">
24.491
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
11.09
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
7.51
</td>
<td style="text-align:right;">
0.910
</td>
<td style="text-align:right;">
0.954
</td>
<td style="text-align:right;">
0.917
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2014
</td>
<td style="text-align:right;">
1.37
</td>
<td style="text-align:right;">
0.377
</td>
<td style="text-align:right;">
0.151
</td>
<td style="text-align:right;">
0.105
</td>
<td style="text-align:right;">
3.595
</td>
<td style="text-align:right;">
0.049
</td>
<td style="text-align:right;">
1.44
</td>
<td style="text-align:right;">
0.244
</td>
<td style="text-align:right;">
23.68
</td>
<td style="text-align:right;">
0.304
</td>
<td style="text-align:right;">
0.552
</td>
<td style="text-align:right;">
0.590
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2015
</td>
<td style="text-align:right;">
1.61
</td>
<td style="text-align:right;">
0.092
</td>
<td style="text-align:right;">
0.320
</td>
<td style="text-align:right;">
0.054
</td>
<td style="text-align:right;">
1.717
</td>
<td style="text-align:right;">
0.208
</td>
<td style="text-align:right;">
5.98
</td>
<td style="text-align:right;">
0.001
</td>
<td style="text-align:right;">
14.38
</td>
<td style="text-align:right;">
0.833
</td>
<td style="text-align:right;">
0.913
</td>
<td style="text-align:right;">
0.857
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2016
</td>
<td style="text-align:right;">
2.91
</td>
<td style="text-align:right;">
0.077
</td>
<td style="text-align:right;">
0.713
</td>
<td style="text-align:right;">
0.099
</td>
<td style="text-align:right;">
0.772
</td>
<td style="text-align:right;">
0.477
</td>
<td style="text-align:right;">
7.18
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
10.83
</td>
<td style="text-align:right;">
0.861
</td>
<td style="text-align:right;">
0.928
</td>
<td style="text-align:right;">
0.878
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2017
</td>
<td style="text-align:right;">
1.78
</td>
<td style="text-align:right;">
0.103
</td>
<td style="text-align:right;">
0.131
</td>
<td style="text-align:right;">
0.075
</td>
<td style="text-align:right;">
1.373
</td>
<td style="text-align:right;">
0.279
</td>
<td style="text-align:right;">
1.73
</td>
<td style="text-align:right;">
0.153
</td>
<td style="text-align:right;">
15.40
</td>
<td style="text-align:right;">
0.423
</td>
<td style="text-align:right;">
0.650
</td>
<td style="text-align:right;">
0.634
</td>
</tr>
<tr>
<td style="text-align:left;">
WF2010
</td>
<td style="text-align:right;">
2.52
</td>
<td style="text-align:right;">
0.065
</td>
<td style="text-align:right;">
0.337
</td>
<td style="text-align:right;">
0.144
</td>
<td style="text-align:right;">
0.453
</td>
<td style="text-align:right;">
0.643
</td>
<td style="text-align:right;">
2.34
</td>
<td style="text-align:right;">
0.059
</td>
<td style="text-align:right;">
15.05
</td>
<td style="text-align:right;">
0.573
</td>
<td style="text-align:right;">
0.757
</td>
<td style="text-align:right;">
0.701
</td>
</tr>
<tr>
<td style="text-align:left;">
WF2011
</td>
<td style="text-align:right;">
3.18
</td>
<td style="text-align:right;">
0.697
</td>
<td style="text-align:right;">
0.207
</td>
<td style="text-align:right;">
0.179
</td>
<td style="text-align:right;">
3.907
</td>
<td style="text-align:right;">
0.039
</td>
<td style="text-align:right;">
1.16
</td>
<td style="text-align:right;">
0.376
</td>
<td style="text-align:right;">
13.29
</td>
<td style="text-align:right;">
0.136
</td>
<td style="text-align:right;">
0.369
</td>
<td style="text-align:right;">
0.536
</td>
</tr>
<tr>
<td style="text-align:left;">
WF2012
</td>
<td style="text-align:right;">
4.06
</td>
<td style="text-align:right;">
0.489
</td>
<td style="text-align:right;">
0.335
</td>
<td style="text-align:right;">
0.179
</td>
<td style="text-align:right;">
2.731
</td>
<td style="text-align:right;">
0.092
</td>
<td style="text-align:right;">
1.87
</td>
<td style="text-align:right;">
0.123
</td>
<td style="text-align:right;">
10.41
</td>
<td style="text-align:right;">
0.466
</td>
<td style="text-align:right;">
0.683
</td>
<td style="text-align:right;">
0.652
</td>
</tr>
<tr>
<td style="text-align:left;">
WF2013
</td>
<td style="text-align:right;">
3.67
</td>
<td style="text-align:right;">
0.116
</td>
<td style="text-align:right;">
0.531
</td>
<td style="text-align:right;">
0.138
</td>
<td style="text-align:right;">
0.846
</td>
<td style="text-align:right;">
0.445
</td>
<td style="text-align:right;">
3.86
</td>
<td style="text-align:right;">
0.007
</td>
<td style="text-align:right;">
10.10
</td>
<td style="text-align:right;">
0.741
</td>
<td style="text-align:right;">
0.861
</td>
<td style="text-align:right;">
0.794
</td>
</tr>
<tr>
<td style="text-align:left;">
WF2014
</td>
<td style="text-align:right;">
2.51
</td>
<td style="text-align:right;">
0.067
</td>
<td style="text-align:right;">
0.408
</td>
<td style="text-align:right;">
0.231
</td>
<td style="text-align:right;">
0.290
</td>
<td style="text-align:right;">
0.751
</td>
<td style="text-align:right;">
1.76
</td>
<td style="text-align:right;">
0.146
</td>
<td style="text-align:right;">
19.17
</td>
<td style="text-align:right;">
0.433
</td>
<td style="text-align:right;">
0.658
</td>
<td style="text-align:right;">
0.638
</td>
</tr>
<tr>
<td style="text-align:left;">
WF2015
</td>
<td style="text-align:right;">
3.11
</td>
<td style="text-align:right;">
0.166
</td>
<td style="text-align:right;">
0.550
</td>
<td style="text-align:right;">
0.067
</td>
<td style="text-align:right;">
2.481
</td>
<td style="text-align:right;">
0.112
</td>
<td style="text-align:right;">
8.24
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
8.32
</td>
<td style="text-align:right;">
0.879
</td>
<td style="text-align:right;">
0.937
</td>
<td style="text-align:right;">
0.892
</td>
</tr>
<tr>
<td style="text-align:left;">
WF2016
</td>
<td style="text-align:right;">
3.91
</td>
<td style="text-align:right;">
0.219
</td>
<td style="text-align:right;">
0.526
</td>
<td style="text-align:right;">
0.066
</td>
<td style="text-align:right;">
3.297
</td>
<td style="text-align:right;">
0.060
</td>
<td style="text-align:right;">
7.93
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
6.59
</td>
<td style="text-align:right;">
0.874
</td>
<td style="text-align:right;">
0.935
</td>
<td style="text-align:right;">
0.888
</td>
</tr>
<tr>
<td style="text-align:left;">
WF2017
</td>
<td style="text-align:right;">
2.66
</td>
<td style="text-align:right;">
0.160
</td>
<td style="text-align:right;">
0.135
</td>
<td style="text-align:right;">
0.059
</td>
<td style="text-align:right;">
2.733
</td>
<td style="text-align:right;">
0.092
</td>
<td style="text-align:right;">
2.30
</td>
<td style="text-align:right;">
0.063
</td>
<td style="text-align:right;">
9.08
</td>
<td style="text-align:right;">
0.565
</td>
<td style="text-align:right;">
0.752
</td>
<td style="text-align:right;">
0.697
</td>
</tr>
</tbody>
</table>
The above table shows the within-environment ANOVA considering a fixed-effect model. For each environment the Mean Squares for block, genotypes and error are shown. Estimated F-value and the probability error are also shown for block and genotype effects. Some measures of experimental precision are calculated, namelly, coefficient of variation, $CV = (\\sqrt{MS\_{res}}/Mean) \\times 100$; the heritability, *h*2 = (*M**S*<sub>*g**e**n*</sub> − *M**S*<sub>*r**e**s*</sub>)/*M**S*<sub>*g**e**n*</sub>; the accuracy of selection, $As = \\sqrt{h2}$; and the coefficient of determination (R2).

``` r
# printing the WAAS object
data = WAAS1$anova
kable(data,  align  = "l", booktabs = T, format = "html", linesep = "") %>%
kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12) %>%
row_spec(9, bold = T) %>%
add_indent(c(5:13))
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:left;">
Df
</th>
<th style="text-align:left;">
Sum Sq
</th>
<th style="text-align:left;">
Mean Sq
</th>
<th style="text-align:left;">
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
<td style="text-align:left;">
15
</td>
<td style="text-align:left;">
285.910
</td>
<td style="text-align:left;">
19.061
</td>
<td style="text-align:left;">
60.24
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
<td style="text-align:left;">
32
</td>
<td style="text-align:left;">
10.125
</td>
<td style="text-align:left;">
0.316
</td>
<td style="text-align:left;">
3.07
</td>
<td style="text-align:left;">
0.0000003177891856598183
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
<td style="text-align:left;">
9
</td>
<td style="text-align:left;">
13.084
</td>
<td style="text-align:left;">
1.454
</td>
<td style="text-align:left;">
14.09
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
<td style="text-align:left;">
135
</td>
<td style="text-align:left;">
39.755
</td>
<td style="text-align:left;">
0.294
</td>
<td style="text-align:left;">
2.85
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
<td style="text-align:left;">
23
</td>
<td style="text-align:left;">
13.062
</td>
<td style="text-align:left;">
0.568
</td>
<td style="text-align:left;">
5.50
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
<td style="text-align:left;">
21
</td>
<td style="text-align:left;">
10.687
</td>
<td style="text-align:left;">
0.509
</td>
<td style="text-align:left;">
4.93
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
<td style="text-align:left;">
19
</td>
<td style="text-align:left;">
6.803
</td>
<td style="text-align:left;">
0.358
</td>
<td style="text-align:left;">
3.47
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
<td style="text-align:left;">
17
</td>
<td style="text-align:left;">
3.685
</td>
<td style="text-align:left;">
0.217
</td>
<td style="text-align:left;">
2.10
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
<td style="text-align:left;font-weight: bold;">
15
</td>
<td style="text-align:left;font-weight: bold;">
2.670
</td>
<td style="text-align:left;font-weight: bold;">
0.178
</td>
<td style="text-align:left;font-weight: bold;">
1.72
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
<td style="text-align:left;">
13
</td>
<td style="text-align:left;">
1.454
</td>
<td style="text-align:left;">
0.112
</td>
<td style="text-align:left;">
1.08
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
<td style="text-align:left;">
11
</td>
<td style="text-align:left;">
0.727
</td>
<td style="text-align:left;">
0.066
</td>
<td style="text-align:left;">
0.64
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
<td style="text-align:left;">
9
</td>
<td style="text-align:left;">
0.465
</td>
<td style="text-align:left;">
0.052
</td>
<td style="text-align:left;">
0.50
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
<td style="text-align:left;">
7
</td>
<td style="text-align:left;">
0.203
</td>
<td style="text-align:left;">
0.029
</td>
<td style="text-align:left;">
0.28
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
<td style="text-align:left;">
288
</td>
<td style="text-align:left;">
29.725
</td>
<td style="text-align:left;">
0.103
</td>
<td style="text-align:left;">
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
<td style="text-align:left;">
479
</td>
<td style="text-align:left;">
378.599
</td>
<td style="text-align:left;">
0.790
</td>
<td style="text-align:left;">
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

> printing the METAAB object

``` r
options(digits = 4)
data = WAAS1$model[, c(1:3,13:17, 21:22)]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)
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
22.043
</td>
<td style="text-align:right;">
100.00
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
61.02
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
G10
</td>
<td style="text-align:right;">
2.506
</td>
<td style="text-align:right;">
0.5003
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.00
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
36.888
</td>
<td style="text-align:right;">
58.28
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
47.58
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
G3
</td>
<td style="text-align:right;">
2.941
</td>
<td style="text-align:right;">
0.2037
</td>
<td style="text-align:right;">
81.560
</td>
<td style="text-align:right;">
91.96
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
86.76
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
35.802
</td>
<td style="text-align:right;">
40.97
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
38.39
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
G5
</td>
<td style="text-align:right;">
2.566
</td>
<td style="text-align:right;">
0.2143
</td>
<td style="text-align:right;">
11.152
</td>
<td style="text-align:right;">
88.68
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
49.92
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
G6
</td>
<td style="text-align:right;">
2.549
</td>
<td style="text-align:right;">
0.2439
</td>
<td style="text-align:right;">
8.108
</td>
<td style="text-align:right;">
79.50
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
43.81
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
G7
</td>
<td style="text-align:right;">
2.710
</td>
<td style="text-align:right;">
0.3843
</td>
<td style="text-align:right;">
38.256
</td>
<td style="text-align:right;">
35.97
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
37.11
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
G8
</td>
<td style="text-align:right;">
3.039
</td>
<td style="text-align:right;">
0.2465
</td>
<td style="text-align:right;">
100.000
</td>
<td style="text-align:right;">
78.69
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
89.34
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
12.676
</td>
<td style="text-align:right;">
16.45
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
14.56
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
NF2010
</td>
<td style="text-align:right;">
1.989
</td>
<td style="text-align:right;">
0.2973
</td>
<td style="text-align:right;">
23.021
</td>
<td style="text-align:right;">
38.98
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
31.00
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
NF2011
</td>
<td style="text-align:right;">
2.536
</td>
<td style="text-align:right;">
0.2263
</td>
<td style="text-align:right;">
43.329
</td>
<td style="text-align:right;">
62.92
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
53.13
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
NF2012
</td>
<td style="text-align:right;">
3.057
</td>
<td style="text-align:right;">
0.4128
</td>
<td style="text-align:right;">
62.623
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
31.31
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
NF2013
</td>
<td style="text-align:right;">
2.175
</td>
<td style="text-align:right;">
0.1623
</td>
<td style="text-align:right;">
29.930
</td>
<td style="text-align:right;">
84.54
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
57.23
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
NF2014
</td>
<td style="text-align:right;">
1.368
</td>
<td style="text-align:right;">
0.1164
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
100.00
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
50.00
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
NF2015
</td>
<td style="text-align:right;">
1.609
</td>
<td style="text-align:right;">
0.1663
</td>
<td style="text-align:right;">
8.912
</td>
<td style="text-align:right;">
83.17
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
46.04
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
NF2016
</td>
<td style="text-align:right;">
2.910
</td>
<td style="text-align:right;">
0.3235
</td>
<td style="text-align:right;">
57.170
</td>
<td style="text-align:right;">
30.15
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
43.66
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
NF2017
</td>
<td style="text-align:right;">
1.782
</td>
<td style="text-align:right;">
0.1588
</td>
<td style="text-align:right;">
15.344
</td>
<td style="text-align:right;">
85.70
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
50.52
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
WF2010
</td>
<td style="text-align:right;">
2.521
</td>
<td style="text-align:right;">
0.1787
</td>
<td style="text-align:right;">
42.746
</td>
<td style="text-align:right;">
78.99
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
60.87
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
WF2011
</td>
<td style="text-align:right;">
3.180
</td>
<td style="text-align:right;">
0.2086
</td>
<td style="text-align:right;">
67.197
</td>
<td style="text-align:right;">
68.90
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
68.05
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
WF2012
</td>
<td style="text-align:right;">
4.064
</td>
<td style="text-align:right;">
0.3627
</td>
<td style="text-align:right;">
100.000
</td>
<td style="text-align:right;">
16.93
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
58.47
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
WF2013
</td>
<td style="text-align:right;">
3.675
</td>
<td style="text-align:right;">
0.2564
</td>
<td style="text-align:right;">
85.572
</td>
<td style="text-align:right;">
52.76
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
69.17
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
WF2014
</td>
<td style="text-align:right;">
2.507
</td>
<td style="text-align:right;">
0.2429
</td>
<td style="text-align:right;">
42.241
</td>
<td style="text-align:right;">
57.35
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
49.80
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
WF2015
</td>
<td style="text-align:right;">
3.107
</td>
<td style="text-align:right;">
0.3486
</td>
<td style="text-align:right;">
64.497
</td>
<td style="text-align:right;">
21.68
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
43.09
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
WF2016
</td>
<td style="text-align:right;">
3.910
</td>
<td style="text-align:right;">
0.2115
</td>
<td style="text-align:right;">
94.296
</td>
<td style="text-align:right;">
67.92
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
81.11
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
WF2017
</td>
<td style="text-align:right;">
2.663
</td>
<td style="text-align:right;">
0.1167
</td>
<td style="text-align:right;">
48.031
</td>
<td style="text-align:right;">
99.90
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
73.97
</td>
<td style="text-align:right;">
2
</td>
</tr>
</tbody>
</table>
In this example, the scores of the nine PCA were not shown. The output generated by the `WAAS.AMMI` function shows the following results: **type**, genotype (GEN) or environment (ENV); **Code**, the code attributed to each level of the factors; **Y**, the response variable (in this case the grain yield); **WAAS** the weighted average of the absolute scores, estimated with all PCA axes with *P*-value ≤ 0.05; **PctWAAS** and **PctResp** that are the percentage values for the WAAS and Y, respectively; **OrResp** and **OrWAAS** that are the ranks attributed to the genotype and environment regarding the Y or WAAS, respectively; **WAASY** is the weighted average of absolute scores and response variable. In this case, considering equal weights for PctResp and PctWAAS, the WAASY for G1 is estimated by: *W**A**A**S*<sub>*G*1</sub> = \[(86.32 × 50)+(98.88 × 50)\]/50 + 50 = 92.60. Then the \*\*OrWAASY\* is the rank for the WAASY value. The genotype (or environment) with the largest WAASY value has the first ranked. See [Estimating the WAASBY index](#estimating-the-waasby-index) for a detailed explanation.

In this example, the scores of the nine PCAs were not shown. The output generated by the `WAAS.AMMI` function shows the following results: **type**, genotype (GEN) or environment (ENV); **Code**, the code attributed to each level of the factors; **Y**, the response variable (in this case the grain yield); **WAAS** the weighted average of the absolute scores, estimated with all significant IPCA; **PctWAAS** and **PctResp**, the rescaled variable for the WAAS and Y, respectively; **OrResp** and **OrWAAS**, the ranks attributed to the genotype and environment regarding the Y or WAAS, respectively; **WAASY** is the simultanteous selection index that weights between response variable and stability. In this case, considering equal weights for PctResp and PctWAAS, the WAASY for G1 would be then: *W**A**A**S**Y*<sub>*G*1</sub> = \[(22.043 × 50)+(100 × 50)\]/50 + 50 = 92.60, which, with equal weights is equivalente to the aritmetic mean of pctResp and PctWAAS. Then the *OrWAASY* is the rank for the WAASY value. The genotype (or environment) with the largest WAASY value has the first ranked. Please, refer to [Estimating the WAASBY index](#estimating-the-waasby-index) for a detailed explanation.

Number of axes declared manually
--------------------------------

The second option to compute the WAAS is by manually declaring a specific number of multiplicative terms. In this case, the number of terms declared is used independently of its significance. Let us, for the moment, assume that after a cross-validation procedure the AMMI7 was the most predictively accurate AMMI model and the researcher will use this model. The additional argument `naxis` in the function `WAAS.AMMI` is then used to overwrite the default chose of significant terms.

``` r
WAAS2 = WAAS.AMMI(dataset,
                  resp = GY,
                  gen = GEN,
                  env = ENV,
                  rep = BLOCK,
                  naxis = 7)
```

The only difference in this output compared to those from [section 5.1](#assuming-a-given-probability-error-for-chosing-the-number-of-axes) is that here we declared that seven PCA axes should be used for computing the WAAS value. Thus, only the values of WAAS, OrWAAS, WAASY and OrWAASY may have significant changes.

Other AMMI-based stability indexes
==================================

The following AMMI-based stability indexes tested in the article may be computed using the function `AMMI_indexes()`:

-   **AMMI stability value, ASV, (Purchase, Hatting, and Deventer [2000](#ref-Purchase2000)).**

$$
ASV = \\sqrt {{{\\left\[ {\\frac{{IPCA{1\_{ss}}}}{{IPCA{2\_{ss}}}} \\times \\left( {IPCA{1\_{score}}} \\right)} \\right\]}^2} + {{\\left( {IPCA{2\_{score}}} \\right)}^2}}
$$

-   **Sums of the absolute value of the IPCA scores**

*S**I**P**C*<sub>*i*</sub> = ∑<sub>*k* = 1</sub><sup>*P*</sup>||*λ*<sub>*k*</sub><sup>0.5</sup>*a*<sub>*i**k*</sub>|

-   **Averages of the squared eigenvector values**

*E**V*<sub>*i*</sub> = ∑<sub>*k* = 1</sub><sup>*P*</sup>a<sub>*i**k*</sub><sup>2</sup>/*P*
 described by Sneller, Kilgore-Norquest, and Dombek ([1997](#ref-Sneller1997)), where *P* is the number of IPCA retained via F-tests;

-   **absolute value of the relative contribution of IPCAs to the interaction (Zali et al. [2012](#ref-Zali2012)).**

*Z**a*<sub>*i*</sub> = ∑<sub>*k* = 1</sub><sup>*P*</sup>*θ*<sub>*k*</sub>*a*<sub>*i**k*</sub>

where *θ*<sub>*k*</sub> is the percentage sum of squares explained by the *k*-th IPCA. Simultaneous selection indexes (ssi), are computed by summation of the ranks of the ASV, SIPC, EV and Za indexes and the ranks of the mean yields (E. Farshadfar [2008](#ref-Farshadfar2008)), which results in ssiASV, ssiSIPC, ssiEV, and ssiZa, respectively.

``` r
stab_indexes = AMMI_indexes(AMMI_model)
kable(stab_indexes, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
Code
</th>
<th style="text-align:right;">
Y
</th>
<th style="text-align:right;">
rY
</th>
<th style="text-align:right;">
ASV
</th>
<th style="text-align:right;">
rASV
</th>
<th style="text-align:right;">
ssiASV
</th>
<th style="text-align:right;">
SIPC
</th>
<th style="text-align:right;">
rSIPC
</th>
<th style="text-align:right;">
ssiSIPC
</th>
<th style="text-align:right;">
EV
</th>
<th style="text-align:right;">
rEV
</th>
<th style="text-align:right;">
ssiEV
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
G1
</td>
<td style="text-align:right;">
2.624
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.3378
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.8731
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0222
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
G10
</td>
<td style="text-align:right;">
2.506
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
1.3139
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
2.3288
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
0.1734
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
20
</td>
</tr>
<tr>
<td style="text-align:left;">
G2
</td>
<td style="text-align:right;">
2.703
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.1557
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1.6045
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.1680
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
13
</td>
</tr>
<tr>
<td style="text-align:left;">
G3
</td>
<td style="text-align:right;">
2.941
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.1693
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1.0824
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0486
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
G4
</td>
<td style="text-align:right;">
2.697
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.6776
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
1.8304
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0.1232
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
11
</td>
</tr>
<tr>
<td style="text-align:left;">
G5
</td>
<td style="text-align:right;">
2.566
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0.4212
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
1.0560
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.0374
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:left;">
G6
</td>
<td style="text-align:right;">
2.549
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.2596
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
1.2745
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.0787
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:left;">
G7
</td>
<td style="text-align:right;">
2.710
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1.0162
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
1.8831
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.1618
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
11
</td>
</tr>
<tr>
<td style="text-align:left;">
G8
</td>
<td style="text-align:right;">
3.039
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.5326
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
1.2419
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0586
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
G9
</td>
<td style="text-align:right;">
2.574
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
1.0446
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
2.1135
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
0.1281
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
14
</td>
</tr>
</tbody>
</table>
Estimating the WAASB (based on SVD of BLUP-interaction effects)
===============================================================

The `WAASB()` function computes the weighted average of absolute scores considering all possible IPCA from the Singular Value Decomposition of the genotype-vs-environment interaction random effects (BLUPs), considering the following equation:

$$
        WAASB\_i  = 
        \\sum\_{k = 1}^{p} |IPCA\_{ik} \\times EP\_k|/ \\sum\_{k = 1}^{p}EP\_k
$$

where *W**A**A**S**B*<sub>*i*</sub> is the weighted average of absolute scores of the *i*th genotype; *I**P**C**A*<sub>*i**k*</sub> is the scores of the *i*th genotype in the *k*th IPCA; and *E**P*<sub>*k*</sub> is the explained variance of the *k*th PCA for *k* = 1, 2, ..,*p*, *p* = *m**i**n*(*G* − 1; *E* − 1).

``` r
# Assuming equal weights for productivity and stability (default)
WAASB = WAASB(dataset,
              resp = GY,
              gen = GEN,
              env = ENV,
              rep = BLOCK)
```

Diagnostic plot for residuals
-----------------------------

The function `autoplot()` is used to generate diagnostic plots of residuals of the model. The normality of the random effects of genotype and interaction effects may be also obtained by using `type = "re"`.

``` r
autoplot(WAASB)
```

<img src="D:\Desktop\METAAB\README_files/figure-markdown_github/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

``` r
autoplot(WAASB, type = "re")
```

<img src="D:\Desktop\METAAB\README_files/figure-markdown_github/unnamed-chunk-15-2.png" style="display: block; margin: auto;" />

Printing the model outputs
--------------------------

### Likelihood Ratio Tests

``` r
options(digits = 5)
data = WAASB$LRT
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
npar
</th>
<th style="text-align:right;">
logLik
</th>
<th style="text-align:right;">
AIC
</th>
<th style="text-align:right;">
LRT
</th>
<th style="text-align:right;">
Df
</th>
<th style="text-align:right;">
Pr(&gt;Chisq)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Complete
</td>
<td style="text-align:right;">
51
</td>
<td style="text-align:right;">
-260.38
</td>
<td style="text-align:right;">
622.77
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
Genotype
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
-269.04
</td>
<td style="text-align:right;">
638.08
</td>
<td style="text-align:right;">
17.305
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
3e-05
</td>
</tr>
<tr>
<td style="text-align:left;">
Gen vs Env
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
-287.89
</td>
<td style="text-align:right;">
675.78
</td>
<td style="text-align:right;">
55.005
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0e+00
</td>
</tr>
</tbody>
</table>
The output `LRT` contains the Likelihood Ratio Tests for genotype and genotype-vs-environment random effects.

### Variance components and genetic parameters

``` r
options(digits = 7)
data = WAASB$ESTIMATES
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)
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
0.191119576855307
</td>
</tr>
<tr>
<td style="text-align:left;">
Heritability
</td>
<td style="text-align:left;">
0.126366948791287
</td>
</tr>
<tr>
<td style="text-align:left;">
GEIr2
</td>
<td style="text-align:left;">
0.333599634106794
</td>
</tr>
<tr>
<td style="text-align:left;">
Heribatility of means
</td>
<td style="text-align:left;">
0.797430712385577
</td>
</tr>
<tr>
<td style="text-align:left;">
Accuracy
</td>
<td style="text-align:left;">
0.892989760515526
</td>
</tr>
<tr>
<td style="text-align:left;">
rge
</td>
<td style="text-align:left;">
0.381853266248619
</td>
</tr>
<tr>
<td style="text-align:left;">
CVg
</td>
<td style="text-align:left;">
5.77536600368121
</td>
</tr>
<tr>
<td style="text-align:left;">
CVr
</td>
<td style="text-align:left;">
11.9391409605518
</td>
</tr>
<tr>
<td style="text-align:left;">
CV ratio
</td>
<td style="text-align:left;">
0.483733798165516
</td>
</tr>
</tbody>
</table>
In the output `ESTIMATES`, beyond the variance components for the declared random effects, some important parameters are also shown. **Heribatility** is the broad-sense heritability, h<sub>*g*</sub><sup>2</sup>, estimated by $\\mathop h\\nolimits\_g^2 = \\mathop {\\hat\\sigma} \\nolimits\_g^2 /\\left( {\\mathop {\\hat\\sigma} \\nolimits\_g^2 + \\mathop {\\hat\\sigma} \\nolimits\_i^2 + \\mathop {\\hat\\sigma} \\nolimits\_e^2 } \\right)$ where $\\mathop {\\hat\\sigma} \\nolimits\_g^2$ is the genotypic variance; $\\mathop {\\hat\\sigma} \\nolimits\_i^2$ is the genotype-by-environment interaction variance; and $\\mathop {\\hat\\sigma} \\nolimits\_e^2$ is the residual variance. **GEIr2** is the coefficient of determination of the interaction effects, r<sub>*i*</sub><sup>2</sup>, estimated by $\\mathop r\\nolimits\_i^2 = \\mathop {\\hat\\sigma} \\nolimits\_i^2 /\\left( {\\mathop {\\hat\\sigma} \\nolimits\_g^2 + \\mathop {\\hat\\sigma} \\nolimits\_i^2 + \\mathop {\\hat\\sigma} \\nolimits\_e^2 } \\right)$; **Heribatility of means** is the heribability on the mean basis, h<sub>*g**m*</sub><sup>2</sup>, estimated by $\\mathop h\\nolimits\_{gm}^2 = \\mathop {\\hat\\sigma} \\nolimits\_g^2 /\\left\[ {\\mathop {\\hat\\sigma} \\nolimits\_g^2 + \\mathop {\\hat\\sigma} \\nolimits\_i^2 /e + \\mathop {\\hat\\sigma} \\nolimits\_e^2 /\\left( {eb} \\right)} \\right\]$, where *e* and *b* are the number of environments and blocks, respectively; **Accuracy** is the accuracy of selection, *Ac*, estimated by $Ac = \\sqrt{\\mathop h\\nolimits\_{gm}^2}$ ; **rge** is the genotype-environment correlation, r<sub>*g**e*</sub>, estimated by $\\mathop r\\nolimits\_{ge} = \\mathop {\\hat\\sigma} \\nolimits\_g^2 /\\left({\\mathop {\\hat\\sigma} \\nolimits\_g^2 + \\mathop {\\hat\\sigma} \\nolimits\_i^2} \\right)$; **CVg** is the the genotypic coefficient of variation, estimated by $\\left( {\\sqrt {\\mathop {\\hat \\sigma }\\nolimits\_g^2 } /\\mu } \\right) \\times 100$ , where *μ* is the grand mean; **CVr** is the residual coefficient of variation, estimated by $\\left( {\\sqrt {\\mathop {\\hat \\sigma }\\nolimits\_e^2 } /\\mu } \\right) \\times 100$ ; **CV ratio** is the ratio between genotypic and residual coefficient of variation.

### Some useful information

``` r
options(digits = 4)
data = WAASB$Details
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)
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
0.899 (Genotype G10 in NF2014 )
</td>
</tr>
<tr>
<td style="text-align:left;">
Max
</td>
<td style="text-align:left;">
4.812 (Genotype G8 in WF2016 )
</td>
</tr>
<tr>
<td style="text-align:left;">
MinENV
</td>
<td style="text-align:left;">
Environment NF2014 (1.3682)
</td>
</tr>
<tr>
<td style="text-align:left;">
MaxENV
</td>
<td style="text-align:left;">
Environment WF2012 (4.0643)
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
The following pieces of information are provided in `Details` output. **WgtResponse** is the weight for the response variable in estimating WAASB; **WgtWAAS** is the weight for stability; **Ngen** is the number of genotypes; **Nenv** is the number of environments; **OVmean** is the overall mean; **Min** is the minimum value observed (returning the genotype and environment); **Max** is the maximum observed; **MinENV** is the environment with the lower mean; **MaxENV** is the environment with the largest mean observed; **MinGEN** is the genotype with the lower mean; **MaxGEN** is the genotype with the largest mean.

### The WAASB object

``` r
options(digits = 4)
data = WAASB$model[, c(1:3,13:17, 21:22)]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)
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
WAASBY
</th>
<th style="text-align:right;">
OrWAASBY
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
22.043
</td>
<td style="text-align:right;">
94.44
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
58.24
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
G10
</td>
<td style="text-align:right;">
2.506
</td>
<td style="text-align:right;">
0.4496
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.00
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
36.888
</td>
<td style="text-align:right;">
76.63
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
56.76
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
G3
</td>
<td style="text-align:right;">
2.941
</td>
<td style="text-align:right;">
0.1415
</td>
<td style="text-align:right;">
81.560
</td>
<td style="text-align:right;">
100.00
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
90.78
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
G4
</td>
<td style="text-align:right;">
2.697
</td>
<td style="text-align:right;">
0.2996
</td>
<td style="text-align:right;">
35.802
</td>
<td style="text-align:right;">
48.68
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
42.24
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
G5
</td>
<td style="text-align:right;">
2.566
</td>
<td style="text-align:right;">
0.1837
</td>
<td style="text-align:right;">
11.152
</td>
<td style="text-align:right;">
86.30
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
48.73
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
G6
</td>
<td style="text-align:right;">
2.549
</td>
<td style="text-align:right;">
0.1627
</td>
<td style="text-align:right;">
8.108
</td>
<td style="text-align:right;">
93.11
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
50.61
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
G7
</td>
<td style="text-align:right;">
2.710
</td>
<td style="text-align:right;">
0.2994
</td>
<td style="text-align:right;">
38.256
</td>
<td style="text-align:right;">
48.73
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
43.49
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
G8
</td>
<td style="text-align:right;">
3.039
</td>
<td style="text-align:right;">
0.2063
</td>
<td style="text-align:right;">
100.000
</td>
<td style="text-align:right;">
78.95
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
89.48
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
G9
</td>
<td style="text-align:right;">
2.574
</td>
<td style="text-align:right;">
0.4017
</td>
<td style="text-align:right;">
12.676
</td>
<td style="text-align:right;">
15.54
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
14.11
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
NF2010
</td>
<td style="text-align:right;">
1.989
</td>
<td style="text-align:right;">
0.1974
</td>
<td style="text-align:right;">
23.021
</td>
<td style="text-align:right;">
60.28
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
41.65
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
NF2011
</td>
<td style="text-align:right;">
2.536
</td>
<td style="text-align:right;">
0.1660
</td>
<td style="text-align:right;">
43.329
</td>
<td style="text-align:right;">
71.70
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
57.52
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
NF2012
</td>
<td style="text-align:right;">
3.057
</td>
<td style="text-align:right;">
0.3628
</td>
<td style="text-align:right;">
62.623
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
31.31
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
NF2013
</td>
<td style="text-align:right;">
2.175
</td>
<td style="text-align:right;">
0.1673
</td>
<td style="text-align:right;">
29.930
</td>
<td style="text-align:right;">
71.25
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
50.59
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
NF2014
</td>
<td style="text-align:right;">
1.368
</td>
<td style="text-align:right;">
0.0985
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
96.32
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
48.16
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
NF2015
</td>
<td style="text-align:right;">
1.609
</td>
<td style="text-align:right;">
0.1776
</td>
<td style="text-align:right;">
8.912
</td>
<td style="text-align:right;">
67.48
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
38.20
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
NF2016
</td>
<td style="text-align:right;">
2.910
</td>
<td style="text-align:right;">
0.3004
</td>
<td style="text-align:right;">
57.170
</td>
<td style="text-align:right;">
22.76
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
39.96
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:left;">
ENV
</td>
<td style="text-align:left;">
NF2017
</td>
<td style="text-align:right;">
1.782
</td>
<td style="text-align:right;">
0.1091
</td>
<td style="text-align:right;">
15.344
</td>
<td style="text-align:right;">
92.45
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
53.90
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
WF2010
</td>
<td style="text-align:right;">
2.521
</td>
<td style="text-align:right;">
0.1750
</td>
<td style="text-align:right;">
42.746
</td>
<td style="text-align:right;">
68.44
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
55.59
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
WF2011
</td>
<td style="text-align:right;">
3.180
</td>
<td style="text-align:right;">
0.1435
</td>
<td style="text-align:right;">
67.197
</td>
<td style="text-align:right;">
79.91
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
73.55
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
WF2012
</td>
<td style="text-align:right;">
4.064
</td>
<td style="text-align:right;">
0.2494
</td>
<td style="text-align:right;">
100.000
</td>
<td style="text-align:right;">
41.33
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
70.66
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
WF2013
</td>
<td style="text-align:right;">
3.675
</td>
<td style="text-align:right;">
0.2617
</td>
<td style="text-align:right;">
85.572
</td>
<td style="text-align:right;">
36.86
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
61.22
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
WF2014
</td>
<td style="text-align:right;">
2.507
</td>
<td style="text-align:right;">
0.2304
</td>
<td style="text-align:right;">
42.241
</td>
<td style="text-align:right;">
48.24
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
45.24
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
WF2015
</td>
<td style="text-align:right;">
3.107
</td>
<td style="text-align:right;">
0.2362
</td>
<td style="text-align:right;">
64.497
</td>
<td style="text-align:right;">
46.13
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
55.32
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
WF2016
</td>
<td style="text-align:right;">
3.910
</td>
<td style="text-align:right;">
0.2109
</td>
<td style="text-align:right;">
94.296
</td>
<td style="text-align:right;">
55.35
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
74.82
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
WF2017
</td>
<td style="text-align:right;">
2.663
</td>
<td style="text-align:right;">
0.0884
</td>
<td style="text-align:right;">
48.031
</td>
<td style="text-align:right;">
100.00
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
74.02
</td>
<td style="text-align:right;">
2
</td>
</tr>
</tbody>
</table>
This output generated by the `WAASB` function is very similar to those shown in the [sections 5.1](#assuming-a-given-probability-error-for-chosing-the-number-of-axes) and [5.2](#declaring-a-specific-number-of-axis-to-be-used). The main difference here, is that the singular value decomposition is based on the BLUP interaction effect matrix and the **WAASB** in this output is the weighted average of the absolute scores, estimated with all estimated PCA axes, where instead **WAAS** that is estimated considering only PCA axes with *P*-value ≤ 0.05.

### BLUP for genotypes

``` r
options(digits = 4)
data = WAASB$BLUPgen[1:10,]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)
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
### Plotting the BLUP for genotypes

``` r
# No file exported
p1 = plot.blup(WAASB)
p2 = plot.blup(WAASB, 
               col.shape  =  c("gray20", "gray80")) + coord_flip()
plot_grid(p1, p2,
          labels = c("p1", "p2"))
```

<img src="D:\Desktop\METAAB\README_files/figure-markdown_github/unnamed-chunk-21-1.png" style="display: block; margin: auto;" />

This output shows the predicted means for genotypes. **BLUPg** is the genotypic effect $(\\hat{g}\_{i})$ estimated by $\\hat{g}\_{i} = h\_g^2(\\bar{y}\_{i.}-\\bar{y}\_{..})$ where *h*<sub>*g*</sub><sup>2</sup> is the shrinkage effect for genotype. **Predicted** is the predicted mean estimated by $\\hat{g}\_{i}+\\mu$ where *μ* is the grand mean. **LL** and **UL** are the lower and upper limits, respectively, estimated by $(\\hat{g}\_{i}+\\mu)\\pm{CI}$. *C**I* is the confidence interval for BLUP prediction assuming a given probability error, where $CI = t\\times\\sqrt{((1-Ac)\\times{\\mathop \\sigma \\nolimits\_g^2)}}$ where *t* is the Student's *t* value for a two-tailed t test at a given probability error; *A**c* is the accuracy of selection and σ<sub>*g*</sub><sup>2</sup> is the genotypic variance.

### BLUP for genotypes X environment combination

``` r
options(digits = 4)
data = WAASB$BLUPgge[1:10,]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)
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
NF2010
</td>
<td style="text-align:left;">
G1
</td>
<td style="text-align:right;">
-0.0244
</td>
<td style="text-align:right;">
-0.0536
</td>
<td style="text-align:right;">
-0.0780
</td>
<td style="text-align:right;">
1.911
</td>
<td style="text-align:right;">
1.810
</td>
<td style="text-align:right;">
2.012
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G10
</td>
<td style="text-align:right;">
0.2612
</td>
<td style="text-align:right;">
-0.1473
</td>
<td style="text-align:right;">
0.1138
</td>
<td style="text-align:right;">
2.103
</td>
<td style="text-align:right;">
2.002
</td>
<td style="text-align:right;">
2.204
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G2
</td>
<td style="text-align:right;">
-0.0070
</td>
<td style="text-align:right;">
0.0095
</td>
<td style="text-align:right;">
0.0025
</td>
<td style="text-align:right;">
1.991
</td>
<td style="text-align:right;">
1.891
</td>
<td style="text-align:right;">
2.092
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G3
</td>
<td style="text-align:right;">
-0.0193
</td>
<td style="text-align:right;">
0.1994
</td>
<td style="text-align:right;">
0.1802
</td>
<td style="text-align:right;">
2.169
</td>
<td style="text-align:right;">
2.068
</td>
<td style="text-align:right;">
2.270
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G4
</td>
<td style="text-align:right;">
-0.0083
</td>
<td style="text-align:right;">
0.0049
</td>
<td style="text-align:right;">
-0.0034
</td>
<td style="text-align:right;">
1.986
</td>
<td style="text-align:right;">
1.885
</td>
<td style="text-align:right;">
2.086
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G5
</td>
<td style="text-align:right;">
-0.1509
</td>
<td style="text-align:right;">
-0.0999
</td>
<td style="text-align:right;">
-0.2508
</td>
<td style="text-align:right;">
1.738
</td>
<td style="text-align:right;">
1.637
</td>
<td style="text-align:right;">
1.839
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G6
</td>
<td style="text-align:right;">
-0.0769
</td>
<td style="text-align:right;">
-0.1128
</td>
<td style="text-align:right;">
-0.1897
</td>
<td style="text-align:right;">
1.799
</td>
<td style="text-align:right;">
1.698
</td>
<td style="text-align:right;">
1.900
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G7
</td>
<td style="text-align:right;">
0.3551
</td>
<td style="text-align:right;">
0.0153
</td>
<td style="text-align:right;">
0.3705
</td>
<td style="text-align:right;">
2.359
</td>
<td style="text-align:right;">
2.259
</td>
<td style="text-align:right;">
2.460
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G8
</td>
<td style="text-align:right;">
-0.0033
</td>
<td style="text-align:right;">
0.2778
</td>
<td style="text-align:right;">
0.2745
</td>
<td style="text-align:right;">
2.263
</td>
<td style="text-align:right;">
2.163
</td>
<td style="text-align:right;">
2.364
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G9
</td>
<td style="text-align:right;">
-0.3262
</td>
<td style="text-align:right;">
-0.0934
</td>
<td style="text-align:right;">
-0.4196
</td>
<td style="text-align:right;">
1.569
</td>
<td style="text-align:right;">
1.468
</td>
<td style="text-align:right;">
1.670
</td>
</tr>
</tbody>
</table>
This output shows the predicted means for each genotype and environment combination. **BLUPg** is the genotypic effect described above. **BLUPge** is the genotypic effect of the *i*th genotype in the *j*th environment $(\\hat{g}\_{ij})$ estimated by $\\hat{g}\_{ij} = h\_g^2(\\bar{y}\_{i.}-\\bar{y}\_{..})+h\_{ge}^2(y\_{ij}-\\bar{y}\_{i.}-\\bar{y}\_{.j}+\\bar{y}\_{..})$, where *h*<sub>*g**e*</sub><sup>2</sup> is the shrinkage effect for the genotype-by-environment interaction; **BLUPg+ge** is *B**L**U**P*<sub>*g*</sub> + *B**L**U**P*<sub>*g**e*</sub>; **Predicted** is the predicted mean ($\\hat{y}\_{ij}$) estimated by $\\hat{y}\_{ij} = \\bar{y}\_{.j}+BLUP\_{g+ge}$.

Eigenvalues from the SVD
------------------------

``` r
options(digits = 4)
data = WAASB$PCA
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)
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
``` r
plot.eigen(WAASB, size.lab = 14, size.tex = 14)
```

<img src="D:\Desktop\METAAB\README_files/figure-markdown_github/unnamed-chunk-24-1.png" style="display: block; margin: auto;" />

The above output shows the eigenvalues and the proportion of variance explained by each principal component axis of the BLUP interaction effects matrix.

### Phenotypic means

``` r
options(digits = 4)
data = WAASB$MeansGxE[1:10,]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)
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
NF2010
</td>
<td style="text-align:left;">
G1
</td>
<td style="text-align:right;">
1.898
</td>
<td style="text-align:right;">
0.1339
</td>
<td style="text-align:right;">
0.1097
</td>
<td style="text-align:right;">
2.638
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G10
</td>
<td style="text-align:right;">
2.244
</td>
<td style="text-align:right;">
0.1339
</td>
<td style="text-align:right;">
-0.7340
</td>
<td style="text-align:right;">
2.408
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G2
</td>
<td style="text-align:right;">
1.988
</td>
<td style="text-align:right;">
0.1339
</td>
<td style="text-align:right;">
0.0630
</td>
<td style="text-align:right;">
2.711
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G3
</td>
<td style="text-align:right;">
2.159
</td>
<td style="text-align:right;">
0.1339
</td>
<td style="text-align:right;">
-0.0830
</td>
<td style="text-align:right;">
2.930
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G4
</td>
<td style="text-align:right;">
1.981
</td>
<td style="text-align:right;">
0.1339
</td>
<td style="text-align:right;">
0.2908
</td>
<td style="text-align:right;">
2.736
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G5
</td>
<td style="text-align:right;">
1.657
</td>
<td style="text-align:right;">
0.1339
</td>
<td style="text-align:right;">
0.0776
</td>
<td style="text-align:right;">
2.576
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G6
</td>
<td style="text-align:right;">
1.758
</td>
<td style="text-align:right;">
0.1339
</td>
<td style="text-align:right;">
0.1195
</td>
<td style="text-align:right;">
2.565
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G7
</td>
<td style="text-align:right;">
2.551
</td>
<td style="text-align:right;">
0.1339
</td>
<td style="text-align:right;">
0.6690
</td>
<td style="text-align:right;">
2.800
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G8
</td>
<td style="text-align:right;">
2.262
</td>
<td style="text-align:right;">
0.1339
</td>
<td style="text-align:right;">
-0.0203
</td>
<td style="text-align:right;">
3.037
</td>
</tr>
<tr>
<td style="text-align:left;">
NF2010
</td>
<td style="text-align:left;">
G9
</td>
<td style="text-align:right;">
1.393
</td>
<td style="text-align:right;">
0.1339
</td>
<td style="text-align:right;">
-0.4923
</td>
<td style="text-align:right;">
2.508
</td>
</tr>
</tbody>
</table>
In this output, *Y* is the phenotypic mean for each genotype and environment combination (*y*<sub>*i**j*</sub>), estimated by *y*<sub>*i**j*</sub> = ∑<sub>*k*</sub>*y*<sub>*i**j*</sub>/*B* with *k* = 1, 2, ...*B*.

Biplots
-------

We will show how biplots may be obtained for both traditional AMMI model, fitted by the function `WAAS.AMMI()` and the mixed-effect model fitted by the function `WAASB()`. Provided that an object of class "WAAS.AMMI" or "WAASB" is available in the global environment, the graphics may be obtained using the function `plot.scores()`. To do that, we will revisit the previusly fitted model `WAASB` . Please, refer to `?plot.scores` for more details. Four types of graphics can be generated: 1 = *P**C*1 × *P**C*2; 2 = *G**Y* × *P**C*1; 3 = *G**Y* × *W**A**A**S**B*; and 4 = a graphic with nominal yield as a function of the environment PCA1 scores.

### biplot type 1: PC1 x PC2

``` r
library(cowplot)
p1 = plot.scores(WAASB, type = 1)
p2 = plot.scores(WAASB,
                 type = 1,
                 polygon = TRUE,
                 col.gen = "black",
                 col.env = "gray70",
                 col.segm.env = "gray70",
                 axis.expand = 1.5)
plot_grid(p1, p2, labels = c("p1","p2"))
```

<img src="D:\Desktop\METAAB\README_files/figure-markdown_github/unnamed-chunk-26-1.png" style="display: block; margin: auto;" />

### biplot type 2: GY x PC1

``` r
p3 = plot.scores(WAASB, type = 2)
p4 = plot.scores(WAASB, type = 2,
                 col.segm.env = "transparent") +
                 theme_gray() +
                 theme(legend.position = c(0.1, 0.9),
                       legend.background = element_rect(fill = NA))

plot_grid(p3, p4, labels = c("p3","p4"))
```

<img src="D:\Desktop\METAAB\README_files/figure-markdown_github/unnamed-chunk-27-1.png" style="display: block; margin: auto;" />

### biplot type 3: GY x WAASB

The quadrants proposed in the following biplot represent the four classifications proposed here regarding the joint interpretation of productivity and stability. The genotypes or environments included in quadrant I can be considered unstable genotypes or environments with high discrimination ability, and with productivity below the grand mean. In quadrant II are included unstable genotypes, although with productivity above the grand mean. The environments included in this quadrant deserve special attention since, in addition to providing high magnitudes of the response variable, they present a good discrimination ability. Genotypes within quadrant III have low productivity, but can be considered stable due to the lower values of WAASB. The lower this value, the more stable the genotype can be considered. The environments included in this quadrant can be considered as poorly productive and with low discrimination ability. The genotypes within the quadrant IV are higly productive and broadly adapted due to the high magnitude of the response variable and high stability performance (lower values of WAASB).

``` r
p5 = plot.scores(WAASB, type = 3)
p6 = plot.scores(WAASB, type = 3,
                 x.lab = "My customized x label",
                 size.shape = 3,
                 size.tex = 2,
                 x.lim = c(1.2, 4.7),
                 x.breaks = seq(1.5, 4.5, by = 0.5)) + 
                 theme(legend.position = c(0.1, 0.9))
plot_grid(p5, p6, labels = c("p5","p6"))
```

<img src="D:\Desktop\METAAB\README_files/figure-markdown_github/unnamed-chunk-28-1.png" style="display: block; margin: auto;" />

### biplot type 4 : nominal yield and environment PCA1

``` r
plot.scores(WAASB,
            type = 4, size.tex = 1.5)
```

<img src="D:\Desktop\METAAB\README_files/figure-markdown_github/unnamed-chunk-29-1.png" style="display: block; margin: auto;" />

Estimating the WAASBY index
===========================

The function `WAASBratio` considers both stability (WAASB) and productivity for genotype ranking considering the following model:

$$
WAASB{Y\_i} = \\frac{{\\left( {r{G\_i} \\times {\\theta \_Y}} \\right) + \\left( {r{W\_i} \\times {\\theta \_S}} \\right)}}{{{\\theta \_Y} + {\\theta \_S}}}
$$

where *W**A**A**S**B**Y*<sub>*i*</sub> is the superiority index for the *i*-th genotype that weights between performance and stability; *r**G*<sub>*i*</sub> and *r**W*<sub>*i*</sub> are the rescaled values (0-100) for GY and WAASB, respectively; *θ*<sub>*Y*</sub> and *θ*<sub>*S*</sub> are the weights for GY and WAASB, respectively.

This function provides the option of attributing weights for stability and productivity in genotype ranking. This is important depending on the goal of a selection strategy or cultivar reccomendation. For example, if the goal is to select a genotype with high yield (independently of the stability performance), that genotype with the first rank in a WAASB/GY = 0/100 ratio should be selected. The reciprocal is true. Aiming at selecting high-stable genotypes (independently of the productivity), that genotype with the first rank in a WAASB/GY = 100/0 ratio should be selected. By default, the increment on the WAASB/GY ratio is equal to 5 and the WAASBY values are saved when the WAASB/GY ratio is equal to 50/50. Thus, twenty one different scenarios are computed, and for each scenario, the genotypes are ranked based on the WAASBY values.

In the following example we will assume that we want to obtain the ranks changing the WAASB/GY ratio in 10% each scenario and to plot the WAASBY values considering a WAASB/GY ratio equal 30/70. **Important! THE ARGUMENT `saveWAASY` MUST BE DIVISIBLE BY `increment`!**

``` r
WAASBYratio = WAASBYratio(dataset,
                          resp = GY,
                          gen = GEN,
                          env = ENV,
                          rep = BLOCK,
                          increment = 10,
                          saveWAASY = 30)
```

This procedure can also be used with the traditional AMMI analysis. This approach is easily implemented using the `WAASratio.AMMI` function shown in the following example.

``` r
WAASBYratio2 = WAASratio.AMMI(dataset,
                              resp = GY,
                              gen = GEN,
                              env = ENV,
                              rep = BLOCK,
                              increment = 10,
                              saveWAASY = 30)
```

Printing the model outputs
--------------------------

-   WAASBY values

``` r
options(digits = 4)
data = WAASBYratio$WAASY
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)
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
70
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
0.00
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
70
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
13.54
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
70
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
33.61
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
70
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
33.70
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
70
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
39.66
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
70
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
41.40
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
70
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
43.76
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
70
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
48.81
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
70
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
87.09
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
70
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
93.69
</td>
<td style="text-align:left;">
above
</td>
</tr>
</tbody>
</table>
In this output, **PesResp** and **PesWAAS** are the weights attributed to response variable and stability, respecively, to compute de WAASY values.

-   Genotype ranking for each scenario of WAASBY/GY ratio.

``` r
options(digits = 4)
data = WAASBYratio$hetcomb
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)
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
8
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
</tr>
<tr>
<td style="text-align:left;">
G5
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
6
</td>
<td style="text-align:right;">
6
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
</tr>
<tr>
<td style="text-align:left;">
G6
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
5
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
6
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
</tr>
<tr>
<td style="text-align:left;">
G7
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
4
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
9
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
7
</td>
</tr>
</tbody>
</table>
-   Genotype ranking depending on the number of multiplicative terms used to estimate the WAASB index.

``` r
options(digits = 4)
data = WAASBYratio$hetdata
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)
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
8
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
Plotting the WAASBY values
--------------------------

``` r
p1 = plot.WAASBY(WAASBYratio)
p2 = plot.WAASBY(WAASBYratio, col.shape = c("gray20", "gray80")) 
plot_grid(p1, p2, labels = c("p1", "p2"))
```

<img src="D:\Desktop\METAAB\README_files/figure-markdown_github/unnamed-chunk-35-1.png" style="display: block; margin: auto;" />

Plotting the heat map graphics
------------------------------

The first type of heatmap shows the genotype ranking depending on the number of principal component axes used for estimating the WAASB index. An euclidean distance-based dendrogram is used for grouping the genotype ranking for both genotypes and principal component axes. The second type of heatmap shows the genotype ranking depending on the WAASB/GY ratio. The ranks obtained with a ratio of 100/0 considers exclusively the stability for genotype ranking. On the other hand, a ratio of 0/100 considers exclusively the productivity for genotype ranking. Four clusters are estimated (1) unproductive and unstable genotypes; (2) productive, but unstable genotypes; (3) stable, but unproductive genotypes; and (4), productive and stable genotypes.

### Ranks of genotypes depending on the number of PCA used to estimate the WAASB

``` r
plot(WAASBYratio, type = 1)
```

<img src="D:\Desktop\METAAB\README_files/figure-markdown_github/unnamed-chunk-36-1.png" style="display: block; margin: auto;" />

### Ranks of genotypes depending on the WAASB/GY ratio

``` r
plot(WAASBYratio, type = 2)
```

<img src="D:\Desktop\METAAB\README_files/figure-markdown_github/unnamed-chunk-37-1.png" style="display: block; margin: auto;" />

Stability indexes proposed by Resende 2007
==========================================

``` r
res_inde = Resende_indexes(WAASB)
kable(res_inde, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
HMGV
</th>
<th style="text-align:right;">
HMGV\_order
</th>
<th style="text-align:right;">
RPGV
</th>
<th style="text-align:right;">
GY\_RPGV
</th>
<th style="text-align:right;">
RPGV\_order
</th>
<th style="text-align:right;">
HMRPGV
</th>
<th style="text-align:right;">
GY\_HMRPGV
</th>
<th style="text-align:right;">
HMRPGV\_order
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
G1
</td>
<td style="text-align:right;">
2.359
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.9706
</td>
<td style="text-align:right;">
2.612
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.9683
</td>
<td style="text-align:right;">
2.606
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
2.177
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.9229
</td>
<td style="text-align:right;">
2.483
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.9068
</td>
<td style="text-align:right;">
2.440
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
2.460
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
1.0062
</td>
<td style="text-align:right;">
2.707
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.9986
</td>
<td style="text-align:right;">
2.687
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
G3
</td>
<td style="text-align:right;">
2.695
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1.0931
</td>
<td style="text-align:right;">
2.941
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1.0914
</td>
<td style="text-align:right;">
2.937
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
2.450
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
1.0024
</td>
<td style="text-align:right;">
2.697
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.9990
</td>
<td style="text-align:right;">
2.688
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
G5
</td>
<td style="text-align:right;">
2.349
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.9582
</td>
<td style="text-align:right;">
2.579
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.9560
</td>
<td style="text-align:right;">
2.573
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
G6
</td>
<td style="text-align:right;">
2.342
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0.9539
</td>
<td style="text-align:right;">
2.567
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0.9519
</td>
<td style="text-align:right;">
2.562
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
G7
</td>
<td style="text-align:right;">
2.502
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1.0173
</td>
<td style="text-align:right;">
2.737
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
1.0085
</td>
<td style="text-align:right;">
2.714
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
2.793
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1.1301
</td>
<td style="text-align:right;">
3.041
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1.1267
</td>
<td style="text-align:right;">
3.032
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
2.260
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.9453
</td>
<td style="text-align:right;">
2.544
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.9358
</td>
<td style="text-align:right;">
2.518
</td>
<td style="text-align:right;">
9
</td>
</tr>
</tbody>
</table>
References
==========

Farshadfar, E. 2008. “Incorporation of AMMI stability value and grain yield in a single non-parametric index (GSI) in bread wheat.” *Pakistan Journal of Biological Sciences* 11 (14): 1791–6. <http://www.ncbi.nlm.nih.gov/pubmed/18817218>.

Gauch, H.G, and Zobel R.W. 1988. “Predictive and Postdictive Success of Statistical Analyses of Yield Trials.” *Theor. Appl. Genet.* 76 (1): 1–10. doi:[10.1007/BF00288824](https://doi.org/10.1007/BF00288824).

Gollob, H. F. 1968. “A statistical model which combines features of factor analytic and analysis of variance techniques.” *Psychometrika* 33 (1): 73–115. doi:[10.1007/BF02289676](https://doi.org/10.1007/BF02289676).

Olivoto, T. 2018. “WAASBdata.” *Mendeley Data* V1. doi:[10.17632/2sjz32k3s3.1](https://doi.org/10.17632/2sjz32k3s3.1).

Piepho, H.P. 1994. “Best Linear Unbiased Prediction (Blup) for Regional Yield Trials: A Comparison to Additive Main Effects and Multiplicative Interaction (Ammi) Analysis.” *Theor. Appl. Genet.* 89 (5): 647–54. doi:[10.1007/BF00222462](https://doi.org/10.1007/BF00222462).

Purchase, J. L., Hesta Hatting, and C. S. van Deventer. 2000. “Genotype × environment interaction of winter wheat ( Triticum aestivum L.) in South Africa: II. Stability analysis of yield performance.” *South African Journal of Plant and Soil* 17 (3). Taylor & Francis Group: 101–7. doi:[10.1080/02571862.2000.10634878](https://doi.org/10.1080/02571862.2000.10634878).

Sneller, C. H., L. Kilgore-Norquest, and D. Dombek. 1997. “Repeatability of Yield Stability Statistics in Soybean.” *Crop Science* 37 (2). Crop Science Society of America: 383–90. doi:[10.2135/cropsci1997.0011183X003700020013x](https://doi.org/10.2135/cropsci1997.0011183X003700020013x).

Zali, H., E. Farshadfar, S.H. Sabaghpour, and R. Karimizadeh. 2012. “Evaluation of genotype × environment interaction in chickpea using measures of stability from AMMI model.” *Annals of Biological Research* 3 (7): 3126–36. <http://eprints.icrisat.ac.in/7173/>.
