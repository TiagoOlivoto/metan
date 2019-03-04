---
title: "Analyzing multienvironment trials using AMMI"
always_allow_html: yes
output:
  prettydoc::html_pretty:
    theme: cayman
    highlight: vignette
    readme: TRUE
    toc: yes
    number_sections: yes
    keep_md: TRUE
fig_caption: yes
link-citations: true
bibliography: METAABref.bib    
vignette: >
  %\VignetteIndexEntry{Analyzing multienvironment trials using AMMI}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


# Getting started

In this section, we will use the data in `data_ge`. For more information, please, see `?data_ge`. Other data sets can be used provided that the following columns are in the dataset: environment, genotype, block/replicate and response variable(s).


```r
library(METAAB)
library(ggplot2)
library(cowplot) # used in this material to arrange the graphics
library(kableExtra) # Used to 
dataset = data_ge
str(dataset)
```

```
## Classes 'tbl_df', 'tbl' and 'data.frame':	420 obs. of  5 variables:
##  $ ENV: Factor w/ 14 levels "E1","E10","E11",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ GEN: Factor w/ 10 levels "G1","G10","G2",..: 1 1 1 3 3 3 4 4 4 5 ...
##  $ REP: Factor w/ 3 levels "1","2","3": 1 2 3 1 2 3 1 2 3 1 ...
##  $ GY : num  2.17 2.5 2.43 3.21 2.93 ...
##  $ HM : num  44.9 46.9 47.8 45.2 45.3 ...
```

## The AMMI model
The estimate of the response for the *i*th genotype in the *j*th environment using The Additive Main Effect and Multiplicative interaction (AMMI) model, is given as follows:

$$
{y_{ij}} = \mu  + {\alpha_i} + {\tau_j} + \sum\limits_{k = 1}^p {{\lambda _k}{a_{ik}}} {t_{jk}} + {\rho _{ij}} + {\varepsilon _{ij}}
$$
where ${\lambda_k}$ is the singular value for the *k*-th interaction principal component axis (IPCA); $a_{ik}$ is the *i*-th element of the *k*-th eigenvector; $t_{jk}$ is the *j*th element of the *k*th eigenvector. A residual $\rho _{ij}$ remains, if not all *p* IPCA are used, where $p \le min(g - 1; e - 1)$.
 
The AMMI model is fitted with the `WAAS.AMMI()` function. The first argument is the data, in our example `dataset`. The second argument (`resp`) is the response variable to be analyzed. The function allow a single variable (in this case GY) or a vector of response variables. The arguments (`gen`, `env`, and `rep`) are the name of the columns that contais the levels for genotypes, environments, and replications, respectively. The last argument (`verbose`) control if the code is run silently or not.


```r
AMMI_model = WAAS.AMMI(data = dataset,
                       resp = GY,
                       gen = GEN,
                       env = ENV,
                       rep = REP,
                       verbose = FALSE)
```

## Within-environment ANOVA

A within-environment ANOVA considering a fixed-effect model is computed. For each environment the Mean Squares for block, genotypes and error are shown. Estimated F-value and the probability error are also shown for block and genotype effects. Some measures of experimental precision are calculated, namelly, coefficient of variation, $CV = (\sqrt{MS_{res}}/Mean) \times 100$; the heritability, $h2 = (MS_{gen} - MS_{res})/MS_{gen}$, and the accuracy of selection, $As = \sqrt{h2}$.



```r
# printing the WAAS object
options(digits = 3)
data = AMMI_model$GY$individual$individual
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                full_width = F, position = "left", font_size = 12)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
 <thead>
  <tr>
   <th style="text-align:left;"> ENV </th>
   <th style="text-align:right;"> Mean </th>
   <th style="text-align:right;"> MSblock </th>
   <th style="text-align:right;"> MSgen </th>
   <th style="text-align:right;"> MSres </th>
   <th style="text-align:right;"> Fcal(Blo) </th>
   <th style="text-align:right;"> Pr&gt;F(Blo) </th>
   <th style="text-align:right;"> Fcal(Gen) </th>
   <th style="text-align:right;"> Pr&gt;F(Gen) </th>
   <th style="text-align:right;"> CV(%) </th>
   <th style="text-align:right;"> h2 </th>
   <th style="text-align:right;"> AS </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> E1 </td>
   <td style="text-align:right;"> 2.52 </td>
   <td style="text-align:right;"> 0.065 </td>
   <td style="text-align:right;"> 0.337 </td>
   <td style="text-align:right;"> 0.144 </td>
   <td style="text-align:right;"> 0.453 </td>
   <td style="text-align:right;"> 0.643 </td>
   <td style="text-align:right;"> 2.34 </td>
   <td style="text-align:right;"> 0.059 </td>
   <td style="text-align:right;"> 15.05 </td>
   <td style="text-align:right;"> 0.573 </td>
   <td style="text-align:right;"> 0.757 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E10 </td>
   <td style="text-align:right;"> 2.17 </td>
   <td style="text-align:right;"> 0.654 </td>
   <td style="text-align:right;"> 0.296 </td>
   <td style="text-align:right;"> 0.027 </td>
   <td style="text-align:right;"> 24.508 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 11.09 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 7.51 </td>
   <td style="text-align:right;"> 0.910 </td>
   <td style="text-align:right;"> 0.954 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E11 </td>
   <td style="text-align:right;"> 1.37 </td>
   <td style="text-align:right;"> 0.377 </td>
   <td style="text-align:right;"> 0.151 </td>
   <td style="text-align:right;"> 0.105 </td>
   <td style="text-align:right;"> 3.594 </td>
   <td style="text-align:right;"> 0.049 </td>
   <td style="text-align:right;"> 1.44 </td>
   <td style="text-align:right;"> 0.244 </td>
   <td style="text-align:right;"> 23.68 </td>
   <td style="text-align:right;"> 0.304 </td>
   <td style="text-align:right;"> 0.552 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E12 </td>
   <td style="text-align:right;"> 1.61 </td>
   <td style="text-align:right;"> 0.092 </td>
   <td style="text-align:right;"> 0.320 </td>
   <td style="text-align:right;"> 0.054 </td>
   <td style="text-align:right;"> 1.717 </td>
   <td style="text-align:right;"> 0.208 </td>
   <td style="text-align:right;"> 5.98 </td>
   <td style="text-align:right;"> 0.001 </td>
   <td style="text-align:right;"> 14.38 </td>
   <td style="text-align:right;"> 0.833 </td>
   <td style="text-align:right;"> 0.913 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E13 </td>
   <td style="text-align:right;"> 2.91 </td>
   <td style="text-align:right;"> 0.077 </td>
   <td style="text-align:right;"> 0.713 </td>
   <td style="text-align:right;"> 0.099 </td>
   <td style="text-align:right;"> 0.772 </td>
   <td style="text-align:right;"> 0.477 </td>
   <td style="text-align:right;"> 7.18 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 10.83 </td>
   <td style="text-align:right;"> 0.861 </td>
   <td style="text-align:right;"> 0.928 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E14 </td>
   <td style="text-align:right;"> 1.78 </td>
   <td style="text-align:right;"> 0.104 </td>
   <td style="text-align:right;"> 0.131 </td>
   <td style="text-align:right;"> 0.075 </td>
   <td style="text-align:right;"> 1.374 </td>
   <td style="text-align:right;"> 0.278 </td>
   <td style="text-align:right;"> 1.73 </td>
   <td style="text-align:right;"> 0.153 </td>
   <td style="text-align:right;"> 15.40 </td>
   <td style="text-align:right;"> 0.423 </td>
   <td style="text-align:right;"> 0.650 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E2 </td>
   <td style="text-align:right;"> 3.18 </td>
   <td style="text-align:right;"> 0.698 </td>
   <td style="text-align:right;"> 0.207 </td>
   <td style="text-align:right;"> 0.179 </td>
   <td style="text-align:right;"> 3.912 </td>
   <td style="text-align:right;"> 0.039 </td>
   <td style="text-align:right;"> 1.16 </td>
   <td style="text-align:right;"> 0.376 </td>
   <td style="text-align:right;"> 13.29 </td>
   <td style="text-align:right;"> 0.136 </td>
   <td style="text-align:right;"> 0.369 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E3 </td>
   <td style="text-align:right;"> 4.06 </td>
   <td style="text-align:right;"> 0.489 </td>
   <td style="text-align:right;"> 0.335 </td>
   <td style="text-align:right;"> 0.179 </td>
   <td style="text-align:right;"> 2.731 </td>
   <td style="text-align:right;"> 0.092 </td>
   <td style="text-align:right;"> 1.87 </td>
   <td style="text-align:right;"> 0.123 </td>
   <td style="text-align:right;"> 10.41 </td>
   <td style="text-align:right;"> 0.466 </td>
   <td style="text-align:right;"> 0.683 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E4 </td>
   <td style="text-align:right;"> 3.67 </td>
   <td style="text-align:right;"> 0.116 </td>
   <td style="text-align:right;"> 0.531 </td>
   <td style="text-align:right;"> 0.138 </td>
   <td style="text-align:right;"> 0.846 </td>
   <td style="text-align:right;"> 0.446 </td>
   <td style="text-align:right;"> 3.86 </td>
   <td style="text-align:right;"> 0.007 </td>
   <td style="text-align:right;"> 10.10 </td>
   <td style="text-align:right;"> 0.741 </td>
   <td style="text-align:right;"> 0.861 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E5 </td>
   <td style="text-align:right;"> 3.91 </td>
   <td style="text-align:right;"> 0.219 </td>
   <td style="text-align:right;"> 0.526 </td>
   <td style="text-align:right;"> 0.066 </td>
   <td style="text-align:right;"> 3.297 </td>
   <td style="text-align:right;"> 0.060 </td>
   <td style="text-align:right;"> 7.93 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 6.59 </td>
   <td style="text-align:right;"> 0.874 </td>
   <td style="text-align:right;"> 0.935 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E6 </td>
   <td style="text-align:right;"> 2.66 </td>
   <td style="text-align:right;"> 0.160 </td>
   <td style="text-align:right;"> 0.135 </td>
   <td style="text-align:right;"> 0.059 </td>
   <td style="text-align:right;"> 2.729 </td>
   <td style="text-align:right;"> 0.092 </td>
   <td style="text-align:right;"> 2.30 </td>
   <td style="text-align:right;"> 0.063 </td>
   <td style="text-align:right;"> 9.09 </td>
   <td style="text-align:right;"> 0.565 </td>
   <td style="text-align:right;"> 0.752 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E7 </td>
   <td style="text-align:right;"> 1.99 </td>
   <td style="text-align:right;"> 0.381 </td>
   <td style="text-align:right;"> 0.337 </td>
   <td style="text-align:right;"> 0.091 </td>
   <td style="text-align:right;"> 4.185 </td>
   <td style="text-align:right;"> 0.032 </td>
   <td style="text-align:right;"> 3.70 </td>
   <td style="text-align:right;"> 0.009 </td>
   <td style="text-align:right;"> 15.17 </td>
   <td style="text-align:right;"> 0.730 </td>
   <td style="text-align:right;"> 0.854 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E8 </td>
   <td style="text-align:right;"> 2.54 </td>
   <td style="text-align:right;"> 0.817 </td>
   <td style="text-align:right;"> 0.215 </td>
   <td style="text-align:right;"> 0.028 </td>
   <td style="text-align:right;"> 29.369 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 7.72 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 6.58 </td>
   <td style="text-align:right;"> 0.870 </td>
   <td style="text-align:right;"> 0.933 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E9 </td>
   <td style="text-align:right;"> 3.06 </td>
   <td style="text-align:right;"> 0.583 </td>
   <td style="text-align:right;"> 0.679 </td>
   <td style="text-align:right;"> 0.111 </td>
   <td style="text-align:right;"> 5.253 </td>
   <td style="text-align:right;"> 0.016 </td>
   <td style="text-align:right;"> 6.12 </td>
   <td style="text-align:right;"> 0.001 </td>
   <td style="text-align:right;"> 10.90 </td>
   <td style="text-align:right;"> 0.837 </td>
   <td style="text-align:right;"> 0.915 </td>
  </tr>
</tbody>
</table>

## The AMMI table



```r
data = AMMI_model$GY$anova
kable(data,  align  = "l", booktabs = T, format = "html", linesep = "") %>%
kable_styling(bootstrap_options = "striped", "condensed",
              position = "left", full_width = F, font_size = 12) %>%
row_spec(5:8, bold = T) %>%
add_indent(c(5:13))
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;"> Df </th>
   <th style="text-align:left;"> Sum Sq </th>
   <th style="text-align:left;"> Mean Sq </th>
   <th style="text-align:left;"> F value </th>
   <th style="text-align:left;"> Pr(&gt;F) </th>
   <th style="text-align:left;"> Percent </th>
   <th style="text-align:left;"> Accumul </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ENV </td>
   <td style="text-align:left;"> 13 </td>
   <td style="text-align:left;"> 279.574 </td>
   <td style="text-align:left;"> 21.506 </td>
   <td style="text-align:left;"> 62.33 </td>
   <td style="text-align:left;"> 0.000 </td>
   <td style="text-align:left;"> . </td>
   <td style="text-align:left;"> . </td>
  </tr>
  <tr>
   <td style="text-align:left;"> REP(ENV) </td>
   <td style="text-align:left;"> 28 </td>
   <td style="text-align:left;"> 9.662 </td>
   <td style="text-align:left;"> 0.345 </td>
   <td style="text-align:left;"> 3.57 </td>
   <td style="text-align:left;"> 0.000 </td>
   <td style="text-align:left;"> . </td>
   <td style="text-align:left;"> . </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GEN </td>
   <td style="text-align:left;"> 9 </td>
   <td style="text-align:left;"> 12.995 </td>
   <td style="text-align:left;"> 1.444 </td>
   <td style="text-align:left;"> 14.93 </td>
   <td style="text-align:left;"> 0.000 </td>
   <td style="text-align:left;"> . </td>
   <td style="text-align:left;"> . </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENV:GEN </td>
   <td style="text-align:left;"> 117 </td>
   <td style="text-align:left;"> 31.220 </td>
   <td style="text-align:left;"> 0.267 </td>
   <td style="text-align:left;"> 2.76 </td>
   <td style="text-align:left;"> 0.000 </td>
   <td style="text-align:left;"> . </td>
   <td style="text-align:left;"> . </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold; padding-left: 2em;" indentlevel="1"> PC1 </td>
   <td style="text-align:left;font-weight: bold;"> 21 </td>
   <td style="text-align:left;font-weight: bold;"> 10.749 </td>
   <td style="text-align:left;font-weight: bold;"> 0.512 </td>
   <td style="text-align:left;font-weight: bold;"> 5.29 </td>
   <td style="text-align:left;font-weight: bold;"> 0.000 </td>
   <td style="text-align:left;font-weight: bold;"> 34.4 </td>
   <td style="text-align:left;font-weight: bold;"> 34.4 </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold; padding-left: 2em;" indentlevel="1"> PC2 </td>
   <td style="text-align:left;font-weight: bold;"> 19 </td>
   <td style="text-align:left;font-weight: bold;"> 9.924 </td>
   <td style="text-align:left;font-weight: bold;"> 0.522 </td>
   <td style="text-align:left;font-weight: bold;"> 5.40 </td>
   <td style="text-align:left;font-weight: bold;"> 0.000 </td>
   <td style="text-align:left;font-weight: bold;"> 31.8 </td>
   <td style="text-align:left;font-weight: bold;"> 66.2 </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold; padding-left: 2em;" indentlevel="1"> PC3 </td>
   <td style="text-align:left;font-weight: bold;"> 17 </td>
   <td style="text-align:left;font-weight: bold;"> 4.039 </td>
   <td style="text-align:left;font-weight: bold;"> 0.238 </td>
   <td style="text-align:left;font-weight: bold;"> 2.46 </td>
   <td style="text-align:left;font-weight: bold;"> 0.001 </td>
   <td style="text-align:left;font-weight: bold;"> 12.9 </td>
   <td style="text-align:left;font-weight: bold;"> 79.2 </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold; padding-left: 2em;" indentlevel="1"> PC4 </td>
   <td style="text-align:left;font-weight: bold;"> 15 </td>
   <td style="text-align:left;font-weight: bold;"> 3.074 </td>
   <td style="text-align:left;font-weight: bold;"> 0.205 </td>
   <td style="text-align:left;font-weight: bold;"> 2.12 </td>
   <td style="text-align:left;font-weight: bold;"> 0.010 </td>
   <td style="text-align:left;font-weight: bold;"> 9.8 </td>
   <td style="text-align:left;font-weight: bold;"> 89 </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> PC5 </td>
   <td style="text-align:left;"> 13 </td>
   <td style="text-align:left;"> 1.446 </td>
   <td style="text-align:left;"> 0.111 </td>
   <td style="text-align:left;"> 1.15 </td>
   <td style="text-align:left;"> 0.318 </td>
   <td style="text-align:left;"> 4.6 </td>
   <td style="text-align:left;"> 93.6 </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> PC6 </td>
   <td style="text-align:left;"> 11 </td>
   <td style="text-align:left;"> 0.932 </td>
   <td style="text-align:left;"> 0.085 </td>
   <td style="text-align:left;"> 0.88 </td>
   <td style="text-align:left;"> 0.561 </td>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:left;"> 96.6 </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> PC7 </td>
   <td style="text-align:left;"> 9 </td>
   <td style="text-align:left;"> 0.567 </td>
   <td style="text-align:left;"> 0.063 </td>
   <td style="text-align:left;"> 0.65 </td>
   <td style="text-align:left;"> 0.754 </td>
   <td style="text-align:left;"> 1.8 </td>
   <td style="text-align:left;"> 98.4 </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> PC8 </td>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:left;"> 0.362 </td>
   <td style="text-align:left;"> 0.052 </td>
   <td style="text-align:left;"> 0.54 </td>
   <td style="text-align:left;"> 0.804 </td>
   <td style="text-align:left;"> 1.2 </td>
   <td style="text-align:left;"> 99.6 </td>
  </tr>
  <tr>
   <td style="text-align:left; padding-left: 2em;" indentlevel="1"> PC9 </td>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> 0.126 </td>
   <td style="text-align:left;"> 0.025 </td>
   <td style="text-align:left;"> 0.26 </td>
   <td style="text-align:left;"> 0.934 </td>
   <td style="text-align:left;"> 0.4 </td>
   <td style="text-align:left;"> 100 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:left;"> 252 </td>
   <td style="text-align:left;"> 24.367 </td>
   <td style="text-align:left;"> 0.097 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> . </td>
   <td style="text-align:left;"> . </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Total </td>
   <td style="text-align:left;"> 419 </td>
   <td style="text-align:left;"> 357.816 </td>
   <td style="text-align:left;"> 0.854 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> . </td>
   <td style="text-align:left;"> . </td>
  </tr>
</tbody>
</table>


Nine interaction principal component axis (IPCA) were fitted and four were significant at 5% probability error. Based on this result, the AMMI4 model would be the best model to predict the yielding of the genotypes in the studied environments.


# Estimating the response variable based on significant IPCA axes

An interesting feature of `METAAB` is the implementation of the S3 method `predict()`. The response variable of a two-way table (for example, the yield of *m* genotypes in *n* environments) may be estimated using the function `predict(object)`, where `object` is an object of class `WAAS.AMMI`. This estimation is based on the number of multiplicative terms declared in the function. If `naxis = 0` is declared, only the main effects (AMMI0) are considered. In this case, the estimated mean will be the estimate from OLS estimation. If `naxis = 1`, the AMMI1 (with one multiplicative term) is used for estimating the response variable. If `naxis = min(g-1; e-1)`, the AMMIF is fitted. A summary of all possible AMMI models is presented below.

| Member of AMMI family  | Espected response of the *i*-th genotype in the *j*th environment|
|:------------------------|:------------------------------------------------------------------------------|
| AMMI0            | $\hat{y}_{ij} = \bar{y}_{i.} + \bar{y}_{.j} - \bar{y}_{..}$                   |
| AMMI1            |$\hat{y}_{ij} = \bar{y}_{i.} + \bar{y}_{.j} - \bar{y}_{..} +\lambda_1 a_{i1}t_{j1}$ |
| AMMI2            |$\hat{y}_{ij} = \bar{y}_{i.} + \bar{y}_{.j} - \bar{y}_{..} +\lambda_1 a_{i1}t_{j1}+\lambda_2 a_{i2}t_{j2}$ |
| ...              |                                                                               |
| AMMIF            |$\hat{y}_{ij} = \bar{y}_{i.} + \bar{y}_{.j} - \bar{y}_{..} +\lambda_1 a_{i1}t_{j1}+\lambda_2 a_{i2}t_{j2}+...+\lambda_p a_{ip}t_{jp}$ 


Procedures based on postdictive success, such as Gollobs's test [@Gollob:1968] or predictive success, such as cross-validation [@Piepho:1994] should be used to define the number of IPCA used for estimating the response variable in AMMI analysis. This package provides both. The `WAAS.AMMI()` function compute traditional AMMI analysis showing the number of significant axes according to Gollobs's test. On the other hand, `validation.AMMIF()` function provides cross-validation of AMMI-model family, considering a completely randomized design (CRD) or a randomized complete block design (RCBD).



```r
predicted = predict(AMMI_model, naxis = 4)
predicted = predicted$GY[1:5,]
kable(predicted, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", full_width = T)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> ENV </th>
   <th style="text-align:left;"> GEN </th>
   <th style="text-align:right;"> Y </th>
   <th style="text-align:right;"> resOLS </th>
   <th style="text-align:right;"> Ypred </th>
   <th style="text-align:right;"> ResAMMI </th>
   <th style="text-align:right;"> YpredAMMI </th>
   <th style="text-align:right;"> AMMI0 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> E1 </td>
   <td style="text-align:left;"> G1 </td>
   <td style="text-align:right;"> 2.37 </td>
   <td style="text-align:right;"> -0.084 </td>
   <td style="text-align:right;"> 2.45 </td>
   <td style="text-align:right;"> 0.0712 </td>
   <td style="text-align:right;"> 2.52 </td>
   <td style="text-align:right;"> 2.45 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E1 </td>
   <td style="text-align:left;"> G10 </td>
   <td style="text-align:right;"> 1.97 </td>
   <td style="text-align:right;"> -0.344 </td>
   <td style="text-align:right;"> 2.32 </td>
   <td style="text-align:right;"> -0.3539 </td>
   <td style="text-align:right;"> 1.96 </td>
   <td style="text-align:right;"> 2.32 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E1 </td>
   <td style="text-align:left;"> G2 </td>
   <td style="text-align:right;"> 2.90 </td>
   <td style="text-align:right;"> 0.311 </td>
   <td style="text-align:right;"> 2.59 </td>
   <td style="text-align:right;"> 0.2904 </td>
   <td style="text-align:right;"> 2.88 </td>
   <td style="text-align:right;"> 2.59 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E1 </td>
   <td style="text-align:left;"> G3 </td>
   <td style="text-align:right;"> 2.89 </td>
   <td style="text-align:right;"> 0.087 </td>
   <td style="text-align:right;"> 2.80 </td>
   <td style="text-align:right;"> -0.0452 </td>
   <td style="text-align:right;"> 2.76 </td>
   <td style="text-align:right;"> 2.80 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> E1 </td>
   <td style="text-align:left;"> G4 </td>
   <td style="text-align:right;"> 2.59 </td>
   <td style="text-align:right;"> 0.100 </td>
   <td style="text-align:right;"> 2.49 </td>
   <td style="text-align:right;"> 0.0494 </td>
   <td style="text-align:right;"> 2.54 </td>
   <td style="text-align:right;"> 2.49 </td>
  </tr>
</tbody>
</table>


Only the first five values are shown. The following values are presented: **ENV** is the environment; **GEN** is the genotype; **Y** is the response variable; **resOLS** is the residual ($\hat{z}_{ij}$) estimated by the Ordinary Least Square (OLS), where $\hat{z}_{ij} = y_{ij} - \bar{y}_{i.} - \bar{y}_{.j} + \bar{y}_{ij}$; **Ypred** is the predicted value by OLS ($\hat{y}_{ij} = y_{ij} -\hat{z}_{ij}$); **ResAMMI** is the residual estimated by the AMMI model ($\hat{a}_{ij}$) considering the number of multiplicative terms informed in the function (in this case 5), where $\hat{a}_{ij} = \lambda_1\alpha_{i1}\tau_{j1}+...+\lambda_5\alpha_{i5}\tau_{j5}$; **YpredAMMI** is the predicted value by AMMI model  $\hat{ya}_{ij} = \bar{y}_{i.} + \bar{y}_{.j} - \bar{y}_{ij}+\hat{a}_{ij}$; and **AMMI0** is the predicted value when no multiplicative terms are used, i.e., $\hat{y}_{ij} = \bar{y}_{i.} + \bar{y}_{.j} - \bar{y}_{ij}$.


# Estimating the WAAS index
The `WAAS.AMMI()` function computes the Weighted Average of Absolute Scores considering (i) all principal component axes that were significant ($p < 0.05$ by default); or (ii) declaring a specific number of axes to be used, according to the following equation:

$$
        WAAS_i  = 
        \sum_{k = 1}^{p} |IPCA_{ik} \times EP_k|/ \sum_{k = 1}^{p}EP_k
$$

where $WAAS_i$ is the weighted average of absolute scores of the *i*th genotype; $PCA_{ik}$ is the score of the *i*th genotype in the *k*th IPCA; and $EP_k$ is the explained variance of the *k*th IPCA for $k = 1,2,..,p$, considering *p* the number of significant PCAs, or a declared number of PCAs. The following functions may be used to do that.

## Number of axes based on F-test
In this example only IPCAs with *P*-value < 0.05 will be considered in the WAAS estimation. This is the default setting and the model was already fitted and stored into AMMI_model>GY>model



```r
data = AMMI_model$GY$model[, c(1:3,13:17, 21:22)]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
 <thead>
  <tr>
   <th style="text-align:left;"> type </th>
   <th style="text-align:left;"> Code </th>
   <th style="text-align:right;"> Y </th>
   <th style="text-align:right;"> WAAS </th>
   <th style="text-align:right;"> PctResp </th>
   <th style="text-align:right;"> PctWAAS </th>
   <th style="text-align:right;"> OrResp </th>
   <th style="text-align:right;"> OrWAAS </th>
   <th style="text-align:right;"> WAASY </th>
   <th style="text-align:right;"> OrWAASY </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> GEN </td>
   <td style="text-align:left;"> G1 </td>
   <td style="text-align:right;"> 2.60 </td>
   <td style="text-align:right;"> 0.126 </td>
   <td style="text-align:right;"> 24.88 </td>
   <td style="text-align:right;"> 100.0 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 62.4 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GEN </td>
   <td style="text-align:left;"> G10 </td>
   <td style="text-align:right;"> 2.47 </td>
   <td style="text-align:right;"> 0.551 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.0 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.0 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GEN </td>
   <td style="text-align:left;"> G2 </td>
   <td style="text-align:right;"> 2.74 </td>
   <td style="text-align:right;"> 0.365 </td>
   <td style="text-align:right;"> 51.26 </td>
   <td style="text-align:right;"> 43.8 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 47.5 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GEN </td>
   <td style="text-align:left;"> G3 </td>
   <td style="text-align:right;"> 2.96 </td>
   <td style="text-align:right;"> 0.131 </td>
   <td style="text-align:right;"> 90.93 </td>
   <td style="text-align:right;"> 99.0 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 95.0 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GEN </td>
   <td style="text-align:left;"> G4 </td>
   <td style="text-align:right;"> 2.64 </td>
   <td style="text-align:right;"> 0.272 </td>
   <td style="text-align:right;"> 32.06 </td>
   <td style="text-align:right;"> 65.7 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 48.9 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GEN </td>
   <td style="text-align:left;"> G5 </td>
   <td style="text-align:right;"> 2.54 </td>
   <td style="text-align:right;"> 0.255 </td>
   <td style="text-align:right;"> 12.42 </td>
   <td style="text-align:right;"> 69.6 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 41.0 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GEN </td>
   <td style="text-align:left;"> G6 </td>
   <td style="text-align:right;"> 2.53 </td>
   <td style="text-align:right;"> 0.275 </td>
   <td style="text-align:right;"> 11.80 </td>
   <td style="text-align:right;"> 65.0 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 38.4 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GEN </td>
   <td style="text-align:left;"> G7 </td>
   <td style="text-align:right;"> 2.74 </td>
   <td style="text-align:right;"> 0.434 </td>
   <td style="text-align:right;"> 50.65 </td>
   <td style="text-align:right;"> 27.4 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 39.0 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GEN </td>
   <td style="text-align:left;"> G8 </td>
   <td style="text-align:right;"> 3.00 </td>
   <td style="text-align:right;"> 0.307 </td>
   <td style="text-align:right;"> 100.00 </td>
   <td style="text-align:right;"> 57.5 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 78.8 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GEN </td>
   <td style="text-align:left;"> G9 </td>
   <td style="text-align:right;"> 2.51 </td>
   <td style="text-align:right;"> 0.401 </td>
   <td style="text-align:right;"> 7.32 </td>
   <td style="text-align:right;"> 35.4 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 21.3 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENV </td>
   <td style="text-align:left;"> E1 </td>
   <td style="text-align:right;"> 2.52 </td>
   <td style="text-align:right;"> 0.199 </td>
   <td style="text-align:right;"> 42.74 </td>
   <td style="text-align:right;"> 73.1 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 57.9 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENV </td>
   <td style="text-align:left;"> E10 </td>
   <td style="text-align:right;"> 2.17 </td>
   <td style="text-align:right;"> 0.193 </td>
   <td style="text-align:right;"> 29.93 </td>
   <td style="text-align:right;"> 74.6 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 52.3 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENV </td>
   <td style="text-align:left;"> E11 </td>
   <td style="text-align:right;"> 1.37 </td>
   <td style="text-align:right;"> 0.151 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 85.3 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 42.7 </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENV </td>
   <td style="text-align:left;"> E12 </td>
   <td style="text-align:right;"> 1.61 </td>
   <td style="text-align:right;"> 0.192 </td>
   <td style="text-align:right;"> 8.91 </td>
   <td style="text-align:right;"> 75.1 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 42.0 </td>
   <td style="text-align:right;"> 12 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENV </td>
   <td style="text-align:left;"> E13 </td>
   <td style="text-align:right;"> 2.91 </td>
   <td style="text-align:right;"> 0.322 </td>
   <td style="text-align:right;"> 57.17 </td>
   <td style="text-align:right;"> 42.1 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 49.6 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENV </td>
   <td style="text-align:left;"> E14 </td>
   <td style="text-align:right;"> 1.78 </td>
   <td style="text-align:right;"> 0.205 </td>
   <td style="text-align:right;"> 15.34 </td>
   <td style="text-align:right;"> 71.7 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 43.5 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENV </td>
   <td style="text-align:left;"> E2 </td>
   <td style="text-align:right;"> 3.18 </td>
   <td style="text-align:right;"> 0.293 </td>
   <td style="text-align:right;"> 67.20 </td>
   <td style="text-align:right;"> 49.4 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 58.3 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENV </td>
   <td style="text-align:left;"> E3 </td>
   <td style="text-align:right;"> 4.06 </td>
   <td style="text-align:right;"> 0.310 </td>
   <td style="text-align:right;"> 100.00 </td>
   <td style="text-align:right;"> 45.2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 72.6 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENV </td>
   <td style="text-align:left;"> E4 </td>
   <td style="text-align:right;"> 3.67 </td>
   <td style="text-align:right;"> 0.345 </td>
   <td style="text-align:right;"> 85.57 </td>
   <td style="text-align:right;"> 36.4 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 61.0 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENV </td>
   <td style="text-align:left;"> E5 </td>
   <td style="text-align:right;"> 3.91 </td>
   <td style="text-align:right;"> 0.255 </td>
   <td style="text-align:right;"> 94.29 </td>
   <td style="text-align:right;"> 59.0 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 76.6 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENV </td>
   <td style="text-align:left;"> E6 </td>
   <td style="text-align:right;"> 2.66 </td>
   <td style="text-align:right;"> 0.093 </td>
   <td style="text-align:right;"> 48.03 </td>
   <td style="text-align:right;"> 100.0 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 74.0 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENV </td>
   <td style="text-align:left;"> E7 </td>
   <td style="text-align:right;"> 1.99 </td>
   <td style="text-align:right;"> 0.300 </td>
   <td style="text-align:right;"> 23.02 </td>
   <td style="text-align:right;"> 47.6 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 35.3 </td>
   <td style="text-align:right;"> 13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENV </td>
   <td style="text-align:left;"> E8 </td>
   <td style="text-align:right;"> 2.54 </td>
   <td style="text-align:right;"> 0.305 </td>
   <td style="text-align:right;"> 43.33 </td>
   <td style="text-align:right;"> 46.5 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 44.9 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ENV </td>
   <td style="text-align:left;"> E9 </td>
   <td style="text-align:right;"> 3.06 </td>
   <td style="text-align:right;"> 0.489 </td>
   <td style="text-align:right;"> 62.62 </td>
   <td style="text-align:right;"> 0.0 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 31.3 </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
</tbody>
</table>

In this example, the scores of the nine PCA were not shown. The output generated by the `WAAS.AMMI()` function shows the following results: **type**, genotype (GEN) or environment (ENV); **Code**, the code attributed to each level of the factors; **Y**, the response variable (in this case the grain yield); **WAAS** the weighted average of the absolute scores, estimated with all PCA axes with *P*-value $\le$ 0.05; **PctWAAS** and **PctResp** that are the percentage values for the WAAS and Y, respectively; **OrResp** and **OrWAAS** that are the ranks attributed to the genotype and environment regarding the Y or WAAS, respectively; **WAASY** is the weighted average of absolute scores and response variable. In this case, considering equal weights for PctResp and PctWAAS, the WAASY for G1 is estimated by: $WAAS_{G1} = [(86.32\times50)+(98.88\times50)]/50+50 = 92.60$. Then the **OrWAASY* is the rank for the WAASY value. The genotype (or environment) with the largest WAASY value has the first ranked. See [Estimating the WAASBY index](#estimating-the-waasby-index) for a detailed explanation.



## Number of axes declared manually
The second option to compute the WAAS is by manually declaring a specific number of multiplicative terms. In this case, the number of terms declared is used independently of its significance. Let us, for the moment, assume that after a cross-validation procedure the AMMI7 was the most predictively accurate AMMI model and the researcher will use this model. The additional argument `naxis` in the function `WAAS.AMMI` is then used to overwrite the default chose of significant terms. 


```r
AMMI_model_2 = WAAS.AMMI(dataset,
                         resp = GY,
                         gen = GEN,
                         env = ENV,
                         rep = REP,
                         naxis = 7, # Use 7 IPCA for computing WAAS
                         verbose = FALSE)
```

The only difference in this output is that here we declared that seven IPCA axes should be used for computing the WAAS value. Thus, only the values of WAAS, OrWAAS, WAASY and OrWAASY may have significant changes.


## Biplots

Provided that an object of class "WAAS.AMMI" is available in the global environment, the graphics may be obtained using the function `plot.scores()`. To do that, we will revisit the previusly fitted model `AMMI_model` . Please, refer to `?plot.scores` for more details. Four types of graphics can be generated: 1 = $PC1 \times PC2$;  2 = $GY \times PC1$; 3 = $GY \times WAASB$; and 4 = a graphic with nominal yield as a function of the environment PCA1 scores.

### biplot type 1: PC1 x PC2


```r
library(cowplot)
p1 = plot.scores(AMMI_model$GY, type = 1)
p2 = plot.scores(AMMI_model$GY,
                 type = 1,
                 polygon = TRUE,
                 col.gen = "black",
                 col.env = "gray70",
                 col.segm.env = "gray70",
                 axis.expand = 1.5)
plot_grid(p1, p2, labels = c("p1","p2"))
```

<img src="vignettes_AMMI_files/figure-html/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />


### biplot type 2: GY x PC1


```r
p3 = plot.scores(AMMI_model$GY, type = 2)
p4 = plot.scores(AMMI_model$GY, type = 2,
                 col.segm.env = "transparent") +
                 theme_gray() +
                 theme(legend.position = c(0.1, 0.9),
                       legend.background = element_rect(fill = NA))

plot_grid(p3, p4, labels = c("p3","p4"))
```

<img src="vignettes_AMMI_files/figure-html/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

### biplot type 3: GY x WAAS

The quadrants in the following biplot represent four classes of genotypes/environments regarding the joint interpretation of mean performance and stability. The genotypes or environments included in quadrant I can be considered unstable genotypes or environments with high discrimination ability, and with productivity below the grand mean. In quadrant II are included unstable genotypes, although with productivity above the grand mean. The environments included in this quadrant deserve special attention since, in addition to providing high magnitudes of the response variable, they present a good discrimination ability. Genotypes within quadrant III have low productivity, but can be considered stable due to the lower values of WAASB. The lower this value, the more stable the genotype can be considered. The environments included in this quadrant can be considered as poorly productive and with low discrimination ability. The genotypes within the quadrant IV are higly productive and broadly adapted due to the high magnitude of the response variable and high stability performance (lower values of WAASB).    


```r
p5 = plot.scores(AMMI_model$GY, type = 3)
p6 = plot.scores(AMMI_model$GY, type = 3,
                 x.lab = "My customized x label",
                 size.shape = 3,
                 size.tex.pa = 2,
                 x.lim = c(1.2, 4.7),
                 x.breaks = seq(1.5, 4.5, by = 0.5)) + 
                 theme(legend.position = c(0.1, 0.9))
plot_grid(p5, p6, labels = c("p5","p6"))
```

<img src="vignettes_AMMI_files/figure-html/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />



### biplot type 4 : nominal yield and environment IPCA1


```r
plot.scores(AMMI_model$GY,
            type = 4, size.tex.pa = 1.5)
```

<img src="vignettes_AMMI_files/figure-html/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />


# Simultaneous selection for mean performance and stability

The WAASY index is used for genotype ranking considering both the stability (WAAS) and mean performance based on the following model:

$$
WAASY{_i} = \frac{{\left( {r{G_i} \times {\theta _Y}} \right) + \left( {r{W_i} \times {\theta _S}} \right)}}{{{\theta _Y} + {\theta _S}}}
$$

where $WAASY_i$ is the superiority index for the *i*-th genotype that weights between performance and stability; $rG_i$ and $rW_i$ are the rescaled values (0-100) for GY and WAASB, respectively;  $\theta _Y$ and $\theta_S$ are the weights for GY and WAASB, respectively.

This index was also already computed and stored into AMMI_model>GY>model. An intuitively plot may be obtained by running


```r
library(ggplot2)
p1 = plot.WAASBY(AMMI_model$GY)
p2 = plot.WAASBY(AMMI_model$GY, col.shape = c("gray20", "gray80"))
plot_grid(p1, p2, labels = c("p1", "p2"))
```

<img src="vignettes_AMMI_files/figure-html/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

The values of WAASY in the plot above were computed considering equal weights for mean performance and stability. Different weights may be assigned using the `wresp` argument of the `WAAS.AMMI()` function.

# Estimating the WAASY in different scenarios

In the following example, we will assume that we want to obtain the ranks considering different scenarios (different weights). Supposing that the WAAS/GY weight ratio is changed by 10% each scenario the following function is used.


```r
WAASratio = WAASratio.AMMI(dataset,
                           resp = GY,
                           gen = GEN,
                           env = ENV,
                           rep = REP,
                           increment = 10)
```


## Printing the model outputs

The genotype ranking for each scenario of WAASY/GY weight ratio is shown bellow


```r
options(digits = 4)
data = WAASratio$hetcomb
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> 100/0 </th>
   <th style="text-align:right;"> 90/10 </th>
   <th style="text-align:right;"> 80/20 </th>
   <th style="text-align:right;"> 70/30 </th>
   <th style="text-align:right;"> 60/40 </th>
   <th style="text-align:right;"> 50/50 </th>
   <th style="text-align:right;"> 40/60 </th>
   <th style="text-align:right;"> 30/70 </th>
   <th style="text-align:right;"> 20/80 </th>
   <th style="text-align:right;"> 10/90 </th>
   <th style="text-align:right;"> 0/100 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> G1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G2 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G3 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G4 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G5 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G6 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G7 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G8 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G9 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
</tbody>
</table>

In addition, the genotype ranking depending on the number of multiplicative terms used to estimate the WAAS index is also computed.


```r
options(digits = 4)
data = WAASratio$hetdata
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)
```

<table class="table table-striped" style="font-size: 12px; width: auto !important; ">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> WAAS </th>
   <th style="text-align:right;"> 9PC </th>
   <th style="text-align:right;"> 8PC </th>
   <th style="text-align:right;"> 7PC </th>
   <th style="text-align:right;"> 6PC </th>
   <th style="text-align:right;"> 5PC </th>
   <th style="text-align:right;"> 4PC </th>
   <th style="text-align:right;"> 3PC </th>
   <th style="text-align:right;"> 2PC </th>
   <th style="text-align:right;"> 1PC </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> G1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G2 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G3 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G4 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G5 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G6 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G7 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G8 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G9 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 8 </td>
  </tr>
</tbody>
</table>


## Plotting the heat map graphics
The first type of heatmap shows the genotype ranking depending on the number of principal component axes used for estimating the WAASB index. An euclidean distance-based dendrogram is used for grouping the genotype ranking for both genotypes and principal component axes. The second type of heatmap shows the genotype ranking depending on the WAASB/GY ratio. The ranks obtained with a ratio of 100/0 considers exclusively the stability for genotype ranking. On the other hand, a ratio of 0/100 considers exclusively the productivity for genotype ranking.


### Ranks of genotypes depending on the number of PCA used to estimate the WAAS

```r
plot(WAASratio, type = 1)
```

<img src="vignettes_AMMI_files/figure-html/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

### Ranks of genotypes depending on the WAAS/GY ratio

```r
plot(WAASratio, type = 2)
```

<img src="vignettes_AMMI_files/figure-html/unnamed-chunk-17-1.png" style="display: block; margin: auto;" />




# Other AMMI-based stability indexes
The following AMMI-based stability indexes may be computed using the function `AMMI_indexes()`:
 
 * **AMMI stability value, ASV, [@Purchase2000].**

$$
ASV = \sqrt {{{\left[ {\frac{{IPCA{1_{ss}}}}{{IPCA{2_{ss}}}} \times \left( {IPCA{1_{score}}} \right)} \right]}^2} + {{\left( {IPCA{2_{score}}} \right)}^2}}
$$

* **Sums of the absolute value of the IPCA scores**

$$
SIP{C_i} = \sum\nolimits_{k = 1}^P {\left| {\mathop {\lambda }\nolimits_k^{0.5} {a_{ik}}} \right|}
$$

* **Averages of the squared eigenvector values**

$$
E{V_i} = \sum\nolimits_{k = 1}^P {\mathop a\nolimits_{ik}^2 } /P
$$
described by @Sneller1997, where *P* is the number of IPCA retained via F-tests;

* **absolute value of the relative contribution of IPCAs to the interaction [@Zali2012].**

$$
Z{a_i} = \sum\nolimits_{k = 1}^P {{\theta _k}{a_{ik}}}
$$

where ${\theta _k}$ is the percentage sum of squares explained by the *k*-th IPCA. Simultaneous selection indexes (ssi), are computed by summation of the ranks of the ASV, SIPC, EV and Za indexes and the ranks of the mean yields [@Farshadfar2008], which results in ssiASV, ssiSIPC, ssiEV, and ssiZa, respectively.

The `AMMI_index()` function has two arguments. The first (x) is the model, which must be an object of the class `WAAS.AMMI`. The second, (order.y) is the order for ranking the response variable. By default, it is set to NULL, which means that the response variable is ordered in descending order. If `x` is a list with more than one variable, `order.y` must be a vector of the same length of x. Each element of the vector must be one of the "h" or "l". If "h" is used, the response variable will be ordered from maximum to minimum. If "l" is used then the response variable will be ordered from minimum to maximum.


```r
stab_indexes = AMMI_indexes(AMMI_model, order.y = NULL)
kable(stab_indexes, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)
```

<table class="kable_wrapper table table-striped" style="font-size: 12px; width: auto !important; ">
<tbody>
  <tr>
   <td> 

<table>
 <thead>
  <tr>
   <th style="text-align:left;"> Code </th>
   <th style="text-align:right;"> Y </th>
   <th style="text-align:right;"> rY </th>
   <th style="text-align:right;"> ASV </th>
   <th style="text-align:right;"> rASV </th>
   <th style="text-align:right;"> ssiASV </th>
   <th style="text-align:right;"> SIPC </th>
   <th style="text-align:right;"> rSIPC </th>
   <th style="text-align:right;"> ssiSIPC </th>
   <th style="text-align:right;"> EV </th>
   <th style="text-align:right;"> rEV </th>
   <th style="text-align:right;"> ssiEV </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> G1 </td>
   <td style="text-align:right;"> 2.604 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.3458 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.4628 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 0.0149 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G10 </td>
   <td style="text-align:right;"> 2.471 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 1.2255 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 2.0685 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 0.2101 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 20 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G2 </td>
   <td style="text-align:right;"> 2.744 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.2493 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 1.5443 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 0.1791 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G3 </td>
   <td style="text-align:right;"> 2.955 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.1131 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.5515 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 0.0207 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G4 </td>
   <td style="text-align:right;"> 2.642 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.5939 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 1.0358 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 0.0521 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G5 </td>
   <td style="text-align:right;"> 2.537 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 0.4304 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.9967 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.0433 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G6 </td>
   <td style="text-align:right;"> 2.534 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 0.2652 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 1.1397 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.0911 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 14 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G7 </td>
   <td style="text-align:right;"> 2.741 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 0.6632 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 1.7873 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.1913 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 13 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G8 </td>
   <td style="text-align:right;"> 3.004 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.5739 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 1.1778 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 0.0669 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G9 </td>
   <td style="text-align:right;"> 2.510 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 0.9827 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 1.4950 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.1306 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
</tbody>
</table>

 </td>
  </tr>
</tbody>
</table>


# References

