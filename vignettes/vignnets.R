## ---- eval = FALSE-------------------------------------------------------
#  # download the package from Github
#  devtools::install_github("TiagoOlivoto/WAASB")
#  

## ----warning = FALSE, message = FALSE------------------------------------
# Importing data
require(WAASB)
dataset = read.csv("https://data.mendeley.com/datasets/2sjz32k3s3/1/files/07764a07-172a-4285-85db-c31bc39ae480/WAASBdata.csv?dl=1")


## ----eval = TRUE---------------------------------------------------------
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


## ----echo = TRUE---------------------------------------------------------
options(digits = 4)
RMSEweat = rbind(AMMIweat$RMSEmean, BLUPweat$RMSEmean)
RMSEweat = dplyr::mutate(RMSEweat, CROP = "Wheat")
RMSEweat = RMSEweat[order(RMSEweat[,2], decreasing = F),]
#print(RMSEweat)

## ----echo = FALSE, warning = FALSE---------------------------------------
library(kableExtra)
options(digits = 4)
kable(RMSEweat, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", full_width = F, position = "left", font_size = 12)


## ----eval = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center"----
# binding AMMI and BLUP RMSEs
RMSEweat = list(RMSE = rbind(AMMIweat$RMSE, 
                           BLUPweat$RMSE))
# Plotting the RMSE values
plot.validation.AMMIF(RMSEweat,
                      violin = FALSE,
                      col.boxplot = "gray75")

## ----echo = TRUE---------------------------------------------------------
# Predicting the yield of 10 oat cultivars in 16 environments using 5 multiplicative terms.
predictoat = predict.AMMI(dataset,
                          resp = "GY",
                          naxis = 5)


## ----echo = FALSE--------------------------------------------------------
library(kableExtra)
options(digits = 4)
predictoat = predictoat[1:10,]
kable(predictoat, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", full_width = F, position = "float_left", font_size = 12)


## ----echo = TRUE---------------------------------------------------------
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


## ----eval = FALSE--------------------------------------------------------
#  options (digits = 4)
#  # printing the WAASB object
#  print(WAAS1$anova)
#  

## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(kableExtra)
options(digits = 4, width = 200)
data = WAAS1$anova
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12) %>%
  row_spec(9, bold = T) %>%
  add_indent(c(5:13))


## ----eval = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  options (digits = 4)
#  # printing the WAASB object
#  print(WAAS1$model[, c(1:3,13:17, 21:22)])
#  

## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(kableExtra)
options(digits = 4)
data = WAAS1$model[, c(1:3,13:17, 21:22)]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


