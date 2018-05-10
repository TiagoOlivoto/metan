## ---- eval = FALSE, collapse=TRUE, comment = "#"-------------------------
#  # download the package from Github
#  devtools::install_github("TiagoOlivoto/WAASB")
#  

## ----warning = FALSE, message = FALSE------------------------------------
# Importing data
require(WAASB)
dataset = read.csv("https://data.mendeley.com/datasets/2sjz32k3s3/1/files/07764a07-172a-4285-85db-c31bc39ae480/WAASBdata.csv?dl=1")


## ----eval = TRUE, collapse=TRUE, comment = "#"---------------------------
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
#  print(WAAS1$model[, c(1:3,13:17, 21:22)])
#  

## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(kableExtra)
options(digits = 4)
data = WAAS1$model[, c(1:3,13:17, 21:22)]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center"-------------------------------------------------------------------------------------------------------------------------------

plot.scores(WAAS1,
            type = 3)

## ----echo = TRUE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

WAAS2 = WAAS.AMMI(dataset,
                 resp = "GY",
                 naxis = 7,
                 weight.response = 50,
                 weight.WAAS = 50)


## ----eval = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  options (digits = 4)
#  # printing the WAASB object
#  print(WAAS2$model[, c(1:3,13:17, 21:22)])
#  

## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(kableExtra)
options(digits = 4)
data = WAAS2$model[, c(1:3,13:17, 21:22)]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center"-------------------------------------------------------------------------------------------------------------------------------
plot.scores(WAAS2,
            type = 3)

## ----echo = TRUE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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

## ----eval = FALSE, warning=F, message=F---------------------------------------------------------------------------------------------------------------------------------------------------------------
#  print(WAASB$LRT)

## ----echo = FALSE, warning=F, message=F---------------------------------------------------------------------------------------------------------------------------------------------------------------
library(kableExtra)
options (digits = 5)
dt_footnote = cbind(WAASB$LRT, WAASB3$LRT)

names(dt_footnote)[1] <- paste0(names(dt_footnote)[1], 
                                footnote_marker_symbol(1))
names(dt_footnote)[2] <- paste0(names(dt_footnote)[2], 
                                footnote_marker_symbol(2))
names(dt_footnote)[3] <- paste0(names(dt_footnote)[3], 
                                footnote_marker_symbol(3))
names(dt_footnote)[7] <- paste0(names(dt_footnote)[7], 
                                footnote_marker_symbol(4))

kable(dt_footnote, "html", align = "c", escape = F) %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12) %>%
  add_header_above(c(" ", "Genotype LRT" = 2,
                     "Interaction LRT" = 2,
                     "Genotype LRT" = 2,
                     "Environment LRT" = 2,
                     "Interaction LRT" = 2)) %>%
  
  add_header_above(c(" ", "Genotype random effect" = 4,
                     "Genotype and environment random effects" = 6)) %>%
  
  footnote(symbol = c("Reduced model without genotype effect; ",
                      "Complete model; ",
                      "Reduced model without genotype-vs-environment interaction effect; ",
                      "Reduced model without environment effect; "),
           footnote_as_chunk = F)

## ----eval = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  options (digits = 4)
#  print(WAASB$ESTIMATES)
#  

## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(kableExtra)
options(digits = 7)
data = WAASB$ESTIMATES
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ----eval = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  options (digits = 4)
#  print(WAASB$Details)
#  

## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(kableExtra)
options(digits = 4)
data = WAASB$Details
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ----eval = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  options (digits = 4)
#  print(WAASB$model[, c(1:3,13:17, 21:22)])
#  

## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(kableExtra)
options(digits = 4)
data = WAASB$model[, c(1:3,13:17, 21:22)]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ----eval = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  options(digits = 4)
#  print(WAASB$BLUPgen[1:10,])
#  

## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(kableExtra)
options(digits = 4)
data = WAASB$BLUPgen[1:10,]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)

  


## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center"-------------------------------------------------------------------------------------------------------------------------------
# No file exported
plot.blup(WAASB)


## ----eval = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  options (digits = 4)
#  print(WAASB$BLUPgge[1:10,])
#  

## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(kableExtra)
options(digits = 4)
data = WAASB$BLUPgge[1:10,]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ----eval = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  options (digits = 4)
#  
#  print(WAASB$PCA)
#  

## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(kableExtra)
options(digits = 4)
data = WAASB$PCA
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center"-------------------------------------------------------------------------------------------------------------------------------

plot.eigen(WAASB, size.lab = 14, size.tex = 14)


## ----eval = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  options (digits = 4)
#  print(WAASB$MeansGxE[1:10,])

## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(kableExtra)
options(digits = 4)
data = WAASB$MeansGxE[1:10,]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center"-------------------------------------------------------------------------------------------------------------------------------
# No file exported
plot.scores(WAASB,
            type = 1)


## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center"-------------------------------------------------------------------------------------------------------------------------------

plot.scores(WAASB,
            type = 2)

# Save to *.pdf file (default)
plot.scores(WAASB,
            type = 2,
            export = TRUE)


## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center", message=F, warning=F---------------------------------------------------------------------------------------------------------

plot.scores(WAASB,
            type = 3)

# Save to a *.tiff file with resolution of 600 dpi.
plot.scores(WAASB,
            type = 3,
            export = TRUE,
            file.type = "tiff",
            resolution = 600)


## ----echo = TRUE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

WAASBYratio = WAASBYratio(dataset,
                          resp = "GY",
                          increment = 10,
                          saveWAASY = 50)


## ----eval = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  
#  
#  print(WAASBYratio$WAASY)
#  

## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(kableExtra)
options(digits = 4)
data = WAASBYratio$WAASY
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ----eval = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  
#  print(WAASBYratio$hetcomb)
#  

## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(kableExtra)
options(digits = 4)
data = WAASBYratio$hetcomb
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ----eval = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  
#  print(WAASBYratio$hetdata)
#  

## ----echo = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(kableExtra)
options(digits = 4)
data = WAASBYratio$hetdata
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ----echo = TRUE, fig.height = 3, fig.width = 5.5, fig.align = "center"-------------------------------------------------------------------------------------------------------------------------------

plot.WAASBY(WAASBYratio,
            legend.pos = c(0.9, 0.2))


## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center"-------------------------------------------------------------------------------------------------------------------------------

plot.WAASBYratio(WAASBYratio,
                 type = 1)

# save to a *.pdf file (default)
plot.WAASBYratio(WAASBYratio,
                 type = 1,
                 export = T)


## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center", message = F, warning = F-----------------------------------------------------------------------------------------------------

plot.WAASBYratio(WAASBYratio,
                 type = 2)

#save to a *tiff file
plot.WAASBYratio(WAASBYratio,
                type = 2,
               export = T,
              file.type = "tiff",
             resolution = 600)


