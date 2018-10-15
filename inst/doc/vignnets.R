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


## ----echo = TRUE---------------------------------------------------------
options(digits = 4)
RMSEweat = rbind(AMMIweat$RMSPDmean, BLUPweat$RMSPDmean)
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
RMSPDweat = list(RMSPD = rbind(AMMIweat$RMSPD, 
                           BLUPweat$RMSPD))
# Plotting the RMSPD values
plot.validation.AMMIF(RMSPDweat,
                      violin = FALSE,
                      col.boxplot = "gray75")

## ----echo = TRUE---------------------------------------------------------
model =  WAAS.AMMI(dataset,
                   resp = GY,
                   gen = GEN,
                   env = ENV,
                   rep = BLOCK)
predictoat = predict(model, naxis = 5)


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
                  resp = GY,
                  gen = GEN,
                  env = ENV,
                  rep = BLOCK,
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


## ----eval = FALSE, warning=F, message=F---------------------------------------------------------------------------------------------------------------------------------------------------------------
#  print(WAASB$LRT)

## ----echo = FALSE, warning=F, message=F---------------------------------------------------------------------------------------------------------------------------------------------------------------
library(kableExtra)
options(digits = 5)
dt_footnote = WAASB$LRT

kable(dt_footnote, "html", align = "c", escape = F) %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)

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

# dafault theme of ggplot2
library(ggplot2)
plot.scores(WAASB,
            type = 2,
            col.gen = "black",
            col.env = "blue",
            theme = theme_gray())


## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center", message=F, warning=F---------------------------------------------------------------------------------------------------------

plot.scores(WAASB,
            type = 3)

# Save to a *.tiff file with resolution of 600 dpi.
plot.scores(WAASB,
            type = 3,
            export = TRUE,
            file.type = "tiff",
            resolution = 600)


## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center", message=F, warning=F---------------------------------------------------------------------------------------------------------

plot.scores(WAASB,
            type = 4)



## ----echo = TRUE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

WAASBYratio = WAASBYratio(dataset,
                          resp = GY,
                          gen = GEN,
                          env = ENV,
                          rep = BLOCK,
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


## ----echo = TRUE, fig.height = 3, fig.width = 6, fig.align = "center"---------------------------------------------------------------------------------------------------------------------------------

plot.WAASBY(WAASBYratio,
            theme = theme_waasb() +
                    theme(aspect.ratio = NULL,
                          legend.position = c(0.85, 0.2)))


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


