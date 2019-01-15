## ---- eval = FALSE, collapse=TRUE, comment = "#"-------------------------
#  # download the package from Github
#  devtools::install_github("TiagoOlivoto/METAAB")
#  

## ----warning = FALSE, message = FALSE------------------------------------
# Importing data
require(METAAB)
require(ggplot2)
require(cowplot) # used in this material to arrange the graphics
dataset = read.csv("https://data.mendeley.com/datasets/2sjz32k3s3/2/files/1561de7f-b8fd-4093-b4d1-bfcef299dd22/WAASBdata.csv?dl=1")


## ----eval = TRUE, collapse=TRUE, comment = "#"---------------------------

AMMIweat = validation.AMMIF(dataset,
                            resp = GY,
                            gen = GEN,
                            env = ENV,
                            rep = BLOCK,
                            nboot = 50,
                            nrepval = 2)

# cross-validation for BLUP model
BLUPweat = validation.blup(dataset,
                            resp = GY,
                            gen = GEN,
                            env = ENV,
                            rep = BLOCK,
                            nboot = 50,
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


## ----eval = TRUE, fig.height = 5, fig.width = 10, fig.align = "center", message = F, warning = F----
# binding AMMI and BLUP RMSPDs
RMSPDweat = list(RMSPD = rbind(AMMIweat$RMSPD, 
                               BLUPweat$RMSPD))
class(RMSPDweat) = "validation.AMMIF"

# Plotting the RMSPD values
p1 = plot(RMSPDweat)
p2 = plot(RMSPDweat, width.boxplot = 0.5, col.boxplot = "transparent")
plot_grid(p1, p2, labels = c("p1", "p2"))

## ----echo = FALSE--------------------------------------------------------
AMMI_model = WAAS.AMMI(dataset,
                       resp = GY,
                       gen = GEN,
                       env = ENV,
                       rep = BLOCK)
predictoat = predict(AMMI_model, naxis = 5)

library(kableExtra)
options(digits = 4)
predictoat = predictoat[1:10,]
kable(predictoat, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", full_width = F, position = "float_left", font_size = 12)


## ----echo = TRUE---------------------------------------------------------
# Assuming equal weights for productivity and stability (default)
WAAS1 = WAAS.AMMI(dataset,
                  resp = GY,
                  gen = GEN,
                  env = ENV,
                  rep = BLOCK)

## ------------------------------------------------------------------------
# printing the WAAS object
options(digits = 3)
data = WAAS1$individual$individual
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", full_width = F, position = "left", font_size = 12)


## ------------------------------------------------------------------------
# printing the WAAS object
data = WAAS1$anova
kable(data,  align  = "l", booktabs = T, format = "html", linesep = "") %>%
kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12) %>%
row_spec(9, bold = T) %>%
add_indent(c(5:13))


## ------------------------------------------------------------------------
options(digits = 4)
data = WAAS1$model[, c(1:3,13:17, 21:22)]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)



## ----echo = TRUE---------------------------------------------------------
WAAS2 = WAAS.AMMI(dataset,
                  resp = GY,
                  gen = GEN,
                  env = ENV,
                  rep = BLOCK,
                  naxis = 7)


## ------------------------------------------------------------------------
stab_indexes = AMMI_indexes(AMMI_model)
kable(stab_indexes, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ----echo = TRUE---------------------------------------------------------
# Assuming equal weights for productivity and stability (default)
WAASB = WAASB(dataset,
              resp = GY,
              gen = GEN,
              env = ENV,
              rep = BLOCK)


## ----echo = TRUE, fig.width=10, fig.height=10, fig.align="center"--------
autoplot(WAASB)
autoplot(WAASB, type = "re")


## ----warning=F, message=F------------------------------------------------
options(digits = 5)
data = WAASB$LRT
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)

## ------------------------------------------------------------------------
options(digits = 7)
data = WAASB$ESTIMATES
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ------------------------------------------------------------------------
options(digits = 4)
data = WAASB$Details
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ------------------------------------------------------------------------
options(digits = 4)
data = WAASB$model[, c(1:3,13:17, 21:22)]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ------------------------------------------------------------------------
options(digits = 4)
data = WAASB$BLUPgen[1:10,]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ----echo = TRUE, fig.height = 5, fig.width = 10, fig.align = "center", message = F, warning = F----
# No file exported
p1 = plot.blup(WAASB)
p2 = plot.blup(WAASB, 
               col.shape  =  c("gray20", "gray80")) + coord_flip()
plot_grid(p1, p2,
          labels = c("p1", "p2"))


## ------------------------------------------------------------------------
options(digits = 4)
data = WAASB$BLUPgge[1:10,]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ------------------------------------------------------------------------
options(digits = 4)
data = WAASB$PCA
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)


## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center", message = F, warning = F----

plot.eigen(WAASB, size.lab = 14, size.tex = 14)


## ------------------------------------------------------------------------
options(digits = 4)
data = WAASB$MeansGxE[1:10,]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", position = "left", full_width = F, font_size = 12)

## ----echo = TRUE, fig.height = 5, fig.width = 10, fig.align = "center", message = F, warning = F----
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

## ----echo = TRUE, fig.height = 5, fig.width = 10, fig.align = "center", message = F, warning = F----

p3 = plot.scores(WAASB, type = 2)
p4 = plot.scores(WAASB, type = 2,
                 col.segm.env = "transparent") +
                 theme_gray() +
                 theme(legend.position = c(0.1, 0.9),
                       legend.background = element_rect(fill = NA))

plot_grid(p3, p4, labels = c("p3","p4"))

