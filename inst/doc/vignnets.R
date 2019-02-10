## ---- eval = FALSE, collapse=TRUE, comment = "#"-------------------------
#  # download the package from Github
#  devtools::install_github("TiagoOlivoto/METAAB")
#  

## ----warning = FALSE, message = FALSE------------------------------------
# Importing data
require(METAAB)
require(ggplot2)
require(cowplot) # used in this material to arrange the graphics
dataset = data_ge

## ----eval = TRUE, collapse=TRUE, comment = "#"---------------------------

AMMIweat = validation.AMMIF(dataset,
                            resp = GY,
                            gen = GEN,
                            env = ENV,
                            rep = REP,
                            nboot = 5,
                            nrepval = 2)

# cross-validation for BLUP model
BLUPweat = validation.blup(dataset,
                            resp = GY,
                            gen = GEN,
                            env = ENV,
                            rep = REP,
                            nboot = 5,
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
  kable_styling(bootstrap_options = "striped", "condensed",
                full_width = F, position = "left", font_size = 12)


## ----eval = TRUE, fig.height = 5, fig.width = 10, fig.align = "center", message = F, warning = F----
# binding AMMI and BLUP RMSPDs
RMSPDweat = list(RMSPD = rbind(AMMIweat$RMSPD, 
                               BLUPweat$RMSPD))
class(RMSPDweat) = "validation.AMMIF"

# Plotting the RMSPD values
p1 = plot(RMSPDweat)
p2 = plot(RMSPDweat, width.boxplot = 0.5, col.boxplot = "transparent")
plot_grid(p1, p2, labels = c("p1", "p2"))

## ----echo = TRUE---------------------------------------------------------
AMMI_model = WAAS.AMMI(dataset,
                       resp = GY,
                       gen = GEN,
                       env = ENV,
                       rep = REP,
                       verbose = FALSE)
predictoat = predict(AMMI_model, naxis = 5)

library(kableExtra)
options(digits = 4)
predictoat = predictoat[[1]][1:10,]
kable(predictoat, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                full_width = F, position = "float_left", font_size = 12)


## ----echo = TRUE---------------------------------------------------------
# Assuming equal weights for productivity and stability (default)
WAAS1 = WAAS.AMMI(dataset,
                  resp = GY,
                  gen = GEN,
                  env = ENV,
                  rep = REP)

## ------------------------------------------------------------------------
# printing the WAAS object
options(digits = 3)
data = WAAS1$GY$individual$individual
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                full_width = F, position = "left", font_size = 12)


## ------------------------------------------------------------------------
# printing the WAAS object
data = WAAS1$GY$anova
kable(data,  align  = "l", booktabs = T, format = "html", linesep = "") %>%
kable_styling(bootstrap_options = "striped", "condensed",
              position = "left", full_width = F, font_size = 12) %>%
row_spec(9, bold = T) %>%
add_indent(c(5:13))


## ------------------------------------------------------------------------
options(digits = 4)
data = WAAS1$GY$model[, c(1:3,13:17, 21:22)]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)



## ----echo = TRUE---------------------------------------------------------
WAAS2 = WAAS.AMMI(dataset,
                  resp = GY,
                  gen = GEN,
                  env = ENV,
                  rep = REP,
                  naxis = 7,
                  verbose = FALSE)


## ------------------------------------------------------------------------
stab_indexes = AMMI_indexes(AMMI_model)
kable(stab_indexes, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ----echo = TRUE---------------------------------------------------------
# Assuming equal weights for productivity and stability (default)
WAASB = WAASB(dataset,
              resp = GY,
              gen = GEN,
              env = ENV,
              rep = REP)


## ----echo = TRUE, fig.width=7, fig.height=7, fig.align="center"----------
library(ggplot2)
autoplot(WAASB$GY)


## ----warning=F, message=F------------------------------------------------
options(digits = 5)
data = WAASB$GY$LRT
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)

## ------------------------------------------------------------------------
options(digits = 7)
data = WAASB$GY$ESTIMATES
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ------------------------------------------------------------------------
options(digits = 4)
data = WAASB$GY$Details
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ------------------------------------------------------------------------
options(digits = 4)
data = WAASB$GY$model[, c(1:3,13:17, 21:22)]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ------------------------------------------------------------------------
options(digits = 4)
data = WAASB$GY$blupGEN[1:10,]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ----echo = TRUE, fig.height = 5, fig.width = 10, fig.align = "center", message = F, warning = F----
# No file exported
p1 = plot.blup(WAASB$GY)
p2 = plot.blup(WAASB$GY, 
               col.shape  =  c("gray20", "gray80")) + coord_flip()
plot_grid(p1, p2,
          labels = c("p1", "p2"))


## ------------------------------------------------------------------------
options(digits = 4)
data = WAASB$GY$BLUPgge[1:10,]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ------------------------------------------------------------------------
options(digits = 4)
data = WAASB$GY$PCA
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center", message = F, warning = F----

plot.eigen(WAASB$GY, size.lab = 14, size.tex.lab = 14)


## ------------------------------------------------------------------------
options(digits = 4)
data = WAASB$GY$MeansGxE[1:10,]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)

## ----echo = TRUE, fig.height = 5, fig.width = 10, fig.align = "center", message = F, warning = F----
library(cowplot)
p1 = plot.scores(WAASB$GY, type = 1)
p2 = plot.scores(WAASB$GY,
                 type = 1,
                 polygon = TRUE,
                 col.gen = "black",
                 col.env = "gray70",
                 col.segm.env = "gray70",
                 axis.expand = 1.5)
plot_grid(p1, p2, labels = c("p1","p2"))

## ----echo = TRUE, fig.height = 5, fig.width = 10, fig.align = "center", message = F, warning = F----

p3 = plot.scores(WAASB$GY, type = 2)
p4 = plot.scores(WAASB$GY, type = 2,
                 col.segm.env = "transparent") +
                 theme_gray() +
                 theme(legend.position = c(0.1, 0.9),
                       legend.background = element_rect(fill = NA))

plot_grid(p3, p4, labels = c("p3","p4"))

## ----echo = TRUE, fig.height = 5, fig.width = 10, fig.align = "center", message=F, warning=F----

p5 = plot.scores(WAASB$GY, type = 3)
p6 = plot.scores(WAASB$GY, type = 3,
                 x.lab = "My customized x label",
                 size.shape = 3,
                 size.tex.pa = 2,
                 x.lim = c(1.2, 4.7),
                 x.breaks = seq(1.5, 4.5, by = 0.5)) + 
                 theme(legend.position = c(0.1, 0.9))
plot_grid(p5, p6, labels = c("p5","p6"))

## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center", message=F, warning=F----

plot.scores(WAASB$GY,
            type = 4, size.tex.pa = 1.5)



## ----echo = TRUE---------------------------------------------------------
WAASBYratio = WAASBYratio(dataset,
                          resp = GY,
                          gen = GEN,
                          env = ENV,
                          rep = REP,
                          increment = 50,
                          saveWAASY = 50)


## ----echo = TRUE, eval=FALSE---------------------------------------------
#  WAASBYratio2 = WAASratio.AMMI(dataset,
#                                resp = GY,
#                                gen = GEN,
#                                env = ENV,
#                                rep = REP,
#                                increment = 50,
#                                saveWAASY = 50)
#  

## ------------------------------------------------------------------------
options(digits = 4)
data = WAASBYratio$WAASY
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ------------------------------------------------------------------------

options(digits = 4)
data = WAASBYratio$hetcomb
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ------------------------------------------------------------------------

options(digits = 4)
data = WAASBYratio$hetdata
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ----echo = TRUE, fig.height = 4, fig.width = 10, fig.align = "center", message = F, warning = F----
library(ggplot2)
p1 = plot.WAASBY(WAASBYratio)
p2 = plot.WAASBY(WAASBYratio, col.shape = c("gray20", "gray80"))
plot_grid(p1, p2, labels = c("p1", "p2"))

## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center", message = F, warning = F----

plot(WAASBYratio, type = 1)



## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center", message = F, warning = F----

plot(WAASBYratio, type = 2)



## ------------------------------------------------------------------------
res_inde = Resende_indexes(WAASB)
kable(res_inde, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ------------------------------------------------------------------------
multivariate = WAASB = WAASB(dataset,
               resp = c(GY, HM),
               gen = GEN,
               env = ENV,
               rep = REP)

FAI = FAI.BLUP(multivariate,
               SI = 10,
               DI = c("max", "max"), 
               UI = c("min", "min"))

plot(FAI)

## ------------------------------------------------------------------------
MTSI_MODEL = WAASB(data_ge,
                   resp = c(GY, HM),
                   gen = GEN,
                   env = ENV,
                   rep = REP)

MTSI_index = MTSI(MTSI_MODEL)

plot(MTSI_index)

## ------------------------------------------------------------------------
MTSI_MODEL = WAASB(data_ge,
                   resp = c(GY, HM),
                   gen = GEN,
                   env = ENV,
                   rep = REP,
                   mresp = c(100, 100), #Default
                   wresp = c(65, 65))

MTSI_index = MTSI(MTSI_MODEL, index = "WAASBY")

plot(MTSI_index)

