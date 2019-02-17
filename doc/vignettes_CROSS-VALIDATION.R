## ----warning = FALSE, message = FALSE------------------------------------
# Importing data
library(METAAB)
library(ggplot2)
library(cowplot) # used in this material to arrange the graphics
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
RMSEweat = RMSEweat[order(RMSEweat$mean, decreasing = F),]

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

