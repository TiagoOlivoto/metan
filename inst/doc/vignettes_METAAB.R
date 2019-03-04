## ----warning = FALSE, message = FALSE------------------------------------
# Importing data
library(METAAB)
dataset = data_ge

## ------------------------------------------------------------------------
multivariate = WAASB(dataset,
                     resp = c(GY, HM),
                     gen = GEN,
                     env = ENV,
                     rep = REP)

FAI = FAI.BLUP(multivariate,
               SI = 10,
               DI = c("max", "max"), 
               UI = c("min", "min"))

plot(FAI)

