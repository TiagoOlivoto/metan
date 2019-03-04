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

