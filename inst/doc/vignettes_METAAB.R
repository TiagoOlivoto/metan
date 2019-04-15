## ---- message=FALSE, warning=FALSE---------------------------------------
library(METAAB)
library(kableExtra)
library(cowplot)
library(magrittr)
str(data_ge)

## ---- fig.height=3.5, fig.width=5----------------------------------------
CVAL = validation.AMMIF(data_ge,
                        resp = GY,
                        gen = GEN,
                        env = ENV,
                        rep = REP,
                        nboot = 50,
                        nrepval = 2)
plot(CVAL)

## ------------------------------------------------------------------------
model <- data_ge %>% WAAS.AMMI(ENV, GEN, REP, GY)

## ---- fig.height=8, fig.width=5,  message=FALSE, warning=FALSE-----------
p1 <- plot.scores(model$GY)
p2 <- plot.scores(model$GY,
                 type = 1,
                 polygon = TRUE,
                 col.gen = "black",
                 col.env = "gray70",
                 col.segm.env = "gray70",
                 axis.expand = 1.5)
plot_grid(p1, p2,
          align = "v",
          labels = c("p1","p2"),
          ncol = 1)

## ------------------------------------------------------------------------
predicted <- predict(model, naxis = 4)
predicted <- head(predicted$GY)
kable(predicted, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", full_width = F)

## ------------------------------------------------------------------------
model2 <- data_ge %>% WAASB(ENV, GEN, REP, GY)

## ------------------------------------------------------------------------
LRT = model2$GY$LRT
kable(LRT, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", full_width = F)

VC = model2$GY$random
kable(VC, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", full_width = F)

ESTIMATES = model2$GY$ESTIMATES
kable(ESTIMATES, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", full_width = F)

BP = model2$GY$blupGEN
kable(BP, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", full_width = F)

## ---- fig.height=8, fig.width=5------------------------------------------
p1 <- plot.blup(model2$GY)
p2 <- plot.blup(model2$GY,
                prob = 0.1,
                col.shape  =  c("gray20", "gray80")) + coord_flip()
plot_grid(p1, p2,
          align = "v",
          labels = c("p1", "p2"),
          ncol = 1)

## ------------------------------------------------------------------------
data <- head(model2$GY$BLUPgge)
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)

