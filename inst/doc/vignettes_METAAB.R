## ---- eval=FALSE---------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("TiagoOlivoto/METAAB")

## ------------------------------------------------------------------------
library(METAAB)
str(data_ge)

## ------------------------------------------------------------------------
model <- WAAS.AMMI(data_ge,
                  resp = GY,
                  gen = GEN,
                  env = ENV,
                  rep = REP)

## ---- fig.height=8, fig.width=5,  message=FALSE, warning=FALSE-----------
library(cowplot)
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
library(kableExtra)
options(digits = 4)
predicted <- predict(model, naxis = 4)
predicted <- predicted$GY[1:5,]
kable(predicted, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", full_width = F)

## ------------------------------------------------------------------------
model2 <- WAASB(data_ge,
                resp = GY,
                gen = GEN,
                env = ENV,
                rep = REP)

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
data <- model2$GY$BLUPgge[1:5,]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)

