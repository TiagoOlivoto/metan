## ----warning = FALSE, message = FALSE------------------------------------
library(METAAB)
library(ggplot2)
library(cowplot)
library(kableExtra)
dataset = data_ge
str(dataset)

## ----echo = TRUE---------------------------------------------------------
WAASB_model = WAASB(dataset,
                    resp = GY,
                    gen = GEN,
                    env = ENV,
                    rep = REP)


## ----echo = TRUE, fig.width=7, fig.height=7, fig.align="center"----------
autoplot(WAASB_model$GY)
autoplot(WAASB_model$GY, type = "re")

## ----warning=F, message=F------------------------------------------------
options(digits = 5)
data = WAASB_model$GY$LRT
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)

## ------------------------------------------------------------------------
data = WAASB_model$GY$ESTIMATES
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ------------------------------------------------------------------------
data = WAASB_model$GY$Details
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ------------------------------------------------------------------------
data = WAASB_model$GY$model[, c(1:3,13:17, 21:22)]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ------------------------------------------------------------------------
data = WAASB_model$GY$blupGEN[1:5,]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ----echo = TRUE, fig.height = 5, fig.width = 10, fig.align = "center", message = F, warning = F----
# No file exported
p1 = plot.blup(WAASB_model$GY)
p2 = plot.blup(WAASB_model$GY, 
               col.shape  =  c("gray20", "gray80")) + coord_flip()
plot_grid(p1, p2,
          labels = c("p1", "p2"))


## ------------------------------------------------------------------------
data = WAASB_model$GY$BLUPgge[1:5,]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ------------------------------------------------------------------------
data = WAASB_model$GY$PCA
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center", message = F, warning = F----

plot.eigen(WAASB_model$GY, size.lab = 14, size.tex.lab = 14)


## ------------------------------------------------------------------------
data = WAASB_model$GY$MeansGxE[1:10,]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)

## ----echo = TRUE, fig.height = 5, fig.width = 10, fig.align = "center", message = F, warning = F----
p1 = plot.scores(WAASB_model$GY, type = 1)
p2 = plot.scores(WAASB_model$GY,
                 type = 1,
                 polygon = TRUE,
                 col.gen = "black",
                 col.env = "gray70",
                 col.segm.env = "gray70",
                 axis.expand = 1.5)
plot_grid(p1, p2, labels = c("p1","p2"))

## ----echo = TRUE, fig.height = 5, fig.width = 10, fig.align = "center", message = F, warning = F----

p3 = plot.scores(WAASB_model$GY, type = 2)
p4 = plot.scores(WAASB_model$GY, type = 2,
                 col.segm.env = "transparent") +
                 theme_gray() +
                 theme(legend.position = c(0.1, 0.9),
                       legend.background = element_rect(fill = NA))

plot_grid(p3, p4, labels = c("p3","p4"))

## ----echo = TRUE, fig.height = 5, fig.width = 10, fig.align = "center", message=F, warning=F----

p5 = plot.scores(WAASB_model$GY, type = 3)
p6 = plot.scores(WAASB_model$GY, type = 3,
                 x.lab = "My customized x label",
                 size.shape = 3,
                 size.tex.pa = 2,
                 x.lim = c(1.2, 4.7),
                 x.breaks = seq(1.5, 4.5, by = 0.5)) + 
                 theme(legend.position = c(0.1, 0.9))
plot_grid(p5, p6, labels = c("p5","p6"))

## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center", message=F, warning=F----

plot.scores(WAASB_model$GY,
            type = 4, size.tex.pa = 1.5)



## ----echo = TRUE, fig.height = 4, fig.width = 10, fig.align = "center", message = F, warning = F----
p1 = plot.WAASBY(WAASB_model$GY)
p2 = plot.WAASBY(WAASB_model$GY, col.shape = c("gray20", "gray80"))
plot_grid(p1, p2, labels = c("p1", "p2"))

## ----echo = TRUE---------------------------------------------------------
WAASBYratio = WAASBYratio(dataset,
                          resp = GY,
                          gen = GEN,
                          env = ENV,
                          rep = REP,
                          increment = 50)


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
p1 = plot.WAASBY(WAASBYratio)
p2 = plot.WAASBY(WAASBYratio, col.shape = c("gray20", "gray80"))
plot_grid(p1, p2, labels = c("p1", "p2"))

## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center", message = F, warning = F----

plot(WAASBYratio, type = 1)



## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center", message = F, warning = F----

plot(WAASBYratio, type = 2)



## ------------------------------------------------------------------------
Resende_model = WAASB(dataset,
                      resp = c(GY, HM),
                      gen = GEN,
                      env = ENV,
                      rep = REP)
res_inde = Resende_indexes(Resende_model)
kable(res_inde, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


