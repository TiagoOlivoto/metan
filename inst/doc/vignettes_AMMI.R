## ----warning = FALSE, message = FALSE------------------------------------
library(METAAB)
library(ggplot2)
library(cowplot) # used in this material to arrange the graphics
library(kableExtra) # Used to 
dataset = data_ge
str(dataset)

## ----echo = TRUE---------------------------------------------------------
AMMI_model = WAAS.AMMI(data = dataset,
                       resp = GY,
                       gen = GEN,
                       env = ENV,
                       rep = REP,
                       verbose = FALSE)


## ------------------------------------------------------------------------
# printing the WAAS object
options(digits = 3)
data = AMMI_model$GY$individual$individual
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                full_width = F, position = "left", font_size = 12)


## ----echo = TRUE---------------------------------------------------------
data = AMMI_model$GY$anova
kable(data,  align  = "l", booktabs = T, format = "html", linesep = "") %>%
kable_styling(bootstrap_options = "striped", "condensed",
              position = "left", full_width = F, font_size = 12) %>%
row_spec(5:8, bold = T) %>%
add_indent(c(5:13))

## ----echo = TRUE---------------------------------------------------------
predicted = predict(AMMI_model, naxis = 4)
predicted = predicted$GY[1:5,]
kable(predicted, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed", full_width = T)

## ------------------------------------------------------------------------
data = AMMI_model$GY$model[, c(1:3,13:17, 21:22)]
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)



## ----echo = TRUE---------------------------------------------------------
AMMI_model_2 = WAAS.AMMI(dataset,
                         resp = GY,
                         gen = GEN,
                         env = ENV,
                         rep = REP,
                         naxis = 7, # Use 7 IPCA for computing WAAS
                         verbose = FALSE)


## ----echo = TRUE, fig.height = 5, fig.width = 10, fig.align = "center", message = F, warning = F----
library(cowplot)
p1 = plot.scores(AMMI_model$GY, type = 1)
p2 = plot.scores(AMMI_model$GY,
                 type = 1,
                 polygon = TRUE,
                 col.gen = "black",
                 col.env = "gray70",
                 col.segm.env = "gray70",
                 axis.expand = 1.5)
plot_grid(p1, p2, labels = c("p1","p2"))

## ----echo = TRUE, fig.height = 5, fig.width = 10, fig.align = "center", message = F, warning = F----

p3 = plot.scores(AMMI_model$GY, type = 2)
p4 = plot.scores(AMMI_model$GY, type = 2,
                 col.segm.env = "transparent") +
                 theme_gray() +
                 theme(legend.position = c(0.1, 0.9),
                       legend.background = element_rect(fill = NA))

plot_grid(p3, p4, labels = c("p3","p4"))

## ----echo = TRUE, fig.height = 5, fig.width = 10, fig.align = "center", message=F, warning=F----

p5 = plot.scores(AMMI_model$GY, type = 3)
p6 = plot.scores(AMMI_model$GY, type = 3,
                 x.lab = "My customized x label",
                 size.shape = 3,
                 size.tex.pa = 2,
                 x.lim = c(1.2, 4.7),
                 x.breaks = seq(1.5, 4.5, by = 0.5)) + 
                 theme(legend.position = c(0.1, 0.9))
plot_grid(p5, p6, labels = c("p5","p6"))

## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center", message=F, warning=F----

plot.scores(AMMI_model$GY,
            type = 4, size.tex.pa = 1.5)



## ----echo = TRUE, fig.height = 4, fig.width = 10, fig.align = "center", message = F, warning = F----
library(ggplot2)
p1 = plot.WAASBY(AMMI_model$GY)
p2 = plot.WAASBY(AMMI_model$GY, col.shape = c("gray20", "gray80"))
plot_grid(p1, p2, labels = c("p1", "p2"))

## ----echo = TRUE---------------------------------------------------------
WAASratio = WAASratio.AMMI(dataset,
                           resp = GY,
                           gen = GEN,
                           env = ENV,
                           rep = REP,
                           increment = 10)


## ------------------------------------------------------------------------

options(digits = 4)
data = WAASratio$hetcomb
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ------------------------------------------------------------------------

options(digits = 4)
data = WAASratio$hetdata
kable(data, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center", message = F, warning = F----

plot(WAASratio, type = 1)



## ----echo = TRUE, fig.height = 5, fig.width = 5.5, fig.align = "center", message = F, warning = F----

plot(WAASratio, type = 2)



## ------------------------------------------------------------------------
stab_indexes = AMMI_indexes(AMMI_model, order.y = NULL)
kable(stab_indexes, "html") %>%
  kable_styling(bootstrap_options = "striped", "condensed",
                position = "left", full_width = F, font_size = 12)


