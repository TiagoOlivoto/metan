## ----global_options, include = FALSE-------------------------------------
knitr::opts_chunk$set(comment = "", cache = TRUE)


## ---- message=FALSE, warning=FALSE---------------------------------------
library(metan)
library(cowplot) # used to arrange the graphics
library(kableExtra) # Used to make the tables

# Function to make HTML tables
print_table = function(table){
  kable(table, "html", digits = 3) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed",        "responsive"), font_size = 12)
}
str(data_ge)

## ---- fig.height=3.5, fig.width=5----------------------------------------
CVAL = cv_ammif(data_ge,
                env = ENV,                
                gen = GEN,
                rep = REP,
                resp = GY,
                nboot = 500,
                nrepval = 2)
plot(CVAL)

## ------------------------------------------------------------------------
model <- waas(data_ge, ENV, GEN, REP, GY)

## ---- fig.height=8, fig.width=5,  message=FALSE, warning=FALSE-----------
p1 <- plot_scores(model$GY, axis.expand = 1.2)
p2 <- plot_scores(model$GY,
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
print_table(head(predicted$GY))

## ------------------------------------------------------------------------
model2 <- waasb(data_ge, ENV, GEN, REP, GY)

## ----fig.height=12, fig.width=4------------------------------------------
library(ggplot2)
autoplot(model2$GY, which = c(1, 2, 7), mfrow = c(3,1))

## ----fig.height=8, fig.width=4-------------------------------------------
autoplot(model2$GY, type = "re", which = c(1, 2), mfrow = c(2,1))

## ------------------------------------------------------------------------
print_table(model2$GY$LRT)

print_table(model2$GY$random)

print_table(model2$GY$ESTIMATES)

print_table(model2$GY$blupGEN)

## ---- fig.height=8, fig.width=4------------------------------------------
p1 <- plot_blup(model2$GY)
p2 <- plot_blup(model2$GY,
                prob = 0.1,
                col.shape  =  c("gray20", "gray80")) +
      coord_flip()
plot_grid(p1, p2,
          align = "v",
          labels = c("p1", "p2"),
          ncol = 1)

## ------------------------------------------------------------------------
print_table(head(model2$GY$BLUPgge))

## ------------------------------------------------------------------------
get_model_data(model2, what = "WAASB") %>% 
  print_table()

## ------------------------------------------------------------------------
index = model2 %>%
  Resende_indexes()
  print_table(index$GY)

## ----echo = TRUE---------------------------------------------------------
# Using a data.frame
gge_model = gge(data_ge, ENV, GEN, GY)

# Using a two-way table
ge_table = make_mat(data_ge, GEN, ENV, GY)
gge_model = gge(ge_table, table = TRUE)


## ----echo = TRUE, fig.width = 4, fig.height=8, message=F, warning=F------
p1 = plot(gge_model)
p2 = plot(gge_model, type = 2)
plot_grid(p1, p2,
          align = "v",
          labels = c("p1", "p2"),
          ncol = 1)

## ------------------------------------------------------------------------
stat_ge <- data_ge %>% ge_stats(ENV, GEN, REP, GY)

## ---- eval=FALSE---------------------------------------------------------
#  summary(stat_ge, export = TRUE)
#  

## ---- eval=FALSE, warning=FALSE------------------------------------------
#  data_ge2 %>%
#    waasb(ENV, GEN, REP,
#          resp = c(KW, NKE, PH, EH, TKW),
#          verbose = FALSE) %>%
#    mtsi(verbose = FALSE) %>%
#    plot()
#  

