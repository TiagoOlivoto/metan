## ---- message=FALSE, warning=FALSE---------------------------------------
library(metan)
library(cowplot) # used to arrange the graphics
library(kableExtra) # Used to make the tables
library(magrittr) # for the forward-pipe operator %>%

# Function to make HTML tables
print_table = function(table){
  kable(table, "html", digits = 3) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), font_size = 12)
}
str(data_ge)

## ---- fig.height=3.5, fig.width=5----------------------------------------
CVAL = cv_ammif(data_ge,
                resp = GY,
                gen = GEN,
                env = ENV,
                rep = REP,
                nboot = 50,
                nrepval = 2)
plot(CVAL)

## ------------------------------------------------------------------------
model <- data_ge %>% waas(ENV, GEN, REP, GY)

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
model2 <- data_ge %>% waasb(ENV, GEN, REP, GY)

## ----fig.height=12, fig.width=4------------------------------------------
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
                col.shape  =  c("gray20", "gray80")) + coord_flip()
plot_grid(p1, p2,
          align = "v",
          labels = c("p1", "p2"),
          ncol = 1)

## ------------------------------------------------------------------------
print_table(head(model2$GY$BLUPgge))

## ------------------------------------------------------------------------
index = model2 %>%
  Resende_indexes()
  print_table(index$GY)

## ------------------------------------------------------------------------
stat_ge <- data_ge %>% ge_stats(ENV, GEN, REP, GY)

## ---- eval=FALSE---------------------------------------------------------
#  summary(stat_ge, export = TRUE)

