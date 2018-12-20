autoplot.WAASB = function(x,
                          type = "res",
                          conf = 0.95,
                          labels = FALSE,
                          theme = theme_waasb(),
                          alpha = 0.2,
                          fill.hist = "gray",
                          col.hist = "black",
                          col.point = "black",
                          col.line = "red",
                          col.lab.out = "red",
                          size.lab.out = 2.5,
                          bins = 30,
                          which = c(1:4),
                          mfrow = c(2 , 2),
                          ...){
if (type == "re" & max(which) >= 5){
  stop("When type =\"re\", 'which' must be a value between 1 and 4" )
}
if (type == "res"){
  df = data.frame(x$residuals)
  df$id = rownames(df)
  df = data.frame(df[order(df$.scresid),])
  P <- ppoints(nrow(df))
  df$z = qnorm(P)
  n <- nrow(df)
  Q.x <- quantile(df$.scresid, c(0.25, 0.75))
  Q.z <- qnorm(c(0.25, 0.75))
  b <- diff(Q.x)/diff(Q.z)
  coef <- c(Q.x[1] - b * Q.z[1], b)
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (coef[2]/dnorm(df$z)) * sqrt(P * (1 - P)/n)
  fit.value <- coef[1] + coef[2] * df$z
  df$upper <- fit.value + zz * SE
  df$lower <- fit.value - zz * SE
  df$label <- ifelse(df$.scresid > df$.scresid | df$.scresid < df$lower, rownames(df),"")
  df$factors = paste(df$ENV, df$GEN)

  # Residuals vs .fitted

p1 = ggplot(df, aes(.fitted, .resid)) +
    geom_point(col = col.point)  +
    geom_smooth(se = F, method = "loess", col = col.line) +
    geom_hline(yintercept = 0, linetype = 2, col = "gray")+
    labs(x = "Fitted values",  y = "Residual")+
    theme
if (labels != FALSE){
  p1 = p1 + ggrepel::geom_text_repel(aes(.fitted, .resid,
                                      label = (label)),
                                  color = col.lab.out,
                                  size = size.lab.out)
} else{p1 = p1}

  # normal qq
p2 = ggplot(df, aes(z, .scresid)) +
    geom_point(col = col.point) +
    geom_abline(intercept = coef[1], slope = coef[2],
                size = 1,
                col = col.line) +
    geom_ribbon(aes_(ymin = ~lower, ymax = ~upper), alpha = 0.2)+
    labs(x = "Theoretical quantiles",
         y = "Sample quantiles") +
    ggtitle("Normal Q-Q")+
    theme
if (labels != FALSE){
p2 = p2 + ggrepel::geom_text_repel(aes(z, .scresid,
                               label = (label)),
                           color = col.lab.out,
                           size = size.lab.out)
} else {p2 = p2}


# scale-location
p3 = ggplot(df, aes(.fitted, sqrt(abs(.resid))))+
  geom_point(col = col.point)  +
  geom_smooth(se = F, method = "loess", col = col.line) +
  labs(x = "Fitted Values",
       y = expression(sqrt("|Standardized residuals|"))) +
  ggtitle("Scale-location") +
  theme
if (labels != FALSE){
p3 = p3 + ggrepel::geom_text_repel(aes(.fitted, sqrt(abs(.resid)),
                               label = (label)),
                           color = col.lab.out,
                           size = size.lab.out)
} else {p3 = p3}

# Residuals vs Factor-levels
p4 = ggplot(df, aes(factors, .scresid))+
  geom_point(col = col.point)  +
  geom_hline(yintercept = 0, linetype = 2, col = "gray")+
  labs(x = "Fitted values",
       y = "Standardized residuals") +
  ggtitle("Residuals vs factor-levels") +
  theme
if (labels != FALSE){
p4 = p4 + ggrepel::geom_text_repel(aes(factors, .scresid,
                               label = (label)),
                           color = col.lab.out,
                           size = size.lab.out)
} else {p4 = p4}

# Histogram of residuals
p5 = ggplot(df, aes(x = .resid)) +
  geom_histogram(bins = bins, colour = col.hist, fill = fill.hist,
                 aes(y = ..density..)) +
  stat_function(fun = dnorm,
                color = col.line,
                size = 1,
                args = list(mean = mean(df$.resid),
                            sd = sd(df$.resid))) +
  labs(x = "Raw residuals", y = "Density")+
  ggtitle("Histogram of residuals") +
  theme

# Residuals vs order
p6 = ggplot(df, aes(as.numeric(id), .scresid, group = 1))+
  geom_point(col = col.point)  +
  geom_line(col = col.line) +
  geom_hline(yintercept = 0, linetype = 2, col = col.line) +
  labs(x = "Observation order",
       y = "Standardized residuals") +
  ggtitle("Residuals vs observation order") +
  theme



p7 = ggplot(df, aes(.fitted, Y)) +
  geom_point(col = col.point) +
  facet_wrap(~GEN) +
  geom_abline(intercept = 0,slope = 1, col = col.line) +
  theme_waasb() +
  theme(panel.spacing = unit(0, "cm")) +
  labs(x = "Fitted values", y = "Observed values")


plots <- list(p1, p2, p3, p4, p5, p6, p7)

}
if (type == "re"){
  df = data.frame(x$BLUPgge)[,1:4]
  df$id = rownames(df)
  df = data.frame(df[order(df$BLUPge),])
  P <- ppoints(nrow(df))
  df$z = qnorm(P)
  n <- nrow(df)
  Q.x <- quantile(df$BLUPge, c(0.25, 0.75))
  Q.z <- qnorm(c(0.25, 0.75))
  b <- diff(Q.x)/diff(Q.z)
  coef <- c(Q.x[1] - b * Q.z[1], b)
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (coef[2]/dnorm(df$z)) * sqrt(P * (1 - P)/n)
  fit.value <- coef[1] + coef[2] * df$z
  df$upper <- fit.value + zz * SE
  df$lower <- fit.value - zz * SE
  df$label <- ifelse(df$BLUPg > df$BLUPg | df$BLUPg < df$lower, rownames(df),"")
  df$factors = paste(df$ENV, df$GEN)

  dfgen = data.frame(x$BLUPgen)[,2:3]
  dfgen$id = rownames(dfgen)
  dfgen = data.frame(dfgen[order(dfgen$BLUPg),])
  P2 <- ppoints(nrow(dfgen))
  dfgen$z = qnorm(P2)
  n2 <- nrow(dfgen)
  Qx <- quantile(dfgen$BLUPg, c(0.25, 0.75))
  Qz <- qnorm(c(0.25, 0.75))
  b2 <- diff(Qx)/diff(Qz)
  coef2 <- c(Qx[1] - b2 * Qz[1], b2)
  zz2 <- qnorm(1 - (1 - conf)/2)
  SE2 <- (coef2[2]/dnorm(dfgen$z)) * sqrt(P2 * (1 - P2)/n2)
  fit.value2 <- coef2[1] + coef2[2] * dfgen$z
  dfgen$upper <- fit.value2 + zz2 * SE2
  dfgen$lower <- fit.value2 - zz2 * SE2
  dfgen$label <- ifelse(dfgen$BLUPg > dfgen$BLUPg | dfgen$BLUPg < dfgen$lower, rownames(dfgen),"")

  # normal qq GEI effects
p1 = ggplot(df, aes(z, BLUPge)) +
    geom_point(col = col.point) +
    geom_abline(intercept = coef[1], slope = coef[2],
                size = 1,
                col = col.line) +
    geom_ribbon(aes_(ymin = ~lower, ymax = ~upper), alpha = 0.2) +
    labs(x = "Theoretical quantiles",
         y = "Sample quantiles") +
    ggtitle("Q-Q | GEI effects") +
    theme
  if (labels != FALSE){
    p1 = p1 + ggrepel::geom_text_repel(aes(z, BLUPge,
                                           label = (label)),
                                       color = col.lab.out,
                                       size = size.lab.out)
  } else {p1 = p1}


# normal qq Genotype effects
p2 = ggplot(dfgen, aes(z, BLUPg)) +
  geom_point(col = col.point) +
  geom_abline(intercept = coef2[1], slope = coef2[2],
              size = 1,
              col = col.line) +
  geom_ribbon(aes_(ymin = ~lower, ymax = ~upper), alpha = 0.2) +
  labs(x = "Theoretical quantiles",
       y = "Sample quantiles") +
  ggtitle("Q-Q | genotype effects") +
  theme
if (labels != FALSE){
  p2 = p2 + ggrepel::geom_text_repel(aes(z, BLUPg,
                                         label = (label)),
                                     color = col.lab.out,
                                     size = size.lab.out)
} else {p2 = p2}

# random effects vs Factor-levels
p3 = ggplot(df, aes(BLUPge, factors))+
  geom_point(col = col.point)  +
  geom_vline(xintercept = 0, linetype = 2, col = "gray")+
  labs(x = "Random effects",
       y = "Factor-levels") +
  ggtitle("Random effects vs factor-levels") +
  theme +
theme(axis.text.y = element_blank())


# random effects vs genotypes
p4 = ggplot(dfgen, aes(BLUPg, GEN))+
  geom_point(col = col.point)  +
  geom_vline(xintercept = 0, linetype = 2, col = "gray")+
  labs(x = "Random effects ",
       y = "Genotypes") +
  ggtitle("Random effects vs Genotypes") +
  theme

plots <- list(p1, p2, p3, p4)

}

  # making the plots
grid::grid.newpage()

  if (prod(mfrow)>1) {
    mypos <- expand.grid(1:mfrow[1], 1:mfrow[2])
    mypos <- mypos[with(mypos, order(Var1)), ]
    grid::pushViewport(viewport(layout = grid.layout(mfrow[1], mfrow[2])))
    formatter <- function(.){}
  } else {
    mypos <- data.frame(matrix(1, length(which), 2))
    grid::pushViewport(viewport(layout = grid.layout(1, 1)))
    formatter <- function(.) {
      .dontcare <- readline("Hit <Return> to see next plot: ")
      grid::grid.newpage()
    }
  }

  j <- 1
  for (i in which){
    formatter()
    print(plots[[i]], vp = viewport(layout.pos.row = mypos[j,][1], layout.pos.col=mypos[j,][2]))
    j <- j + 1
  }
}



