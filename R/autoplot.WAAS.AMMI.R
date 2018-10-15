autoplot.WAAS.AMMI = function(x,
                              conf = 0.95,
                              theme = theme_waasb(),
                              alpha = 0.2,
                              fill.hist = "gray",
                              col.hist = "black",
                              col.point = "black",
                              col.line = "red",
                              col.lab.out = "red",
                              bins = 30,
                              which = c(1:4),
                              mfrow = c(2 , 2),
                              ...){

  df = x$residuals
  df$id = rownames(df)
  df = data.frame(df[order(df$stdres),])
  P <- ppoints(nrow(df))
  df$z = qnorm(P)
  n <- nrow(df)
  Q.x <- quantile(df$stdres, c(0.25, 0.75))
  Q.z <- qnorm(c(0.25, 0.75))
  b <- diff(Q.x)/diff(Q.z)
  coef <- c(Q.x[1] - b * Q.z[1], b)
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (coef[2]/dnorm(df$z)) * sqrt(P * (1 - P)/n)
  fit.value <- coef[1] + coef[2] * df$z
  df$upper <- fit.value + zz * SE
  df$lower <- fit.value - zz * SE
  df$label <- ifelse(df$stdres > df$upper | df$stdres < df$lower, rownames(df),"")


  # residuals vs fitted

p1 = ggplot(df, aes(fitted, resid)) +
    geom_point(col = col.point)  +
    geom_smooth(se = F, method = "loess", col = col.line) +
    geom_hline(yintercept = 0, linetype = 2, col = "gray")+
    labs(x = "Fitted Values",
         y = "Residual") +
    geom_text(aes(label = label), hjust = "inward", col = col.lab.out)+
    ggtitle("Residuals vs Fitted") +
    theme



  # normal qq
p2 = ggplot(df, aes(z, stdres)) +
    geom_point(col = col.point) +
    geom_abline(intercept = coef[1], slope = coef[2],
                col = col.line,
                size = 1) +
    geom_ribbon(aes_(ymin = ~lower, ymax = ~upper), alpha = 0.2) +
    labs(x = "Theoretical Quantiles",
         y = "Standardized Residuals") +
    ggtitle("Normal Q-Q") +
    geom_text(aes(label = label), hjust = "inward", col = col.lab.out)+
  theme



# scale-location
p3 = ggplot(df, aes(fitted, sqrt(abs(resid))))+
  geom_point(col = col.point)  +
  geom_smooth(se = F, method = "loess", col = col.line) +
  labs(x = "Fitted Values",
       y = expression(sqrt("|Standardized Residuals|"))) +
  geom_text(aes(label = label), hjust = "inward", col = col.lab.out)+
  ggtitle("Scale-location") +
  theme

# Residuals vs Factor-levels
p4 = ggplot(df, aes(factors, stdres))+
  geom_point(col = col.point)  +
  geom_hline(yintercept = 0, linetype = 2, col = "gray")+
  labs(x = "Fitted Values",
       y = "Standardized Residuals") +
  geom_text(aes(label = label), hjust = "inward", col = col.lab.out)+
  ggtitle("Residuals vs Factor-levels") +
  theme


# Histogram of residuals
p5 = ggplot(df, aes(x = resid)) +
  geom_histogram(bins = bins, colour = col.hist, fill = fill.hist,
                 aes(y = ..density..)) +
  stat_function(fun = dnorm,
                color = "red",
                size = 1,
                args = list(mean = mean(df$resid),
                            sd = sd(df$resid))) +
  labs(x = "Raw residuals", y = "Density")+
  theme

# Residuals vs order
p6 = ggplot(df, aes(as.numeric(id), stdres, group = 1))+
  geom_point(col = col.point)  +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 2, col = col.line) +
  labs(x = "Observation order",
       y = "Standardized Residuals") +
  ggtitle("Residuals vs Observation order") +
  theme



plots <- list(p1, p2, p3, p4, p5, p6)

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



