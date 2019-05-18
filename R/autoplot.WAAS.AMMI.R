## autoplot WAAS.AMMI.R
#' Several types of residual plots
#'
#' Residual plots for a output model of class "WAAS.AMMI" and "WAASB". Six sets
#' of plots are produced: (1) Residuals vs fitted, (2) normal Q-Q plot for the
#' residuals, (3) scale-location plot (standardized residuals vs Fitted
#' Values), (4) standardized residuals vs Factor-levels, (5) Histogram of raw
#' residuals and (6) standardized residuals vs observation order. For a "WAASB"
#' object, normal Q-Q plot for random effects may also be obtained declaring
#' \code{type = "re"}
#'
#'
#' @param x An object of class \code{WAAS.AMMI} or \code{WAASB}.
#' @param conf Level of confidence interval to use in the Q-Q plot (0.95 by
#' default).
#' @param labels Logical argument. If \code{TRUE} labels the points outside
#' confidence interval limits.
#' @param theme The theme to use in the graphics. Default is \code{theme =
#' theme_waasb()}.
#' @param alpha The transparency of confidence band in the Q-Q plot. Must be a
#' number between 0 (opaque) and 1 (full transparency).
#' @param fill.hist The color to fill the histogram. Default is "gray".
#' @param col.hist The color of the border of the the histogram. Default is
#' "black".
#' @param col.point The color of the points in the graphic. Default is "black".
#' @param col.line The color of the lines in the graphic. Default is "red".
#' @param col.lab.out The color of the labels for the "outlying" points.
#' @param size.lab.out The size of the labels for the "outlying" points.
#' @param size.tex.lab The size of the text in axis text and labels.
#' @param size.shape The size of the shape in the plots.
#' @param bins The number of bins to use in the histogram. Default is 30.
#' @param which Which graphics should be plotted. Default is \code{which =
#' c(1:4)} that means that the first four graphics will be plotted.
#' @param mfrow Allow creating a simple multi-paneled plot. Must be a vector of
#' length 2. The first argument is the number of rows and the second the number
#' of columns. Default is \code{mfrow = c(2, 2)}.
#' @param ... Additional parameter for the function
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method autoplot WAAS.AMMI
#' @export
#' @templateVar fun autoplot.WAAS.AMMI
#' @template template-depr_fun
NULL

#' @templateVar old autoplot.WAAS.AMMI
#' @templateVar new autoplot.waas
#' @template template-depr_pkg
#'
#' @export
autoplot.WAAS.AMMI <- function(x, conf = 0.95, labels = FALSE, theme = theme_waasb(),
    alpha = 0.2, fill.hist = "gray", col.hist = "black", col.point = "black", col.line = "red",
    col.lab.out = "red", size.lab.out = 2.5, size.tex.lab = 10, size.shape = 1.5,
    bins = 30, which = c(1:4), mfrow = c(2, 2), ...) {

    df <- x$residuals
    df$id <- rownames(df)
    df <- data.frame(df[order(df$stdres), ])
    P <- ppoints(nrow(df))
    df$z <- qnorm(P)
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
    df$label <- ifelse(df$stdres > df$upper | df$stdres < df$lower, rownames(df),
        "")

    # residuals vs fitted

    p1 <- ggplot(df, aes(fitted, resid)) + geom_point(col = col.point) + geom_smooth(se = F,
        method = "loess", col = col.line) + geom_hline(yintercept = 0, linetype = 2,
        col = "gray") + labs(x = "Fitted values", y = "Residual") + ggtitle("Residuals vs fitted") +
        theme %+replace% theme(axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"), plot.title = element_text(size = size.tex.lab,
                hjust = 0, vjust = 1))
    if (labels != FALSE) {
        p1 <- p1 + geom_text(aes(label = label), size = size.lab.out, hjust = "inward",
            col = col.lab.out)
    } else {
        p1 <- p1
    }

    # normal qq
    p2 <- ggplot(df, aes(z, stdres)) + geom_point(col = col.point) + geom_abline(intercept = coef[1],
        slope = coef[2], col = col.line, size = 1) + geom_ribbon(aes_(ymin = ~lower,
        ymax = ~upper), alpha = alpha) + labs(x = "Theoretical quantiles", y = "Sample quantiles") +
        ggtitle("Normal Q-Q") + theme %+replace% theme(axis.text = element_text(size = size.tex.lab,
        colour = "black"), axis.title = element_text(size = size.tex.lab, colour = "black"),
        plot.title = element_text(size = size.tex.lab, hjust = 0, vjust = 1))
    if (labels != FALSE) {
        p2 <- p2 + geom_text(aes(label = label), size = size.lab.out, hjust = "inward",
            col = col.lab.out)
    } else {
        p2 <- p2
    }

    # scale-location
    p3 <- ggplot(df, aes(fitted, sqrt(abs(resid)))) + geom_point(col = col.point) +
        geom_smooth(se = F, method = "loess", col = col.line) + labs(x = "Fitted values",
        y = expression(sqrt("|Standardized residuals|"))) + ggtitle("Scale-location") +
        theme %+replace% theme(axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"), plot.title = element_text(size = size.tex.lab,
                hjust = 0, vjust = 1))
    if (labels != FALSE) {
        p3 <- p3 + geom_text(aes(label = label), size = size.lab.out, hjust = "inward",
            col = col.lab.out)
    } else {
        p3 <- p3
    }

    # Residuals vs Factor-levels
    p4 <- ggplot(df, aes(factors, stdres)) + geom_point(col = col.point) + geom_hline(yintercept = 0,
        linetype = 2, col = "gray") + labs(x = "Fitted values", y = "Standardized residuals") +
        ggtitle("Residuals vs factor-levels") + theme %+replace% theme(axis.text = element_text(size = size.tex.lab,
        colour = "black"), axis.title = element_text(size = size.tex.lab, colour = "black"),
        plot.title = element_text(size = size.tex.lab, hjust = 0, vjust = 1))
    if (labels != FALSE) {
        p4 <- p4 + geom_text(aes(label = label), size = size.lab.out, hjust = "inward",
            col = col.lab.out)
    } else {
        p4 <- p4
    }


    # Histogram of residuals
    p5 <- ggplot(df, aes(x = resid)) + geom_histogram(bins = bins, colour = col.hist,
        fill = fill.hist, aes(y = ..density..)) + stat_function(fun = dnorm, color = col.line,
        size = 1, args = list(mean = mean(df$resid), sd = sd(df$resid))) + labs(x = "Raw residuals",
        y = "Density") + ggtitle("Histogram of residuals") + theme

    # Residuals vs order
    p6 <- ggplot(df, aes(as.numeric(id), stdres, group = 1)) + geom_point(col = col.point) +
        geom_line() + geom_hline(yintercept = 0, linetype = 2, col = col.line) + labs(x = "Observation order",
        y = "Standardized Residuals") + ggtitle("Residuals vs observation order") +
        theme %+replace% theme(axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"), plot.title = element_text(size = size.tex.lab,
                hjust = 0, vjust = 1))



    p7 <- ggplot(df, aes(.fitted, Y)) + geom_point(col = col.point) + facet_wrap(~GEN) +
        geom_abline(intercept = 0, slope = 1, col = col.line) + labs(x = "Fitted values",
        y = "Observed values") + ggtitle("1:1 line plot") + theme %+replace% theme(axis.text = element_text(size = size.tex.lab,
        colour = "black"), axis.title = element_text(size = size.tex.lab, colour = "black"),
        plot.title = element_text(size = size.tex.lab, hjust = 0, vjust = 1), panel.spacing = unit(0,
            "cm"))

    plots <- list(p1, p2, p3, p4, p5, p6, p7)

    # making the plots
    grid::grid.newpage()

    if (prod(mfrow) > 1) {
        mypos <- expand.grid(1:mfrow[1], 1:mfrow[2])
        mypos <- mypos[with(mypos, order(Var1)), ]
        grid::pushViewport(viewport(layout = grid.layout(mfrow[1], mfrow[2])))
        formatter <- function(.) {
        }
    } else {
        mypos <- data.frame(matrix(1, length(which), 2))
        grid::pushViewport(viewport(layout = grid.layout(1, 1)))
        formatter <- function(.) {
            .dontcare <- readline("Hit <Return> to see next plot: ")
            grid::grid.newpage()
        }
    }

    j <- 1
    for (i in which) {
        formatter()
        print(plots[[i]], vp = viewport(layout.pos.row = mypos[j, ][1], layout.pos.col = mypos[j,
            ][2]))
        j <- j + 1
    }
}



