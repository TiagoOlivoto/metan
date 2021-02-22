#' @name lineplots
#' @title Fast way to create line plots
#' @description
#' `r badge('stable')`
#'
#' * `plot_lines()` Creates a line plot based on one quantitative factor
#' and one numeric variable.  It can be used to show the results of a one-way
#' trial with **quantitative treatments**.
#' * `plot_factlines()` Creates a line plot based on: one categorical and one quantitative factor and
#' one numeric variable.  It can be used to show the results of a
#' two-way trial with **qualitative-quantitative treatment structure**.
#' @param .data The data set
#' @param x,y The variables to be mapped to the `x` and `y` axes,
#'   respectively.
#' @param group The grouping variable. Valid for `plot_factlines()` only.
#' @param fit The polynomial degree to use. It must be between 1 (linear fit) to
#'   4 (fourth-order polynomial regression.). In `plot_factlines()`, if
#'   `fit` is a lenth 1 vector, i.e., 1, the fitted curves of all levels in
#'   `group` will be fitted with polynomial degree `fit`. To use a
#'   different polynomial degree for each level in `group`, use a numeric
#'   vector with the same length of the variable in `group`.
#' @param level The fonfidence level. Defaults to `0.05`.
#' @param confidence Display confidence interval around smooth? (TRUE by
#'   default)
#' @param xlab,ylab The labels of the axes x and y, respectively. Defaults to
#'   `NULL`.
#' @param n.dodge The number of rows that should be used to render the x labels.
#'   This is useful for displaying labels that would otherwise overlap.
#' @param check.overlap Silently remove overlapping labels, (recursively)
#'   prioritizing the first, last, and middle labels.
#' @param legend.position Valid argument for `plot_factlines`. The position
#'   of the legend. Defaults to 'bottom'.
#' @param grid Valid argument for `plot_factlines`. Logical argument. If
#'   `TRUE` then a grid will be created.
#' @param scales Valid argument for `plot_factlines`. If `grid = TRUE`
#'   scales controls how the scales are in the plot. Possible values are
#'   `'free'` (default), `'fixed'`, `'free_x'` or
#'   `'free_y'`.
#' @param col The colour to be used in the line plot and points.
#' @param alpha The alpha for the color in confidence band
#' @param size.shape The size for the shape in plot
#' @param size.line The size for the line in the plot
#' @param size.text The size of the text
#' @param fontfam The family of the font text.
#' @param plot_theme The graphical theme of the plot. Default is
#'   `plot_theme = theme_metan()`. For more details, see
#'   [ggplot2::theme()].
#' @return An object of class `gg, ggplot`.
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso [plot_bars()] and [plot_factbars()]
#' @examples
#' \donttest{
#' library(metan)
#' # One-way line plot
#' df1 <- data.frame(group = "A",
#'                   x = c(0, 100, 200, 300, 400),
#'                   y = c(3.2, 3.3, 4.0, 3.8, 3.4))
#' plot_lines(df1, x, y, fit = 2)
#'
#' # Two-way line plot
#' df2 <- data.frame(group = "B",
#'                   x = c(0, 100, 200, 300, 400),
#'                   y = c(3.2, 3.3, 3.7, 3.9, 4.1))
#' facts <- rbind(df1, df2)
#'
#' p1 <- plot_factlines(facts, x, y, group = group, fit = 1)
#' p2 <- plot_factlines(facts,
#'                      x = x,
#'                      y = y,
#'                      group = group,
#'                      fit = c(2, 1),
#'                      confidence = FALSE)
#' arrange_ggplot(p1, p2)
#'}
plot_lines <- function(.data,
                       x,
                       y,
                       fit,
                       level = 0.95,
                       confidence = TRUE,
                       xlab = NULL,
                       ylab = NULL,
                       n.dodge = 1,
                       check.overlap = FALSE,
                       col = "red",
                       alpha = 0.2,
                       size.shape = 1.5,
                       size.line = 1,
                       size.text = 12,
                       fontfam = "sans",
                       plot_theme = theme_metan()) {
  cl <- match.call()
  if (col == TRUE) {
    stop(paste0("The argument col = ", cl$col, " is invalid. Please, informe one valid color or FALSE to a black and white plot"))
  }
  data2 <- .data %>% select(x = !!enquo(x), y = !!enquo(y))
  if (fit == 1) {
    formula <- as.formula("y ~ x")
  }
  if (fit == 2) {
    formula <- as.formula("y ~ poly(x, 2)")
  }
  if (fit == 3) {
    formula <- as.formula("y ~ poly(x, 3)")
  }
  if (fit == 4) {
    formula <- as.formula("y ~ poly(x, 4)")
  }
  if (col == FALSE) {
    linetype <- 1
    p_smooth <- stat_smooth(method = "lm",
                            formula = formula,
                            data = data2,
                            level = level,
                            alpha = alpha,
                            col = "black",
                            se = confidence,
                            size = size.line)
  } else {
    linetype <- 1
    p_smooth <- stat_smooth(method = "lm",
                            formula = formula,
                            data = data2,
                            level = level,
                            alpha = alpha,
                            col = col,
                            se = confidence,
                            size = size.line)
  }
  if (is.null(ylab) == TRUE) {
    ylab <- cl$y
  } else {
    ylab <- ylab
  }
  if (is.null(xlab) == TRUE) {
    xlab <- cl$x
  } else {
    xlab <- xlab
  }
  if (col == FALSE) {
    p <- ggplot(data2, aes(x = x, y = y)) +
      geom_point(size = size.shape)
  } else {
    p <- ggplot(data2, aes(x = x, y = y)) +
      geom_point(size = size.shape, col = col)
  }
  p <- p +
    p_smooth +
    plot_theme %+replace%
    theme(axis.ticks.length = unit(0.2, "cm"),
          axis.text = element_text(size = size.text, family = fontfam, colour = "black"),
          axis.title = element_text(size = size.text, family = fontfam, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")) +
    labs(y = ylab, x = xlab) +
    scale_x_continuous(guide = guide_axis(n.dodge = n.dodge, check.overlap = check.overlap))
  return(p)
}

#' @name lineplots
#' @export
plot_factlines <- function(.data,
                           x,
                           y,
                           group,
                           fit,
                           level = 0.95,
                           confidence = TRUE,
                           xlab = NULL,
                           ylab = NULL,
                           n.dodge = 1,
                           check.overlap = FALSE,
                           legend.position = "bottom",
                           grid = FALSE,
                           scales = "free",
                           col = TRUE,
                           alpha = 0.2,
                           size.shape = 1.5,
                           size.line = 1,
                           size.text = 12,
                           fontfam = "sans",
                           plot_theme = theme_metan()) {
  if (length(fit) == 1 & grid == TRUE) {
    stop("When grid is TRUE the argument fit must have the same length of the grouping variable.")
  }
  if (max(fit) >= 5) {
    stop("The maximum polynomial degree is 4.")
  }
  cl <- match.call()
  data2 <- .data %>% select(factors = !!enquo(group), x = !!enquo(x),
                            y = !!enquo(y))
  group <- as.factor(data2$factors)
  p_smooth <- list()
  levels <- levels(group)
  if (length(fit) > 1) {
    for (i in 1:length(levels)) {
      levelname <- levels[i]
      mycond <- quote(group == levelname)
      if (fit[i] == 1) {
        formula <- as.formula("y ~ x")
      }
      if (fit[i] == 2) {
        formula <- as.formula("y ~ poly(x, 2)")
      }
      if (fit[i] == 3) {
        formula <- as.formula("y ~ poly(x, 3)")
      }
      if (fit[i] == 4) {
        formula <- as.formula("y ~ poly(x, 4)")
      }
      if (col == FALSE) {
        linetype <- i
        p_smooth[[paste(levels[i])]] <- stat_smooth(method = "lm",
                                                    formula = formula,
                                                    data = subset(data2, eval(mycond)),
                                                    level = level, linetype = linetype,
                                                    alpha = alpha,
                                                    col = "black",
                                                    size = size.line,
                                                    se = confidence)
      } else {
        linetype <- 1
        p_smooth[[paste(levels[i])]] <- stat_smooth(method = "lm",
                                                    formula = formula,
                                                    data = subset(data2, eval(mycond)),
                                                    level = level,
                                                    linetype = linetype,
                                                    alpha = alpha,
                                                    size = size.line,
                                                    se = confidence)
      }
    }
  } else {
    if (fit == 1) {
      formula <- as.formula("y ~ x")
    }
    if (fit == 2) {
      formula <- as.formula("y ~ poly(x, 2)")
    }
    if (fit == 3) {
      formula <- as.formula("y ~ poly(x, 3)")
    }
    if (fit == 4) {
      formula <- as.formula("y ~ poly(x, 4)")
    }
    if (col == TRUE) {
      p_smooth <- stat_smooth(method = "lm",
                              formula = formula,
                              data = data2,
                              level = level,
                              linetype = 1,
                              size = size.line,
                              se = confidence)
    } else {
      p_smooth <- stat_smooth(method = "lm",
                              formula = formula,
                              data = data2,
                              level = level,
                              linetype = 1,
                              col = "black",
                              size = size.line,
                              se = confidence)
    }
  }
  if (grid == TRUE) {
    legend.position <- "none"
  } else {
    legend.position <- legend.position
  }
  if (is.null(ylab) == TRUE) {
    ylab <- cl$y
  } else {
    ylab <- ylab
  }
  if (is.null(xlab) == TRUE) {
    xlab <- cl$x
  } else {
    xlab <- xlab
  }
  if (col == FALSE) {
    p <- ggplot(data2, aes(x = x, y = y)) +
      geom_point(aes(shape = factors), size = size.shape)
  } else {
    if (length(fit) > 1) {
      p <- ggplot(data2, aes(x = x, y = y, colour = factors)) +
        geom_point(aes(colour = factors), size = size.shape)
    } else {
      p <- ggplot(data2, aes(x = x, y = y)) +
        geom_point(aes(colour = factors), size = size.shape)
    }
  }
  p <- p +
    p_smooth +
    plot_theme %+replace%
    theme(axis.ticks.length = unit(0.2, "cm"),
          axis.text = element_text(size = size.text, family = fontfam, colour = "black"),
          axis.title = element_text(size = size.text, family = fontfam, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
          legend.title = element_blank(),
          legend.position = legend.position,
          legend.text = element_text(size = size.text, family = fontfam)) +
    labs(y = ylab, x = xlab) +
    scale_x_continuous(guide = guide_axis(n.dodge = n.dodge, check.overlap = check.overlap))
  if (grid == TRUE) {
    p <- p + facet_wrap(~factors, scales = scales)
  } else {
    p <- p
  }
  return(p)
}

