#' @title Fast way to create a bar plot
#' @description Create a bar plot based on categorical (one or two) variables
#' and one numeric variable.
#' @param .data The data set
#' @param ... A comma-separated list of unquoted variable names. Must be up to two
#' variables.
#' @param resp The response variable
#' @param y.expand A multiplication factor to expand the y axis.. Defaults to 1.
#' @param y.breaks The breaks to be plotted in the y-axis. Defaults to waiver().
#' \code{authomatic breaks}. The same arguments than \code{x.breaks} can be
#' used.
#' @param xlab The x label
#' @param ylab The y label
#' @param lab.bar A vector of characters to show in each bar. Defaults to NULL.
#' @param lab.bar.vjust The vertical adjust for the labels in the bar. Defaults to -0.2
#' @param size.text.bar The size of the text in the bar labels.
#' @param errorbar Logical argument, set to TRUE. In this case, an error bar is shown.
#' @param stat.erbar The statistic to be shown in the errorbar. Must be one of the 'se'
#' (standard error, default), 'sd' (standard deviation), or 'ci' confidence interval, based
#' on the confidence level
#' @param width.erbar The width of the error bar.
#' @param level The confidence level
#' @param invert Logical argument. If \code{TRUE}, the order of the factors entered
#' in changes in the graph
#' @param col Logical argument. If \code{FALSE}, a gray scale is used.
#' @param palette The color palette to be used. For more details, see
#' \code{?scale_colour_brewer}
#' @param width.bar The width of the bars in the graph. Defaults to 0.9
#' possible values [0-1].
#' @param lab.x.angle The angle of the caption text. Default is 0.
#' @param lab.x.hjust The horizontal adjustment of the caption text. Defaults to 0.5.
#' @param lab.x.vjust The vertical adjustment of the caption text. Defaults to 1.
#' Use this argument to adjust the text when the angle of the text is different from 0.
#' @param legend.position The legend position.
#' @param size.text The size of the text
#' @param fontfam The family of the font text
#' @param na.rm Should 'NA' values be removed to compute the statistics?
#' Defaults to true
#' @param verbose Logical argument. If TRUE a tibble containing the mean, N,
#' standard deviation, standard error of mean and confidence interval is returned.
#' @export
#' @seealso \code{\link{plot_lines}}, \code{\link{plot_factlines}}
#'
#' @examples
#' library(metan)
#' plot_factbars(data_ge2,
#'               GEN,
#'               ENV,
#'               resp = PH,
#'               palette = 'Greys')
plot_factbars <- function(.data, ..., resp, y.expand = 1, y.breaks = waiver(),
                          xlab = NULL, ylab = NULL, lab.bar = NULL, lab.bar.vjust = -0.2,
                          size.text.bar = 5, errorbar = TRUE, stat.erbar = "se", width.erbar = 0.3,
                          level = 0.95, invert = FALSE, col = TRUE, palette = "Spectral",
                          width.bar = 0.9, lab.x.angle = 0, lab.x.hjust = 0.5, lab.x.vjust = 1,
                          legend.position = "bottom", size.text = 12, fontfam = "sans",
                          na.rm = TRUE, verbose = FALSE) {
  cl <- match.call()
  datac <- .data %>% mutate_at(quos(...), as.factor) %>% select(...,
                                                                Y = !!enquo(resp)) %>% group_by(...) %>% summarise(N = n(),
                                                                                                                   mean_var = mean(Y, na.rm = na.rm), sd = sd(Y, na.rm = na.rm),
                                                                                                                   se = sd/sqrt(n()), ci = se * qt(level/2 + 0.5, n() -
                                                                                                                                                     1))
  nam <- names(select(.data, ...))
  if (length(nam) > 1) {
    names(datac) <- c("x", "y", "N", "mean_var", "sd", "se",
                      "ci")
  } else {
    names(datac) <- c("x", "N", "mean_var", "sd", "se", "ci")
  }
  if (is.null(ylab) == T) {
    ylab <- cl$resp
  } else {
    ylab <- ylab
  }
  if (invert == FALSE) {
    if (is.null(xlab) == T) {
      xlab <- nam[1]
    } else {
      xlab <- xlab
    }
  } else {
    if (is.null(xlab) == T) {
      xlab <- nam[2]
    } else {
      xlab <- xlab
    }
  }
  y.lim <- c(0, (max(datac$mean_var) + max(datac$ci)) * y.expand)
  pd <- ggplot2::position_dodge(width.bar)
  if (length(nam) > 1) {
    if (invert == FALSE) {
      p <- ggplot2::ggplot(data = datac, aes(x = x, y = mean_var,
                                             fill = y)) + geom_bar(aes(fill = y), colour = "black",
                                                                   stat = "identity", position = position_dodge(),
                                                                   width = width.bar) + scale_fill_brewer(palette = palette)
    } else {
      p <- ggplot2::ggplot(data = datac, aes(x = y, y = mean_var,
                                             fill = x)) + geom_bar(aes(fill = x), colour = "black",
                                                                   stat = "identity", position = position_dodge(),
                                                                   width = width.bar) + scale_fill_brewer(palette = palette)
    }
  } else {
    p <- ggplot2::ggplot(data = datac, aes(x = x, y = mean_var)) +
      geom_bar(stat = "identity", position = position_dodge(),
               width = width.bar)
  }
  p <- p + ggplot2::theme_bw()
  if (col == FALSE) {
    p <- p + scale_fill_grey(start = 0, end = 0.9)
  } else {
    p <- p
  }
  if (errorbar == TRUE) {
    if (stat.erbar == "ci") {
      p <- p + geom_errorbar(aes(ymin = mean_var - ci,
                                 ymax = mean_var + ci), width = width.erbar, position = pd)
    }
    if (stat.erbar == "sd") {
      p <- p + geom_errorbar(aes(ymin = mean_var - sd,
                                 ymax = mean_var + sd), width = width.erbar, position = pd)
    }
    if (stat.erbar == "se") {
      p <- p + geom_errorbar(aes(ymin = mean_var - se,
                                 ymax = mean_var + se), width = width.erbar, position = pd)
    }
  }
  if (!missing(lab.bar)) {
    if (length(lab.bar) > 1 & length(lab.bar) != nrow(datac)) {
      stop("The labels must be either length 1 or the same as the levels of ",
           paste(quos(...)), " (", nrow(datac), ")")
    }
    p <- p + geom_text(aes(label = lab.bar), position = pd,
                       vjust = lab.bar.vjust, size = size.text.bar)
  }
  p <- p + ggplot2::theme(axis.ticks.length = unit(0.2, "cm"),
                          axis.text = element_text(size = size.text, family = fontfam,
                                                   colour = "black"), axis.title = element_text(size = size.text,
                                                                                                family = fontfam, colour = "black"), axis.text.x = element_text(angle = lab.x.angle,
                                                                                                                                                                hjust = lab.x.hjust, vjust = lab.x.vjust, size = size.text,
                                                                                                                                                                colour = "black"), axis.ticks = element_line(colour = "black"),
                          plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"), legend.title = element_blank(),
                          legend.position = legend.position, legend.text = element_text(size = size.text,
                                                                                        family = fontfam), panel.border = element_rect(colour = "black",
                                                                                                                                       fill = NA, size = 1), panel.grid.major.x = element_blank(),
                          panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(),
                          panel.grid.minor.y = element_blank()) + ggplot2::labs(y = ylab,
                                                                                x = xlab) + scale_y_continuous(limits = y.lim, breaks = y.breaks,
                                                                                                               expand = expand_scale(mult = c(0, 0.1)))
  if (verbose == TRUE) {
    print(datac)
  }
  return(p)
}
