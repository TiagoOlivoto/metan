#' @name barplots
#' @title Fast way to create bar plots
#' @description
#' * \code{plot_bars()} Creates a bar plot based on one categorical variable and
#' one numeric variable.  It can be used to show the results of a one-way trial
#' with \strong{qualitative treatments}.
#' * \code{plot_factbars()} Creates a bar plot based on two categorical
#' variables and one numeric variable.  It can be used to show the results of a
#' two-way trial with \strong{qualitative-qualitative treatment structure}.
#' @param .data The data set.
#' @param ... Argument valid for \code{plot_factbars()}. A comma-separated list
#'   of unquoted variable names. Sets the two variables to be mapped to the
#'   \code{x} axis.
#' @param resp Argument valid for \code{plot_factbars()}. The response variable
#'   to be mapped to the y axis.
#' @param x,y Argument valid for \code{plot_bars()} The variables to be mapped
#'   to the \code{x} and \code{y} axes, respectively.
#' @param order Argument valid for \code{plot_bars()}. Controls the order of the
#'   factor in the \code{x} axis. Defaults to the order of the factors in
#'   \code{.data}. Use \code{order = "asce"} or \code{order = "desc"} to reorder
#'   the labels to ascending or descending order, respectively, based on the
#'   values of the variable \code{y}.
#' @param y.lim The range of y axis. Defaults to \code{NULL} (maximum and
#'   minimum values of the data set). New values can be inserted as \code{y.lim
#'   = c(y.min, y.max)}.
#' @param y.breaks The breaks to be plotted in the y-axis. Defaults to waiver().
#'   \code{authomatic breaks}. The same arguments than \code{x.breaks} can be
#'   used.
#' @param y.expand,y.contract A multiplication range expansion/contraction
#'   factor. \code{y.expand} expands the upper limit of the y escale, while
#'   \code{y.contract} contracts the lower limit of the y scale. By default
#'   \code{y.expand = 0.05} and \code{y.contract = 0} produces a plot without
#'   spacing in the lower y limit and an expansion in the upper y limit.
#' @param xlab,ylab The labels of the axes x and y, respectively. Defaults to
#'   \code{NULL}.
#' @param n.dodge The number of rows that should be used to render the x labels.
#'   This is useful for displaying labels that would otherwise overlap.
#' @param check.overlap Silently remove overlapping labels, (recursively)
#'   prioritizing the first, last, and middle labels.
#' @param lab.bar A vector of characters to show in each bar. Defaults to NULL.
#' @param lab.bar.hjust,lab.bar.vjust The horizontal and vertical adjust for the
#'   labels in the bar. Defaults to 0.5 and -0.5, respectively.
#' @param lab.bar.angle The angle for the labels in the plot. Defaults to 0. Use
#'   in combination with \code{lab.bar.hjust} and \code{lab.bar.vjust} to best
#'   fit the labels in the plot.
#' @param size.text.bar The size of the text in the bar labels.
#' @param values Logical argument. Shows the values in the plot bar?
#'   Defaults to \code{FALSE}
#' @param values.hjust,values.vjust The horizontal and vertical adjust
#'   for the values in the bar. Defaults to 0.5 and 1.5, respectively. If
#'   \code{values = TRUE} the values are shown bellow the error bar.
#' @param values.angle The angle for the labels in the plot. Defaults to 0.
#'   Use in combination with \code{values.hjust} and \code{values.vjust}
#'   to best fit the values in the plot bar.
#' @param values.digits The significant digits to show if \code{values
#'   = TRUE}. Defaults to \code{2}.
#' @param values.size The size of the text for values shown in the bars.
#'   Defaults to \code{3}.
#' @param lab.x.hjust,lab.x.vjust The horizontal and vertical adjust for the
#'   labels in the bar. Defaults to 0.5 and 1, respectively.
#' @param lab.x.angle The angle for the labels in x axis. Defaults to 0. Use
#'   in combination with \code{lab.x.hjust} and \code{lab.x.vjust} to best
#'   fit the labels in the axis.
#' @param errorbar Logical argument, set to TRUE. In this case, an error bar is
#'   shown.
#' @param stat.erbar The statistic to be shown in the errorbar. Must be one of
#'   the \code{stat.erbar = "se"} (standard error, default), \code{stat.erbar =
#'   "sd"} (standard deviation), or \code{stat.erbar = "ci"} (confidence
#'   interval), based on the confidence level in the argument \code{level}.
#' @param width.erbar The width of the error bar. Defaults to 25% of
#'   \code{width.bar}.
#' @param level The confidence level
#' @param invert Logical argument. If \code{TRUE}, rotate the plot in
#'   \code{plot_bars()} and invert the order of the factors in
#'   \code{plot_factbars()}.
#' @param color.bar,fill.bar Argument valid for \code{plot_bars()}. The color and
#'   fill values of the bars.
#' @param col Logical argument valid for \code{plot_factbars()}. If
#'   \code{FALSE}, a gray scale is used.
#' @param palette Argument valid for \code{plot_factbars()} The color palette to
#'   be used. For more details, see \code{?scale_colour_brewer}
#' @param width.bar The width of the bars in the graph. Defaults to 0.9.
#'   Possible values are in the range 0-1.
#' @param legend.position The position of the legend in the plot.
#' @param size.line The size of the line in the bars. Default to \code{0.5}.
#' @param size.text The size of the text. Default to \code{12}.
#' @param fontfam The family of the font text. Defaults to \code{"sans"}.
#' @param na.rm Should 'NA' values be removed to compute the statistics?
#'   Defaults to true
#' @param verbose Logical argument. If TRUE a tibble containing the mean, N,
#'   standard deviation, standard error of mean and confidence interval is
#'   returned.
#' @param plot_theme The graphical theme of the plot. Default is
#'   \code{plot_theme = theme_metan()}. For more details, see
#'   \code{\link[ggplot2]{theme}}.
#' @return An object of class \code{gg, ggplot}.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @seealso \code{\link{plot_lines}}, \code{\link{plot_factlines}}
#'
#' @examples
#' \donttest{
#' library(metan)
#' # two categorical variables
#' plot_factbars(data_ge2,
#'               GEN,
#'               ENV,
#'               resp = PH)
#'
#' # one categorical variable
#' p1 <- plot_bars(data_g, GEN, PH)
#' p2 <- plot_bars(data_g, GEN, PH,
#'                 n.dodge = 2, # two rows for x labels
#'                 y.expand = 0.1, # expand y scale
#'                 y.contract = -0.75, # contract the lower limit
#'                 errorbar = FALSE, # remove errorbar
#'                 color.bar = "red", # color of bars
#'                 fill.bar = alpha_color("cyan", 75), # create a transparent color
#'                 lab.bar = letters[1:13]) # add labels
#' arrange_ggplot(p1, p2)
#'}
plot_bars <- function(.data,
                      x,
                      y,
                      order = NULL,
                      y.lim = NULL,
                      y.breaks = waiver(),
                      y.expand = 0.05,
                      y.contract = 0,
                      xlab = NULL,
                      ylab = NULL,
                      n.dodge = 1,
                      check.overlap = FALSE,
                      color.bar = "black",
                      fill.bar = "gray",
                      lab.bar = NULL,
                      lab.bar.hjust = 0.5,
                      lab.bar.vjust = -0.5,
                      lab.bar.angle = 0,
                      size.text.bar = 5,
                      values = FALSE,
                      values.hjust = 0.5,
                      values.vjust = 1.5,
                      values.angle = 0,
                      values.digits = 2,
                      values.size = 4,
                      lab.x.hjust = 0.5,
                      lab.x.vjust = 1,
                      lab.x.angle = 0,
                      errorbar = TRUE,
                      stat.erbar = "se",
                      width.erbar = NULL,
                      level = 0.95,
                      invert = FALSE,
                      width.bar = 0.9,
                      size.line = 0.5,
                      size.text = 12,
                      fontfam = "sans",
                      na.rm = TRUE,
                      verbose = FALSE,
                      plot_theme = theme_metan()) {
  if(!missing(order) && !order %in% c("asce", "desc")){
    stop("Argument order must be one of 'asce' or 'desc'", call. = FALSE)
  }
  width.erbar <- ifelse(missing(width.erbar), width.bar/4, width.erbar)
  cl <- match.call()
  datac <-
    .data %>%
    as_factor({{x}}) %>%
    select({{x}}, {{y}}) %>%
    group_by({{x}}) %>%
    desc_stat({{y}}, stats = c("n, mean, sd.amo, ci, se"), level = level)
  if(errorbar == TRUE){
    if(stat.erbar == "ci"){
      datac %<>% add_cols(max = mean + ci,
                          min = mean - ci)
    }
    if(stat.erbar == "sd"){
      datac %<>% add_cols(max = mean + sd.amo,
                          min = mean - sd.amo)
    }
    if(stat.erbar == "se"){
      datac %<>% add_cols(max = mean + se,
                          min = mean - se)
    }
  } else{
    datac %<>% add_cols(max = mean,
                        min = mean)
  }
  ylab <- ifelse(is.null(ylab), cl$y, ylab)
  xlab <- ifelse(is.null(xlab), cl$x, xlab)
  if(!missing(order)){
    if(order == "asce"){
      p <- ggplot(datac, aes(reorder({{x}}, mean), mean))
    }
    if(order == "desc"){
      p <- ggplot(datac, aes(reorder({{x}}, -mean), mean))
    }
  } else{
    p <- ggplot(datac, aes(x = {{x}}, y = mean))
  }
  p <- p +
    geom_bar(stat = "identity",
             width = width.bar,
             color = color.bar,
             size = size.line,
             fill = fill.bar)
  if (errorbar == TRUE) {
    if (stat.erbar == "ci") {
      p <- p + geom_errorbar(aes(ymin = mean - ci,
                                 ymax = mean + ci),
                             size = size.line,
                             width = width.erbar)
    }
    if (stat.erbar == "sd") {
      p <- p + geom_errorbar(aes(ymin = mean - sd.amo,
                                 ymax = mean + sd.amo),
                             size = size.line,
                             width = width.erbar)
    }
    if (stat.erbar == "se") {
      p <- p + geom_errorbar(aes(ymin = mean - se,
                                 ymax = mean + se),
                             size = size.line,
                             width = width.erbar)
    }
  }
  if (!missing(lab.bar)) {
    if (length(lab.bar) > 1 & length(lab.bar) != nrow(datac)) {
      stop("The labels must be either length 1 or the same as the levels of ",
           paste(xlab), " (", nrow(datac), ")", call. = FALSE)
    }
    p <- p + geom_text(aes(label = lab.bar, y = max),
                       vjust = lab.bar.vjust,
                       hjust = lab.bar.hjust,
                       size = size.text.bar,
                       family = fontfam,
                       angle = lab.bar.angle)
  }
  if(values == TRUE){
    p <- p + geom_text(aes(label = round(mean, values.digits), y = min),
                       vjust = values.vjust,
                       hjust = values.hjust,
                       size = values.size,
                       family = fontfam,
                       angle = values.angle)
  }
  p <- p +
    plot_theme %+replace%
    theme(axis.ticks.length = unit(0.2, "cm"),
          axis.text = element_text(size = size.text, family = fontfam, colour = "black"),
          axis.text.x = element_text(angle = lab.x.angle, vjust = lab.x.vjust, hjust = lab.x.hjust),
          axis.title = element_text(size = size.text, family = fontfam, colour = "black"),
          axis.ticks = element_line(colour = "black", size = size.line),
          plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")) +
    theme(panel.border = element_rect(size = size.line)) +
    labs(y = ylab, x = xlab) +
    scale_y_continuous(limits = y.lim,
                       breaks = y.breaks,
                       expand = expansion(mult = c(y.contract, y.expand))) +
    scale_x_discrete(guide = guide_axis(n.dodge = n.dodge, check.overlap = check.overlap))

  if (verbose == TRUE) {
    print(datac)
  }
  if(invert == TRUE){
    return(p + coord_flip())
  }
  return(p)
}



#' @name barplots
#' @export

plot_factbars <- function(.data,
                          ...,
                          resp,
                          y.lim = NULL,
                          y.breaks = waiver(),
                          y.expand = 0.05,
                          y.contract = 0,
                          xlab = NULL,
                          ylab = NULL,
                          n.dodge = 1,
                          check.overlap = FALSE,
                          lab.bar = NULL,
                          lab.bar.hjust = 0.5,
                          lab.bar.vjust = -0.5,
                          lab.bar.angle = 0,
                          size.text.bar = 5,
                          values = FALSE,
                          values.hjust = 0.5,
                          values.vjust = 1.5,
                          values.angle = 0,
                          values.digits = 2,
                          values.size = 4,
                          lab.x.hjust = 0.5,
                          lab.x.vjust = 1,
                          lab.x.angle = 0,
                          errorbar = TRUE,
                          stat.erbar = "se",
                          width.erbar = NULL,
                          level = 0.95,
                          invert = FALSE,
                          col = TRUE,
                          palette = "Spectral",
                          width.bar = 0.9,
                          legend.position = "bottom",
                          size.line = 0.5,
                          size.text = 12,
                          fontfam = "sans",
                          na.rm = TRUE,
                          verbose = FALSE,
                          plot_theme = theme_metan()) {
  width.erbar <- ifelse(missing(width.erbar), width.bar/4, width.erbar)
  cl <- match.call()
  datac <-
    .data %>%
    mutate(across(c(...), as.factor)) %>%
    select(..., Y = {{resp}}) %>%
    group_by(...) %>%
    summarise(N = n(),
              mean_var = mean(Y, na.rm = na.rm),
              sd = sd(Y, na.rm = na.rm), se = sd/sqrt(n()),
              ci = se * qt(level/2 + 0.5, n() - 1),
              .groups = "drop")
  nam <- names(select(.data, ...))
  if(errorbar == TRUE){
    if(stat.erbar == "ci"){
      datac %<>% add_cols(max = mean_var + ci,
                          min = mean_var - ci)
    }
    if(stat.erbar == "sd"){
      datac %<>% add_cols(max = mean_var + sd.amo,
                          min = mean_var - sd.amo)
    }
    if(stat.erbar == "se"){
      datac %<>% add_cols(max = mean_var + se,
                          min = mean_var - se)
    }
  } else{
    datac %<>% add_cols(max = mean_var,
                        min = mean_var)
  }
  if (length(nam) > 1) {
    names(datac) <- c("x", "y", "N", "mean_var", "sd", "se", "ci", "max", "min")
  } else {
    names(datac) <- c("x", "N", "mean_var", "sd", "se", "ci", "max", "min")
  }
  if (is.null(ylab) == TRUE) {
    ylab <- cl$resp
  } else {
    ylab <- ylab
  }
  if (invert == FALSE) {
    if (is.null(xlab) == TRUE) {
      xlab <- nam[1]
    } else {
      xlab <- xlab
    }
  } else {
    if (is.null(xlab) == TRUE) {
      xlab <- nam[2]
    } else {
      xlab <- xlab
    }
  }
  pd <- position_dodge(width.bar)
  if (length(nam) > 1) {
    if (invert == FALSE) {
      p <- ggplot(data = datac, aes(x = x, y = mean_var, fill = y)) +
        geom_bar(aes(fill = y),
                 colour = "black",
                 stat = "identity",
                 position = position_dodge(),
                 width = width.bar,
                 size = size.line) +
        scale_fill_brewer(palette = palette)
    } else {
      p <- ggplot(data = datac, aes(x = y, y = mean_var, fill = x)) +
        geom_bar(aes(fill = x),
                 colour = "black",
                 stat = "identity",
                 position = position_dodge(),
                 width = width.bar,
                 size = size.line) +
        scale_fill_brewer(palette = palette)
    }
  } else {
    p <- ggplot(data = datac, aes(x = x, y = mean_var)) +
      geom_bar(stat = "identity",
               position = position_dodge(),
               width = width.bar,
               size = size.line)
  }
  if (col == FALSE) {
    p <- p + scale_fill_grey(start = 0, end = 0.9)
  } else {
    p <- p
  }
  if (errorbar == TRUE) {
    if (stat.erbar == "ci") {
      p <- p + geom_errorbar(aes(ymin = mean_var - ci,
                                 ymax = mean_var + ci),
                             width = width.erbar,
                             size = size.line,
                             position = pd)
    }
    if (stat.erbar == "sd") {
      p <- p + geom_errorbar(aes(ymin = mean_var - sd,
                                 ymax = mean_var + sd),
                             width = width.erbar,
                             size = size.line,
                             position = pd)
    }
    if (stat.erbar == "se") {
      p <- p + geom_errorbar(aes(ymin = mean_var - se,
                                 ymax = mean_var + se),
                             width = width.erbar,
                             size = size.line,
                             position = pd)
    }
  }
  if (!missing(lab.bar)) {
    if (length(lab.bar) > 1 & length(lab.bar) != nrow(datac)) {
      stop("The labels must be either length 1 or the same as the levels of ",
           paste(quos(...)), " (", nrow(datac), ")")
    }
    p <- p + geom_text(aes(label = lab.bar, y = max),
                       position = pd,
                       vjust = lab.bar.vjust,
                       hjust = lab.bar.hjust,
                       size = size.text.bar,
                       family = fontfam,
                       angle = lab.bar.angle)
  }
  if(values == TRUE){
    p <- p + geom_text(aes(label = round(mean_var, values.digits), y = min),
                       position = pd,
                       vjust = values.vjust,
                       hjust = values.hjust,
                       size = values.size,
                       family = fontfam,
                       angle = values.angle)
  }
  p <- p +
    plot_theme %+replace%
    theme(axis.ticks.length = unit(0.2, "cm"),
          axis.text = element_text(size = size.text, family = fontfam, colour = "black"),
          axis.text.x = element_text(angle = lab.x.angle, vjust = lab.x.vjust, hjust = lab.x.hjust),
          axis.title = element_text(size = size.text, family = fontfam, colour = "black"),
          axis.ticks = element_line(colour = "black", size = size.line),
          plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
          legend.title = element_blank(),
          legend.position = legend.position,
          legend.text = element_text(size = size.text, family = fontfam)) +
    theme(panel.border = element_rect(size = size.line)) +
    labs(y = ylab, x = xlab) +
    scale_y_continuous(limits = y.lim,
                       breaks = y.breaks,
                       expand = expansion(mult = c(y.contract, y.expand))) +
    scale_x_discrete(guide = guide_axis(n.dodge = n.dodge, check.overlap = check.overlap))
  if (verbose == TRUE) {
    print(datac)
  }
  return(p)
}

