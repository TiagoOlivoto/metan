#' Visualization of a correlation matrix
#' @description
#' `r badge('stable')`
#'
#' Graphical and numerical visualization of a correlation matrix
#'
#'
#' @param .data The data. Should, preferentially, contain numeric variables
#'   only. If `.data` has factor-columns, these columns will be deleted
#'   with a warning message.
#' @param ... Variables to use in the correlation. If no variable is informed
#'   all the numeric variables from `.data` are used.
#' @param col.by A categorical variable to map the color of the points by.
#'   Defaults to `NULL`.
#' @param upper The visualization method for the upper triangular correlation
#'   matrix. Must be one of `'corr'` (numeric values), `'scatter'`
#'   (the scatterplot for each pairwise combination), or `NULL` to set a
#'   blank diagonal.
#' @param lower The visualization method for the lower triangular correlation
#'   matrix. Must be one of `'corr'` (numeric values), `'scatter'`
#'   (the scatterplot for each pairwise combination), or `NULL` to set a
#'   blank diagonal.
#' @param decimal.mark The decimal mark. Defaults to `"."`.
#' @param axis.labels Should the axis labels be shown in the plot? Set to
#'   `FALSE`.
#' @param show.labels.in Where to show the axis labels. Defaults to "show"
#'   bottom and left. Use "diag" to show the labels on the diagonal. In this
#'   case, the diagonal layer (boxplot, density or histogram) will be
#'   overwritten.
#' @param size.axis.label The size of the text for axis labels if
#'   `axis.labels = TRUE`. Defaults to 12.
#' @param size.varnames The size of the text for variable names. Defaults to 12.
#' @param col.varnames The color of the text for variable names. Defaults to
#'   "black".
#' @param diag Should the diagonal be shown?
#' @param diag.type The type of plot to show in the diagonal if `diag
#'   TRUE`. It must be one of the 'histogram' (to show an histogram), 'density'
#'   to show the Kernel density, or 'boxplot' (to show a boxplot).
#' @param  bins The number of bins, Defaults to 20.
#' @param col.diag If `diag = TRUE` then `diagcol` is the color for
#'   the distribution. Set to gray.
#' @param alpha.diag Alpha-transparency scale (0-1) to make the diagonal plot
#'   transparent. 0 = fully transparent; 1 = full color. Set to 0.15
#' @param col.up.panel,col.lw.panel,col.dia.panel The color for the upper,
#'   lower, and diagonal panels, respectively. Set to 'gray'.
#' @param prob The probability of error. Significant correlations will be
#'   highlighted with '*', '**', and '***' (0.05, 0.01, and 0.001,
#'   respectively). Scatterplots with significant correlations may be
#'   color-highlighted.
#' @param col.sign The color that will highlight the significant correlations.
#'   Set to 'green'.
#' @param alpha.sign Alpha-transparency scale (0-1) to make the plot area
#'   transparent. 0 = fully transparent; 1 = full color. Set to 0.15
#' @param lab.position The position that the labels will appear. Set to
#'   `'tr'`, i.e., the legends will appear in the top and right of the
#'   plot. Other allowed options are `'tl'` (top and left), `'br'`
#'   (bottom and right), `'bl'` (bottom and left).
#' @param progress `NULL` (default) for a progress bar in interactive
#'   sessions with more than 15 plots, `TRUE` for a progress bar,
#'   `FALSE` for no progress bar.
#' @param smooth Should a linear smooth line be shown in the scatterplots? Set
#'   to `FALSE`.
#' @param col.smooth The color for the smooth line.
#' @param confint Should a confidence band be shown with the smooth line? Set to
#'   `TRUE`.
#' @param size.point The size of the points in the plot. Set to `0.5`.
#' @param shape.point The shape of the point, set to `1`.
#' @param alpha.point Alpha-transparency scale (0-1) to make the points
#'   transparent. 0 = fully transparent; 1 = full color. Set to 0.7
#' @param fill.point The color to fill the points. Valid argument if points are
#'   between 21 and 25.
#' @param col.point The color for the edge of the point, set to `black`.
#' @param size.line The size of the line (smooth and diagonal).
#' @param minsize The size of the letter that will represent the smallest
#'   correlation coefficient.
#' @param maxsize The size of the letter that will represent the largest
#'   correlation coefficient.
#' @param pan.spacing The space between the panels. Set to 0.15.
#' @param digits The number of digits to show in the plot.
#' @param export Logical argument. If `TRUE`, then the plot is exported to
#'   the current directory.
#' @param file.type The format of the file if `export = TRUE`.  Set to
#'   `'pdf'`. Other possible values are `*.tiff` using `file.type
#'   = 'tiff'`.
#' @param file.name The name of the plot when exported. Set to `NULL`,
#'   i.e., automatically.
#' @param width The width of the plot, set to `8`.
#' @param height The height of the plot, set to `7`.
#' @param resolution The resolution of the plot if `file.type = 'tiff'` is
#'   used. Set to `300` (300 dpi).
#' @return An object of class `gg, ggmatrix`.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' dataset <- data_ge2 %>% select_cols(1:7)
#'
#' # Default plot setting
#' corr_plot(dataset)
#'
#' # Chosing variables to be correlated
#' corr_plot(dataset, PH, EH, EL)
#'
#'
#' # Axis labels, similar to the function pairs()
#' # Gray scale
#' corr_plot(dataset, PH, EH, EL,
#'           shape.point = 19,
#'           size.point = 2,
#'           alpha.point = 0.5,
#'           alpha.diag = 0,
#'           pan.spacing = 0,
#'           col.sign = 'gray',
#'           alpha.sign = 0.3,
#'           axis.labels = TRUE)
#'
#' corr_plot(dataset, PH, EH, EL,
#'           prob = 0.01,
#'           shape.point = 21,
#'           col.point = 'black',
#'           fill.point = 'orange',
#'           size.point = 2,
#'           alpha.point = 0.6,
#'           maxsize = 4,
#'           minsize = 2,
#'           smooth = TRUE,
#'           size.line = 1,
#'           col.smooth = 'black',
#'           col.sign = 'cyan',
#'           col.up.panel = 'black',
#'           col.lw.panel = 'black',
#'           col.dia.panel = 'black',
#'           pan.spacing = 0,
#'           lab.position = 'tl')
#'}
corr_plot <- function(.data, ...,
                      col.by = NULL,
                      upper = "corr",
                      lower = "scatter",
                      decimal.mark = ".",
                      axis.labels = FALSE,
                      show.labels.in = "show",
                      size.axis.label = 12,
                      size.varnames = 12,
                      col.varnames = "black",
                      diag = TRUE,
                      diag.type = "histogram",
                      bins = 20,
                      col.diag = "gray",
                      alpha.diag = 1,
                      col.up.panel = "gray",
                      col.lw.panel = "gray",
                      col.dia.panel = "gray",
                      prob = 0.05,
                      col.sign = "green",
                      alpha.sign = 0.15,
                      lab.position = "tr",
                      progress = NULL,
                      smooth = FALSE,
                      col.smooth = "red",
                      confint = TRUE,
                      size.point = 1,
                      shape.point = 19,
                      alpha.point = 0.7,
                      fill.point = NULL,
                      col.point = "black",
                      size.line = 0.5,
                      minsize = 2,
                      maxsize = 3,
                      pan.spacing = 0.15,
                      digits = 2,
                      export = FALSE,
                      file.type = "pdf",
                      file.name = NULL,
                      width = 8,
                      height = 7,
                      resolution = 300) {
  if(!show.labels.in %in% c("show", "internal", "none")){
    stop("The argument 'show.labels.in' must be one of the 'show', 'internal', or 'none'. ")
  }
  if (!lab.position %in% c("tr", "tl", "br", "bl")) {
    stop("The argument 'lab.position' must be one of the 'tr', 'tl', 'br', or 'bl'.")
  }
  if (!diag.type %in% c("histogram", "density", "boxplot")) {
    stop("The argument 'diag.type' must be one of the 'boxplot', 'histogram' or 'density'.")
  }
  if (!is.null(upper)) {
    if (!upper %in% c("corr", "scatter", NULL)) {
      stop("The argument 'upper' must be one of the 'corr', 'scatter' or 'NULL'.")
    }
  }
  if (!is.null(lower)) {
    if (!lower %in% c("corr", "scatter", NULL)) {
      stop("The argument 'lower' must be one of the 'corr', 'scatter' or 'NULL'.")
    }
  }
  col.by.test <- missing(col.by)
  if(!col.by.test && !missing(fill.point) || !col.by.test && !missing(col.point)){
    message("Arguments 'fill.point' and 'col.point' overwritten by 'col.by'")
  }
  if (missing(...)) {
    data <- select(.data, where(is.numeric), {{col.by}})
    if(has_na(data)){
      data <- remove_rows_na(data)
      has_text_in_num(data)
    }
  }
  if (!missing(...)) {
    data <- select(.data, ..., {{col.by}})
    if(has_na(data)){
      data <- remove_rows_na(data)
      has_text_in_num(data)
    }
  }
  w <- c(21:25)
  if (is.null(fill.point) == TRUE && any(w == shape.point) && col.by.test) {
    stop("If 'shape.point' is a value between 21 and 25, you must provide a color to fill the shape using the argument 'fill.point.'")
  }
  if (!is.null(fill.point) && !shape.point %in% w) {
    stop("If 'fill.point' is informed, then declare 'shape.point' between 21 and 25.", call. = FALSE)
  }
  my_custom_cor <- function(data, mapping, color = I("black"),
                            sizeRange = c(minsize, maxsize), ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    ct <- cor.test(x, y)
    sig <- symnum(ct$p.value, corr = FALSE, na = FALSE, cutpoints = c(0,
                                                                      0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**",
                                                                                                              "*", ".", " "))
    r <- unname(ct$estimate)
    rt <- format(r, digits = digits, nsmall = digits, decimal.mark = decimal.mark,  scientific = FALSE)[1]
    cex <- max(sizeRange)
    percent_of_range <- function(percent, range) {
      percent * diff(range) + min(range, na.rm = TRUE)
    }
    GGally::ggally_text(label =  rt,
                        mapping = aes(),
                        xP = 0.5, yP = 0.5, size = I(percent_of_range(cex *
                                                                        abs(r), sizeRange)), color = color, ...) + geom_text(aes_string(x = 0.8,
                                                                                                                                        y = 0.8), label = sig, size = I(cex), color = color,
                                                                                                                             ...) + theme_classic() + theme(panel.background = ggplot2::element_rect(color = col.up.panel),
                                                                                                                                                            axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                                                                                                                                                            axis.text.y = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank())
  }
  my_custom_smooth <- function(data, mapping, ...) {
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    ct <- cor.test(x, y)
    r <- unname(ct$p.value)
    rt <- format(r, digits = digits, nsmall = digits)[1]
    tt <- as.character(rt)
    p <- ggplot2::ggplot(data = data, mapping = mapping)
    if (!is.null(fill.point) && col.by.test) {
      p <- p + geom_point(color = I(col.point),
                          fill = fill.point,
                          shape = shape.point,
                          size = size.point,
                          alpha = alpha.point)
    } else {
      if(col.by.test){
        p <- p + geom_point(color = col.point,
                            shape = shape.point,
                            size = size.point,
                            alpha = alpha.point)
      } else{
        p <- p + geom_point(aes(color = {{col.by}},
                                fill = {{col.by}}),
                            shape = shape.point,
                            size = size.point,
                            alpha = alpha.point)
      }
    }
    p <- p + theme_classic() + theme(panel.background = ggplot2::element_rect(fill = "white",
                                                                              color = col.lw.panel))
    if (smooth == TRUE) {
      p <- p + geom_smooth(method = "lm",
                           formula = y ~ x,
                           se = confint,
                           size = size.line,
                           color = col.smooth, ...)
    }
    if (r < prob) {
      p <- p + theme(panel.background = ggplot2::element_rect(fill = ggplot2::alpha(col.sign,
                                                                                    alpha.sign)))
    }
    p
  }
  ggally_mysmooth <- function(data, mapping, ...) {

    dia <- ggplot(data = data, mapping = mapping)
    if (diag.type == "density") {
      dia <- dia + geom_density(fill = ggplot2::alpha(col.diag, alpha.diag), col = "black")
    }
    if (diag.type == "boxplot") {
      y <- GGally::eval_data_col(data, mapping$x)
      x = mean(y)
      dia <- dia + stat_boxplot(aes(y = y, x = x, group = 1),
                                geom ='errorbar',
                                color = transparent_color(),
                                size = size.line,
                                width = 0.9) +
        stat_boxplot(aes(y = y, x = x, group = 1),
                     geom ='errorbar',
                     size = size.line,
                     width = 0.2) +
        ggplot2::geom_boxplot(aes(y = y, x = x, group = 1),
                              fill = ggplot2::alpha(col.diag,
                                                    alpha.diag),
                              width = 0.4,
                              size = size.line,
                              col = "black")

    }
    if (diag.type == "histogram") {
      dia <- dia + geom_histogram(fill = ggplot2::alpha(col.diag, alpha.diag),
                                  col = "black",
                                  size = size.line,
                                  bins = bins)
    }
    dia <- dia + theme_classic() +
      theme(panel.background = element_rect(fill = ggplot2::alpha("white", 1), color = col.dia.panel))
  }
  if (diag == TRUE) {
    diag <- list(continuous = ggally_mysmooth)
  } else {
    diag <- NULL
  }
  if (!is.null(upper)) {
    if (upper == "corr") {
      upper <- list(continuous = my_custom_cor)
    }
    if (upper == "scatter") {
      upper <- list(continuous = my_custom_smooth)
    }
  }
  if (is.null(upper)) {
    upper <- NULL
  }
  if (!is.null(lower)) {
    if (lower == "corr") {
      lower <- list(continuous = my_custom_cor)
    }
    if (lower == "scatter") {
      lower <- list(continuous = my_custom_smooth)
    }
  }
  if (is.null(lower)) {
    lower <- NULL
  }
  if (lab.position == "tr") {
    switch <- NULL
  }
  if (lab.position == "tl") {
    switch <- "y"
  }
  if (lab.position == "br") {
    switch <- "x"
  }
  if (lab.position == "bl") {
    switch <- "both"
  }
  if (axis.labels == TRUE) {
    axis.labels <- show.labels.in
  } else {
    axis.labels <- "none"
  }
  p1 <- GGally::ggpairs(data,
                        columns =  which( unlist(lapply(data, is.numeric))),
                        upper = upper, lower = lower,
                        switch = switch, diag = diag, progress = progress, axisLabels = axis.labels)+
    theme(panel.spacing = grid::unit(pan.spacing, "lines"),
          axis.text = element_text(size = size.axis.label, color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          strip.text = element_text(size = size.varnames, colour = col.varnames))
  if (export == FALSE) {
    return(p1)
  } else if (file.type == "pdf") {
    if (is.null(file.name)) {
      pdf("Scatterplot Correlation.pdf", width = width,
          height = height)
    } else pdf(paste0(file.name, ".pdf"), width = width, height = height)
    print(p1)
    dev.off()
  }
  if (file.type == "tiff") {
    if (is.null(file.name)) {
      tiff(filename = "Scatterplot Correlation.tiff", width = width,
           height = height, units = "in", compression = "lzw",
           res = resolution)
    } else tiff(filename = paste0(file.name, ".tiff"), width = width,
                height = height, units = "in", compression = "lzw",
                res = resolution)
    print(p1)
    dev.off()
  }
}

