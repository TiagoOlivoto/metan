#' Arrange separate ggplots into the same graphic
#'
#' This is a wraper function around \code{\link[patchwork]{wrap_plots}()} and
#' \code{\link[patchwork]{plot_annotation}()} to arrange ggplot2 objects.
#' @param ... multiple \code{ggplot}s or a list containing \code{ggplot}
#'   objects.
#' @param nrow,ncol The number of rows and columns, respectively.
#' @param widths,heights The relative widths and heights of each column and row
#'   in the grid. Will get repeated to match the dimensions of the grid.
#' @param guides A string specifying how guides should be treated in the layout.
#'   Defaults to `'auto'`. Other possible values are `'keep'` and `'collect'`.
#'   In this case, will collect guides below to the given nesting level,
#'   removing duplicates.
#' @param design Specification of the location of areas in the layout.
#' @param legend.position The position of the legends in the plot if
#'   \code{guides = "collect"} Default to `'bottom'`.
#' @param title,subtitle,caption Text strings to use for the various plot
#' annotations.
#' @param tag_levels A character vector defining the enumeration format to use
#' at each level. Possible values are `'a'` for lowercase letters, `'A'` for
#' uppercase letters, `'1'` for numbers, `'i'` for lowercase Roman numerals, and
#' `'I'` for uppercase Roman numerals. It can also be a list containing
#' character vectors defining arbitrary tag sequences. If any element in the
#' list is a scalar and one of `'a'`, `'A'`, `'1'`, `'i`, or `'I'`, this level
#' will be expanded to the expected sequence.
#' @param tag_prefix,tag_suffix Strings that should appear before or after the
#' tag.
#' @param tag_sep A separator between different tag levels.
#' @param theme A ggplot theme specification to use for the plot. Only elements
#'   related to the titles as well as plot margin and background is used.
#' @return A `patchwork` object
#' @import patchwork
#' @export
#'
#' @examples
#' \donttest{
#' library(ggplot2)
#' library(metan)
#' p1 <- ggplot(mtcars, aes(wt, mpg)) +
#'       geom_point()
#'
#' p2 <- ggplot(mpg, aes(class, hwy)) +
#'              geom_boxplot()
#'
#' # Default plot
#' arrange_ggplot(p1, p2)
#'
#' # Insert plot annotation, titles and subtitles
#' arrange_ggplot(p1, p2,
#'                ncol = 1,
#'                tag_levels = list(c("(P1)", "(P2)")),
#'                title = "My grouped ggplot",
#'                subtitle = "Made with arrange_ggplot()",
#'                caption = "P1 = scatter plot\nP2 = boxplot",
#'                theme = theme(plot.title = element_text(size = 20,
#'                                                        face = "bold"),
#'                              plot.subtitle = element_text(size = 10,
#'                                                           face = "italic"),
#'                              plot.caption  = element_text(size = 10,
#'                                                           face = "italic")))
#' }
#'
arrange_ggplot <- function(...,
                           nrow = NULL,
                           ncol = NULL,
                           widths = NULL,
                           heights = NULL,
                           guides = NULL,
                           design = NULL,
                           legend.position = "bottom",
                           title = NULL,
                           subtitle = NULL,
                           caption = NULL,
                           tag_levels = NULL,
                           tag_prefix = NULL,
                           tag_suffix = NULL,
                           tag_sep = NULL,
                           theme = NULL) {
  p <-
  wrap_plots(...,
             ncol = ncol,
             nrow = nrow,
             widths = widths,
             heights = heights,
             design = design,
             guides = guides) +
    plot_annotation(title = title,
                    subtitle = subtitle,
                    caption = caption,
                    tag_levels = tag_levels,
                    tag_prefix = tag_prefix,
                    tag_suffix = tag_suffix,
                    tag_sep = tag_sep,
                    theme = theme)
  if(!missing(guides)){
    p <-
      p +
      plot_layout(guides = guides) &
      theme(legend.position = legend.position)
  }
  return(p)
}
