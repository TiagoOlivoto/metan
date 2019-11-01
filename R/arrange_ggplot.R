#' Arrange multiple ggplot2 graphics in a single image window
#'
#' This is a helper function to arrange ggplot2 objects in the metan package
#' @param ... An object of class \code{gg}
#' @param nrow The number of rows
#' @param ncol The number of columns
#'
#' @return None.
#' @export
#'
#' @examples
#' library(ggplot2)
#' p1 <- ggplot(mtcars, aes(wt, mpg)) +
#'       geom_point()
#'
#' p2 <- ggplot(mpg, aes(class, hwy)) +
#'              geom_boxplot()
#'
#' arrange_ggplot(p1, p2)
#'
arrange_ggplot <- function(..., nrow = NULL, ncol = NULL) {
plot_grid(plotlist = ..., nrow = nrow, ncol = ncol)
}

