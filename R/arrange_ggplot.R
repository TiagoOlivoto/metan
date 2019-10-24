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
  dots <- list(...)
  n <- length(dots)
  if(missing(nrow) & missing(ncol)) {
    nrow = floor(n/2)
    ncol = ceiling(n/nrow)
  }
  if(missing(nrow)){
    nrow = ceiling(n/ncol)
  }
  if(missing(ncol)){
    ncol = ceiling(n/nrow)
  }
  grid.newpage()
  viewport(layout = grid.layout(nrow, ncol)) %>%
    pushViewport()
  p <- 1
  for(i in seq(1, nrow)){
    tbrow <- i
    for(j in seq(1, ncol)){
      tble <- p
      if(p > n) break
      vpl <- function(x, y){
        viewport(layout.pos.row = x, layout.pos.col = y)
      }
      print(dots[[tble]], vp = vpl(tbrow, j))
      p <- p + 1
    }
  }
}
