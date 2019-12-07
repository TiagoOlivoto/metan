#' Harmonic mean
#'
#' Helper function to compute the harmonic mean. The harmonic mean is the
#' reciprocal of the arithmetic mean of the reciprocals.
#' @param x A numeric vector or a data frame.
#' @param ... Variables to compute the harmonic mean. If no variable is informed
#'   and \code{x} is a data frame, all the numeric variables will be used.
#' @param na.rm A logical value indicating whether \code{NA} values should be
#'   stripped before computation proceeds.
#' @seealso \code{\link{gm_mean}}
#' @note Not useful if there are elements with values \code{<=0}.
#' @return The harmonic mean(s) of \code{x}. If \code{x} is a numeric vector, the
#'   function returns a numeric value. If a data frame is used then a numeric
#'   vector with the harmonic mean for each variable is returned.
#' @export
#'
#' @examples
#' num <- c(1:10, 50)
#' hm_mean(num)
#'
#' num_df <- make_mat(data_ge, ENV, GEN, GY)
#' hm_mean(num_df)
hm_mean <- function(x, ..., na.rm = TRUE) {
  if(is.null(nrow(x))){
    1 / mean(1 / x, na.rm = na.rm)
  } else{
    if(missing(...)){
      df <- select_if(x, is.numeric)
    } else{
      df <- select(x, ...) %>%
        select_if(is.numeric)
    }
    1 / (apply( 1 / df, 2, mean, na.rm = na.rm))
  }
}
