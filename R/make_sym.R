#' Make a symmetric matrix on a triangular matrix
#'
#' This function help users to easily make a symmetric matrix using a lower or
#' an upper triangular matrix.
#'
#'
#' @param .matrix The upper or lower triangular matrix.
#' @param make The triangular to built. Default is \code{"upper"}. In this
#' case, a symetric matrix will be built based on the values of a lower
#' triangular values.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'
#' m <- matrix(NA,4,4)
#' m[lower.tri(m)] <- 1:16
#' m
#' make_sym(m)
#'
make_sym <- function(.matrix, make = "upper") {
  if(make == "upper"){
    .matrix[upper.tri(.matrix)] <- t(.matrix)[upper.tri(.matrix)]
  }
  if (make == "lower"){
    .matrix[lower.tri(.matrix)] <- t(.matrix)[lower.tri(.matrix)]
  }
  return(.matrix)
}
