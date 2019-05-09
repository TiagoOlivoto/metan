#' Coerce to an object of class group_factors
#'
#' Functions to check if an object is of class \code{group_factors}
#'
#' A dataframe may be easily coerced to be split into named subsets based on
#' each combination of factors existing in the original dataframe. For example,
#' if the original data has two columns, namely ENV (four levels) and HIB (ten
#' levels), and ten numeric columns, then using \code{as.group_factors} will
#' split the data into 40 10-columns subsets, corresponding to each combination
#' of ENV x HIB.
#'
#' @param x The input data.
#' @importFrom methods is
#' @export
#' @examples
#'
#' library(METAAB)
#' spdata = as.group_factors(iris)
#' is.group_factors(spdata)
#'
is.group_factors = function(x){
  return((class(x) == "group_factors"))
}
