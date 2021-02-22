#' Split a data frame by factors
#' @description
#' `r badge('stable')`
#'
#' Split a data frame into subsets grouping by one or more factors.
#'
#' This function is used to split a data frame into a named list where each
#' element is a level of the grouping variable (or combination of grouping
#' variables).
#' @name split_factors
#' @param .data The data that will be split. Must contain at least one grouping
#'   variable.
#' @param x An object to check for class `split_factors`.
#' @param ... Comma-separated list of unquoted variable names that will be used
#'   to split the data.
#' @param keep_factors Should the grouping columns be kept?
#' @return A list where each element is a named level of the grouping factors.
#'   If more than one grouping variable is used, then each element is the
#'   combination of the grouping variables.
#' @details
#' * `split_factors()` Split a data frame by factors.
#' * `as.splict_factors()` coerce to an object of class `split_factors`
#' * `is.splict_factors()` check if an object is of class `split_factors`
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @importFrom rlang  eval_bare  expr set_names
#' @importFrom dplyr group_vars
#' @examples
#' \donttest{
#' library(metan)
#'
#' g1 <- split_factors(iris, Species)
#' g2 <- split_factors(data_ge, ENV, keep_factors = TRUE)
#'
#' spdata <- as.split_factors(iris)
#'
#' is.split_factors(spdata)
#' }
#'
split_factors <- function(.data, ..., keep_factors = FALSE) {
  grouped <- group_by(.data, ...)
  names <- eval_bare(expr(paste(!!!group_keys(grouped), sep = " | ")))
  gd <- grouped %>%
    group_split(.keep = keep_factors) %>%
    set_names(names)
  return(gd %>% set_class("split_factors"))
}
#' @name split_factors
#' @export
as.split_factors <- function(.data, keep_factors = FALSE) {
  grouped <- .data %>% group_by(across(where(is.factor)))
  names <- eval_bare(expr(paste(!!!group_keys(grouped), sep = " | ")))
  gd <- grouped %>%
    group_split(.keep = keep_factors) %>%
    set_names(names)
  return(gd %>% set_class("split_factors"))
}
#' @name split_factors
#' @export
is.split_factors <- function(x){
  if(has_class(x, "split_factors")){
    return(TRUE)
  } else{
    return(FALSE)
  }
}
