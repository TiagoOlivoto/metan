#' @title Utilities for set operations for many sets
#'
#' @description
#' `r badge('stable')`
#'
#'  Provides alternative function to [base::union()], [base::intersect()], and
#'  [base::setdiff()].
#' * `set_union()`: Returns the union of the sets in `...`.
#' * `set_intersect()`: Returns the intersect of the sets in `...`.
#' * `set_difference()`: Returns the difference of the sets in `...`.
#'
#' @name utils_sets
#' @param ... A list or a comma-separated list of vectors in the same class. If
#'   vector contains duplicates they will be discarded. If the list doesn't have
#'   names the sets will be named as `"set_1"`, "`Set_2"`, `"Set_3"` and so on.
#'   If vectors are given in `...`, the set names will be named with the names
#'   of the objects provided. Data frames are also allowed, provided that common
#'   column names exist.
#' @param pairs Returns the pairwise unions of the sets? Defaults to `FALSE`.
#' @return A vector showing the desired operation of the sets. If `pairs =
#'   TRUE`, returns a list showing the pairwise operation of the sets.
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(metan)
#' (A <- letters[1:4])
#' (B <- letters[2:5])
#' (C <- letters[3:7])
#'
#' set_union(A, B)
#' set_intersect(A, B, C)
#' set_difference(B, C)
#'
#' # Operations with data frames
#' # Add a row id for better understanding
#' sets <- data_ge %>% add_row_id()
#'
#' set_1 <- sets[1:5,]
#' set_2 <- sets[2:6,]
#' set_3 <- sets[3:7,]
#'
#' set_intersect(set_1, set_2, set_3)
#' set_difference(set_2, set_3)
#' set_union(set_1, set_2, set_3)
#' }
#'

set_intersect <- function(..., pairs = FALSE){
  set_helper(..., pairs = pairs, fun = intersect)
}
#' @name utils_sets
#' @export
set_union <- function(..., pairs = FALSE){
  set_helper(..., pairs = pairs, fun = union)
}
#' @name utils_sets
#' @export
set_difference <- function(..., pairs = FALSE){
  set_helper(..., pairs = pairs, fun = setdiff)
}

set_helper <- function(..., pairs = FALSE, fun = intersect){
  set_name <-
    as.character(
      sapply(quos(...), function(x){
        rlang::quo_get_expr(x)
      })
    )
  if(has_df_in_list(list(...))){
    sets <- list(...)
  } else{
    if(is.list(c(...))){
      sets <- as.list(...)
      if(is.null(names(sets))){
        names(sets) <- paste("set", seq_len(length(sets)), sep = "_")
      }
    }else{
      sets <- list(...)
      if(is.null(names(sets))){
        names(sets) <- set_name
      }
    }
  }
  if (length(sets) <= 1) {
    stop("The list must contain at least 2 vectors.", call. = FALSE)
  }
  if (!all_df_in_list(sets) & length(unique(sapply(sets, class))) != 1) {
    stop("Vectors must be in the same class.", call. = FALSE)
  }
  if(is.null(names(sets))){
    names(sets) <- paste("set", seq_len(length(sets)), sep = "_")
  }
  if(pairs == FALSE){
    reduce(sets, fun)
  } else{
    pairs = combn(length(sets), 2)
    inter = vector("list", ncol(pairs))
    for(i in 1:ncol(pairs)) {
      inter[[i]] = fun(sets[[pairs[1, i]]],
                       sets[[pairs[2, i]]])
      names(inter)[i] = paste(names(sets[pairs[1, i]]),
                              names(sets[pairs[2, i]]),
                              sep = "...")
    }
    inter
  }
}
