#' @title Encode variables to a specific format
#'
#' @description
#' `r badge('stable')`
#'
#' Function to quick encode vector or columns to a specific format.
#' * `as_numeric()`: Encode columns to numeric using
#' [as.numeric()].
#' * `as_integer()`: Encode columns to integer using
#' [as.integer()].
#' * `as_logical()`: Encode columns to logical using
#' [as.logical()].
#' * `as_character()`: Encode columns to character using
#' [as.character()].
#' * `as_factor()`: Encode columns to factor using
#' [as.factor()].
#' @name utils_as
#' @param .data A data frame or a vector.
#' @param ... <[`tidy-select`][dplyr_tidy_select]>. If `.data` is a data
#'   frame, then `...` are the variable(s) to encode to a format.
#' @return An object of the same class of `.data` with the variables in
#'   `...` encoded to the specified format.
#' @param .keep Allows you to control which columns from `.data` are
#'   retained in the output.
#' * `"all"` (default) retains all variables.
#' * `"used"` keeps any variables used to make new variables.
#' @param .pull Allows you to pull out the last column of the output. It is
#'   useful in combination with `.keep = "used"`. In this case, a vector
#'   will be created with the used column.
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(metan)
#' library(tibble)
#' df <-
#'   tibble(y = rnorm(5),
#'          x1 = c(1:5),
#'          x2 = c(TRUE, TRUE, FALSE, FALSE, FALSE),
#'          x3 = letters[1:5],
#'          x4 = as.factor(x3))
#' df
#'
#' # Convert y to integer
#' as_integer(df, y)
#' as_integer(df$y)
#'
#' # convert x3 to factor
#' as_factor(df, x3)
#'
#' # Convert all columns to character
#' as_character(df, everything())
#'
#' # Convert x2 to numeric and coerce to a vector
#' as_numeric(df, x2, .keep = "used", .pull = TRUE)
#' }
#'
as_numeric <- function(.data, ..., .keep = "all", .pull = FALSE){
  if (has_class(.data, c("data.frame","tbl_df", "data.table"))){
    res <- mutate(.data, across(c(...), as.numeric), .keep = .keep)
    if (.pull == TRUE){
      res <- pull(res)
    }
    return(res)
  } else{
    return(as.numeric(.data))
  }
}
#' @name utils_as
#' @export
as_integer <- function(.data, ..., .keep = "all", .pull = FALSE){
  if (has_class(.data, c("data.frame","tbl_df", "data.table"))){
    res <- mutate(.data, across(c(...), as.integer), .keep = .keep)
    if (.pull == TRUE){
      res <- pull(res)
    }
    return(res)
  } else{
    return(as.integer(.data))
  }
}
#' @name utils_as
#' @export
as_logical <- function(.data, ..., .keep = "all", .pull = FALSE){
  if (has_class(.data, c("data.frame","tbl_df", "data.table"))){
    res <- mutate(.data, across(c(...), as.logical), .keep = .keep)
    if (.pull == TRUE){
      res <- pull(res)
    }
    return(res)
  } else{
    return(as.logical(.data))
  }
}
#' @name utils_as
#' @export
as_character <-function(.data, ..., .keep = "all", .pull = FALSE){
  if (has_class(.data, c("data.frame","tbl_df", "data.table"))){
    res <- mutate(.data, across(c(...), as.character), .keep = .keep)
    if (.pull == TRUE){
      res <- pull(res)
    }
    return(res)
  } else{
    return(as.character(.data))
  }
}
#' @name utils_as
#' @export
as_factor <- function(.data, ..., .keep = "all", .pull = FALSE){
  if (has_class(.data, c("data.frame","tbl_df", "data.table"))){
    res <- mutate(.data, across(c(...), as.factor), .keep = .keep)
    if (.pull == TRUE){
      res <- pull(res)
    }
    return(res)
  } else{
    return(as.factor(.data))
  }
}
