#' @title Encode variables to a specific format
#'
#' @description Function to quick encode columns to a specific format.
#' * \code{as_numeric()}: Encode columns to numeric using
#' \code{\link{as.numeric}()}.
#' * \code{as_integer()}: Encode columns to integer using
#' \code{\link{as.integer}()}.
#' * \code{as_logical()}: Encode columns to logical using
#' \code{\link{as.logical}()}.
#' * \code{as_character()}: Encode columns to character using
#' \code{\link{as.character}()}.
#' * \code{as_factor()}: Encode columns to factor using
#' \code{\link{as.factor}()}.
#'
#' @name utils_as
#' @param .data A data frame
#' @param ... <[`tidy-select`][dplyr_tidy_select]>. The variable(s) to encode to
#'   a format.
#' @return An object of the same class of \code{.data} with the variables in
#'   \code{...} encoded to the specified format.
#' @param .keep Allows you to control which columns from \code{.data} are
#'   retained in the output.
#' * \code{"all"} (default) retains all variables.
#' * \code{"used"} keeps any variables used to make new variables.
#' @param .pull Allows you to pull out the last column of the output. It is
#'   useful in combination with \code{.keep = "used"}. In this case, a vector
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
    res <- mutate(.data, across(c(...), as.numeric), .keep = .keep)
    if (.pull == TRUE){
      res <- pull(res)
    }
return(res)
}
#' @name utils_as
#' @export
as_integer <- function(.data, ..., .keep = "all", .pull = FALSE){
  res <- mutate(.data, across(c(...), as.integer), .keep = .keep)
  if (.pull == TRUE){
    res <- pull(res)
  }
  return(res)
}
#' @name utils_as
#' @export
as_logical <- function(.data, ..., .keep = "all", .pull = FALSE){
  res <- mutate(.data, across(c(...), as.logical), .keep = .keep)
  if (.pull == TRUE){
    res <- pull(res)
  }
  return(res)
}
#' @name utils_as
#' @export
as_character <-function(.data, ..., .keep = "all", .pull = FALSE){
  res <- mutate(.data, across(c(...), as.character), .keep = .keep)
  if (.pull == TRUE){
    res <- pull(res)
  }
  return(res)
}
#' @name utils_as
#' @export
to_factor <- function(.data, ...){
  warning("'to_factor()' is deprecated as of metan 1.9.0. Use 'as_factor()' instead.", call. = FALSE)
  return(mutate(.data, across(c(...), as.factor)))
}
#' @name utils_as
#' @export
as_factor <- function(.data, ..., .keep = "all", .pull = FALSE){
  res <- mutate(.data, across(c(...), as.factor), .keep = .keep)
  if (.pull == TRUE){
    res <- pull(res)
  }
  return(res)
}

