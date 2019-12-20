#' @title Utilities for manipulating data
#' @name utilities
#' @param x A numeric value or a numeric vector.
#' @param .data A data frame
#' @param ... Variables to round. If no variable is informed, all the numeric
#'   variables from \code{data} are used.
#' @param digits The number of significant figures.
#' @param cols A quoted variable name.
#'
#' @description Utilities for manipulating data.
#' * \code{round_value()}: Round a numeric value or a numeric vector to significant figures.
#'
#' * \code{round_column()}: Round a selected column or a whole data frame to significant figures.
#'
#' * \code{column_exists()}: Check if a column exists in a data frame.
#' @md
round_value <- function(x, digits = 2){
  sapply(x, function(x, digits){
      if(is.na(x)){
        x
      } else if (abs(x) > 10^-digits){
        round(x, digits)
      } else {
        signif(x, digits)
      }
    }, digits
  )
}

#' @name utilities
#' @export
round_column <- function(.data, ...,  digits = 2){
  if (missing(...)){
    .data %<>% dplyr::mutate_if(is.numeric, round_value, digits = digits)
  } else{
    .data %<>% dplyr::mutate_at(vars(...), round_value, digits = digits)
  }
  return(.data)
}

#' @name utilities
#' @export
column_exists <-function(.data, cols){
  if(length(setdiff(cols, colnames(.data))) != 0){
    stop("Can't find the variables ",
         paste(cols, collapse = ", "),
         " in the .data", call. = FALSE)
  } else{
    message("Variable ", paste(cols, collapse = ", "), " is in the data.")
  }
}


#' @title Means by one or more factors
#' @description Computes the mean values for all numeric variables of a data by a factor
#' @param .data A data frame
#' @param ... One or more categorical variables for grouping the data.
#' @export
#' @examples
#' library(metan)
#' data_ge2 %>%
#' means_by(ENV)
means_by <- function(.data, ...){
  .data %>%
    group_by(...) %>%
    summarise_if(is.numeric, mean)
}
