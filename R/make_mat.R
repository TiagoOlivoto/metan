#' Make a two-way table based on categorical and numerical arguments
#'
#' This function help users to easely make a two-way table from a "long format"
#' data.
#'
#'
#' @param .data The dataset. Must contains at least two categorical columns.
#' @param row The column of data in which the mean of each level will
#' correspond to \strong{one line} in the output.
#' @param col The column of data in which the mean of each level will
#' correspond to \strong{one column} in the output.
#' @param value The column of data that contains the values to fill the two-way
#' table.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'
#' library(METAAB)
#' library(dplyr)
#' matrix = data_ge %>% make_mat(row = GEN, col = ENV, val = GY)
#' matrix
#'
make_mat <- function(.data, row, col, value) {
  data = data.frame(.data %>% dplyr::select(!!dplyr::enquo(row),
                                            !!dplyr::enquo(col),
                                            !!dplyr::enquo(value)))
  return(data.frame(tapply(data[, 3], data[, c(1, 2)], mean)))
}
