#' Make a two-way table
#' @description
#' `r badge('stable')`
#'
#' This function help users to easily make a two-way table from a "long format"
#' data.
#'
#'
#' @param .data The dataset. Must contains at least two categorical columns.
#' @param row The column of data in which the mean of each level will
#' correspond to **one line** in the output.
#' @param col The column of data in which the mean of each level will
#' correspond to **one column** in the output.
#' @param value The column of data that contains the values to fill the two-way
#' table.
#' @param fun The function to apply. Defaults to `mean`, i.e., the two-way table
#' will show the mean values for each genotype-environment combination. Other R base functions
#' such as `max`, `min`, `sd`, `var`, or an own function that return
#'  a single numeric value can be used.
#'
#' @return A two-way table with the argument `row` in the rows, `col`
#'   in the columns, filled by the argument `value`.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' matrix <- data_ge %>% make_mat(row = GEN, col = ENV, val = GY)
#' matrix
#'
#' # standart error of mean
#'
#' data_ge %>% make_mat(GEN, ENV, GY, sem)
#'}

make_mat <- function(.data, row, col, value, fun = mean) {
  data <- .data %>%
    select({{row}},
           {{col}},
           {{value}}) %>%
    group_by({{row}}, {{col}}) %>%
    summarise(across(where(is.numeric), fun, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = {{col}}, values_from = {{value}})
  data %<>% column_to_rownames(var = names(data[1]))
  return(data)
}
