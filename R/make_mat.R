make_mat <- function(.data, row, col, value) {
  data = data.frame(.data %>% dplyr::select(!!dplyr::enquo(row),
                                            !!dplyr::enquo(col),
                                            !!dplyr::enquo(value)))
  return(data.frame(tapply(data[, 3], data[, c(1, 2)], mean)))
}
