#' @title Utilities for handling with NA values
#' @name utils_na
#' @param .data A data frame or tibble
#' @param ... Variables to replace \code{NAs}. If \code{...} is null then all
#'   variables with \code{NA} will be replaced. It must be a single variable
#'   name or a comma-separated list of unquoted variables names. Select helpers
#'   are also allowed.
#' @param prop The proportion (percentage) of \code{NA} values to generate in \code{.data}.
#' @param replace The value used for replacement. Defaults to \code{0}. Use
#'   \code{replace = "colmeans"} to replace missing values with colum means.
#' @param verbose Logical argument. If \code{TRUE} (default) shows in console
#'   the rows or columns deleted.
#' @export
#' @return A data frame with rows or columns with \code{NA} values deleted.
#' @description
#' * \code{has_na()}: Check for \code{NA} values in the data and return a logical value.
#' * \code{random_na()}: Generate random \code{NA} values in a two-way table
#' based on a desired proportion.
#' * \code{remove_rows_na()}: Remove rows with \code{NA} values.
#' * \code{remove_cols_na()}: Remove columns with \code{NA} values.
#' * \code{replace_na()} Replace missing values
#' @md
#' @examples
#' \donttest{
#' library(metan)
#' data_with_na <- data_g
#' data_with_na[c(1, 5, 10), c(3:5, 10:15)] <- NA
#' data_with_na
#' has_na(data_with_na)
#' remove_cols_na(data_with_na)
#' remove_rows_na(data_with_na)
#' replace_na(data_with_na)
#' }
remove_rows_na <- function(.data, verbose = TRUE){
row_with_na <- which(complete.cases(.data) ==  FALSE)
if(verbose == TRUE){
warning("Row(s) ", paste(row_with_na, collapse = ", "), " with NA values deleted.", call. = FALSE)
}
return(na.omit(.data))
}

#' @name utils_na
#' @export
remove_cols_na <- function(.data, verbose = TRUE){
  cols_with_na <- names(which(sapply(.data, anyNA)))
  if(verbose == TRUE){
    warning("Column(s) ", paste(cols_with_na, collapse = ", "), " with NA values deleted.", call. = FALSE)
  }
  return(select(.data, -cols_with_na))
}

#' @name utils_na
#' @export
has_na <- function(.data){
   any(complete.cases(.data) ==  FALSE)
}
#' @name utils_na
#' @export
replace_na <- function(.data, ..., replace = 0){
  if (has_class(.data, c("data.frame","tbl_df", "data.table"))){
    if(has_rownames(.data)){
      rnames <- rownames(.data)
    }
  if(missing(...)){
    cols_with_na <- names(which(sapply(.data, anyNA)))
    df <- mutate_at(.data, vars(cols_with_na), ~replace(., is.na(.), replace))
    if(replace == "colmeans"){
    df <- mutate_at(.data, vars(cols_with_na), ~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))
    }
  } else{
    df <- mutate_at(.data, vars(...), ~replace(., is.na(.), replace))
    df <- mutate_at(.data, vars(...), ~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))
  }
  } else{
    df <- replace(.data, is.na(.data), replace)
  }
  if(has_rownames(.data)){
    rownames(df) <- rnames
  }
  return(df)
}
#' @name utils_na
#' @export
random_na <- function(.data, prop){
  if( prop < 1 | prop > 100){
    stop("Argument prob must have a 1-100 interval.")
  }
  .data <- as.matrix(.data)
  cells <- length(.data)
  miss <- (cells * (prop / 100))
  coord_miss <- sample(1:cells, miss, replace = FALSE)
  .data[coord_miss] <- NA
  return(as.matrix(.data))
}
