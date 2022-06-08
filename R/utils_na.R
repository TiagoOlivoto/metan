#' @title Utilities for handling with NA and zero values
#'
#' @description
#' `r badge('stable')`
#'
#'   NAs and zeros can increase the noise in multi-environment trial analysis.
#'   This collection of functions will make it easier to deal with them.
#' @name utils_na_zero
#' @param .data A data frame.
#' @param ... Variables to fill  `NAs` in `fill_na()`, replace
#'   `NAs` in `replace_na()` or zeros in `replace_zero()`. If
#'   `...` is null then all variables in `.data` will be evaluated. It
#'   must be a single variable name or a comma-separated list of unquoted
#'   variables names. Select helpers are also allowed.
#' @param direction Direction in which to fill missing values. Currently either
#'   "down" (the default), "up", "downup" (i.e. first down and then up) or
#'   "updown" (first up and then down).
#' @param prop The proportion (percentage) of `NA` values to generate in
#'   `.data`.
#' @param replacement The value used for replacement. Defaults to `0`. Other
#'   possible values are Use `"colmean"`, `"colmin"`, and `"colmax"` to replace
#'   missing values with column mean, minimum and maximum values, respectively.
#' @param verbose Logical argument. If `TRUE` (default) shows in console
#'   the rows or columns deleted.
#' @export
#' @return A data frame with rows or columns with `NA` values deleted.
#' @description
#' * `fill_na()`: Fills `NA` in selected columns using the next or
#' previous entry.
#' * `has_na(), has_zero()`: Check for `NAs` and `0s` in the
#' data and return a logical value.
#' * `prop_na()` returns the proportion of `NAs` in each column of a data frame.
#' * `random_na()`: Generate random `NA` values in a two-way table
#' based on a desired proportion.
#' * `remove_cols_na()`, `remove_rows_na()`: Remove columns and rows that
#' contains at least one `NA` value.
#' * `remove_cols_all_na()`, `remove_rows_all_na()`: Remove columns and rows
#' where all values are `NAs`.
#' * `remove_cols_zero()`, `remove_rows_zero()`: Remove columns and rows that
#' contains at least one `0` value, respectively.
#' * `select_cols_na(), select_cols_zero()`: Select columns with `NAs`
#' and `0s`, respectively.
#' * `select_rows_na(), select_rows_zero()`: Select rows with `NAs`
#' and `0s`, respectively.
#' * `replace_na(), replace_zero()`: Replace `NAs` and `0s`,
#' respectively, with a `replacement` value.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(metan)
#' data_naz <- iris %>%
#'               group_by(Species) %>%
#'               doo(~head(., n = 3)) %>%
#'               as_character(Species)
#' data_naz
#' data_naz[c(2:3, 6, 8), c(1:2, 4, 5)] <- NA
#' data_naz[c(2, 7, 9), c(2, 3, 4)] <- 0
#' has_na(data_naz)
#' has_zero(data_naz)
#'
#' # Fill NA values of column GEN
#' fill_na(data_naz, Species)
#'
#' # Remove columns
#' remove_cols_na(data_naz)
#' remove_cols_zero(data_naz)
#' remove_rows_na(data_naz)
#' remove_rows_zero(data_naz)
#'
#' # Select columns
#' select_cols_na(data_naz)
#' select_cols_zero(data_naz)
#' select_rows_na(data_naz)
#' select_rows_zero(data_naz)
#'
#' # Replace values
#' replace_na(data_naz)
#' replace_zero(data_naz)
#' }
#'

# Dealing with NAs

#' @importFrom tidyr fill
fill_na <- function(.data, ..., direction = "down"){
  if(!direction %in% c("down", "up", "downup", "updown")){
    stop("Argument 'direction' must be one of 'down', 'up', 'downup', or 'updown'. ")
  }
  if(missing(...)){
    return(fill(.data, everything(), .direction = direction))
  } else{
    return(fill(.data, ..., .direction = direction))
  }
}
#' @name utils_na_zero
#' @export
has_na <- function(.data){
  any(complete.cases(.data) ==  FALSE)
}
#' @name utils_na_zero
#' @export
prop_na <- function(.data, ...){
  if(missing(...)){
  df <- .data %>% select(everything())
  } else{
  df <- .data %>% select(...)
  }
  df <- apply(df, 2, function(x){
    round(length(which(is.na(x))) / length(x), digits = 4)
  }) %>%
    as.data.frame() %>%
    rownames_to_column("variable") %>%
    set_names(c("variable", "prop"))
  return(df)
}

#' @name utils_na_zero
#' @export
remove_rows_na <- function(.data, verbose = TRUE){
row_with_na <- which(complete.cases(.data) ==  FALSE)
if(verbose == TRUE){
warning("Row(s) ", paste(row_with_na, collapse = ", "), " with NA values deleted.", call. = FALSE)
}
return(na.omit(.data))
}

#' @name utils_na_zero
#' @export
remove_rows_all_na <- function(.data, verbose = TRUE){
  row_with_na <- which(apply(.data, 1, function(x){
    all(is.na(x))
  }))
  if(verbose == TRUE){
    warning("Row(s) ", paste(row_with_na, collapse = ", "), " with all NA values deleted.", call. = FALSE)
  }
  return(.data[-row_with_na, ])
}

#' @name utils_na_zero
#' @export
remove_cols_na <- function(.data, verbose = TRUE){
  cols_with_na <- names(which(sapply(.data, anyNA)))
  if(verbose == TRUE){
    warning("Column(s) ", paste(cols_with_na, collapse = ", "), " with NA values deleted.", call. = FALSE)
  }
  return(select(.data, -cols_with_na))
}

#' @name utils_na_zero
#' @export
remove_cols_all_na <- function(.data, verbose = TRUE){
  cols_with_na <- which(sapply(.data, function(x){
    all(is.na(x))
  }))
  if(verbose == TRUE){
    warning("Column(s) ", paste(names(cols_with_na), collapse = ", "), " with all NA values deleted.", call. = FALSE)
  }
  return(.data[, -cols_with_na])
}


#' @name utils_na_zero
#' @export
select_cols_na <- function(.data, verbose = TRUE){
  cols_with_na <- names(which(sapply(.data, anyNA)))
  if(verbose == TRUE){
    warning("Column(s) with NAs: ", paste(cols_with_na, collapse = ", "), call. = FALSE)
  }
  return(select(.data, cols_with_na))
}

#' @name utils_na_zero
#' @export
select_rows_na <- function(.data, verbose = TRUE){
  rows_with_na <- which(complete.cases(.data) == FALSE)
  if(verbose == TRUE){
    warning("Rows(s) with NAs: ", paste(rows_with_na, collapse = ", "), call. = FALSE)
  }
  return(.data[rows_with_na, ])
}
#' @name utils_na_zero
#' @export
replace_na <- function(.data, ..., replacement = 0){
  col_mean <- !is.na(replacement) && replacement == "colmean"
  col_min <- !is.na(replacement) && replacement == "colmin"
  col_max <- !is.na(replacement) && replacement == "colmax"
  if (has_class(.data, c("data.frame","tbl_df", "data.table"))){
    if(has_rownames(.data)){
      rnames <- rownames(.data)
    }
  if(missing(...)){
    cols_with_na <- names(which(sapply(.data, anyNA)))
    df <- mutate(.data, across(all_of(cols_with_na), ~replace(., is.na(.), replacement)))
    if(col_mean){
    df <- mutate(.data, across(all_of(cols_with_na), ~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)))
    }
    if(col_min){
    df <- mutate(.data, across(all_of(cols_with_na), ~ifelse(is.na(.x), min(.x, na.rm = TRUE), .x)))
    }
    if(col_max){
    df <- mutate(.data, across(all_of(cols_with_na), ~ifelse(is.na(.x), max(.x, na.rm = TRUE), .x)))
    }
  } else{
    df <- mutate(.data, across(c(...), ~replace(., is.na(.), replacement)))
    if(col_mean){
    df <- mutate(.data, across(c(...), ~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)))
    }
    if(col_min){
    df <- mutate(.data, across(c(...), ~ifelse(is.na(.x), min(.x, na.rm = TRUE), .x)))
    }
    if(col_max){
    df <- mutate(.data, across(c(...), ~ifelse(is.na(.x), max(.x, na.rm = TRUE), .x)))
    }
  }
  } else{
    df <- replace(.data, is.na(.data), replacement)
  }
  if(has_rownames(.data)){
    rownames(df) <- rnames
  }
  return(df)
}

#' @name utils_na_zero
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






# Dealing with zeros

#' @name utils_na_zero
#' @export
has_zero <- function(.data){
  if(is.vector(.data)){
    any(.data == 0)
  } else{
    any(apply(.data %>% replace_na(replacement = 1) == 0, 2, any) == TRUE)
  }
}

#' @name utils_na_zero
#' @export
remove_rows_zero <- function(.data, verbose = TRUE){
  row_with_zeros <- which(apply(.data == 0, 1, any) ==  TRUE)
  if(verbose == TRUE){
    warning("Row(s) ", paste(row_with_zeros, collapse = ", "), " with 0s deleted.", call. = FALSE)
  }
  return(remove_rows(.data, all_of(row_with_zeros)))
}

#' @name utils_na_zero
#' @export
remove_cols_zero <- function(.data, verbose = TRUE){
  cols_with_zero <- which(apply(.data == 0, 2, any) ==  TRUE)
  if(verbose == TRUE){
    warning("Column(s) ", paste(names(.data[cols_with_zero]), collapse = ", "), " with 0s deleted.", call. = FALSE)
  }
  return(remove_cols(.data, all_of(cols_with_zero)))
}

#' @name utils_na_zero
#' @export
select_cols_zero <- function(.data, verbose = TRUE){
  cols_with_zero <- which(apply(.data == 0, 2, any) ==  TRUE)
  if(verbose == TRUE){
    warning("Column(s) with 0s: ", paste(names(.data[cols_with_zero]), collapse = ", "), call. = FALSE)
  }
  return(select(.data, cols_with_zero))
}

#' @name utils_na_zero
#' @export
select_rows_zero <- function(.data, verbose = TRUE){
  row_with_zeros <- which(apply(.data == 0, 1, any) ==  TRUE)
  if(verbose == TRUE){
    warning("Rows(s) with 0s: ", paste(rownames(.data[row_with_zeros,]), collapse = ", "), call. = FALSE)
  }
  return(.data[row_with_zeros, ])
}




#' @name utils_na_zero
#' @export
replace_zero <- function(.data, ..., replacement = NA){
  if (has_class(.data, c("data.frame","tbl_df", "data.table"))){
    if(has_rownames(.data)){
      rnames <- rownames(.data)
    }
    test <- !is.na(replacement) && replacement == "colmean"
    if(missing(...)){
      cols_with_zero <- which(apply(.data == 0, 2, any) ==  TRUE)
      df <- mutate(.data, across(all_of(cols_with_zero), ~replace(., . == 0, replacement)))
      if(test){
        df <- mutate(.data, across(all_of(cols_with_zero), ~ifelse(.x == 0, mean(.x, na.rm = TRUE), .x)))
      }
    } else{
      df <- mutate(.data, across(c(...), ~replace(., . == 0, replacement)))
      if(test){
        df <- mutate(.data, across(c(...), ~ifelse(.x == 0, mean(.x, na.rm = TRUE), .x)))
      }
    }
  } else{
    df <- replace(.data, which(.data == 0), replacement)
  }
  if(has_rownames(.data)){
    rownames(df) <- rnames
  }
  return(df)
}
