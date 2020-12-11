#' @title Utilities for handling with NA and zero values
#'
#' @description NAs and zeros can increase the noise in multi-environment trial
#'   analysis. This collection of functions will make it easier to deal with
#'   them.
#' @name utils_na_zero
#' @param .data A data frame.
#' @param ... Variables to fill  \code{NAs} in \code{fill_na()}, replace
#'   \code{NAs} in \code{replace_na()} or zeros in \code{replace_zero()}. If
#'   \code{...} is null then all variables in \code{.data} will be evaluated. It
#'   must be a single variable name or a comma-separated list of unquoted
#'   variables names. Select helpers are also allowed.
#' @param direction Direction in which to fill missing values. Currently either
#'   "down" (the default), "up", "downup" (i.e. first down and then up) or
#'   "updown" (first up and then down).
#' @param prop The proportion (percentage) of \code{NA} values to generate in
#'   \code{.data}.
#' @param replacement The value used for replacement. Defaults to \code{0}. Use
#'   \code{replacement. = "colmean"} to replace missing values with column mean.
#' @param verbose Logical argument. If \code{TRUE} (default) shows in console
#'   the rows or columns deleted.
#' @export
#' @return A data frame with rows or columns with \code{NA} values deleted.
#' @description
#' * \code{fill_na()}: Fills \code{NA} in selected columns using the next or
#' previous entry.
#' * \code{has_na(), has_zero()}: Check for \code{NAs} and \code{0s} in the
#' data and return a logical value.
#' * \code{random_na()}: Generate random \code{NA} values in a two-way table
#' based on a desired proportion.
#' * \code{remove_cols_na(), remove_cols_zero()}: Remove columns with \code{NAs}
#' and \code{0s}, respectively.
#' * \code{remove_rows_na(), remove_rows_zero()}: Remove rows with \code{NAs}
#' and \code{0s}, respectively.
#' * \code{select_cols_na(), select_cols_zero()}: Select columns with \code{NAs}
#' and \code{0s}, respectively.
#' * \code{select_rows_na(), select_rows_zero()}: Select rows with \code{NAs}
#' and \code{0s}, respectively.
#' * \code{replace_na(), replace_zero()}: Replace \code{NAs} and \code{0s},
#' respectively, with a \code{replacement} value.
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
#' @importFrom tidyr fill
has_na <- function(.data){
  any(complete.cases(.data) ==  FALSE)
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
remove_cols_na <- function(.data, verbose = TRUE){
  cols_with_na <- names(which(sapply(.data, anyNA)))
  if(verbose == TRUE){
    warning("Column(s) ", paste(cols_with_na, collapse = ", "), " with NA values deleted.", call. = FALSE)
  }
  return(select(.data, -cols_with_na))
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
  test <- !is.na(replacement) && replacement == "colmean"
  if (has_class(.data, c("data.frame","tbl_df", "data.table"))){
    if(has_rownames(.data)){
      rnames <- rownames(.data)
    }
  if(missing(...)){
    cols_with_na <- names(which(sapply(.data, anyNA)))
    df <- mutate(.data, across(all_of(cols_with_na), ~replace(., is.na(.), replacement)))
    if(test){
    df <- mutate(.data, across(all_of(cols_with_na), ~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)))
    }
  } else{
    df <- mutate(.data, across(c(...), ~replace(., is.na(.), replacement)))
    if(test){
    df <- mutate(.data, across(c(...), ~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)))
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
    print(test)
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
