#' @importFrom dplyr mutate
#' @export
dplyr::mutate
#' @importFrom dplyr mutate_if
#' @export
dplyr::mutate_if
#' @importFrom dplyr filter
#' @export
dplyr::filter


#' @title Utilities for handling with numbers and strings
#' @name utils-num-str
#' @param .data A data frame
#' @param ... Variables to round. If no variable is informed, all the numeric
#'   variables from \code{data} are used.
#' @param digits The number of significant figures.
#' @param x A numeric value or a string. Vectors are also allowed.
#' @param pattern A string to be matched. Regular Expression Syntax is also allowed.
#' @param replacement A string for replacement.
#' @description
#' * \code{round_column()}: Round a selected column or a whole data frame to significant figures.
#' * \code{extract_number()}: Extract the number(s) of a string.
#' * \code{replace_number()}: Replace numbers with a replacement.
#' * \code{extract_string()}: Extract all strings, ignoring case.
#' * \code{replace_string()}: Replace all strings with a replacement, ignoring case.
#' @md
#' @examples
#' library(metan)
#' ################ Rounding numbers ###############
#' # All numeric columns
#' round_column(data_ge2, digits = 1)
#'
#' # Round specific columns
#' round_column(data_ge2, EP, digits = 1)
#'
#' ######## Extract or replace numbers ########
#' # Extract numbers
#' extract_number("GEN_10")
#' extract_number(data_ge2$GEN)
#'
#' mutate(data_ge,
#'        GEN_NUMERIC = extract_number(GEN))
#'
#' # Replace numbers
#' mutate(data_ge,
#'        GEN_NUMERIC = replace_number(GEN))
#'
#' mutate(data_ge,
#'        GEN_NUMERIC =
#'          replace_number(GEN,
#'                         pattern = "2",
#'                         with = "_two"))
#'
#' ######## Extract or replace strings ########
#'
#' # Extract strings
#' extract_string("GEN_10")
#' extract_string(data_ge2$GEN)
#'
#' mutate(data_ge,
#'        GEN_STRING = extract_string(GEN))
#'
#' # Replace strings
#' mutate(data_ge,
#'        GEN_STRING = replace_string(GEN))
#'
#' mutate(data_ge,
#'        GEN_CODE =
#'          replace_string(GEN,
#'                         pattern = "G",
#'                         with = "GENOTYPE_"))
#' @export
round_column <- function(.data, ...,  digits = 2){
  if (missing(...)){
    .data %<>% dplyr::mutate_if(is.numeric, round, digits = digits)
  } else{
    .data %<>% dplyr::mutate_at(vars(...), round, digits = digits)
  }
  return(.data)
}

#' @name utils-num-str
#' @export
extract_number <- function(x){
  as.numeric(gsub("[^0-9.-]+", "", as.character(x)))
}


#' @name utils-num-str
#' @export
replace_number <- function(x, pattern = NULL, replacement = ""){
  if(missing(pattern)){
    pattern <- "[0-9]"
  } else{
    pattern <- pattern
  }
  gsub(pattern, replacement, as.character(x))
}


#' @name utils-num-str
#' @export
extract_string <- function(x){
  as.character(gsub("[^A-z.-]+", "", as.character(x)))
}


#' @name utils-num-str
#' @export
replace_string <- function(x, pattern = NULL, replacement = ""){
  if(missing(pattern)){
    pattern <- "[A-z]"
  } else{
    pattern <- pattern
  }
  gsub(pattern, replacement, as.character(x))
}




#' @title Utilities for handling with columns and rows
#' @name utils-rows-cols
#' @description
#' * \code{add_cols()}: Add one or more columns to an existing data frame. If
#' specified .before or .after columns does not exist, columns are appended at
#' the end of the data. Return a data frame with all the original columns in
#' \code{.data} plus the columns declared in \code{...}.
#'
#' * \code{concatenate()}: Concatenate either two columns of a data frame or a
#' column and specified values. Return a data frame with all the original
#' columns in \code{.data} plus the concatenated variable, after the last
#' column.
#'
#' * \code{column_exists()}: Checks if a column exists in a data frame. Return a logical value.
#'
#' @param .data A data frame
#' @param ... Name-value pairs. All values must have one element for each row in
#'   \code{.data} when using \code{add_cols} or one element for each column in
#'   \code{.data} when using \code{add_rows}. Values of length 1 will be
#'   recycled when using \code{add_cols}.
#' @param before,after For \code{add_cols}, one-based column index or column
#'   name where to add the new columns, default: after last column. For
#'   \code{add_rows}, one-based row index where to add the new rows, default:
#'   after last row.
#' @param col_1,col_2 An unquoted variable name for the first and second
#'   columns to concatenate, respectively. If both columns are in \code{.data},
#'   the values of those columns will be concatenated. Else, \code{col_1} or
#'   \code{col_2} must be a numeric or a character vector of length equals to
#'   the number of rows in \code{.data}. Vectors of length 1 will be recycled.
#' @param var_name The name of the new variable containing the concatenated
#'   values from \code{col_1} and \code{col_2}. Defaults to \code{new_var}.
#' @param sep The separator to be between values of \code{col_1} and
#'   \code{col_2}. Defaults to "_".
#' @param cols A quoted variable name to check if it exists in \code{.data}
#' @md
#' @importFrom  tibble add_column add_row
#' @export
#' @examples
#' library(metan)
#'
#' ################# Adding columns #################
#' # Variables x and y after last column
#' data_ge2 %>%
#'   add_columns(x = 10,
#'               y = 30)
#' # Variables x and y before the variable GEN
#' data_ge2 %>%
#'   add_columns(x = 10,
#'               y = 30,
#'               before = "GEN")
#'
#'
#' ########### Concatenating columns ################
#' # Both variables in data
#' concatenate(data_ge, ENV, GEN)
#'
#' # One variable in data
#' concatenate(data_ge, ENV, rep(1:14, each=30))
#'
#' # One variable in data, recycling the second one.
#' concatenate(data_ge, ENV, 1)
#'
#' # None variables in data
#' # Set name of the concatenate variable
#' # Set the separator symbol
#' concatenate(data_ge,
#'             col_1 = rep("ENV", 420),
#'             col_2 = rep(1:14, each=30),
#'             var_name = ENVIRONMENT,
#'             sep = "-")
#' ################### Adding rows ##################
#' data_ge %>%
#'   add_rows(ENV = "TEST",
#'            GEN = "G_TEST",
#'            REP = 3,
#'            GY = 10.3,
#'            HM = 100.11,
#'            after = 1)
#' ########## checking if a column exists ###########
#'
add_columns <- function(.data, ..., before = NULL, after = NULL){
  if(is.character(before)){
    if(!(before %in% colnames(.data))){
      before <- NULL
    }
  }
  if(is.character(after)){
    if(!(.after %in% colnames(.data))){
      after <- NULL
    }
  }
  add_column(.data, ..., .before = before, .after = after)
}

#' @name utils-rows-cols
#' @export
add_rows <- function(.data, ..., before = NULL, after = NULL){
  if(is.character(before)){
    if(!(before %in% rownames(.data))){
      before <- NULL
    }
  }
  if(is.character(after)){
    if(!(.after %in% rownames(.data))){
      after <- NULL
    }
  }
  add_row(.data, ..., .before = before, .after = after)
}

#' @name utils-rows-cols
#' @export
concatenate <- function(.data, col_1, col_2, var_name = new_var,  sep = "_"){
  .data %>%
    mutate(
      {{var_name}} := paste({{col_1}}, sep, {{col_2}}, sep = "")
    )
}


#' @name utils-rows-cols
#' @export
column_exists <-function(.data, cols){
  if(length(setdiff(cols, colnames(.data))) != 0){
    FALSE
  } else{
    TRUE
  }
}


#' @title Means by one or more factors
#' @description Computes the mean for all numeric variables of a data frame,
#'   grouping by one or more factors.
#' @param .data A data frame
#' @param ... One or more categorical variables for grouping the data.
#' @return An object of class tbl_df with the computed means by each level of
#'   the factor(s) declared in \code{...}.
#' @export
#' @examples
#' library(metan)
#' means_by(data_ge2, ENV)
#' means_by(data_ge2, GEN, ENV)
means_by <- function(.data, ...){
  .data %>%
    group_by(...) %>%
    summarise_if(is.numeric, mean)
}


