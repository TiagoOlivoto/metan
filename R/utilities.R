#' @title Utilities for handling with numbers and strings
#' @name utils-num-str
#' @param .data A data frame
#' @param ... Variables to round. If no variable is informed, all the numeric
#'   variables from \code{data} are used.
#' @param digits The number of significant figures.
#' @param var The variable to extract or replace numbers or strings.
#' @param new_var The name of the new variable containing the numbers or
#'   strings extracted or replaced. Defaults to \code{new_var}.
#' @param drop Logical argument. If \code{TRUE} keeps the new variable
#'   \code{new_var} and drops the existing ones. Defaults to \code{FALSE}.
#' @param pattern A string to be matched. Regular Expression Syntax is also
#'   allowed.
#' @param replacement A string for replacement.
#' @param pull Logical argument. If \code{TRUE}, returns the last column (on the
#'   assumption that's the column you've created most recently), as a vector.
#' @description
#' * \code{round_column()}: Round a selected column or a whole data frame to
#' significant figures.
#' * \code{extract_number()}: Extract the number(s) of a string.
#' * \code{replace_number()}: Replace numbers with a replacement.
#' * \code{extract_string()}: Extract all strings, ignoring case.
#' * \code{replace_string()}: Replace all strings with a replacement, ignoring
#' case.
#' * \code{all_lower_case()}: Translate all non-numeric strings of a data frame
#' to lower case. A character vector is also allowed as im
#' * \code{all_upper_case()}: Translate all non-numeric strings of a data frame
#' to upper case.
#' @md
#' @examples
#' \donttest{
#' library(metan)
#'
#' ################ Rounding numbers ###############
#' # All numeric columns
#' round_column(data_ge2, digits = 1)
#'
#' # Round specific columns
#' round_column(data_ge2, EP, digits = 1)
#'
#' ########### Extract or replace numbers ##########
#' # Extract numbers
#' extract_number(data_ge, GEN)
#' extract_number(data_ge,
#'                var = GEN,
#'                drop = TRUE,
#'                new_var = g_number)
#'
#' # Replace numbers
#'
#' replace_number(data_ge, GEN)
#' replace_number(data_ge,
#'                var = GEN,
#'                pattern = "1",
#'                replacement = "_one",
#'                pull = TRUE)
#'
#' ########## Extract or replace strings ##########
#' # Extract strings
#' extract_string(data_ge, GEN)
#' extract_string(data_ge,
#'                var = GEN,
#'                drop = TRUE,
#'                new_var = g_name)
#'
#' # Replace strings
#' replace_string(data_ge, GEN)
#' replace_string(data_ge,
#'                var = GEN,
#'                new_var = GENOTYPE,
#'                pattern = "G",
#'                replacement = "GENOTYPE_")
#' ############# upper and lower cases ############
#' all_lower_case("GENOTYPE")
#' lc <- all_lower_case(data_ge)
#' lc
#' all_lower_case("GENOTYPE")
#'
#' all_upper_case("Genotype")
#' all_upper_case(lc)
#' }
#' @export
#' @importFrom dplyr mutate_if transmute
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
extract_number <- function(.data,
                           var,
                           new_var = new_var,
                           drop = FALSE,
                           pull = FALSE){
  if (drop == FALSE){
    results <-
    mutate(.data,
           {{new_var}} :=   as.numeric(gsub("[^0-9.-]+", "", as.character({{var}})))
    )
    if(pull == TRUE){
      results <- pull(results)
    }
  } else{
    results <-
    transmute(.data,
              {{new_var}} :=   as.numeric(gsub("[^0-9.-]+", "", as.character({{var}})))
    )
    if(pull == TRUE){
      results <- pull(results)
    }
  }
  return(results)
}
#' @name utils-num-str
#' @export
replace_number <- function(.data,
                           var,
                           new_var = new_var,
                           pattern = NULL,
                           replacement = "",
                           drop = FALSE,
                           pull = FALSE){
  if(missing(pattern)){
    pattern <- "[0-9]"
  } else{
    pattern <- pattern
  }
  if (drop == FALSE){
    results <-
    mutate(.data,
           {{new_var}} := gsub(pattern, replacement, as.character({{var}}))
    )
    if(pull == TRUE){
      results <- pull(results)
    }
  } else{
    results <-
    transmute(.data,
              {{new_var}} := gsub(pattern, replacement, as.character({{var}}))
    )
    if(pull == TRUE){
      results <- pull(results)
    }
  }
  return(results)
}
#' @name utils-num-str
#' @export
extract_string <- function(.data,
                           var,
                           new_var = new_var,
                           drop = FALSE,
                           pull = FALSE){
  if (drop == FALSE){
    results <-
    mutate(.data,
           {{new_var}} := as.character(gsub("[^A-z.-]+", "", as.character({{var}})))
    )
    if(pull == TRUE){
      results <- pull(results)
    }
  } else{
    results <-
    transmute(.data,
              {{new_var}} := as.character(gsub("[^A-z.-]+", "", as.character({{var}})))
    )
    if(pull == TRUE){
      results <- pull(results)
    }
  }
  return(results)
}
#' @name utils-num-str
#' @export
replace_string <- function(.data,
                           var,
                           new_var = new_var,
                           pattern = NULL,
                           replacement = "",
                           drop = FALSE,
                           pull = FALSE){
  if(missing(pattern)){
    pattern <- "[A-z]"
  } else{
    pattern <- pattern
  }
  if (drop == FALSE){
    results <-
      mutate(.data,
             {{new_var}} := gsub(pattern, replacement, as.character({{var}}))
      )
    if(pull == TRUE){
      results <- pull(results)
    }
  } else{
    results <-
      transmute(.data,
                {{new_var}} := gsub(pattern, replacement, as.character({{var}}))
      )
    if(pull == TRUE){
      results <- pull(results)
    }
  }
  return(results)
}
#' @name utils-num-str
#' @export
all_upper_case <- function(.data){
  if (any(class(.data) %in% c("data.frame","tbl_df", "data.table"))){
    .data <- as.data.frame(.data)
    mutate_if(.data, ~!is.numeric(.x), toupper) %>%
      as_tibble(rownames = NA)
  } else{
    toupper(.data)
  }
}
#' @name utils-num-str
#' @export
all_lower_case <- function(.data){
  if (any(class(.data) %in% c("data.frame","tbl_df", "data.table"))){
    .data <- as.data.frame(.data)
    mutate_if(.data, ~!is.numeric(.x), tolower) %>%
      as_tibble(rownames = NA)
  } else{
    tolower(.data)
  }
}
NULL

#' @title Utilities for handling with rows and columns
#' @name utils-rows-cols
#' @description
#' * \code{add_cols()}: Add one or more columns to an existing data frame. If
#' specified \code{.before} or \code{.after} columns does not exist, columns are
#' appended at the end of the data. Return a data frame with all the original
#' columns in \code{.data} plus the columns declared in \code{...}. In
#' \code{add_cols()} columns in \code{.data} are available for the expressions.
#' So, it is possible to add a column based on existing data.
#' * \code{add_rows()}: Add one or more rows to an existing data frame. If
#' specified \code{.before} or \code{.after} rows does not exist, rows are
#' appended at the end of the data. Return a data frame with all the original
#' rows in \code{.data} plus the rows declared in \code{...}.
#' * \code{remove_cols()}: Remove one or more columns from a data frame.
#' * \code{remove_rows()}: Remove one or more rows from a data frame.
#' * \code{select_cols()}: Select one or more columns from a data frame.
#' * \code{select_rows()}: Select one or more rows from a data frame.
#' * \code{concatenate()}: Concatenate either two columns of a data frame or a
#' column and specified values. Return a data frame with all the original
#' columns in \code{.data} plus the concatenated variable, .after the last
#' column.
#' * \code{column_exists()}: Checks if a column exists in a data frame. Return a
#' logical value.
#' * \code{get_levels()}: Get the levels of a factor variable.
#' * \code{get_level_size()}: Get the size of each level of a factor variable.
#' * \code{get_all_pairs()}: Get all the possible pairs between the levels of a
#' factor.
#' * \code{select_numeric_cols()}: Select all the numeric columns of a data
#' frame.
#' * \code{select_non_numeric_cols()}: Select all the non-numeric columns of a
#' data frame.
#' @param .data A data frame
#' @param ... Name-value pairs. All values must have one element for each row in
#'   \code{.data} when using \code{add_cols} or one element for each column in
#'   \code{.data} when using \code{add_rows}. Values of length 1 will be
#'   recycled when using \code{add_cols}. For \code{remove_cols()} and
#'   \code{select_cols()},  \code{...} is the column name or column index of the
#'   variable(s) to be dropped. For \code{remove_rows()} and
#'   \code{select_rows()}, \code{...} are the integer row values.
#' @param .before,.after For \code{add_cols}, one-based column index or column
#'   name where to add the new columns, default: .after last column. For
#'   \code{add_rows}, one-based row index where to add the new rows, default:
#'   .after last row.
#' @param col_1,col_2 An unquoted variable name for the first and second
#'   columns to concatenate, respectively. If both columns are in \code{.data},
#'   the values of those columns will be concatenated. Else, \code{col_1} or
#'   \code{col_2} must be a numeric or a character vector of length equals to
#'   the number of rows in \code{.data}. Vectors of length 1 will be recycled.
#' @param new_var The name of the new variable containing the concatenated
#'   values from \code{col_1} and \code{col_2}. Defaults to \code{new_var}.
#' @param sep The separator to be between values of \code{col_1} and
#'   \code{col_2}. Defaults to "_".
#' @param drop Logical argument. If \code{TRUE} keeps the new variable
#'   \code{new_var} and drops the existing ones. Defaults to \code{FALSE}.
#' @param pull Logical argument. If \code{TRUE}, returns the last column (on the
#'   assumption that's the column you've created most recently), as a vector.
#' @param cols A quoted variable name to check if it exists in \code{.data}
#' @param group A factor variable to get the levels.
#' @param levels The levels of a factor or a numeric vector.
#' @md
#' @importFrom  tibble add_column add_row
#' @export
#' @examples
#' \donttest{
#' library(metan)
#'
#' ################# Adding columns #################
#' # Variables x and y .after last column
#' data_ge %>%
#'   add_cols(x = 10,
#'            y = 30)
#' # Variables x and y .before the variable GEN
#' data_ge %>%
#'   add_cols(x = 10,
#'            y = 30,
#'            .before = "GEN")
#'
#' # Creating a new variable based on the existing ones.
#' data_ge %>%
#'   add_cols(GY2 = GY^2,
#'            GY2_HM = GY2 + HM,
#'            .after = "GY")
#'
#' ####### Selecting and removing columns ##########
#' select_cols(data_ge2, GEN, REP)
#' select_cols(data_ge2, 2:3)
#' remove_cols(data_ge2, GEN, REP)
#' remove_cols(data_ge2, 2:3)
#'
#' ######## Selecting and removing rows ###########
#' select_rows(data_ge2, GEN, REP)
#' select_rows(data_ge2, 2:3)
#' remove_rows(data_ge2, GEN, REP)
#' remove_rows(data_ge2, 2:3)
#'
#' ########### Concatenating columns ################
#' # Both variables in data
#' concatenate(data_ge, ENV, GEN)
#'
#' # Combine with add_cols()
#' data_ge2 %>%
#'   add_cols(ENV_GEN = concatenate(., ENV, GEN, pull = TRUE),
#'            .after = "GEN")
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
#'             new_var = ENVIRONMENT,
#'             sep = "-")
#'
#' ################### Adding rows ##################
#' data_ge %>%
#'   add_rows(ENV = "E_TEST",
#'            GEN = "G_TEST",
#'            REP = 3,
#'            GY = 10.3,
#'            HM = 100.11,
#'            .after = 1)
#' ########## checking if a column exists ###########
#' column_exists(data_g, "GEN")
#'
#' ####### get the levels and size of levels ########
#' get_levels(data_g, GEN)
#' get_level_size(data_g, GEN)
#'
#' ############## all possible pairs ################
#' all_pairs(data_g, GEN)
#'
#' ########## select numeric variables only #########
#' select_numeric_cols(data_g)
#' select_non_numeric_cols(data_g)
#' }
add_cols <- function(.data, ..., .before = NULL, .after = NULL){
  if (!missing(.after) && .after == names(.data[ncol(.data)])){
    message("Putting variables .after the last column. Setting '.after' to NULL.", call. = FALSE)
    .after = NULL
  }
   if (!missing(.before)){
    if(is.character(.before)){
      if(!(.before %in% colnames(.data))){
        stop("Column '.before' not in .data")
      }
    } else{
      .before <- colnames(data_ge[.before])
    }
    bfr <- .data[,1:which(colnames(.data) ==  .before)-1]
    aft <- select(.data, -!!colnames(bfr))
    df2 <-
      .data %>%
      mutate(...) %>%
      select(-!!colnames(bfr), -!!colnames(aft))
    results <- cbind(bfr, df2, aft) %>%
      as_tibble()
  } else if (!is.null(.after)) {
    if(is.character(.after)){
      if(!(.after %in% colnames(.data))){
        stop("Column '.after' not in .data")
      }
    } else{
      .after <- colnames(data_ge[.after])
    }
    aft <- .data[,(which(colnames(.data) ==  .after)+1):ncol(.data)]
    bfr <- select(.data, -!!colnames(aft))
    df2 <-
      .data %>%
      mutate(...) %>%
      select(-!!colnames(bfr), -!!colnames(aft))
    results <- cbind(bfr, df2, aft) %>%
      as_tibble()
  } else{
    results <- mutate(.data, ...)
  }
  return(results)
  }
#' @name utils-rows-cols
#' @export
add_rows <- function(.data, ..., .before = NULL, .after = NULL){
  if(is.character(.before)){
    if(!(.before %in% rownames(.data))){
      .before <- NULL
    }
  }
  if(is.character(.after)){
    if(!(.after %in% rownames(.data))){
      .after <- NULL
    }
  }
  add_row(.data, ..., .before = .before, .after = .after)
}
#' @name utils-rows-cols
#' @importFrom rlang quos
#' @export
remove_cols <- function(.data, ...){
  select(.data, -c(!!!quos(...)))
}
#' @name utils-rows-cols
#' @export
remove_rows <- function(.data, ...){
  slice(.data, -c(!!!quos(...)))
}
#' @name utils-rows-cols
#' @export
select_cols <- function(.data, ...){
  select(.data, ...)
}
#' @name utils-rows-cols
#' @export
select_rows <- function(.data, ...){
  slice(.data, ...)
}
#' @name utils-rows-cols
#' @export
concatenate <- function(.data,
                        col_1,
                        col_2,
                        new_var = new_var,
                        sep = "_",
                        drop = FALSE,
                        pull = FALSE){
  if (drop == FALSE){
    results <-
      .data %>%
      mutate({{new_var}} := paste({{col_1}}, sep, {{col_2}}, sep = ""))
    if (pull == TRUE){
      results <- pull(results)
    }
  } else{
    results <-
      .data %>%
      transmute({{new_var}} := paste({{col_1}}, sep, {{col_2}}, sep = ""))
    if (pull == TRUE){
      results <- pull(results)
    }
  }
  return(results)
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
#' @name utils-rows-cols
#' @export
get_levels <- function(.data, group){
  .data %>%
    mutate_if(~!is.numeric(.x), as.factor) %>%
    pull({{group}}) %>%
    levels()
}
#' @name utils-rows-cols
#' @importFrom dplyr count
#' @export
get_level_size <- function(.data, group){
  result <- .data %>%
    group_by({{group}}) %>%
    count()
  n <- result$n
  names(n) <- result %>% pull(1)
  return(n)
}
#' @name utils-rows-cols
#' @export
all_pairs <- function(.data, levels){
  levels <-
    get_levels(.data, {{levels}})
  combn(levels, 2) %>%
    as.data.frame() %>%
    t()
}
#' @name utils-rows-cols
#' @export
select_numeric_cols <- function(.data){
  if(is_grouped_df(.data)){
    .data <- ungroup(.data)
  }
  select_if(.data, is.numeric)
}
#' @name utils-rows-cols
#' @export
select_non_numeric_cols <- function(.data){
  if(is_grouped_df(.data)){
    .data <- ungroup(.data)
  }
  select_if(.data, ~!is.numeric(.x))
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
#' \donttest{
#' library(metan)
#' means_by(data_ge2, ENV)
#' means_by(data_ge2, GEN, ENV)
#'}
means_by <- function(.data, ...){
  .data %>%
    group_by(...) %>%
    summarise_if(is.numeric, mean)
}


