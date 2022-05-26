#' @title Utilities for handling with numbers and strings
#' @name utils_num_str
#' @param .data A data frame
#' @param ... The argument depends on the function used.
#'  * For `round_cols()` `...` are the variables to round. If no
#'  variable is informed, all the numeric variables from `data` are used.
#'  * For `all_lower_case()`, `all_upper_case()`,
#'   `all_title_case()`, `stract_number()`, `stract_string()`,
#'   `remove_strings()`, and `tidy_strings()` `...` are the
#'   variables to apply the function. If no variable is informed, the function
#'   will be applied to all non-numeric variables in `.data`.
#' @param digits The number of significant figures.
#' @param pattern A string to be matched. Regular Expression Syntax is also
#'   allowed.
#' @param replacement A string for replacement.
#' @param ignore_case If `FALSE` (default), the pattern matching is case
#'   sensitive and if `TRUE`, case is ignored during matching.
#' @param sep A character string to separate the terms. Defaults to "_".
#' @description
#' `r badge('stable')`
#'
#' * `all_lower_case()`: Translate all non-numeric strings of a data frame
#' to lower case.
#' * `all_upper_case()`: Translate all non-numeric strings of a data frame
#' to upper case.
#' * `all_title_case()`: Translate all non-numeric strings of a data frame
#' to title case.
#' * `first_upper_case`: Translate the first word of a string to upper
#' case.
#' * `extract_number()`: Extract the number(s) of a string.
#' * `extract_string()`: Extract all strings, ignoring case.
#' * `find_text_in_num()`: Find text characters in a numeric sequence and
#' return the row index.
#' * `has_text_in_num()`: Inspect columns looking for text in numeric
#' sequence and return a warning if text is found.
#' * `remove_space()`: Remove all blank spaces of a string.
#' * `remove_strings()`: Remove all strings of a variable.
#' * `replace_number()`: Replace numbers with a replacement.
#' * `replace_string()`: Replace all strings with a replacement, ignoring
#' case.
#' * `round_cols()`: Round a selected column or a whole data frame to
#' significant figures.
#' * `tidy_strings()`: Tidy up characters strings, non-numeric columns, or
#' any selected columns in a data frame by putting all word in upper case,
#' replacing any space, tabulation, punctuation characters by `'_'`, and
#' putting `'_'` between lower and upper case. Suppose that `str =
#' c("Env1", "env 1", "env.1")` (which by definition should represent a unique
#' level in plant breeding trials, e.g., environment 1) is subjected to
#' `tidy_strings(str)`: the result will be then `c("ENV_1", "ENV_1",
#' "ENV_1")`. See Examples section for more examples.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(metan)
#'
#' ################ Rounding numbers ###############
#' # All numeric columns
#' round_cols(data_ge2, digits = 1)
#'
#' # Round specific columns
#' round_cols(data_ge2, EP, digits = 1)
#'
#' ########### Extract or replace numbers ##########
#' # Extract numbers
#' extract_number(data_ge, GEN)

#' # Replace numbers
#' replace_number(data_ge, GEN)
#' replace_number(data_ge,
#'                GEN,
#'                pattern = 1,
#'                replacement = "_one")
#'
#' ########## Extract, replace or remove strings ##########
#' # Extract strings
#' extract_string(data_ge, GEN)

#'
#' # Replace strings
#' replace_string(data_ge, GEN)
#' replace_string(data_ge,
#'                GEN,
#'                pattern = "G",
#'                replacement = "GENOTYPE_")
#'
#' # Remove strings
#' remove_strings(data_ge)
#' remove_strings(data_ge, ENV)
#'
#'
#' ############ Find text in numeric sequences ###########
#' mixed_text <- data.frame(data_ge)
#' mixed_text[2, 4] <- "2..503"
#' mixed_text[3, 4] <- "3.2o75"
#' find_text_in_num(mixed_text, GY)
#'
#' ############# upper, lower and title cases ############
#'gen_text <- c("This is the first string.", "this is the second one")
#'all_lower_case(gen_text)
#'all_upper_case(gen_text)
#'all_title_case(gen_text)
#'first_upper_case(gen_text)
#'
#'# A whole data frame
#'all_lower_case(data_ge)
#'
#'
#' ############### Tidy up messy text string ##############
#' messy_env <- c("ENV 1", "Env   1", "Env1", "env1", "Env.1", "Env_1")
#' tidy_strings(messy_env)
#'
#' messy_gen <- c("GEN1", "gen 2", "Gen.3", "gen-4", "Gen_5", "GEN_6")
#' tidy_strings(messy_gen)
#'
#' messy_int <- c("EnvGen", "Env_Gen", "env gen", "Env Gen", "ENV.GEN", "ENV_GEN")
#' tidy_strings(messy_int)
#'
#' library(tibble)
#' # Or a whole data frame
#' df <- tibble(Env = messy_env,
#'              gen = messy_gen,
#'              Env_GEN = interaction(Env, gen),
#'              y = rnorm(6, 300, 10))
#' df
#' tidy_strings(df)
#' }
#' @export
all_upper_case <- function(.data, ...){
  helper_case(.data, toupper, ...)
}
#' @name utils_num_str
#' @export
all_lower_case <- function(.data, ...){
  helper_case(.data, tolower, ...)
}
#' @name utils_num_str
#' @export
all_title_case <- function(.data, ...){
  helper_case(.data, to_title, ...)
}
#' @name utils_num_str
#' @export
first_upper_case <- function(.data, ...){
  helper_case(.data, first_upper, ...)
}
#' @name utils_num_str
#' @export
#'
extract_number <- function(.data,
                           ...,
                           pattern = NULL){
  if(missing(pattern)){
    pattern <- "[^0-9.-]+"
  }
  if (has_class(.data, c("data.frame","tbl_df", "data.table"))){
    if(missing(...)){
      results <-
        mutate(.data, across(where(~!is.numeric(.)), gsub, pattern = pattern, replacement = "")) %>%
        as_numeric(everything())
    } else{
      results <-
        mutate(.data, across(c(...), gsub, pattern = pattern, replacement = "")) %>%
        as_numeric(...)
    }
    return(results)
  } else{
    return(as.numeric(gsub("[^0-9.-]+", "", .data)))
  }
}
#' @name utils_num_str
#' @export
extract_string <- function(.data,
                           ...,
                           pattern = NULL){
  if(missing(pattern)){
    pattern <- "[^A-z.-]+"
  }
  if (has_class(.data, c("data.frame","tbl_df", "data.table"))){
    if(missing(...)){
      results <-
        mutate(.data, across(where(~!is.numeric(.)), gsub, pattern = pattern, replacement = ""))
    } else{
      results <-
        mutate(.data, across(c(...), gsub, pattern = pattern, replacement = ""))
    }
    return(results)
  } else{
    return(as.character(gsub("[^A-z.-]+", "", .data)))
  }
}
#' @name utils_num_str
#' @export
find_text_in_num <- function(.data, ...){
  help_text <- function(x){
    which(is.na(suppressWarnings(as.numeric(x))))
  }
  if(!missing(...)){
    .data <- select_cols(.data, ...)
    apply(.data, 2, help_text)

  } else{
    paste0("Line(s): ", paste0(help_text(.data), collapse = ","))
  }
}

#' @name utils_num_str
#' @export
has_text_in_num <- function(.data){
  if(any(sapply(sapply(.data, find_text_in_num), length)> 0)){
    var_nam <- names(which(sapply(sapply(.data, find_text_in_num), length)>0))
    warning("Non-numeric variable(s) ", paste(var_nam, collapse = ", "),
            " will not be selected.\nUse 'find_text_in_num()' to see rows with non-numeric characteres.", call. = FALSE)
  }
}
#' @name utils_num_str
#' @export
remove_space <- function(.data, ...){
  if (has_class(.data, c("data.frame","tbl_df", "data.table"))){
    if(missing(...)){
      results <-
        mutate(.data,
               across(where(~!is.numeric(.x)), gsub, pattern = "[[:space:]]", replacement = ""))
    } else{
      results <-
        mutate(.data,
               across(c(...), gsub, pattern = "[[:space:]]", replacement = ""))
    }

    return(results)
  } else{
    return(gsub(pattern = "[[:space:]]", replacement = "", .data))
  }
}

#' @name utils_num_str
#' @export
remove_strings <- function(.data, ...){
  if (has_class(.data, c("data.frame","tbl_df", "data.table"))){
    if(missing(...)){
      results <-
        mutate(.data, across(everything(), gsub, pattern = "[^0-9.-]", replacement = "")) %>%
        as_numeric(everything())
    } else{
      results <-
        mutate(.data, across(c(...), gsub, pattern = "[^0-9.-]", replacement = "")) %>%
        as_numeric(...)
    }
    return(results)
  } else{
    return(gsub(pattern = "[^0-9.-]", replacement = "", .data)) %>%
      remove_space() %>%
      as.numeric()
  }
}
#' @name utils_num_str
#' @export
replace_number <- function(.data,
                           ...,
                           pattern = NULL,
                           replacement = "",
                           ignore_case = FALSE){
  if(missing(pattern)){
    pattern <- "[0-9]"
  } else{
    pattern <- pattern
  }
  if (has_class(.data, c("data.frame","tbl_df", "data.table"))){
    results <-
      .data %>%
      mutate(across(c(...), ~gsub(pattern, replacement, ., ignore.case = ignore_case)))
    return(results)
  } else{
    return(gsub(pattern, replacement, .data, ignore.case = ignore_case))
  }
}
#' @name utils_num_str
#' @export
replace_string <- function(.data,
                           ...,
                           pattern = NULL,
                           replacement = "",
                           ignore_case = FALSE){
  if(missing(pattern)){
    pattern <- "[A-z]"
  } else {
    pattern <- pattern
  }
  if (has_class(.data, c("data.frame","tbl_df", "data.table"))){
    results <-
      .data %>%
      mutate(across(c(...), ~gsub(pattern, replacement, ., ignore.case = ignore_case)))
    return(results)
  } else{
    return(gsub(pattern, replacement, .data, ignore.case = ignore_case))
  }
}

#' @name utils_num_str
#' @export
#' @importFrom dplyr transmute
round_cols <- function(.data, ...,  digits = 2){
  is_mat <- is.matrix(.data)
  if(is_mat ==  TRUE){
    .data <-
      .data %>%
      as.data.frame() %>%
      rownames_to_column()
  }
  rn_test <- has_rownames(.data)
  if(rn_test == TRUE){
    rnames <- rownames(.data)
  }
  if (missing(...)){
    .data %<>% mutate(across(where(is.numeric), round, digits = digits))
  } else{
    .data %<>% mutate(across(c(...), round, digits = digits))
  }
  if(rn_test == TRUE){
    rownames(.data) <- rnames
  }
  if(is_mat ==  TRUE){
    .data <-
      .data %>%
      column_to_rownames() %>%
      as.matrix()
  }
  return(.data)
}
#' @name utils_num_str
#' @export
tidy_strings <- function(.data, ..., sep = "_"){
  fstr <- function(x){
    str <- toupper(gsub("(?<=[a-z0-9])(?=[A-Z])|[[:space:][:punct:]]+", sep, x, perl = TRUE))
    str <- gsub("(?<=[A-Z])(?=[0-9])|(?<=[0-9])(?=[A-Z])", sep, str, perl = TRUE)
    return(str)
  }
  if (has_class(.data, c("data.frame","tbl_df", "data.table"))){
    if(missing(...)){
      results <-
        mutate(.data, across(where(~!is.numeric(.x)),  fstr))
    } else{
      results <-
        mutate(.data, across(c(...), fstr))
    }
    return(results)
  } else{
    return(fstr(.data))
  }
}



#' @title Utilities for handling with rows and columns
#' @name utils_rows_cols
#' @description
#' `r badge('stable')`
#'
#' * `add_cols()`: Add one or more columns to an existing data frame. If
#' specified `.before` or `.after` columns does not exist, columns are
#' appended at the end of the data. Return a data frame with all the original
#' columns in `.data` plus the columns declared in `...`. In
#' `add_cols()` columns in `.data` are available for the expressions.
#' So, it is possible to add a column based on existing data.
#' * `add_rows()`: Add one or more rows to an existing data frame. If
#' specified `.before` or `.after` rows does not exist, rows are
#' appended at the end of the data. Return a data frame with all the original
#' rows in `.data` plus the rows declared in `...` argument.
#' * `add_row_id()`: Add a column with the row id as the first column in
#' `.data`.
#' * `add_prefix()` and `add_suffix()` add prefixes and suffixes,
#' respectively, in variable names selected in `...` argument.
#' * `all_pairs()`: Get all the possible pairs between the levels of a
#' factor.
#' * `colnames_to_lower()`: Translate all column names to lower case.
#' * `colnames_to_upper()`: Translate all column names to upper case.
#' * `colnames_to_title()`: Translate all column names to title case.
#' * `column_exists()`: Checks if a column exists in a data frame. Return a
#' logical value.
#' * `columns_to_first()`: Move columns to first positions in `.data`.
#' * `columns_to_last()`: Move columns to last positions in `.data`.
#' * `columns_to_rownames()`: Move a column of `.data` to its row
#' names.
#' * `rownames_to_column()`: Move the row names of `.data` to a new
#' column.
#' * `remove_rownames()`: Remove the row names of `.data`.
#' * `concatenate()`: Concatenate columns of a data frame. If `drop =
#' TRUE` then the existing variables are dropped. If `pull = TRUE` then the
#' concatenated variable is pull out to a vector. This is specially useful when
#' using `concatenate` to add columns to a data frame with `add_cols()`.
#' * `get_levels()`: Get the levels of a factor variable.
#' * `get_levels_comb()`: Get the combination of the levels of a factor.
#' * `get_level_size()`: Get the size of each level of a factor variable.
#' * `remove_cols()`: Remove one or more columns from a data frame.
#' * `remove_rows()`: Remove one or more rows from a data frame.
#' * `reorder_cols()`: Reorder columns in a data frame.
#' * `select_cols()`: Select one or more columns from a data frame.
#' * `select_first_col()`: Select first variable, possibly with an offset.
#' * `select_last_col()`: Select last variable, possibly with an offset.
#' * `select_numeric_cols()`: Select all the numeric columns of a data
#' frame.
#' * `select_non_numeric_cols()`: Select all the non-numeric columns of a
#' data frame.
#' * `select_rows()`: Select one or more rows from a data frame.
#' * `tidy_colnames()`: Tidy up column names with
#' [tidy_strings()].
#'
#'
#' @param .data A data frame
#' @param ... The argument depends on the function used.
#' * For `add_cols()` and `add_rows()` is name-value pairs. All values
#' must have one element for each row in `.data` when using
#' `add_cols()` or one element for each column in `.data` when using
#' `add_rows()`. Values of length 1 will be recycled when using
#' `add_cols()`.
#'
#' * For `remove_cols()` and `select_cols()`,  `...` is the
#' column name or column index of the variable(s) to be dropped.
#'
#' * For `add_prefix()` and `add_suffix()`,  `...` is the column
#' name to add the prefix or suffix, respectively. Select helpers are allowed.
#'
#' * For `columns_to_first()` and `columns_to_last()`,  `...` is
#' the column name or column index of the variable(s) to be moved to first or
#' last in `.data`.
#'
#' * For `remove_rows()` and `select_rows()`, `...` is an integer
#' row value.
#'
#' * For `concatenate()`, `...` is the unquoted variable names to be
#' concatenated.
#'
#' * For `get_levels()`, `get_level_comb()`, and
#' `get_level_size()` `...` is the unquoted variable names to get the
#' levels, levels combinations and levels size, respectively.
#'
#' @param .before,.after For `add_cols()`, `concatenate()`, and
#'   `reorder_cols()`, one-based column index or column name where to add
#'   the new columns, default: .after last column. For `add_rows()`,
#'   one-based row index where to add the new rows, default: .after last row.
#' @param prefix,suffix The prefix and suffix used in `add_prefix()` and
#'   `add_suffix()`, respectively.
#' @param new_var The name of the new variable containing the concatenated
#'   values. Defaults to `new_var`.
#' @param var Name of column to use for rownames.
#' @param sep The separator to appear when using `concatenate()`,
#'   `add_prefix()`, or `add_suffix()`. Defaults to to `"_"`.
#' @param drop Logical argument. If `TRUE` keeps the new variable
#'   `new_var` and drops the existing ones. Defaults to `FALSE`.
#' @param pull Logical argument. If `TRUE`, returns the last column (on the
#'   assumption that's the column you've created most recently), as a vector.
#' @param cols A quoted variable name to check if it exists in `.data`.
#' @param levels The levels of a factor or a numeric vector.
#' @param offset Set it to *n* to select the *n*th variable from the
#'   end (for `select_last_col()`) of from the begin (for
#'   `select_first_col()`)
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @importFrom  tibble add_row
#' @importFrom dplyr relocate across
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
#'            .before = GEN)
#'
#' # Creating a new variable based on the existing ones.
#' data_ge %>%
#'   add_cols(GY2 = GY^2,
#'            GY2_HM = GY2 + HM,
#'            .after = GY)
#'
#' ############### Reordering columns ###############
#' reorder_cols(data_ge2, NKR, .before = ENV)
#' reorder_cols(data_ge2, where(is.factor), .after = last_col())
#'
#' ######## Selecting and removing columns ##########
#' select_cols(data_ge2, GEN, REP)
#' remove_cols(data_ge2, GEN, REP)
#'
#' ########## Selecting and removing rows ###########
#' select_rows(data_ge2, 2:3)
#' remove_rows(data_ge2, 2:3)
#'
#' ########### Concatenating columns ################
#' concatenate(data_ge, ENV, GEN, REP)
#' concatenate(data_ge, ENV, GEN, REP, drop = TRUE)
#'
#' # Combine with add_cols() and replace_string()
#'data_ge2 %>%
#'  add_cols(ENV_GEN = concatenate(., ENV, GEN, pull = TRUE),
#'           .after = GEN) %>%
#'  replace_string(ENV_GEN,
#'                 pattern = "H",
#'                 replacement = "HYB_")
#'
#' # Use prefixes and suffixes
#' concatenate(data_ge2, REP, prefix = "REP", new_var = REP)
#'
#' # Use prefixes and suffixes (the ear traits EH, EP, EL, and ED)
#' add_prefix(data_ge2, PH, EH, EP, EL, prefix = "EAR")
#' add_suffix(data_ge2, PH, EH, EP, EL, suffix = "EAR", sep = ".")
#'
#' # Use prefixes and suffixes (colnames)
#' concatenate(data_ge2, REP, prefix = "REP", new_var = REP)
#'
#'
#' ########### formating column names ###############
#' # Creating data with messy column names
#' df <- head(data_ge, 3)
#' colnames(df) <- c("Env", "gen", "Rep", "GY", "hm")
#' df
#' colnames_to_lower(df)
#' colnames_to_upper(df)
#' colnames_to_title(df)
#'
#'
#' ################### Adding rows ##################
#' data_ge %>%
#'   add_rows(GY = 10.3,
#'            HM = 100.11,
#'            .after = 1)
#'
#' ########## checking if a column exists ###########
#' column_exists(data_g, "GEN")
#'
#' ####### get the levels, level combinations and size of levels ########
#' get_levels(data_g, GEN)
#' get_levels_comb(data_ge, ENV, GEN)
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
  results <- mutate(.data, ...)
  pos_dots <- (ncol(results)) - ((ncol(results) - ncol(.data)) - 1)
  if (!missing(.before) | !missing(.after)){
    results <- reorder_cols(results,
                            all_of(pos_dots):ncol(results),
                            .before = {{.before}},
                            .after = {{.after}})
  }
  return(as_tibble(results, rownames = NA))
}
#' @name utils_rows_cols
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
#' @name utils_rows_cols
#' @export
add_row_id <- function(.data, var = "row_id"){
  if(!has_class(.data, "data.frame")){
    stop("The object '",  match.call()[[".data"]], "' must be a data frame.")
  }
  df <- .data
  add_cols(df, {{var}} := seq_len(nrow(df)), .before = 1)
}
#' @name utils_rows_cols
#' @export
all_pairs <- function(.data, levels){
  levels <-
    get_levels(.data, {{levels}})
  combn(levels, 2) %>%
    t() %>%
    as.data.frame()
}
#' @name utils_rows_cols
#' @export
#' @importFrom dplyr rename_with
add_prefix <- function(.data, ..., prefix, sep = "_"){
  rename_with(.data, ~paste(prefix, ., sep = sep), c(...))
}
#' @name utils_rows_cols
#' @export
add_suffix <- function(.data, ..., suffix, sep = "_"){
  rename_with(.data, ~paste(., suffix, sep = sep), c(...))
}
#' @name utils_rows_cols
#' @export
colnames_to_lower <- function(.data){
  colnames(.data) <- all_lower_case(colnames(.data))
  return(.data)
}
#' @name utils_rows_cols
#' @export
colnames_to_upper <- function(.data){
  colnames(.data) <- all_upper_case(colnames(.data))
  return(.data)
}
#' @name utils_rows_cols
#' @export
colnames_to_title <- function(.data){
  colnames(.data) <- all_title_case(colnames(.data))
  return(.data)
}
#' @name utils_rows_cols
#' @export
column_to_first <-function(.data, ...){
  select_cols(.data, ..., everything())
}
#' @name utils_rows_cols
#' @export
column_to_last <-function(.data, ...){
  select_cols(.data, -c(!!!quos(...)), everything())
}
#' @name utils_rows_cols
#' @export
column_to_rownames <- function(.data, var = "rowname"){
  df <-
    as.data.frame(.data) %>%
    remove_rownames()
  if(!var %in% colnames(df)){
    stop("Variable '", var, "' not in data.", call. = FALSE)
  }
  rownames(df) <- df[[var]]
  df[[var]] <- NULL
  df
}
#' @name utils_rows_cols
#' @export
rownames_to_column <-  function(.data, var = "rowname"){
  col_names <- colnames(.data)
  if (var %in% col_names) {
    stop("Column `", var, "` already exists in `.data`.")
  }
  .data[, var] <- rownames(.data)
  rownames(.data) <- NULL
  .data[, c(var, setdiff(col_names, var))]
}
#' @name utils_rows_cols
#' @export
remove_rownames <- function(.data, ...){
  rownames(.data) <- NULL
  .data
}
#' @name utils_rows_cols
#' @export
column_exists <-function(.data, cols){
  if(length(setdiff(cols, colnames(.data))) != 0){
    FALSE
  } else{
    TRUE
  }
}
#' @name utils_rows_cols
#' @export
concatenate <- function(.data,
                        ...,
                        prefix = NULL,
                        suffix = NULL,
                        new_var = new_var,
                        sep = "_",
                        drop = FALSE,
                        pull = FALSE,
                        .before = NULL,
                        .after  = NULL){

  if (drop == FALSE){
    conc <- select(.data, ...)
    results <- mutate(.data,
                      {{new_var}} := paste(prefix,
                                           ifelse(is.null(prefix), "", sep),
                                           apply(conc, 1, paste, collapse = sep),
                                           ifelse(is.null(suffix), "", sep),
                                           suffix, sep = ""))
    if (pull == TRUE){
      results <- pull(results)
    }
    if (!missing(.before) | !missing(.after)){
      results <- reorder_cols(results, {{new_var}}, .before = {{.before}}, .after = {{.after}})
    }
  } else{
    conc <- select(.data, ...)
    results <- transmute(.data,
                         {{new_var}} := paste(prefix,
                                              ifelse(is.null(prefix), "", sep),
                                              apply(conc, 1, paste, collapse = sep),
                                              ifelse(is.null(suffix), "", sep),
                                              suffix, sep = ""))
    if (pull == TRUE){
      results <- pull(results)
    }
  }
  return(results)
}
#' @name utils_rows_cols
#' @export
get_levels <- function(.data, ...){
  if(missing(...)){
    df <-
      select(.data, everything()) %>% select_non_numeric_cols()
  } else{
    df <- select(.data, ...) %>% select_non_numeric_cols()
  }
  results <-
    df %>%
    as_factor(everything()) %>%
    map(levels)
  if(length(results) == 1){
    return(results[[1]])
  } else{
    return(results)
  }
}
#' @name utils_rows_cols
#' @export
get_levels_comb <- function(.data, ...){
  if(missing(...)){
    df <-
      select(.data, everything()) %>% select_non_numeric_cols()
  } else{
    df <- select(.data, ...) %>% select_non_numeric_cols()
  }
  df %>%
    as_factor() %>%
    map(levels) %>%
    expand.grid() %>%
    as_tibble()
}
#' @name utils_rows_cols
#' @importFrom dplyr count
#' @export
get_level_size <- function(.data, ...){
  return(n_by(.data, ...))
}
#' @name utils_rows_cols
#' @importFrom rlang quos
#' @export
reorder_cols <- function(.data, ..., .before = NULL, .after = NULL){
  return(relocate(.data, ..., .before = {{.before}}, .after = {{.after}}))
}
#' @name utils_rows_cols
#' @export
remove_cols <- function(.data, ...){
  select(.data, -c(!!!quos(...)))
}
#' @name utils_rows_cols
#' @export
remove_rows <- function(.data, ...){
  slice(.data, -c(!!!quos(...)))
}
#' @name utils_rows_cols
#' @export
select_first_col <- function(.data, offset = NULL){
  if(!missing(offset) && offset == 0){
    stop("`offset` must be greater than the 0", call. = FALSE)
  }
  offset <- ifelse(missing(offset), ncol(.data) - 1, ncol(.data) - offset)
  select(.data, last_col(offset = offset))
}
#' @name utils_rows_cols
#' @export
select_last_col <- function(.data, offset = NULL){
  offset <- ifelse(missing(offset), 0, offset)
  select(.data, last_col(offset = offset))
}
#' @name utils_rows_cols
#' @export
select_numeric_cols <- function(.data){
  if(is_grouped_df(.data)){
    .data <- ungroup(.data)
  }
  select_if(.data, is.numeric)
}
#' @name utils_rows_cols
#' @export
select_non_numeric_cols <- function(.data){
  if(is_grouped_df(.data)){
    .data <- ungroup(.data)
  }
  select_if(.data, ~!is.numeric(.x))
}
#' @name utils_rows_cols
#' @export
select_cols <- function(.data, ...){
  select(.data, ...)
}
#' @name utils_rows_cols
#' @export
select_rows <- function(.data, ...){
  slice(.data, ...)
}
#' @name utils_rows_cols
#' @export
tidy_colnames <- function(.data, sep = "_"){
  colnames(.data) <- tidy_strings(colnames(.data), sep = sep)
  .data
}





#' @title Useful functions for computing descriptive statistics
#' @name utils_stats
#' @description
#' `r badge('stable')`
#'
#' * **The following functions compute descriptive statistics by levels of
#' a factor or combination of factors quickly.**
#'    - `cv_by()` For computing coefficient of variation.
#'    - `max_by()` For computing maximum values.
#'    - `means_by()` For computing arithmetic means.
#'    - `min_by()` For compuing minimum values.
#'    - `n_by()` For getting the length.
#'    - `sd_by()` For computing sample standard deviation.
#'    - `var_by()` For computing sample variance.
#'    - `sem_by()` For computing standard error of the mean.
#'
#' * **Useful functions for descriptive statistics. All of them work
#' naturally with `\%>\%`, handle grouped data and multiple variables (all
#' numeric variables from `.data` by default).**
#'    - `av_dev()` computes the average absolute deviation.
#'    - `ci_mean_t()` computes the t-interval for the mean.
#'    - `ci_mean_z()` computes the z-interval for the mean.
#'    - `cv()` computes the coefficient of variation.
#'    - `freq_table()` Computes a frequency table for either numeric and
#'    categorical/discrete data. For numeric data, it is possible to define the
#'    number of classes to be generated.
#'    - `hmean(), gmean()` computes the harmonic and geometric means,
#' respectively. The harmonic mean is the reciprocal of the arithmetic mean of
#' the reciprocals. The geometric mean is the *n*th root of *n*
#' products.
#'    - `kurt()` computes the kurtosis like used in SAS and SPSS.
#'    - `range_data()` Computes the range of the values.
#'    - `n_valid()` The valid (not `NA`) length of a data.
#'    - `n_unique()` Number of unique values.
#'    - `n_missing()` Number of missing values.
#'    - `row_col_mean(), row_col_sum()` Adds a row with the mean/sum of
#'    each variable and a column with the the mean/sum for each row of the data.
#'    - `sd_amo(), sd_pop()` Computes sample and populational standard
#' deviation, respectively.
#'    - `sem()` computes the standard error of the mean.
#'    - `skew()` computes the skewness like used in SAS and SPSS.
#'    - `ave_dev()` computes the average of the absolute deviations.
#'    - `sum_dev()` computes the sum of the absolute deviations.
#'    - `sum_sq()` computes the sum of the squared values.
#'    - `sum_sq_dev()` computes the sum of the squared deviations.
#'    - `var_amo(), var_pop()` computes sample and populational variance.
#'
#' [desc_stat()] is wrapper function around the above ones and can be
#' used to compute quickly all these statistics at once.
#'
#' @details  The function `freq_table()` computes a frequency table for either
#'   numerical or categorical variables. If a variable is categorical or
#'   discrete (integer values), the number of classes will be the number of
#'   levels that the variable contains.
#'
#'If a variable (say, data) is continuous, the number of classes (k) is given by
#'the square root of the number of samples (n) if `n =< 100` or `5 * log10(n)`
#'if `n > 100`.
#'
#'The amplitude (\mjseqn{A}) of the data is used to define the size of the class (\mjseqn{c}),
#'given by
#'
#' \loadmathjax
#' \mjsdeqn{c = \frac{A}{n - 1}}
#'
#' The lower limit of the first class (LL1) is given by min(data) - c / 2. The
#' upper limit is given by LL1 + c. The limits of the other classes are given in
#' the same way. After the creation of the classes, the absolute and relative
#' frequencies within each class are computed.
#'
#' @param .data A data frame or a numeric vector.
#' @param ... The argument depends on the function used.
#' * For `*_by` functions, `...` is one or more categorical variables
#'  for grouping the data. Then the statistic required will be computed for all
#'  numeric variables in the data. If no variables are informed in `...`,
#'  the statistic will be computed ignoring all non-numeric variables in
#'  `.data`.
#' * For the other statistics, `...` is a comma-separated of unquoted
#'  variable names to compute the statistics. If no variables are informed in n
#'  `...`, the statistic will be computed for all numeric variables in
#'  `.data`.
#' @param .vars Used to select variables in the `*_by()` functions. One or more
#'   unquoted expressions separated by commas. Variable names can be used as if
#'   they were positions in the data frame, so expressions like `x:y` can be
#'   used to select a range of variables. Defaults to `everything()`.
#' @param na.rm If `FALSE`, the default, missing values are removed with a
#'   warning. If `TRUE`, missing values are silently removed.
#' @param var The variable to compute the frequency table. See `Details` for
#'   more details.
#' @param k The number of classes to be created. See `Details` for
#'   more details.
#' @param digits The number of significant figures to show. Defaults to 2.
#' @param table A frequency table computed with [freq_table()].
#' @param xlab,ylab The `x` and `y` labels.
#' @param fill,color The color to fill the bars and color the border of the bar,
#'   respectively.
#' @param ygrid Shows a grid line on the `y` axis? Defaults to `TRUE`.
#' freq_hist <- function(table,
#' @param level The confidence level for the confidence interval of the mean.
#'   Defaults to 0.95.
#' @return
#'  * Functions `*_by()` returns a `tbl_df` with the computed statistics by
#'  each level of the factor(s) declared in `...`.
#'  * All other functions return a named integer if the input is a data frame
#'  or a numeric value if the input is a numeric vector.
#'  * `freq_table()` Returns a list with the frequency table and the breaks used
#'  for class definition. These breaks can be used to construct an histogram of
#'  the variable.
#' @md
#' @references Ferreira, Daniel Furtado. 2009. Estatistica Basica. 2 ed. Vicosa,
#'   MG: UFLA.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(metan)
#' # means of all numeric variables by ENV
#' means_by(data_ge2, GEN, ENV)
#'
#' # Coefficient of variation for all numeric variables
#' # by GEN and ENV
#' cv_by(data_ge2, GEN, ENV)
#'
#' # Skewness of a numeric vector
#' set.seed(1)
#' nvec <- rnorm(200, 10, 1)
#' skew(nvec)
#'
#' # Confidence interval 0.95 for the mean
#' # All numeric variables
#' # Grouped by levels of ENV
#' data_ge2 %>%
#'   group_by(ENV) %>%
#'   ci_mean_t()
#'
#' # standard error of the mean
#' # Variable PH and EH
#' sem(data_ge2, PH, EH)
#'
#' # Frequency table for variable NR
#' data_ge2 %>%
#'   freq_table(NR)
#'}
#'
#' @export
av_dev <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    sum(abs(df - mean(df, na.rm = na.rm)), na.rm = na.rm) / length(which(!is.na(df)))
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
ci_mean_t <- function(.data, ..., na.rm = FALSE, level = 0.95) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    qt((0.5 + level/2), (length(which(!is.na(df))) - 1)) * sd(df, na.rm = na.rm)/sqrt(length(which(!is.na(df))))
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
ci_mean_z <- function(.data, ..., na.rm = FALSE, level = 0.95) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    qnorm((0.5 + level/2)) * sd(df, na.rm = na.rm)/sqrt(length(which(!is.na(df))))
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
cv <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    sd(df, na.rm = na.rm)/mean(df, na.rm = na.rm) * 100
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
freq_table <- function(.data, var, k = NULL, digits = 2){
  if(is_grouped_df(.data)){
    res <-
      metan::doo(.data,
                 ~freq_table(., {{var}}, k = k, digits = digits))
    freqs <-
      res %>%
      mutate(freqs = map(data, ~.x %>% .[["freqs"]])) |>
      unnest(freqs) |>
      remove_cols(data)
    breaks <-
      res |>
      mutate(freqs = map(data, ~.x %>% .[["breaks"]])) |>
      remove_cols(data)
    list_breaks <- breaks$freqs
    names(list_breaks) <- breaks$cor_grao
    return(list(freqs = freqs,
                breaks = breaks))

  } else{
    # function to create a frequence table with continuous variable
    # adapted from https://bendeivide.github.io/book-epaec/book-epaec.pdf
    freq_quant <- function(data, k = NULL, digits = digits){
      # the number of observations
      n <- length(data)

      # check the number of classes
      if(is.null(k)){
        if (n > 100) {
          k <- round(5 * log10(n), 0)
        } else{
          k <- round(sqrt(n), 0)
        }
      } else{
        k <- k
      }
      # data range
      rang <- range(data)
      # amplitude
      A <- diff(rang)
      # the size of the class
      c <- round(A / (k - 1), digits = digits)

      # lower and upper limit of the first class
      LI1 <- min(rang) - c / 2
      vi <- c(LI1,     rep(0, k - 1))
      vs <- c(LI1 + c, rep(0, k - 1))

      # build the other classes
      for (i in 2:k) {
        vi[i] <- vi[i - 1] + c
        vs[i] <- vs[i - 1] + c
      }
      vi <- round(vi, digits = digits)
      vs <- round(vs, digits = digits)
      # Find the frequency of each class
      freq <- function(x, vi, vs, k) {
        freq <- rep(0, k)
        for (i in 1:(k - 1)) {
          freq[i] <- length(x[x >= vi[i] & x < vs[i]])
        }
        freq[k] <- length(x[x >= vi[k] & x <= vs[k]])
        return(freq)
      }

      # absolute frequency
      fi <- freq(data, vi, vs, k)
      # check if any class is empty
      if(any(fi == 0)){
        warning("An empty class is not advised. Try to reduce the number of classes with the `k` argument", call. = FALSE)
      }
      # building the classes
      classe <- paste(vi, "|--- ", vs)
      classe[k] <- paste(vi[k], "|---|", vs[k])
      freqs <-
        data.frame(class = classe,
                   abs_freq = fi) |>
        mutate(abs_freq_ac = cumsum(abs_freq),
               rel_freq = abs_freq / sum(abs_freq),
               rel_freq_ac = cumsum(rel_freq))
      freqs[nrow(freqs) + 1, ] <- c("Total", sum(freqs[, 2]), sum(freqs[, 2]), 1, 1)
      freqs <-
        freqs |>
        as_numeric(2:5) |>
        round_cols(digits = digits)

      breaks <- sort(c(vi, vs))
      return(
        structure(
          list(freqs = freqs,
               LL = vi,
               UL = vs,
               vartype = "continuous"),
          class = "freq_table"
        )
      )
    }

    #check the class of the variable
    class_data <- .data |> pull({{var}}) |> class()
    # if variable is discrete or categorical
    if(class_data %in% c("character", "factor", "integer")){
      df <-
        .data %>%
        count({{var}}) |>
        as_character(1) |>
        mutate(abs_freq = n,
               abs_freq_ac = cumsum(abs_freq),
               rel_freq = abs_freq / sum(abs_freq),
               rel_freq_ac = cumsum(rel_freq)) |>
        remove_cols(n) |>
        as.data.frame()
      df[nrow(df) + 1, ] <- c("Total", sum(df[, 2]), sum(df[, 2]), 1, 1)
      df <- df |> as_numeric(2:5)
      return(
        structure(
          list(freqs = df,
               vartype = "categorical"),
          class = "freq_table"
        )
      )
    }
    # if variable is numeric
    if(class_data == "numeric"){
      data <- .data |> pull({{var}})
      # apply the function freq_quant in the numeric vector
      freq_quant(data, k = k, digits = digits)
    }
  }
}

#' @name utils_stats
#' @export
#' @importFrom graphics axis barplot grid plot.new plot.window rect title
freq_hist <- function(table,
                      xlab = NULL,
                      ylab = NULL,
                      fill = "gray",
                      color = "black",
                      ygrid = TRUE) {
  if (class(table) != "freq_table"){
    stop("Class of object `table` is not valid. Please use `freq_table()` to create a valid object.")
  }
  if (table$vartype == "categorical") {

    classes <- as.character(table$freqs[[1]])
    classes <- classes[-length(classes)]
    freqs <- table$freqs[[2]]
    freqs <- freqs[-length(freqs)]
    plot.new()

    if (is.null(xlab)) {
      xlab <- gettext("Groups")
    }
    if (is.null(ylab)) {
      ylab <- gettext("Frequency")
    }
    barplot(freqs ~ classes, xaxt = "n", xlab = "", ylab = "")

    if(isTRUE(ygrid)){
      grid(nx = NA, ny = NULL, col = "gray")
      opar <- par(new = TRUE)
      on.exit(par(opar))
    }
    barplot(freqs ~ classes,
            xlab = xlab,
            ylab = ylab)
  }
  if (table$vartype == "continuous") {
    xvar1 <- table$LL
    xvar2 <- table$UL
    freqs <- table$freqs$abs_freq
    yvar <- freqs[-length(freqs)]
    # Limiares
    xlim <- c(min(xvar1), max(xvar2))
    ylim <- c(0, 1.2 * max(yvar))

    # Area de plotagem
    plot.new()
    plot.window(xlim, ylim)

    # Labels
    if (is.null(xlab)) {
      xlab <- gettext("Classes")
    }
    if (is.null(ylab)) {
      ylab <- gettext("Frequency")
    }

    title(xlab = xlab, ylab = ylab)

    if(isTRUE(ygrid)){
      grid(nx = NA, ny = NULL, col = "gray")
      opar <- par(new = TRUE)
    }

    rect(xvar1,
         0,
         xvar2,
         yvar,
         col = fill,
         border = color)
    xvar <- c(xvar1, max(xvar2))
    axis(1, at = xvar)
    axis(2)
  }
}
#' @name utils_stats
#' @export
hmean <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    1 / mean(1 / df, na.rm = na.rm)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
gmean <- function(.data, ..., na.rm = FALSE){
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    exp(sum(log(df[df > 0]), na.rm = na.rm) / length(df))
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
kurt <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    n <- length(which(!is.na(df)))
    tmp <- sum((df - mean(df, na.rm = na.rm))^4, na.rm = na.rm)/(var(df, na.rm = na.rm))^2
    n * (n + 1) * tmp/((n - 1) * (n - 2) * (n - 3)) - 3 * (n - 1)^2/((n - 2) * (n - 3))
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
n_missing <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    sum(is.na(df))
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(everything(), funct)) %>%
        ungroup()
    } else{
      .data %>%
        summarise(across(c(...), funct))
    }
  }
}
#' @name utils_stats
#' @export
n_unique <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    length(unique(df))
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(everything(), funct)) %>%
        ungroup()
    } else{
      .data %>%
        summarise(across(c(...), funct))
    }
  }
}
#' @name utils_stats
#' @export
n_valid <- function(.data, ..., na.rm = FALSE){
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    length(which(!is.na(df)))
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
pseudo_sigma <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    IQR(df, na.rm = na.rm) / 1.35
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
range_data <- function(.data, ..., na.rm = FALSE){
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    max(df, na.rm = na.rm) - min(df, na.rm = na.rm)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
row_col_mean <- function(.data, na.rm = FALSE) {
  if(na.rm == FALSE & has_na(.data)){
    warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
    message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
    na.rm <- TRUE
  }
  if(!any(sapply(.data, is.numeric))){
    stop("All columns in '.data' must be numeric")
  }
  mat <- as.matrix(.data)
  row_means <- rowMeans(mat, na.rm = na.rm)
  col_means <- colMeans(mat, na.rm = na.rm)
  cmeans <- suppressWarnings(cbind(mat,  row_means) %>% rbind(col_means))
  rownames(cmeans) <- c(rownames(mat), "col_means")
  cmeans[nrow(cmeans), ncol(cmeans)] <- mean(mat, na.rm = na.rm)
  return(cmeans)
}
#' @name utils_stats
#' @export
row_col_sum <- function(.data, na.rm = FALSE) {
  if(na.rm == FALSE & has_na(.data)){
    warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
    message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
    na.rm <- TRUE
  }
  if(!any(sapply(.data, is.numeric))){
    stop("All columns in '.data' must be numeric")
  }
  mat <- as.matrix(.data)
  row_sums <- rowSums(mat, na.rm = na.rm)
  col_sums <- colSums(mat, na.rm = na.rm)
  cmeans <- suppressWarnings(cbind(mat,  row_sums) %>% rbind(col_sums))
  rownames(cmeans) <- c(rownames(mat), "col_sums")
  cmeans[nrow(cmeans), ncol(cmeans)] <- sum(mat, na.rm = na.rm)
  return(cmeans)
}
#' @name utils_stats
#' @export
sd_amo <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    sqrt(sum((df - mean(df, na.rm = na.rm))^2, na.rm = na.rm) / (length(which(!is.na(df))) - 1))
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
sd_pop <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    sqrt(sum((df - mean(df, na.rm = na.rm))^2, na.rm = na.rm) / length(which(!is.na(df))))
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
sem <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    sd(df, na.rm = na.rm) / sqrt(length(which(!is.na(df))))
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
skew <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    n <- length(which(!is.na(df)))
    sum((df - mean(df, na.rm = na.rm))^3, na.rm = na.rm)/sd(df, na.rm = na.rm)^3 * n / ((n - 1) * (n - 2))
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
sum_dev <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    sum(abs(df - mean(df, na.rm = na.rm)), na.rm = na.rm)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
ave_dev <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    mean(abs(df - mean(df, na.rm = na.rm)), na.rm = na.rm)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
sum_sq_dev <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    sum((df - mean(df, na.rm = na.rm))^2, na.rm = na.rm)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}

#' @name utils_stats
#' @export
sum_sq <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    sum(df ^ 2, na.rm = na.rm)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
var_pop <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    sd_pop(df, na.rm = na.rm)^2
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
var_amo <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    if(na.rm == FALSE & has_na(df)){
      warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
      message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
      na.rm <- TRUE
    }
    sd_amo(df, na.rm = na.rm)^2
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    } else{
      .data %>%
        select_cols(group_vars(.), ...) %>%
        summarise(across(where(is.numeric), funct)) %>%
        ungroup()
    }
  }
}

# main statistics, possible by one or more factors
#' @name utils_stats
#' @export
cv_by <- function(.data,
                  ...,
                  .vars = everything(),
                  na.rm = FALSE){
  if(na.rm == FALSE & has_na(.data)){
    warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
    message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
    na.rm <- TRUE
  }
  dfg <- group_by(.data, ...)
  dfg |>
    summarise(across(where(is.numeric), cv, na.rm = na.rm), .groups = "drop") |>
    select(group_vars(dfg), {{.vars}}) |>
    ungroup()
}
#' @name utils_stats
#' @export
max_by <- function(.data,
                   ...,
                   .vars = everything(),
                   na.rm = FALSE){
  if(na.rm == FALSE & has_na(.data)){
    warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
    message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
    na.rm <- TRUE
  }
  dfg <- group_by(.data, ...)
  dfg |>
    summarise(across(where(is.numeric), max, na.rm = na.rm), .groups = "drop") |>
    select(group_vars(dfg), {{.vars}}) |>
    ungroup()
}
#' @name utils_stats
#' @export
means_by <- function(.data,
                     ...,
                     .vars = everything(),
                     na.rm = FALSE){
  if(na.rm == FALSE & has_na(.data)){
    warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
    message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
    na.rm <- TRUE
  }
  dfg <- group_by(.data, ...)
  dfg |>
    summarise(across(where(is.numeric), mean, na.rm = na.rm), .groups = "drop") |>
    select(group_vars(dfg), {{.vars}}) |>
    ungroup()
}
#' @name utils_stats
#' @export
min_by <- function(.data,
                   ...,
                   .vars = everything(),
                   na.rm = FALSE){
  if(na.rm == FALSE & has_na(.data)){
    warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
    message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
    na.rm <- TRUE
  }
  dfg <- group_by(.data, ...)
  dfg |>
    summarise(across(where(is.numeric), min, na.rm = na.rm), .groups = "drop") |>
    select(group_vars(dfg), {{.vars}}) |>
    ungroup()
}
#' @name utils_stats
#' @export
n_by <- function(.data,
                 ...,
                 .vars = everything(),
                 na.rm = FALSE){
  if(na.rm == FALSE & has_na(.data)){
    warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
    message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
    na.rm <- TRUE
  }
  dfg <- group_by(.data, ...)
  dfg |>
    summarise(across(everything(), ~sum(!is.na(.))), .groups = "drop") |>
    select(group_vars(dfg), {{.vars}}) |>
    ungroup()
}
#' @name utils_stats
#' @export
sd_by <- function(.data,
                  ...,
                  .vars = everything(),
                  na.rm = FALSE){
  if(na.rm == FALSE & has_na(.data)){
    warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
    message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
    na.rm <- TRUE
  }
  dfg <- group_by(.data, ...)
  dfg |>
    summarise(across(where(is.numeric), sd, na.rm = na.rm), .groups = "drop") |>
    select(group_vars(dfg), {{.vars}}) |>
    ungroup()
}
#' @name utils_stats
#' @export
var_by <- function(.data,
                   ...,
                   .vars = everything(),
                   na.rm = FALSE){
  if(na.rm == FALSE & has_na(.data)){
    warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
    message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
    na.rm <- TRUE
  }
  dfg <- group_by(.data, ...)
  dfg |>
    summarise(across(where(is.numeric), var, na.rm = na.rm), .groups = "drop") |>
    select(group_vars(dfg), {{.vars}}) |>
    ungroup()
}
#' @name utils_stats
#' @export
sem_by <- function(.data,
                   ...,
                   .vars = everything(),
                   na.rm = FALSE){
  if(na.rm == FALSE & has_na(.data)){
    warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
    message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
    na.rm <- TRUE
  }
  dfg <- group_by(.data, ...)
  dfg |>
    summarise(across(where(is.numeric), sem, na.rm = na.rm), .groups = "drop") |>
    select(group_vars(dfg), {{.vars}}) |>
    ungroup()
}
#' @name utils_stats
#' @export
sum_by <- function(.data,
                   ...,
                   .vars = everything(),
                   na.rm = FALSE){
  if(na.rm == FALSE & has_na(.data)){
    warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
    message("To remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.")
    na.rm <- TRUE
  }
  dfg <- group_by(.data, ...)
  dfg |>
    summarise(across(where(is.numeric), sum, na.rm = na.rm), .groups = "drop") |>
    select(group_vars(dfg), {{.vars}}) |>
    ungroup()
}

#' @title Utilities for handling with matrices
#' @description
#' `r badge('stable')`
#'
#'  These functions help users to make upper, lower, or symmetric matrices
#'  easily.
#'
#' @name utils_mat
#' @param x A matrix to apply the function. It must be a symmetric (square)
#'   matrix in `make_upper_tri()` and `make_lower_tri()` or a
#'   triangular matrix in `make_sym()`. `tidy_sym()` accepts both
#'   symmetrical or triangular matrices.
#' @param lower A square matrix to fill the lower diagonal of the new matrix.
#' @param upper A square matrix to fill the upper diagonal of the new matrix.
#' @param diag What show in the diagonal of the matrix. Default to `NA`.
#' @param make The triangular to built. Default is `"upper"`. In this case,
#'   a symmetric matrix will be built based on the values of a lower triangular
#'   matrix.
#' @param keep_diag Keep diagonal values in the tidy data frame? Defaults to
#'   `TRUE`.
#' @details
#' * `make_upper_tri()` makes an upper triangular matrix using a symmetric
#' matrix.
#' * `make_lower_tri()` makes a lower triangular matrix using a symmetric
#' matrix.
#' * `make_sym()` makes a lower triangular matrix using a symmetric matrix.
#' * `tidy_sym()` transform a symmetric matrix into tidy data frame.
#' @md
#' @return An upper, lower, or symmetric matrix, or a tidy data frame.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' m <- cor(select_cols(data_ge2, 5:10))
#' make_upper_tri(m)
#' make_lower_tri(m)
#' make_lower_tri(m) %>%
#' make_sym(diag = 0)
#' tidy_sym(m)
#' tidy_sym(make_lower_tri(m))
#'
#' }
#'
make_upper_tri<-function(x, diag = NA){
  x[lower.tri(x)] <- NA
  diag(x) <- diag
  return(x)
}
#' @name utils_mat
#' @export
make_lower_tri<-function(x, diag = NA){
  x[upper.tri(x)] <- NA
  diag(x) <- diag
  return(x)
}
#' @name utils_mat
#' @export
make_lower_upper<-function(lower, upper, diag = NA){
  if(any(rownames(upper) == rownames(lower)) == FALSE){
    stop("Row names of 'upper' and 'lower' don't match.")
  }
  if(any(colnames(upper) == colnames(lower)) == FALSE){
    stop("Column names of 'upper' and 'lower' don't match.")
  }
  if (nrow(lower) != ncol(lower)) {
    stop("lower matrix must be square")
  }
  if (nrow(upper) != ncol(upper)) {
    stop("upper matrix must be square")
  }
  if (nrow(lower) != ncol(upper)) {
    stop("lower and upper matrices must have the same dimensions")
  }
  result <- matrix(NA, nrow = nrow(upper), ncol = ncol(upper))
  result[lower.tri(result)] <- t(lower)[lower.tri(t(lower))]
  result[upper.tri(result)] <- t(upper)[upper.tri(t(upper))]
  diag(result) <- diag
  rownames(result) <- colnames(result) <- rownames(lower)
  return(result)
}
#' @name utils_mat
#' @export
make_sym <- function(x, make = "upper", diag = NA) {
  if(make == "upper"){
    x[upper.tri(x)] <- t(x)[upper.tri(x)]
    diag(x) <- diag
  }
  if (make == "lower"){
    x[lower.tri(x)] <- t(x)[lower.tri(x)]
    diag(x) <- diag
  }
  return(x)
}
#' @name utils_mat
#' @export
tidy_sym <- function(x, keep_diag = TRUE){
  res <-
    x %>%
    as_tibble(rownames = "group1") %>%
    pivot_longer(names_to = "group2",
                 values_to = "value",
                 -group1) %>%
    remove_rows_na(verbose = FALSE) %>%
    arrange(group1)
  if(keep_diag == FALSE){
    res %<>%  remove_rows(which(res[1] == res[2]))
  }
  return(res)
}
#'Alternative to dplyr::do for doing anything
#'
#'
#' Provides an alternative to the `dplyr:do()` using `nest()`,
#' `mutate()` and `map()` to apply a function to a grouped data frame.
#'
#'  If the applied function returns a data frame, then the output will be
#'  automatically unnested. Otherwise, the output includes the grouping
#'  variables and a column named "data" , which is a "list-columns" containing
#'  the results for group combinations.
#'
#'@param .data a (grouped) data frame
#'@param .fun A function, formula, or atomic vector.
#'@param ... Additional arguments passed on to `.fun`
#'@param unnest Logical argument defaults to `TRUE` to control if results
#'  of `.fun` should be unnested. Valid only if the result is of class
#'  `data.frame` or `tbl_df`.
#' @export
#'@return a data frame
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'@importFrom purrr map
#'@importFrom tidyr nest unnest
#'@examples
#'\donttest{
#'library(metan)
#'# Head the first two lines of each environment
#'data_ge2 %>%
#'  group_by(ENV) %>%
#'  doo(~head(., 2))
#'
#' # Genotype analysis for each environment using 'gafem()'
#' # variable PH
#' data_ge2 %>%
#'   group_by(ENV) %>%
#'   doo(~gafem(., GEN, REP, PH, verbose = FALSE))
#'}
doo <- function(.data, .fun, ..., unnest = TRUE){
  if(is_grouped_df(.data)){
    out <- nest(.data)
  } else{
    out <- nest(.data, data = everything())
  }
  out %<>%
    ungroup() %>%
    mutate(data = purrr::map(.data$data, droplevels)) %>%
    mutate(data = purrr::map(.data$data, .fun, ...))
  if(has_class(out$data[[1]], c("tbl_df", "data.frame")) & unnest == TRUE){
    out <- unnest(out, cols = c(data))
  }
  if(is_grouped_df(.data)){
    out <- arrange(out, !!!syms(group_vars(.data)))
  }
  return(out)
}

#' @title Utilities for handling with classes
#' @name utils_class
#' @param x An object
#' @param class The class to add or remove
#' @details
#' * `add_class()`: add a class to the object `x` keeping all the other class(es).
#' * `has_class()`: Check if a class exists in object `x` and returns a logical value.
#' * `set_class()`: set a class to the object `x`.
#' * `remove_class()`: remove a class from the object `x`.
#' @md
#' @return The object `x` with the class added or removed.
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'@examples
#'\donttest{
#'library(metan)
#'df <-
#' data_ge2 %>%
#' add_class("my_class")
#'class(df)
#'has_class(df, "my_class")
#'remove_class(df, "my_class") %>% class()
#'set_class(df, "data_frame") %>% class()
#'}
#'
add_class <- function(x, class){
  class(x) <- unique(c(class(x), class))
  return(x)
}
#' @name utils_class
#' @export
has_class <- function(x, class){
  any(class(x)  %in%  class)
}
#' @name utils_class
#' @export
remove_class <- function(x, class){
  if(!class %in% class(x)){
    stop("Class not found in object ", match.call()[["x"]], call. = FALSE)
  }
  class(x) <- setdiff(class(x), class)
  return(x)
}
#' @name utils_class
#' @export
set_class <- function(x, class){
  class(x) <- class
  return(x)
}

#' Generate significance stars from p-values
#'
#' Generate significance stars from p-values using R's standard definitions.
#'
#' Mapping from p_value ranges to symbols:
#' * **0 - 0.0001**: '****'
#' * **0.0001 - 0.001**: '***'
#' * **0.001 - 0.01**: '**'
#' * **0.01 - 0.05**: '*'
#' * **0.05 - 1.0**: 'ns'
#' @md
#' @param p_value A numeric vector of p-values
#'
#' @return A character vector containing the same number of elements as p-value,
#'   with an attribute "legend" providing the conversion pattern.
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#'\donttest{
#' p_vals <- c(0.01, 0.043, 0.1, 0.0023, 0.000012)
#' stars_pval(p_vals)
#' }
#'
stars_pval <- function(p_value){
  unclass(
    symnum(
      p_value,
      corr = FALSE,
      na = FALSE,
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05,  1),
      symbols = c("****", "***", "**", "*",  "ns")
    )
  )
}



#' Check if a data set is balanced
#'
#' Check if a data set coming from multi-environment trials is balanced, i.e.,
#' all genotypes are in all environments.
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param resp The response variable.
#' @return A logical value
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#'\donttest{
#' unb <- data_ge %>%
#'         remove_rows(1:3) %>%
#'         droplevels()
#' is_balanced_trial(data_ge, ENV, GEN, GY)
#' is_balanced_trial(unb, ENV, GEN, GY)
#' }
#'

is_balanced_trial <- function(.data, env, gen, resp){
  mat <-
    .data %>%
    make_mat({{env}}, {{gen}}, {{resp}})
  if(has_na(mat) == TRUE){
    return(FALSE)
  } else{
    return(TRUE)
  }
}



#' Utilities for data Copy-Pasta
#' @name utils_data
#' @description
#' `r badge('stable')`
#'
#' These functions allows interacting with the system clipboard. It is possible
#' read from the clipboard or write a data frame or matrix to the clipboard.
#' * `clip_read()` read data from the clipboard.
#' * `clip_write()` write data to the clipboard.
#'
#' @param .data The data that should be copied to the clipboard. Only data frames and matrices are allowed
#' @param header If the copied data has a header row for dataFrame, defaults to
#'   `TRUE`.
#' @param sep The separator which should be used in the copied output.
#' @param row_names Decides if the output should keep row names or not, defaults
#'   to `FALSE`.
#' @param col_names Decides if the output should keep column names or not,
#'   defaults to `TRUE`.
#' @md
#' @param ... Further arguments to be passed to [utils::read.table()].
#' @export
#' @importFrom utils read.table write.table
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return Nothing
#'
clip_read <- function(header = TRUE, sep = "\t", ...){
  os <- get_os()
  if(os == "windows"){
    df <- read.table("clipboard", sep = sep, header = header, ...)
  }
  if(os == "osx"){
    df <- read.table(pipe("pbpaste"), sep = sep, header = header, ...)
  }
  if(os == "linux"){
    if (!file.exists(Sys.which("xclip")[1L])){
      stop("Cannot find xclip. Try installing it with 'sudo apt-get install xclip'")
    }
    df <- read.table(pipe("xclip -selection clipboard -o", open = "r"), ...)
  }
  return(as_tibble(df))
}
#' @name utils_data
#' @export
clip_write <- function(.data,
                       sep = "\t",
                       row_names = FALSE,
                       col_names = TRUE, ...){
  os <- get_os()
  if(os == "windows"){
    if(!has_class(.data , c("data.frame", "matrix"))){
      stop("Only data frames/tibbles/matrices allowed.")
    }
    write.table(.data,
                "clipboard",
                sep = sep,
                row.names = row_names,
                col.names = col_names,
                ...)
    message("Object '", match.call()[".data"], "' copied to the clipboard")
  }
  if(os == "osx"){
    if(!has_class(.data , c("data.frame", "matrix"))){
      stop("Only data frames/tibbles/matrices allowed.")
    }
    write.table(.data,
                file = pipe("pbcopy"),
                sep = sep,
                row.names = row_names,
                col.names = col_names,
                ...)
    message("Object '", match.call()[".data"], "' copied to the clipboard")
  }
  if(os == "linux"){
    if(!has_class(.data , c("data.frame", "matrix"))){
      stop("Only data frames/tibbles/matrices allowed.")
    }
    if (!file.exists(Sys.which("xclip")[1L])){
      stop("Cannot find xclip. Try installing it with 'sudo apt-get install xclip'")
    }
    write.table(.data,
                file = pipe("xclip -selection clipboard -i", open="w"),
                sep = sep,
                row.names = row_names,
                col.names = col_names,
                ...)
  }
}



### For internal use only ##
# Check labels
check_labels <- function(.data){
  if(any(sapply(.data, grepl, pattern = ":"))){
    stop("Using ':' in genotype or environment labels is not allowed. Use '_' instead.\ne.g., replace_string(data, ENV, pattern = ':', replacement = '_', new_var = ENV)", call. = FALSE)
  }
}
# Case
helper_case <- function(.data, fun, ...){
  if (has_class(.data, c("data.frame","tbl_df", "data.table"))){
    .data <- as.data.frame(.data)
    if(!missing(...)){
      mutate(.data, across(c(...), fun)) %>%
        as_tibble(rownames = NA)
    } else{
      mutate(.data, across(where(~!is.numeric(.x)), fun)) %>%
        as_tibble(rownames = NA)
    }
  } else{
    fun(.data)
  }
}
first_upper <- function(x){
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}
to_title <- function(x){
  gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2",
       all_lower_case(x),
       perl = TRUE)
}


# Get the OS
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else {
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os)){
      os <- "osx"
    }
    if (grepl("linux-gnu", R.version$os)){
      os <- "linux"
    }
  }
  return(all_lower_case(os))
}
# Coord radar for mtsi and mgidi indexes
coord_radar <- function (theta = "x", start = 0, direction = 1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- ifelse(theta == "x", "y", "x")
  ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start,
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}
# check for data frames in lists
has_df_in_list <- function(x){
  !any(
    is.na(
      map(x, class) %>%
        lapply(function(x){
          match('data.frame', x)
        })
    )
  )
}
# check if all are dataframes
all_df_in_list <- function(x){
  !all(
    is.na(
      map(x, class) %>%
        lapply(function(x){
          match('data.frame', x)
        })
    )
  )
}
# convert seconds to a hh:mm:ss format.
sec_to_hms <- function(t){
  paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0"),
        formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"),
        formatC(t %% 60, width = 2, format = "d", flag = "0"),
        sep = ":"
  )
}
# Embed a lifecycle badge in documentation
# Adapted from
badge <- function(stage){
  stages <- c("experimental", "stable", "superseded", "deprecated")
  if(!stage %in% stages){
    stop("'stage' must be one of the ", stages)
  }
  url <- paste0("https://lifecycle.r-lib.org/articles/stages.html#",
                stage)
  html <- sprintf("\\href{%s}{\\figure{%s}{options: alt='[%s]'}}",
                  url, file.path(sprintf("lifecycle-%s.svg", stage)),
                  first_upper_case(stage))
  text <- sprintf("\\strong{[%s]}", first_upper_case(stage))
  sprintf("\\ifelse{html}{%s}{%s}", html, text)
}


# Tools for deprecating functions and arguments
build_msg <- function(when, what, with, message){
  if(!grepl("\\(", what) && !grepl("\\)", what)){
    stop("must have function call syntax")
  }
  get_arg <- function(x){
    gsub("[\\(\\)]", "", regmatches(x, gregexpr("\\(.*?\\)", x))[[1]])
  }
  pkg <- ifelse(grepl("::", what), sub("::.*", "", what), ".")
  funct <- paste(sub("\\(.*", "", sub(".*::", "", what)), "()", sep = "")
  arg <- ifelse(!grepl("\\()\\b", what), get_arg(what) %>% remove_space(), NA)
  if(!is.na(arg)){
    msge <- paste0("Argument `", arg, "` of `", funct, "` is deprecated as of ", pkg, " ", when, ".")
    if(!is.null(with)){
      if(!grepl("\\(", with) && !grepl("\\)", with)){
        stop("'with' must have function call syntax", call. = FALSE)
      }
      pkg_with <- ifelse(grepl("::", with), sub("::.*", "", with), ".")
      funct_with <- paste(sub("\\(.*", "", sub(".*::", "", with)), "()", sep = "")
      with_mesg <-
        case_when(pkg_with == pkg & funct_with ==  funct ~ paste0("Please use the `", get_arg(with), "` argument instead."),
                  pkg_with == pkg & funct_with !=  funct ~ paste0("Please use the `", get_arg(with), "` argument of `", funct_with, "` instead." ),
                  pkg_with != pkg  & pkg_with == "." ~ paste0("Please use the `", get_arg(with), "` argument of `", funct_with, "` instead."),
                  pkg_with != pkg ~ paste0("Please use the `", get_arg(with), "` argument of `", paste0(pkg_with, "::", funct_with), "` instead."))
      msge <- paste0(msge, "\n", with_mesg)
    }
    msge <- ifelse(is.null(message), msge, paste0(msge, "\n", message))
  } else{
    msge <-  paste0("`",funct, "` is deprecated as of ", pkg, " ", when)
    if(!is.null(with)){
      if(!grepl("\\(", with) && !grepl("\\)", with)){
        stop("'with' must have function call syntax", call. = FALSE)
      }
      pkg_with <- ifelse(grepl("::", with), sub("::.*", "", with), NA)
      funct_with <- paste(sub("\\(.*", "", sub(".*::", "", with)), "()", sep = "")
      with_mesg <-
        case_when(pkg_with == pkg ~ paste0("Please use `", funct_with, "` instead"),
                  pkg_with != pkg ~ paste0("Please use the function `", funct_with, "` of the `",pkg_with, "` package."))
      with_mesg <- ifelse(is.null(message), with_mesg, paste(message))
      msge <- paste(msge, "\n", with_mesg, sep = "")
    }
    msge <- ifelse(is.null(message), msge, paste0(msge, "\n", message))
  }
  return(msge)
}
deprecated_warning <- function(when, what, with = NULL, message = NULL){
  warning(build_msg(when, what, with, message), call. = FALSE)
}
deprecated_error <- function(when, what, with = NULL, message = NULL){
  stop(build_msg(when, what, with, message), call. = FALSE)
}
deprecated <- function(){
  "missing_arg"
}
is_present <- function(x){
  x != "missing_arg"
}
