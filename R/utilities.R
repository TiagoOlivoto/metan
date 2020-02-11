#' @title Utilities for handling with numbers and strings
#' @name utils_num_str
#' @param .data A data frame
#' @param ... The argument depends on the function used.
#' * For \code{round_cols()} \code{...} are the variables to round. If no
#' variable is informed, all the numeric variables from \code{data} are used.
#' * For \code{all_lower_case()}, \code{all_upper_case()},
#' \code{all_title_case()}, \code{remove_strings()}, and \code{tidy_strings()}
#' \code{...} are the variables to apply the function. If no variable is
#' informed, the function will be applied to all non-numeric variables in
#' \code{.data}.
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
#' @param .before,.after For \code{replace_sting()}, \code{replace_number()},
#'   \code{extract_string()}, ,and  \code{extract_number()} one-based column
#'   index or column name where to add the new columns.
#' @param sep A character string to separate the terms. Defaults to "_".
#' @description
#' * \code{all_lower_case()}: Translate all non-numeric strings of a data frame
#' to lower case (
#'  \code{"Env"} to \code{"env"}).
#' * \code{all_upper_case()}: Translate all non-numeric strings of a data frame
#' to upper case (e.g., \code{"Env"} to \code{"ENV"}).
#' * \code{all_title_case()}: Translate all non-numeric strings of a data frame
#' to title case (e.g., \code{"ENV"} to \code{"Env"}).
#' * \code{extract_number()}: Extract the number(s) of a string.
#' * \code{extract_string()}: Extract all strings, ignoring case.
#' * \code{find_text_in_num()}: Find text characters in a numeric sequence and
#' return the row index.
#' * \code{has_text_in_num()}: Inspect columns looking for text in numeric
#' sequence and return a warning if text is found.
#' * \code{remove_strings()}: Remove all strings of a variable.
#' * \code{replace_number()}: Replace numbers with a replacement.
#' * \code{replace_string()}: Replace all strings with a replacement, ignoring
#' case.
#' * \code{round_cols()}: Round a selected column or a whole data frame to
#' significant figures.
#' * \code{tidy_strings()}: Tidy up characters strings, non-numeric columns, or
#' any selected columns in a data frame by putting all word in upper case,
#' replacing any space, tabulation, punctuation characters by \code{'_'}, and
#' putting \code{'_'} between lower and upper case. Suppose that \code{str =
#' c("Env1", "env 1", "env.1")} (which by definition should represent a unique
#' level in plant breeding trials, e.g., environment 1) is subjected to
#' \code{tidy_strings(str)}: the result will be then \code{c("ENV_1", "ENV_1",
#' "ENV_1")}. See Examples section for more examples.
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
#' ########## Extract, replace or remove strings ##########
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
#'
#' # Remove strings
#' remove_strings(data_ge)
#' remove_strings(data_ge, ENV)
#'
#'
#' ############ Find text in numeric sequences ###########
#' mixed_text <- data_ge
#' mixed_text[2, 4] <- "2..503"
#' mixed_text[3, 4] <- "3.2o75"
#' find_text_in_num(mixed_text, GY)
#'
#' ############# upper, lower and title cases ############
#'gen_text <- c("GEN 1", "Gen 1", "gen 1")
#'all_lower_case(gen_text)
#'all_upper_case(gen_text)
#'all_title_case(gen_text)
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
#' # Or a whole data frame
#' library(tibble)
#' df <- tibble(Env = messy_env,
#'              gen = messy_gen,
#'              Env_GEN = interaction(Env, gen),
#'              y = rnorm(6, 300, 10))
#' df
#' tidy_strings(df)
#' }
#' @export
all_upper_case <- function(.data, ...){
  if (any(class(.data) %in% c("data.frame","tbl_df", "data.table"))){
    .data <- as.data.frame(.data)
    if(!missing(...)){
      mutate_at(.data, vars(...), toupper) %>%
        as_tibble(rownames = NA)
    } else{
      mutate_if(.data, ~!is.numeric(.x), toupper) %>%
        as_tibble(rownames = NA)
    }
  } else{
    toupper(.data)
  }
}
#' @name utils_num_str
#' @export
all_lower_case <- function(.data, ...){
  if (any(class(.data) %in% c("data.frame","tbl_df", "data.table"))){
    .data <- as.data.frame(.data)
    if(!missing(...)){
      mutate_at(.data, vars(...), tolower) %>%
        as_tibble(rownames = NA)
    } else{
    mutate_if(.data, ~!is.numeric(.x), tolower) %>%
      as_tibble(rownames = NA)
    }
  } else{
    tolower(.data)
  }
}
#' @name utils_num_str
#' @export
all_title_case <- function(.data, ...){
  to_title <- function(x){
    x <- tolower(x)
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    return(x)
  }
  if (any(class(.data) %in% c("data.frame","tbl_df", "data.table"))){
    .data <- as.data.frame(.data)
    if(!missing(...)){
     mutate_at(.data, vars(...), to_title) %>%
      as_tibble(rownames = NA)
    } else{
      mutate_if(.data, ~!is.numeric(.x), to_title) %>%
        as_tibble(rownames = NA)
    }
  } else{
    return(to_title(.data))
  }
}

#' @name utils_num_str
#' @export
extract_number <- function(.data,
                           var,
                           new_var = new_var,
                           drop = FALSE,
                           pull = FALSE,
                           .before = NULL,
                           .after  = NULL){
  if (any(class(.data) %in% c("data.frame","tbl_df", "data.table"))){
  if (drop == FALSE){
    results <- .data %>%
      mutate({{new_var}} :=   as.numeric(gsub("[^0-9.-]+", "", as.character({{var}}))))
    if (pull == TRUE){
      results <- pull(results)
    }
    if (!is.null(.before) | !is.null(.after)){
      results <- reorder_cols(results, {{new_var}}, .before = .before, .after = .after)
    }
  } else{
    results <- .data %>%
      transmute({{new_var}} :=   as.numeric(gsub("[^0-9.-]+", "", as.character({{var}}))))
    if(pull == TRUE){
      results <- pull(results)
    }
  }
  return(results)
  } else{
    as.numeric(gsub("[^0-9.-]+", "", .data))
  }
}
#' @name utils_num_str
#' @export
extract_string <- function(.data,
                           var,
                           new_var = new_var,
                           drop = FALSE,
                           pull = FALSE,
                           .before = NULL,
                           .after  = NULL){
  if (any(class(.data) %in% c("data.frame","tbl_df", "data.table"))){
  if (drop == FALSE){
    results <- .data %>%
      mutate({{new_var}} := as.character(gsub("[^A-z.-]+", "", as.character({{var}}))))
    if(pull == TRUE){
      results <- pull(results)
    }
    if (!is.null(.before) | !is.null(.after)){
      results <- reorder_cols(results, {{new_var}}, .before = .before, .after = .after)
    }
  } else{
    results <- .data %>%
      transmute({{new_var}} := as.character(gsub("[^A-z.-]+", "", as.character({{var}}))))
    if(pull == TRUE){
      results <- pull(results)
    }
  }
  return(results)
  } else{
    as.character(gsub("[^A-z.-]+", "", .data))
  }
}
#' @name utils_num_str
#' @export
find_text_in_num <- function(.data, ...){
  if(!missing(...)){
    .data <- select_cols(.data, ...)
    if(ncol(.data)>1){
      stop("Only one variable is accepted.", call. = FALSE)
    }
    .data %<>% pull()
  }
  which(is.na(suppressWarnings(as.numeric(.data))))
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
remove_strings <- function(.data, ...){
  if (any(class(.data) %in% c("data.frame","tbl_df", "data.table"))){
  if(missing(...)){
    vars <- vars(everything())
  } else{
    vars <- vars(...)
  }
  results <-
    mutate_at(.data,
              .vars = vars,
              .funs = gsub,
              pattern = "[A-z]",
              replacement = "") %>%
    mutate_at(.vars = vars,
              .funs = as.numeric)
  return(results)
  } else{
    return(gsub(pattern, replacement, .data))
  }
}
#' @name utils_num_str
#' @export
replace_number <- function(.data,
                           var,
                           new_var = new_var,
                           pattern = NULL,
                           replacement = "",
                           drop = FALSE,
                           pull = FALSE,
                           .before = NULL,
                           .after  = NULL){
  if(missing(pattern)){
    pattern <- "[0-9]"
  } else{
    pattern <- pattern
  }
  if (any(class(.data) %in% c("data.frame","tbl_df", "data.table"))){
  if (drop == FALSE){
    results <- .data %>%
      mutate({{new_var}} := gsub(pattern, replacement, as.character({{var}})))
    if(pull == TRUE){
      results <- pull(results)
    }
    if (!is.null(.before) | !is.null(.after)){
      results <- reorder_cols(results, {{new_var}}, .before = .before, .after = .after)
    }
  } else{
    results <- .data %>%
      transmute({{new_var}} := gsub(pattern, replacement, as.character({{var}})))
    if(pull == TRUE){
      results <- pull(results)
    }
  }
  return(results)
  } else{
    return(gsub(pattern, replacement, .data))
  }
}
#' @name utils_num_str
#' @export
replace_string <- function(.data,
                           var,
                           new_var = new_var,
                           pattern = NULL,
                           replacement = "",
                           drop = FALSE,
                           pull = FALSE,
                           .before = NULL,
                           .after  = NULL){
  if(missing(pattern)){
    pattern <- "[A-z]"
  } else {
    pattern <- pattern
  }
  if (any(class(.data) %in% c("data.frame","tbl_df", "data.table"))){
  if (drop == FALSE){
    results <- .data %>%
      mutate({{new_var}} := gsub(pattern, replacement, as.character({{var}})))
    if(pull == TRUE){
      results <- pull(results)
    }
    if (!is.null(.before) | !is.null(.after)){
      results <- reorder_cols(results, {{new_var}}, .before = .before, .after = .after)
    }
  } else{
    results <- .data %>%
      transmute({{new_var}} := gsub(pattern, replacement, as.character({{var}})))
    if(pull == TRUE){
      results <- pull(results)
    }
  }
  return(results)
  } else{
    return(gsub(pattern, replacement, .data))
  }
}
#' @name utils_num_str
#' @export
#' @importFrom dplyr mutate_if transmute
round_cols <- function(.data, ...,  digits = 2){
  rn_test <- has_rownames(.data)
  if(rn_test == TRUE){
    rnames <- rownames(.data)
  }
  if (missing(...)){
    .data %<>% mutate_if(is.numeric, round, digits = digits)
  } else{
    .data %<>% mutate_at(vars(...), round, digits = digits)
  }
  if(rn_test == TRUE){
    rownames(.data) <- rnames
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
  if (any(class(.data) %in% c("data.frame","tbl_df", "data.table"))){
    if(missing(...)){
      results <-
        mutate_if(.data, ~!is.numeric(.x),  fstr)
    } else{
      results <-
        mutate_at(.data, vars(...), fstr)
    }
    return(results)
  } else{
    return(fstr(.data))
  }
}



#' @title Utilities for handling with rows and columns
#' @name utils_rows_cols
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
#' * \code{all_pairs()}: Get all the possible pairs between the levels of a
#' factor.
#' * \code{colnames_to_lower}: Translate all column names to lower case.
#' * \code{colnames_to_upper}: Translate all column names to upper case.
#' * \code{colnames_to_title}: Translate all column names to title case.
#' * \code{column_exists()}: Checks if a column exists in a data frame. Return a
#' logical value.
#' * \code{columns_to_first()}: Move columns to first positions in \code{.data}.
#' * \code{columns_to_last()}: Move columns to last positions in \code{.data}.
#' * \code{concatenate()}: Concatenate columns of a data frame. If \code{drop =
#' TRUE} then the existing variables are dropped. If \code{pull = TRUE} then the
#' concatenated variable is pull out to a vector. This is specially useful when
#' using \code{concatenate} to add columns to a data frame with \code{add_cols}.
#' * \code{get_levels()}: Get the levels of a factor variable.
#' * \code{get_level_size()}: Get the size of each level of a factor variable.
#' * \code{remove_cols()}: Remove one or more columns from a data frame.
#' * \code{remove_rows()}: Remove one or more rows from a data frame.
#' * \code{reorder_cols()}: Reorder columns in a data frame.
#' * \code{select_cols()}: Select one or more columns from a data frame.
#' * \code{select_first_col()}: Select first variable, possibly with an offset.
#' * \code{select_last_col()}: Select last variable, possibly with an offset.
#' * \code{select_numeric_cols()}: Select all the numeric columns of a data
#' frame.
#' * \code{select_non_numeric_cols()}: Select all the non-numeric columns of a
#' data frame.
#' * \code{select_rows()}: Select one or more rows from a data frame.
#'
#' @param .data A data frame
#' @param ... The argument depends on the function used.
#' * For \code{add_cols()} and \code{add_rows()} is name-value pairs. All values
#' must have one element for each row in \code{.data} when using
#' \code{add_cols()} or one element for each column in \code{.data} when using
#' \code{add_rows()}. Values of length 1 will be recycled when using
#' \code{add_cols()}.
#'
#' * For \code{remove_cols()} and \code{select_cols()},  \code{...} is the
#' column name or column index of the variable(s) to be dropped.
#'
#' * For \code{columns_to_first()} and \code{columns_to_last()},  \code{...} is
#' the column name or column index of the variable(s) to be moved to first or
#' last in \code{.data}.
#'
#' * For \code{remove_rows()} and \code{select_rows()}, \code{...} is an integer
#' row value.
#'
#' * For \code{concatenate()}, \code{...} is the unquoted variable names to be
#' concatenated.
#' @param .before,.after For \code{add_cols()}, \code{concatenate()}, and
#'   \code{reorder_cols()}, one-based column index or column name where to add
#'   the new columns, default: .after last column. For \code{add_rows()},
#'   one-based row index where to add the new rows, default: .after last row.
#' @param new_var The name of the new variable containing the concatenated
#'   values. Defaults to \code{new_var}.
#' @param sep The separator to appear between concatenated variables. Defaults
#'   to "_".
#' @param drop Logical argument. If \code{TRUE} keeps the new variable
#'   \code{new_var} and drops the existing ones. Defaults to \code{FALSE}.
#' @param pull Logical argument. If \code{TRUE}, returns the last column (on the
#'   assumption that's the column you've created most recently), as a vector.
#' @param cols A quoted variable name to check if it exists in \code{.data}.
#' @param group A factor variable to get the levels.
#' @param levels The levels of a factor or a numeric vector.
#' @param offset Set it to \emph{n} to select the \emph{n}th variable from the
#'   end (for \code{select_last_col()}) of from the begin (for
#'   \code{select_first_col()})
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
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
#' ############### Reordering columns ###############
#' reorder_cols(data_ge2, NKR, .before = "ENV")
#' reorder_cols(data_ge2, ENV, GEN, .after = "ED")
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
#'           .after = "GEN") %>%
#'  replace_string(ENV_GEN,
#'                 pattern = "H",
#'                 replacement = "HYB_",
#'                 .after = "ENV_GEN")
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
#'   add_rows(ENV = "E_TEST",
#'            GEN = "G_TEST",
#'            REP = 3,
#'            GY = 10.3,
#'            HM = 100.11,
#'            .after = 1)
#'
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
  results <- mutate(.data, ...)
  pos_dots <- (ncol(results)) - ((ncol(results) - ncol(.data)) - 1)
  if (!is.null(.before) | !is.null(.after)){
    results <- reorder_cols(results,
                            pos_dots:ncol(results),
                            .before = .before,
                            .after = .after)
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
all_pairs <- function(.data, levels){
  levels <-
    get_levels(.data, {{levels}})
  combn(levels, 2) %>%
    t() %>%
    as.data.frame()
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
                        new_var = new_var,
                        sep = "_",
                        drop = FALSE,
                        pull = FALSE,
                        .before = NULL,
                        .after  = NULL){
  if (drop == FALSE){
    conc <- select(.data, ...)
    results <- mutate(.data, {{new_var}} := apply(conc, 1, paste, collapse = sep))
    if (pull == TRUE){
      results <- pull(results)
    }
    if (!is.null(.before) | !is.null(.after)){
      results <- reorder_cols(results, {{new_var}}, .before = .before, .after = .after)
    }
  } else{
    conc <- select(.data, ...)
    results <- transmute(.data, {{new_var}} := apply(conc, 1, paste, collapse = sep))
    if (pull == TRUE){
      results <- pull(results)
    }
  }
  return(results)
}
#' @name utils_rows_cols
#' @export
get_levels <- function(.data, group){
  .data %>%
    mutate_if(~!is.numeric(.x), as.factor) %>%
    pull({{group}}) %>%
    levels()
}
#' @name utils_rows_cols
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
#' @name utils_rows_cols
#' @importFrom rlang quos
#' @export
reorder_cols <- function(.data, ..., .before = NULL, .after = NULL){
  args <- match.call()
  if (is.null(.after) && is.null(.before)){
    stop("At least '.before' or '.after' must be informed", call. = FALSE)
  }
  if (!is.null(.after) && !is.null(.before)){
    stop("'.before' and '.after' cannot be evaluated together.", call. = FALSE)
  }
  if (!is.null(.before)){
    if(is.character(.before)){
      if(!(.before %in% colnames(.data))){
        stop(paste("Column", args[[".before"]], "not present in .data"), call. = FALSE)
      }
    } else{
      .before <- colnames(.data[.before])
    }
    sel <- select(.data, ...)
    bfr <- .data[,1:which(colnames(.data) ==  .before)-1]
    if(any(colnames(sel) %in% colnames(bfr))){
      bfr <- bfr[, -which(names(bfr)  %in% names(sel))]
    }
    aft <- select(.data, -!!colnames(bfr), -!!colnames(sel))
    results <- cbind(bfr, sel, aft)
  }
  if(!is.null(.after)) {
    if(is.character(.after)){
      if(!(.after %in% colnames(.data))){
        stop(paste("Column", args[[".after"]], "not present in .data"), call. = FALSE)
      }
    } else{
      .after <- colnames(.data[.after])
    }
    if(which(colnames(.data) ==  .after)+1 > ncol(.data)){
      results <- select(.data, -c(!!!quos(...)), everything())
    } else{
      sel <- select(.data, ...)
      aft <- .data[,(which(colnames(.data) ==  .after)+1):ncol(.data)]
      if(any(colnames(sel) %in% colnames(aft))){
        aft <- aft[, -which(names(aft)  %in% names(sel))]
      }
      bfr <- select(.data, -!!colnames(aft), -!!colnames(sel))
      results <- cbind(bfr, sel, aft)
    }
  }
  return(as_tibble(results))
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





#' @title Useful functions for computing descriptive statistics
#' @name utils_stats
#' @description
#' * \strong{The following functions compute descriptive statistics by levels of
#' a factor or combination of factors quickly.}
#'    - \code{cv_by()} For computing coefficient of variation.
#'    - \code{max_by()} For computing maximum values.
#'    - \code{means_by()} For computing arithmetic means.
#'    - \code{min_by()} For compuing minimum values.
#'    - \code{n_by()} For getting the length.
#'    - \code{sd_by()} For computing sample standard deviation.
#'    - \code{sem_by()} For computing standard error of the mean.
#'
#' * \strong{Useful functions for descriptive statistics. All of them work
#' naturally with \code{\%>\%}, handle grouped data and multiple variables (all
#' numeric variables from \code{.data} by default).}
#'    - \code{av_dev()} computes the average absolute deviation.
#'    - \code{ci_mean()} computes the confidence interval for the mean.
#'    - \code{cv()} computes the coefficient of variation.
#'    - \code{freq_table()} Computes frequency fable. Handles grouped data.
#'    - \code{hm_mean(), gm_mean()} computes the harmonic and geometric means,
#' respectively. The harmonic mean is the reciprocal of the arithmetic mean of
#' the reciprocals. The geometric mean is the \emph{n}th root of \emph{n}
#' products.
#'    - \code{kurt()} computes the kurtosis like used in SAS and SPSS.
#'    - \code{range_data()} Computes the range of the values.
#'    - \code{sd_amo(), sd_pop()} Computes sample and populational standard
#' deviation, respectively.
#'    - \code{sem()} computes the standard error of the mean.
#'    - \code{skew()} computes the skewness like used in SAS and SPSS.
#'    - \code{sum_dev()} computes the sum of the absolute deviations.
#'    - \code{sum_sq_dev()} computes the sum of the squared deviations.
#'    - \code{var_amo(), var_pop()} computes sample and populational variance.
#'    - \code{valid_n()} Return the valid (not \code{NA}) length of a data.
#'
#' \code{\link{desc_stat}} is wrapper function around the above ones and can be
#' used to compute quickly all these statistics at once.
#'
#' @param .data A data frame or a numeric vector.
#' @param ... The argument depends on the function used.
#'  * For \code{*_by} functions, \code{...} is one or more categorical variables
#'  for grouping the data. Then the statistic required will be computed for all
#'  numeric variables in the data. If no variables are informed in \code{...},
#'  the statistic will be computed ignoring all non-numeric variables in
#'  \code{.data}.
#'  * For the other statistics, \code{...} is a comma-separated of unquoted
#'  variable names to compute the statistics. If no variables are informed in n
#'  \code{...}, the statistic will be computed for all numeric variables in
#'  \code{.data}.
#'
#' @param na.rm A logical value indicating whether \code{NA} values should be
#'    stripped before the computation proceeds. Defaults to \code{FALSE}.
#' @param level The confidence level for the confidence interval of the mean.
#'   Defaults to 0.95.
#' @return
#'  * Functions \code{*_by()} returns a tbl_df with the computed statistics by
#'  each level of the factor(s) declared in \code{...}.
#'  * All other functions return a nammed integer if the input is a data frame
#'  or a numeric value if the input is a numeric vector.
#' @md
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
#'   ci_mean()
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
    sum(abs(df - mean(df, na.rm = na.rm)), na.rm = na.rm) / length(which(!is.na(df)))
  }
  if(has_na(.data) && na.rm == FALSE){
    stop("NA values in data. Use 'na.rm = TRUE' to remove NAs from analysis.\nTo remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.", call. = FALSE)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    } else{
      .data %>%
       select_cols(group_vars(.), ...) %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
ci_mean <- function(.data, ..., na.rm = FALSE, level = 0.95) {
  funct <- function(df){
    qt((0.5 + level/2), (length(which(!is.na(df))) - 1)) * sd(df, na.rm = na.rm)/sqrt(length(which(!is.na(df))))
  }
  if(has_na(.data) && na.rm == FALSE){
    stop("NA values in data. Use 'na.rm = TRUE' to remove NAs from analysis.\nTo remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.", call. = FALSE)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    } else{
      .data %>%
       select_cols(group_vars(.), ...) %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
cv <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    sd(df, na.rm = na.rm)/mean(df, na.rm = na.rm) * 100
  }
  if(has_na(.data) && na.rm == FALSE){
    stop("NA values in data. Use 'na.rm = TRUE' to remove NAs from analysis.\nTo remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.", call. = FALSE)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    } else{
      .data %>%
       select_cols(group_vars(.), ...) %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
freq_table <- function(.data, ...){
  if(is_grouped_df(.data)){
    dplyr::do(.data, freq_table(., ...))
  }
  .data %>%
    count(...) %>%
    mutate(rel_freq = n / sum(n),
           cum_freq = cumsum(rel_freq))
}
#' @name utils_stats
#' @export
hm_mean <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    1 / mean(1 / df, na.rm = na.rm)
  }
  if(has_na(.data) && na.rm == FALSE){
    stop("NA values in data. Use 'na.rm = TRUE' to remove NAs from analysis.\nTo remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.", call. = FALSE)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    } else{
      .data %>%
       select_cols(group_vars(.), ...) %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
gm_mean <- function(.data, ..., na.rm = FALSE){
  funct <- function(df){
    exp(sum(log(df[df > 0]), na.rm = na.rm) / length(df))
  }
  if(has_na(.data) && na.rm == FALSE){
    stop("NA values in data. Use 'na.rm = TRUE' to remove NAs from analysis.\nTo remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.", call. = FALSE)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    } else{
      .data %>%
       select_cols(group_vars(.), ...) %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
kurt <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    n <- length(which(!is.na(df)))
    tmp <- sum((df - mean(df, na.rm = na.rm))^4, na.rm = na.rm)/(var(df, na.rm = na.rm))^2
    n * (n + 1) * tmp/((n - 1) * (n - 2) * (n - 3)) - 3 * (n - 1)^2/((n - 2) * (n - 3))
  }
  if(has_na(.data) && na.rm == FALSE){
    stop("NA values in data. Use 'na.rm = TRUE' to remove NAs from analysis.\nTo remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.", call. = FALSE)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    } else{
      .data %>%
       select_cols(group_vars(.), ...) %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
range_data <- function(.data, ..., na.rm = FALSE){
  funct <- function(df){
    max(df, na.rm = na.rm) - min(df, na.rm = na.rm)
  }
  if(has_na(.data) && na.rm == FALSE){
    stop("NA values in data. Use 'na.rm = TRUE' to remove NAs from analysis.\nTo remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.", call. = FALSE)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    } else{
      .data %>%
       select_cols(group_vars(.), ...) %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
sd_amo <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    sqrt(sum((df - mean(df, na.rm = na.rm))^2, na.rm = na.rm) / (length(which(!is.na(df))) - 1))
  }
  if(has_na(.data) && na.rm == FALSE){
    stop("NA values in data. Use 'na.rm = TRUE' to remove NAs from analysis.\nTo remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.", call. = FALSE)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    } else{
      .data %>%
       select_cols(group_vars(.), ...) %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
sd_pop <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    sqrt(sum((df - mean(df, na.rm = na.rm))^2, na.rm = na.rm) / length(which(!is.na(df))))
  }
  if(has_na(.data) && na.rm == FALSE){
    stop("NA values in data. Use 'na.rm = TRUE' to remove NAs from analysis.\nTo remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.", call. = FALSE)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    } else{
      .data %>%
       select_cols(group_vars(.), ...) %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
sem <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    sd(df, na.rm = na.rm) / sqrt(length(which(!is.na(df))))
  }
  if(has_na(.data) && na.rm == FALSE){
    stop("NA values in data. Use 'na.rm = TRUE' to remove NAs from analysis.\nTo remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.", call. = FALSE)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    } else{
      .data %>%
       select_cols(group_vars(.), ...) %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
skew <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    n <- length(which(!is.na(df)))
    sum((df - mean(df, na.rm = na.rm))^3, na.rm = na.rm)/sd(df, na.rm = na.rm)^3 * n / ((n - 1) * (n - 2))
  }
  if(has_na(.data) && na.rm == FALSE){
    stop("NA values in data. Use 'na.rm = TRUE' to remove NAs from analysis.\nTo remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.", call. = FALSE)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    } else{
      .data %>%
       select_cols(group_vars(.), ...) %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
sum_dev <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    sum(abs(df - mean(df, na.rm = na.rm)), na.rm = na.rm)
  }
  if(has_na(.data) && na.rm == FALSE){
    stop("NA values in data. Use 'na.rm = TRUE' to remove NAs from analysis.\nTo remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.", call. = FALSE)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    } else{
      .data %>%
       select_cols(group_vars(.), ...) %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
sum_sq_dev <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    sum((df - mean(df, na.rm = na.rm))^2, na.rm = na.rm)
  }
  if(has_na(.data) && na.rm == FALSE){
    stop("NA values in data. Use 'na.rm = TRUE' to remove NAs from analysis.\nTo remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.", call. = FALSE)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    } else{
      .data %>%
       select_cols(group_vars(.), ...) %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
var_pop <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    sd_pop(df, na.rm = na.rm)^2
  }
  if(has_na(.data) && na.rm == FALSE){
    stop("NA values in data. Use 'na.rm = TRUE' to remove NAs from analysis.\nTo remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.", call. = FALSE)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    } else{
      .data %>%
       select_cols(group_vars(.), ...) %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
var_amo <- function(.data, ..., na.rm = FALSE) {
  funct <- function(df){
    sd_amo(df, na.rm = na.rm)^2
  }
  if(has_na(.data) && na.rm == FALSE){
    stop("NA values in data. Use 'na.rm = TRUE' to remove NAs from analysis.\nTo remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.", call. = FALSE)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    } else{
      .data %>%
       select_cols(group_vars(.), ...) %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    }
  }
}
#' @name utils_stats
#' @export
valid_n <- function(.data, ..., na.rm = FALSE){
  funct <- function(df){
    length(which(!is.na(df)))
  }
  if(has_na(.data) && na.rm == FALSE){
    stop("NA values in data. Use 'na.rm = TRUE' to remove NAs from analysis.\nTo remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.", call. = FALSE)
  }
  if(is.null(nrow(.data))){
    funct(.data)
  } else{
    if(missing(...)){
      .data %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    } else{
      .data %>%
       select_cols(group_vars(.), ...) %>%
        summarise_if(is.numeric, funct) %>%
        ungroup()
    }
  }
}
# main statistics, possible by one or more factors
#' @name utils_stats
#' @export
cv_by <- function(.data, ...){
  group_by(.data, ...) %>%
    summarise_if(is.numeric, cv) %>%
    ungroup()
}
#' @name utils_stats
#' @export
max_by <- function(.data, ...){
  group_by(.data, ...) %>%
    summarise_if(is.numeric, max) %>%
    ungroup()
}
#' @name utils_stats
#' @export
means_by <- function(.data, ...){
  group_by(.data, ...) %>%
    summarise_if(is.numeric, mean) %>%
    ungroup()
}
#' @name utils_stats
#' @export
min_by <- function(.data, ...){
  group_by(.data, ...) %>%
    summarise_if(is.numeric, min) %>%
    ungroup()
}
#' @name utils_stats
#' @export
n_by <- function(.data, ...){
    group_by(.data, ...) %>%
    summarise_all(list(~sum(!is.na(.)))) %>%
    ungroup()
}
#' @name utils_stats
#' @export
sd_by <- function(.data, ...){
  group_by(.data, ...) %>%
    summarise_if(is.numeric, sd) %>%
    ungroup()
}
#' @name utils_stats
#' @export
sem_by <- function(.data, ...){
  group_by(.data, ...) %>%
    summarise_if(is.numeric, sem) %>%
    ungroup()
}
#' @title Utilities for handling with matrices
#'
#' @description These functions help users to make upper, lower, or symmetric
#'   matrices easily.
#'
#' @name utils_mat
#' @param x A matrix to apply the function. It must be a symmetric (square)
#'   matrix in \code{make_upper_tri()} and \code{make_lower_tri()} or a
#'   triangular matrix in \code{make_sym()}. \code{tidy_sym()} accepts both
#'   symmetrical or triangular matrices.
#' @param diag What show in the diagonal of the matrix. Default to \code{NA}.
#' @param make The triangular to built. Default is \code{"upper"}. In this case,
#'   a symmetric matrix will be built based on the values of a lower triangular
#'   matrix.
#' @param keep_diag Keep diagonal values in the tidy data frame? Defaults to
#'   \code{TRUE}.
#' @details
#' * \code{make_upper_tri()} makes an upper triangular matrix using a symmetric
#' matrix.
#' * \code{make_lower_tri()} makes a lower triangular matrix using a symmetric
#' matrix.
#' * \code{make_sym()} makes a lower triangular matrix using a symmetric matrix.
#' * \code{tidy_sym()} transform a symmetric matrix into tidy data frame.
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
#' Provides an alternative to the \code{dplyr:do()} using \code{nest()},
#' \code{mutate()} and \code{map()} to apply a function to a grouped data frame.
#'
#'  If the applied function returns a data frame, then the output will be
#'  automatically unnested. Otherwise, the output includes the grouping
#'  variables and a column named "data" , which is a "list-columns" containing
#'  the results for group combinations.
#'
#'@param .data a (grouped) data frame
#'@param .fun A function, formula, or atomic vector.
#'@param ... Additional arguments passed on to \code{.fun}
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'@return a data frame
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
doo <- function(.data, .fun, ...){
  if(is_grouped_df(.data)){
    out <- nest(.data)
  } else{
    out <- nest(.data, data = everything())
  }
  out %<>%
    ungroup() %>%
    mutate(data = purrr::map(.data$data, droplevels)) %>%
    mutate(data = purrr::map(.data$data, .fun, ...))
  if(inherits(out$data[[1]], c("tbl_df"))){
    out <- suppressWarnings(unnest(out))
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
#' * \code{add_class()}: add a class to the object \code{x} keeping all the other class(es).
#' * \code{set_class()}: set a class to the object \code{x}.
#' * \code{remove_class()}: remove a class from the object \code{x}.
#' @md
#' @return The object \code{x} with the class added or removed.
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'@examples
#'\donttest{
#'library(metan)
#'df <-
#' data_ge2 %>%
#' add_class("my_class")
#'class(df)
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
#' * \strong{0 - 0.0001}: '****'
#' * \strong{0.0001 - 0.001}: '***'
#' * \strong{0.001 - 0.01}: '**'
#' * \strong{0.01 - 0.05}: '*'
#' * \strong{0.05 - 1.0}: 'ns'
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
