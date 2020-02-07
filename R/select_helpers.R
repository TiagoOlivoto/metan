#' @aliases select_helpers
#' @importFrom tidyselect contains
#' @export
tidyselect::contains
#' @importFrom tidyselect ends_with
#' @export
tidyselect::ends_with
#' @importFrom tidyselect everything
#' @export
tidyselect::everything
#' @importFrom tidyselect matches
#' @export
tidyselect::matches
#' @importFrom tidyselect num_range
#' @export
tidyselect::num_range
#' @importFrom tidyselect one_of
#' @export
tidyselect::one_of
#' @importFrom tidyselect starts_with
#' @export
tidyselect::starts_with
#' @importFrom tidyselect last_col
#' @export
tidyselect::last_col
#' @importFrom dplyr select
#' @export
dplyr::select
#' @importFrom dplyr mutate
#' @export
dplyr::mutate
#' @importFrom dplyr group_by
#' @export
dplyr::group_by
#' @importFrom dplyr rename
#' @export
dplyr::rename
#' @importFrom tibble column_to_rownames
#' @export
tibble::column_to_rownames
#' @importFrom tibble rownames_to_column
#' @export
tibble::rownames_to_column
NULL

#' @title Select helper
#' @name Select_helper
#' @param prefix A prefix that start the variable name.
#' @param suffix A suffix that end the variable name.
#' @param n The length of variable names to select. For \code{width_of()} the
#'   selected variables contains \code{n} characters. For
#'   \code{width_greater_than()} and \code{width_less_than()} the selected
#'   variables contains greater and less characteres than \code{n},
#'   respectively.
#' @param vars A character vector of variable names. When called from inside
#'   selecting functions like \code{\link{select_cols}} these are automatically
#'   set to the names of the table.
#' @importFrom dplyr intersect
#' @export
#' @description These functions allow you to select variables based operations
#'   with prefixes and suffixes and length of names.
#' * \code{difference_var()}: Select variables that start with a prefix \strong{AND
#' NOT} end wiht a suffix.
#'
#' * \code{intersect_var()}: Select variables that start with a prefix
#' \strong{AND} end wiht a suffix.
#'
#' * \code{union_var()}: Select variables that start with a prefix \strong{OR}
#' end wiht a suffix.
#' * \code{width_of()}: Select variables with width of \code{n}.
#' * \code{width_greater_than()}: Select variables with width greater than \code{n}.
#' * \code{width_less_than()}: Select variables with width less than \code{n}.
#' * \code{lower_case_only()}: Select variables that contains lower case only
#' (e.g., "env").
#' * \code{upper_case_only()}: Select variables that contains upper case only
#' (e.g., "ENV").
#' * \code{title_case_only()}: Select variables that contains upper case in the first
#' character only (e.g., "Env").
#' @md
#' @examples
#' \donttest{
#' library(metan)
#'
#'
#' # Select variables that start with "C" and not end with "D".
#' data_ge2 %>%
#' select_cols(difference_var("C", "D"))
#'
#' # Select variables that start with "C" and end with "D".
#' data_ge2 %>%
#' select_cols(intersect_var("C", "D"))
#'
#' # Select variables that start with "C" or end with "D".
#' data_ge2 %>%
#' select_cols(union_var("C", "D"))
#'
#' # Select variables with width name of 4
#' data_ge2 %>%
#' select_cols(width_of(4))
#'
#' # Select variables with width name greater than 2
#' data_ge2 %>%
#' select_cols(width_greater_than(2))
#'
#' # Select variables with width name less than 3
#' data_ge2 %>%
#' select_cols(width_less_than(3))
#'
#' # Creating data with messy column names
#' df <- head(data_ge, 3)
#' colnames(df) <- c("Env", "gen", "Rep", "GY", "hm")
#' select_cols(df, lower_case_only())
#' select_cols(df, upper_case_only())
#' select_cols(df, title_case_only())
#' }
#'
difference_var <- function(prefix, suffix) {
  setdiff(starts_with(prefix), ends_with(suffix))
}
#' @name Select_helper
#' @importFrom dplyr union
#' @export
intersect_var <- function(prefix, suffix) {
  intersect(starts_with(prefix), ends_with(suffix))
}
#' @name Select_helper
#' @importFrom dplyr setdiff
#' @export
union_var <- function(prefix, suffix) {
  union(starts_with(prefix), ends_with(suffix))
}
#' @name Select_helper
#' @export
#' @importFrom tidyselect peek_vars
width_of <- function(n, vars = peek_vars(fn = "width_of")) {
  which(nchar(vars) == n)
}
#' @name Select_helper
#' @export
width_greater_than <- function(n, vars = peek_vars(fn = "width_greater_than")) {
  which(nchar(vars) > n)
}
#' @name Select_helper
#' @export
width_less_than <- function(n, vars = peek_vars(fn = "width_less_than")) {
  which(nchar(vars) < n)
}
#' @name Select_helper
#' @export
lower_case_only <- function(vars = peek_vars(fn = "lower_case_only")) {
  vars[!grepl("[A-Z]", vars)]
}
#' @name Select_helper
#' @export
upper_case_only <- function(vars = peek_vars(fn = "upper_case_only")) {
  vars[!grepl("[a-z]", vars)]
}
#' @name Select_helper
#' @export
title_case_only <- function(vars = peek_vars(fn = "title_case_only")) {
  vars[grep("^[A-Z].[a-z]", vars)]
}
