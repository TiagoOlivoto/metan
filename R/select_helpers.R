#' @aliases select_helpers
#' @importFrom dplyr contains
#' @export
dplyr::contains
#' @importFrom dplyr ends_with
#' @export
dplyr::ends_with
#' @importFrom dplyr everything
#' @export
dplyr::everything
#' @importFrom dplyr matches
#' @export
dplyr::matches
#' @importFrom dplyr num_range
#' @export
dplyr::num_range
#' @importFrom dplyr one_of
#' @export
dplyr::one_of
#' @importFrom dplyr starts_with
#' @export
dplyr::starts_with
#' @importFrom dplyr last_col
#' @export
dplyr::last_col
#' @importFrom dplyr select
#' @export
dplyr::select
NULL

#' @title Select helper
#' @name Select_helper
#' @param prefix A prefix that starts the variable name.
#' @param suffix A suffix that ends the variable name.
#' @importFrom dplyr intersect
#' @export
#' @description These functions allow you to select variables based prefixes and suffixes.
#' * \code{intersect_var()}: Select variables that starts with a prefix
#' \strong{AND} ends wiht a suffix.
#' * \code{union_var()}: Select variables that starts with a prefix \strong{OR}
#' ends wiht a suffix.
#' * \code{difference_var()}: Select variables that starts with a prefix \strong{AND
#' NOT} ends wiht a suffix.
#' @md
#' @examples
#' library(metan)
#' # Select variables that starts with "C" and ends with "D".
#' data_ge2 %>%
#' select(intersect_var("C", "D"))
#'
#' # Select variables that starts with "D" or ends with "D".
#' data_ge2 %>%
#' select(union_var("C", "D"))
#'
#' # Select variables that starts with "D" and not ends with "D".
#' data_ge2 %>%
#' select(difference_var("C", "D"))
#'
intersect_var <- function(prefix, suffix) {
  intersect(starts_with(prefix), ends_with(suffix))
}
NULL

#' @title Select helper
#' @name Select_helper
#' @keywords internal
#' @importFrom dplyr union
#' @export
union_var <- function(prefix, suffix) {
  union(starts_with(prefix), ends_with(suffix))
}

#' @title Select helper
#' @name Select_helper
#' @keywords internal
#' @importFrom dplyr setdiff
#' @export
difference_var <- function(prefix, suffix) {
  setdiff(starts_with(prefix), ends_with(suffix))
}
