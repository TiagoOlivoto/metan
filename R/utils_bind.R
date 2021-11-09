#' Helper function for binding rows
#' @description
#' * [rbind_fill_id()] `r badge('stable')` Implements the common pattern of `do.call(rbind, dfs)`
#' with data frame identifier and filling of missing values.
#'
#' @name utils_bind
#' @param ... The dataframes. Either a list of data frames, or a comma-separated
#'   list of dataframes.
#' @param .fill When row-binding, columns are matched by name, and any missing
#'   columns will be filled with `NA` Defaults to `NA`.
#' @param .id Data frame identifier. If a comma-separated list of data frames is
#'   supplied, the labels are taken from the names of the objects. When a list
#'   of data frames is supplied, the labels are taken from the names of the
#'   list. If no names are found, a numeric sequence is used instead.
#' @return A data frame.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#' \donttest{
#' (df1 <- data.frame(v1 = c(1, 2), v2 = c(2, 3)))
#' (df2 <- data.frame(v3 = c(4, 5)))
#' rbind_fill_id(df1, df2)
#' rbind_fill_id(df1, df2,
#'              .fill = ".",
#'              .id = "dfs")
#'
#' # Named list
#' list <- list(a = df1, b = df2)
#' rbind_fill_id(list, .id = "dfs")
#'
#' # Unnamed list
#' list <- list(df1, df2)
#' rbind_fill_id(list, .id = "dfs")
#'
#' }
#'
rbind_fill_id <- function(..., .id = NULL, .fill = NA){
    set_name <-
        as.character(
            sapply(quos(...), function(x){
                rlang::quo_get_expr(x)
            })
        )
    dfs <- list(...)
    if(length(dfs) == 1){
        dfs <- reduce(dfs, unnest) %>% set_class("list")
    }
    if(is.null(names(dfs))){
        names(dfs) <- set_name
    }
    bind <-
        bind_rows(dfs, .id = .id) %>%
        replace_na(everything(), replacement = .fill)
    return(bind)
}
