#' Select a set of predictors with minimal multicollinearity
#' @description
#' `r badge('stable')`
#'
#' Select a set of predictors with minimal multicollinearity using the variance
#' inflation factor (VIF) as criteria to remove collinear variables. The
#' algorithm will: **(i)** compute the VIF value of the correlation matrix
#' containing the variables selected in `...`; **(ii)** arrange the
#' VIF values and delete the variable with the highest VIF; and **(iii)**
#' iterate step **ii** until VIF value is less than or equal to
#' `max_vif`.
#'
#' @param .data The data set containing the variables.
#' @param ... Variables to be submitted to selection. If `...` is null then
#'   all the numeric variables from `.data` are used. It must be a single
#'   variable name or a comma-separated list of unquoted variables names.
#' @param max_vif The maximum value for the Variance Inflation Factor
#'   (threshold) that will be accepted in the set of selected predictors.
#' @param missingval How to deal with missing values. For more information,
#'   please see [stats::cor()].
#'
#' @return A data frame showing the number of selected predictors, maximum VIF
#'   value, condition number, determinant value, selected predictors and removed
#'   predictors from the original set of variables.
#' @export
#'
#' @examples
#' \donttest{
#' library(metan)
#' # All numeric variables
#' non_collinear_vars(data_ge2)
#'
#' # Select variables and choose a VIF threshold to 5
#' non_collinear_vars(data_ge2, EH, CL, CW, KW, NKE, max_vif = 5)
#' }
non_collinear_vars <- function(.data,
                               ...,
                               max_vif = 10,
                               missingval = "pairwise.complete.obs"){
  if(is_grouped_df(.data)){
    result <-
      .data |>
      doo(non_collinear_vars,
          ...,
          max_vif = max_vif,
          missingval = missingval)
    return(result)
  }
  if(!missing(...)){
    xxx <-  select_cols(.data, ...)
    xxx %<>% select_numeric_cols()
    if(has_na(xxx)){
      xxx <- remove_rows_na(xxx)
      has_text_in_num(xxx)
    }
  } else{
    xxx <- select_numeric_cols(.data)
    if(has_na(xxx)){
      xxx <- remove_rows_na(xxx)
      has_text_in_num(xxx)
    }
  }
  varin <- xxx
  cor.xx <- cor(xxx, use = missingval)
  VIF <- data.frame(VIF = diag(solve_svd(cor.xx))) %>%
    rownames_to_column("VAR") %>%
    arrange(VIF)
  repeat {
    VIF2 <- slice(VIF, -n())
    xxx2 <- .data[VIF2$VAR]
    VIF3 <- data.frame(VIF = diag(solve_svd(cor(xxx2, use = missingval))))%>%
      rownames_to_column("VAR") %>%
      arrange(VIF)
    if (max(VIF3$VIF) <= max_vif)
      break
    VIF <- VIF3
  }
  xxx <- .data[VIF3$VAR] %>% as.data.frame()
  selectedpred <- VIF3$VAR
  deleted <- varin %>%
    select(-selectedpred)
  eval <- eigen(cor(xxx2))$value
  cn <- max(eval) / min(eval)
  dtm <- det(cor(xxx2))
  df <- data.frame(Parameter = c("Predictors",
                                 "VIF",
                                 "Condition Number",
                                 "Determinant",
                                 "Selected",
                                 "Removed"),
                   values = c(nrow(VIF3),
                              round(max(VIF3$VIF), 3),
                              round(cn, 3),
                              round(dtm, 10),
                              paste0(selectedpred, collapse = ", "),
                              paste0(names(deleted), collapse = ", ")))
  return(df)
}
