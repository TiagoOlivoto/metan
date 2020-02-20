#' Select a set of predictors with minimal multicollinearity
#'
#' Select a set of predictors with minimal multicollinearity using the variance
#' inflation factor (VIF) as criteria to remove collinear variables. The
#' algorithm will: \strong{(i)} compute the VIF value of the correlation matrix
#' containing the variables selected in \code{...}; \strong{(ii)} arrange the
#' VIF values and delete the variable with the highest VIF; and \strong{(iii)}
#' iterate step \strong{ii} until VIF value is less than or equal to
#' \code{max_vif}.
#'
#' @param .data The data set containing the variables.
#' @param ... Variables to be submitted to selection. If \code{...} is null then
#'   all the numeric variables from \code{.data} are used. It must be a single
#'   variable name or a comma-separated list of unquoted variables names.
#' @param max_vif The maximum value for the Variance Inflation Factor
#'   (threshold) that will be accepted in the set of selected predictors.
#' @param missingval How to deal with missing values. For more information,
#'   please see \code{\link[stats]{cor}()}.
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
  if(!missing(...)){
    xxx <-  select_cols(.data, ...)
    has_text_in_num(xxx)
    xxx %<>% select_numeric_cols()
  } else{
    xxx <- select_numeric_cols(.data)
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
