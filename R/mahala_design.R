#' Mahalanobis distance from designed experiments
#' @description
#' `r badge('stable')`
#'
#' Compute the Mahalanobis distance using data from an experiment conducted in a
#' randomized complete block design or completely randomized design.
#'
#'
#' @param .data The dataset containing the columns related to Genotypes,
#'   replication/block and response variables, possible with grouped data passed
#'   from [dplyr::group_by()].
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#'   replications/blocks.
#' @param resp The response variables. For example `resp = c(var1, var2,
#'   var3)`.
#' @param design The experimental design. Must be RCBD or CRD.
#'@param by One variable (factor) to compute the function by. It is a shortcut
#'  to [dplyr::group_by()]. To compute the statistics by more than
#'  one grouping variable use that function.
#' @param return What the function return? Default is 'distance', i.e., the
#'   Mahalanobis distance. Alternatively, it is possible to return the matrix of
#'   means `return = 'means'`, or the variance-covariance matrix of
#'   residuals `return = 'covmat'`.
#' @return A symmetric matrix with the Mahalanobis' distance. If `.data` is
#'   a grouped data passed from [dplyr::group_by()] then the results
#'   will be returned into a list-column of data frames.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' maha <- mahala_design(data_g,
#'                       gen = GEN,
#'                       rep = REP,
#'                       resp = everything(),
#'                       return = "covmat")
#'
#' # Compute one distance for each environment (all numeric variables)
#'maha_group <- mahala_design(data_ge,
#'                            gen = GEN,
#'                            rep = REP,
#'                            resp = everything(),
#'                            by = ENV)
#'
#' # Return the variance-covariance matrix of residuals
#'cov_mat <- mahala_design(data_ge,
#'                         gen = GEN,
#'                         rep = REP,
#'                         resp = c(GY, HM),
#'                         return = 'covmat')
#'}
mahala_design <- function(.data,
                          gen,
                          rep,
                          resp,
                          design = "RCBD",
                          by = NULL,
                          return = "distance") {
  if (!design %in% c("RCBD", "CRD")) {
    stop("The experimental design must be RCBD or CRD.")
  }
  if (!missing(by)){
    if(length(as.list(substitute(by))[-1L]) != 0){
      stop("Only one grouping variable can be used in the argument 'by'.\nUse 'group_by()' to pass '.data' grouped by more than one variable.", call. = FALSE)
    }
    .data <- group_by(.data, {{by}})
  }
  if(is_grouped_df(.data)){
    results <-
      .data %>%
      doo(mahala_design,
          gen = {{gen}},
          rep = {{rep}},
          resp = {{resp}},
          design = design,
          return = return)
    return(add_class(results, "mahala_group"))
  }
  factors <- select(.data,
                    GEN = {{gen}},
                    REP = {{rep}}) %>%
    as_factor(1:2)
  GEN <- factors$GEN
  REP <- factors$REP
    vars <- .data %>%
      select({{resp}}, -{{gen}}, -{{rep}}) %>%
      select_numeric_cols()
    nvar <- ncol(vars)
    mat <- matrix(nrow = nvar, ncol = nvar)
    covdata <- data.frame(matrix(nrow = nrow(.data), ncol = nvar))
    vin <- 0
    for (var in 1:nvar) {
      vin <- vin + 1
      Y <- vars[[var]]
      covdata[, vin] <- Y
      if (design == "RCBD") {
        model <- anova(aov(Y ~ GEN + REP))
        diag(mat)[vin] <- model[3, 3]
      } else {
        model <- anova(aov(Y ~ GEN))
        diag(mat)[vin] <- model[2, 3]
      }
      colnames(covdata)[[vin]] <- paste(names(vars[var]))
    }
    means <-
      data.frame(cbind(GEN, covdata)) %>%
      mean_by(GEN) %>%
      column_to_rownames("GEN")
    covdata2 <- comb_vars(data.frame(covdata), order = "first")
    index <- data.frame(t(combn(ncol(mat), 2)))
    temp <- NULL
    for (i in 1:ncol(covdata2)) {
      if (design == "RCBD") {
        model <- anova(aov(covdata2[[i]] ~ GEN + REP))
        temp[i] <- (model[3, 3] - diag(mat)[[index[i, 1]]] - diag(mat)[[index[i, 2]]])/2
      } else {
        model <- anova(aov(covdata2[[i]] ~ GEN))
        temp[i] <- (model[2, 3] - diag(mat)[[index[i,
                                                   1]]] - diag(mat)[[index[i, 2]]])/2
      }
    }
    mat[lower.tri(mat, diag = FALSE)] <- temp
    rownames(mat) <- colnames(mat) <- colnames(means)
    dist <- mahala(.means = means, covar = make_sym(mat, diag = diag(mat)), inverted = FALSE)
    if (return == "distance") {
      return(as.data.frame(dist))
    }
    if (return == "covmat") {
      return(as.data.frame(make_sym(mat, diag = diag(mat))))
    }
    if (return == "means") {
      return(as_tibble(means))
    }
}
