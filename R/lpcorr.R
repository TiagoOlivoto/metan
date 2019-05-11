#' Linear and Partial Correlation Coefficients
#'
#' Estimates the linear and partial correlation coefficients using as input a
#' data frame or a correlation matrix.
#'
#'
#' @param .data The data to be analyzed. Must be a symmetric correlation matrix
#' or, a dataframe containing the predictor variables, or an object of class
#' \code{split_factors}.
#' @param n If a correlation matrix is provided, then \code{n} is the number of
#' objects used to compute the correlation coefficients.
#' @param verbose If \code{verbose = TRUE} then some results are shown in the
#' console.
#' @return If a grouping factor is used then a list is returned with the
#' following values. \item{linear.mat}{The matrix of linear correlation.}
#' \item{partial.mat}{The matrix of partial correlations.}
#' \item{results}{Hypotesis testing for each pairwise comparation.}
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'
#' library(METAAB)
#' library(dplyr)
#' partial1 = lpcor(iris)
#'
#' # Alternatively using the pipe operator %>%
#' partial2 = iris %>% lpcor()
#'
#' partial3 = cor(iris[1:4]) %>%
#'            lpcor(n = nrow(iris),
#'                  verbose = FALSE)
#'
#' partial4 = iris %>%
#'            split_factors(Species) %>%
#'            lpcor()
#'
lpcor <- function(.data, n = NULL, verbose = TRUE) {

  if(!is.matrix(.data) && !is.data.frame(.data) && !is.split_factors(.data)){
    stop("The object 'x' must be a correlation matrix, a data.frame or an object of class split_factors")
  }
  if(is.matrix(.data) && is.null(n)){
    stop("You have a matrix but the sample size used to compute the correlations (n) was not declared.")
  }
  if(is.data.frame(.data) && !is.null(n)){
    stop("You cannot informe the sample size because a data frame was used as input.")
  }

  if(is.split_factors(.data) && !is.null(n)){
    stop("You cannot informe the sample size because a data frame was used as input.")
  }

    internal = function(x){
      if(is.matrix(x)){
        cor.x = x
      }
      if(is.data.frame(x)){
        cor.x = cor(x)
        n = nrow(x)
      }
      nvar <- ncol(cor.x)
      df <- n - nvar
      if(df < 0){
        warning("The number of variables is higher than the number of individuals. Hypothesis testing will not be made.",
                call. = FALSE)
      }
    m <- as.matrix(cor.x)
    X.resid <- -(solve(m))
    diag(X.resid) <- 1/(1 - (1 - 1/diag(solve(m))))
    X.resid <- cov2cor(X.resid)
    results <- data.frame(linear = as.vector(t(m)[lower.tri(m, diag = F)]))
    results <- suppressWarnings(dplyr::mutate(results,
                             partial = as.vector(t(X.resid)[lower.tri(X.resid, diag = F)]),
                             t = partial/(sqrt(1 - partial^2)) * sqrt(n - nvar), prob = 2 *
                               (1 - pt(abs(t), df = df))))
    names <- colnames(x)
    combnam <- combn(names, 2, paste, collapse = " x ")
    rownames(results) <- names(sapply(combnam, names))
    if(verbose == TRUE){
      print.data.frame(results)
    }
    return(list(linear.mat = m,
                partial.mat = as.matrix(X.resid),
                results = results))
    }

    if (is.matrix(.data)) {
       out = internal(.data)
    }

    if (any(class(.data) == "split_factors")) {
      out = lapply(seq_along(.data), function(x){
        if(verbose == TRUE){
          cat("\n----------------------------------------------------------------------------\n")
          cat("Level:", names(.data)[[x]], "\n")
          cat("----------------------------------------------------------------------------\n")
        }
        internal(.data[[x]])
      })
      names(out) = names(.data)
    }

    if (is.data.frame(.data)) {
        if(sum(lapply(.data, is.factor)==TRUE)>0){
        }
        data = .data[ , unlist(lapply(.data, is.numeric))]
        out = internal(data)
        if(verbose == TRUE){
          message("The factors ", paste0(collapse = " ", names(.data[ , unlist(lapply(.data, is.factor)) ])),
                  " where excluded to perform the analysis. If you want to perform an analysis for each level of a factor, use the function 'split_factors() before.' ")
        }
    }

    if (class(.data) == "split_factors") {
    invisible(structure(out, class = "lpcor_group"))
    } else {
    invisible(structure(out, class = "lpcor"))
    }
}
