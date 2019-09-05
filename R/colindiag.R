#' Collinearity Diagnostics
#'
#' Perform a (multi)collinearity diagnostic of a correlation matrix of
#' predictor variables based several indicators, as shown by Olivoto et al.
#' (2017).
#'
#'
#' @param .data The data to be analyzed. Must be a symmetric correlation matrix
#' or, a dataframe containing the predictor variables, or an object of class
#' \code{split_factors}.
#' @param n If a correlation matrix is provided, then \code{n} is the number of
#' objects used to compute the correlation coefficients.
#' @param verbose If \code{verbose = TRUE} then some results are shown in the
#' console.
#' @return
#'
#' The following values are returned. Please, note that if a grouping variable
#' is used, then the results are returned into a list.
#'
#' \item{cormat}{A symmetric Pearson's coefficient correlation matrix between
#' the variables}
#'
#' \item{corlist}{A hypothesis testing for each of the correlation
#' coefficients}
#'
#' \item{evalevet}{The eigenvalues with associated eigenvectors of the
#' correlation matrix}
#'
#' \item{VIF}{The Variance Inflation Factors, being the diagonal elements of
#' the inverse of the correlation matrix.}
#'
#' \item{CN}{The Condition Number of the correlation matrix, given by the ratio
#' between the largest and smallest eigenvalue.}
#'
#' \item{det}{The determinant of the correlation matrix.}
#'
#' \item{largest_corr}{The largest correlation (in absolute value) observed.}
#'
#' \item{smallest_corr}{The smallest correlation (in absolute value) observed.}
#'
#' \item{weight_var}{The variables with largest eigenvector (largest weight) in
#' the eigenvalue of smallest value, sorted in decreasing order.}
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Olivoto, T., V.Q. Souza, M. Nardino, I.R. Carvalho, M. Ferrari,
#' A.J. Pelegrin, V.J. Szareski, and D. Schmidt. 2017. Multicollinearity in
#' path analysis: a simple method to reduce its effects. Agron. J. 109:131-142.
#' doi:10.2134/agronj2016.04.0196.
#' \href{https://dl.sciencesocieties.org/publications/aj/abstracts/109/1/131}{doi:10.2134/agronj2016.04.0196}
#' @export
#' @examples
#'
#' # Using the correlation matrix
#' library(metan)
#' cor_iris = cor(iris[,1:4])
#' n = nrow(iris)
#' colindiag(cor_iris, n = n, verbose = FALSE)
#'
#' # Using the pipe operator %>%
#' cor_iris %>% colindiag(n = n, verbose = FALSE)
#'
#' # Diagnostic by species, storing into an object
#' col_diag_spec = iris %>%
#'                 split_factors(Species) %>%
#'                 colindiag()
#'
#'
colindiag <- function(.data, n = NULL, verbose = TRUE) {
  if (!class(.data) %in% c("matrix", "data.frame", "split_factors",
                           "covcor_design")) {
    stop("The object 'x' must be a correlation matrix, a data.frame or an object of class split_factors")
  }
  if (is.matrix(.data) && is.null(n)) {
    stop("You have a matrix but the sample size used to compute the correlations (n) was not declared.")
  }
  if (class(.data) == "covcor_design" && is.null(n)) {
    stop("You have a list of matrices but the sample size used to compute the correlations (n) was not declared.")
  }
  if (is.data.frame(.data) && !is.null(n)) {
    stop("You cannot informe the sample size because a data frame was used as input.")
  }
  if (is.split_factors(.data) && !is.null(n)) {
    stop("You cannot informe the sample size because a data frame was used as input.")
  }
  internal <- function(x) {
    if (is.matrix(x)) {
      cor.x <- x
    }
    if (is.data.frame(x)) {
      cor.x <- cor(x)
      n <- nrow(x)
    }
    eigen <- eigen(cor.x)
    Det <- det(cor.x)
    NC <- max(eigen$values)/min(eigen$values)
    Aval <- data.frame(eigen$values)
    names(Aval) <- "Eigenvalues"
    Avet <- data.frame(t(eigen$vectors))
    names(Avet) <- colnames(x)
    AvAvet <- cbind(Aval, Avet)
    VIF <- data.frame(diag(solve_svd(cor.x)))
    names(VIF) <- "VIF"
    results <- data.frame(linear = as.vector(t(cor.x)[lower.tri(cor.x,
                                                                diag = F)]))
    results <- dplyr::mutate(results, t = linear * (sqrt(n -
                                                           2))/(sqrt(1 - linear^2)), prob = 2 * (1 - pt(abs(t),
                                                                                                        df = n - 2)))
    names <- colnames(x)
    combnam <- combn(names, 2, paste, collapse = " x ")
    rownames(results) <- names(sapply(combnam, names))
    largest_corr <- paste0(rownames(results)[which.max(abs(results$linear))],
                           " = ", round(results[which.max(abs(results$linear)),
                                                1], 3))
    smallest_corr <- paste0(rownames(results)[which.min(abs(results$linear))],
                            " = ", round(results[which.min(abs(results$linear)),
                                                 1], 3))
    ncorhigh <- sum(results$linear >= abs(0.8))
    if (verbose == TRUE) {
      if (NC > 1000) {
        cat(paste0("Severe multicollinearity in the matrix! Pay attention on the variables listed bellow\n",
                   "CN = ", round(NC, 3), "\n"))
      }
      if (NC < 100) {
        cat(paste0("Weak multicollinearity in the matrix\n",
                   "NC = ", round(NC, 3), "\n"))
      }
      if (NC > 100 & NC < 1000) {
        cat(paste0("The multicollinearity in the matrix should be investigated.\n",
                   "NC = ", round(NC, 3), "\n", "Largest VIF = ",
                   max(VIF), "\n"))
      }
    }
    ultimo <- data.frame(Peso = t(AvAvet[c(nrow(AvAvet)),
                                         ])[-c(1), ])
    abs <- data.frame(Peso = abs(ultimo[, "Peso"]))
    rownames(abs) <- rownames(ultimo)
    ultimo <- abs[order(abs[, "Peso"], decreasing = T), ,
                  drop = FALSE]
    pesovarname <- paste(rownames(ultimo), collapse = " > ")
    final <- list(cormat = data.frame(cor.x), corlist = results,
                  evalevet = AvAvet, VIF = VIF, CN = NC, det = Det,
                  largest_corr = largest_corr, smallest_corr = smallest_corr,
                  weight_var = pesovarname)
    if (verbose == TRUE) {
      cat(paste0("Matrix determinant: ", round(Det, 7)),
          "\n")
      cat(paste0("Largest correlation: ", largest_corr),
          "\n")
      cat(paste0("Smallest correlation: ", smallest_corr),
          "\n")
      cat(paste0("Number of VIFs > 10: ", length(VIF >
                                                   10)), "\n")
      cat(paste0("Number of correlations with r >= |0.8|: ",
                 ncorhigh), "\n")
      cat(paste0("Variables with largest weight in the last eigenvalues: ",
                 "\n", pesovarname), "\n")
    }
    return(invisible(final))
  }
  if (is.matrix(.data)) {
    out <- internal(.data)
  }
  if (any(class(.data) %in% c("split_factors", "covcor_design"))) {
    out <- lapply(seq_along(.data), function(x) {
      if (verbose == TRUE) {
        cat("\n----------------------------------------------------------------------------\n")
        cat("Level:", names(.data)[[x]], "\n")
        cat("----------------------------------------------------------------------------\n")
      }
      internal(.data[[x]])
    })
    names(out) <- names(.data)
  }
  if (is.data.frame(.data)) {
    data <- .data[, unlist(lapply(.data, is.numeric))]
    out <- internal(data)
    if (verbose == TRUE) {
      if (sum(lapply(.data, is.factor) == TRUE) > 0) {
        message("The factors ", paste0(collapse = " ",
                                       names(.data[, unlist(lapply(.data, is.factor))])),
                " where excluded to perform the analysis. If you want to perform an analysis for each level of a factor, use the function 'split_factors() before.' ")
      }
    }
  }
  invisible(out)
}
