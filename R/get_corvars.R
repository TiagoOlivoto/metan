#' Generate normal, correlated variables
#'
#' Given the mean and desired correlations, generate normal, correlated
#' variables.
#'
#' @param n The number of samples required.
#' @param mu A vector with the means for the variables.
#' @param sigma A symmetric, positive-definite matrix with the (co)variance or
#'   correlation matrix of the variables.
#' @param tol Tolerance (relative to largest variance) for numerical lack of
#'   positive-definiteness in sigma.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @return
#' A tibble containing the simulated data.
#' @export
#'
#' @examples
#' \donttest{
#' sigma <- matrix(c(1,  .3,  0,
#'                   .3,   1, .9,
#'                   0,   .9,  1),3,3)
#' mu <- c(6,50,5)
#'
#' set.seed(04182019)
#' df <- get_corvars(n = 10000, mu = mu, sigma = sigma)
#'means_by(df)
#' cor(df)
#' }
get_corvars <- function (n = 10, mu, sigma, tol = 1e-06) {
  # Adapted from MASS::mvrnorm()
  p <- length(mu)
  if (!all(dim(sigma) == c(p, p))) {
    stop("incompatible arguments")
  }
  eS <- eigen(sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) {
    stop("'sigma' is not positive definite")
  }
  X <- matrix(rnorm(p * n), n)
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(sigma))) {
    nm <- dn[[1L]]
  }
  dimnames(X) <- list(nm, NULL)
  if (n == 1) {
    res <- drop(X)
  } else{
    res <- t(X)
  }
  return(as_tibble(data.frame(res)))
}
