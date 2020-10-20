#' Generate a covariance matrix
#'
#' Given the variances and desired correlations, generate a covariance matrix
#'
#' @param cormat A symmetric matrix with desired correlations.
#' @param var A numeric vector with variances. It must have length equal to the
#'   number of elements in the diagonal of \code{cormat}.
#'
#' @return A (co)variance matrix
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#'
#' @examples
#' \donttest{
#' cormat <-
#' matrix(c(1,  0.9, -0.4,
#'          0.9,  1,  0.6,
#'         -0.4, 0.6, 1),
#'       nrow = 3,
#'       ncol = 3)
#' get_covmat(cormat, var =  c(16, 25, 9))
#' }
#'
get_covmat <- function(cormat, var){
  vars <- diag(var)
  covariances <- matrix(NA, ncol = ncol(cormat), nrow = nrow(cormat))
  diag(covariances) <- var
  for (i in 1:nrow(vars)) {
    for (j in 1:ncol(vars)){
      if(i == j){
        next
      } else{
        covariances[i, j] <- cormat[i, j] * sqrt(vars[i, i] * vars[j, j])
      }
    }
  }
  rownames(covariances) <- colnames(covariances) <- paste("V", 1:ncol(covariances), sep = "")
  return(covariances)
}
