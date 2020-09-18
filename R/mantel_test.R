#' Mantel test
#'
#' Performs a Mantel test between two correlation/distance matrices. The
#' function calculates the correlation between two matrices, the Z-score that is
#' is the sum of the products of the corresponding elements of the matrices and
#' a two-tailed p-value (null hypothesis: \code{r = 0}).
#'
#' @param mat1,mat2 A correlation matrix or an object of class \code{dist}.
#' @param nboot The number of permutations to be used. Defaults to \code{1000}.
#' @param plot if \code{plot = TRUE}, plots the density estimate of the
#'   permutation distribution along with the observed Z-score as a vertical
#'   line.
#' @return
#' * \code{mantel_r} The correlation between the two matrices.
#' * \code{z_score} The Z-score.
#' * \code{p-value} The quantile of the observed Z-score. in the permutation
#' distribution.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' # Test if the correlation of traits (data_ge2 dataset)
#' # changes between A1 and A2 levels of factor ENV
#' A1 <- corr_coef(data_ge2 %>% subset(ENV == "A1"))[["cor"]]
#' A2 <- corr_coef(data_ge2 %>% subset(ENV == "A2"))[["cor"]]
#' mantel_test(A1, A2, plot = TRUE)
#'
#'}
mantel_test <- function(mat1, mat2, nboot = 1000, plot = FALSE){
  mat1 <- as.matrix(mat1)
  mat2 <- as.matrix(mat2)
  if(!identical(dim(mat1), dim(mat2))){
    stop("Matrices with different dimension are not allowed.", call. = FALSE)
  }
  diag(mat1) <- diag(mat2) <- 0
  permute_mat <- function(mat1, n){
    s <- sample(1:n)
    mat1[s, s]
  }
  stat_z <- function(mat1, mat2) {
    sum(mat1 * mat2) / 2
  }
  n <- nrow(mat1)
  real <- stat_z(mat1, mat2)
  null_st <- replicate(nboot, stat_z(mat1, permute_mat(mat2, n)))
  pval <- (2 * min(sum(null_st >= real), sum(null_st <= real)) + 1) / (nboot + 1)
  pval <- ifelse(pval > 1, 1, pval)
  mantel_r <- cor(as.vector(mat1[lower.tri(mat1)]), as.vector(mat2[lower.tri(mat2)]))
  if(plot == TRUE){
    p1 <-
    ggplot(data.frame(null_st), aes(x = null_st)) +
      geom_density() +
      geom_vline(xintercept = real,
                 linetype = "dashed") +
      theme_metan(color.background = "white") +
      labs(x = NULL)
    plot(p1)
  }
  return(list(mantel_r = mantel_r, z_score = real, p_value = pval))
}
