#' Missing value imputation
#'
#' Given a matrix with missing values, impute the missing entries using
#' different algorithms. See \strong{Details} section.
#' @details
#' \strong{\code{EM-AMMI} algorithm}
#'
#' The \code{EM-AMMI} algorithm completes a data set with missing values according to both
#' main and interaction effects. The algorithm works as follows (Gauch and
#' Zobel, 1990):
#'  1. The initial values are calculated as the grand mean increased by main
#'  effects of rows and main effects of columns. That way, the matrix of
#'  observations is pre-filled in.
#'  2. The parameters of the AMMI model are estimated.
#'  3. The adjusted means are calculated based on the AMMI model with
#'  \code{naxis} principal components.
#'  4. The missing cells are filled with the adjusted means.
#'  5. If the maximum distance between the missing value estimations in the two
#'  successive iteration steps ((\code{d0 - d1})) is greater than the assumed
#'  \code{tol}, the steps 2 through 5 are repeated. Declare convergence if
#'  \code{(d0 - d1) < tol}. If \code{max_iter} is achieved without convergence,
#'  the algorithm will stop with a warning.
#'
#' \strong{\code{EM-SVD} algorithm}
#'
#' The \code{EM-SVD} algorithm impute the missing entries using a low-rank Singular
#' Value Decomposition approximation estimated by the Expectation-Maximization
#' algorithm. The algorithm works as follows:
#'  1. Initialize all \code{NA} values to the column means.
#'  2. Compute the first \code{naxis} terms of the SVD of the completed matrix
#'  3. Replace the previously missing values with their approximations from the SVD
#'  4. If the maximum distance between the missing value estimations in the two
#'  successive iteration steps ((\code{d0 - d1})) is greater than the assumed
#'  \code{tol}, the steps 2 through 3 are repeated. Declare convergence if
#'  \code{(d0 - d1) < tol}. If \code{max_iter} is achieved without convergence,
#'  the algorithm will stop with a warning.
#'
#' \strong{\code{colmeans} algorithm}
#' The \code{colmeans} algorithm simply impute the missing entires using the
#' column mean of the respective entire. Thus, there is no iteractive process.
#'
#'
#' @md
#'
#' @param .data A matrix to impute the missing entries. Frequently a two-way
#'   table of genotype means in each environment.
#' @param naxis The rank of the Singular Value Approximation. Defaults to \code{1}.
#' @param algorithm The algorithm to impute missing values. Defaults to
#'   \code{"EM-SVD"}. Other possible values are \code{"EM-AMMI"} and
#'   \code{"colmeans"}. See \strong{Details} section.
#' @param tol The convergence tolerance for the algorithm.
#' @param max_iter The maximum number of steps to take. If \code{max_iter} is
#'   achieved without convergence, the algorithm will stop with a warning.
#' @param simplified Valid argument when \code{algorithm = "EM-AMMI"}. IF
#'   \code{FALSE} (default), the current effects of rows and columns change from
#'   iteration to iteration. If \code{TRUE}, the general mean and effects of
#'   rows and columns are computed in the first iteration only, and in next
#'   iterations uses these values.
#' @param verbose Logical argument. If \code{verbose = FALSE} the code will run
#'   silently.
#'
#' @return An object of class \code{imv} with the following values:
#'
#'  * \strong{.data} The imputed matrix
#'  * \strong{pc_ss} The sum of squares representing variation explained by the
#'  principal components
#'  * \strong{iter} The final number of iterations.
#'  * \strong{final_tol} The maximum change of the estimated values for missing cells in the last step of iteration.
#'  * \strong{final_axis} The final number of principal component axis.
#'  * \strong{convergence} Logical value indicating whether the modern converged.
#' @export
#'
#' @examples
#' \donttest{
#' library(metan)
#' tb <- (1:10) %*% t(1:5) %>%
#'   random_na(prop = 20)
#' tb
#' mod <- impute_miss_val(tb)
#' mod$.data
#'
#' }
impute_miss_val <- function(.data,
                            naxis = 1,
                            algorithm = "EM-SVD",
                            tol = 1e-10,
                            max_iter = 1000,
                            simplified = FALSE,
                            verbose = TRUE){
  if(!has_na(.data)){
    stop("No missing values found in data.")
  }
  if(!algorithm %in% c("EM-AMMI", "EM-SVD", "colmeans")){
    stop("The algorithm must be of of 'EM-AMMI', 'EM-SVD', and 'colmeans")
  }
  if(algorithm == "EM-AMMI"){
    .data <- as.matrix(.data)
    missing <- matrix(1, nrow(.data), ncol(.data))
    missing[is.na(.data)] <- 0
    max_ipc <- min(c(rowSums(missing), colSums(missing))) - 1
    axis_used <- ifelse(max_ipc < naxis, max_ipc, naxis)
    data_in <- .data
    data_mean <- mean(c(.data), na.rm = TRUE)
    row_m <- matrix(rowMeans(.data, na.rm = TRUE), nrow(.data),ncol(.data))
    col_m <- t(matrix(colMeans(.data,na.rm = TRUE), ncol(.data),nrow(.data)))
    estimated <- (-data_mean) + row_m + col_m
    data_in[is.na(data_in)] <- estimated[is.na(.data)]
    iter <- 1
    new_data <- data_in
    change <- tol + 1
    while((change > tol) & (iter < max_iter)){
      if (iter == 1 | !simplified){
        grand_mean <- mean(new_data)
        int_mat <- new_data - grand_mean
        int_mat <- scale(int_mat, scale = FALSE)
        col_center <- attr(int_mat,"scaled:center")
        int_mat <- t(scale(t(int_mat), scale = FALSE))
        row_center <- attr(int_mat,"scaled:center")
      } else {
        int_mat <- new_data - grand_mean - row_eff - col_eff
      }
      if (axis_used >=1){
        SVD <- La.svd(int_mat)
        diag_l <- diag(SVD$d, nrow = axis_used)
        d <- SVD$d[1:axis_used]
        u <- SVD$u[,1:axis_used]
        v <- SVD$v[1:axis_used,]
        adj_val <- u %*% diag_l %*% v
      } else {
        adj_val <- 0
      }
      if(iter == 1 | !simplified){
        row_eff <- matrix(row_center, nrow(.data), ncol(.data))
        col_eff <- t(matrix(col_center,ncol(.data), nrow(.data)))
      }
      new_data2 <- new_data
      new_data2[is.na(.data)] <- (grand_mean + row_eff + col_eff + adj_val)[is.na(.data)]
      change <- max(abs(c((new_data2 - new_data)[is.na(.data)])))
      iter <- iter + 1
      new_data <- new_data2
    }
  }
  if(algorithm == "EM-SVD"){
    data_in <- as.matrix(replace_na(data.frame(.data),  replace = "colmeans"))
    max_ipc <- min(nrow(data_in), ncol(data_in)) - 1
    axis_used <- ifelse(max_ipc < naxis, max_ipc, naxis)
    iter <- 1
    new_data <- data_in
    change <- tol + 1
    while((change > tol) & (iter < max_iter)){
      SVD <- La.svd(new_data)
      diag_l <- diag(SVD$d, nrow = axis_used)
      d <- SVD$d[1:axis_used]
      u <- SVD$u[,1:axis_used]
      v <- SVD$v[1:axis_used,]
      adj_val <- u %*% diag_l %*% v
      new_data2 <- new_data
      new_data2[is.na(.data)] <- adj_val[is.na(.data)]
      change <- max(abs(c((new_data2 - new_data)[is.na(.data)])))
      iter <- iter + 1
      new_data <- new_data2
    }
  }
  if(algorithm == "colmeans"){
    new_data <- as.matrix(replace_na(data.frame(.data),  replace = "colmeans"))
    pc_ss <- NULL
    iter <- NULL
    final_tol <- NULL
    final_naxis <- NULL
    converged <- NULL
    change <- NULL
    axis_used <- NULL
  }
  if(algorithm %in% c("EM-AMMI", 'EM-SVD')){
    pc_ss <- d^2
    converged <- ifelse(change <= tol, TRUE, FALSE)
    if (axis_used < naxis){
      converged = FALSE
    }
    if (axis_used == 0){
      SVD <- list(d = 0)
    }
    if(verbose == TRUE){
      cat("----------------------------------------------\n")
      cat("Convergence information\n")
      cat("----------------------------------------------\n")
      cat("Number of iterations:", iter)
      cat("\nFinal tolerance:", change)
      cat("\nNumber of axis:", axis_used)
      cat("\nConvergence:", converged)
      cat("\n----------------------------------------------\n")
    }
    if(!converged){
      warning("Maximum number of iterations achieved without convergence.\nTo increase the number of iterations use 'max_iter'.", call. = FALSE)
    }
  }
  return(list(
    .data = as.matrix(new_data),
    pc_ss = pc_ss,
    iter = iter,
    final_tol = change,
    final_naxis = axis_used,
    convergence = converged) %>% set_class("imv"))
}
