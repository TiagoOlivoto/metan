#' Collinearity Diagnostics
#'
#' Perform a (multi)collinearity diagnostic of a correlation matrix of predictor
#' variables using several indicators, as shown by Olivoto et al. (2017).
#'
#'
#' @param .data The data to be analyzed. It must be a symmetric correlation
#'   matrix, or a data frame, possible with grouped data passed from
#'   \code{\link[dplyr]{group_by}()}.
#' @param ... Variables to use in the correlation. If \code{...} is null then
#'   all the numeric variables from \code{.data} are used. It must be a single
#'   variable name or a comma-separated list of unquoted variables names.
#'@param by One variable (factor) to compute the function by. It is a shortcut
#'  to \code{\link[dplyr]{group_by}()}. To compute the statistics by more than
#'  one grouping variable use that function.
#' @param n If a correlation matrix is provided, then \code{n} is the number of
#'   objects used to compute the correlation coefficients.
#' @return
#'
#' If \code{.data} is a grouped data passed from \code{\link[dplyr]{group_by}()}
#' then the results will be returned into a list-column of data frames.
#'
#' * \strong{cormat} A symmetric Pearson's coefficient correlation matrix
#' between the variables
#'
#' * \strong{corlist} A hypothesis testing for each of the correlation
#' coefficients
#'
#' * \strong{evalevet} The eigenvalues with associated eigenvectors of the
#' correlation matrix
#'
#' * \strong{VIF} The Variance Inflation Factors, being the diagonal elements of
#' the inverse of the correlation matrix.
#'
#' * \strong{CN} The Condition Number of the correlation matrix, given by the
#' ratio between the largest and smallest eigenvalue.
#'
#' * \strong{det} The determinant of the correlation matrix.
#'
#' * \strong{ncorhigh} Number of correlation greather than |0.8|.
#'
#' * \strong{largest_corr} The largest correlation (in absolute value) observed.
#'
#' * \strong{smallest_corr} The smallest correlation (in absolute value)
#' observed.
#'
#' * \strong{weight_var} The variables with largest eigenvector (largest weight)
#' in the eigenvalue of smallest value, sorted in decreasing order.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Olivoto, T., V.Q. Souza, M. Nardino, I.R. Carvalho, M. Ferrari,
#'   A.J. Pelegrin, V.J. Szareski, and D. Schmidt. 2017. Multicollinearity in
#'   path analysis: a simple method to reduce its effects. Agron. J.
#'   109:131-142. doi:10.2134/agronj2016.04.0196.
#'   \href{https://dl.sciencesocieties.org/publications/aj/abstracts/109/1/131}{doi:10.2134/agronj2016.04.0196}
#'
#' @references Olivoto, T., M. Nardino, I.I.R. Carvalho, D.N. Follmann, M.
#'   Ferrari, A.J. de Pelegrin, V.J. Szareski, A.C. de Oliveira, B.O. Caron, and
#'   V.Q. de Souza. 2017. Optimal sample size and data arrangement method in
#'   estimating correlation matrices with lesser collinearity: A statistical
#'   focus in maize breeding. African J. Agric. Res. 12:93-103.
#'   \href{https://academicjournals.org/journal/AJAR/article-abstract/5DB554B62359}{doi:10.5897/AJAR2016.11799}.
#'
#' @export
#' @examples
#'\donttest{
#' # Using the correlation matrix
#' library(metan)
#'
#' cor_iris <- cor(iris[,1:4])
#' n <- nrow(iris)
#'
#' col_diag <- colindiag(cor_iris, n = n)
#'
#'
#' # Using a data frame
#' col_diag_gen <- data_ge2 %>%
#'                 group_by(GEN) %>%
#'                 colindiag()
#'
#' # Diagnostic by levels of a factor
#' # For variables with "N" in variable name
#' col_diag_gen <- data_ge2 %>%
#'                 group_by(GEN) %>%
#'                 colindiag(contains("N"))
#'}
colindiag <- function(.data,
                      ...,
                      by = NULL,
                      n = NULL) {

  if (!has_class(.data, c("matrix", "data.frame", "grouped_df", "covcor_design", "tbl_df"))) {
    stop("The object 'x' must be a correlation matrix, a data.frame or a grouped data.frame")
  }
  if (is.matrix(.data) && is.null(n)) {
    stop("You have a matrix but the sample size used to compute the correlations (n) was not declared.")
  }
  if (is.data.frame(.data) && !is.null(n)) {
    stop("You cannot informe the sample size because a data frame was used as input.")
  }
  if (!missing(by)){
    if(length(as.list(substitute(by))[-1L]) != 0){
      stop("Only one grouping variable can be used in the argument 'by'.\nUse 'group_by()' to pass '.data' grouped by more than one variable.", call. = FALSE)
    }
    .data <- group_by(.data, {{by}})
  }
  if(is_grouped_df(.data)){
    results <- .data %>%
      doo(colindiag,
          ...,
          n = n)
    return(results %>% set_class(c("tbl_df", "tbl",  "data.frame", "colingroup", "colindiag")))
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
    ultimo <- data.frame(Peso = t(AvAvet[c(nrow(AvAvet)), ])[-c(1), ])
    abs <- data.frame(Peso = abs(ultimo[, "Peso"]))
    rownames(abs) <- rownames(ultimo)
    ultimo <- abs[order(abs[, "Peso"], decreasing = T), , drop = FALSE]
    pesovarname <- paste(rownames(ultimo), collapse = " > ")
    final <- list(cormat = data.frame(cor.x),
                  corlist = results,
                  evalevet = AvAvet,
                  VIF = VIF,
                  CN = NC,
                  det = Det,
                  ncorhigh = ncorhigh,
                  largest_corr = largest_corr,
                  smallest_corr = smallest_corr,
                  weight_var = pesovarname)
    invisible(final)
  }
  if (is.matrix(.data)) {
    out <- internal(.data)
  }
  if (is.data.frame(.data)) {
    if(!missing(...)){
      dfs <-  select_cols(.data, ...) %>%
        select_numeric_cols()
      if(has_na(dfs)){
        dfs <- remove_rows_na(dfs)
        has_text_in_num(dfs)
      }
    } else{
      dfs <- select_numeric_cols(.data)
      if(has_na(dfs)){
        dfs <- remove_rows_na(dfs)
        has_text_in_num(dfs)
      }
    }
    out <- internal(dfs)
  }
  invisible(out %>% set_class("colindiag"))
}


#' Print an object of class colindiag
#'
#' Print the \code{colindiag} object in two ways. By default, the results
#' are shown in the R console. The results can also be exported to the directory
#' into a *.txt file.
#'
#'
#' @param x The object of class \code{colindiag}
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print colindiag
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' col <- colindiag(data_ge2)
#' print(col)
#' }
print.colindiag <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "colindiag print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  if(has_class(x, "colingroup")){
    for (i in 1:nrow(x)){
      df <- x[i,]
      names <-
        select(df, -data) %>%
        concatenate(everything(), pull = TRUE)
      cat("Level:", names, "\n")
      cat("---------------------------------------------------------------------------\n")
      dat <- df[["data"]][[1]]
      CN <- dat$CN
      VIF <- dat$VIF
    if (CN > 1000) {
      cat(paste0("Severe multicollinearity in the matrix! Pay attention on the variables listed bellow\n",
                 "CN = ", round(CN, digits), "\n"))
    }
    if (CN < 100) {
      cat(paste0("Weak multicollinearity in the matrix\n",
                 "CN = ", round(CN, digits), "\n"))
    }
    if (CN > 100 & CN < 1000) {
      cat(paste0("The multicollinearity in the matrix should be investigated.\n",
                 "CN = ", round(CN, digits), "\n", "Largest VIF = ",
                 max(VIF), "\n"))
    }
      cat(paste0("Matrix determinant: ", round(dat$det, 7)), "\n")
      cat(paste0("Largest correlation: ", dat$largest_corr), "\n")
      cat(paste0("Smallest correlation: ", dat$smallest_corr), "\n")
      cat(paste0("Number of VIFs > 10: ", length(which(VIF > 10))), "\n")
      cat(paste0("Number of correlations with r >= |0.8|: ",dat$ncorhigh), "\n")
      cat(paste0("Variables with largest weight in the last eigenvalues: \n", dat$weight_var), "\n")
      cat("---------------------------------------------------------------------------\n")
    }
  } else{
    CN <- x$CN
    VIF <- x$VIF
    if (CN > 1000) {
      cat(paste0("Severe multicollinearity in the matrix! Pay attention on the variables listed bellow\n",
                 "CN = ", round(CN, digits), "\n"))
    }
    if (CN < 100) {
      cat(paste0("Weak multicollinearity in the matrix\n",
                 "CN = ", round(CN, digits), "\n"))
    }
    if (CN > 100 & CN < 1000) {
      cat(paste0("The multicollinearity in the matrix should be investigated.\n",
                 "CN = ", round(CN, digits), "\n", "Largest VIF = ",
                 max(VIF), "\n"))
    }
    cat(paste0("Matrix determinant: ", round(x$det, 7)), "\n")
    cat(paste0("Largest correlation: ", x$largest_corr), "\n")
    cat(paste0("Smallest correlation: ", x$smallest_corr), "\n")
    cat(paste0("Number of VIFs > 10: ", length(which(VIF > 10))), "\n")
    cat(paste0("Number of correlations with r >= |0.8|: ",x$ncorhigh), "\n")
    cat(paste0("Variables with largest weight in the last eigenvalues: \n", x$weight_var), "\n")
  }
  if (export == TRUE) {
    sink()
  }
}
