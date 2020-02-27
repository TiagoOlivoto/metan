#' Linear and Partial Correlation Coefficients
#'
#' Estimates the linear and partial correlation coefficients using as input a
#' data frame or a correlation matrix.
#'
#'
#' @param .data The data to be analyzed. It must be a symmetric correlation
#'   matrix or a data frame, possible with grouped data passed from
#'   \code{\link[dplyr]{group_by}()}.
#' @param ... Variables to use in the correlation. If \code{...} is null
#'   (Default) then all the numeric variables from \code{.data} are used. It
#'   must be a single variable name or a comma-separated list of unquoted
#'   variables names.
#'@param by One variable (factor) to compute the function by. It is a shortcut
#'  to \code{\link[dplyr]{group_by}()}. To compute the statistics by more than
#'  one grouping variable use that function.
#' @param n If a correlation matrix is provided, then \code{n} is the number of
#'   objects used to compute the correlation coefficients.
#' @param method a character string indicating which correlation coefficient is
#'   to be computed. One of 'pearson' (default), 'kendall', or 'spearman'.
#' @return If \code{.data} is a grouped data passed from
#'   \code{\link[dplyr]{group_by}()} then the results will be returned into a
#'   list-column of data frames, containing:
#' * \strong{linear.mat} The matrix of linear correlation.
#' * \strong{partial.mat} The matrix of partial correlations.
#' * \strong{results} Hypothesis testing for each pairwise comparison.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' partial1 <- lpcor(iris)
#'
#' # Alternatively using the pipe operator %>%
#' partial2 <- iris %>% lpcor()
#'
#' # Using a correlation matrix
#' partial3 <- cor(iris[1:4]) %>%
#'             lpcor(n = nrow(iris),
#'                   verbose = FALSE)
#'
#' # Select all numeric variables and compute the partial correlation
#' # For each level of \code{Species}
#'
#' partial4 <- lpcor(iris, everithig(), by = Species)
#'}
lpcor <- function(.data,
                  ...,
                  by = NULL,
                  n = NULL,
                  method = "pearson",
                  verbose = TRUE) {
  if (!missing(by)){
    if(length(as.list(substitute(by))[-1L]) != 0){
      stop("Only one grouping variable can be used in the argument 'by'.\nUse 'group_by()' to pass '.data' grouped by more than one variable.", call. = FALSE)
    }
    .data <- group_by(.data, {{by}})
  }
  if(is_grouped_df(.data)){
    results <- .data %>%
      doo(lpcor,
          ...,
          n = n,
          method = method,
          verbose = verbose)
    return(set_class(results, c("lpcor", "lpcor_group", "tbl_df", "tbl",  "data.frame")))
  }
  if (!is.matrix(.data) && !is.data.frame(.data)) {
    stop("The object 'x' must be a correlation matrix or a data.frame.")
  }
  if (is.matrix(.data) && is.null(n)) {
    stop("You have a matrix but the sample size used to compute the correlations (n) was not declared.")
  }
  if (is.data.frame(.data) && !is.null(n)) {
    stop("You cannot informe the sample size because a data frame was used as input.")
  }
  internal <- function(x) {
    if (is.matrix(x)) {
      cor.x <- x
    }
    if (is.data.frame(x)) {
      cor.x <- cor(x, method = method)
      n <- nrow(x)
    }
    nvar <- ncol(cor.x)
    df <- n - nvar
    if (df < 0) {
      warning("The number of variables is higher than the number of individuals. Hypothesis testing will not be made.",
              call. = FALSE)
    }
    m <- as.matrix(cor.x)
    X.resid <- -(solve_svd(m))
    diag(X.resid) <- 1/(1 - (1 - 1/diag(solve_svd(m))))
    X.resid <- cov2cor(X.resid)
    results <- data.frame(linear = as.vector(t(m)[lower.tri(m, diag = F)]))
    results <- suppressWarnings(
      results %>%
        mutate(partial = as.vector(t(X.resid)[lower.tri(X.resid, diag = F)]),
               t = partial/(sqrt(1 - partial^2)) * sqrt(n - nvar),
               prob = 2 * (1 - pt(abs(t), df = df)))
    )
    names <- colnames(x)
    combnam <- combn(names, 2, paste, collapse = " x ")
    results <- mutate(results,
                      Pairs = names(sapply(combnam, names))) %>%
      select(Pairs, everything())
    return(list(linear.mat = m,
                partial.mat = as.matrix(X.resid),
                results = as_tibble(results)))
  }
  if (is.matrix(.data)) {
    out <- internal(.data)
  }
  if (is.data.frame(.data)) {
    if (missing(...)){
      dfs <- select_numeric_cols(.data)
    } else{
      dfs <- select(.data, ...) %>%
        select_numeric_cols()
    }
    out <- internal(dfs)
  }
  invisible(add_class(out, class = "lpcor"))
}


#' Print the partial correlation coefficients
#'
#' Print an object of class \code{lpcor} or or \code{lpcor_group} in two ways.
#' By default, the results are shown in the R console. The results can also be
#' exported to the directory.
#'
#'
#' @param x An object of class \code{lpcor} or \code{lpcor_group}.
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported
#'   to the working directory
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print lpcor
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' pcor <- lpcor(data_ge2, NR, NKR, NKE)
#' print(pcor)
#'
#' # Compute the correlations for each level of the factor ENV
#' lpc2 <- lpcor(data_ge2,
#'               NR, NKR, NKE,
#'               by = ENV,
#'               verbose = FALSE)
#' print(lpc2)
#' }
print.lpcor <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "Lpcor print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  if(any(class(x) == "lpcor_group")){
    x %>%
    mutate(af = map(data,
                    ~.x %>%
                      .[[3]])) %>%
      unnest(cols = af) %>%
      remove_cols(data) %>%
      print()
  } else{
    print(x[[3]])
  }
  if (export == TRUE) {
    sink()
  }
}
