#' Generate correlated variables
#' @description
#' `r badge('experimental')`
#'
#' Generate correlated variables using a vector of know values and desired
#' maximum and minimum correlations
#'
#' @param y A vector to generate variables correlated with.
#' @param min_cor The minimum desired correlation.
#' @param max_cor The maximum desired correlation.
#' @param nvars The number of variables.
#' @param constant A constant. Use `operation` to define which operation is
#'   used.
#' @param operation The operation to be applied to the `constant` value.
#' @param x An optional vector of the same length of `y`. If not informed
#'   (default) then a normally distributed variable (mean = 0, sd = 1) will be
#'   used.
#'
#' @return A data frame with the `y` variable and the correlated variables.
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(metan)
#' y <- rnorm(n = 10)
#' cor_vars <- correlated_vars(y, nvar = 6)
#' plot(cor_vars)
#' }
#'
#'
correlated_vars <- function(y,
                            min_cor = -1,
                            max_cor = 1,
                            nvars,
                            constant = NULL,
                            operation = "*",
                            x = NULL){
  rho <- round(seq(min_cor, max_cor, length.out = nvars), digits = 2)
  if (missing(x)) x <- rnorm(length(y))
  y_res <- residuals(lm(x ~ y))
  df <- cbind(y,
              sapply(rho, function(rho){
                rho * sd(y_res) * y + y_res * sd(y) * sqrt(1 - rho ^ 2)
              })
  ) %>%
    as.data.frame()
  names(df) <- paste(c("y", paste("r", rho, sep = "")))
  if(!is.null(constant)){
    if(length(constant) > 1 & length(constant) != ncol(df)){
      stop("Leng of 'constant' not valid")
    }
    df <- sweep(df, 1, STATS = constant, FUN = operation)
  }
  return(list(df = df) %>% set_class("correlated_vars"))
}

#' Plot an object of class correlated_vars
#'
#' @param x An object of class correlated_vars.
#' @param ... Currently not used.
#'
#' @return An object of class gg.
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' y <- rnorm(n = 10)
#' cor_vars <- correlated_vars(y, nvar = 6)
#' plot(cor_vars)
#' }

plot.correlated_vars <- function(x, ...){
  x[[1]] %>%
    pivot_longer(-y) %>%
    ggplot(aes(y, value, group=name)) +
    geom_smooth(method="lm",
                formula = 'y ~ x',
                color="Black") +
    geom_rug(sides="b") +
    geom_point(alpha=1/2) +
    facet_wrap(~ name, scales="free")
}
