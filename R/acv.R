#' Adjusted Coefficient of Variation
#' \loadmathjax
#' @description
#' `r badge('stable')`
#'
#'  Computes the scale-adjusted coefficient of variation,
#'   *acv*, (Doring and Reckling, 2018) to account for the systematic
#'   dependence of \mjseqn{\sigma^2} from \mjseqn{\mu}. The *acv* is
#'   computed as follows:
#'   \mjsdeqn{acv = \frac{\sqrt{10^{\tilde v_i}}}{\mu_i}\times 100}
#'   where \mjseqn{\tilde v_i} is the adjusted logarithm of the variance
#'   computed as:
#' \mjsdeqn{\tilde v_i = a + (b - 2)\frac{1}{n}\sum m_i + 2m_i + e_i}
#' being \mjseqn{a} and \mjseqn{b} the coefficients of the linear regression for
#' \mjseqn{log_{10}} of the variance over the \mjseqn{log_{10}} of the mean;
#' \mjseqn{ m_i} is the \mjseqn{log_{10}} of the mean, and \mjseqn{ e_i} is the
#' Power Law Residuals (POLAR), i.e., the residuals for the previously described
#' regression.
#' @param mean A numeric vector with mean values.
#' @param var A numeric vector with variance values.
#' @param na.rm If `FALSE`, the default, missing values are removed with a
#'   warning. If `TRUE`, missing values are silently removed.
#'
#' @return A tibble with the following columns
#' * **mean** The mean values.
#' * **var** The variance values;
#' * **log10_mean** The base 10 logarithm of mean;
#' * **log10_var** The base 10 logarithm of variance;
#' * **POLAR** The Power Law Residuals;
#' * **cv** The standard coefficient of variation;
#' * **acv** Adjusted coefficient of variation.
#' @references Doring, T.F., and M. Reckling. 2018. Detecting global trends of
#'   cereal yield stability by adjusting the coefficient of variation. Eur. J.
#'   Agron. 99: 30-36. \doi{10.1016/j.eja.2018.06.007}
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @md
#' @export
#'
#' @examples
#' \donttest{
#' ################# Table 1 from Doring and Reckling (2018)  ###########
#'
#' # Mean values
#' u <- c(0.5891, 0.6169, 0.7944, 1.0310, 1.5032, 3.8610, 4.6969, 6.1148,
#'        7.1526, 7.5348, 1.2229, 1.6321, 2.4293, 2.5011, 3.0161)
#'
#' # Variances
#' v <- c(0.0064, 0.0141, 0.0218, 0.0318, 0.0314, 0.0766, 0.0620, 0.0822,
#'        0.1605, 0.1986, 0.0157, 0.0593, 0.0565, 0.1997, 0.2715)
#'
#' library(metan)
#' acv(u, v)
#' }
acv <- function(mean, var, na.rm = FALSE) {
  if(na.rm == FALSE & has_na(mean) | na.rm == FALSE & has_na(var)){
    warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
    mean <- na.omit(mean)
    var <- na.omit(var)
  }
  if(na.rm == TRUE){
    mean <- na.omit(mean)
    var <- na.omit(var)
  }
  if(!is.numeric(mean) | !is.numeric(var)){
    stop("Argument 'mean' and 'var' must be numeric.")
  }
  if(length(mean) != length(var)){
    stop("'mean' and 'var' must have the same length.")
  }
  mi <- log10(mean)
  m <- mean(mi)
  vi <- log10(var)
  mod <- lm(vi ~ mi)
  a <- coef(mod)[[1]]
  b <- coef(mod)[[2]]
  ui <- residuals(mod)
  a_adj <- (a + (b - 2) * mean(mi))
  vi_adj <- a_adj + 2 * mi + ui
  acv <- sqrt(10 ^ vi_adj) / mean * 100
  results <-
    tibble(
      mean = mean,
      var = var,
      log10_mean = mi,
      log10_var = vi,
      POLAR = ui,
      cv = sqrt(var) /  mean  * 100,
      acv = acv
    )
  return(results)
}
