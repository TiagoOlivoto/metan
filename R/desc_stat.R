#' Descriptive statistics
#'
#' Compute the most used measures of central tendency, position, and dispersion.
#'
#'
#' @param .data The data to be analyzed. Must be a dataframe or an object of
#' class \code{split_factors}.
#' @param ... A single variable name or a comma-separated list of unquoted variables names.
#' @param values An alternative way to pass the data to the function. It must be a numeric
#' vector.
#' @param stats The descriptive statistics to show. It must be one of the \code{"AV.dev"} (average deviation),
#' \code{"CI.mean"} (confidence interval for the mean), \code{"CV"} (coefficient of variation),
#' \code{"IQR"} (interquartile range), \code{"Kurt"} (kurtosis), \code{"max"} (maximum value),
#' \code{"mean"} (arithmetic mean), \code{"median"} (median), \code{"min"} (minimum value),
#' \code{"n"} (the length of the data), \code{"norm.pval"} (the p-value for the Shapiro-Wilk test),
#'  \code{"norm.stat"} (the statistic for the Shapiro-Wilk test), \code{"Q2.5"} (the percentile 2.5%),
#'  \code{"Q25"} (the first quartile, Q1), \code{"Q75"} (the third quartile, Q3), \code{"Q97.5"} (the percentile 97.5%),
#'  \code{"SD.amo"} (the sample standard deviation), \code{"SD.pop"} (the population standard deviation),
#'  \code{"SE.mean"} (the standard error of the mean), \code{"skew"} (the skewness), \code{"var.amo"} (the sample variance),
#'  \code{"var.pop"} (the population variance).
#' @param level The confidence level to compute the confidence interval of mean. Defaults to 0.95.
#' @param digits The number of significant figures.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @importFrom tidyr spread gather separate
#' @examples
#'
#' library(metan)
#' library(dplyr)
#'
#' data_ge2 %>%
#'   desc_stat(TKW)
#'
#' vect = data_ge2 %>%
#'   select(TKW) %>%
#'   pull()
#'
#' desc_stat(values = vect)
#'
#' data_ge2 %>%
#'   split_factors(ENV) %>%
#'   desc_stat(EP, EL, PH, CL, CW, NR, NKR,
#'   stats = c("mean", "SE.mean", "CV"))
#'

desc_stat = function(.data = NULL,
                     ...,
                     values = NULL,
                     stats = c("AV.dev", "CI.mean", "CV", "IQR", "Kurt", "max", "mean",
                               "median", "min", "n", "norm.pval", "norm.stat", "Q2.5",
                               "Q25", "Q75", "Q97.5", "SD.amo", "SD.pop", "SE.mean",
                               "skew", "var.amo", "var.pop"),
                     level = 0.95,
                     digits = 4){

if(any(!stats %in% c("AV.dev", "CI.mean", "CV", "IQR", "Kurt", "max", "mean",
                     "median", "min", "n", "norm.pval", "norm.stat", "Q2.5",
                     "Q25", "Q75", "Q97.5", "SD.amo", "SD.pop", "SE.mean",
                     "skew", "var.amo", "var.pop")) == TRUE){
  stop("Invalid value for the argument 'stat'. Allowed values are one of the AV.dev, CI.mean, CV, IQR, Kurt, max, mean, median, min, n, norm.pval, norm.stat, Q2.5, Q25, Q75, Q97.5, SD.amo, SD.pop, SE.mean, skew, var.amo, and var.pop.")
}
  if(!missing(.data) & !missing(values)){
    stop("You can not inform a vector of values if a data frame is used as imput.")
  }
  if(!missing(.data) & missing(...)){
    stop("At least one variable must be informed when using a data frame as input.")
  }
  if(missing(.data) & missing(values)){
    stop("Invalid input. Please, use argument 'values' to inform a numeric vector, or the argument '.data' to inform a dataset.")
  }
  # Helper functions
  var_p = function(x){
    sum((x - mean(x))^2) / length(x)
  }
  var_a = function(x){
    sum((x - mean(x))^2) / (length(x) - 1)
  }
  sd_p = function(x){
    sqrt(sum((x - mean(x))^2) / length(x))
  }
  AV_dev = function(x){
    sum(abs(x - mean(x)))/length(x)
  }
  SE_mean = function(x){
    sd(x)/sqrt(length(x))
  }
  CI_mean = function(x, level = 0.95){
    qt((0.5 + level/2), (length(x) - 1)) * sd(x)/sqrt(length(x))
  }
  CV = function(x){
    sd(x)/mean(x) * 100
  }
  skew <- function(x){
    sum((x - mean(x))^3)/(length(x) * sqrt(var_a(x))^3)
  }
  Kurt <- function(x){
    sum((x - mean(x))^4)/(length(x) * var_a(x)^2) - 3
  }
  norm_st = function(x){
    shapiro.test(x)[[1]]
  }
  norm_pv = function(x){
    shapiro.test(x)[[2]]
  }
  if (any(class(.data) == "split_factors")) {
    dfs = list()
    datain = .data
    for (k in 1:length(.data)) {
      .data = datain[[k]]
      nam = names(datain[k])
      data <- dplyr::select(.data, ...)
      cat("---------------------------------------------------------------------------\n")
      cat("Level", nam, "\n")
      cat("---------------------------------------------------------------------------\n")
      data %<>%
        summarise_all(funs(n = n(),
                           mean = mean,
                           min = min,
                           Q2.5 = quantile(., 0.025),
                           Q25 = quantile(., 0.25),
                           median = median,
                           Q75 = quantile(., 0.75),
                           Q97.5 = quantile(., 0.975),
                           max = max,
                           IQR = IQR,
                           AV.dev = AV_dev,
                           var.pop = var_p,
                           var.amo = var_a,
                           SD.pop = sd_p,
                           SD.amo = sd,
                           SE.mean = SE_mean,
                           CI.mean = CI_mean,
                           skew = skew,
                           Kurt = Kurt,
                           norm.stat = norm_st,
                           norm.pval = norm_pv,
                           CV = CV))
      if(ncol(data) > 22){
        statistics = suppressWarnings(data %>% gather(stat, val) %>%
                                        separate(stat, into = c("var", "stat"), sep = "_") %>%
                                        make_mat(stat, var, val) %>%
                                        as.data.frame() %>%
                                        rownames_to_column("stat") %>%
                                        dplyr::filter(stat %in% stats))
        invisible(statistics)
        print(statistics, digits = digits, row.names = FALSE)
      }
      if(ncol(data) == 22){
        statistics = t(data) %>%
          as.data.frame() %>%
          rownames_to_column("stat") %>%
          dplyr::filter(stat %in% stats)
        invisible(statistics)
        print(statistics, digits = digits, row.names = FALSE)
      }

    }
  } else {
    if(is.null(values) ==  FALSE){
      data = data.frame(values)
    } else {
      data <- dplyr::select(.data, ...)
    }
    data %<>%
      summarise_all(funs(n = n(),
                         mean = mean,
                         min = min,
                         Q2.5 = quantile(., 0.025),
                         Q25 = quantile(., 0.25),
                         median = median,
                         Q75 = quantile(., 0.75),
                         Q97.5 = quantile(., 0.975),
                         max = max,
                         IQR = IQR,
                         AV.dev = AV_dev,
                         var.pop = var_p,
                         var.amo = var_a,
                         SD.pop = sd_p,
                         SD.amo = sd,
                         SE.mean = SE_mean,
                         CI.mean = CI_mean,
                         skew = skew,
                         Kurt = Kurt,
                         norm.stat = norm_st,
                         norm.pval = norm_pv,
                         CV = CV))

    if(ncol(data) > 22){
      statistics = suppressWarnings(data %>% gather(stat, val) %>%
                                      separate(stat, into = c("var", "stat"), sep = "_") %>%
                                      make_mat(stat, var, val) %>%
                                      as.data.frame() %>%
                                      rownames_to_column("stat") %>%
                                      dplyr::filter(stat %in% stats))
      invisible(statistics)
      print(statistics, digits = digits, row.names = FALSE)
    }
    if(ncol(data) == 22){
      statistics = t(data) %>%
        as.data.frame() %>%
        rownames_to_column("stat") %>%
      dplyr::filter(stat %in% stats)
      invisible(statistics)
      print(statistics, digits = digits, row.names = FALSE)
    }
  }
}
