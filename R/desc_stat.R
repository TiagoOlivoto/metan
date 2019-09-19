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
#' @param stats The descriptive statistics to show. Defaults to \code{NULL} (all statistics shown).
#'  It must be one of the \code{'AV.dev'} (average deviation),
#' \code{'CI.mean'} (confidence interval for the mean), \code{'CV'} (coefficient of variation),
#' \code{'IQR'} (interquartile range), \code{'Kurt'} (kurtosis), \code{'max'} (maximum value),
#' \code{'mean'} (arithmetic mean), \code{'median'} (median), \code{'min'} (minimum value),
#' \code{'n'} (the length of the data), \code{'norm.pval'} (the p-value for the Shapiro-Wilk test),
#'  \code{'norm.stat'} (the statistic for the Shapiro-Wilk test), \code{'Q2.5'} (the percentile 2.5%),
#'  \code{'Q25'} (the first quartile, Q1), \code{'Q75'} (the third quartile, Q3), \code{'Q97.5'} (the percentile 97.5%),
#'  \code{range} (The range of data), \code{'SD.amo'} (the sample standard deviation), \code{'SD.pop'} (the population standard deviation),
#'  \code{'SE.mean'} (the standard error of the mean), \code{'skew'} (the skewness), \code{'var.amo'} (the sample variance),
#'  \code{'var.pop'} (the population variance). Use a comma-separated vector of names to select the statistics.
#'  For example, \code{stats = c("media, mean, CV, n")}.
#' @param level The confidence level to compute the confidence interval of mean. Defaults to 0.95.
#' @param digits The number of significant figures.
#' @param verbose Logical argument. If \code{verbose = FALSE} the code is run silently.
#' @return A tibble with the statistics in the lines and variables in columns. If
#' \code{.data} is an object of class \code{split_factors}, then the statistics will be
#' shown for each level of the grouping variable in the function \code{split_factors()}.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @importFrom tidyr spread gather separate
#' @importFrom rlang quo quo_text
#' @examples
#' library(metan)
#' library(dplyr)
#' desc1 <- data_ge2 %>% desc_stat(TKW)
#'
#' vect <- data_ge2 %>%
#'   select(TKW) %>%
#'   pull()
#'
#' desc2 <- desc_stat(values = vect)
#'
#' desc3 <-
#' data_ge2 %>%
#'   split_factors(ENV) %>%
#'   desc_stat(EP, EL, PH, CL, CW, NR, NKR,
#'   stats = c('mean, SE.mean, CV'),
#'   verbose = FALSE)
#'
desc_stat <- function(.data = NULL, ..., values = NULL, stats = NULL,
                      level = 0.95, digits = 4, verbose = TRUE) {
  test = c("AV.dev", "CI.mean", "CV", "IQR", "Kurt", "max", "mean", "median", "min", "range", "n", "norm.pval", "norm.stat", "Q2.5", "Q25", "Q75", "Q97.5", "SD.amo", "SD.pop", "SE.mean", "skew", "var.amo", "var.pop")
  if(!missing(stats)){
    stats = unlist(strsplit(stats, split=", "))
  } else {
    stats = unlist(strsplit(test, split=", "))
  }
  if (any(!stats %in% test) == TRUE) {
    stop("Invalid value for the argument 'stat'. Allowed values are one of the AV.dev, CI.mean, CV, IQR, Kurt, max, mean, median, min, range, n, norm.pval, norm.stat, Q2.5, Q25, Q75, Q97.5, SD.amo, SD.pop, SE.mean, skew, var.amo, and var.pop. Did you accidentally omit the space between the comma and the following word?")
  }
  if (!missing(.data) & !missing(values)) {
    stop("You can not inform a vector of values if a data frame is used as imput.")
  }
  if (!missing(.data) & missing(...)) {
    stop("At least one variable must be informed when using a data frame as input.")
  }
  if (missing(.data) & missing(values)) {
    stop("Invalid input. Please, use argument 'values' to inform a numeric vector, or the argument '.data' to inform a dataset.")
  }
  # Helper functions
  var_p <- function(x) {
    sum((x - mean(x))^2)/length(x)
  }
  range_data <- function(x){
    max(x) - min(x)
  }
  var_a <- function(x) {
    sum((x - mean(x))^2)/(length(x) - 1)
  }
  sd_p <- function(x) {
    sqrt(sum((x - mean(x))^2)/length(x))
  }
  AV_dev <- function(x) {
    sum(abs(x - mean(x)))/length(x)
  }
  SE_mean <- function(x) {
    sd(x)/sqrt(length(x))
  }
  CI_mean <- function(x, level = 0.95) {
    qt((0.5 + level/2), (length(x) - 1)) * sd(x)/sqrt(length(x))
  }
  CV <- function(x) {
    sd(x)/mean(x) * 100
  }
  skew <- function(x) {
    sum((x - mean(x))^3)/(length(x) * sqrt(var_a(x))^3)
  }
  Kurt <- function(x) {
    sum((x - mean(x))^4)/(length(x) * var_a(x)^2) - 3
  }
  norm_st <- function(x) {
    shapiro.test(x)[[1]]
  }
  norm_pv <- function(x) {
    shapiro.test(x)[[2]]
  }
  if (any(class(.data) == "split_factors")) {
    dfs <- list()
    datain <- .data
    for (k in 1:length(.data)) {
      .data <- datain[[k]]
      nam <- names(datain[k])
      data <- dplyr::select(.data, ...)
      if(verbose == TRUE){
        cat("---------------------------------------------------------------------------\n")
        cat(nam, "\n")
        cat("---------------------------------------------------------------------------\n")
      }
      data %<>% summarise_all(funs(n = n(), mean = mean, range = range_data, min = min,
                                   Q2.5 = quantile(., 0.025), Q25 = quantile(., 0.25), median = median,
                                   Q75 = quantile(., 0.75), Q97.5 = quantile(., 0.975), max = max,
                                   IQR = IQR, AV.dev = AV_dev, var.pop = var_p, var.amo = var_a,
                                   SD.pop = sd_p, SD.amo = sd, SE.mean = SE_mean, CI.mean = CI_mean,
                                   skew = skew, Kurt = Kurt, norm.stat = norm_st, norm.pval = norm_pv,
                                   CV = CV))
      if (ncol(data) > 23) {
        statistics <- suppressWarnings(data %>% gather(stat, val) %>%
                                         separate(stat, into = c("var", "stat"), sep = "_") %>%
                                         make_mat(stat, var, val) %>%
                                         as_tibble(rownames = NA) %>%
                                         rownames_to_column("Statistic") %>%
                                         dplyr::filter(Statistic %in% stats))
        dfs[[paste(nam)]] <- statistics
        if(verbose == TRUE){
          print(statistics, digits = digits, row.names = FALSE)
        }
      }
      if (ncol(data) == 23) {
        statistics <- t(data) %>%
          as_tibble(rownames = NA) %>%
          rownames_to_column("Statistic") %>%
          dplyr::filter(Statistic %in% stats)
        dfs[[paste(nam)]] <- statistics
        if(verbose == TRUE){
          print(statistics, digits = digits, row.names = FALSE)
        }
      }
    }
    df = do.call(rbind, lapply(dfs, function(x){
      x
    })) %>%
      rownames_to_column("LEVEL") %>%
      arrange(Statistic) %>%
      mutate(LEVEL = paste(rep(names(dfs), nrow(dfs[[1]])))) %>%
      arrange(LEVEL)
    invisible(df)
  } else {
    if (is.null(values) == FALSE) {
      data <- data.frame(values)
    } else {
      data <- dplyr::select(.data, ...)
    }
    data %<>% summarise_all(funs(n = n(), mean = mean, range = range_data, min = min, Q2.5 = quantile(.,
                                                                                  0.025), Q25 = quantile(., 0.25), median = median, Q75 = quantile(.,
                                                                                                                                                   0.75), Q97.5 = quantile(., 0.975), max = max, IQR = IQR, AV.dev = AV_dev,
                                 var.pop = var_p, var.amo = var_a, SD.pop = sd_p, SD.amo = sd,
                                 SE.mean = SE_mean, CI.mean = CI_mean, skew = skew, Kurt = Kurt,
                                 norm.stat = norm_st, norm.pval = norm_pv, CV = CV))
    if (ncol(data) > 23) {
      statistics <- suppressWarnings(data %>% gather(stat, val) %>%
                                       separate(stat, into = c("var", "stat"), sep = "_") %>%
                                       make_mat(stat, var, val) %>%
                                       as_tibble(rownames = NA) %>%
                                       rownames_to_column("Statistic") %>%
                                       dplyr::filter(Statistic %in% stats))
    }
    if (ncol(data) == 23) {
      statistics <- t(data) %>%
        as_tibble(rownames = NA) %>%
        rownames_to_column("Statistic") %>%
        dplyr::filter(Statistic %in% stats)
      names(statistics)[2] <- quo_text(quo(...))

    }
    if(verbose == TRUE){
    print(statistics, digits = digits, row.names = FALSE)
    }
    invisible(statistics)
  }
}
