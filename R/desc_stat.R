#'Descriptive statistics
#'
#'Compute the most used measures of central tendency, position, and dispersion.
#'
#'
#'@param .data The data to be analyzed. Must be a dataframe or an object of
#'  class \code{split_factors}.
#'@param ... A single variable name or a comma-separated list of unquoted
#'  variables names. If no variable is informed, all the numeric from
#'  \code{.data} variables will be used.
#'@param values An alternative way to pass the data to the function. It must be
#'  a numeric vector.
#'@param stats The descriptive statistics to show. Defaults to \code{"main"}
#'  (main statistics). Set to \code{"all"} to compute all the statistics bellow
#'  or chose one (or more) of the following: \code{'AV.dev'} (average
#'  deviation), \code{'CI.mean'} (confidence interval for the mean), \code{'CV'}
#'  (coefficient of variation), \code{'IQR'} (interquartile range),
#'  \code{'gm.mean'} (geometric mean), \code{'hm.mean'} (harmonic mean),
#'  \code{'Kurt'} (kurtosis), \code{'mad'} (median absolute deviation),
#'  \code{'max'} (maximum value), \code{'mean'} (arithmetic mean),
#'  \code{'median'} (median), \code{'min'} (minimum value), \code{'n'} (the
#'  length of the data), \code{'norm.pval'} (the p-value for the Shapiro-Wilk
#'  test), \code{'norm.stat'} (the statistic for the Shapiro-Wilk test),
#'  \code{'Q2.5'} (the percentile 2.5%), \code{'Q25'} (the first quartile, Q1),
#'  \code{'Q75'} (the third quartile, Q3), \code{'Q97.5'} (the percentile
#'  97.5%), \code{range} (The range of data), \code{'SD.amo'} (the sample
#'  standard deviation), \code{'SD.pop'} (the population standard deviation),
#'  \code{'SE.mean'} (the standard error of the mean), \code{'skew'} (the
#'  skewness), \code{sum} (the sum of the values), \code{sum.dev} (the sum of
#'  the absolute deviations), \code{sum.sq.dev} (the sum of the squared
#'  deviations), \code{valid.n} (The size of sample with valid number (not NA),
#'  \code{'var.amo'} (the sample variance), \code{'var.pop'} (the population
#'  variance). Use a comma-separated vector of names to select the statistics.
#'  For example, \code{stats = c("median, mean, CV, n")}.
#'@param hist Logical argument defaults to \code{FALSE}. If \code{hist = TRUE}
#'  then a histogram is created for each selected variable.
#'@param level The confidence level to compute the confidence interval of mean.
#'  Defaults to 0.95.
#'@param digits The number of significant digits.
#'@param na.rm Logical. Should missing values be removed?
#'@param verbose Logical argument. If \code{verbose = FALSE} the code is run
#'  silently.
#'@return A tibble with the statistics in the lines and variables in columns. If
#'  \code{.data} is an object of class \code{split_factors}, then the statistics
#'  will be shown for each level of the grouping variable in the function
#'  \code{\link{split_factors}}.
#'@details In cases when the statistics are computed for more than two variables
#'  with data coming from the function \code{\link{split_factors}} the results
#'  are returned in a \emph{long} format. Thus, use the function
#'  \code{\link{desc_wider}} to convert it into a \emph{wide} format (levels of
#'  the factors in the rows and statistics in the columns).
#'@author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'@export
#'@importFrom tidyr spread gather separate pivot_wider pivot_longer
#'@importFrom rlang quo as_label
#' @examples
#' library(metan)
#'
#' desc_stat(data_ge2, TKW)
#'
#' # Compute all statistics
#' # Use a numeric vector as input data
#' vect <- data_ge2$TKW
#' desc_stat(values = vect, stats = "all")
#'
#' # Select specific statistics
#' desc_stat(values = c(12, 13, 19, 21, 8, NA, 23, NA),
#'           na.rm = TRUE,
#'           stats = c('mean, SE.mean, CV, n, valid.n'))
#'
#' # Compute the statistics for each level of "ENV"
#' stats <-
#' data_ge2 %>%
#'   split_factors(GEN) %>%
#'   desc_stat(EP, EL, EH, ED, PH, CD,
#'   stats = c('mean, min, max, median, SE.mean, CI.mean, n, CV'),
#'   verbose = FALSE)
#'
#' # To get a 'wide' format with the statistics of the variable EP above.
#' desc_wider(stats, PH)
#'
desc_stat <- function(.data = NULL, ..., values = NULL, stats = "main", hist = FALSE,
                      level = 0.95, digits = 4, na.rm = FALSE, verbose = TRUE) {
  all_f = c("AV.dev", "CI.mean", "CV", "gm.mean", "hm.mean", "IQR", "Kurt", "mad", "max", "mean", "median", "min", "n", "norm.pval", "norm.stat", "Q2.5", "Q25", "Q75", "Q97.5", "range", "SD.amo", "SD.pop", "SE.mean", "skew", "sum", "sum.dev", "sum.sq.dev", "valid.n", "var.amo", "var.pop")
  main_f = c("CV", "max", "mean", "median", "min", "SE.mean", "var.amo")
  if(!stats %in% c("main", "all")){
    stats = unlist(strsplit(stats, split=", "))
  } else {
    if(stats == "main"){
    stats = main_f
    } else{
    stats = unlist(strsplit(all_f, split=", "))
    }
  }
  if (any(!stats %in% c(all_f, main_f)) == TRUE) {
    stop("Invalid value for the argument 'stat'. Allowed values are one of the AV.dev, CI.mean, CV, IQR, Kurt, mad, max, mean, median, min, n, norm.pval, norm.stat, Q2.5, Q25, Q75, Q97.5, range, SD.amo, SD.pop, SE.mean, skew, sum, sum.dev, sum.sq.dev, valid.n, var.amo, and var.pop. Did you accidentally omit the space between the comma and the following word?")
  }
  if (!missing(.data) & !missing(values)) {
    stop("You can not inform a vector of values if a data frame is used as imput.")
  }
  # if (!missing(.data) & missing(...)) {
  #   stop("At least one variable must be informed when using a data frame as input.")
  # }
  if (missing(.data) & missing(values)) {
    stop("Invalid input. Please, use argument 'values' to inform a numeric vector, or the argument '.data' to inform a dataset.")
  }
  # Helper functions
  var_p <- function(x) {
    sum((x - mean(x, na.rm = na.rm))^2, na.rm = na.rm) / length(which(!is.na(x)))
  }
  range_data <- function(x){
    max(x, na.rm = na.rm) - min(x, na.rm = na.rm)
  }
  var_a <- function(x) {
    sum((x - mean(x, na.rm = na.rm))^2, na.rm = na.rm) / (length(which(!is.na(x))) - 1)
  }
  sd_p <- function(x) {
    sqrt(sum((x - mean(x, na.rm = na.rm))^2, na.rm = na.rm) / length(which(!is.na(x))))
  }
  AV_dev <- function(x) {
    sum(abs(x - mean(x, na.rm = na.rm)), na.rm = na.rm) / length(which(!is.na(x)))
  }
  sum_dev <- function(x) {
    sum(abs(x - mean(x, na.rm = na.rm)), na.rm = na.rm)
  }
  sum_sq_dev <- function(x) {
    sum((x - mean(x, na.rm = na.rm))^2, na.rm = na.rm)
  }
  SE_mean <- function(x) {
    sd(x, na.rm = na.rm) / sqrt(length(which(!is.na(x))))
  }
  CI_mean <- function(x, level = 0.95) {
    qt((0.5 + level/2), (length(which(!is.na(x))) - 1)) * sd(x, na.rm = na.rm)/sqrt(length(which(!is.na(x))))
  }
  CV <- function(x) {
    sd(x, na.rm = na.rm)/mean(x, na.rm = na.rm) * 100
  }
  skew <- function(x) {
    sum((x - mean(x, na.rm = na.rm))^3, na.rm = na.rm)/(length(which(!is.na(x))) * sqrt(var_a(x))^3)
  }
  Kurt <- function(x) {
    sum((x - mean(x, na.rm = na.rm))^4, na.rm = na.rm)/(length(which(!is.na(x))) * var_a(x)^2) - 3
  }
  norm_st <- function(x) {
    shapiro.test(x)[[1]]
  }
  norm_pv <- function(x) {
    shapiro.test(x)[[2]]
  }
  valid_n <- function(x){
    length(which(!is.na(x)))
  }
  if (any(class(.data) == "split_factors")) {
    dfs <- list()
    datain <- .data[[1]]
    for (k in 1:length(datain)) {
      data <- datain[[k]]
      nam <- names(datain[k])
      if (missing(...)){
        data <- select_if(data, is.numeric)
      } else{
        data <- dplyr::select(data, ...)
      }
      if(verbose == TRUE){
        cat("---------------------------------------------------------------------------\n")
        cat(nam, "\n")
        cat("---------------------------------------------------------------------------\n")
      }
      if(any(sapply(data, is.na) == TRUE) & na.rm == FALSE){
        stop("NAs values in data. Set the argument 'na.rm' to TRUE to compute the statistics removing missing values.")
      }
      if(hist == TRUE){
        stats_facet <- data %>% pivot_longer(1:ncol(data))
        nbins <- round(1 + 3.322 * log(nrow(data)), 0)
        plt <-
          ggplot(stats_facet, aes(value)) +
          geom_histogram(color = "black", fill = "gray70", bins = nbins)+
          facet_wrap(.~name, scales = "free")+
          labs(x = "Observed value",
               y = "Count")+
          scale_y_continuous(expand = expand_scale(mult = c(0, .1)))+
          theme(axis.text = element_text(size = 12, colour = "black"),
                axis.title = element_text(size = 14, color = "black"),
                axis.ticks.length = unit(0.25, "cm"),
                axis.ticks = element_line(color = "black"),
                panel.grid.minor = element_blank(),
                strip.background = element_rect(color = "black", fill = NA),
                panel.border = element_rect(color = "black", fill = NA))+
          ggtitle(paste("Histogram for ", nam))
        print(plt)
      }
      data %<>%  summarise_all(funs(n = n(),
                                    valid.n = valid_n,
                                    mean = mean(., na.rm = na.rm),
                                    gm.mean = gm_mean(., na.rm = na.rm),
                                    hm.mean = hm_mean(., na.rm = na.rm),
                                    range = range_data,
                                    min = min(., na.rm = na.rm),
                                    Q2.5 = quantile(., 0.025, na.rm = na.rm),
                                    Q25 = quantile(., 0.25, na.rm = na.rm),
                                    median = median(., na.rm = na.rm),
                                    Q75 = quantile(., 0.75, na.rm = na.rm),
                                    Q97.5 = quantile(., 0.975, na.rm = na.rm),
                                    max = max(., na.rm = na.rm),
                                    IQR = IQR(., na.rm = na.rm),
                                    AV.dev = AV_dev,
                                    mad = mad(., na.rm = na.rm),
                                    var.pop = var_p,
                                    var.amo = var_a,
                                    SD.pop = sd_p,
                                    SD.amo = sd(., na.rm = na.rm),
                                    SE.mean = SE_mean,
                                    CI.mean = CI_mean,
                                    skew = skew,
                                    Kurt = Kurt,
                                    norm.stat = norm_st,
                                    norm.pval = norm_pv,
                                    CV = CV,
                                    sum = sum(., na.rm = na.rm),
                                    sum.dev = sum_dev,
                                    sum.sq.dev = sum_sq_dev))
      if (ncol(data) > 30) {
        statistics <- suppressWarnings(data %>% gather(stat, val) %>%
                                         separate(stat,
                                                  into = c("var", "stat"),
                                                  sep = "_(?=[^_]*$)") %>%
                                         make_mat(stat, var, val) %>%
                                         as_tibble(rownames = NA) %>%
                                         rownames_to_column("Statistic") %>%
                                         dplyr::filter(Statistic %in% stats))
        dfs[[paste(nam)]] <- statistics
        if(verbose == TRUE){
          print(statistics, digits = digits, row.names = FALSE)
        }
      }
      if (ncol(data) == 30) {
        statistics <- t(data) %>%
          as_tibble(rownames = NA) %>%
          rownames_to_column("Statistic") %>%
          dplyr::filter(Statistic %in% stats)
        names(statistics)[ncol(statistics)] <- as_label(quo(...))
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
      arrange(LEVEL) %>%
      separate(LEVEL, into = .data[[2]], sep = "[/|]")

    if(ncol(df) - length(.data[[2]]) == 2){
      invisible(df %>% pivot_wider(names_from = Statistic, values_from = as_label(quo(...))))
    } else{
      invisible(df)
    }

  } else {
    if (is.null(values) == FALSE) {
      if(!is.vector(values)){
        stop("The data in 'values' must be a numeric vector. Please check and fix.")
      }
      data <- data.frame(values)
    } else {
      if (!missing(.data) & missing(...)){
        data <- select_if(.data, is.numeric)
      } else{
        data <- dplyr::select(.data, ...)
      }
    }
    if(any(sapply(data, is.na) == TRUE) & na.rm == FALSE){
      stop("NAs values in data. Set the argument 'na.rm' to TRUE to compute the statistics removing missing values.")
    }
    if(hist == TRUE){
      stats_facet <- data %>% pivot_longer(1:ncol(data))
      nbins <- round(1 + 3.322 * log(nrow(data)), 0)
      plt <-
        ggplot(stats_facet, aes(value)) +
        geom_histogram(color = "black", fill = "gray70", bins = nbins)+
        facet_wrap(.~name, scales = "free")+
        labs(x = "Observed value",
             y = "Count")+
        scale_y_continuous(expand = expand_scale(mult = c(0, .1)))+
        theme(axis.text = element_text(size = 12, colour = "black"),
              axis.title = element_text(size = 14, color = "black"),
              axis.ticks.length = unit(0.25, "cm"),
              axis.ticks = element_line(color = "black"),
              panel.grid.minor = element_blank(),
              strip.background = element_rect(color = "black", fill = NA),
              panel.border = element_rect(color = "black", fill = NA))
      print(plt)
    }
    data %<>%  summarise_all(funs(n = n(),
                                  valid.n = valid_n,
                                  mean = mean(., na.rm = na.rm),
                                  gm.mean = gm_mean(., na.rm = na.rm),
                                  hm.mean = hm_mean(., na.rm = na.rm),
                                  range = range_data,
                                  min = min(., na.rm = na.rm),
                                  Q2.5 = quantile(., 0.025, na.rm = na.rm),
                                  Q25 = quantile(., 0.25, na.rm = na.rm),
                                  median = median(., na.rm = na.rm),
                                  Q75 = quantile(., 0.75, na.rm = na.rm),
                                  Q97.5 = quantile(., 0.975, na.rm = na.rm),
                                  max = max(., na.rm = na.rm),
                                  IQR = IQR(., na.rm = na.rm),
                                  AV.dev = AV_dev,
                                  mad = mad(., na.rm = na.rm),
                                  var.pop = var_p,
                                  var.amo = var_a,
                                  SD.pop = sd_p,
                                  SD.amo = sd(., na.rm = na.rm),
                                  SE.mean = SE_mean,
                                  CI.mean = CI_mean,
                                  skew = skew,
                                  Kurt = Kurt,
                                  norm.stat = norm_st,
                                  norm.pval = norm_pv,
                                  CV = CV,
                                  sum = sum(., na.rm = na.rm),
                                  sum.dev = sum_dev,
                                  sum.sq.dev = sum_sq_dev))

    if (ncol(data) > 30) {
      statistics <- suppressWarnings(data %>% gather(stat, val) %>%
                                       separate(stat,
                                                into = c("var", "stat"),
                                                sep = "_(?=[^_]*$)") %>%
                                       make_mat(stat, var, val) %>%
                                       as_tibble(rownames = NA) %>%
                                       rownames_to_column("Statistic") %>%
                                       dplyr::filter(Statistic %in% stats))
    }
    if (ncol(data) == 30) {
      statistics <- t(data) %>%
        as_tibble(rownames = NA) %>%
        rownames_to_column("Statistic") %>%
        dplyr::filter(Statistic %in% stats)
      names(statistics)[2] <- ifelse(missing(...), "value", as_label(quo(...)))
      if(verbose == TRUE){
        print(statistics, digits = digits, row.names = FALSE)
      }
    }
    invisible(statistics)
  }
}
