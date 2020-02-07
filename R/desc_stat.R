#'Descriptive statistics
#'
#'* \code{desc_stat()} Computes the most used measures of central tendency,
#'position, and dispersion.
#'* \code{desc_wider()} is useful to put the variables in columns and grouping
#'variables in rows. The table is filled with a statistic chosen with the
#'argument \code{stat}.
#'
#'@name desc_stat
#'@param .data The data to be analyzed. It can be a data frame (possible with
#'  grouped data passed from \code{\link[dplyr]{group_by}()} or a numeric
#'  vector. For \code{desc_wider()} \code{.data} is an object of class
#'  \code{desc_stat}.
#'@param ... A single variable name or a comma-separated list of unquoted
#'  variables names. If no variable is informed, all the numeric variables from
#'  \code{.data} will be used. Select helpers are allowed.
#'@param by One variable (factor) to compute the function by. It is a shortcut
#'  to \code{\link[dplyr]{group_by}()}. To compute the statistics by more than
#'  one grouping variable use that function.
#'@param values Deprecated argument. It will be retired in the next release.
#'@param stats The descriptive statistics to show. This is used to filter the
#'  output after computation. Defaults to \code{"main"} (cv, max, mean median,
#'  min, sd.amo, se, ci ). Other allowed values are \code{"all"} to
#'  show all the statistics, \code{"robust"} to show robust statistics,
#'  \code{"quantile"} to show quantile statistics, or chose one (or more) of the
#'  following:
#'  * \code{"av.dev"}: average deviation.
#'  * \code{"ci"}: 95 percent confidence interval of the mean.
#'  * \code{"cv"}: coefficient of variation.
#'  * \code{"iqr"}: interquartile range.
#'  * \code{"gm.mean"}: geometric mean.
#'  * \code{"hm.mean"}: harmonic mean.
#'  * \code{"Kurt"}: kurtosis.
#'  * \code{"mad"}: median absolute deviation.
#'  * \code{"max"}: maximum value.
#'  * \code{"mean"}: arithmetic mean.
#'  * \code{"median"}: median.
#'  * \code{"min"}: minimum value.
#'  * \code{"n"}: the length of the data.
#'  * \code{"q2.5", "q25", "q75", "q97.5"}: the percentile 2.5\%, first
#'  quartile, third quartile, and percentile 97.5\%, respectively.
#'  * \code{range}: The range of data).
#'  * \code{"sd.amo", "sd.pop"}: the sample and population standard deviation.
#'  * \code{"se"}: the standard error of the mean.
#'  * \code{"skew"}: skewness.
#'  * \code{"sum"}. the sum of the values.
#'  * \code{"sum.dev"}: the sum of the absolute deviations.
#'  * \code{"sum.sq.dev"}: the sum of the squared deviations.
#'  * \code{"valid.n"}: The size of sample with valid number (not NA).
#'  * \code{"var.amo", "var.pop"}: the sample and population variance.
#'
#'  Use a names to select the statistics. For example, \code{stats = c("median,
#'  mean, cv, n")}. Note that the statistic names \strong{are not}
#'  case-sensitive. Both comma or space can be used as separator.
#'@param hist Logical argument defaults to \code{FALSE}. If \code{hist = TRUE}
#'  then a histogram is created for each selected variable.
#'@param level The confidence level to compute the confidence interval of mean.
#'  Defaults to 0.95.
#'@param digits The number of significant digits.
#'@param na.rm Logical. Should missing values be removed? Defaults to \code{FALSE}.
#'@param verbose Logical argument. If \code{verbose = FALSE} the code is run
#'  silently.
#' @param plot_theme The graphical theme of the plot. Default is
#'   \code{plot_theme = theme_metan()}. For more details, see
#'   \code{\link[ggplot2]{theme}}.
#'@param var Deprecated argument. It will be retired in the next release.
#'@param which A statistic to fill the table.
#'@return
#' * \code{desc_stats()} returns a tibble with the statistics in the columns and
#' variables (with possible grouping factors) in rows.
#' * \code{desc_wider()} returns a tibble with variables in columns and grouping
#' factors in rows.
#'@md
#'@author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'@export
#'@importFrom tidyr separate pivot_wider pivot_longer
#'@importFrom rlang quo as_label
#'@importFrom dplyr group_cols
#' @examples
#' \donttest{
#' library(metan)
#' #===============================================================#
#' # Example 1: main statistics (coefficient of variation, maximum,#
#' # mean, median, minimum, sample standard deviation, standard    #
#' # error and confidence interval of the mean) for all numeric    #
#' # variables in data                                             #
#' #===============================================================#
#'
#' desc_stat(data_ge2)
#'
#' #===============================================================#
#' #Example 2: robust statistics using a numeric vector as input   #
#' # data
#' #===============================================================#
#' vect <- data_ge2$TKW
#' desc_stat(vect, stats = "robust")
#'
#' #===============================================================#
#' # Example 3: Select specific statistics. In this example, NAs   #
#' # are removed before analysis with a warning message            #
#' #===============================================================#
#' desc_stat(c(12, 13, 19, 21, 8, NA, 23, NA),
#'           stats = c('mean, se, cv, n, valid.n'),
#'           na.rm = TRUE)
#'
#' #===============================================================#
#' # Example 4: Select specific variables and compute statistics by#
#' # levels of a factor variable (GEN)                             #
#' #===============================================================#
#' stats <-
#'   desc_stat(data_ge2,
#'             EP, EL, EH, ED, PH, CD,
#'             by = GEN)
#' stats
#'
#' # To get a 'wide' format with the maximum values for all variables
#' desc_wider(stats, max)
#'
#' #===============================================================#
#' # Example 5: Compute all statistics for all numeric variables   #
#' # by two or more factors. Note that group_by() was used to pass #
#' # grouped data to the function desc_stat()                      #
#' #===============================================================#
#'
#' data_ge2 %>%
#'   group_by(ENV, GEN) %>%
#'   desc_stat()
#'
#'}
desc_stat <- function(.data = NULL,
                      ...,
                      by = NULL,
                      values = "deprecated",
                      stats = "main",
                      hist = FALSE,
                      level = 0.95,
                      digits = 4,
                      na.rm = FALSE,
                      verbose = TRUE,
                      plot_theme = theme_metan()) {
  if (!missing(by)){
    if(length(as.list(substitute(by))[-1L]) != 0){
      stop("Only one grouping variable can be used in the argument 'by'.\nUse 'group_by()' to pass '.data' grouped by more than one variable.", call. = FALSE)
    }
    .data <- group_by(.data, {{by}})
  }
  if(is_grouped_df(.data)){
    results <- .data %>%
      doo(desc_stat, ...,
          stats = stats,
          hist = hist,
          level = level,
          digits = digits,
          na.rm = na.rm,
          verbose = verbose,
          plot_theme = plot_theme)
    return(results)
  }
  all <- c("av.dev", "ci", "cv", "gm.mean", "hm.mean", "iqr", "kurt", "mad", "max", "mean", "median", "min", "n", "q2.5", "q25", "q75", "q97.5", "range", "sd.amo", "sd.pop", "se", "skew", "sum", "sum.dev", "sum.sq.dev", "valid.n", "var.amo", "var.pop")
  stats <- strsplit(
    case_when(
      all_lower_case(stats) == "main" ~ c("cv, max, mean, median, min, sd.amo, se, ci"),
      all_lower_case(stats) == "all" ~ c("av.dev, ci, cv, gm.mean, hm.mean, iqr, kurt, mad, max, mean, median, min, n, q2.5, q25, q75, q97.5, range, sd.amo, sd.pop, se, skew, sum, sum.dev, sum.sq.dev, valid.n, var.amo, var.pop"),
      all_lower_case(stats) == "robust" ~ c("n, median, iqr"),
      all_lower_case(stats) == "quantile" ~ c("n, min, q25, median, q75, max"),
      TRUE ~ all_lower_case(stats)
    ), "\\s*(\\s|,)\\s*")[[1]]

  if (!any(stats %in% c("all", "main", "robust", "quantile", all))) {
    stop("Invalid value for the argument 'stat'. Allowed values are:\nav.dev, ci, cv, iqr, kurt, mad, max, mean, median, min, n, q2.5, q25, q75, q97.5, range, sd.amo, sd.pop, se, skew, sum, sum.dev, sum.sq.dev, valid.n, var.amo, and var.pop.\nAlternatively, you can set the following groups of statistics:\n'main', 'all', 'robust', or 'quantile'.", call. = FALSE)
  }
  if(any(class(.data) == "numeric")){
    .data <- data.frame(val = .data)
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))

  if (!missing(.data) & missing(...)){
    data <- select_numeric_cols(.data)
  } else{
    data <- select(.data, ...) %>%
      select_numeric_cols()
  }
  if(has_na(data) == TRUE && na.rm == FALSE){
    stop("NA values in data. Use 'na.rm = TRUE' to remove NAs from analysis.\nTo remove rows with NA use `remove_rows_na()'. \nTo remove columns with NA use `remove_cols_na()'.", call. = FALSE)
  }
  if(hist == TRUE){
    stats_facet <- data %>% select_numeric_cols() %>%  pivot_longer(everything())
    nbins <- round(1 + 3.322 * log(nrow(data)), 0)
    plt <-
      ggplot(stats_facet, aes(value)) +
      geom_histogram(color = "black", fill = "gray70", bins = nbins)+
      facet_wrap(.~name, scales = "free")+
      labs(x = "Observed value",
           y = "Count")+
      scale_y_continuous(expand = expand_scale(mult = c(0, .1)))+
      plot_theme
    print(plt)
  }
  results <-
    data %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    group_by(variable) %>%
    summarise_all(list(n = ~n(),
                       valid.n = ~valid_n(., na.rm = na.rm),
                       mean = ~mean(., na.rm = na.rm),
                       gm.mean = ~gm_mean(., na.rm = na.rm),
                       hm.mean = ~hm_mean(., na.rm = na.rm),
                       range = ~range_data(., na.rm = na.rm),
                       min = ~min(., na.rm = na.rm),
                       q2.5 = ~quantile(., 0.025, na.rm = na.rm),
                       q25 = ~quantile(., 0.25, na.rm = na.rm),
                       median = ~median(., na.rm = na.rm),
                       q75 = ~quantile(., 0.75, na.rm = na.rm),
                       q97.5 = ~quantile(., 0.975, na.rm = na.rm),
                       max = ~max(., na.rm = na.rm),
                       iqr = ~IQR(., na.rm = na.rm),
                       av.dev = ~av_dev(., na.rm = na.rm),
                       mad = ~mad(., na.rm = na.rm),
                       var.pop = ~var_pop(., na.rm = na.rm),
                       var.amo = ~var_amo(., na.rm = na.rm),
                       sd.pop = ~sd_pop(., na.rm = na.rm),
                       sd.amo = ~sd_amo(., na.rm = na.rm),
                       se = ~sem(., na.rm = na.rm),
                       ci = ~ci_mean(., na.rm = na.rm),
                       skew = ~skew(., na.rm = na.rm),
                       kurt = ~kurt(., na.rm = na.rm),
                       cv = ~cv(., na.rm = na.rm),
                       sum = ~sum(., na.rm = na.rm),
                       sum.dev = ~sum_dev(., na.rm = na.rm),
                       sum.sq.dev = ~sum_sq_dev(., na.rm = na.rm))) %>%
    select(group_cols(), variable, {{stats}}) %>%
    mutate_if(is.numeric, round, digits = digits)
  return(results)
}
#'@name desc_stat
#'@export
desc_wider <- function(.data, which, var = "deprecated") {
  factors = .data %>% select_non_numeric_cols()
  numeric = .data %>% select({{which}})
  cbind(factors, numeric) %>%
    pivot_wider(values_from = {{which}}, names_from = variable)
}
