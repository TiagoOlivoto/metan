#'Descriptive statistics
#' @description
#' `r badge('stable')`
#'
#'* `desc_stat()` Computes the most used measures of central tendency,
#'position, and dispersion.
#'* `desc_wider()` is useful to put the variables in columns and grouping
#'variables in rows. The table is filled with a statistic chosen with the
#'argument `stat`.
#'
#'@name desc_stat
#'@param .data The data to be analyzed. It can be a data frame (possible with
#'  grouped data passed from [dplyr::group_by()] or a numeric
#'  vector. For `desc_wider()` `.data` is an object of class
#'  `desc_stat`.
#'@param ... A single variable name or a comma-separated list of unquoted
#'  variables names. If no variable is informed, all the numeric variables from
#'  `.data` will be used. Select helpers are allowed.
#'@param by One variable (factor) to compute the function by. It is a shortcut
#'  to [dplyr::group_by()]. To compute the statistics by more than
#'  one grouping variable use that function.
#'@param stats The descriptive statistics to show. This is used to filter the
#'  output after computation. Defaults to `"main"` (cv, max, mean median,
#'  min, sd.amo, se, ci ). Other allowed values are `"all"` to
#'  show all the statistics, `"robust"` to show robust statistics,
#'  `"quantile"` to show quantile statistics, or chose one (or more) of the
#'  following:
#'  * `"av.dev"`: average deviation.
#'  * `"ci.t"`: t-interval (95% confidence interval) of the mean.
#'  * `"ci.z"`: z-interval (95% confidence interval) of the mean.
#'  * `"cv"`: coefficient of variation.
#'  * `"iqr"`: interquartile range.
#'  * `"gmean"`: geometric mean.
#'  * `"hmean"`: harmonic mean.
#'  * `"Kurt"`: kurtosis.
#'  * `"mad"`: median absolute deviation.
#'  * `"max"`: maximum value.
#'  * `"mean"`: arithmetic mean.
#'  * `"median"`: median.
#'  * `"min"`: minimum value.
#'  * `"n"`: the length of the data.
#'  * `"n.valid"`: The valid (Not `NA`) number of elements
#'  * `"n.missing"`: The number of missing values
#'  * `"n.unique"`: The length of unique elements.
#'  * `"ps"`: the pseudo-sigma (iqr / 1.35).
#'  * `"q2.5", "q25", "q75", "q97.5"`: the percentile 2.5\%, first
#'  quartile, third quartile, and percentile 97.5\%, respectively.
#'  * `range`: The range of data).
#'  * `"sd.amo", "sd.pop"`: the sample and population standard deviation.
#'  * `"se"`: the standard error of the mean.
#'  * `"skew"`: skewness.
#'  * `"sum"`. the sum of the values.
#'  * `"sum.dev"`: the sum of the absolute deviations.
#'  * `"ave.sq.dev"`: the average of the squared deviations.
#'  * `"sum.sq.dev"`: the sum of the squared deviations.
#'  * `"n.valid"`: The size of sample with valid number (not NA).
#'  * `"var.amo", "var.pop"`: the sample and population variance.
#'
#'  Use a names to select the statistics. For example, `stats = c("median,
#'  mean, cv, n")`. Note that the statistic names **are not**
#'  case-sensitive. Both comma or space can be used as separator.
#'@param hist Logical argument defaults to `FALSE`. If `hist = TRUE`
#'  then a histogram is created for each selected variable.
#'@param level The confidence level to compute the confidence interval of mean.
#'  Defaults to 0.95.
#'@param digits The number of significant digits.
#'@param na.rm Logical. Should missing values be removed? Defaults to `FALSE`.
#'@param verbose Logical argument. If `verbose = FALSE` the code is run
#'  silently.
#' @param plot_theme The graphical theme of the plot. Default is
#'   `plot_theme = theme_metan()`. For more details, see
#'   [ggplot2::theme()].
#'@param which A statistic to fill the table.
#'@return
#' * `desc_stats()` returns a tibble with the statistics in the columns and
#' variables (with possible grouping factors) in rows.
#' * `desc_wider()` returns a tibble with variables in columns and grouping
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
#'           stats = c('mean, se, cv, n, n.valid'),
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
  all <- c("av.dev", "ci.t", "ci.z", "cv", "gmean", "hmean", "iqr", "kurt", "mad", "max", "mean", "median", "min", "n", "n.valid", "n.missing", "n.unique", "ps", "q2.5", "q25", "q75", "q97.5", "range", "sd.amo", "sd.pop", "se", "skew", "sum", "sum.dev", "ave.dev", "sum.sq.dev",  "var.amo", "var.pop")
  stats <- strsplit(
    case_when(
      all_lower_case(stats) == "main" ~ c("cv, max, mean, median, min, sd.amo, se, ci.t"),
      all_lower_case(stats) == "all" ~ c("av.dev, ci.t, ci.z, cv, gmean, hmean, iqr, kurt, mad, max, mean, median, min, n, n.valid, n.missing, n.unique, ps, q2.5, q25, q75, q97.5, range, sd.amo, sd.pop, se, skew, sum, sum.dev, ave.dev, sum.sq.dev, n.valid, var.amo, var.pop"),
      all_lower_case(stats) == "robust" ~ c("n, median, iqr, ps"),
      all_lower_case(stats) == "quantile" ~ c("n, min, q25, median, q75, max"),
      TRUE ~ all_lower_case(stats)
    ), "\\s*(\\s|,)\\s*")[[1]]

  if (!any(stats %in% c("all", "main", "robust", "quantile", all))) {
    stop("Invalid value for the argument 'stat'. Allowed values are:\nav.dev, ci.t, ci.z, cv, iqr, kurt, mad, max, mean, median, min, n, n.valid, n.missing, n.unique, q2.5, q25, q75, q97.5, range, sd.amo, sd.pop, se, skew, sum, sum.dev, ave.dev, sum.sq.dev, n.valid, var.amo, and var.pop.\nAlternatively, you can set the following groups of statistics:\n'main', 'all', 'robust', or 'quantile'.", call. = FALSE)
  }
  if(has_class(.data, "numeric")){
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
  if(na.rm == FALSE & has_na(data)){
    warning("NA values removed to compute the function. Use 'na.rm = TRUE' to suppress this warning.", call. = FALSE)
    na.rm <- TRUE
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
      scale_y_continuous(expand = expansion(mult = c(0, .1)))+
      plot_theme
    print(plt)
  }
  results <-
    data %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    group_by(variable) %>%
    summarise(across(value,
                     list(n = ~n(),
                          n.valid = ~n_valid(., na.rm = na.rm),
                          n.missing = ~n_missing(., na.rm = na.rm),
                          n.unique = ~n_unique(., na.rm = na.rm),
                          mean = ~mean(., na.rm = na.rm),
                          gmean = ~gmean(., na.rm = na.rm),
                          hmean = ~hmean(., na.rm = na.rm),
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
                          ps = ~pseudo_sigma(., na.rm = na.rm),
                          var.pop = ~var_pop(., na.rm = na.rm),
                          var.amo = ~var_amo(., na.rm = na.rm),
                          sd.pop = ~sd_pop(., na.rm = na.rm),
                          sd.amo = ~sd_amo(., na.rm = na.rm),
                          se = ~sem(., na.rm = na.rm),
                          ci.t = ~ci_mean_t(., na.rm = na.rm, level = level),
                          ci.z = ~ci_mean_z(., na.rm = na.rm, level = level),
                          skew = ~skew(., na.rm = na.rm),
                          kurt = ~kurt(., na.rm = na.rm),
                          cv = ~cv(., na.rm = na.rm),
                          sum = ~sum(., na.rm = na.rm),
                          ave.dev = ~ave_dev(., na.rm = na.rm),
                          sum.dev = ~sum_dev(., na.rm = na.rm),
                          sum.sq.dev = ~sum_sq_dev(., na.rm = na.rm)),
                     .names = "{fn}"),
              .groups = "drop") %>%
    select(group_cols(), variable, {{stats}}) %>%
    round_cols(digits = digits)
  return(results)
}
#'@name desc_stat
#'@export
desc_wider <- function(.data, which) {
  factors <- .data %>% select_non_numeric_cols()
  numeric <- .data %>% select({{which}})
  cbind(factors, numeric) %>%
    pivot_wider(values_from = {{which}}, names_from = variable)
}
