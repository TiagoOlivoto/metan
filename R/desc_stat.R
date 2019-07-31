#' Descriptive statistics
#'
#' Compute the most used measures of central tendency, position, and dispersion.
#'
#'
#' @param .data The data to be analyzed. Must be a dataframe or an object of
#' class \code{split_factors}. If
#' @param var The variable to be analyzed.
#' @param values An alternative way to pass the data to the function. It must be a numeric
#' vector.
#' @param level The confidence level to compute the confidence interval of mean. Defaults to 0.95.
#' @param digits The number of significant figures.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#' \dontrun{
#'
#' library(metan)
#' library(dplyr)
#'
#' data_ge2 %>%
#'   desc_stat(var = TKW)
#'
#' vect = data_ge2 %>%
#'   select(TKW) %>%
#'   pull()
#'
#' desc_stat(values = vect)
#'
#' data_ge2 %>%
#'   split_factors(ENV) %>%
#'   desc_stat(var = c(EP, EL, PH, CL, CW, NR, NKR))
#'}
#'

desc_stat = function(.data = NULL, var = NULL, values = NULL, level = 0.95, digits = 4){
  if(!missing(.data) & !missing(values)){
    stop("You can not inform a vector of values if a data frame is used as imput.")
  }
  if(!missing(.data) & missing(var)){
    stop("At least one variable must be informed when using a data frame as input.")
  }
  if(missing(.data) & missing(values)){
    stop("Invalid input. Please, use argument 'values' to inform a numeric vector, or the argument '.data' to inform a dataset.")
  }
  # Functions
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
      data <- as.data.frame(dplyr::select(.data, !!!dplyr::quos(!!dplyr::enquo(var))))
      cat("---------------------------------------------------------------------------\n")
      cat("Level", nam, "\n")
      cat("---------------------------------------------------------------------------\n")
      data %<>%
        summarise_all(funs(n = n(),
                           mean = mean,
                           min = min,
                           q25 = quantile(., 0.25),
                           median = median,
                           q75 = quantile(., 0.75),
                           max = max,
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

      if(ncol(data) > 15){
        print(
          suppressWarnings(data %>% gather(stat, val) %>%
                             separate(stat, into = c("var", "stat"), sep = "_") %>%
                             make_mat(stat, var, val)),
          digits = digits
        )
      }
      if(ncol(data) == 15){
        print(t(data), digits = digits)
      }

    }
  } else {
    if(is.null(values) ==  FALSE){
      data = data.frame(values)
    } else {
      data <- dplyr::select(.data, !!!dplyr::quos(!!dplyr::enquo(var)))
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
      print(
        suppressWarnings(data %>% gather(stat, val) %>%
                           separate(stat, into = c("var", "stat"), sep = "_") %>%
                           make_mat(stat, var, val)),
        digits = digits
      )
    }
    if(ncol(data) == 22){
      print(t(data), digits = digits)
    }
  }
}
