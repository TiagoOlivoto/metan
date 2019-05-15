#' Rescale a continuous vector to have specified minimum and maximum values
#'
#' Helper function used in the WAASB package. It rescales a continuous variable
#' to have specified minimum and maximum var. Missing values are not
#' allowed.
#'
#' The function rescale a continuous variable as follows:
#' \deqn{Rv_i = (Nmax - Nmin)/(Omax - Omin) * (O_i - Omax) + Nmax}
#' Where \eqn{Rv_i} is the rescaled value of the ith position of the variable/
#' vector; \eqn{Nmax} and \eqn{Nmin} are the new maximum and minimum values;
#' \eqn{Omax and Omin} are the maximum and minimum values of the original data,
#' and \eqn{O_i} is the ith value of the original data.
#'
#' There are basically two options to use \code{resca} to reescale a variable.
#' The first is passing a data frame to \code{.data} argument and selecting a
#' continuous variable to be scaled using \code{var}. The function will return
#' the original data frame with an additional variable (rescaled) which is the
#' rescaled value for the variable input. The second option is pass a numeric
#' vector in the argument \code{values}. The output, of course, will be a numeric
#' vector of rescaled values.
#'
#' @param .data The dataset
#' @param var The continuous variable to manipulate
#' @param values Optional vector of values to manipulate.
#' @param new_min The minimum value of the new scale. Default is 0.
#' @param new_max The maximum value of the new scale. Default is 100
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'
#' resca(values = c(1:10))
#'
#'
resca <- function(.data = NULL, var = NULL, values = NULL, new_min = 0, new_max = 100) {
  if(!missing(.data) && missing(var)){
    stop("You must inform a continous variable in '.data' to rescale.")
  }
  if(!missing(.data) && !missing(values)){
    stop("You can not inform a vector of values if a data frame is used as imput.")
  }
  new_v <- function(v) {
    (new_max - new_min)/(max(values) - min(values)) * (v - max(values)) + new_max
  }
  if(!missing(.data)){
    values =  unlist(dplyr::select(.data, !!enquo(var))[,1])
  } else {
    values = values
  }
  rescaled = sapply(values, new_v)
  if(!missing(.data)){
    out = .data %>% mutate(rescaled = rescaled)
    return(out)
  } else{
    return(rescaled)
  }
}
