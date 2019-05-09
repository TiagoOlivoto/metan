#' Rescale a continuous vector to have specified minimum and maximum values
#'
#' Helper function used in the WAASB package. It rescales a continuous vector
#' to have specified minimum and maximum values. Missing values are not
#' allowed.
#'
#'
#' @param values continuous vector of values to manipulate.
#' @param new_min The minimum value of the new scale. Default is 0.
#' @param new_max The maximum value of the new scale. Default is 100
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'
#' resca(1:10)
#' resca(c(20, 30, 40), new_min = 5, new_max = 10)
#'
resca <- function(values, new_min = 0, new_max = 100) {
    new_v <- function(v) {
        (new_max - new_min)/(max(values) - min(values)) * (v - max(values)) + new_max
    }
    return(sapply(values, new_v))
}
