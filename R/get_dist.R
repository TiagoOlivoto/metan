#' Get a distance matrix
#'
#' Get the distance matrices from objects fitted with the function
#' \code{\link{clustering}}. This is especially useful to get distance matrices
#' from several objects to be further analyzed using \code{\link{pairs_mantel}}.
#'
#' @param ... Object(s) of class \code{clustering}.]
#' @param digits The number of significant figures. Defaults to \code{2}.
#' @return A list of class \code{clustering}.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#' \donttest{
#' library(metan)
#'d <- data_ge2 %>%
#'       means_by(GEN) %>%
#'       column_to_rownames("GEN") %>%
#'       clustering()
#'get_dist(d)
#'}
#'
get_dist <- function(..., digits = 2) {
  df <- list(...)
  if(!any(sapply(df, class) == "clustering")){
    stop("Only objects of class 'clustering' are allowed.")
  }
  model_names <- unlist(lapply(substitute(list(...))[-1], deparse))
  res <- lapply(df, function(x) {
    x[["distance"]] %>% round_cols(digits = digits)
  })
  names(res) <- model_names
  return(set_class(res, "clustering"))
}
