#' Plot an object of class clustering
#'
#' Plot the multitrait stability index using ggplot-generated graphics. .
#'
#'
#' @param x An object of class \code{clustering}
#' @param type The type of plot. Must be one of the 'dendrogram' or
#' 'cophenetic'.
#' @param ... Other arguments passed from the function \code{plot.dendrogram}.
#' @method plot clustering
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'
plot.clustering = function(x, type = "dendrogram", ...){
  if (type == "dendrogram"){
    plot(as.dendrogram(x$hc), ...)
  }
  if (type == "cophenetic"){
    return(x$cofgrap)
  }
}
