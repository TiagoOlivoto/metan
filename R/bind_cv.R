#' Bind cross-validation objects.
#'
#' Helper function that combines objects of class \code{cv_ammi}, \code{cv_ammif}
#'  or \code{cv_blup}. It is useful when looking for a boxplot containing
#' the RMSPD values of those cross-validation procedures.
#'
#'
#' @param ... Input objects of class \code{cv_ammi}, \code{cv_ammif} or
#' \code{cv_blup}.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'
#' \dontrun{
#' library(METAAB)
#' AMMI = cv_ammif(data_ge,
#'                         resp = GY,
#'                         gen = GEN,
#'                         env = ENV,
#'                         rep = REP,
#'                         nboot = 100,
#'                         nrepval = 2)
#' BLUP = cv_blup(data_ge,
#'                         resp = GY,
#'                         gen = GEN,
#'                         env = ENV,
#'                         rep = REP,
#'                         nboot = 100,
#'                         nrepval = 2)
#' bind_data = bind_cv(AMMI, BLUP)
#' plot(bind_data)
#' }
#'
bind_cv = function(...){
  class = list(...)
  if(sum(lapply(class, function(x) !class(x) %in% c("cv_ammi", "cv_ammif", "cv_blup") == TRUE)>0)){
    stop("The object must be of the class 'cv_ammi', 'cv_ammif', or 'cv_blup'.")
  }
  dots = list(...)
  data = do.call(rbind, lapply(dots, function(x){
    x$RMSPD
    })
  )
  return(structure(list(RMSPD = data),
                   class = "cv_ammif"))
}
