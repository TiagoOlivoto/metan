#' Bind cross-validation objects
#'
#' Helper function that combines objects of class \code{cv_ammi}, \code{cv_ammif}
#'  or \code{cv_blup}. It is useful when looking for a boxplot containing
#' the RMSPD values of those cross-validation procedures.
#'
#'
#' @param ... Input objects of class \code{cv_ammi}, \code{cv_ammif} or
#' \code{cv_blup}.
#' @param bind What data should be used? To plot the RMSPD, use 'boot' (default).
#' Use \code{bind = "means"} to return the RMSPD mean for each model.
#' @param sort Used to sort the RMSPD mean in ascending order.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'
#' \dontrun{
#' library(metan)
#' AMMI = cv_ammif(data_ge,
#'                 resp = GY,
#'                 gen = GEN,
#'                 env = ENV,
#'                 rep = REP,
#'                 nboot = 100,
#'                 nrepval = 2)
#' BLUP = cv_blup(data_ge,
#'                resp = GY,
#'                gen = GEN,
#'                env = ENV,
#'                rep = REP,
#'                nboot = 100,
#'                nrepval = 2)
#' bind_data = bind_cv(AMMI, BLUP)
#' plot(bind_data)
#'
#' print(bind_cv(AMMI, BLUP, bind = "means"))
#' }
#'
bind_cv = function(..., bind = "boot", sort = TRUE){
  class = list(...)
  if(sum(lapply(class, function(x) !class(x) %in% c("cv_ammi", "cv_ammif", "cv_blup") == TRUE)>0)){
    stop("The object must be of the class 'cv_ammi', 'cv_ammif', or 'cv_blup'.")
  }
  dots = list(...)
  if (bind == "boot"){
  data = do.call(rbind, lapply(dots, function(x){
    x$RMSPD
    })
  )
  if (sort == TRUE){
    data = data %>% arrange(mean)
  }
  }
  if (bind == "means"){
    data = do.call(rbind, lapply(dots, function(x){
      x$RMSPDmean
    })
    )
  }
  if (sort == TRUE){
    data = data %>% arrange(mean)
  }
  return(structure(list(RMSPD = data),
                   class = "cv_ammif"))
}
