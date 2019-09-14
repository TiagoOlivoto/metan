#' Summarize an object of class can_cor
#'
#' Summarize an object of class \code{can_cor} object in two ways. By default,
#' the results are shown in the R console. The results can also be exported to
#' the directory.
#'
#'
#' @param object The \code{can_cor} object
#' @param export A logical argument. If \code{TRUE|T}, a *.txt file is exported
#' to the working directory
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Other arguments of the function
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method summary can_cor
#' @export
#' @examples
#'
#' library(metan)
#' cc <- can_corr(data_ge2,
#'                FG = c(PH, EH, EP),
#'                SG = c(EL, CL, CD, CW, KW, NR, TKW),
#'                verbose = FALSE)
#' summary(cc)
#' summary(cc, export = TRUE,
#'         file.name = "canonical results")
#'
summary.can_cor <- function(object, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (!class(object) == "can_cor") {
    stop("The object must be of class 'can_cor'")
  }
  if (export == TRUE) {
    if (is.null(file.name) == TRUE) {
      file.name <- "Canonical summary"
    } else {
      file.name <- file.name
    }
    sink(paste0(file.name, ".txt"))
    cat("---------------------------------------------------------------------------\n")
    cat("Matrix (correlation/covariance) between variables of first group (FG)\n")
    cat("---------------------------------------------------------------------------\n")
    print(object$MFG, digits = digits)
    cat("\n---------------------------------------------------------------------------\n")
    cat("Collinearity diagnostic between first group\n")
    cat("---------------------------------------------------------------------------\n")
    colindiag(object$MFG, n = nrow(object$Score_FG))
    cat("\n---------------------------------------------------------------------------\n")
    cat("Matrix (correlation/covariance) between variables of second group (SG)\n")
    cat("---------------------------------------------------------------------------\n")
    print(object$MSG, digits = digits)
    cat("\n---------------------------------------------------------------------------\n")
    cat("Collinearity diagnostic between second group\n")
    cat("---------------------------------------------------------------------------\n")
    colindiag(object$MSG, n = nrow(object$Score_SG))
    cat("\n---------------------------------------------------------------------------\n")
    cat("Matrix (correlation/covariance) between FG and SG)\n")
    cat("---------------------------------------------------------------------------\n")
    print(object$MFG_SG, digits = digits)
    cat("\n---------------------------------------------------------------------------\n")
    cat("Correlation of the canonical pairs and hipotesis testing \n")
    cat("---------------------------------------------------------------------------\n")
    print(object$SigTest, digits = digits)
    cat("\n---------------------------------------------------------------------------\n")
    cat("Canonical coefficients of the first group \n")
    cat("---------------------------------------------------------------------------\n")
    print(object$Coef_FG, digits = digits)
    cat("\n---------------------------------------------------------------------------\n")
    cat("Canonical coefficients of the second group \n")
    cat("---------------------------------------------------------------------------\n")
    print(object$Coef_SG, digits = digits)
    cat("\n---------------------------------------------------------------------------\n")
    cat("Canonical loads of the first group \n")
    cat("---------------------------------------------------------------------------\n")
    print(object$Loads_FG, digits = digits)
    cat("\n---------------------------------------------------------------------------\n")
    cat("Canonical loads of the second group \n")
    cat("---------------------------------------------------------------------------\n")
    print(object$Loads_SG, digits = digits)
    cat("\n-------------------------- End of data ------------------------------------\n\n\n\n")
    sink()
  } else {
      cat("---------------------------------------------------------------------------\n")
      cat("Matrix (correlation/covariance) between variables of first group (FG)\n")
      cat("---------------------------------------------------------------------------\n")
      print(object$MFG, digits = digits)
      cat("\n---------------------------------------------------------------------------\n")
      cat("Collinearity diagnostic between first group\n")
      cat("---------------------------------------------------------------------------\n")
      colindiag(object$MFG, n = nrow(object$Score_FG))
      cat("\n---------------------------------------------------------------------------\n")
      cat("Matrix (correlation/covariance) between variables of second group (SG)\n")
      cat("---------------------------------------------------------------------------\n")
      print(object$MSG, digits = digits)
      cat("\n---------------------------------------------------------------------------\n")
      cat("Collinearity diagnostic between second group\n")
      cat("---------------------------------------------------------------------------\n")
      colindiag(object$MSG, n = nrow(object$Score_SG))
      cat("\n---------------------------------------------------------------------------\n")
      cat("Matrix (correlation/covariance) between FG and SG)\n")
      cat("---------------------------------------------------------------------------\n")
      print(object$MFG_SG, digits = digits)
      cat("\n---------------------------------------------------------------------------\n")
      cat("Correlation of the canonical pairs and hipotesis testing \n")
      cat("---------------------------------------------------------------------------\n")
      print(object$SigTest, digits = digits)
      cat("\n---------------------------------------------------------------------------\n")
      cat("Canonical coefficients of the first group \n")
      cat("---------------------------------------------------------------------------\n")
      print(object$Coef_FG, digits = digits)
      cat("\n---------------------------------------------------------------------------\n")
      cat("Canonical coefficients of the second group \n")
      cat("---------------------------------------------------------------------------\n")
      print(object$Coef_SG, digits = digits)
      cat("\n---------------------------------------------------------------------------\n")
      cat("Canonical loads of the first group \n")
      cat("---------------------------------------------------------------------------\n")
      print(object$Loads_FG, digits = digits)
      cat("\n---------------------------------------------------------------------------\n")
      cat("Canonical loads of the second group \n")
      cat("---------------------------------------------------------------------------\n")
      print(object$Loads_SG, digits = digits)
      cat("\n-------------------------- End of data ------------------------------------\n\n\n\n")  }
}
