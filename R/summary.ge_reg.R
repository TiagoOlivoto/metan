#' Summary a ge_reg object
#'
#' Summary the \code{ge_reg} object in two ways. By default, the results
#' are shown in the R console. The results can also be exported to the
#' directory into a *.txt file.
#'
#'
#' @param object The \code{ge_reg} object
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported
#' to the working directory.
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Other arguments of the function
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method summary ge_reg
#' @export
#' @examples
#'
#' library(metan)
#' model = ge_reg(data_ge2, ENV, GEN, REP, PH)
#' summary(model)
#'
summary.ge_reg <- function(object, export = FALSE, file.name = NULL, digits = 3,
                              ...) {
  if (!class(object) == "ge_reg") {
    stop("The object must be of class 'ge_reg'")
  }

  if (export == TRUE) {
    if (is.null(file.name) == T) {
      file.name <- "ge_reg summary"
    } else {
      file.name <- file.name
    }
    sink(paste0(file.name, ".txt"))
    for (i in 1:length(object)) {
      var <- object[[i]]
      cat("Variable", names(object)[i], "\n")
      cat("---------------------------------------------------------------------------\n")
      cat("Joint-regression Analysis of variance\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$anova, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Regression parameters\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$regression, row.names = FALSE, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("\n\n\n")
    }
    sink()
  } else {
    for (i in 1:length(object)) {
      var <- object[[i]]
      cat("Variable", names(object)[i], "\n")
      cat("---------------------------------------------------------------------------\n")
      cat("Joint-regression Analysis of variance\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$anova, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Regression parameters\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$regression, row.names = FALSE, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("\n\n\n")
      }
  }
}
