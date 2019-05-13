#' Summary an Annicchiarico object
#'
#' Summary the \code{Annicchiarico} object in two ways. By default, the results
#' are shown in the R console. The results can also be exported to the
#' directory into a *.txt file.
#'
#'
#' @param object The \code{Annicchiarico} object
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported
#' to the working directory.
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Other arguments of the function
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method summary Annicchiarico
#' @export
#' @examples
#'
#' library(metan)
#' Ann = Annicchiarico(data_ge2,
#'                     env = ENV,
#'                     gen = GEN,
#'                     rep = REP,
#'                     resp = PH)
#' summary(Ann)
#'
summary.Annicchiarico <- function(object, export = FALSE, file.name = NULL, digits = 3,
                              ...) {
  if (!class(object) == "Annicchiarico") {
    stop("The object must be of class 'Annicchiarico'")
  }

  if (export == TRUE) {
    if (is.null(file.name) == T) {
      file.name <- "Annicchiarico summary"
    } else {
      file.name <- file.name
    }
    sink(paste0(file.name, ".txt"))
    for (i in 1:length(object)) {
      var <- object[[i]]
      cat("Variable", names(object)[i], "\n")
      cat("---------------------------------------------------------------------------\n")
      cat("Environmental index\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$environments, row.names = FALSE, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Analysis for all environments\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$general, row.names = FALSE, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Analysis for favorable environments\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$favorable, row.names = FALSE, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Analysis for unfavorable environments\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$unfavorable, row.names = FALSE, digits = digits)
      cat("\n\n\n")
    }
    sink()
  } else {
    for (i in 1:length(object)) {
      var <- object[[i]]
      cat("Variable", names(object)[i], "\n")
      cat("---------------------------------------------------------------------------\n")
      cat("Environmental index\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$environments, row.names = FALSE, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Analysis for all environments\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$general, row.names = FALSE, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Analysis for favorable environments\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$favorable, row.names = FALSE, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Analysis for unfavorable environments\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$unfavorable, row.names = FALSE, digits = digits)
      cat("\n\n\n")
    }
  }
}




