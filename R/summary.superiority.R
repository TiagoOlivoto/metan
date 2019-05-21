#' Summary a superiority object
#'
#' Summary the \code{superiority} object in two ways. By default, the results
#' are shown in the R console. The results can also be exported to the
#' directory into a *.txt file.
#'
#'
#' @param object The \code{superiority} object
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported
#' to the working directory.
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Other arguments of the function
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method summary superiority
#' @export
#' @examples
#'
#' library(metan)
#' model = superiority(data_ge2, ENV, GEN, REP, PH)
#' summary(model)
#'
summary.superiority <- function(object, export = FALSE, file.name = NULL, digits = 3,
                           ...) {
  if (!class(object) == "superiority") {
    stop("The object must be of class 'superiority'")
  }

  if (export == TRUE) {
    if (is.null(file.name) == T) {
      file.name <- "superiority summary"
    } else {
      file.name <- file.name
    }
    sink(paste0(file.name, ".txt"))
    for (i in 1:length(object)) {
      var <- object[[i]]
      cat("Variable", names(object)[i], "\n")
      cat("---------------------------------------------------------------------------\n")
      cat("Superiority index considering all, favorable and unfavorable environments\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$index, row.names = FALSE, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("\n\n\n")
    }
    sink()
  } else {
    for (i in 1:length(object)) {
      var <- object[[i]]
      cat("Variable", names(object)[i], "\n")
      cat("---------------------------------------------------------------------------\n")
      cat("Superiority index considering all, favorable and unfavorable environments\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$index, row.names = FALSE, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("\n\n\n")
    }
  }
}
