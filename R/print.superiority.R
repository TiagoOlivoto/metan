#' Print an object ofclass \code{superiority}
#'
#' Print the \code{superiority} object in two ways. By default, the results are
#' shown in the R console. The results can also be exported to the directory
#' into a *.txt file.
#'
#'
#' @param x An object of class \code{superiority}.
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Other arguments of the function
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print superiority
#' @export
#' @examples
#'
#' library(metan)
#' model <- superiority(data_ge2, ENV, GEN, REP, PH)
#' print(model)
print.superiority <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (!class(x) == "superiority") {
    stop("The object must be of class 'superiority'")
  }
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "superiority summary", file.name)
    sink(paste0(file.name, ".txt"))
  }
  backup_options <- options()
  options(pillar.sigfig = digits, ...)
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    cat("Superiority index considering all, favorable and unfavorable environments\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$index)
    cat("---------------------------------------------------------------------------\n")
    if (export == TRUE) {
      sink()
    }
  }
  options(backup_options)
}
