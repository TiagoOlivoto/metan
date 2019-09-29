#' Print an object of class Annicchiarico
#'
#' Print the \code{Annicchiarico} object in two ways. By default, the results
#' are shown in the R console. The results can also be exported to the
#' directory into a *.txt file.
#'
#'
#' @param object The \code{Annicchiarico} object
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported
#' to the working directory.
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output.
#' See \code{\link[tibble]{trunc_mat}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print Annicchiarico
#' @export
#' @examples
#'
#' library(metan)
#' Ann = Annicchiarico(data_ge2,
#'                     env = ENV,
#'                     gen = GEN,
#'                     rep = REP,
#'                     resp = PH)
#' print(Ann)
#'
print.Annicchiarico <- function(object, export = FALSE, file.name = NULL, digits = 3,
                                  ...) {
  if (!class(object) == "Annicchiarico") {
    stop("The object must be of class 'Annicchiarico'")
  }
  backup_options <- options()
  options(pillar.sigfig = digits, ...)
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "Annicchiarico print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  for (i in 1:length(object)) {
    var <- object[[i]]
    cat("Variable", names(object)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    cat("Environmental index\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$environments)
    cat("---------------------------------------------------------------------------\n")
    cat("Analysis for all environments\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$general)
    cat("---------------------------------------------------------------------------\n")
    cat("Analysis for favorable environments\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$favorable)
    cat("---------------------------------------------------------------------------\n")
    cat("Analysis for unfavorable environments\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$unfavorable)
    if (export == TRUE) {
      sink()
    }
  }
  options(backup_options)
}
