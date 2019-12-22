#' Print an object of class ecovalence
#'
#' Print the \code{ecovalence} object in two ways. By default, the results
#' are shown in the R console. The results can also be exported to the directory
#' into a *.txt file.
#'
#'
#' @param x The \code{ecovalence} x
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print ecovalence
#' @export
#' @examples
#'
#' library(metan)
#' eco <- ecovalence(data_ge2,
#'   env = ENV,
#'   gen = GEN,
#'   rep = REP,
#'   resp = PH
#' )
#' print(eco)
print.ecovalence <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (!class(x) == "ecovalence") {
    stop("The object must be of class 'ecovalence'")
  }
  on.exit(options(options()))
  options(pillar.sigfig = digits, ...)
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "ecovalence print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    cat("Genotypic confidence index\n")
    cat("---------------------------------------------------------------------------\n")
    print(var)
  }
  cat("\n\n\n")
  if (export == TRUE) {
    sink()
  }
}
