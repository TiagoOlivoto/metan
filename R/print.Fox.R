#' Print an object of class Fox
#'
#' Print the \code{Fox} object in two ways. By default, the results
#' are shown in the R console. The results can also be exported to the directory
#' into a *.txt file.
#'
#'
#' @param x The \code{Fox} x
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print Fox
#' @export
#' @examples
#'
#' library(metan)
#' fox = Fox(data_ge2,
#'           env = ENV,
#'           gen = GEN,
#'           rep = REP,
#'           resp = PH)
#' print(fox)
print.Fox <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (!class(x) == "Fox") {
    stop("The object must be of class 'Fox'")
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "Fox print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    cat("Fox TOP third criteria\n")
    cat("---------------------------------------------------------------------------\n")
    print(var)
  }
  cat("\n\n\n")
  if (export == TRUE) {
    sink()
  }
}
