#' Print an object of class ge_reg
#'
#' Print the \code{ge_reg} object in two ways. By default, the results are shown
#' in the R console. The results can also be exported to the directory into a
#' *.txt file.
#'
#' @param x An object of class \code{ge_reg}.
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print ge_reg
#' @export
#' @examples
#'
#' library(metan)
#' model <- ge_reg(data_ge2, ENV, GEN, REP, PH)
#' print(model)
print.ge_reg <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (!class(x) == "ge_reg") {
    stop("The object must be of class 'ge_reg'")
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "ge_reg print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    cat("Joint-regression Analysis of variance\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$anova)
    cat("---------------------------------------------------------------------------\n")
    cat("Regression parameters\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$regression)
    cat("---------------------------------------------------------------------------\n")
    cat("\n\n\n")
  }
  if (export == TRUE) {
    sink()
  }
}
