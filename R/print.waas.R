#' Print an object of class waas
#'
#' Print the \code{waas} object in two ways. By default, the results are shown
#' in the R console. The results can also be exported to the directory.
#'
#'
#' @param x An object of class \code{waas}.
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported to
#'   the working directory
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print waas
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' model <- waas(data_ge,
#'   resp = c(GY, HM),
#'   gen = GEN,
#'   env = ENV,
#'   rep = REP
#' )
#' print(model)
#' }
print.waas <- function(x, export = FALSE, file.name = NULL, digits = 4, ...) {
  if (!class(x) == "waas") {
    stop("The object must be of class 'waas'")
  }
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "waas print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    cat("Individual analysis of variance\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$individual$individual, ...)
    cat("---------------------------------------------------------------------------\n")
    cat("AMMI analysis table\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$anova, ...)
    cat("---------------------------------------------------------------------------\n")
    cat("Weighted average of the absolute scores\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$model, ...)
    cat("---------------------------------------------------------------------------\n")
    cat("Some information regarding the analysis\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$Details, ...)
    cat("\n\n\n")
  }
  if (export == TRUE) {
    sink()
  }
}
