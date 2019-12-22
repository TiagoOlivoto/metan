#' Print an object of class performs_ammi
#'
#' Print the \code{performs_ammi} object in two ways. By default, the results are shown
#' in the R console. The results can also be exported to the directory.
#'
#'
#' @param x An object of class \code{performs_ammi}.
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported to
#'   the working directory
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble]{formatting}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print performs_ammi
#' @export
#' @examples
#'
#' library(metan)
#' model = performs_ammi(data_ge, ENV, GEN, REP,
#'                       resp = c(GY, HM))
#' print(model)
print.performs_ammi <- function(x, export = FALSE, file.name = NULL, digits = 4, ...) {
  if (!class(x) == "performs_ammi") {
    stop("The object must be of class 'performs_ammi'")
  }
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "performs_ammi print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  on.exit(options(options()))
  options(pillar.sigfig = digits, ...)
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    cat("AMMI analysis table\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$ANOVA)
    cat("---------------------------------------------------------------------------\n")
    cat("Scores for genotypes and environments\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$model)
    cat("\n\n\n")
  }
  if (export == TRUE) {
    sink()
  }
}
