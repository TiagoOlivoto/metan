#' Print an object of class ge_stats
#'
#' Print the \code{ge_stats} object in two ways. By default, the results are
#' shown in the R console. The results can also be exported to the directory
#' into a *.txt file.
#'
#'
#' @param x An object of class \code{ge_stats}.
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported
#' to the working directory.
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output.
#' See \code{\link[tibble]{trunc_mat}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print ge_stats
#' @export
#' @examples
#'
#' library(metan)
#' model = ge_stats(.data = data_ge2,
#'                  env = ENV,
#'                  gen = GEN,
#'                  rep = REP,
#'                  resp = PH)
#' print(model)
#'
print.ge_stats <- function(x, export = FALSE, file.name = NULL, digits = 3,
                              ...) {
  if (!class(x) == "ge_stats") {
    stop("The object must be of class 'ge_stats'")
  }
    if (export == TRUE) {
      file.name <- ifelse(is.null(file.name) == TRUE, "ge_stats print", file.name)
      sink(paste0(file.name, ".txt"))
    }
  backup_options <- options()
  options(pillar.sigfig = digits, ...)
    for (i in 1:length(x)) {
      var <- x[[i]]
      cat("Variable", names(x)[i], "\n")
      cat("---------------------------------------------------------------------------\n")
      cat("Individual analysis of variance\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$individual)
      cat("---------------------------------------------------------------------------\n")
      cat("Genotype-vs-environment mean\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$ge_mean)
      cat("---------------------------------------------------------------------------\n")
      cat("Genotype-vs-environment effects\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$ge_effect)
      cat("---------------------------------------------------------------------------\n")
      cat("Genotype + Genotype-vs-environment effects\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$gge_effect)
      cat("---------------------------------------------------------------------------\n")
      cat("Genotype-vs-environment statistics\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$ge_stats)
      cat("---------------------------------------------------------------------------\n")
      cat("Regression-based stability (anova)\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$ER$anova)
      cat("---------------------------------------------------------------------------\n")
      cat("Regression-based stability (parameters)\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$ER$regression)
      cat("---------------------------------------------------------------------------\n")
      cat("Genotypic confidence index (all environments)\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$ANN$general)
      cat("---------------------------------------------------------------------------\n")
      cat("Wricke's Ecovalence\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$Ecoval)
      cat("---------------------------------------------------------------------------\n")
      cat("Nonparametric Superiority index\n")
      cat("---------------------------------------------------------------------------\n")
      print(var$Superiority$index)
      if (export == TRUE) {
        sink()
      }
    }
  options(backup_options)
}
