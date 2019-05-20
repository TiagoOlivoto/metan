#' Summary a ge_stats object
#'
#' Summary the \code{ge_stats} object in two ways. By default, the results are
#' shown in the R console. The results can also be exported to the directory
#' into a *.txt file.
#'
#'
#' @param object The \code{ge_stats} object
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported
#' to the working directory.
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Other arguments of the function
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method summary ge_stats
#' @export
#' @examples
#'
#' library(metan)
#' model = ge_stats(.data = data_ge2,
#'                  env = ENV,
#'                  gen = GEN,
#'                  rep = REP,
#'                  resp = PH)
#' summary(model)
#'
summary.ge_stats <- function(object, export = FALSE, file.name = NULL, digits = 3,
                              ...) {
  if (!class(object) == "ge_stats") {
    stop("The object must be of class 'ge_stats'")
  }

  if (export == TRUE) {
    if (is.null(file.name) == T) {
      file.name <- "ge_stats summary"
    } else {
      file.name <- file.name
    }
    sink(paste0(file.name, ".txt"))
    for (i in 1:length(object)) {
      var <- object[[i]]
      cat("Variable", names(object)[i], "\n")
      cat("---------------------------------------------------------------------------\n")
      cat("Individual analysis of variance\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$individual$individual, row.names = FALSE, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Genotype-vs-environment mean\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$ge_mean, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Genotype-vs-environment effects\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$ge_effect, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Genotype + Genotype-vs-environment effects\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$gge_effect, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Genotype-vs-environment statistics\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$ge_stats, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Regression-based stability (anova)\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$ER$anova, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Regression-based stability (parameters)\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$ER$regression, digits = digits, row.names = FALSE)
      cat("---------------------------------------------------------------------------\n")
      cat("Genotypic confidence index (all environments)\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$ANN$general, digits = digits, row.names = FALSE)
      cat("---------------------------------------------------------------------------\n")
      cat("Wricke's Ecovalence\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$Ecoval, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Nonparametric Superiority index\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$Superiority$index, digits = digits, row.names = FALSE)
      cat("------------------------------- End of data -----------------------------\n\n\n\n")

    }
    sink()
  } else {
    for (i in 1:length(object)) {
      var <- object[[i]]
      cat("Variable", names(object)[i], "\n")
      cat("---------------------------------------------------------------------------\n")
      cat("Individual analysis of variance\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$individual$individual, row.names = FALSE, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Genotype-vs-environment mean\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$ge_mean, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Genotype-vs-environment effects\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$ge_effect, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Genotype + Genotype-vs-environment effects\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$gge_effect, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Genotype-vs-environment statistics\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$ge_stats, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Regression-based stability (anova)\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$ER$anova, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Regression-based stability (parameters)\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$ER$regression, digits = digits, row.names = FALSE)
      cat("---------------------------------------------------------------------------\n")
      cat("Genotypic confidence index (all environments)\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$ANN$general, digits = digits, row.names = FALSE)
      cat("---------------------------------------------------------------------------\n")
      cat("Wricke's Ecovalence\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$Ecoval, digits = digits)
      cat("---------------------------------------------------------------------------\n")
      cat("Nonparametric Superiority index\n")
      cat("---------------------------------------------------------------------------\n")
      print.data.frame(var$Superiority$index, digits = digits, row.names = FALSE)
      cat("------------------------------- End of data -----------------------------\n\n\n\n")
    }
  }
}
