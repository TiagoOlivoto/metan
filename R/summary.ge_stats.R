summary.ge_stats <- function(object, export = FALSE, file.name = NULL, digits = 4,
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
    options(max.print = 99999999, width = 90, digits = digits)
    for (i in 1:length(object)) {
      var <- object[[i]]
      cat("Variable", names(object)[i], "\n")
      cat("------------------------------ Individual analysis of variance--- ---------------------\n")
      print.data.frame(var$individual$individual, row.names = FALSE)
      cat("\n")
      cat("------------------------------- Genotype-vs-environment mean --------------------------\n")
      print.data.frame(var$ge_mean)
      cat("\n")
      cat("------------------------------ Genotype-vs-environment effects -------------------------\n")
      print.data.frame(var$ge_effect)
      cat("\n")
      cat("------------------------ Genotype + Genotype-vs-environment effects --------------------\n")
      print.data.frame(var$gge_effect)
      cat("\n")
      cat("----------------------------- Genotype-vs-environment statistics -----------------------\n")
      print.data.frame(var$ge_stats)
      cat("\n")
      cat("----------------------------- Regression-based stability (anova) -----------------------\n")
      print.data.frame(var$ER$anova)
      cat("\n")
      cat("--------------------------- Regression-based stability (parameters) ---------------------\n")
      print.data.frame(var$ER$regression)
      cat("\n")
      cat("-------------------------------------- End of data -------------------------------------\n\n\n\n")

    }
    sink()
  } else {
    for (i in 1:length(object)) {
      var <- object[[i]]
      cat("Variable", names(object)[i], "\n")
      cat("------------------------------ Individual analysis of variance--- ---------------------\n")
      print.data.frame(var$individual$individual, row.names = FALSE)
      cat("\n")
      cat("------------------------------- Genotype-vs-environment mean --------------------------\n")
      print.data.frame(var$ge_mean)
      cat("\n")
      cat("------------------------------ Genotype-vs-environment effects -------------------------\n")
      print.data.frame(var$ge_effect)
      cat("\n")
      cat("------------------------ Genotype + Genotype-vs-environment effects --------------------\n")
      print.data.frame(var$gge_effect)
      cat("\n")
      cat("----------------------------- Genotype-vs-environment statistics -----------------------\n")
      print.data.frame(var$ge_stats)
      cat("\n")
      cat("----------------------------- Regression-based stability (anova) -----------------------\n")
      print.data.frame(var$ER$anova)
      cat("\n")
      cat("--------------------------- Regression-based stability (parameters) ---------------------\n")
      print.data.frame(var$ER$regression)
      cat("\n")
      cat("-------------------------------------- End of data -------------------------------------\n\n\n\n")
    }
  }
}




