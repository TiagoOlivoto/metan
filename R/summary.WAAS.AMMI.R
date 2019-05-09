#' Summary a WAAS.AMMI object
#'
#' Summary the \code{WAAS.AMMI} object in two ways. By default, the results are
#' shown in the R console. The results can also be exported to the directory.
#'
#'
#' @param object The \code{WAAS.AMMI} object
#' @param export A logical argument. If \code{TRUE|T}, a *.txt file is exported
#' to the working directory
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Other arguments of the function
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method summary WAAS.AMMI
#' @export
#' @examples
#'
#' library(METAAB)
#' model = WAAS.AMMI(data_ge,
#'                   resp = c(GY, HM),
#'                   gen = GEN,
#'                   env = ENV,
#'                   rep = REP)
#' summary(model)
#'
summary.WAAS.AMMI <- function(object, export = FALSE, file.name = NULL, digits = 4,
    ...) {

    class <- class(object)
    if (!class == "WAAS.AMMI") {
        stop("The object must be of class 'WAAS.AMMI'")
    }

    if (export == TRUE) {
        if (is.null(file.name) == T) {
            file.name <- "WAAS.AMMI Summary"
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
            cat("-------------------------------------- AMMI analysis -----------------------------------\n")
            print.data.frame(var$anova)
            cat("\n")
            cat("------------------------- Weighted average of the absolute scores ----------------------\n")
            print.data.frame(var$model, row.names = FALSE)
            cat("\n")
            cat("---------------------------- Means for genotypes-vs-environments -----------------------\n")
            print.data.frame(var$MeansGxE, row.names = FALSE)
            cat("\n")
            cat("------------------------- Some information regarding the analysis ----------------------\n")
            print.data.frame(var$Details, row.names = FALSE)
            cat("-------------------------------------- End of data -------------------------------------\n\n\n\n")

        }
        sink()
    } else {
        for (i in 1:length(object)) {
            var <- object[[i]]
            cat("Variable", names(object)[i], "\n")
            cat("----------------------------- Individual analysis of variance---------------------------\n")
            print.data.frame(var$individual$individual, row.names = FALSE)
            cat("\n")
            cat("-------------------------------------- AMMI analysis -----------------------------------\n")
            print.data.frame(var$anova)
            cat("\n")
            cat("------------------------- Weighted average of the absolute scores ----------------------\n")
            print.data.frame(var$model, row.names = FALSE)
            cat("\n")
            cat("---------------------------- Means for genotypes-vs-environments -----------------------\n")
            print.data.frame(var$MeansGxE, row.names = FALSE)
            cat("\n")
            cat("------------------------- Some information regarding the analysis ----------------------\n")
            print.data.frame(var$Details, row.names = FALSE)
            cat("------------------------------------ End of data ---------------------------------------\n")
        }
    }
}




