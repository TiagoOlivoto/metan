#' Summary a WAASB object
#'
#' Summary the \code{WAASB} object in two ways. By default, the results are
#' shown in the R console. The results can also be exported to the directory.
#'
#'
#' @param object The \code{WAASB} object
#' @param export A logical argument. If \code{TRUE|T}, a *.txt file is exported
#' to the working directory
#' @param blup A logical argument. If \code{TRUE|T}, the blups are shown.
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Other arguments of the function
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method summary WAASB
#' @export
#' @examples
#'
#' library(METAAB)
#' model = WAASB(data_ge,
#'               resp = c(GY, HM),
#'               gen = GEN,
#'               env = ENV,
#'               rep = REP)
#' summary(model)
#'
summary.WAASB <- function(object, export = FALSE, blup = FALSE, file.name = NULL,
    digits = 4, ...) {

    class <- class(object)
    if (!class == "WAASB") {
        stop("The object must be of class 'WAASB'")
    }

    if (export == TRUE) {
        if (is.null(file.name) == T) {
            file.name <- "WAASB summary"
        } else {
            file.name <- file.name
        }
        sink(paste0(file.name, ".txt"))
        options(max.print = 99999999, width = 110)
        for (i in 1:length(object)) {
            var <- object[[i]]
            cat("Variable", names(object)[i], "\n")
            cat("----------------------- Individual fixed-model analysis of variance ---------------------\n")
            print.data.frame(var$individual$individual, row.names = FALSE)
            cat("\n")
            cat("---------------------------------- Fixed effects --------------------------------------\n")
            print.data.frame(var$fixed)
            cat("\n")
            cat("--------------------------------- Random effects --------------------------------------\n")
            print.data.frame(var$random, row.names = FALSE)
            cat("\n")
            cat("------------------------------ Likelihood ratio test ----------------------------------\n")
            print.data.frame(var$LRT)
            cat("\n")
            cat("-----------------Variance components and genetic parameters --------------------------\n")
            print.data.frame(var$ESTIMATES, row.names = FALSE)
            cat("\n")
            cat("------------- Principal component analysis of the G x E interaction matrix -------------\n")
            print.data.frame(var$PCA, row.names = FALSE)
            cat("\n")
            if (blup == TRUE) {
                cat("---------------------------------- BLUPs for genotypes ---------------------------------\n")
                print.data.frame(var$BLUPgen, row.names = FALSE)
                cat("\n")
                cat("--------------------------- BLUPs for genotypes-vs-environments ------------------------\n")
                print.data.frame(var$BLUPgge, row.names = FALSE)
                cat("\n")
            }
            cat("------------------------- Some information regarding the analysis ----------------------\n")
            print.data.frame(var$Details, row.names = FALSE)
            cat("-------------------------------------- End of data -------------------------------------\n\n\n\n")
        }
        sink()
    } else {
        for (i in 1:length(object)) {
            var <- object[[i]]
            options(digits = digits)
            cat("Variable", names(object)[i], "\n")
            cat("-----------------------------------------------------------------------------------------\n")
            cat("----------------------- Individual fixed-model analysis of variance ---------------------\n")
            print.data.frame(var$individual$individual, row.names = FALSE)
            cat("\n")
            cat("---------------------------------- Fixed effects --------------------------------------\n")
            print.data.frame(var$fixed)
            cat("\n")
            cat("--------------------------------- Random effects --------------------------------------\n")
            print.data.frame(var$random, row.names = FALSE)
            cat("\n")
            cat("------------------------------ Likelihood ratio test ----------------------------------\n")
            print.data.frame(var$LRT)
            cat("\n")
            cat(" -----------------Variance components and genetic parameters --------------------------\n")
            print.data.frame(var$ESTIMATES, row.names = FALSE)
            cat("\n")
            cat("------------- Principal component analysis of the G x E interaction matrix -------------\n")
            print.data.frame(var$PCA, row.names = FALSE)
            cat("\n")
            if (blup == TRUE) {
                cat("---------------------------------- BLUPs for genotypes ---------------------------------\n")
                print.data.frame(var$BLUPgen, row.names = FALSE)
                cat("\n")
                cat("--------------------------- BLUPs for genotypes-vs-environments ------------------------\n")
                print.data.frame(var$BLUPgge, row.names = FALSE)
                cat("\n")
            }
            cat("------------------------- Some information regarding the analysis ----------------------\n")
            print.data.frame(var$Details, row.names = FALSE)
            cat("-------------------------------------- End of data -------------------------------------\n")

        }
    }
}



