#' Summarise a waas object
#'
#' Summary the \code{waas} object in two ways. By default, the results are
#' shown in the R console. The results can also be exported to the directory.
#'
#'
#' @param object The \code{waas} object
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported
#' to the working directory
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Other arguments of the function
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method summary waas
#' @export
#' @examples
#'
#' library(metan)
#' model = waas(data_ge,
#'              resp = c(GY, HM),
#'              gen = GEN,
#'              env = ENV,
#'              rep = REP)
#' summary(model)
#'
summary.waas <- function(object, export = FALSE, file.name = NULL, digits = 4,
    ...) {

    class <- class(object)
    if (!class == "waas") {
        stop("The object must be of class 'waas'")
    }

    if (export == TRUE) {
        if (is.null(file.name) == TRUE) {
            file.name <- "waas summary"
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
            cat("AMMI analysis table\n")
            cat("---------------------------------------------------------------------------\n")
            print.data.frame(var$anova, digits = digits)
            cat("---------------------------------------------------------------------------\n")
            cat("Weighted average of the absolute scores\n")
            cat("---------------------------------------------------------------------------\n")
            print.data.frame(var$model, row.names = FALSE, digits = digits)
            cat("---------------------------------------------------------------------------\n")
            cat("Some information regarding the analysis\n")
            cat("---------------------------------------------------------------------------\n")
            print.data.frame(var$Details, row.names = FALSE, digits = digits)
            cat("------------------------------- End of data -------------------------------\n\n\n\n")

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
            cat("AMMI analysis table\n")
            cat("---------------------------------------------------------------------------\n")
            print.data.frame(var$anova, digits = digits)
            cat("---------------------------------------------------------------------------\n")
            cat("Weighted average of the absolute scores\n")
            cat("---------------------------------------------------------------------------\n")
            print.data.frame(var$model, row.names = FALSE, digits = digits)
            cat("---------------------------------------------------------------------------\n")
            cat("Some information regarding the analysis\n")
            cat("---------------------------------------------------------------------------\n")
            print.data.frame(var$Details, row.names = FALSE, digits = digits)
            cat("------------------------------- End of data -------------------------------\n\n\n\n")
        }
    }
}




