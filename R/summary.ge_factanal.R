#' Summary a ge_factanal object
#'
#' Summary the \code{ge_factanal} object in two ways. By default, the results are
#' shown in the R console. The results can also be exported to the directory.
#'
#'
#' @param object The \code{ge_factanal} object
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported
#' to the working directory
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Other arguments of the function
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method summary ge_factanal
#' @export
#' @examples
#'
#' model = ge_factanal(data_ge2,
#'                     env = ENV,
#'                     gen = GEN,
#'                     rep = REP,
#'                     resp = PH)
#' summary(model)
#'
summary.ge_factanal <- function(object, export = FALSE, file.name = NULL, digits = 4,
    ...) {

    class <- class(object)
    if (!class == "ge_factanal") {
        stop("The object must be of class 'ge_factanal'")
    }
    if (export == TRUE) {
        if (is.null(file.name) == T) {
            file.name <- "ge_factanal Summary"
        } else {
            file.name <- file.name
        }
        sink(paste0(file.name, ".txt"))
        for (i in 1:length(object)) {
            var <- object[[i]]
            cat("Variable", names(object)[i], "\n")
            cat("------------------------------------------------------------------------------------\n")
            cat("Correlation matrix among environments\n")
            cat("------------------------------------------------------------------------------------\n")
            print(var$cormat, digits = digits)
            cat("------------------------------------------------------------------------------------\n")
            cat("Eigenvalues and explained variance\n")
            cat("------------------------------------------------------------------------------------\n")
            print(var$PCA, digits = digits, row.names = FALSE)
            cat("------------------------------------------------------------------------------------\n")
            cat("Initial loadings\n")
            cat("------------------------------------------------------------------------------------\n")
            print(var$initial.loadings, digits = digits, row.names = FALSE)
            cat("------------------------------------------------------------------------------------\n")
            cat("Loadings after varimax rotation and commonalities\n")
            cat("------------------------------------------------------------------------------------\n")
            print(var$FA, digits = digits, row.names = FALSE)
            cat("------------------------------------------------------------------------------------\n")
            cat("Environmental stratification based on factor\n")
            cat("------------------------------------------------------------------------------------\n")
            print(var$env_strat, digits = digits, row.names = FALSE)
            cat("------------------------------------------------------------------------------------\n")
            cat("Mean = mean; Min = minimum; Max = maximum; CV = coefficient of variation (%)\n")
            cat("The summary statistics are based on the men values of ", length(unique(var$data$REP)), "replicates\n")
            cat("------------------------------------------------------------------------------------\n")
        }
        sink()
    } else {
        for (i in 1:length(object)) {
            var <- object[[i]]
            cat("Variable", names(object)[i], "\n")
            cat("------------------------------------------------------------------------------------\n")
            cat("Correlation matrix among environments\n")
            cat("------------------------------------------------------------------------------------\n")
            print(var$cormat, digits = digits)
            cat("------------------------------------------------------------------------------------\n")
            cat("Eigenvalues and explained variance\n")
            cat("------------------------------------------------------------------------------------\n")
            print(var$PCA, digits = digits, row.names = FALSE)
            cat("------------------------------------------------------------------------------------\n")
            cat("Initial loadings\n")
            cat("------------------------------------------------------------------------------------\n")
            print(var$initial.loadings, digits = digits, row.names = FALSE)
            cat("------------------------------------------------------------------------------------\n")
            cat("Loadings after varimax rotation and commonalities\n")
            cat("------------------------------------------------------------------------------------\n")
            print(var$FA, digits = digits, row.names = FALSE)
            cat("------------------------------------------------------------------------------------\n")
            cat("Environmental stratification based on factor\n")
            cat("------------------------------------------------------------------------------------\n")
            print(var$env_strat, digits = digits, row.names = FALSE)
            cat("------------------------------------------------------------------------------------\n")
            cat("Mean = mean; Min = minimum; Max = maximum; CV = coefficient of variation (%)\n")
            cat("The summary statistics are based on the men values of ", length(unique(var$data$REP)), "replicates\n")
            cat("------------------------------------------------------------------------------------\n")
        }
    }
}
