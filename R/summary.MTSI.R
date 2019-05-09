#' Summary a MTSI object
#'
#' Summary a \code{MTSI} object in two ways. By default, the results are shown
#' in the R console. The results can also be exported to the directory.
#'
#'
#' @param object The \code{MTSI} object
#' @param export A logical argument. If \code{TRUE|T}, a *.txt file is exported
#' to the working directory
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Other arguments of the function
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method summary MTSI
#' @export
#' @examples
#' \dontrun{
#' library(METAAB)
#' # Based on stability only
#' MTSI_MODEL = WAASB(data_ge,
#'                    resp = c(GY, HM),
#'                    gen = GEN,
#'                    env = ENV,
#'                    rep = REP)
#'
#' MTSI_index = MTSI(MTSI_MODEL)
#' summary(MTSI_index)
#' summary(MTSI_index,
#'         export = TRUE,
#'         file.name = "my results")
#' }
#'
summary.MTSI <- function(object, export = FALSE, file.name = NULL, digits = 4, ...) {
    class <- class(object)
    if (!class == "MTSI") {
        stop("The object must be of class 'MTSI'")
    }

    if (export == TRUE) {
        if (is.null(file.name) == T) {
            file.name <- "MTSI summary"
        } else {
            file.name <- file.name
        }
        sink(paste0(file.name, ".txt"))
        options(max.print = 99999999, width = 110)

        cat("-------------------- Correlation matrix used used in factor analysis -----------------\n")
        print(object$cormat)
        cat("\n")
        cat("---------------------------- Principal component analysis -----------------------------\n")
        print(object$PCA)
        cat("\n")
        cat("--------------------------------- Initial loadings -----------------------------------\n")
        print(object$initial.loadings, row.names = FALSE)
        cat("\n")
        cat("-------------------------- Loadings after varimax rotation ---------------------------\n")
        print(object$finish.loadings)
        cat("\n")
        cat("--------------------------- Scores for genotypes-ideotype -----------------------------\n")
        print(rbind(object$scores.gen, object$scores.ide))
        cat("\n")
        cat("---------------------------- Multitrait stability index ------------------------------\n")
        print(object$MTSI)
        cat("\n")
        cat("------------------------------ Selection differential ---------------------------------\n")
        print(object$selection.diferential)
        cat("\n")
        cat("-------------------------- Mean of Selection differential -----------------------------\n")
        print(object$selec.dif.mean)
        cat("\n")
        cat("-------------------------------- Selected genotypes -----------------------------------\n")
        cat(object$Selected)
        cat("\n")
        cat("-------------------------------------- End of data ------------------------------------\n\n\n\n")
        sink()
    } else {
        options(digits = digits)
        cat("-------------------- Correlation matrix used used in factor analysis -----------------\n")
        print(object$cormat)
        cat("\n")
        cat("---------------------------- Principal component analysis -----------------------------\n")
        print(object$PCA)
        cat("\n")
        cat("--------------------------------- Initial loadings -----------------------------------\n")
        print(object$initial.loadings, row.names = FALSE)
        cat("\n")
        cat("-------------------------- Loadings after varimax rotation ---------------------------\n")
        print(object$finish.loadings)
        cat("\n")
        cat("--------------------------- Scores for genotypes-ideotype -----------------------------\n")
        print(rbind(object$scores.gen, object$scores.ide))
        cat("\n")
        cat("---------------------------- Multitrait stability index ------------------------------\n")
        print(object$MTSI)
        cat("\n")
        cat("------------------------------ Selection differential ---------------------------------\n")
        print(object$selection.diferential)
        cat("\n")
        cat("-------------------------- Mean of Selection differential -----------------------------\n")
        print(object$selec.dif.mean)
        cat("\n")
        cat("-------------------------------- Selected genotypes ----------------------------------\n")
        cat(object$Selected)
        cat("\n")
        cat("------------------------------------ End of data -------------------------------------\n\n\n\n")

    }
}



