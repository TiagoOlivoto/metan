#' Summarise a mtsi object
#'
#' Summary a \code{mtsi} object in two ways. By default, the results are shown
#' in the R console. The results can also be exported to the directory.
#'
#'
#' @param object The \code{mtsi} object
#' @param export A logical argument. If \code{TRUE|T}, a *.txt file is exported
#' to the working directory
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Other arguments of the function
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method summary mtsi
#' @export
#' @examples
#' \dontrun{
#' library(metan)
#' # Based on stability only
#' MTSI_MODEL = waasb(data_ge,
#'                    resp = c(GY, HM),
#'                    gen = GEN,
#'                    env = ENV,
#'                    rep = REP)
#'
#' MTSI_index = mtsi(MTSI_MODEL)
#' summary(MTSI_index)
#' summary(MTSI_index,
#'         export = TRUE,
#'         file.name = "my results")
#' }
#'
summary.mtsi <- function(object, export = FALSE, file.name = NULL, digits = 4, ...) {
    class <- class(object)
    if (!class == "mtsi") {
        stop("The object must be of class 'mtsi'")
    }

    if (export == TRUE) {
        if (is.null(file.name) == T) {
            file.name <- "mtsi summary"
        } else {
            file.name <- file.name
        }
        sink(paste0(file.name, ".txt"))
        cat("-------------------- Correlation matrix used used in factor analysis -----------------\n")
        print(object$cormat, digits = digits)
        cat("\n")
        cat("---------------------------- Principal component analysis -----------------------------\n")
        print(object$PCA, digits = digits)
        cat("\n")
        cat("--------------------------------- Initial loadings -----------------------------------\n")
        print(object$initial.loadings, row.names = FALSE, digits = digits)
        cat("\n")
        cat("-------------------------- Loadings after varimax rotation ---------------------------\n")
        print(object$finish.loadings, digits = digits)
        cat("\n")
        cat("--------------------------- Scores for genotypes-ideotype -----------------------------\n")
        print(rbind(object$scores.gen, object$scores.ide), digits = digits)
        cat("\n")
        cat("---------------------------- Multitrait stability index ------------------------------\n")
        print(object$MTSI, digits = digits)
        cat("\n")
        cat("------------------------------ Selection differential ---------------------------------\n")
        print(object$selection.diferential, digits = digits)
        cat("\n")
        cat("-------------------------- Mean of Selection differential -----------------------------\n")
        print(object$selec.dif.mean, digits = digits)
        cat("\n")
        cat("-------------------------------- Selected genotypes -----------------------------------\n")
        cat(object$Selected, digits = digits)
        cat("\n")
        cat("-------------------------------------- End of data ------------------------------------\n\n\n\n")
        sink()
    } else {
        cat("-------------------- Correlation matrix used used in factor analysis -----------------\n")
        print(object$cormat, digits = digits)
        cat("\n")
        cat("---------------------------- Principal component analysis -----------------------------\n")
        print(object$PCA, digits = digits)
        cat("\n")
        cat("--------------------------------- Initial loadings -----------------------------------\n")
        print(object$initial.loadings, row.names = FALSE, digits = digits)
        cat("\n")
        cat("-------------------------- Loadings after varimax rotation ---------------------------\n")
        print(object$finish.loadings, digits = digits)
        cat("\n")
        cat("--------------------------- Scores for genotypes-ideotype -----------------------------\n")
        print(rbind(object$scores.gen, object$scores.ide), digits = digits)
        cat("\n")
        cat("---------------------------- Multitrait stability index ------------------------------\n")
        print(object$MTSI, digits = digits)
        cat("\n")
        cat("------------------------------ Selection differential ---------------------------------\n")
        print(object$selection.diferential, digits = digits)
        cat("\n")
        cat("-------------------------- Mean of Selection differential -----------------------------\n")
        print(object$selec.dif.mean, digits = digits)
        cat("\n")
        cat("-------------------------------- Selected genotypes ----------------------------------\n")
        cat(object$Selected, digits = digits)
        cat("\n")
        cat("------------------------------------ End of data -------------------------------------\n\n\n\n")

    }
}



