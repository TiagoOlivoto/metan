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



