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




