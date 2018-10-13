summary.WAAS.AMMI = function(object,
                         export = FALSE,
                         file.name = NULL,
                         digits = 4,
                         ...){

  class = class(object)
  if(!class == "WAAS.AMMI"){
    stop("The object must be of class 'WAAS.AMMI'")
  }

  options(digits = digits)

  cat("------------------------ Some informations related to the analysis ---------------------\n")
  cat("WAAS = Weighted average of Absolute Scores ","\n")
  cat("PctResp = Percentage of the i-th observation in relation to the maximum productivity","\n")
  cat("PctWAAS = Percentage of the i-th observation in relation to the minimum WAAS ","\n")
  cat("WAASY = Weighted average of Absolute Scores and Yield (this index allow weighting the productivity and the stability) ","\n")
  cat("------------------------------ Individual analysis of variance--- ---------------------\n")
  print.data.frame(object$individual$individual, row.names = FALSE)
  cat("\n")
  cat("-------------------------------------- AMMI analysis -----------------------------------\n")
  print.data.frame(object$anova)
  cat("\n")
  cat("------------------------- Weighted average of the absolute scores ----------------------\n")
  print.data.frame(object$model, row.names = FALSE)
  cat("\n")
  cat("---------------------------- Means for genotypes-vs-environments -----------------------\n")
  print.data.frame(object$MeansGxE, row.names = FALSE)
  cat("\n")
  cat("------------------------- Some information regarding the analysis ----------------------\n")
  print.data.frame(object$Details, row.names = FALSE)
  cat("-------------------------------------- End of data -------------------------------------\n")

  if(export  ==  TRUE){
    if(is.null(file.name)==T){
      file.name = "WAAS.AMMI Summary"
    } else{file.name = file.name}
    sink(paste0(file.name,".txt"))
    options(max.print = 99999999, width = 90)

    cat("------------------------ Some informations related to the analysis ---------------------\n")
    cat("WAAS = Weighted average of Absolute Scores ","\n")
    cat("PctResp = Percentage of the i-th observation in relation to the maximum productivity","\n")
    cat("PctWAAS = Percentage of the i-th observation in relation to the minimum WAAS ","\n")
    cat("WAASY = Weighted average of Absolute Scores and Yield (this index allow weighting the productivity and the stability) ","\n")
    cat("----------------------------- Individual analysis of variance---------------------------\n")
    print.data.frame(object$individual$individual, row.names = FALSE)
    cat("\n")
    cat("-------------------------------------- AMMI analysis -----------------------------------\n")
    print.data.frame(object$anova)
    cat("\n")
    cat("------------------------- Weighted average of the absolute scores ----------------------\n")
    print.data.frame(object$model, row.names = FALSE)
    cat("\n")
    cat("---------------------------- Means for genotypes-vs-environments -----------------------\n")
    print.data.frame(object$MeansGxE, row.names = FALSE)
    cat("\n")
    cat("------------------------- Some information regarding the analysis ----------------------\n")
    print.data.frame(object$Details, row.names = FALSE)
    cat("------------------------------------ End of data ---------------------------------------\n")
    sink()
  }
}




