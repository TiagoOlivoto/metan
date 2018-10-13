summary.WAASB = function(object,
                         export = FALSE,
                         file.name = NULL,
                         digits = 4,
                         ...){

  class = class(object)
  if(!class == "WAASB"){
    stop("The object must be of class 'WAASB'")
  }

  options(digits = digits)

  cat("------------------------ Some informations related to the analysis ----------------------\n")
  cat("WAAS = Weighted average of Absolute Scores ","\n")
  cat("PctResp = Percentage of the i-th observation in relation to the maximum productivity","\n")
  cat("PctWAAS = Percentage of the i-th observation in relation to the minimum WAAS ","\n")
  cat("WAASY = Weighted average of Absolute Scores and Yield (this index allow weighting the productivity and the stability) ","\n")
  cat("----------------------- Individual fixed-model analysis of variance ---------------------\n")
  print.data.frame(object$individual$individual, row.names = FALSE)
  cat("\n")
  cat("\nFixed effects\n")
  print.data.frame(object$REML$fixed, row.names = FALSE)
  cat("\n")
  cat("\nRandom effects\n")
  print.data.frame(object$REML$random, row.names = FALSE)
  cat("\n")
  cat("\nStatistics\n")
  print.data.frame(object$REML$statistics, row.names = FALSE)
  cat("\n")
  cat("\nSome estimates\n")
  print.data.frame(object$ESTIMATES, row.names = FALSE)
  cat("\n")
  cat("------------------------- Weighted average of the absolute scores ----------------------\n")
  print.data.frame(object$model, row.names = FALSE)
  cat("\n")
  cat("------------- Principal component analysis of the G x E interaction matrix -------------\n")
  print.data.frame(object$PCA, row.names = FALSE)
  cat("\n")
  cat("---------------------------------- BLUPs for genotypes ---------------------------------\n")
  print.data.frame(object$BLUPgen, row.names = FALSE)
  cat("\n")
  cat("--------------------------- BLUPs for genotypes-vs-environments ------------------------\n")
  print.data.frame(object$BLUPgge, row.names = FALSE)
  cat("\n")
  cat("------------------------ Average of Genotypes in the environments ----------------------\n")
  print.data.frame(object$MeansGxE, row.names = FALSE)
  cat("\n")
  cat("------------------------- Some information regarding the analysis ----------------------\n")
  print.data.frame(object$Details, row.names = FALSE)
  cat("-------------------------------------- End of data -------------------------------------\n")

  if(export  ==  TRUE){
    if(is.null(file.name) == T){
      file.name = "WAASB summary"
    } else{file.name = file.name}
    sink(paste0(file.name,".txt"))
    options(max.print = 99999999, width = 90)

    cat("--------------------- Some informations related to the analysis -----------------------------------------\n")
    cat("WAAS = Weighted average of Absolute Scores ","\n")
    cat("PctResp = Percentage of the i-th observation in relation to the maximum productivity","\n")
    cat("PctWAAS = Percentage of the i-th observation in relation to the minimum WAAS ","\n")
    cat("WAASY = Weighted average of Absolute Scores and Yield (this index allow weighting the productivity and the stability) ","\n")
    cat("----------------------- Individual fixed-model analysis of variance ---------------------\n")
    print.data.frame(object$individual$individual, row.names = FALSE)
    cat("\n")
    cat("\nFixed effects\n")
    print.data.frame(object$REML$fixed, row.names = FALSE)
    cat("\n")
    cat("\nRandom effects\n")
    print.data.frame(object$REML$random, row.names = FALSE)
    cat("\n")
    cat("\nStatistics\n")
    print.data.frame(object$REML$statistics, row.names = FALSE)
    cat("\n")
    cat("\nSome estimates\n")
    print.data.frame(object$ESTIMATES, row.names = FALSE)
    cat("\n")
    cat("------------------------- Weighted average of the absolute scores ----------------------\n")
    print.data.frame(object$model, row.names = FALSE)
    cat("\n")
    cat("------------- Principal component analysis of the G x E interaction matrix -------------\n")
    print.data.frame(object$PCA, row.names = FALSE)
    cat("\n")
    cat("---------------------------------- BLUPs for genotypes ---------------------------------\n")
    print.data.frame(object$BLUPgen, row.names = FALSE)
    cat("\n")
    cat("--------------------------- BLUPs for genotypes-vs-environments ------------------------\n")
    print.data.frame(object$BLUPgge, row.names = FALSE)
    cat("\n")
    cat("------------------------ Average of Genotypes in the environments ----------------------\n")
    print.data.frame(object$MeansGxE, row.names = FALSE)
    cat("\n")
    cat("------------------------- Some information regarding the analysis ----------------------\n")
    print.data.frame(object$Details, row.names = FALSE)
    cat("-------------------------------------- End of data -------------------------------------\n")
    sink()
  }

}




