summary.WAASB = function(data, export = FALSE, file.name = "WAASB"){

  cat("----------------------------------------- Some informations related to the analysis -----------------------------------------\n")
  cat("WAAS = Weighted average of Absolute Scores ","\n")
  cat("PctResp = Percentage of the i-th observation in relation to the maximum productivity","\n")
  cat("PctWAAS = Percentage of the i-th observation in relation to the minimum WAAS ","\n")
  cat("WAASY = Weighted average of Absolute Scores and Yield (this index allow weighting the productivity and the stability) ","\n")
  cat("-----------------------------------------------------------------------------------------------------------------------------\n")
  cat("------------------------------------------------- Summary of analysis of variance--- ----------------------------------------\n")
  cat("\nFixed effects")
  print.data.frame(data$REML$fixed, row.names = FALSE)
  cat("\n")
  cat("\nRandom effects")

  print.data.frame(data$REML$random, row.names = FALSE)
  cat("\n")
  cat("\nStatistics")

  print.data.frame(data$REML$statistics, row.names = FALSE)
  cat("\n")
  cat("--------------------------------------------- Weighted average of the absolute scores ----------------------------------------\n")
  print.data.frame(data$WAASB, row.names = FALSE)
  cat("\n")
  cat("---------------------------------- Principal component analysis of the G x E interaction matrix-------------------------------\n")
  print.data.frame(data$PCA, row.names = FALSE)
  cat("\n")
  cat("------------------------------------------------------- BLUPs for genotypes ---------------------------------------------------\n")
  print.data.frame(data$BLUPgen, row.names = FALSE)
  cat("\n")
  cat("----------------------------------------------- BLUPs for genotypes-vs-environments-------------------------------------------\n")
  print.data.frame(data$BLUPgge, row.names = FALSE)
  cat("\n")
  cat("-------------------------------------------- Average of Genotypes in the environments -----------------------------------------\n")
  print.data.frame(data$MeansGxE, row.names = FALSE)
  cat("\n")
  cat("\n")
  cat("The analysis consisted of ", data$Ngen, " genotypes and ", data$Nenv, " environments. The Weighted Average of Absolute scores was computed considering", length(data$PCA$PC),
      "PCs. The values of WAASY were estimated considering the weights for Response variable and stability equal to ", data$WgtResponse, "% and ", data$WgtWAAS, "%, respectively. The overall mean was ", round(data$OVmean, 4),
      ". The minimum and maximum values were", data$Min," and ", data$Max, ", respectively. Regarding the environments, that one with the lower mean was the", data$MinENV,
      ". The Environment with the largest mean was", data$MaxENV,". Regarding the Genotypes, that one  with the lower mean across the environmnets was the", data$MinGEN,
      ". The Genotype with the largest mean across, the tested environment was the",  data$MaxGEN, ".")
  cat("\n")
  cat("\n")
  cat("------------------------------------------------------ End of data -----------------------------------------------------------\n")

  if(export == TRUE|T){

    sink(paste0(file.name,".txt"))
    options(max.print=99999999)

    cat("----------------------------------------- Some informations related to the analysis -----------------------------------------\n")
    cat("WAAS = Weighted average of Absolute Scores ","\n")
    cat("PctResp = Percentage of the i-th observation in relation to the maximum productivity","\n")
    cat("PctWAAS = Percentage of the i-th observation in relation to the minimum WAAS ","\n")
    cat("WAASY = Weighted average of Absolute Scores and Yield (this index allow weighting the productivity and the stability) ","\n")
    cat("-----------------------------------------------------------------------------------------------------------------------------\n")
    cat("------------------------------------------------- Summary of analysis of variance--- ----------------------------------------\n")
    cat("\nFixed effects")
    print.data.frame(data$REML$fixed, row.names = FALSE)
    cat("\n")
    cat("\nRandom effects")

    print.data.frame(data$REML$random, row.names = FALSE)
    cat("\n")
    cat("\nStatistics")

    print.data.frame(data$REML$statistics, row.names = FALSE)
    cat("\n")
    cat("--------------------------------------------- Weighted average of the absolute scores ----------------------------------------\n")
    print.data.frame(data$WAASB, row.names = FALSE)
    cat("\n")
    cat("---------------------------------- Principal component analysis of the G x E interaction matrix-------------------------------\n")
    print.data.frame(data$PCA, row.names = FALSE)
    cat("\n")
    cat("------------------------------------------------------- BLUPs for genotypes ---------------------------------------------------\n")
    print.data.frame(data$BLUPgen, row.names = FALSE)
    cat("\n")
    cat("----------------------------------------------- BLUPs for genotypes-vs-environments-------------------------------------------\n")
    print.data.frame(data$BLUPgge, row.names = FALSE)
    cat("\n")
    cat("-------------------------------------------- Average of Genotypes in the environments -----------------------------------------\n")
    print.data.frame(data$MeansGxE, row.names = FALSE)
    cat("\n")
    cat("\n")
    cat("The analysis consisted of ", data$Ngen, " genotypes and ", data$Nenv, " environments. The Weighted Average of Absolute scores was computed considering", length(data$PCA$PC),
        "PCs. The values of WAASY were estimated considering the weights for Response variable and stability equal to ", data$WgtResponse, "% and ", data$WgtWAAS, "%, respectively. The overall mean was ", round(data$OVmean, 4),
        ". The minimum and maximum values were", data$Min," and ", data$Max, ", respectively. Regarding the environments, that one with the lower mean was the", data$MinENV,
        ". The Environment with the largest mean was", data$MaxENV,". Regarding the Genotypes, that one  with the lower mean across the environmnets was the", data$MinGEN,
        ". The Genotype with the largest mean across, the tested environment was the",  data$MaxGEN, ".")
    cat("\n")
    cat("\n")
    cat("------------------------------------------------------ End of data -----------------------------------------------------------\n")
    sink()
  }
}
