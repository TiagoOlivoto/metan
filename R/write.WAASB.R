write.WAASB = function(data,
                      file.type = "xlsx",
                      file.name = "Results WAAS"){

  if(file.type == "xlsx"){

  xlsx::write.xlsx(data, file = paste0(file.name,".",file.type), sheetName = "WAAS",
                   append=F, row.names = T)

  }

  if(file.type == "txt"){

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
    print.data.frame(data$WAAS, row.names = FALSE)
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
    cat("------------------------------------------------------ End of data -----------------------------------------------------------\n")
    cat("\n")
    cat("The analysis consisted of ", data$Ngen, " genotypes and ", data$Nenv, " environments. Considering the p-value informed in the analysis, the number \n",
        "of significant Principal Components were", data$NSPC, ". Thus, the Weighted Average of Absolute scores was computed considering \n",
        "this number of PCs. The values of WAASY were estimated considering the weights for Response variable and stability\n",
        "equal to ", data$WgtResponse, "% and ", data$WgtWAAS, "%, respectively. The overall mean was ", round(data$OVmean, 4)," The minimum and maximum values were \n",
        data$Min," and ", data$Max, ", respectively. Regarding the environments, that one with the \n",
        "lower mean was the", data$MinENV,". The Environment with the largest mean was", data$MaxENV,". Regarding the \n",
        "Genotypes, that one  with the lower mean across the environmnets was the", data$MinGEN,". The Genotype with the largest  \n",
        "mean across, the tested environment was the",  data$MaxGEN, ".")
    cat("\n")
    cat("\n")
    cat("------------------------------------------------------ End of data -----------------------------------------------------------\n")
    sink()




  }



}
