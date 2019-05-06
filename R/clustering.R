clustering = function(.data,
                      ... = NULL,
                      means_by = NULL,
                      scale = FALSE,
                      selvar = FALSE,
                      verbose = TRUE,
                      distmethod = "euclidean",
                      clustmethod = "average",
                      nclust = NULL){
  if (scale == TRUE && selvar == TRUE){
    stop("It is not possible to execute the algorithm for variable selection when 'scale = TRUE'. Please, verify.")
  }
  if (!distmethod %in% c("euclidean", "maximum", "manhattan", "canberra",
                         "binary", "minkowski", "pearson", "spearman", "kendall")){
    stop("The argument 'distmethod' is incorrect. It should be one of the 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski', 'pearson', 'spearman', or 'kendall'.")
  }
  if (!clustmethod %in% c("complete", "ward.D", "ward.D2", "single",
                          "average", "mcquitty", "median", "centroid")){
    stop("The argument 'distmethod' is incorrect. It should be one of the 'ward.D', 'ward.D2', 'single', 'average', 'mcquitty', 'median' or 'centroid'.")
  }

if(any(class(.data) == "group_factors")){
    dfs = list()
    datain = .data
for (k in 1:length(.data)){
  .data = datain[[k]]
  nam  = names(datain[k])
  if (!missing(...) && !missing(means_by)){
    data = suppressWarnings(dplyr::select(.data, !!dplyr::enquo(means_by), ...) %>%
                              dplyr::group_by(!!dplyr::enquo(means_by)) %>%
                              dplyr::summarise_all(mean) %>%
                              as.data.frame())
    rownames(data) = data[,1]
    data[,1] = NULL
    data = data %>% select_if(function(x) any(!is.na(x)))
  }
  if (!missing(...) && missing(means_by)){
    data = dplyr::select(.data, ...)
  }
  if (missing(...) && !missing(means_by)){
    data = suppressWarnings(dplyr::group_by(.data, !!dplyr::enquo(means_by)) %>%
                              dplyr::summarise_all(mean) %>%
                              as.data.frame())
    rownames(data) = data[,1]
    data[,1] = NULL
    data = data %>% select_if(function(x) any(!is.na(x)))
  }
  if (missing(...) && missing(means_by)){
    data = .data[ , unlist(lapply(.data, is.numeric))]
    if (verbose == TRUE){
      if (sum(lapply(.data, is.factor)==TRUE)>0){
        message("The columns ", paste0(collapse = " ", names(.data[ , unlist(lapply(.data, is.factor)) ])),
                " where excluded. Use 'means_by' to compute the distances using the means of a factor. If you want to compute the distances for each level of a factor, use the function 'group_factors() before.' ")
      }
    }
  }
  if (scale == TRUE){
    data = data.frame(scale(data, center = FALSE, scale = apply(data, 2, sd, na.rm = TRUE)))
  } else {
    data = data
    }
  if (selvar == TRUE){
    n = (ncol(data)-1)
    statistics = data.frame(matrix(nrow = n, ncol = 6))
    ModelEstimates = list()
    modelcode = 1
    namesv = "-"
    original = data
    if (distmethod %in% c("pearson", "spearman", "kendall")){
      dein = as.dist(cor(t(data), method = distmethod))
    } else{
      dein = dist(data, method = distmethod, diag = T, upper = T)
    }
    if(verbose == TRUE){
    cat("\n\n\n\n----------------------------------------------------------------------------\n")
    cat("Level", nam, "\n")
    cat("----------------------------------------------------------------------------\n")
    }
    for (i in 1:n){
      if (distmethod %in% c("pearson", "spearman", "kendall")){
        de = as.dist(cor(t(data), method = distmethod))
      } else{
        de = dist(data, method = distmethod, diag = T, upper = T)
      }
      hc = hclust(de, method = clustmethod)
      d2 = cophenetic(hc)
      cof = cor(d2, de)
      mant = ade4::mantel.rtest(de, dein, nrepet = 1000)
      mantc = mant$obs
      mantp = mant$pvalue
      evect = data.frame(t(prcomp(data)$rotation))
      var = abs(evect)[nrow(evect),]
      names = apply(var, 1, function(x) which(x == max(x)))
      npred = ncol(data)
      statistics[i,1] = paste("Model",modelcode)
      statistics[i,2] = namesv
      statistics[i,3] = cof
      statistics[i,4] = npred
      statistics[i,5] = mantc
      statistics[i,6] = mantp
      mat = as.matrix(de)
      mat = as.data.frame(mat)
      Results = list(nvars = npred,
                     excluded = namesv,
                     namevars = names(data),
                     distance = mat,
                     cormantel = mantc,
                     pvmant = mantp)
      namesv = names(data[names])
      data2 = data.frame(data[-(match(c(namesv), names(data)))])
      data = data2
      ModelEstimates[[paste("Model",modelcode)]] = Results
      names(statistics) = c("Model", "excluded", "cophenetic", "remaining", "cormantel", "pvmantel")
      if(verbose == TRUE){
        cat(paste("Calculating model ",modelcode, " with ", npred,
                  " variables. ",namesv, " excluded in this step (",
                  round(modelcode/n*100,1),"%).\n", sep = ""))
      }
      modelcode = modelcode + 1
    }
    cat("Done!","\n")
    cat("--------------------------------------------------------------------------","\n")
    cat("\nSummary of the adjusted models","\n")
    cat("--------------------------------------------------------------------------","\n")
    print(statistics, row.names = F)
    cat("--------------------------------------------------------------------------")
    cofgrap = ggplot2::ggplot(statistics, ggplot2::aes(x = remaining, y = cophenetic))+
      ggplot2::geom_point(size = 3)+
      ggplot2::theme_bw()+
      ggplot2::geom_line(size = 1)+
      ggplot2::theme(axis.ticks.length = unit(.2, "cm"),
                     axis.text = ggplot2::element_text(size = 12, colour = "black"),
                     axis.title = ggplot2::element_text(size = 12, colour = "black"),
                     axis.ticks = ggplot2::element_line(colour = "black"),
                     plot.margin = margin(0.5, 0.5, 0.2, 0.6, "cm"),
                     axis.title.y = ggplot2::element_text(margin = margin(r=16)),
                     legend.title = ggplot2::element_blank(),
                     legend.text = ggplot2::element_text(size=12),
                     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
                     panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank(),
                     panel.grid.minor.y = ggplot2::element_blank())+
      ggplot2::labs(x = "Number of variables", y = "Cophenetic correlation")
    model = statistics$Model[which.max(statistics$cophenetic)]
    predvar = ModelEstimates[[model]]$namevars
    data = data.frame(original[(match(c(predvar), names(original)))])
    if (verbose == TRUE){
    cat("\nSuggested variables to be used in the analysis","\n")
    cat("--------------------------------------------------------------------------","\n")
    cat("The clustering was calculated with the ", model,
        "\nThe variables included in this model were...\n",
        predvar,"\n")
    cat("--------------------------------------------------------------------------")
    }
  } else {data = data
          cofgrap = NULL
          statistics = NULL
         }
  if (distmethod %in% c("pearson", "spearman", "kendall")){
    de = as.dist(1- cor(t(data), method = distmethod))
  } else {
    de = dist(data, method = distmethod, diag = T, upper = T)
  }
  mat = as.matrix(de)
  hc = hclust(de, method = clustmethod)
  d2 = cophenetic(hc)
  cof = cor(d2, de)
  k = 1.25
  pcorte = mean(hc$height) + k * sd(hc$height)
  if (!missing(nclust)){
    groups <- cutree(hc, k = nclust)
    Mgroups <- cbind(data, groups)
  distance <- hc$height[(length(hc$height) - nclust):length(hc$height)]
  Sim <- (1 - distance/max(de))
  Passos <- 1:length(Sim)
  Simgroups <- length(Sim):1
  similarity <- Sim * 100
  Tab <- cbind(Passos, Simgroups, round(similarity, 3), round(distance, 2))
  colnames(Tab) <- c("Steps", "Groups", "Similarity", "Distance")
    TabResgroups <- NULL
    MGr <- cbind(data, groups)
    for (i in 1:nclust) {
      NewGroups <- subset(MGr, groups == i)
      GrupCalc <- NewGroups[, 1:(ncol(NewGroups) - 1)]
      Qtd.Elementos <- nrow(NewGroups)
      if (Qtd.Elementos == 1)
        Media <- GrupCalc
      else Media <- apply(GrupCalc, 2, mean)
      if (Qtd.Elementos == 1)
        SqG <- 0
      else SqG <- sum(sweep(GrupCalc, 2, Media)^2)
      TabResgroups <- rbind(TabResgroups, c(i, Qtd.Elementos,
                                            SqG, Media))
    }
    colnames(TabResgroups) <- c("Cluster", "Number of Elements",
                                "Sum_sq", paste(colnames(TabResgroups[, 4:(ncol(TabResgroups))])))
  } else{
    TabResgroups = NULL
    Tab = NULL
  }
  Sqt <- sum(sweep(data, 2, apply(data, 2, mean))^2)
  temp = structure(list(data = data,
              cutpoint = pcorte,
              distance = mat,
              de = de,
              hc = hc,
              cophenetic = cof,
              Sqt = Sqt,
              tab = as.data.frame(Tab),
              clusters = as.data.frame(TabResgroups),
              cofgrap = cofgrap,
              statistics = statistics),
              class = "clustering")
  dfs[[paste(nam)]] = temp
}
return(structure(dfs, class = "group_clustering"))
} else {
  if (!missing(...) && !missing(means_by)){
    data = suppressWarnings(dplyr::select(.data, !!dplyr::enquo(means_by), ...) %>%
                              dplyr::group_by(!!dplyr::enquo(means_by)) %>%
                              dplyr::summarise_all(mean) %>%
                              as.data.frame())
    rownames(data) = data[,1]
    data[,1] = NULL
    data = data %>% select_if(function(x) any(!is.na(x)))
  }
  if (!missing(...) && missing(means_by)){
    data = dplyr::select(.data, ...)
  }
  if (missing(...) && !missing(means_by)){
    data = suppressWarnings(dplyr::group_by(.data, !!dplyr::enquo(means_by)) %>%
                              dplyr::summarise_all(mean) %>%
                              as.data.frame())
    rownames(data) = data[,1]
    data[,1] = NULL
    data = data %>% select_if(function(x) any(!is.na(x)))
  }
  if (missing(...) && missing(means_by)){
    data = .data[ , unlist(lapply(.data, is.numeric))]
    if (verbose == TRUE){
      if (sum(lapply(.data, is.factor)==TRUE)>0){
        message("The columns ", paste0(collapse = " ", names(.data[ , unlist(lapply(.data, is.factor)) ])),
                " where excluded. Use 'means_by' to compute the distances using the means of a factor. If you want to compute the distances for each level of a factor, use the function 'group_factors() before.' ")
      }
    }
  }
  if (scale == TRUE){
    data = data.frame(scale(data, center = FALSE, scale = apply(data, 2, sd, na.rm = TRUE)))
  } else{data = data}

  if (selvar == TRUE){
    n = (ncol(data)-1)
    statistics = data.frame(matrix(nrow = n, ncol = 6))
    ModelEstimates = list()
    modelcode = 1
    namesv = "-"
    original = data
    if (distmethod %in% c("pearson", "spearman", "kendall")){
      dein = as.dist(cor(t(data), method = distmethod))
    } else{
      dein = dist(data, method = distmethod, diag = T, upper = T)
    }
    for (i in 1:n){
      if (distmethod %in% c("pearson", "spearman", "kendall")){
        de = as.dist(cor(t(data), method = distmethod))
      } else{
        de = dist(data, method = distmethod, diag = T, upper = T)
      }
      hc = hclust(de, method = clustmethod)
      d2 = cophenetic(hc)
      cof = cor(d2, de)
      mant = ade4::mantel.rtest(de, dein, nrepet = 1000)
      mantc = mant$obs
      mantp = mant$pvalue
      evect = data.frame(t(prcomp(data)$rotation))
      var = abs(evect)[nrow(evect),]
      names = apply(var, 1, function(x) which(x == max(x)))
      npred = ncol(data)
      statistics[i,1] = paste("Model",modelcode)
      statistics[i,2] = namesv
      statistics[i,3] = cof
      statistics[i,4] = npred
      statistics[i,5] = mantc
      statistics[i,6] = mantp
      mat = as.matrix(de)
      mat = as.data.frame(mat)
      Results = list(nvars = npred,
                     excluded = namesv,
                     namevars = names(data),
                     distance = mat,
                     cormantel = mantc,
                     pvmant = mantp)
      namesv = names(data[names])
      data2 = data.frame(data[-(match(c(namesv), names(data)))])
      data = data2
      ModelEstimates[[paste("Model",modelcode)]] = Results
      names(statistics) = c("Model", "excluded", "cophenetic", "remaining", "cormantel", "pvmantel")
      if(verbose == TRUE){
        cat(paste("Calculating model ",modelcode, " with ", npred,
                  " variables. ",namesv, " excluded in this step (",
                  round(modelcode/n*100,1),"%).\n", sep = ""))
      }
      modelcode = modelcode + 1
    }
    cat("Done!","\n")
    cat("--------------------------------------------------------------------------","\n")
    cat("\nSummary of the adjusted models","\n")
    cat("--------------------------------------------------------------------------","\n")
    print(statistics, row.names = F)
    cat("--------------------------------------------------------------------------")
    cofgrap = ggplot2::ggplot(statistics, ggplot2::aes(x = remaining, y = cophenetic))+
      ggplot2::geom_point(size = 3)+
      ggplot2::theme_bw()+
      ggplot2::geom_line(size = 1)+
      ggplot2::theme(axis.ticks.length = unit(.2, "cm"),
                     axis.text = ggplot2::element_text(size = 12, colour = "black"),
                     axis.title = ggplot2::element_text(size = 12, colour = "black"),
                     axis.ticks = ggplot2::element_line(colour = "black"),
                     plot.margin = margin(0.5, 0.5, 0.2, 0.6, "cm"),
                     axis.title.y = ggplot2::element_text(margin = margin(r=16)),
                     legend.title = ggplot2::element_blank(),
                     legend.text = ggplot2::element_text(size=12),
                     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
                     panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank(),
                     panel.grid.minor.y = ggplot2::element_blank())+
      ggplot2::labs(x = "Number of variables", y = "Cophenetic correlation")
    model = statistics$Model[which.max(statistics$cophenetic)]
    predvar = ModelEstimates[[model]]$namevars
    data = data.frame(original[(match(c(predvar), names(original)))])
    if (verbose == TRUE){
      cat("\nSuggested variables to be used in the analysis","\n")
      cat("--------------------------------------------------------------------------","\n")
      cat("The clustering was calculated with the ", model,
          "\nThe variables included in this model were...\n",
          predvar,"\n")
      cat("--------------------------------------------------------------------------")
    }
  } else {data = data
  cofgrap = NULL
  statistics = NULL
  }
  if (distmethod %in% c("pearson", "spearman", "kendall")){
    de = as.dist(1- cor(t(data), method = distmethod))
  } else {
    de = dist(data, method = distmethod, diag = T, upper = T)
  }
  mat = as.matrix(de)
  hc = hclust(de, method = clustmethod)
  d2 = cophenetic(hc)
  cof = cor(d2, de)
  k = 1.25
  pcorte = mean(hc$height) + k * sd(hc$height)
  if (!missing(nclust)){
    groups <- cutree(hc, k = nclust)
    Mgroups <- cbind(data, groups)
    distance <- hc$height[(length(hc$height) - nclust):length(hc$height)]
    Sim <- (1 - distance/max(de))
    Passos <- 1:length(Sim)
    Simgroups <- length(Sim):1
    similarity <- Sim * 100
    Tab <- cbind(Passos, Simgroups, round(similarity, 3), round(distance, 2))
    colnames(Tab) <- c("Steps", "Groups", "Similarity", "Distance")
    TabResgroups <- NULL
    MGr <- cbind(data, groups)
    for (i in 1:nclust) {
      NewGroups <- subset(MGr, groups == i)
      GrupCalc <- NewGroups[, 1:(ncol(NewGroups) - 1)]
      Qtd.Elementos <- nrow(NewGroups)
      if (Qtd.Elementos == 1)
        Media <- GrupCalc
      else Media <- apply(GrupCalc, 2, mean)
      if (Qtd.Elementos == 1)
        SqG <- 0
      else SqG <- sum(sweep(GrupCalc, 2, Media)^2)
      TabResgroups <- rbind(TabResgroups, c(i, Qtd.Elementos,
                                            SqG, Media))
    }
    colnames(TabResgroups) <- c("Cluster", "Number of Elements",
                                "Sum_sq", paste(colnames(TabResgroups[, 4:(ncol(TabResgroups))])))
  } else{
    TabResgroups = NULL
    Tab = NULL
  }
  Sqt <- sum(sweep(data, 2, apply(data, 2, mean))^2)
  return(structure(list(data = data,
                        cutpoint = pcorte,
                        distance = mat,
                        de = de,
                        hc = hc,
                        cophenetic = cof,
                        Sqt = Sqt,
                        tab = as.data.frame(Tab),
                        clusters = as.data.frame(TabResgroups),
                        cofgrap = cofgrap,
                        statistics = statistics),
                   class = "clustering"))
}
}
