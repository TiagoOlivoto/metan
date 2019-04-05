path_coeff = function(.data,
                      resp,
                      pred = NULL,
                      exclude = FALSE,
                      correction = NULL,
                      knumber = 50,
                      brutstep = FALSE,
                      maxvif = 10,
                      missingval = "pairwise.complete.obs",
                      verbose = TRUE){
if (missing(resp) == TRUE){
  stop("A response variable (dependent) must be declared.")
}
kincrement = 1 / knumber
  if(any(class(.data) == "group_factors")){
    dfs = list()
    for (k in 1:length(.data)){
      data = .data[[k]]
      nam  = names(.data[k])
      resp = dplyr::enquo(resp)
      if(!missing(pred)){
        prcod <- dplyr::quos(!!dplyr::enquo(pred))
      }

      if (brutstep == FALSE){
        if (missing(pred)){
          pr = data %>% dplyr::select(-!!resp)

        } else {
          if (exclude == TRUE){
            pr = data %>% dplyr::select(-c(!!!prcod), -!!resp)
          } else{
            pr = data %>% dplyr::select(!!!prcod)
          }
        }
        pv = pr
        names = colnames(pv)
        y = data %>% select(!!resp)
        cor.y = cor(pv, y, use = missingval)
        cor.x = cor(pv, use = missingval)
        if (is.null(correction) == FALSE){
          diag(cor.x) = diag(cor.x) + correction

        } else {
          cor.x = cor(pv, use = missingval)
        }


        if (is.null(correction) == TRUE){
          betas = data.frame(matrix(nrow = knumber, ncol = length(pv)+1))
          cc = 0
          nvar = length(pv) + 1
          for (i in 1:knumber){
            cor.x2 = cor.x
            diag(cor.x2) = diag(cor.x2) + cc
            betas[i,1] = cc
            betas[i,2:nvar] = t(solve(cor.x2, cor.y))
            cc = cc + kincrement
          }
          names(betas) = paste0(c("K", names(pv)))
          xx <- betas$K
          yy <- colnames(betas)
          fila <- knumber
          col <- length(yy)
          total <- fila * (col-1)
          x <- character(length = total)
          y <- character(length = total)
          z <- numeric(length = total)
          k <- 0
          for (i in 1:fila) {
            for (j in 2:col) {
              k <- k  +  1
              x[k] <- xx[i]
              y[k] <- yy[j]
              z[k] <- betas[i, j]
            }
          }
          x = as.numeric(x)
          betas = data.frame(K = x, VAR = y, direct = z)


          p1 = ggplot2::ggplot(betas, ggplot2::aes(K, direct, col = VAR)) +
            ggplot2::geom_line(size = 1)+
            ggplot2::theme_bw() +
            ggplot2::theme(axis.ticks.length = unit(.2, "cm"),
                  axis.text = element_text(size = 12, colour = "black"),
                  axis.title = element_text(size = 12, colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  legend.position = "bottom",
                  plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                  legend.title = element_blank(),
                  panel.border = element_rect(colour = "black", fill = NA, size = 1),
                  panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
                  panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
            ggplot2::labs(x = "k values", y = "Beta values")+
            ggplot2::geom_abline(intercept = 0, slope = 0)+
            ggplot2::scale_x_continuous(limits = c(0, 1), breaks=seq(0,1, by = 0.1))
        } else {
          p1 = "No graphic generated due to correction value"
        }
        eigen = eigen(cor.x)
        Det = det(cor.x)
        NC = max(eigen$values)/min(eigen$values)
        Aval = data.frame(eigen$values)
        names(Aval) = "Eigenvalues"
        Avet = data.frame(t(eigen$vectors))
        names(Avet) = names
        AvAvet = cbind(Aval, Avet)
        Direct = solve(cor.x, cor.y)
        n = ncol(cor.x)
        Coeff = data.frame(cor.x)
        for (i in 1:n) {
          for (j in 1:n) {
            Coeff[i, j] = Direct[j] * cor.x[i, j]
          }
        }
        Residual = 1 - t(Direct) %*% cor.y
        R2 = t(Direct) %*% cor.y
        VIF = data.frame(diag(solve(cor.x)))
        if(verbose == TRUE){
        names(VIF) = "VIF"
        cat("\n----------------------------------------------------------------------------\n")
        cat("Level", nam, "\n")
        cat("----------------------------------------------------------------------------\n")
        if (NC > 1000){
          cat(paste0("Severe multicollinearity. \n",
                     "Condition Number = ", round(NC,3), "\n",
                     "Please, consider using a correction factor, or use 'brutstep = TRUE'. \n"))
        }

        if (NC < 100){
          cat(paste0("Weak multicollinearity. \n",
                    "Condition Number = ", round(NC,3), "\n",
                    "You will probably have path coefficients close to being unbiased. \n"))
        }

        if (NC > 100 & NC < 1000 ){
          cat(paste0("Moderate multicollinearity! \n",
                     "Condition Number = ", round(NC,3), "\n",
                     "Please, cautiosely evaluate the VIF and matrix determinant.", "\n"))
        }
        }
        last = data.frame(weight = t(AvAvet[c(nrow(AvAvet)),])[-c(1),])
        abs = data.frame(weight = abs(last[,"weight"]))
        rownames(abs) = rownames(last)
        last = abs[order(abs[,"weight"], decreasing = T), , drop = FALSE]
        weightvarname = paste(rownames(last), collapse = " > ")

        temp = structure(list(Corr.x = data.frame(cor.x),
                              Corr.y = data.frame(cor.y),
                              Coefficients = data.frame(t(Coeff)),
                              Eigen = AvAvet,
                              VIF = VIF,
                              plot = p1,
                              Predictors = pred,
                              CN = NC,
                              Det = Det,
                              R2 = R2,
                              Residual = Residual,
                              Response = resp,
                              weightvar = weightvarname),
                         class = "path_coeff")
        dfs[[paste(nam)]] = temp
      }


      if (brutstep == TRUE){
        yyy = data %>% dplyr::select(!!resp)
        xxx = data %>% dplyr::select(-c(!!resp))
        cor.xx = cor(xxx, use = missingval)
        VIF = data.frame(diag(solve(cor.xx)))
        names(VIF) = "VIF"
        VIF =  VIF[order(VIF[,"VIF"], decreasing = F), , drop = FALSE]

        repeat{
          VIF2 = VIF[order(VIF[-1,], decreasing = F), , drop = FALSE]
          pred2 = rownames(VIF2)
          xxx2 = data[rownames(VIF2)]
          VIF3 = data.frame(VIF = diag(solve(cor(xxx2, use = missingval))))
          VIF3 = VIF3[order(VIF3[,"VIF"], decreasing = F), , drop = FALSE]
          if (max(VIF3$VIF) <= maxvif) break
          VIF = VIF3
        }
        xxx = data[rownames(VIF3)]

        selectedpred = rownames(VIF3)
        npred = ncol(xxx)-1
        statistics = data.frame(matrix(nrow = npred-1, ncol = 8))
        ModelEstimates = list()
        modelcode = 1
        nproced = npred - 1
        if (verbose == TRUE){
        cat("\n----------------------------------------------------------------------------\n")
        cat("Level", nam, "\n")
        cat("----------------------------------------------------------------------------\n")
        cat(paste("The algorithm has selected a set of ",nrow(VIF3),
                  "predictors with largest VIF = ", round(max(VIF3$VIF),3),".", sep = ""), "\n")
        cat("Selected predictors:",paste0(selectedpred), "\n")
        cat(paste("A forward stepwise-based selection procedure will fit",nproced, "models.", "\n"))
        cat("--------------------------------------------------------------------------","\n")
        }
        for (i in 1:nproced){
          FDSel =  FWDselect::selection(x = xxx,
                                        y = unlist(yyy),
                                        q = npred,
                                        method = "lm",
                                        criterion = "aic",
                                        cluster = F)
          predstw = FDSel$Variable_names
          x = data[,c(predstw)]
          names = colnames(x)
          y = data %>% dplyr::select(!!resp)
          cor.y = cor(x, y, use = missingval)
          cor.x = cor(x, use = missingval)
          if (is.null(correction) == F){
            diag(cor.x) = diag(cor.x) + correction
          } else {
            cor.x = cor(x, use = missingval)
          }
          if (is.null(correction) == T){
            betas = data.frame(matrix(nrow = knumber, ncol = length(predstw)+1))
            cc = 0
            nvar = length(predstw) + 1
            for (i in 1:knumber){
              cor.x2 = cor.x
              diag(cor.x2) = diag(cor.x2) + cc
              betas[i,1] = cc
              betas[i,2:nvar] = t(solve(cor.x2, cor.y))
              cc = cc + kincrement
            }
            names(betas) = paste0(c("K", predstw))
            xx <- betas$K
            yy <- colnames(betas)
            fila <- knumber
            col <- length(yy)
            total <- fila * (col-1)
            x <- character(length = total)
            y <- character(length = total)
            z <- numeric(length = total)
            k <- 0
            for (i in 1:fila) {
              for (j in 2:col) {
                k <- k  +  1
                x[k] <- xx[i]
                y[k] <- yy[j]
                z[k] <- betas[i, j]
              }
            }
            x = as.numeric(x)
            betas = data.frame(K = x, VAR = y, direct = z)
            p1 = ggplot2::ggplot(betas, ggplot2::aes(K, direct, col = VAR)) +
              ggplot2::geom_line(size = 1)+
              ggplot2::theme_bw() +
              ggplot2::theme(axis.ticks.length = unit(.2, "cm"),
                    axis.text = element_text(size = 12, colour = "black"),
                    axis.title = element_text(size = 12, colour = "black"),
                    axis.ticks = element_line(colour = "black"),
                    legend.position = "bottom",
                    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                    legend.title = element_blank(),
                    panel.border = element_rect(colour = "black", fill = NA, size = 1),
                    panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
                    panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
              ggplot2::labs(x = "k values", y = expression(paste( beta, " values")))+
              ggplot2::geom_abline(intercept = 0, slope = 0)+
              ggplot2::scale_x_continuous(breaks=seq(0,1,by=0.1))

          } else {
            p1 = "No graphic generated due to correction value"
          }
          eigen = eigen(cor.x)
          Det = det(cor.x)
          NC = max(eigen$values)/min(eigen$values)
          Aval = as.data.frame(eigen$values)
          names(Aval) = "Eigenvalues"
          Avet = as.data.frame(eigen$vectors)
          names(Avet) = names
          AvAvet = cbind(Aval, Avet)
          Direct = solve(cor.x, cor.y)
          n = ncol(cor.x)
          Coeff = data.frame(cor.x)
          for (i in 1:n) {
            for (j in 1:n) {
              Coeff[i, j] = Direct[j] * cor.x[i, j]
            }
          }
          Residual = 1 - t(Direct) %*% cor.y
          R2 = t(Direct) %*% cor.y
          VIF = data.frame(diag(solve(cor.x)))
          names(VIF) = "VIF"
          last = data.frame(weight=t(AvAvet[c(nrow(AvAvet)),])[-c(1),])
          abs = data.frame(weight = abs(last[,"weight"]))
          rownames(abs) = rownames(last)
          last = abs[order(abs[,"weight"], decreasing = T), , drop = FALSE]
          weightvarname = paste(rownames(last), collapse = ' > ')
          cor.y = data.frame(cor.y)
          Results = list(Corr.x = data.frame(cor.x),
                         Corr.y = data.frame(cor.y),
                         Coefficients = data.frame(t(Coeff)),
                         Eigen = AvAvet,
                         VIF = VIF,
                         plot = p1,
                         Predictors = pred,
                         CN = NC,
                         Det = Det,
                         R2 = R2,
                         Residual = Residual,
                         Response = resp,
                         weightvar = weightvarname)
          ModelEstimates[[paste("Model", modelcode)]] = Results
          statistics[i, 1] = paste("Model", modelcode)
          statistics[i, 2] = FDSel$Information_Criterion
          statistics[i, 3] = npred
          statistics[i, 4] = NC
          statistics[i, 5] = Det
          statistics[i, 6] = R2
          statistics[i, 7] = Residual
          statistics[i, 8] = max(VIF)
          cat(paste("Adjusting the model ",modelcode, " with ", npred,
                    " predictor variables (",round(modelcode/nproced*100,2),"% concluded)","\n"), sep = "")
          npred = npred - 1
          modelcode = modelcode + 1
        }
        statistics = statistics[-c(1),]
        names(statistics) = c("Model", "AIC", "Numpred", "CN", "Determinant", "R2", "Residual", "maxVIF")
        if(verbose == TRUE){
        cat("Done!","\n")
        cat("\n\n")
        cat("--------------------------------------------------------------------------","\n")
        cat("Summary of the adjusted models","\n")
        cat("--------------------------------------------------------------------------","\n")
        print(statistics)
        cat("--------------------------------------------------------------------------")
        }
        temp = list(Models = ModelEstimates,
                    Summary = statistics,
                    Selectedpred = selectedpred)
        dfs[[paste(nam)]] = temp
      }
    }
    return(dfs)
  } else{

    if(sum(lapply(.data, is.factor)==TRUE)>0){
      if (verbose == TRUE){
      message("The factors ", paste0(collapse = " ", names(.data[ , unlist(lapply(.data, is.factor)) ])),
              " where excluded to perform the analysis. If you want to perform an analysis for each level of a factor, use the function 'group_factors() before.' ")
      }
    }
    data = .data[ , unlist(lapply(.data, is.numeric))]
    nam  = names(.data)
    resp = dplyr::enquo(resp)
    if(!missing(pred)){
      prcod <- dplyr::quos(!!dplyr::enquo(pred))
    }

    if (brutstep == FALSE){
      if (missing(pred)){
        pr = data %>% dplyr::select(-!!resp)

      } else {
        if (exclude == TRUE){
          pr = data %>% dplyr::select(-c(!!!prcod), -!!resp)
        } else{
          pr = data %>% dplyr::select(!!!prcod)
        }
      }

      names = colnames(pr)
      y = data %>% select(!!resp)
      cor.y = cor(pr, y, use = missingval)
      cor.x = cor(pr, use = missingval)
      if (is.null(correction) == FALSE){
        diag(cor.x) = diag(cor.x) + correction

      } else {
        cor.x = cor(pr, use = missingval)
      }


      if (is.null(correction) == TRUE){
        betas = data.frame(matrix(nrow = knumber, ncol = length(pr)+1))
        cc = 0
        nvar = length(pr) + 1
        for (i in 1:knumber){
          cor.x2 = cor.x
          diag(cor.x2) = diag(cor.x2) + cc
          betas[i,1] = cc
          betas[i,2:nvar] = t(solve(cor.x2, cor.y))
          cc = cc + kincrement
        }
        names(betas) = paste0(c("K", names(pr)))
        xx <- betas$K
        yy <- colnames(betas)
        fila <- knumber
        col <- length(yy)
        total <- fila * (col-1)
        x <- character(length = total)
        y <- character(length = total)
        z <- numeric(length = total)
        k <- 0
        for (i in 1:fila) {
          for (j in 2:col) {
            k <- k  +  1
            x[k] <- xx[i]
            y[k] <- yy[j]
            z[k] <- betas[i, j]
          }
        }
        x = as.numeric(x)
        betas = data.frame(K = x, VAR = y, direct = z)


        p1 = ggplot2::ggplot(betas, ggplot2::aes(K, direct, col = VAR)) +
          ggplot2::geom_line(size = 1)+
          ggplot2::theme_bw() +
          ggplot2::theme(axis.ticks.length = unit(.2, "cm"),
                axis.text = element_text(size = 12, colour = "black"),
                axis.title = element_text(size = 12, colour = "black"),
                axis.ticks = element_line(colour = "black"),
                legend.position = "bottom",
                plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                legend.title = element_blank(),
                panel.border = element_rect(colour = "black", fill = NA, size = 1),
                panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
                panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
          ggplot2::labs(x = "k values", y = "Beta values")+
          ggplot2::geom_abline(intercept = 0, slope = 0)+
          ggplot2::scale_x_continuous(breaks=seq(0,1,by=0.1))
      } else {
        p1 = "No graphic generated due to correction value"
      }
      eigen = eigen(cor.x)
      Det = det(cor.x)
      NC = max(eigen$values)/min(eigen$values)
      Aval = data.frame(eigen$values)
      names(Aval) = "Eigenvalues"
      Avet = data.frame(t(eigen$vectors))
      names(Avet) = names
      AvAvet = cbind(Aval, Avet)
      Direct = solve(cor.x, cor.y)
      n = ncol(cor.x)
      Coeff = data.frame(cor.x)
      for (i in 1:n) {
        for (j in 1:n) {
          Coeff[i, j] = Direct[j] * cor.x[i, j]
        }
      }
      Residual = 1 - t(Direct) %*% cor.y
      R2 = t(Direct) %*% cor.y
      VIF = data.frame(diag(solve(cor.x)))
      names(VIF) = "VIF"
      if (verbose == TRUE){
        if (NC > 1000){
          cat(paste0("Severe multicollinearity. \n",
                     "Condition Number = ", round(NC,3), "\n",
                     "Please, consider using a correction factor, or use 'brutstep = TRUE'. \n"))
        }

        if (NC < 100){
          cat(paste0("Weak multicollinearity. \n",
                     "Condition Number = ", round(NC,3), "\n",
                     "You will probably have path coefficients close to being unbiased. \n"))
        }

        if (NC > 100 & NC < 1000 ){
          cat(paste0("Moderate multicollinearity! \n",
                     "Condition Number = ", round(NC,3), "\n",
                     "Please, cautiosely evaluate the VIF and matrix determinant.", "\n"))
        }
      }

      last = data.frame(weight = t(AvAvet[c(nrow(AvAvet)),])[-c(1),])
      abs = data.frame(weight = abs(last[,"weight"]))
      rownames(abs) = rownames(last)
      last = abs[order(abs[,"weight"], decreasing = T), , drop = FALSE]
      weightvarname = paste(rownames(last), collapse = " > ")

      temp = structure(list(Corr.x = data.frame(cor.x),
                            Corr.y = data.frame(cor.y),
                            Coefficients = data.frame(t(Coeff)),
                            Eigen = AvAvet,
                            VIF = VIF,
                            plot = p1,
                            Predictors = names(pr),
                            CN = NC,
                            Det = Det,
                            R2 = R2,
                            Residual = Residual,
                            Response = resp,
                            weightvar = weightvarname),
                       class = "path_coeff")
      return(temp)
    }


    if (brutstep == TRUE){
      yyy = data %>% dplyr::select(!!resp) %>% as.data.frame()
      xxx = data %>% dplyr::select(-c(!!resp)) %>% as.data.frame()
      cor.xx = cor(xxx, use = missingval)
      VIF = data.frame(diag(solve(cor.xx)))
      names(VIF) = "VIF"
      VIF =  VIF[order(VIF[,"VIF"], decreasing = F), , drop = FALSE]

      repeat{
        VIF2 = VIF[order(VIF[-1,], decreasing = F), , drop = FALSE]
        pred2 = rownames(VIF2)
        xxx2 = data[rownames(VIF2)]
        VIF3 = data.frame(VIF = diag(solve(cor(xxx2, use = missingval))))
        VIF3 = VIF3[order(VIF3[,"VIF"], decreasing = F), , drop = FALSE]
        if (max(VIF3$VIF) <= maxvif) break
        VIF = VIF3
      }
      xxx = data[rownames(VIF3)]

      selectedpred = rownames(VIF3)
      npred = ncol(xxx)-1
      statistics = data.frame(matrix(nrow = npred-1, ncol = 8))
      ModelEstimates = list()
      modelcode = 1
      nproced = npred - 1
      if (verbose == TRUE){
      cat("--------------------------------------------------------------------------","\n")
      cat(paste("The algorithm has selected a set of ",nrow(VIF3),
                "predictors with largest VIF = ", round(max(VIF3$VIF),3),".", sep = ""), "\n")
      cat("Selected predictors:",paste0(selectedpred), "\n")
      cat(paste("A forward stepwise-based selection procedure will fit",nproced, "models.", "\n"))
      cat("--------------------------------------------------------------------------","\n")
      }
      for (i in 1:nproced){
        FDSel =  FWDselect::selection(x = xxx,
                                      y = unlist(yyy),
                                      q = npred,
                                      method = "lm",
                                      criterion = "aic",
                                      cluster = F)
        predstw = FDSel$Variable_names
        x = data[,c(predstw)]
        names = colnames(x)
        y = data %>% dplyr::select(!!resp)

        cor.y = cor(x, y, use = missingval)
        cor.x = cor(x, use = missingval)
        if (is.null(correction) == F){
          diag(cor.x) = diag(cor.x) + correction
        } else {
          cor.x = cor(x, use = missingval)
        }
        if (is.null(correction) == T){
          betas = data.frame(matrix(nrow = knumber, ncol = length(predstw)+1))
          cc = 0
          nvar = length(predstw) + 1
          for (i in 1:knumber){
            cor.x2 = cor.x
            diag(cor.x2) = diag(cor.x2) + cc
            betas[i,1] = cc
            betas[i,2:nvar] = t(solve(cor.x2, cor.y))
            cc = cc + kincrement
          }
          names(betas) = paste0(c("K", predstw))
          xx <- betas$K
          yy <- colnames(betas)
          fila <- knumber
          col <- length(yy)
          total <- fila * (col-1)
          x <- character(length = total)
          y <- character(length = total)
          z <- numeric(length = total)
          k <- 0
          for (i in 1:fila) {
            for (j in 2:col) {
              k <- k  +  1
              x[k] <- xx[i]
              y[k] <- yy[j]
              z[k] <- betas[i, j]
            }
          }
          x = as.numeric(x)
          betas = data.frame(K = x, VAR = y, direct = z)

          p1 = ggplot2::ggplot(betas, ggplot2::aes(K, direct, col = VAR)) +
            ggplot2::geom_line(size = 1)+
            ggplot2::theme_bw() +
            ggplot2::theme(axis.ticks.length = unit(.2, "cm"),
                  axis.text = element_text(size = 12, colour = "black"),
                  axis.title = element_text(size = 12, colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  legend.position = "bottom",
                  plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                  legend.title = element_blank(),
                  panel.border = element_rect(colour = "black", fill = NA, size = 1),
                  panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
                  panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
            ggplot2::labs(x = "k values", y = expression(paste( beta, " values")))+
            ggplot2::geom_abline(intercept = 0, slope = 0)+
            ggplot2::scale_x_continuous(breaks=seq(0,1,by=0.1))

        } else {
          p1 = "No graphic generated due to correction value"
        }
        eigen = eigen(cor.x)
        Det = det(cor.x)
        NC = max(eigen$values)/min(eigen$values)
        Aval = as.data.frame(eigen$values)
        names(Aval) = "Eigenvalues"
        Avet = as.data.frame(eigen$vectors)
        names(Avet) = names
        AvAvet = cbind(Aval, Avet)
        Direct = solve(cor.x, cor.y)
        n = ncol(cor.x)
        Coeff = data.frame(cor.x)
        for (i in 1:n) {
          for (j in 1:n) {
            Coeff[i, j] = Direct[j] * cor.x[i, j]
          }
        }
        Residual = 1 - t(Direct) %*% cor.y
        R2 = t(Direct) %*% cor.y
        VIF = data.frame(diag(solve(cor.x)))
        names(VIF) = "VIF"
        last = data.frame(weight=t(AvAvet[c(nrow(AvAvet)),])[-c(1),])
        abs = data.frame(weight = abs(last[,"weight"]))
        rownames(abs) = rownames(last)
        last = abs[order(abs[,"weight"], decreasing = T), , drop = FALSE]
        weightvarname = paste(rownames(last), collapse = ' > ')
        cor.y = data.frame(cor.y)
        Results = list(Corr.x = data.frame(cor.x),
                       Corr.y = data.frame(cor.y),
                       Coefficients = data.frame(t(Coeff)),
                       Eigen = AvAvet,
                       VIF = VIF,
                       plot = p1,
                       Predictors = predstw,
                       CN = NC,
                       Det = Det,
                       R2 = R2,
                       Residual = Residual,
                       Response = resp,
                       weightvar = weightvarname)
        ModelEstimates[[paste("Model", modelcode)]] = Results
        statistics[i, 1] = paste("Model", modelcode)
        statistics[i, 2] = FDSel$Information_Criterion
        statistics[i, 3] = npred
        statistics[i, 4] = NC
        statistics[i, 5] = Det
        statistics[i, 6] = R2
        statistics[i, 7] = Residual
        statistics[i, 8] = max(VIF)
        if (verbose == TRUE){
        cat(paste("Adjusting the model ",modelcode, " with ", npred,
                  " predictors (",round(modelcode/nproced*100,2),"% concluded)","\n", sep = ""))
        }
        npred = npred - 1
        modelcode = modelcode + 1
      }
      statistics = statistics[-c(1),]
      names(statistics) = c("Model", "AIC", "Numpred", "CN", "Determinant", "R2", "Residual", "maxVIF")
      if (verbose == TRUE){
      cat("Done!","\n")
      cat("--------------------------------------------------------------------------","\n")
      cat("Summary of the adjusted models","\n")
      cat("--------------------------------------------------------------------------","\n")
      print(statistics)
      cat("--------------------------------------------------------------------------")
      }
      temp = list(Models = ModelEstimates,
                  Summary = statistics,
                  Selectedpred = selectedpred)
      return(temp)
    }
  }
}

