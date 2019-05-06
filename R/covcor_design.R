covcor_design = function(.data, gen, rep, resp, design = "RCBD", type = NULL){
  if (!design %in% c("RCBD", "CRD")) {
    stop("The experimental design must be RCBD or CRD.")
  }
  if (!is.null(type)){
    if (!type %in% c(c("pcor", "gcor", "rcor", "pcov", "gcov", "rcov", "means"))){
      stop("The type must be one of the 'pcor', 'gcor', 'rcor', 'pcov', 'gcov', 'rcov', or 'means'. ")
    }
  }
  if(any(class(.data) == "group_factors")){
    dfs = list()
    for (k in 1:length(.data)){
      datain <- .data[[k]]
      nam  = names(.data[k])
      GEN <- factor(eval(substitute(gen), eval(datain)))
      REP <- factor(eval(substitute(rep), eval(datain)))
      NREP <- length(unique(REP))
      d <- match.call()
      nvar <- as.numeric(ifelse(length(d$resp) > 1, length(d$resp) - 1, length(d$resp)))
      covdata = data.frame(matrix(nrow = nrow(datain), ncol = nvar))
      vin <- 0
      mst = NULL
      msr = NULL
      for (var in 2:length(d$resp)) {
        vin <- vin + 1
        if (length(d$resp) > 1) {
          Y <- eval(substitute(resp)[[var]], eval(datain))
        } else {
          Y <- eval(substitute(resp), eval(datain))
        }
        covdata[, vin] = Y
        if (design == "RCBD"){
          model = anova(aov(Y ~ GEN + REP))
          mst[vin] = model[1, 3]
          msr[vin] = model[3, 3]
        } else {
          model = anova(aov(Y ~ GEN))
          mst[vin] = model[1, 3]
          msr[vin] = model[2, 3]
        }
        colnames(covdata)[[vin]] = paste(d$resp[var])
      }
      ms = data.frame(mst = mst, msr = msr) %>% dplyr::mutate(tr = mst - msr)
      vres = diag(ms[,2])
      vfen = diag(ms[,1] / 3)
      vgen = (diag(ms[, 1]) - diag(ms[, 2])) / 3
      means = data.frame(cbind(GEN, covdata)) %>%
        dplyr::group_by(GEN) %>%
        dplyr::summarise_all(mean) %>%
        dplyr::ungroup() %>%
        dplyr::select(-GEN) %>%
        as.data.frame()
      covdata2 = comb_vars(data.frame(covdata), order = "second")
      index = data.frame(t(combn(5, 2)))
      index = index[with(index, order(X2)), ]
      covres = NULL
      covfen = NULL
      covgen = NULL
      cores = NULL
      corfen = NULL
      corgen = NULL
      for (i in 1:ncol(covdata2)){
        if (design == "RCBD"){
          model = anova(aov(covdata2[[i]] ~ GEN + REP))
          tcovres = (model[3, 3] - ms[index[i, 1], 2] - ms[index[i, 2], 2])/2
          tcovfen = ((model[1, 3] - ms[index[i, 1], 1] - ms[index[i, 2], 1])/2) / NREP
          tcovgen = (tcovfen * NREP - tcovres) / NREP
          covres[i] = tcovres
          covfen[i] = tcovfen
          covgen[i] = tcovgen
          corfen[i] = tcovfen / sqrt((ms[index[i, 1], 1] / NREP) * (ms[index[i, 2], 1] / NREP))
          corgen[i] = tcovgen / sqrt((ms[index[i, 1], 3] / NREP) * (ms[index[i, 2], 3] / NREP))
          cores[i] = tcovres / sqrt((ms[index[i, 1], 2]) * (ms[index[i, 2], 2]))
        } else {
          model = anova(aov(covdata2[[i]] ~ GEN))
          tcovres = (model[2, 3] - ms[index[i, 1], 2] - ms[index[i, 2], 2])/2
          tcovfen = ((model[1, 3] - ms[index[i, 1], 1] - ms[index[i, 2], 1])/2) / NREP
          tcovgen = (tcovfen * NREP - tcovres) / NREP
          covres[i] = tcovres
          covfen[i] = tcovfen
          covgen[i] = tcovgen
          corfen[i] = tcovfen / sqrt((ms[index[i, 1], 1] / NREP) * (ms[index[i, 2], 1] / NREP))
          corgen[i] = tcovgen / sqrt((ms[index[i, 1], 3] / NREP) * (ms[index[i, 2], 3] / NREP))
          cores[i] = tcovres / sqrt((ms[index[i, 1], 2]) * (ms[index[i, 2], 2]))
        }
      }
      corres = matrix(1, nvar, nvar)
      corrgen = matrix(1, nvar, nvar)
      corrfen = matrix(1, nvar, nvar)
      vres[lower.tri(vres, diag = F)] <- covres
      vfen[lower.tri(vfen, diag = F)] <- covfen
      vgen[lower.tri(vgen, diag = F)] <- covgen
      corres[lower.tri(corres, diag = F)] <- cores
      corrfen[lower.tri(corrfen, diag = F)] <- corfen
      corrgen[lower.tri(corrgen, diag = F)] <- corgen
      colnames(vres) <- rownames(vres) <- names(means)
      colnames(vfen) <- rownames(vfen) <- names(means)
      colnames(vgen) <- rownames(vgen) <- names(means)
      colnames(corres) <- rownames(corres) <- names(means)
      colnames(corrfen) <- rownames(corrfen) <- names(means)
      colnames(corrgen) <- rownames(corrgen) <- names(means)

      if (is.null(type)){
        tmp = list(geno_cov = as.matrix(make_sym(vgen)),
                   phen_cov = as.matrix(make_sym(vfen)),
                   resi_cov = as.matrix(make_sym(vres)),
                   geno_cor = as.matrix(make_sym(corrgen)),
                   phen_cor = as.matrix(make_sym(corrfen)),
                   resi_cor = as.matrix(make_sym(corres)))
        dfs[[paste(nam)]] = tmp
      } else {
        if(type == "pcor"){
          dfs[[paste(nam)]] =  as.matrix(corrfen)
        }
        if(type == "gcor"){
          dfs[[paste(nam)]] = as.matrix(corrgen)
        }
        if(type == "rcor"){
          dfs[[paste(nam)]] = as.matrix(corres)
        }
        if(type == "pcov"){
          dfs[[paste(nam)]] = as.matrix(vfen)
        }
        if(type == "gcov"){
          dfs[[paste(nam)]] = as.matrix(vgen)
        }
        if(type == "rcov"){
          dfs[[paste(nam)]] = as.matrix(vres)
        }
        if(type == "means"){
          return(means)
        }
      }
    }
    return(structure(dfs, class = "covcor_design"))
  } else {
  datain <- .data
  GEN <- factor(eval(substitute(gen), eval(datain)))
  REP <- factor(eval(substitute(rep), eval(datain)))
  NREP <- length(unique(REP))
  d <- match.call()
  nvar <- as.numeric(ifelse(length(d$resp) > 1, length(d$resp) - 1, length(d$resp)))
  covdata = data.frame(matrix(nrow = nrow(.data), ncol = nvar))
  vin <- 0
  mst = NULL
  msr = NULL
  for (var in 2:length(d$resp)) {
    vin <- vin + 1
    if (length(d$resp) > 1) {
      Y <- eval(substitute(resp)[[var]], eval(datain))
    } else {
      Y <- eval(substitute(resp), eval(datain))
    }
    covdata[, vin] = Y
    if (design == "RCBD"){
      model = anova(aov(Y ~ GEN + REP))
      mst[vin] = model[1, 3]
      msr[vin] = model[3, 3]
    } else {
      model = anova(aov(Y ~ GEN))
      mst[vin] = model[1, 3]
      msr[vin] = model[2, 3]
    }
    colnames(covdata)[[vin]] = paste(d$resp[var])
  }
    ms = data.frame(mst = mst, msr = msr) %>% dplyr::mutate(tr = mst - msr)
    vres = diag(ms[,2])
    vfen = diag(ms[,1] / 3)
    vgen = (diag(ms[, 1]) - diag(ms[, 2])) / 3
    means = data.frame(cbind(GEN, covdata)) %>%
      dplyr::group_by(GEN) %>%
      dplyr::summarise_all(mean) %>%
      dplyr::ungroup() %>%
      dplyr::select(-GEN) %>%
      as.data.frame()
  covdata2 = comb_vars(data.frame(covdata), order = "second")
  index = data.frame(t(combn(5, 2)))
  index = index[with(index, order(X2)), ]
  covres = NULL
  covfen = NULL
  covgen = NULL
  cores = NULL
  corfen = NULL
  corgen = NULL
  for (i in 1:ncol(covdata2)){
    if (design == "RCBD"){
      model = anova(aov(covdata2[[i]] ~ GEN + REP))
      tcovres = (model[3, 3] - ms[index[i, 1], 2] - ms[index[i, 2], 2])/2
      tcovfen = ((model[1, 3] - ms[index[i, 1], 1] - ms[index[i, 2], 1])/2) / NREP
      tcovgen = (tcovfen * NREP - tcovres) / NREP
      covres[i] = tcovres
      covfen[i] = tcovfen
      covgen[i] = tcovgen
      corfen[i] = tcovfen / sqrt((ms[index[i, 1], 1] / NREP) * (ms[index[i, 2], 1] / NREP))
      corgen[i] = tcovgen / sqrt((ms[index[i, 1], 3] / NREP) * (ms[index[i, 2], 3] / NREP))
      cores[i] = tcovres / sqrt((ms[index[i, 1], 2]) * (ms[index[i, 2], 2]))
    } else {
      model = anova(aov(covdata2[[i]] ~ GEN))
      tcovres = (model[2, 3] - ms[index[i, 1], 2] - ms[index[i, 2], 2])/2
      tcovfen = ((model[1, 3] - ms[index[i, 1], 1] - ms[index[i, 2], 1])/2) / NREP
      tcovgen = (tcovfen * NREP - tcovres) / NREP
      covres[i] = tcovres
      covfen[i] = tcovfen
      covgen[i] = tcovgen
      corfen[i] = tcovfen / sqrt((ms[index[i, 1], 1] / NREP) * (ms[index[i, 2], 1] / NREP))
      corgen[i] = tcovgen / sqrt((ms[index[i, 1], 3] / NREP) * (ms[index[i, 2], 3] / NREP))
      cores[i] = tcovres / sqrt((ms[index[i, 1], 2]) * (ms[index[i, 2], 2]))
    }
  }
  corres = matrix(1, nvar, nvar)
  corrgen = matrix(1, nvar, nvar)
  corrfen = matrix(1, nvar, nvar)
  vres[lower.tri(vres, diag = F)] <- covres
  vfen[lower.tri(vfen, diag = F)] <- covfen
  vgen[lower.tri(vgen, diag = F)] <- covgen
  corres[lower.tri(corres, diag = F)] <- cores
  corrfen[lower.tri(corrfen, diag = F)] <- corfen
  corrgen[lower.tri(corrgen, diag = F)] <- corgen
  colnames(vres) <- rownames(vres) <- names(means)
  colnames(vfen) <- rownames(vfen) <- names(means)
  colnames(vgen) <- rownames(vgen) <- names(means)
  colnames(corres) <- rownames(corres) <- names(means)
  colnames(corrfen) <- rownames(corrfen) <- names(means)
  colnames(corrgen) <- rownames(corrgen) <- names(means)

  if (is.null(type)){
  return(structure(list(geno_cov = as.matrix(make_sym(vgen)),
              phen_cov = as.matrix(make_sym(vfen)),
              resi_cov = as.matrix(make_sym(vres)),
              geno_cor = as.matrix(make_sym(corrgen)),
              phen_cor = as.matrix(make_sym(corrfen)),
              resi_cor = as.matrix(make_sym(corres)),
              means = means),
              class = "covcor_design"))
  }
  if(type == "pcor"){
    return(as.matrix(corrfen))
  }
  if(type == "gcor"){
    return(as.matrix(corrgen))
  }
  if(type == "rcor"){
    return(as.matrix(corres))
  }
  if(type == "pcov"){
    return(as.matrix(vfen))
  }
  if(type == "gcov"){
    return(as.matrix(vgen))
  }
  if(type == "rcov"){
    return(as.matrix(vres))
  }
  if(type == "means"){
    return(means)
  }
  }
}
