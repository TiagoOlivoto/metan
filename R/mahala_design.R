mahala_design = function(.data, gen, rep, resp, design = "RCBD", return = "distance"){
  if (!design %in% c("RCBD", "CRD")) {
    stop("The experimental design must be RCBD or CRD.")
  }
  if(any(class(.data) == "group_factors")){
    dfs = list()
    for (k in 1:length(.data)){
      datain <- .data[[k]]
      nam  = names(.data[k])
      GEN <- factor(eval(substitute(gen), eval(datain)))
      REP <- factor(eval(substitute(rep), eval(datain)))
      d <- match.call()
      nvar <- as.numeric(ifelse(length(d$resp) > 1, length(d$resp) - 1, length(d$resp)))
      mat = matrix(nrow = nvar, ncol = nvar)
      covdata = data.frame(matrix(nrow = nrow(datain), ncol = nvar))
      vin <- 0
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
          diag(mat)[vin] = model[3, 3]
        } else {
          model = anova(aov(Y ~ GEN))
          diag(mat)[vin] = model[2, 3]
        }
        colnames(covdata)[[vin]] = paste(d$resp[var])
      }
      means = data.frame(cbind(GEN, covdata)) %>%
        dplyr::group_by(GEN) %>%
        dplyr::summarise_all(mean) %>%
        dplyr::select(-GEN)
      covdata2 = comb_vars(data.frame(covdata), order = "second")
      index = data.frame(t(combn(ncol(mat), 2)))
      index = index[with(index, order(X2)), ]
      temp = NULL
      for (i in 1:ncol(covdata2)){
        if (design == "RCBD"){
          model = anova(aov(covdata2[[i]] ~ GEN + REP))
          temp[i] = (model[3, 3] - diag(mat)[[index[i, 1]]] - diag(mat)[[index[i, 2]]])/2
        } else {
          model = anova(aov(covdata2[[i]] ~ GEN))
          temp[i] = (model[2, 3] - diag(mat)[[index[i, 1]]] - diag(mat)[[index[i, 2]]])/2
        }
      }
      mat[upper.tri(mat, diag = F)] = temp
      mat[lower.tri(mat, diag = F)] = temp
      rownames(mat) = colnames(means)
      colnames(mat) = colnames(means)
      dist = mahala(.means = means, covar = mat, inverted = FALSE)
      if (return == "distance"){
        dfs[[paste(nam)]] = dist
      }
      if (return == "covmat"){
        dfs[[paste(nam)]] = mat
      }
      if (return == "means"){
        dfs[[paste(nam)]] = means
      }
    }
    return(structure(dfs, class = "mahala_group"))
  } else {
    datain <- .data
    GEN <- factor(eval(substitute(gen), eval(datain)))
    REP <- factor(eval(substitute(rep), eval(datain)))
    d <- match.call()
    nvar <- as.numeric(ifelse(length(d$resp) > 1, length(d$resp) - 1, length(d$resp)))
    mat = matrix(nrow = nvar, ncol = nvar)
    covdata = data.frame(matrix(nrow = nrow(.data), ncol = nvar))
    vin <- 0
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
        diag(mat)[vin] = model[3, 3]
      } else {
        model = anova(aov(Y ~ GEN))
        diag(mat)[vin] = model[2, 3]
      }
      colnames(covdata)[[vin]] = paste(d$resp[var])
    }
    means = data.frame(cbind(GEN, covdata)) %>%
            dplyr::group_by(GEN) %>%
            dplyr::summarise_all(mean) %>%
            dplyr::select(-GEN)
    covdata2 = comb_vars(data.frame(covdata), order = "second")
    index = data.frame(t(combn(ncol(mat), 2)))
    index = index[with(index, order(X2)), ]
    temp = NULL
    for (i in 1:ncol(covdata2)){
      if (design == "RCBD"){
        model = anova(aov(covdata2[[i]] ~ GEN + REP))
        temp[i] = (model[3, 3] - diag(mat)[[index[i, 1]]] - diag(mat)[[index[i, 2]]])/2
      } else {
        model = anova(aov(covdata2[[i]] ~ GEN))
        temp[i] = (model[2, 3] - diag(mat)[[index[i, 1]]] - diag(mat)[[index[i, 2]]])/2
      }
    }
    mat[upper.tri(mat, diag = F)] = temp
    mat[lower.tri(mat, diag = F)] = temp
    rownames(mat) = colnames(means)
    colnames(mat) = colnames(means)
    dist = mahala(.means = means, covar = mat, inverted = FALSE)
    if (return == "distance"){
      return(dist)
    }
    if (return == "covmat"){
      return(mat)
    }
    if (return == "means"){
      return(means)
    }
  }
}
