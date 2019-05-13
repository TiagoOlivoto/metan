#' Mahalanobis distance using data from a designed experiment
#'
#' Compute the Mahalanobis distance using data from an experiment conducted in
#' a randomized complete block design or completely randomized design.
#'
#'
#' @param .data The dataset containing the columns related to Genotypes,
#' replication/block and response variables. Alternatively, it is possible to
#' use an object of class 'split_factors' to compute the distance for each
#' level of the grouping factor. See \code{?split_factors}.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#' replications/blocks.
#' @param resp The response variables. For example \code{resp = c(var1, var2,
#' var3)}.
#' @param design The experimental design. Must be RCBD or CRD.
#' @param return What the function return? Default is "distance", i.e., the
#' Mahalanobis distance. Alternatively, it is possible to return the matrix of
#' means \code{return = "means"}, or the variance-covariance matrix of
#' residuals \code{return = "covmat"}.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'
#' library(metan)
#' library(dplyr)
#' maha_group = mahala_design(data_ge,
#'                            gen = GEN,
#'                            rep = REP,
#'                            resp = c(GY, HM))
#'
#' # Compute one distance for each environment
#' maha_group = data_ge %>%
#'              split_factors(ENV, keep_factors = TRUE) %>%
#'              mahala_design(GEN, REP, c(GY, HM))
#'
#' # Return the variance-covariance matrix of residuals
#' cov_mat = mahala_design(data_ge,
#'                            gen = GEN,
#'                            rep = REP,
#'                            resp = c(GY, HM),
#'                            return = "covmat")
#'
mahala_design = function(.data, gen, rep, resp, design = "RCBD", return = "distance"){
  if (!design %in% c("RCBD", "CRD")) {
    stop("The experimental design must be RCBD or CRD.")
  }
  if(any(class(.data) == "split_factors")){
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
