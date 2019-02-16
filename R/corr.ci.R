corr.ci = function(data = NA, r = NULL, n = NULL) {
  if(any(!is.na(data))){

    if(is.matrix(data) == TRUE & is.null(n) == TRUE){
      stop("You have a matrix as entire but do not declared the sample size used to estimate the matrix")
    }

    if(is.matrix(data) != TRUE & is.null(n) != TRUE){
      stop("You have a dataframe as entire but declared the sample size used. When a dataframe is used, the sample size is automatically caltulated.")
    }
    if(is.null(n)==TRUE){
      n = nrow(data)
    }else{n = n}

    if(any(!is.na(data))) {
      data = cor(data, use = "complete", method = "pearson")
    }
    m = as.matrix(data)

  results = data.frame(Corr = as.vector(t(m)[lower.tri(m, diag=F)]))
  CI =  dplyr::mutate(results,
                           CI = (0.45304^abs(Corr))*2.25152*(n^-0.50089),
                           LL = Corr - CI,
                           UL = Corr + CI)
  combnam = combn(colnames(data), 2, paste, collapse = " x ")
  rownames(CI) = names(sapply(combnam, names))

  return(CI)
  } else{
  CI=(0.45304^r)*2.25152*(n^-0.50089)
  UP=r+CI
  LP=r-CI
  cat("-------------------------------------------------","\n")
  cat("Nonparametric 95% half-width confidence interval", "\n")
  cat("-------------------------------------------------","\n")
  cat(paste0("Level of significance: 5%","\n",
                                    "Correlation coefficient: ", r,
                                    "\nSample size: ", n,
                                    "\nConfidence interval: ", round(CI,4),
                                    "\nTrue parameter range from: ", round(LP,4)," to ",
                                    round(UP,4)),"\n")
  cat("-------------------------------------------------","\n")
  }
}
