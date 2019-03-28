colindiag = function(x, group_var = NULL, n = NULL, verbose = TRUE){
  if(!is.matrix(x) && !is.data.frame(x)){
    stop("The object 'x' must be a correlation matrix or a data.frame")
  }
  if(is.matrix(x) && is.null(n)){
    stop("You have a matrix but the sample size used to compute the correlations (n) was not declared.")
  }
  
  internal = function(x){
    if(is.matrix(x)){
      cor.x = x
    }
    if(is.data.frame(x)){
      cor.x = cor(x)
      n = nrow(x)
    }
    eigen = eigen(cor.x)
    Det = det(cor.x)
    NC = max(eigen$values)/min(eigen$values)
    Aval = data.frame(eigen$values)
    names(Aval) = "Eigenvalues"
    Avet = data.frame(t(eigen$vectors))
    names(Avet) = colnames(x) 
    AvAvet = cbind(Aval, Avet)
    VIF = data.frame(diag(solve(cor.x)))
    names(VIF) = "VIF"
    results = data.frame(linear = as.vector(t(cor.x)[lower.tri(cor.x,diag=F)]))
    results = dplyr::mutate(results,
                            t = linear*(sqrt(n-2))/(sqrt(1-linear^2)),
                            prob = 2*(1-pt(abs(t), df = n-2)))
    names = colnames(x)
    combnam = combn(names,2, paste, collapse = " x ")
    rownames(results) = names(sapply(combnam, names))
    largest_corr = paste0(rownames(results)[which.max(abs(results$linear))], " = ",
                          round(results[which.max(abs(results$linear)),1],3))
    smallest_corr = paste0(rownames(results)[which.min(abs(results$linear))], " = ",
                           round(results[which.min(abs(results$linear)),1],3))
    ncorhigh = sum(results$linear >= abs(0.8))
    if(verbose == TRUE){
      if (NC > 1000){
      cat(paste0("Severe multicollinearity in the matrix! Pay attention on the variables listed bellow\n",
                 "CN = ", round(NC,3), "\n"))
    }
    if (NC < 100){
      cat(paste0("Weak multicollinearity in the matrix\n",
                 "NC = ", round(NC,3), "\n"))
    }
    if (NC > 100 & NC < 1000 ){
      cat(paste0("The multicollinearity in the matrix should be investigated.\n",
                 "NC = ", round(NC,3), "\n",
                 "Largest VIF = ", max(VIF), "\n"))
    }
  }
    ultimo = data.frame(Peso=t(AvAvet[c(nrow(AvAvet)),])[-c(1),])
    abs = data.frame(Peso = abs(ultimo[,"Peso"]))
    rownames(abs) = rownames(ultimo)
    ultimo = abs[order(abs[,"Peso"], decreasing = T), , drop = FALSE]
    pesovarname = paste(rownames(ultimo), collapse = ' > ')
    final = list(cormat = data.frame(cor.x),
                 corlist = results,
                 evalevet = AvAvet,
                 VIF = VIF,
                 CN = NC,
                 det = Det,
                 largest_corr = largest_corr,
                 smallest_corr = smallest_corr, 
                 weight_var = pesovarname)
    if(verbose == TRUE){
    cat(paste0("Matrix determinant: ", round(Det,7)),  "\n")
    cat(paste0("Largest correlation: ", largest_corr),  "\n")
    cat(paste0("Smallest correlation: ", smallest_corr),  "\n")
    cat(paste0("Number of correlations with r >= |0.8|: ", ncorhigh),  "\n")
    cat(paste0("Variables with largest weight in the last eigenvalues: ","\n",
               pesovarname),  "\n")
    }
    
    return(invisible(final))
  }
  
  if(is.matrix(x)){
    out = internal(x)
  }
  
  if(is.data.frame(x)){
  if(!missing(group_var)){
    group_var <- dplyr::enquo(group_var)
  re =  x  %>%
    split(dplyr::pull(., !!group_var))
data = lapply(re, function(x){
     unlist(lapply(x, is.factor))
    message("The factors ", paste0(collapse = " ", names(x[ , unlist(lapply(x, is.factor)) ])),
            " where excluded from the data")
  x[ , unlist(lapply(x, is.numeric))]
  }
  )
out = lapply(seq_along(data), function(x){
  if(verbose == TRUE){
  cat("\n----------------------------------------------------------------------------\n")
  cat("Level", names(data)[[x]], "\n")
  cat("----------------------------------------------------------------------------\n")
  }
  internal(data[[x]])
})
names(out) = names(data)
  } else{
    if(sum(lapply(x, is.factor)==TRUE)>0){
    message("The factors ", paste0(collapse = " ", names(x[ , unlist(lapply(x, is.factor)) ])),
            " where excluded from the data")
    }
    data = x[ , unlist(lapply(x, is.numeric))]
    out = internal(data)
  }
}
  invisible(out)
}
