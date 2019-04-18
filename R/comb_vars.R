comb_vars = function(.data, order = "second", FUN = "+", verbose = TRUE){
  FUN = match.fun(FUN)
  if (!order %in% c("first", "second")) {
    stop("The orde must be one of 'first' or 'second'.")
  }
  if (verbose == TRUE){
    if (sum(lapply(.data, is.factor) == TRUE)>0){
    message("The columns factors ", paste0(collapse = " ", names(.data[ , unlist(lapply(.data, is.factor))])),
            " where deleted. Only the numeric variables were used.")
    }
  }
  x = .data[ , unlist(lapply(.data, is.numeric))] %>% as.data.frame()
  cb = data.frame(t(combn(ncol(x), 2)))
  if (order == "second"){
  cb = cb[with(cb, order(X2)), ]
  }
  nvars = data.frame(matrix(nrow = nrow(x), ncol = nrow(cb)))
  for (i in 1:nrow(cb)){
    nvars[, i] = FUN(x[[cb[i, 1]]], x[[cb[i, 2]]])
    colnames(nvars)[[i]] =  paste(colnames(x[cb[i, 1]]), colnames(x[cb[i, 2]]))
  }
  return(data.frame(nvars))
}
