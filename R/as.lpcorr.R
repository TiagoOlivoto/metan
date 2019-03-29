as.lpcor = function(x, y, ...){
  data = list(x, y, ...)
  if(length(data)== 1){
    stop(call. = FALSE, "At least two matrices must be used")
  }
  if(length(which(sapply(data, function(x) identical(nrow(x), ncol(x)))==TRUE)) != length(data)){
    stop(call. = FALSE, "All matrices in the list must be a square matrix. Please, check and fix.")
  }
  if (length(unique(unique(sapply(data, function(x) dim(x)))[1,1:length(data)])) != 1){
    stop(call. = FALSE, "All matrices in the list must be the same dimension. Please, check and fix.")
  }
  invisible(structure(data, class = "lpcor"))
}


