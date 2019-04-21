as.lpcor = function(...){
  data = list(...)
  if(length(which(sapply(data, function(x) identical(nrow(x), ncol(x)))==TRUE)) != length(data)){
    stop(call. = FALSE, "All matrices in the list must be a square matrix. Please, check and fix.")
  }
  if (length(unique(unique(sapply(data, function(x) dim(x)))[1,1:length(data)])) != 1){
    stop(call. = FALSE, "All matrices in the list must be the same dimension. Please, check and fix.")
  }
  data =  lapply(data, function(x) as.matrix(x))
  names(data) = paste("mat", 1:length(data), sep = "_")
  invisible(structure(data, class = "lpcor"))
}



