plot.clustering = function(x, type = "dendrogram", ...){
  if (type == "dendrogram"){
    plot(as.dendrogram(x$hc), ...)
  }
  if (type == "cophenetic"){
    return(x$cofgrap)
  }
}
