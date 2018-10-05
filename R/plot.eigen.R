plot.eigen = function(x,
                     export = FALSE,
                     theme = theme_waasb(),
                     file.type = "pdf",
                     file.name = NULL,
                     width = 6,
                     height = 6,
                     size.shape = 3.5,
                     size.line = 1,
                     col.shape = "blue",
                     col.line = "blue",
                     y.lab = "Eigenvalue",
                     x.lab = "Number of multiplicative terms",
                     resolution = 300,
                     ...){
  class = class(x)
  if(class != "WAASB"){
    stop("The object 'x' must be a 'WAASB' object.")
  }
eigen = x$PCA
p = ggplot2::ggplot(eigen, aes(x = PC, y = Eigenvalue, group = 1))  +
    geom_point(stat = 'identity', col = col.shape,  size = size.shape)   +
    geom_line(col = col.line, size = size.line) +
    theme +
    labs(x = x.lab, y = y.lab)

  if (export  ==  F|FALSE) {
    return(p)
  } else

    if(file.type == "pdf"){
      if (is.null(file.name)){
        pdf("Eigenvalues.pdf",width = width, height = height)
      } else
      pdf(paste0(file.name, ".pdf"), width = width, height = height)
      plot(p)
      dev.off()
    }

  if (file.type == "tiff"){
    if (is.null(file.name)){
    tiff(filename = "Eigenvalues.tiff",width = width, height = height, units = "in", compression = "lzw", res = resolution)
    } else
    tiff(filename = paste0(file.name, ".tiff"), width = width, height = height, units = "in", compression = "lzw", res = resolution)
      plot(p)
    dev.off()
  }
}
