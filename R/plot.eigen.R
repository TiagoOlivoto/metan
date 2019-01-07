plot.eigen = function(x,
                     export = FALSE,
                     theme = theme_waasb(),
                     file.type = "pdf",
                     file.name = NULL,
                     width = 6,
                     height = 6,
                     size.shape = 3.5,
                     size.line = 1,
                     size.tex.lab = 12,
                     y.lab = "Eigenvalue",
                     y2.lab = "Accumulated variance",
                     x.lab = "Number of multiplicative terms",
                     resolution = 300,
                     ...){
  class = class(x)
  if(class != "WAASB"){
    stop("The object 'x' must be a 'WAASB' object.")
  }
eigen = x$PCA
eigen$PC = factor(eigen$PC, levels = eigen$PC)
scaleFactor <- 100/max(eigen$Eigenvalue)

p = ggplot2::ggplot(eigen, aes(x =  PC, group = 1)) +
  geom_line(aes(y = Eigenvalue, col = "Eigenvalue"), size = size.line) +
  geom_point(aes(y = Eigenvalue, col = "Eigenvalue"), size = size.shape) +
  geom_line(aes(y = Accumulated/scaleFactor, col = "Percentage"), size = size.line) +
  geom_point(aes(y = Accumulated/scaleFactor, col = "Percentage"), size = size.shape) +
  scale_y_continuous(sec.axis = sec_axis(~.*scaleFactor,
                                         name = y2.lab)) +
  labs(x = x.lab, y = y.lab) +
  theme %+replace%
  theme(axis.text = element_text(size = size.tex.lab, colour = "black"),
        axis.title = element_text(size = size.tex.lab, colour = "black"),
        legend.position = c(0.15, 0.1))


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
