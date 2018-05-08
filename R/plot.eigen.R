plot.eigen = function(data,
                     export = FALSE,
                     file.type = "pdf",
                     width = 6,
                     height = 6,
                     size.lab = 12,
                     size.tex = 12,
                     size.shape = 3.5,
                     size.line = 1,
                     col.shape = c("blue"),
                     col.line = c("blue"),
                     y.lab = "Eigenvalue",
                     x.lab = "Number of multiplicative terms",
                     resolution = 300){
eigen = data$PCA
p = ggplot2::ggplot(eigen, aes(x = PC, y = Eigenvalue, group = 1))  +
    geom_point(stat = 'identity', col = col.shape,  size = size.shape)   +
    geom_line(col = col.line, size = size.line) +
    theme_bw() +
    theme(axis.ticks.length = unit(.2, "cm"),
          axis.text = element_text(size = size.tex, colour = "black"),
          axis.text.x = element_text(size = size.tex, colour = "black"),
          axis.title = element_text(size = size.lab, colour = "black"),
          axis.ticks = element_line(colour = "black"),
          plot.margin = margin(0.2, 0.2, 0.2, 0.7, "cm"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    labs(x = x.lab, y = y.lab)

  if (export  ==  F|FALSE) {
    plot(p)
  } else

    if(file.type == "pdf"){
      pdf("Eigenvalues.pdf",width = width, height = height)
      plot(p)
      dev.off()
    }

  if (file.type == "tiff"){
    tiff(filename = "Eigenvalues.tiff",width = width, height = height, units = "in", compression = "lzw", res = resolution)
    plot(p)
    dev.off()
  }
}
