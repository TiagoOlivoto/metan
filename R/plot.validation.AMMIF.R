plot.validation.AMMIF = function(x,
                                violin = TRUE,
                                export = FALSE,
                                file.type = "pdf",
                                file.name = NULL,
                                theme = theme_waasb(),
                                width = 6,
                                height = 6,
                                resolution = 300,
                                col.violin = "gray90",
                                col.boxplot = "gray70",
                                width.boxplot = 0.2,
                                x.lim = NULL,
                                x.breaks = waiver(),
                                ...){

if (class(x) != "validation.AMMIF" & class(x) != "validation.blup"){
  warning("The object 'x' should be a 'validation.AMMIF' or a 'validation.blup' object.")
}

if (violin == T|FALSE){
dodge = position_dodge(width = 1)
p1 = ggplot2::ggplot(x$RMSPD, aes(x = MODEL, y = RMSPD)) +
     geom_violin(position = dodge, fill = col.violin) +
     geom_boxplot(width = width.boxplot, position = dodge, fill = col.boxplot) +
     theme +
     coord_flip() +
     scale_y_continuous(limits = x.lim, breaks = x.breaks) +
     labs(x = "\nTested models",
          y = expression(paste("Root mean square prediction difference (Mg ha"^-1,")")) )
  }else{dodge = position_dodge(width = 1)
p1 = ggplot2::ggplot(x$RMSPD, aes(x = MODEL, y = RMSPD)) +
      geom_boxplot(width = 0.2, position = dodge, fill = col.boxplot) +
      theme +
      coord_flip() +
      scale_y_continuous(limits = x.lim, breaks = x.breaks) +
      labs(x = "\nTested models",
           y = expression(paste("Root mean square prediction difference (Mg ha"^-1,")")) )
  }

  if(export  ==  F|FALSE) {
    return(p1)
  } else

    if(file.type == "pdf"){
      if (is.null(file.name)){
      pdf("RMSPD validation.pdf",width = width, height = height)
      } else
      pdf(paste0(file.name, ".pdf"), width = width, height = height)
      plot(p1)
      dev.off()
    }

  if (file.type == "tiff"){
    if (is.null(file.name)){
    tiff(filename = "RMSPD validation.tiff",width = width, height = height, units = "in", compression = "lzw", res = resolution)
    } else
    tiff(filename = paste0(file.name, ".tiff"), width = width, height = height, units = "in", compression = "lzw", res = resolution)
    plot(p1)
    dev.off()
  }

}
