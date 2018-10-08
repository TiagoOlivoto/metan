plot.validation.AMMIF = function(x,
                                violin = FALSE,
                                export = FALSE,
                                x.lab = NULL,
                                y.lab = NULL,
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

  if (is.null(y.lab) == F){
    y.lab = y.lab
  } else
    y.lab = expression(paste("Root mean square prediction difference (Mg ha"^-1,")"))

  if (is.null(x.lab) == F){
    x.lab = x.lab
  } else
    x.lab = "AMMI family models"

if (violin == TRUE){
dodge = position_dodge(width = 1)
p1 = ggplot2::ggplot(x$RMSPD, aes(x = MODEL, y = RMSPD)) +
     geom_violin(position = dodge, fill = col.violin) +
     geom_boxplot(width = width.boxplot, position = dodge, fill = col.boxplot) +
     stat_summary(fun.y = mean,
                  geom = "point",
                  shape = 23,
                  fill = "black")+
     theme +
     coord_flip() +
     scale_y_continuous(limits = x.lim, breaks = x.breaks) +
     labs(x = x.lab,
          y = y.lab)
  }else{dodge = position_dodge(width = 1)
p1 = ggplot2::ggplot(x$RMSPD, aes(x = MODEL, y = RMSPD)) +
      geom_boxplot(width = width.boxplot, position = dodge, fill = col.boxplot) +
      stat_summary(fun.y = mean,
                   geom = "point",
                   shape = 23,
                   fill = "black")+
      theme +
      coord_flip() +
      scale_y_continuous(limits = x.lim, breaks = x.breaks) +
      labs(x = x.lab,
           y = y.lab)
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
