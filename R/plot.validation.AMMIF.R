plot.validation.AMMIF = function(data,
                                violin = T,
                                export = F,
                                file.type = "pdf",
                                width = 6,
                                height = 6,
                                resolution = 300,
                                col.violin = "gray90",
                                col.boxplot = "gray70",
                                size.text = 12,
                                size.title = 12,
                                width.boxplot = 0.2,
                                x.lim = NULL,
                                x.breaks = waiver()){

if (violin == T|FALSE){
dodge = position_dodge(width = 1)
p1 = ggplot2::ggplot(data$RMSE, aes(x = MODEL, y = RMSE)) +
     geom_violin(position = dodge, fill = col.violin) +
     geom_boxplot(width = width.boxplot, position = dodge, fill = col.boxplot) +
  theme_bw() +
  theme(axis.ticks.length = unit(.2, "cm"),
        axis.text = element_text(size = size.text, colour = "black"),
        axis.title = element_text(size = size.title, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.background = element_blank(),
        legend.title = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  coord_flip() +
  scale_y_continuous(limits = x.lim, breaks = x.breaks) +
  labs(x = "\nTested models",
       y = expression(paste("Root mean square error (Mg ha"^-1,")")) )
  }else{dodge = position_dodge(width = 1)
p1 = ggplot2::ggplot(data$RMSE, aes(x = MODEL, y = RMSE)) +
      geom_boxplot(width = 0.2, position = dodge, fill = col.boxplot) +
      theme_bw() +
      theme(axis.ticks.length = unit(.2, "cm"),
            axis.text = element_text(size = size.text, colour = "black"),
            axis.title = element_text(size = size.title, colour = "black"),
            axis.ticks = element_line(colour = "black"),
            legend.background = element_blank(),
            legend.title = element_blank(),
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            panel.border = element_rect(colour = "black", fill = NA, size = 1),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank()) +
      coord_flip() +
      scale_y_continuous(limits = x.lim, breaks = x.breaks) +
      labs(x = "\nTested models",
           y = expression(paste("Root mean square error (Mg ha"^-1,")")) )
  }

  if(export  ==  F|FALSE) {
    plot(p1)
  } else

    if(file.type == "pdf"){
      pdf("RMSE validation.pdf",width = width, height = height)
      plot(p1)
      dev.off()
    }

  if (file.type == "tiff"){
    tiff(filename = "RMSE validation.tiff",width = width, height = height, units = "in", compression = "lzw", res = resolution)
    plot(p1)
    dev.off()
  }

}
