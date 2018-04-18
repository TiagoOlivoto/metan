plot.WAASBY = function(data,
                      export = F,
                      file.type = "pdf",
                      width = 6,
                      height = 6,
                      size.lab = 12,
                      size.shape = 3.5,
                      legend.pos = c(0.9, 0.1),
                      col.shape = c("blue", "red"),
                      y.lab = "Genotypes",
                      x.breaks = waiver(),
                      resolution = 300)
 {

  p1 = ggplot2::ggplot(data$WAASY, aes(x=Code, y=WAASY)) +
    geom_point(stat='identity', aes(col=Mean), size=size.shape)  +
    geom_segment(aes(y = min(data$WAASY$WAASY),
                     x = `Code`,
                     yend = WAASY,
                     xend = `Code`),
                 color = "black") +
    coord_flip()+
    scale_color_manual(name="Average",
                       values = col.shape,
                       labels = c("Above", "Below")) +
    theme_bw()+
    theme(axis.ticks.length = unit(.2, "cm"),
          axis.text = element_text(size = size.lab, colour = "black"),
          axis.title = element_text(size = size.lab, colour = "black"),
          legend.position = legend.pos,
          axis.ticks = element_line(colour = "black"),
          plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())+
    scale_y_continuous(limits = c(min(data$WAASY$WAASY), 100), breaks = x.breaks)+
    labs(x = y.lab, y = "\nWAASBY (%)")

  if(export == F|FALSE) {
    plot(p1)
  } else

  if(file.type=="pdf"){
      pdf("WAASY values.pdf",  width = width, height = height)
      plot(p1)
      dev.off()
    }

  if (file.type=="tiff"){
      tiff(filename="WAASY values.tiff",width = width, height = height,
           units = "in", compression = "lzw", res = resolution)
      plot(p1)
      dev.off()
    }

  if (file.type=="pptx"){
  doc <- pptx()
  doc <- addSlide(doc, "Two Content" )
  doc <- addPlot(doc, function() print(p1), vector.graphic = TRUE )
  writeDoc(doc, file = "WAASY values.pptx")

  }


}
