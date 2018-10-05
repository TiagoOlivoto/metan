plot.WAASBY = function(x,
                      export = F,
                      file.type = "pdf",
                      file.name = NULL,
                      theme = theme_waasb(),
                      width = 6,
                      height = 6,
                      size.shape = 3.5,
                      col.shape = c("blue", "red"),
                      y.lab = "Genotypes",
                      x.breaks = waiver(),
                      resolution = 300,
                      ...){
  class = class(x)
  if (class != "WAASratio.AMMI" & class != "WAASBYratio") {
    stop("The object 'x' must be a 'WAASratio.AMMI' or a 'WAASBYratio' object.")
  }
  data = x
p1 = ggplot2::ggplot(data$WAASY, aes(x = Code, y = WAASY)) +
    geom_point(stat = 'identity', aes(col = Mean), size = size.shape)  +
    geom_segment(aes(y = min(data$WAASY$WAASY),
                     x = `Code`,
                     yend = WAASY,
                     xend = `Code`),
                     color = "black") +
    coord_flip() +
    scale_color_manual(name = "Average",
                       values = col.shape,
                       labels = c("Above", "Below")) +
    theme +
    scale_y_continuous(limits = c(min(data$WAASY$WAASY), 100), breaks = x.breaks) +
    labs(x = y.lab, y = "\nWAASBY (%)")

  if(export  ==  F|FALSE) {
    return(p1)
  } else

  if(file.type == "pdf"){
    if (is.null(file.name)){
      pdf("WAASY values.pdf",  width = width, height = height)
    } else
      pdf(paste0(file.name, ".pdf"),  width = width, height = height)
      plot(p1)
      dev.off()
    }

  if (file.type == "tiff"){
    if (is.null(file.name)){
      tiff(filename = "WAASY values.tiff",width = width, height = height,
           units = "in", compression = "lzw", res = resolution)
    } else
      tiff(filename = paste0(file.name, ".tiff"), width = width, height = height,
           units = "in", compression = "lzw", res = resolution)
      plot(p1)
      dev.off()
    }
}
