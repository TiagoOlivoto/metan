plot.scores = function(x,
                     type,
                     file.type = "pdf",
                     export = FALSE,
                     file.name = NULL,
                     width = 8,
                     height = 7,
                     x.lim = NULL,
                     x.breaks = waiver(),
                     x.lab = NULL,
                     y.lab = NULL,
                     y.lim = NULL,
                     y.breaks = waiver(),
                     shape.gen = 21,
                     shape.env = 23,
                     size.shape = 3.5,
                     size.lab = 18,
                     size.tex = 3.5,
                     size.line = 0.5,
                     size.leg = 14,
                     size.segm.line = 0.5,
                     leg.position = "tr",
                     leg.lab = c("Gen", "Env"),
                     line.type = "dashed",
                     line.alpha = 0.9,
                     col.line = "gray",
                     col.gen = "gray75",
                     col.env = "red",
                     col.alpha = 0.9,
                     col.segm.gen = "transparent",
                     col.segm.env = "grey50",
                     resolution = 300,
                     ...){

  if (leg.position  ==  "tr"){
    leg.pos = c(.87, 0.9) #lt
  }

  if (leg.position  ==  "tl"){
    leg.pos = c(.12, 0.9) #lt
  }

  if (leg.position  ==  "br"){
    leg.pos = c(0.87, 0.12)
  }


  if (leg.position  ==  "bl"){
    leg.pos = c(0.12, 0.12)
  }

  class = class(x)


if (type == 1){

  if (is.null(y.lab) == F){
    y.lab = y.lab
  } else
    y.lab = paste0("PC2 (", round(x$PCA[2,3],2), "%)")

  if (is.null(x.lab) == F){
    x.lab = x.lab
  } else
    x.lab = paste0("\nPC1 (", round(x$PCA[1,3],2), "%)")



p1 = ggplot2::ggplot(x$model, aes(PC1, PC2, shape = type, fill = type))  +
      geom_point(size = size.shape, aes(fill = type), alpha = col.alpha)  +
      scale_shape_manual(labels = leg.lab, values = c(shape.gen,  shape.env))  +
      scale_fill_manual(labels = leg.lab, values = c( col.gen, col.env)) +
      ggrepel::geom_text_repel(aes(PC1, PC2, label = (Code)), size = size.tex)  +
      theme_bw() +
      theme(axis.ticks.length = unit(.2, "cm"),
            axis.text = element_text(size = size.lab, colour = "black"),
            axis.title = element_text(size = size.lab, colour = "black"),
            axis.ticks = element_line(colour = "black"),
            legend.position = leg.pos,
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            legend.title = element_blank(),
            legend.text = element_text(size = size.leg),
            panel.border = element_rect(colour = "black", fill = NA, size = 1),
            panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
            panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
      labs(x = paste("\n", x.lab), y = paste(y.lab)) +
      scale_x_continuous(limits = x.lim, breaks = x.breaks) +
      scale_y_continuous(limits = y.lim, breaks = y.breaks) +
      geom_vline(xintercept = 0, linetype = line.type, color = col.line, size = size.line, alpha = line.alpha) +
      geom_hline(yintercept = 0, linetype = line.type, color = col.line, size = size.line, alpha = line.alpha) +
      geom_segment(data = x$WAAS, aes(x = 0, y = 0, xend = PC1, yend = PC2, size =  type, color = type, group = type),
                   arrow = arrow(length = unit(0.15, 'cm'))) +
      scale_color_manual(name = "", values = c( col.segm.gen, col.segm.env), theme(legend.position = "none")) +
      scale_size_manual(name = "", values = c(size.segm.line, size.segm.line),theme(legend.position = "none"))

    if (export  ==  F|FALSE) {
      plot(p1)
    } else

      if (file.type == "pdf"){
        if (is.null(file.name)){
        pdf("PC1 x PC2.pdf",width = width, height = height)
        } else
        pdf(paste0(file.name, ".pdf"), width = width, height = height)
        plot(p1)
        dev.off()
      }

    if (file.type == "tiff"){
      if (is.null(file.name)){
      tiff(filename = "PC1 x PC2.tiff",width = width, height = height, units = "in", compression = "lzw", res = resolution)
      } else
      tiff(filename = paste0(file.name, ".tiff"), width = width, height = height, units = "in", compression = "lzw", res = resolution)
      plot(p1)
      dev.off()
    }

  }

if (type == 2){

  if (is.null(y.lab) == F){
    y.lab = y.lab
  } else
    y.lab = paste0("PC1 (", round(x$PCA[1,3],2), "%)")

  if (is.null(x.lab) == F){
    x.lab = x.lab
  } else
    x.lab = paste0("Grain yield")

     mean = mean(x$model$Y)
p2 = ggplot2::ggplot(x$model, aes(Y, PC1, shape = type, fill = type))  +
      geom_point(size = size.shape, aes(fill = type), alpha = col.alpha)  +
      scale_shape_manual(labels = leg.lab, values = c(shape.gen,  shape.env))  +
      scale_fill_manual(labels = leg.lab, values = c( col.gen, col.env)) +
      ggrepel::geom_text_repel(aes(Y, PC1, label = (Code)), size = size.tex)  +
      theme_bw() +
      theme(axis.ticks.length = unit(.2, "cm"),
            axis.text = element_text(size = size.lab, colour = "black"),
            axis.title = element_text(size = size.lab, colour = "black"),
            axis.ticks = element_line(colour = "black"),
            legend.position = leg.pos,
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            legend.title = element_blank(),
            legend.text = element_text(size = size.leg),
            panel.border = element_rect(colour = "black", fill = NA, size = 1),
            panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
            panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
      labs(x = paste("\n", x.lab), y = paste(y.lab)) +
      scale_x_continuous(limits = x.lim, breaks = x.breaks) +
      scale_y_continuous(limits = y.lim, breaks = y.breaks) +
      geom_vline(xintercept = mean(x$model$Y), linetype = line.type, color = col.line, size = size.line, alpha = line.alpha) +
      geom_hline(yintercept = 0, linetype = line.type, size = size.line, color = col.line, alpha = line.alpha) +
      geom_segment(data = x$model, aes(x = mean, y = 0, xend = Y, yend = PC1, size = type, color = type, group = type),
                   arrow = arrow(length = unit(0.15, 'cm'))) +
      scale_color_manual(name = "", values = c( col.segm.gen, col.segm.env), theme(legend.position = "none")) +
      scale_size_manual(name = "", values = c(size.segm.line, size.segm.line),theme(legend.position = "none"))

    if (export  ==  F|FALSE) {
      plot(p2)
    } else

      if (file.type == "pdf"){
        if (is.null(file.name)){
          pdf("GY x PC1.pdf",width = width, height = height)
        } else
          pdf(paste0(file.name, ".pdf"), width = width, height = height)
        plot(p2)
        dev.off()

      }

    if (file.type == "tiff"){
      if (is.null(file.name)){
        tiff(filename = "GY x PC1.tiff",width = width, height = height, units = "in", compression = "lzw", res = resolution)
      } else
        tiff(filename = paste0(file.name, ".tiff"), width = width, height = height, units = "in", compression = "lzw", res = resolution)
      plot(p2)
      dev.off()
    }
}

if (type == 3){

  if (is.null(y.lab) == F){
    y.lab = y.lab
  } else
    y.lab = paste0("Weighted average of the absolute scores")

  if (is.null(x.lab) == F){
    x.lab = x.lab
  } else
    x.lab = paste0("Grain yield")

    if (class  ==  "WAASB"){
    m1 = mean(x$model$Y)
    m2 = mean(x$model$WAASB)
    I = grid::grobTree(textGrob("I", x = 0.02,  y = 0.98, hjust = 0))
    II = grid::grobTree(textGrob("II", x = 0.97,  y = 0.97, hjust = 0))
    III = grid::grobTree(textGrob("III", x = 0.01,  y = 0.03, hjust = 0))
    IV = grid::grobTree(textGrob("IV", x = 0.96,  y = 0.03, hjust = 0))
p3 = ggplot2::ggplot(x$model, aes(Y, WAASB, shape = type, fill = type))  +
      geom_point(size = size.shape, aes(fill = type), alpha = col.alpha)  +
      scale_shape_manual(labels = leg.lab, values = c(shape.gen,  shape.env))  +
      scale_fill_manual(labels = leg.lab, values = c(col.gen, col.env)) +
ggrepel::geom_text_repel(aes(Y, WAASB, label = (Code)), size = size.tex)  +
      theme_bw() +
      theme(axis.ticks.length = unit(.2, "cm"),
            axis.text = element_text(size = size.lab, colour = "black"),
            axis.title = element_text(size = size.lab, colour = "black"),
            axis.ticks = element_line(colour = "black"),
            legend.position = leg.pos,
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            axis.title.y = element_text(margin = margin(r = 16)),
            legend.title = element_blank(),
            legend.text = element_text(size = size.leg),
            panel.border = element_rect(colour = "black", fill = NA, size = 1),
            panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
            panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
      labs(x = paste("\n", x.lab), y = paste(y.lab)) +
      scale_x_continuous(limits = x.lim, breaks = x.breaks) +
      scale_y_continuous(limits = y.lim, breaks = y.breaks) +
      geom_vline(xintercept = m1, linetype = line.type, color = col.line, size = size.line, alpha = line.alpha) +
      geom_hline(yintercept = m2, linetype = line.type, color = col.line, size = size.line, alpha = line.alpha) +
      annotation_custom(I) +
      annotation_custom(II) +
      annotation_custom(III) +
      annotation_custom(IV)
    }

if (class  ==  "WAAS.AMMI"){
      m1 = mean(x$model$Y)
      m2 = mean(x$model$WAAS)
      I = grid::grobTree(textGrob("I", x = 0.02,  y = 0.98, hjust = 0))
      II = grid::grobTree(textGrob("II", x = 0.97,  y = 0.97, hjust = 0))
      III = grid::grobTree(textGrob("III", x = 0.01,  y = 0.03, hjust = 0))
      IV = grid::grobTree(textGrob("IV", x = 0.96,  y = 0.03, hjust = 0))
p3 = ggplot2::ggplot(x$model, aes(Y, WAAS, shape = type, fill = type))  +
        geom_point(size = size.shape, aes(fill = type), alpha = col.alpha)  +
        scale_shape_manual(labels = leg.lab, values = c(shape.gen,  shape.env))  +
        scale_fill_manual(labels = leg.lab, values = c(col.gen, col.env)) +
        ggrepel::geom_text_repel(aes(Y, WAAS, label = (Code)), size = size.tex)  +
        theme_bw() +
        theme(axis.ticks.length = unit(.2, "cm"),
              axis.text = element_text(size = size.lab, colour = "black"),
              axis.title = element_text(size = size.lab, colour = "black"),
              axis.ticks = element_line(colour = "black"),
              legend.position = leg.pos,
              plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
              axis.title.y = element_text(margin = margin(r = 16)),
              legend.title = element_blank(),
              legend.text = element_text(size = size.leg),
              panel.border = element_rect(colour = "black", fill = NA, size = 1),
              panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
              panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
        labs(x = paste("\n", x.lab), y = paste(y.lab)) +
        scale_x_continuous(limits = x.lim, breaks = x.breaks) +
        scale_y_continuous(limits = y.lim, breaks = y.breaks) +
        geom_vline(xintercept = m1, linetype = line.type, color = col.line, size = size.line, alpha = line.alpha) +
        geom_hline(yintercept = m2, linetype = line.type, color = col.line, size = size.line, alpha = line.alpha) +
        annotation_custom(I) +
        annotation_custom(II) +
        annotation_custom(III) +
        annotation_custom(IV)
    }

    if (export  ==  F|FALSE) {
      plot(p3)
    } else

      if (file.type == "pdf"){
        if (is.null(file.name)){
          pdf("GY x WAASB.pdf",width = width, height = height)
        } else
          pdf(paste0(file.name, ".pdf"), width = width, height = height)
        plot(p3)
        dev.off()
      }

    if (file.type == "tiff"){
      if (is.null(file.name)){
        tiff(filename = "GY x WAAS.tiff",width = width, height = height, units = "in", compression = "lzw", res = resolution)
      } else
        tiff(filename = paste0(file.name, ".tiff"), width = width, height = height, units = "in", compression = "lzw", res = resolution)
      plot(p3)
      dev.off()
     }
}

  if (type == 4){

    if (is.null(y.lab) == F){
      y.lab = y.lab
    } else
      y.lab = paste0("Nominal Yield (Mg/ha)")

    if (is.null(x.lab) == F){
      x.lab = x.lab
    } else
      x.lab = paste0("Environment PC1 [square root of  (Mg/ha)]")


    min = min(x$MeansGxE$nominal)

    p4 = ggplot2::ggplot(x$MeansGxE, aes(x = envPC1, y = nominal, group = GEN))  +
      geom_line(size = 1, aes(colour = GEN),
                data = subset(x$MeansGxE, envPC1 %in% c(max(envPC1), min(envPC1))))+
      geom_point(aes(x = envPC1, y = min),
                 data = subset(x$MeansGxE, GEN == x$MeansGxE[1,2]))+
      ggrepel::geom_label_repel(data=subset(x$MeansGxE, envPC1 == min(envPC1)),
                                aes(label=GEN, fill=GEN),
                                size=3, color='white',
                                force=5, segment.color='#bbbbbb') +
      ggrepel::geom_text_repel(aes(x = envPC1,
                                   y = min,
                                   label = ENV),
                               force = 5,
                               data = subset(x$MeansGxE, GEN == x$MeansGxE[1,2])) +
      theme_bw() +
      theme(axis.ticks.length = unit(.2, "cm"),
            axis.text = element_text(size = size.lab, colour = "black"),
            axis.title = element_text(size = size.lab, colour = "black"),
            axis.ticks = element_line(colour = "black"),
            legend.position = "none",
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            panel.border = element_rect(colour = "black", fill = NA, size = 1),
            panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
            panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
      scale_x_continuous(limits = x.lim, breaks = x.breaks) +
      scale_y_continuous(limits = y.lim, breaks = y.breaks) +
      labs(x = x.lab, y = y.lab)

    if (export  ==  F|FALSE) {
      plot(p4)
    } else

      if (file.type == "pdf"){
        if (is.null(file.name)){
          pdf("Adaptative reponses (AMMI1 model).pdf",width = width, height = height)
        } else
          pdf(paste0(file.name, ".pdf"), width = width, height = height)
        plot(p4)
        dev.off()

      }

    if (file.type == "tiff"){
      if (is.null(file.name)){
        tiff(filename = "Adaptative reponses (AMMI1 model).tiff",width = width, height = height, units = "in", compression = "lzw", res = resolution)
      } else
        tiff(filename = paste0(file.name, ".tiff"), width = width, height = height, units = "in", compression = "lzw", res = resolution)
      plot(p4)
      dev.off()
    }
  }

}

