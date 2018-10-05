plot.scores = function(x,
                     type,
                     file.type = "pdf",
                     export = FALSE,
                     file.name = NULL,
                     theme = theme_waasb(),
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
                     size.tex = 3.5,
                     size.line = 0.5,
                     size.segm.line = 0.5,
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
      theme +
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
      return(p1)
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
      theme +
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
      return(p2)
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
    I = grid::grobTree(grid::textGrob("I", x = 0.02,  y = 0.98, hjust = 0))
    II = grid::grobTree(grid::textGrob("II", x = 0.97,  y = 0.97, hjust = 0))
    III = grid::grobTree(grid::textGrob("III", x = 0.01,  y = 0.03, hjust = 0))
    IV = grid::grobTree(grid::textGrob("IV", x = 0.96,  y = 0.03, hjust = 0))
p3 = ggplot2::ggplot(x$model, aes(Y, WAASB, shape = type, fill = type))  +
      geom_point(size = size.shape, aes(fill = type), alpha = col.alpha)  +
      scale_shape_manual(labels = leg.lab, values = c(shape.gen,  shape.env))  +
      scale_fill_manual(labels = leg.lab, values = c(col.gen, col.env)) +
ggrepel::geom_text_repel(aes(Y, WAASB, label = (Code)), size = size.tex)  +
     theme +
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
      I = grid::grobTree(grid::textGrob("I", x = 0.02,  y = 0.98, hjust = 0))
      II = grid::grobTree(grid::textGrob("II", x = 0.97,  y = 0.97, hjust = 0))
      III = grid::grobTree(grid::textGrob("III", x = 0.01,  y = 0.03, hjust = 0))
      IV = grid::grobTree(grid::textGrob("IV", x = 0.96,  y = 0.03, hjust = 0))
p3 = ggplot2::ggplot(x$model, aes(Y, WAAS, shape = type, fill = type))  +
        geom_point(size = size.shape, aes(fill = type), alpha = col.alpha)  +
        scale_shape_manual(labels = leg.lab, values = c(shape.gen,  shape.env))  +
        scale_fill_manual(labels = leg.lab, values = c(col.gen, col.env)) +
        ggrepel::geom_text_repel(aes(Y, WAAS, label = (Code)), size = size.tex)  +
        theme +
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
      return(p3)
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
                 data = subset(x$MeansGxE, GEN == x$MeansGxE[1,2])) +
      ggrepel::geom_label_repel(data=subset(x$MeansGxE, envPC1 == min(envPC1)),
                                aes(label = GEN, fill = GEN),
                                size = 3, color = 'white',
                                force = 5, segment.color = '#bbbbbb') +
      ggrepel::geom_text_repel(aes(x = envPC1,
                                   y = min,
                                   label = ENV),
                               force = 5,
                               data = subset(x$MeansGxE, GEN == x$MeansGxE[1,2])) +
      theme +
      theme(legend.position = "none") +
      scale_x_continuous(limits = x.lim, breaks = x.breaks) +
      scale_y_continuous(limits = y.lim, breaks = y.breaks) +
      labs(x = paste("\n", x.lab), y = y.lab)

    if (export  ==  F|FALSE) {
      return(p4)
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

