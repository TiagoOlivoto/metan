plot.scores = function(x,
                     type = 1,
                     polygon = FALSE,
                     file.type = "pdf",
                     export = FALSE,
                     file.name = NULL,
                     theme = theme_waasb(),
                     axis.expand = 1.1,
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
                     size.shape = 2.5,
                     size.tex.lab = 12,
                     size.tex.pa = 3.5,
                     size.line = 0.5,
                     size.segm.line = 0.5,
                     leg.lab = c("Gen", "Env"),
                     line.type = "solid",
                     line.alpha = 0.9,
                     col.line = "gray",
                     col.gen = "blue",
                     col.env = "darkgreen",
                     col.alpha = 0.9,
                     col.segm.gen = "transparent",
                     col.segm.env = "darkgreen",
                     resolution = 300,
                     ...){


  if(polygon == TRUE & type != 1){
    stop("The polygon can be drawn with type 1 graphic only.")
  }

  size.tex.leg = size.tex.pa/0.2917
  class = class(x)
 nenv = nrow(subset(x$model, type == "ENV"))
 ngen = nrow(subset(x$model, type == "GEN"))

if (type == 1){

  if (is.null(y.lab) == F){
    y.lab = y.lab
  } else
    y.lab = paste0("PC2 (", round(x$PCA[2,3],2), "%)")

  if (is.null(x.lab) == F){
    x.lab = x.lab
  } else
    x.lab = paste0("PC1 (", round(x$PCA[1,3],2), "%)")


    if (is.null(x.lim) == F){
    x.lim = x.lim
  } else {
    x.lim = c(min(x$model$PC1 * axis.expand),
              max(x$model$PC1 * axis.expand))
  }

  if (is.null(y.lim) == F){
    y.lim = y.lim
  } else {
    y.lim = c(min(x$model$PC2 * axis.expand),
              max(x$model$PC2 * axis.expand))
  }

p1 = ggplot(x$model, aes(PC1, PC2, shape = type, fill = type))  +
      geom_point(size = size.shape, aes(fill = type), alpha = col.alpha)  +
      scale_shape_manual(labels = leg.lab, values = c(shape.gen,  shape.env))  +
      scale_fill_manual(labels = leg.lab, values = c(col.gen, col.env)) +
      ggrepel::geom_text_repel(aes(PC1, PC2, label = (Code)),
                               size = size.tex.pa, col = c(rep(col.gen, ngen),
                                                        rep(col.env, nenv))) +
      theme %+replace%
      theme(aspect.ratio = 1,
            axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            legend.text = element_text(size = size.tex.leg),
            plot.title = element_text(size = size.tex.lab, hjust = 0, vjust = 1)) +
      labs(x = paste(x.lab), y = paste(y.lab)) +
      scale_x_continuous(limits = x.lim, breaks = x.breaks) +
      scale_y_continuous(limits = y.lim, breaks = y.breaks) +
      geom_vline(xintercept = 0, linetype = line.type, color = col.line, size = size.line, alpha = line.alpha) +
      geom_hline(yintercept = 0, linetype = line.type, color = col.line, size = size.line, alpha = line.alpha) +
      geom_segment(data = x$model, aes(x = 0, y = 0, xend = PC1, yend = PC2, size =  type, color = type, group = type),
                   arrow = arrow(length = unit(0.15, 'cm'))) +
      scale_color_manual(name = "", values = c(col.segm.gen, col.segm.env), theme(legend.position = "none")) +
      scale_size_manual(name = "", values = c(size.segm.line, size.segm.line),theme(legend.position = "none"))

if(polygon == TRUE){

  gen = data.frame(subset(x$model, type == "GEN"))
  coordgenotype = data.frame(subset(x$model, type == "GEN"))[,4:5]
  coordenviroment = data.frame(subset(x$model, type == "ENV"))[,4:5]

  hull <- chull(gen[,4:5])
  indice <- c(hull, hull[1])
  segs <- NULL
  limx <- x.lim
  limy <- y.lim
  i <- 1
  while (is.na(indice[i + 1]) == FALSE) {
    m <- (coordgenotype[indice[i], 2] - coordgenotype[indice[i +
                                                               1], 2])/(coordgenotype[indice[i], 1] - coordgenotype[indice[i +
                                                                                                                             1], 1])
    mperp <- -1/m
    c2 <- coordgenotype[indice[i + 1], 2] - m * coordgenotype[indice[i +
                                                                       1], 1]
    xint <- -c2/(m - mperp)
    xint <- ifelse(xint < 0, min(coordenviroment[, 1],
                                 coordgenotype[, 1]), max(coordenviroment[, 1],
                                                          coordgenotype[, 1]))
    yint <- mperp * xint
    xprop <- ifelse(xint < 0, xint/limx[1], xint/limx[2])
    yprop <- ifelse(yint < 0, yint/limy[1], yint/limy[2])
    m3 <- which(c(xprop, yprop) == max(c(xprop, yprop)))
    m2 <- abs(c(xint, yint)[m3])
    if (m3 == 1 & xint < 0)
      sl1 <- (c(xint, yint)/m2) * abs(limx[1])
    if (m3 == 1 & xint > 0)
      sl1 <- (c(xint, yint)/m2) * limx[2]
    if (m3 == 2 & yint < 0)
      sl1 <- (c(xint, yint)/m2) * abs(limy[1])
    if (m3 == 2 & yint > 0)
      sl1 <- (c(xint, yint)/m2) * limy[2]
    segs <- rbind(segs, sl1)
    i <- i + 1
  }
  rownames(segs) <- NULL
  colnames(segs) <- NULL
  segs <- data.frame(segs)

p1 = p1 + geom_segment(aes(x = X1, y = X2),
               xend = 0,
               yend = 0,
               linetype = 2,
               size = size.segm.line,
               color = "black",
               data = segs,
               inherit.aes = FALSE) +
          geom_polygon(data = gen[indice,],
                       fill = NA,
                       col = col.gen,
                       linetype = 2)


}


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

  if (is.null(x.lim) == F){
    x.lim = x.lim
  } else {
    x.lim = c(min(x$model$Y) - (min(x$model$Y) * axis.expand - min(x$model$Y)),
              max(x$model$Y) + (max(x$model$Y) * axis.expand - max(x$model$Y)))
  }

  if (is.null(y.lim) == F){
    y.lim = y.lim
  } else {
    y.lim = c(min(x$model$PC1 * axis.expand),
              max(x$model$PC1 * axis.expand))
  }
     mean = mean(x$model$Y)
p2 = ggplot2::ggplot(x$model, aes(Y, PC1, shape = type, fill = type))  +
      geom_point(size = size.shape, aes(fill = type), alpha = col.alpha)  +
      scale_shape_manual(labels = leg.lab, values = c(shape.gen,  shape.env))  +
      scale_fill_manual(labels = leg.lab, values = c( col.gen, col.env)) +
      ggrepel::geom_text_repel(aes(Y, PC1, label = (Code)),
                               size = size.tex.pa, col = c(rep(col.gen, ngen),
                                                        rep(col.env, nenv)))  +
      theme %+replace%
      theme(aspect.ratio = 1,
            axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            legend.text = element_text(size = size.tex.leg),
            plot.title = element_text(size = size.tex.lab, hjust = 0, vjust = 1)) +
      labs(x = paste(x.lab), y = paste(y.lab)) +
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

      if (is.null(x.lim) == F){
        x.lim = x.lim
      } else {
        x.lim = c(min(x$model$Y) - (min(x$model$Y) * axis.expand - min(x$model$Y)),
                  max(x$model$Y) + (max(x$model$Y) * axis.expand - max(x$model$Y)))
      }

      if (is.null(y.lim) == F){
        y.lim = y.lim
      } else {
        y.lim = c(min(x$model$WAASB) - (min(x$model$WAASB) * axis.expand - min(x$model$WAASB)),
                  max(x$model$WAASB) + (max(x$model$WAASB) * axis.expand - max(x$model$WAASB)))
      }
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
ggrepel::geom_text_repel(aes(Y, WAASB, label = (Code)),
                         size = size.tex.pa, col = c(rep(col.gen, ngen),
                                                  rep(col.env, nenv)))  +
      theme %+replace%
      theme(aspect.ratio = 1,
            axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            legend.text = element_text(size = size.tex.leg),
            plot.title = element_text(size = size.tex.lab, hjust = 0, vjust = 1)) +
      labs(x = paste(x.lab), y = paste(y.lab)) +
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
  if (is.null(x.lim) == F){
    x.lim = x.lim
  } else {
    x.lim = c(min(x$model$Y) - (min(x$model$Y) * axis.expand - min(x$model$Y)),
              max(x$model$Y) + (max(x$model$Y) * axis.expand - max(x$model$Y)))
  }

  if (is.null(y.lim) == F){
    y.lim = y.lim
  } else {
    y.lim = c(min(x$model$WAAS) - (min(x$model$WAAS) * axis.expand - min(x$model$WAAS)),
              max(x$model$WAAS) + (max(x$model$WAAS) * axis.expand - max(x$model$WAAS)))
  }
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
        ggrepel::geom_text_repel(aes(Y, WAAS, label = (Code)),
                                 size = size.tex.pa, col = c(rep(col.gen, ngen),
                                                          rep(col.env, nenv)))  +
        theme %+replace%
        theme(aspect.ratio = 1,
              axis.text = element_text(size = size.tex.lab, colour = "black"),
              axis.title = element_text(size = size.tex.lab, colour = "black"),
              legend.text = element_text(size = size.tex.leg),
              plot.title = element_text(size = size.tex.lab, hjust = 0, vjust = 1)) +
        labs(x = paste(x.lab), y = paste(y.lab)) +
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

    if (is.null(x.lim) == F){
      x.lim = x.lim
    } else {
      x.lim = c(min(x$MeansGxE$envPC1),
                max(x$MeansGxE$envPC1))
    }

    if (is.null(y.lim) == F){
      y.lim = y.lim
    } else {
      y.lim = c(min(x$MeansGxE$nominal),
                max(x$MeansGxE$nominal))

    }

    min = min(x$MeansGxE$nominal)

    p4 = ggplot2::ggplot(x$MeansGxE, aes(x = envPC1, y = nominal, group = GEN))  +
      geom_line(size = size.line, aes(colour = GEN),
                data = subset(x$MeansGxE, envPC1 %in% c(max(envPC1), min(envPC1))))+
      geom_point(aes(x = envPC1, y = min),
                 data = subset(x$MeansGxE, GEN == x$MeansGxE[1,2])) +
      ggrepel::geom_label_repel(data=subset(x$MeansGxE, envPC1 == min(envPC1)),
                                aes(label = GEN, fill = GEN),
                                size = size.tex.pa, color = 'white',
                                force = 5, segment.color = '#bbbbbb') +
      ggrepel::geom_text_repel(aes(x = envPC1,
                                   y = min,
                                   label = ENV),
                               size = size.tex.pa,
                               force = 5,
                               data = subset(x$MeansGxE, GEN == x$MeansGxE[1,2])) +
      theme %+replace%
      theme(legend.position = "none",
            axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            plot.title = element_text(size = size.tex.lab, hjust = 0, vjust = 1)) +
      scale_x_continuous(limits = x.lim, breaks = x.breaks) +
      scale_y_continuous(limits = y.lim, breaks = y.breaks) +
      labs(x = paste(x.lab), y = y.lab)

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

