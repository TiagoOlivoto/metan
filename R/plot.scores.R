#' @method
#' @export
plot.scores = function(data,
                     type,
                     file.type = "pdf",
                     export = FALSE,
                     width = 8,
                     height = 7,
                     x.lim = NULL,
                     x.breaks = waiver(),
                     x.lab = "Grain yield",
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
                     resolution = 300){


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


if (type == 1){

p1 = ggplot2::ggplot(data$model, aes(PC1, PC2, shape = type, fill = type))  +
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
      labs(x = paste0("\nPC1 (", round(data$PCA[1,3],2), "%)"),
           y = paste0("PC2 (", round(data$PCA[2,3],2), "%)")) +
      scale_x_continuous(limits = x.lim, breaks = x.breaks) +
      scale_y_continuous(limits = y.lim, breaks = y.breaks) +
      geom_vline(xintercept = 0, linetype = line.type, color = col.line, size = size.line, alpha = line.alpha) +
      geom_hline(yintercept = 0, linetype = line.type, color = col.line, size = size.line, alpha = line.alpha) +
      geom_segment(data = data$WAAS, aes(x = 0, y = 0, xend = PC1, yend = PC2, size =  type, color = type, group = type),
                   arrow = arrow(length = unit(0.15, 'cm'))) +
      scale_color_manual(name = "", values = c( col.segm.gen, col.segm.env), theme(legend.position = "none")) +
      scale_size_manual(name = "", values = c(size.segm.line, size.segm.line),theme(legend.position = "none"))

    if (export  ==  F|FALSE) {
      plot(p1)
    } else

      if (file.type == "pdf"){
        pdf("PC1 x PC2.pdf",width = width, height = height)
        plot(p1)
        dev.off()
      }

    if (file.type == "tiff"){
      tiff(filename = "PC1 x PC2.tiff",width = width, height = height, units = "in", compression = "lzw", res = resolution)
      plot(p1)
      dev.off()
    }

  }

if (type == 2){
    mean = mean(data$model$Y)
p2 = ggplot2::ggplot(data$model, aes(Y, PC1, shape = type, fill = type))  +
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
      labs(x = paste("\n", x.lab), y = paste0("PC1 (", round(data$PCA[1,3],2), "%)")) +
      scale_x_continuous(limits = x.lim, breaks = x.breaks) +
      scale_y_continuous(limits = y.lim, breaks = y.breaks) +
      geom_vline(xintercept = mean(data$model$Y), linetype = line.type, color = col.line, size = size.line, alpha = line.alpha) +
      geom_hline(yintercept = 0, linetype = line.type, size = size.line, color = col.line, alpha = line.alpha) +
      geom_segment(data = data$model, aes(x = mean, y = 0, xend = Y, yend = PC1, size = type, color = type, group = type),
                   arrow = arrow(length = unit(0.15, 'cm'))) +
      scale_color_manual(name = "", values = c( col.segm.gen, col.segm.env), theme(legend.position = "none")) +
      scale_size_manual(name = "", values = c(size.segm.line, size.segm.line),theme(legend.position = "none"))

    if (export  ==  F|FALSE) {
      plot(p2)
    } else

      if (file.type == "pdf"){
        pdf("GY x PC1.pdf",width = width, height = height)
        plot(p2)
        dev.off()
      }

    if (file.type == "tiff"){
      tiff(filename = "GY x PC1.tiff",width = width, height = height, units = "in", compression = "lzw", res = resolution)
      plot(p2)
      dev.off()
    }

  }


if (type == 3){
    if (data$object  ==  "WAASB"){
    m1 = mean(data$model$Y)
    m2 = mean(data$model$WAASB)
    I = grid::grobTree(textGrob("I", x = 0.02,  y = 0.98, hjust = 0))
    II = grid::grobTree(textGrob("II", x = 0.97,  y = 0.97, hjust = 0))
    III = grid::grobTree(textGrob("III", x = 0.01,  y = 0.03, hjust = 0))
    IV = grid::grobTree(textGrob("IV", x = 0.96,  y = 0.03, hjust = 0))
    p3 = ggplot2::ggplot(data$model, aes(Y, WAASB, shape = type, fill = type))  +
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
      labs(x = paste("\n", x.lab), y = "Weighted average of the absolute scores") +
      scale_x_continuous(limits = x.lim, breaks = x.breaks) +
      scale_y_continuous(limits = y.lim, breaks = y.breaks) +
      geom_vline(xintercept = m1, linetype = line.type, color = col.line, size = size.line, alpha = line.alpha) +
      geom_hline(yintercept = m2, linetype = line.type, color = col.line, size = size.line, alpha = line.alpha) +
      annotation_custom(I) +
      annotation_custom(II) +
      annotation_custom(III) +
      annotation_custom(IV)
    }

if (data$object  ==  "WAAS"){
      m1 = mean(data$model$Y)
      m2 = mean(data$model$WAAS)
      I = grid::grobTree(textGrob("I", x = 0.02,  y = 0.98, hjust = 0))
      II = grid::grobTree(textGrob("II", x = 0.97,  y = 0.97, hjust = 0))
      III = grid::grobTree(textGrob("III", x = 0.01,  y = 0.03, hjust = 0))
      IV = grid::grobTree(textGrob("IV", x = 0.96,  y = 0.03, hjust = 0))
p3 = ggplot2::ggplot(data$model, aes(Y, WAAS, shape = type, fill = type))  +
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
        labs(x = paste("\n", x.lab), y = "Weighted average of the absolute scores") +
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
        pdf("GY x WAASB.pdf",width = width, height = height)
        plot(p3)
        dev.off()
      }

    if (file.type == "tiff"){
      tiff(filename = "GY x WAAS.tiff",width = width, height = height, units = "in", compression = "lzw", res = resolution)
      plot(p3)
      dev.off()
    }
  }
}

