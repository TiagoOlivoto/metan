plot.blup  =  function(x,
                      prob  =  0.95,
                      export  =  FALSE,
                      file.type  =  "pdf",
                      width  =  6,
                      height  =  6,
                      size.lab  =  12,
                      size.tex  =  12,
                      size.leg  =  12,
                      size.err.bar  =  0.5,
                      size.shape  =  3.5,
                      height.err.bar  =  0.3,
                      x.lim  =  NULL,
                      x.breaks  =  waiver(),
                      leg.pos  =  c(0.9, 0.1),
                      col.shape  =  c("blue", "red"),
                      y.lab  =  "Genotypes",
                      x.lab  =  "Predicted Grain Yield",
                      resolution  =  300,
                      ...){

  PROB  =  ((1 - prob)/2) + prob
  t  =  qt(PROB, 100)
  GV  =  as.numeric(substr(x$ESTIMATES[2,2], start  =  1, stop  =  9))
  AccuGen  =  as.numeric(substr(x$ESTIMATES[8,2], start  =  1, stop  =  9))
  Limits  =  t * sqrt(((1 - AccuGen) * GV))
  blup  =  x$BLUPgen
  blup  =  dplyr::mutate(blup,
                 LL  =  Predicted - Limits,
                 UL  =  Predicted + Limits)
  blup = blup[order(blup$BLUPg), ]
  blup$GEN = factor(blup$GEN, levels  =  blup$GEN)
  blup$Mean = ifelse(blup$Predicted < mean(blup$Predicted), "below", "above")
  blup$CI = blup$Predicted - blup$LL
  mean = mean(blup$Predicted)
p1  =  ggplot2::ggplot(blup, aes(x = Predicted, y = GEN)) +
              geom_point(stat = 'identity', aes(col = Mean),  size  =  size.shape)  +
              geom_errorbarh(aes(
                       xmin = Predicted - CI,
                       xmax = Predicted + CI),
                       size  =  size.err.bar,
                       height  =  height.err.bar) +
    scale_color_manual(name = "Average",
                       values  =  col.shape,
                       labels  =  c("Above", "Below")) +
    theme_bw() +
    theme(axis.ticks.length  =  unit(.2, "cm"),
          axis.text  =  element_text(size  =  size.tex, colour  =  "black"),
          axis.text.x  =  element_text(size  =  size.tex, colour  =  "black"),
          axis.title  =  element_text(size  =  size.lab, colour  =  "black"),
          legend.position  =  leg.pos,
          legend.text  =  element_text(size = size.leg),
          legend.title  =  element_text(size = size.leg),
          axis.ticks  =  element_line(colour  =  "black"),
          plot.margin  =  margin(0.2, 0.2, 0.2, 0.7, "cm"),
          panel.border  =  element_rect(colour  =  "black", fill = NA, size = 1),
          panel.grid.major.x  =  element_blank(), panel.grid.major.y  =  element_blank(),
          panel.grid.minor.x  =  element_blank(), panel.grid.minor.y  =  element_blank()) +
    labs(x  =  x.lab, y  =  y.lab) +
    geom_vline(xintercept  =  mean(blup$Predicted)) +
    scale_x_continuous(limits  =  x.lim, breaks  =  x.breaks)

if (export  ==  F|FALSE) {
  plot(p1)
} else

if (file.type == "pdf"){
  pdf("BLUPs genotypes.pdf",width = width, height = height)
  plot(p1)
  dev.off()
  }

if (file.type == "tiff"){
    tiff(filename = "BLUPs genotypes.tiff",width = width, height = height, units  =  "in", compression  =  "lzw", res = resolution)
    plot(p1)
    dev.off()
 }
}
