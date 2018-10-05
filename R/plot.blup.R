plot.blup  =  function(x,
                      prob  =  0.95,
                      export  =  FALSE,
                      file.type  =  "pdf",
                      file.name = NULL,
                      theme = theme_waasb(),
                      width  =  6,
                      height  =  6,
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
    theme +
    labs(x  =  x.lab, y  =  y.lab) +
    geom_vline(xintercept  =  mean(blup$Predicted)) +
    scale_x_continuous(limits  =  x.lim, breaks  =  x.breaks)

if (export  ==  F|FALSE) {
  return(p1)
} else

if (file.type == "pdf"){
  if (is.null(file.name)){
    pdf("BLUPs genotypes.pdf",width = width, height = height)
  } else
    pdf(paste0(file.name, ".pdf"), width = width, height = height)
  plot(p1)
  dev.off()
  }

if (file.type == "tiff"){
  if (is.null(file.name)){
    tiff(filename = "BLUPs genotypes.tiff", width = width, height = height, units  =  "in", compression  =  "lzw", res = resolution)
  } else
    tiff(filename = paste0(file.name,".tiff"), width = width, height = height, units  =  "in", compression  =  "lzw", res = resolution)
    plot(p1)
    dev.off()
 }
}
