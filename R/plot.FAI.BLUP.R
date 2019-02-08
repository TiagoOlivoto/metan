plot.FAI.BLUP = function(x, ideotype = 1, SI = 15, radar = TRUE, ...){
  if (!class(x) == "WAASB") {
    stop("The object 'x' is not of class 'WAASB'")
  }
  data = data.frame(FAI = x$FAI[[ideotype]],
                    Genotype = names(x$FAI[[ideotype]]))
  ngsel = round(nrow(data)*(SI/100),0)
  data$sel = "Selected"
  data$sel[(ngsel + 1):nrow(data)] = "Nonselected"
  cutpoint = min(subset(data, sel == "Selected")$MTSI)
  p = ggplot(data = data,  aes(x = reorder(Genotype, FAI), y = FAI)) +
    geom_hline(yintercept = cutpoint, col = "red")+
    geom_path(colour = "black", group = 1) +
    geom_point(size = 2, aes(colour = sel)) +
    scale_x_discrete() +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_line(colour = "gray80"),
          axis.text = element_text(colour = "black")) +
    labs(x = "", y = "FAI-BLUP")
  if (radar == TRUE) {
    p = p + coord_polar()
  }
  return(p)
}
