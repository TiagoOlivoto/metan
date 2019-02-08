plot.MTSI = function(x, SI = 15,
                     radar = TRUE,
                     size.point = 2,
                     col.sel = "red",
                     col.nonsel = "black",
                     ...){
  if (!class(x) == "MTSI") {
    stop("The object 'x' is not of class 'MTSI'")
  }
  data = data.frame(MTSI = x$MTSI,
                     Genotype = names(x$MTSI))
  ngsel = round(nrow(data)*(SI/100),0)
  data$sel = "Selected"
  data$sel[(ngsel + 1):nrow(data)] = "Nonselected"
  cutpoint = max(subset(data, sel == "Selected")$MTSI)
  p = ggplot(data = data,  aes(x= reorder(Genotype, -MTSI), y = MTSI)) +
    geom_hline(yintercept = cutpoint, col = col.sel)+
    geom_path(colour = "black", group = 1) +
    geom_point(size = size.point, aes(colour = sel)) +
    scale_x_discrete() +
    scale_y_reverse()+
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          panel.border = element_blank(),
          axis.text = element_text(colour = "black")) +
    labs(y = "Multitrait stability index") +
    scale_color_manual(values = c(col.nonsel, col.sel))
  if (radar == TRUE) {
    p = p + coord_polar(theta = 'x', start = 0, direction = 1)
  }
  return(p)
}
