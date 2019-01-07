theme_waasb = function () {
  theme_bw() %+replace% # allows the entered values to be overwritten
    theme(axis.ticks.length = unit(.2, "cm"),
          axis.ticks = element_line(colour = "black"),
          legend.position = c(0.85, 0.1), # bottom right
          plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
          legend.title = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())
}
