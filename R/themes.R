#' Personalized theme for ggplot2-based graphics
#' @rdname themes
#' @description Two themes that provide plots with a gray background and major
#'   grids (\code{theme_metan}) or minimalistic theme with half-open frame, white
#'   background, and no grid (\code{theme_metan_minimal}). For more details see
#'   \code{\link[ggplot2]{theme}}.
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'
theme_metan = function () {
  theme_gray() %+replace% # allows the entered values to be overwritten
    theme(axis.ticks.length = unit(.2, "cm"),
          plot.title = element_text(face = "bold", hjust = 0, vjust = 3),
          plot.subtitle = element_text(face = "italic", hjust = 0, vjust = 2, size = 8),
          axis.ticks = element_line(colour = "black"),
          legend.position = c(0.85, 0.1),
          legend.key = element_rect(fill = NA, colour = "transparent"),
          legend.background = element_rect(fill = NA, colour = "transparent"),
          plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm"),
          legend.title = element_blank(),
          axis.text = element_text(colour = "black"),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
}

#' @rdname themes
#' @export
#'
theme_metan_minimal = function () {
  theme_gray() %+replace% # allows the entered values to be overwritten
    theme(axis.ticks.length = unit(.2, "cm"),
          plot.title = element_text(face = "bold", hjust = 0, vjust = 3),
          plot.subtitle = element_text(face = "italic", hjust = 0, vjust = 2, size = 8),
          axis.ticks = element_line(colour = "black"),
          legend.position = c(0.85, 0.1),
          legend.key = element_rect(fill = NA, colour = "transparent"),
          legend.background = element_rect(fill = NA, colour = "transparent"),
          plot.margin = margin(0.3, 0.1, 0.1, 0.1, "cm"),
          legend.title = element_blank(),
          axis.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line.x.bottom = element_line(),
          axis.line.y.left = element_line(),
          panel.border = element_blank())
}
