#' Plot the multi-trait stability index
#'
#' Plot the multitrait stability index using ggplot-generated graphics. .
#'
#'
#' @param x An object of class \code{MTSI}
#' @param SI An integer [0-100]. The selection intensity in percentage of the
#' total number of genotypes.
#' @param radar Logical argument. If true (default) a radar plot is generated
#' after using \code{coord_polar()}.
#' @param size.point The size of the point in graphic.
#' @param col.sel The colour for selected genotypes.
#' @param col.nonsel The colour for nonselected genotypes.
#' @param ... Other arguments of the function.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method plot MTSI
#' @export
#' @examples
#'
#' library(metan)
#' MTSI_MODEL = WAASB(data_ge,
#'                    env = ENV,
#'                    gen = GEN,
#'                    rep = REP,
#'                    resp = c(GY, HM))
#'
#' MTSI_index = MTSI(MTSI_MODEL)
#' plot(MTSI_index)
#'
#' # Alternatively using the pipe operator %>%
#' library(dplyr)
#' MTSI_index2 = data_ge %>%
#'               WAASB(ENV, GEN, REP, c(GY, HM)) %>%
#'               MTSI()
#' plot(MTSI_index2)
#'
#'
#'
plot.MTSI <- function(x, SI = 15, radar = TRUE, size.point = 2, col.sel = "red", col.nonsel = "black",
    ...) {

    if (!class(x) == "MTSI") {
        stop("The object 'x' is not of class 'MTSI'")
    }
    data <- data.frame(MTSI = x$MTSI, Genotype = names(x$MTSI))
    ngsel <- round(nrow(data) * (SI/100), 0)
    data$sel <- "Selected"
    data$sel[(ngsel + 1):nrow(data)] <- "Nonselected"
    cutpoint <- max(subset(data, sel == "Selected")$MTSI)
    p <- ggplot(data = data, aes(x = reorder(Genotype, -MTSI), y = MTSI)) + geom_hline(yintercept = cutpoint,
        col = col.sel) + geom_path(colour = "black", group = 1) + geom_point(size = size.point,
        aes(colour = sel)) + scale_x_discrete() + scale_y_reverse() + theme_minimal() +
        theme(legend.position = "bottom", legend.title = element_blank(), axis.title.x = element_blank(),
            panel.border = element_blank(), axis.text = element_text(colour = "black")) +
        labs(y = "Multitrait stability index") + scale_color_manual(values = c(col.nonsel,
        col.sel))
    if (radar == TRUE) {
        p <- p + coord_polar(theta = "x", start = 0, direction = 1)
    }
    return(p)
}
