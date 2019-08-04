#' Plot the multi-trait stability index
#'
#' Plot the multitrait stability index using ggplot-generated graphics. .
#'
#'
#' @param x An object of class \code{mtsi}
#' @param SI An integer [0-100]. The selection intensity in percentage of the
#' total number of genotypes.
#' @param radar Logical argument. If true (default) a radar plot is generated
#' after using \code{coord_polar()}.
#' @param size.point The size of the point in graphic. Defaults to 2.5.
#' @param col.sel The colour for selected genotypes.
#' @param col.nonsel The colour for nonselected genotypes.
#' @param ... Other arguments of the function.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method plot mtsi
#' @export
#' @examples
#' \dontrun{
#' library(metan)
#' mtsi_model = waasb(data_ge,
#'                    env = ENV,
#'                    gen = GEN,
#'                    rep = REP,
#'                    resp = c(GY, HM))
#'
#' mtsi_index = mtsi(mtsi_model)
#' plot(mtsi_index)
#'
#' # Alternatively using the pipe operator %>%
#' library(dplyr)
#' mtsi_index2 = data_ge %>%
#'               waasb(ENV, GEN, REP, c(GY, HM)) %>%
#'               mtsi()
#' plot(mtsi_index2)
#'}
#'
#'
plot.mtsi <- function(x, SI = 15, radar = TRUE, size.point = 2.5, col.sel = "red", col.nonsel = "black",
    ...) {

    if (!class(x) == "mtsi") {
        stop("The object 'x' is not of class 'mtsi'")
    }
    data <- tibble(MTSI = x$MTSI,
                   Genotype = names(x$MTSI),
                   sel = "Selected")
    data$sel[(round(nrow(data) * (SI/100), 0) + 1):nrow(data)] <- "Nonselected"
    cutpoint <- max(subset(data, sel == "Selected")$MTSI)
    p <- ggplot(data = data, aes(x = reorder(Genotype, -MTSI), y = MTSI)) +
        geom_hline(yintercept = cutpoint, col = col.sel) +
        geom_path(colour = "black", group = 1) +
        geom_point(size = size.point, aes(fill = sel), shape = 21,  colour = "black") +
        scale_x_discrete() +
        scale_y_reverse() +
        theme_minimal()+
        theme(legend.position = "bottom",
              legend.title = element_blank(),
              axis.title.x = element_blank(),
              panel.border = element_blank(),
              axis.text = element_text(colour = "black")) +
        labs(y = "Multitrait stability index") +
        scale_fill_manual(values = c(col.nonsel, col.sel))

    if (radar == TRUE) {

        sequence_length = length(unique(data$Genotype))
        first_sequence = c(1:(sequence_length%/%2))
        second_sequence = c((sequence_length%/%2+1):sequence_length)
        first_angles = c(90 - 180/length(first_sequence) * first_sequence)
        second_angles = c(-90 - 180/length(second_sequence) * second_sequence)
        p <- p + coord_polar() +
            theme(axis.text.x = element_text(angle= c(first_angles, second_angles),
                                             vjust = 1),
                  legend.margin = margin(-120,0,0,0))
    }
    return(p)
}
