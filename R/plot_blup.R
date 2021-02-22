#' Plot the BLUPs for genotypes
#' @description
#' `r badge('stable')`
#'
#' Plot the predicted BLUP of the genotypes.
#'
#'
#' @param x The `waasb object`
#' @param var The variable to plot. Defaults to `var = 1` the first
#'   variable of `x`.
#' @param which Which plot to shown. If `which = "gen"` (default) plots the
#'   BLUPs for genotypes. To create a plot showing the BLUPs for
#'   genotype-environment combinations, used `which = "ge"`.
#' @param prob The probability error for constructing confidence interval.
#' @param export Export (or not) the plot. Default is `TRUE`.
#' @param file.type If `export = TRUE`, define the type of file to be
#'   exported. Default is `pdf`, Graphic can also be exported in
#'   `*.tiff` format by declaring `file.type = "tiff"`.
#' @param file.name The name of the file for exportation, default is
#'   `NULL`, i.e. the files are automatically named.
#' @param plot_theme The graphical theme of the plot. Default is
#'   `plot_theme = theme_metan()`. For more details, see
#'   [ggplot2::theme()].
#' @param width The width "inch" of the plot. Default is `6`.
#' @param height The height "inch" of the plot. Default is `6`.
#' @param size.shape The size of the shape (both for genotypes and
#'   environments). Default is `3.5`.
#' @param size.tex.lab The size of the text in axis text and labels.
#' @param err.bar Logical value to indicate if an error bar is shown. Defaults
#'   to `TRUE`.
#' @param size.err.bar The size of the error bar for the plot. Default is
#'   `0.5`.
#' @param height.err.bar The height for error bar. Default is `0.3`.
#' @param x.lim The range of x-axis. Default is `NULL` (maximum and minimum
#'   values of the data set). New arguments can be inserted as `x.lim =
#'   c(x.min, x.max)`.
#' @param x.breaks The breaks to be plotted in the x-axis. Default is
#'   `authomatic breaks`. New arguments can be inserted as `x.breaks =
#'   c(breaks)`
#' @param col.shape A vector of length 2 that contains the color of shapes for
#'   genotypes above and below of the mean, respectively. Default is
#'   `c("blue", "red")`.
#' @param x.lab The label of the x-axis in the plot. Default is `NULL`,
#'   i.e., the name of the selected variable.
#' @param y.lab The label of the y-axis in the plot. Default is
#'   `"Genotypes"`.
#' @param n.dodge The number of rows that should be used to render the Y labels.
#'   This is useful for displaying labels that would otherwise overlap.
#' @param check.overlap Silently remove overlapping labels, (recursively)
#'   prioritizing the first, last, and middle labels.
#' @param panel.spacing Defines the spacing between panels when `which =
#'   "ge"`.
#' @param resolution The resolution of the plot. Parameter valid if
#'   `file.type = "tiff"` is used. Default is `300` (300 dpi)
#' @param ... Currently not used.
#' @return An object of class `gg, ggplot`.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso [plot_scores()], [plot_waasby()]
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' BLUP <- waasb(data_ge2,
#'               resp = PH,
#'               gen = GEN,
#'               env = ENV,
#'               rep = REP)
#' plot_blup(BLUP)
#' plot_blup(BLUP, which = "ge")
#'
#'}
#'
#'
plot_blup <- function(x,
                      var = 1,
                      which = "gen",
                      prob = 0.05,
                      export = FALSE,
                      file.type = "pdf",
                      file.name = NULL,
                      plot_theme = theme_metan(),
                      width = 6,
                      height = 6,
                      err.bar = TRUE,
                      size.err.bar = 0.5,
                      size.shape = 3.5,
                      size.tex.lab = 12,
                      height.err.bar = 0.3,
                      x.lim = NULL,
                      x.breaks = waiver(),
                      col.shape = c("blue", "red"),
                      y.lab = "Genotypes",
                      x.lab = NULL,
                      n.dodge = 1,
                      check.overlap = FALSE,
                      panel.spacing = 0.15,
                      resolution = 300, ...) {
    x.lab <- ifelse(missing(x.lab), names(x[var]), x.lab)
    if(!which %in% c("gen", "ge")){
        stop("Argument 'which' must be one of 'gen' or 'ge'.", call. = FALSE)
    }
    x <- x[[var]]
    if(!class(x)  %in% c("waasb", "gamem")){
        stop("The object 'x' must be of class 'waasb' or 'gamem'.")
    }
    if(class(x) == "gamem"){
        PROB <- ((1 - (1 - prob))/2) + (1 - prob)
        t <- qt(PROB, nlevels(x$residuals$REP))
        GV <- as.numeric(x$ESTIMATES[1, 2])
        AccuGen <- as.numeric(x$ESTIMATES[8, 2])
        Limits <- t * sqrt(((1 - AccuGen) * GV))
        blup <- x$BLUPgen %>%
            mutate(LL = Predicted - Limits,
                   UL = Predicted + Limits,
                   Mean = ifelse(Predicted < mean(Predicted), "below", "above")) %>%
            arrange(Predicted)
    }
    if(class(x) ==  "waasb"){
        PROB <- ((1 - (1 - prob))/2) + (1 - prob)
        t <- qt(PROB, nlevels(x$residuals$REP))
        GV <- as.numeric(x$random[which(x$random$Group == "GEN"), 2])
        AccuGen <- as.numeric(x$ESTIMATES[which(x$ESTIMATES$Parameters == "Accuracy"), 2])
        Limits <- t * sqrt(((1 - AccuGen) * GV))
        if(which == "gen"){
        blup <-
            x$BLUPgen %>%
            mutate(LL = Predicted - Limits,
                   UL = Predicted + Limits,
                   Mean = ifelse(Predicted < mean(Predicted), "below", "above")) %>%
            arrange(Predicted)
        } else{
            blup <-
                x$BLUPint %>%
                means_by(ENV, GEN) %>%
                mutate(LL = Predicted - Limits,
                       UL = Predicted + Limits,
                       Mean = ifelse(Predicted < mean(Predicted), "below", "above"))
        }
    }
    p1 <-
        ggplot(blup, aes(x = Predicted, y = reorder(GEN, Predicted))) +
        geom_vline(xintercept = mean(blup$Predicted), linetype = 2) +
        scale_fill_manual(name = "Average", values = col.shape, labels = c("Above", "Below")) +
        labs(x = x.lab, y = y.lab) +
        scale_x_continuous(limits = x.lim, breaks = x.breaks) +
        scale_y_discrete(guide = guide_axis(n.dodge = n.dodge, check.overlap = check.overlap)) +
        plot_theme %+replace% theme(axis.text = element_text(size = size.tex.lab,colour = "black"),
                               axis.title = element_text(size = size.tex.lab, colour = "black"))
    if(err.bar == TRUE){
        p1 <- p1 +
            geom_errorbarh(aes(xmin = LL, xmax = UL), size = size.err.bar, height = height.err.bar) +
            geom_point(stat = "identity", aes(fill = Mean), shape = 21, size = size.shape)
    }
    if(which == "ge"){
    p1 <- p1 +
        facet_wrap(~ENV) +
        theme(panel.spacing = unit(panel.spacing, "cm"),
              legend.position = "bottom")
    }

    if (export == FALSE) {
        return(p1)
    } else if (file.type == "pdf") {
        if (is.null(file.name)) {
            pdf("BLUPs genotypes.pdf", width = width, height = height)
        } else pdf(paste0(file.name, ".pdf"), width = width, height = height)
        plot(p1)
        dev.off()
    }
    if (file.type == "tiff") {
        if (is.null(file.name)) {
            tiff(filename = "BLUPs genotypes.tiff", width = width, height = height,
                 units = "in", compression = "lzw", res = resolution)
        } else tiff(filename = paste0(file.name, ".tiff"), width = width, height = height,
                    units = "in", compression = "lzw", res = resolution)
        plot(p1)
        dev.off()
    }
}
