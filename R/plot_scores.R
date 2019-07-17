#' Plot scores in different graphical interpretations
#'
#'
#' Plot scores of genotypes and environments in different graphics. \code{1 =
#' PC1 x PC2}, \code{2 = GY x PC1}, \code{3 = GY x WAASB}, and \code{4 =
#' Nominal yield x EPCA1}.
#'
#' The plots type 1 and 2 have the same interpretation than those used in
#' traditional-usage AMMI analysis (well know as AMMI2 and AMMI1,
#' respectively). In the plot type 3, the scores of both genotypes and
#' environments are plotted considering the response variable and the WAASB
#' (stability index that considers all significant principal component axis of
#' traditional AMMI models or all principal component axis estimated with
#' BLUP-interaction effects. Plot type 4 may be used to better understand the
#' well known 'which-won-where' pattern, facilitating the recommendation of
#' appropriate genotypes targeted for specific environments, thus allowing the
#' exploitation of narrow adaptations.
#'
#' @param x The object \code{WAASB} or \code{WAAS.AMMI}
#' @param type Three types of graphics can be generated: \code{1 = PC1 x PC2},
#' default, to make inferences related to the interaction effects; \code{2 = GY
#' x PC1} to make inferences related to stability and productivity; \code{3 =
#' GY x WAASB}, and \code{4 = Nominal yield x Environment PC1}.
#' @param polygon Logical argument. If \code{TRUE}, a polygon is drawn when
#' \code{type 1}.
#' @param file.type The type of file to be exported. Valid parameter if
#' \code{export = T|TRUE}.  Default is \code{'pdf'}. The graphic can also be
#' exported in \code{*.tiff} format by declaring \code{file.type = 'tiff'}.
#' @param export Export (or not) the plot. Default is \code{FALSE}.
#' @param file.name The name of the file for exportation, default is
#' \code{NULL}, i.e. the files are automatically named.
#' @param theme The graphical theme of the plot. Default is `theme =
#' theme_waasb()`. Please, see `?metan::theme_waasb`. An own theme can be
#' applied using the arguments: `theme = theme_waasb() + theme(some stuff
#' here)`. For more details, please, see `?ggplot2::theme`
#' @param axis.expand Multiplication factor to expand the axis limits by to
#' enable fitting of labels. Default is \code{1.1}.
#' @param width The width 'inch' of the plot. Default is \code{8}.
#' @param height The height 'inch' of the plot. Default is \code{7}.
#' @param x.lim The range of x-axis. Default is \code{NULL} (maximum and
#' minimum values of the data set). New arguments can be inserted as
#' \code{x.lim = c(x.min, x.max)}.
#' @param x.breaks The breaks to be plotted in the x-axis. Default is
#' \code{authomatic breaks}. New arguments can be inserted as \code{x.breaks =
#' c(breaks)}
#' @param x.lab The label of x-axis. Each plot has a default value. New
#' arguments can be inserted as \code{x.lab = 'my label'}.
#' @param y.lab The label of y-axis. Each plot has a default value. New
#' arguments can be inserted as \code{y.lab = 'my label'}.
#' @param y.lim The range of y-axis. Default is \code{NULL}. The same arguments
#' than \code{x.lim} can be used.
#' @param y.breaks The breaks to be plotted in the x-axis. Default is
#' \code{authomatic breaks}. The same arguments than \code{x.breaks} can be
#' used.
#' @param shape.gen The shape for genotype indication in the biplot. Default is
#' \code{21} (circle). Values must be between \code{21-25}: \code{21} (circle),
#' \code{22} (square), \code{23} (diamond), \code{24} (up triangle), and
#' \code{25} (low triangle).
#' @param shape.env The shape for environment indication in the biplot. Default
#' is \code{23} (diamond). The same arguments than \code{'shape.gen'}.
#' @param size.shape The size of the shape (both for genotypes and
#' environments). Default is \code{2.2}.
#' @param size.bor.tick The size of tick of shape. Default is \code{0.3}. The
#' size of the shape will be \code{size.shape + size.bor.tick}
#' @param size.tex.lab The size of the text in the axes text and labels.
#' Default is \code{12}.
#' @param size.tex.pa The size of the text of the plot area. Default is
#' \code{3.5}.
#' @param repulsion Force of repulsion between overlapping text labels. Defaults to 1.
#' @param size.line The size of the line that indicate the means in the biplot.
#' Default is \code{0.5}.
#' @param size.segm.line The size of the segment that start in the origin of
#' the biplot and end in the scores values. Default is \code{0.5}.
#' @param leg.lab The labs of legend. Default is \code{Gen} and \code{Env}.
#' @param line.type The type of the line that indicate the means in the biplot.
#' Default is \code{'solid'}. Other values that can be attributed are:
#' \code{'blank'}, no lines in the biplot, \code{'dashed', 'dotted', 'dotdash',
#' 'longdash', and 'twodash'}.
#' @param line.alpha The alpha value that combine the line with the background
#' to create the appearance of partial or full transparency. Default is
#' \code{0.4}. Values must be between '0' (full transparency) to '1' (full
#' color).
#' @param col.line The color of the line that indicate the means in the biplot.
#' Default is \code{'gray'}
#' @param col.gen The shape color for genotypes. Must be one value or a vector
#' of colours with the same length of the number of genotypes. Default is
#' \code{'orange'}. Other values can be attributed. For example,
#' \code{'transparent'}, will make a plot with only an outline around the shape
#' area.
#' @param col.env The shape color for environments. Default is \code{'forestgreen'}.
#' The same usability than \code{'col.gen'}.
#' @param col.alpha The alpha value for the color. Default is \code{0.9}.
#' Values must be between \code{0} (full transparency) to \code{1} (full
#' color).
#' @param col.segm.gen The color of segment for genotypes.Default is
#' \code{'transparent'}. Parameter valid for \code{type = 1} and \code{type =
#' 2} graphics. This segment start in the origin of the biplot and end in the
#' scores values.
#' @param col.segm.env The color of segment for environments. Default is
#' \code{'forestgreen'}. The same usability than \code{'col.segm.gen'}
#' @param resolution The resolution of the plot. Parameter valid if
#' \code{file.type = 'tiff'} is used. Default is \code{300} (300 dpi)
#' @param ... Other arguments of the function
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso \code{\link{plot_eigen}}
#' @importFrom ggforce geom_arc
#' @export
#' @examples
#'
#' library(metan)
#' library(ggplot2)
#' scores = waasb(data_ge,
#'                resp = GY,
#'                gen = GEN,
#'                env = ENV,
#'                rep = REP)
#' # PC1 x PC2
#' plot_scores(scores$GY,
#'             type = 1,
#'             polygon = TRUE)
#'
#' # GY x PC1
#' plot_scores(scores$GY,
#'             type = 2,
#'             col.env = 'olivedrab',
#'             col.gen = 'orange2',
#'             x.lab = 'My own x label')
#'
#' # GY x WAASB
#' plot_scores(scores$GY,
#'             type = 3,
#'             size.tex.pa = 2,
#'             size.tex.lab = 16)
#'
#'
plot_scores <- function(x, type = 1, polygon = FALSE, file.type = "pdf",
export = FALSE, file.name = NULL, theme = theme_waasb(),
axis.expand = 1.1, width = 8, height = 7, x.lim = NULL, x.breaks = waiver(),
x.lab = NULL, y.lab = NULL, y.lim = NULL, y.breaks = waiver(),
shape.gen = 21, shape.env = 23, size.shape = 2.2, size.bor.tick = 0.3,
size.tex.lab = 12, size.tex.pa = 3.5, repulsion = 1, size.line = 0.5,
size.segm.line = 0.5, leg.lab = c("Env", "Gen"), line.type = "solid",
line.alpha = 0.9, col.line = "black", col.gen = "orange",
col.env = "forestgreen", col.alpha = 0.9, col.segm.gen = "transparent",
col.segm.env = "forestgreen", resolution = 300, ...) {

if (polygon == TRUE & type != 1) {
stop("The polygon can be drawn with type 1 graphic only.")
}

size.tex.leg <- size.tex.pa/0.2917
class <- class(x)
nenv <- nrow(subset(x$model, type == "ENV"))
ngen <- nrow(subset(x$model, type == "GEN"))

if (type == 1) {
y.lab = ifelse((is.null(y.lab) == F), y.lab, ifelse(class(x) ==
"waas", paste0("PC2 (", round(x$PCA[2, 6], 2), "%)"),
paste0("PC2 (", round(x$PCA[2, 3], 2), "%)")))
x.lab = ifelse((is.null(x.lab) == F), x.lab, ifelse(class(x) ==
"waas", paste0("PC1 (", round(x$PCA[1, 6], 2), "%)"),
paste0("PC1 (", round(x$PCA[1, 3], 2), "%)")))
if (is.null(x.lim) == F) {
x.lim <- x.lim
} else {
x.lim <- c(min(x$model$PC1 * axis.expand), max(x$model$PC1 *
axis.expand))
}

if (is.null(y.lim) == F) {
y.lim <- y.lim
} else {
y.lim <- c(min(x$model$PC2 * axis.expand), max(x$model$PC2 *
axis.expand))
}

p1 <- ggplot(x$model, aes(PC1, PC2, shape = type, fill = type)) +
geom_vline(xintercept = 0, linetype = line.type,
color = col.line, size = size.line, alpha = line.alpha) +
geom_hline(yintercept = 0, linetype = line.type,
color = col.line, size = size.line, alpha = line.alpha) +
geom_segment(data = x$model, aes(x = 0, y = 0, xend = PC1,
yend = PC2, size = type, color = type, group = type)) +
geom_point(size = size.shape, stroke = size.bor.tick,
aes(fill = type), alpha = col.alpha) + scale_shape_manual(labels = leg.lab,
values = c(shape.env, shape.gen)) + scale_fill_manual(labels = leg.lab,
values = c(col.env, col.gen)) + ggrepel::geom_text_repel(aes(PC1,
PC2, label = (Code)), size = size.tex.pa, col = c(rep(col.gen,
ngen), rep(col.env, nenv)), force = repulsion) +
theme %+replace% theme(aspect.ratio = 1, axis.text = element_text(size = size.tex.lab,
colour = "black"), axis.title = element_text(size = size.tex.lab,
colour = "black"), legend.text = element_text(size = size.tex.leg),
plot.title = element_text(size = size.tex.lab,
  hjust = 0, vjust = 1)) + labs(x = paste(x.lab),
y = paste(y.lab)) + scale_x_continuous(limits = x.lim,
breaks = x.breaks) + scale_y_continuous(limits = y.lim,
breaks = y.breaks) + scale_color_manual(name = "",
values = c(col.segm.env, col.segm.gen), theme(legend.position = "none")) +
scale_size_manual(name = "", values = c(size.segm.line,
size.segm.line), theme(legend.position = "none"))

if (polygon == TRUE) {

gen <- data.frame(subset(x$model, type == "GEN"))
coordgenotype <- data.frame(subset(x$model, type ==
"GEN"))[, 4:5]
coordenviroment <- data.frame(subset(x$model, type ==
"ENV"))[, 4:5]

hull <- chull(gen[, 4:5])
indice <- c(hull, hull[1])
segs <- NULL
limx <- x.lim
limy <- y.lim
i <- 1
while (is.na(indice[i + 1]) == FALSE) {
m <- (coordgenotype[indice[i], 2] - coordgenotype[indice[i +
  1], 2])/(coordgenotype[indice[i], 1] - coordgenotype[indice[i +
  1], 1])
mperp <- -1/m
c2 <- coordgenotype[indice[i + 1], 2] - m * coordgenotype[indice[i +
  1], 1]
xint <- -c2/(m - mperp)
xint <- ifelse(xint < 0, min(coordenviroment[,
  1], coordgenotype[, 1]), max(coordenviroment[,
  1], coordgenotype[, 1]))
yint <- mperp * xint
xprop <- ifelse(xint < 0, xint/limx[1], xint/limx[2])
yprop <- ifelse(yint < 0, yint/limy[1], yint/limy[2])
m3 <- which(c(xprop, yprop) == max(c(xprop, yprop)))
m2 <- abs(c(xint, yint)[m3])
if (m3 == 1 & xint < 0)
  sl1 <- (c(xint, yint)/m2) * abs(limx[1])
if (m3 == 1 & xint > 0)
  sl1 <- (c(xint, yint)/m2) * limx[2]
if (m3 == 2 & yint < 0)
  sl1 <- (c(xint, yint)/m2) * abs(limy[1])
if (m3 == 2 & yint > 0)
  sl1 <- (c(xint, yint)/m2) * limy[2]
segs <- rbind(segs, sl1)
i <- i + 1
}
rownames(segs) <- NULL
colnames(segs) <- NULL
segs <- data.frame(segs)

p1 <- p1 + geom_segment(aes(x = X1, y = X2), xend = 0,
yend = 0, linetype = 2, size = size.segm.line,
color = col.gen, data = segs, inherit.aes = FALSE) +
geom_polygon(data = gen[indice, ], fill = NA,
  col = col.gen, linetype = 2)
}

if (export == F | FALSE) {
return(p1)
} else if (file.type == "pdf") {
if (is.null(file.name)) {
pdf("PC1 x PC2.pdf", width = width, height = height)
} else pdf(paste0(file.name, ".pdf"), width = width,
height = height)
plot(p1)
dev.off()
}

if (file.type == "tiff") {
if (is.null(file.name)) {
tiff(filename = "PC1 x PC2.tiff", width = width,
  height = height, units = "in", compression = "lzw",
  res = resolution)
} else tiff(filename = paste0(file.name, ".tiff"),
width = width, height = height, units = "in",
compression = "lzw", res = resolution)
plot(p1)
dev.off()
}

}

if (type == 2) {
y.lab = ifelse(is.null(y.lab) == F, y.lab, ifelse(class(x) ==
"waas", paste0("PC1 (", round(x$PCA[1, 6], 2), "%)"),
paste0("PC1 (", round(x$PCA[1, 3], 2), "%)")))
x.lab = ifelse(is.null(x.lab) == F, x.lab, paste0("Grain yield"))
if (is.null(x.lim) == F) {
x.lim <- x.lim
} else {
x.lim <- c(min(x$model$Y) - (min(x$model$Y) * axis.expand -
min(x$model$Y)), max(x$model$Y) + (max(x$model$Y) *
axis.expand - max(x$model$Y)))
}

if (is.null(y.lim) == F) {
y.lim <- y.lim
} else {
y.lim <- c(min(x$model$PC1 * axis.expand), max(x$model$PC1 *
axis.expand))
}
mean <- mean(x$model$Y)
p2 <- ggplot2::ggplot(x$model, aes(Y, PC1, shape = type,
fill = type)) + geom_vline(xintercept = mean(x$model$Y),
linetype = line.type, color = col.line, size = size.line,
alpha = line.alpha) + geom_hline(yintercept = 0,
linetype = line.type, size = size.line, color = col.line,
alpha = line.alpha) + geom_segment(data = x$model,
aes(x = mean, y = 0, xend = Y, yend = PC1, size = type,
color = type, group = type)) + geom_point(size = size.shape,
stroke = size.bor.tick, aes(fill = type), alpha = col.alpha) +
scale_shape_manual(labels = leg.lab, values = c(shape.env,
shape.gen)) + scale_fill_manual(labels = leg.lab,
values = c(col.env, col.gen)) + ggrepel::geom_text_repel(aes(Y,
PC1, label = (Code)), size = size.tex.pa, col = c(rep(col.gen,
ngen), rep(col.env, nenv)), force = repulsion) +
theme %+replace% theme(aspect.ratio = 1, axis.text = element_text(size = size.tex.lab,
colour = "black"), axis.title = element_text(size = size.tex.lab,
colour = "black"), legend.text = element_text(size = size.tex.leg),
plot.title = element_text(size = size.tex.lab,
  hjust = 0, vjust = 1)) + labs(x = paste(x.lab),
y = paste(y.lab)) + scale_x_continuous(limits = x.lim,
breaks = x.breaks) + scale_y_continuous(limits = y.lim,
breaks = y.breaks) + scale_color_manual(name = "",
values = c(col.segm.env, col.segm.gen), theme(legend.position = "none")) +
scale_size_manual(name = "", values = c(size.segm.line,
size.segm.line), theme(legend.position = "none"))

if (export == F | FALSE) {
return(p2)
} else if (file.type == "pdf") {
if (is.null(file.name)) {
pdf("GY x PC1.pdf", width = width, height = height)
} else pdf(paste0(file.name, ".pdf"), width = width,
height = height)
plot(p2)
dev.off()

}

if (file.type == "tiff") {
if (is.null(file.name)) {
tiff(filename = "GY x PC1.tiff", width = width,
  height = height, units = "in", compression = "lzw",
  res = resolution)
} else tiff(filename = paste0(file.name, ".tiff"),
width = width, height = height, units = "in",
compression = "lzw", res = resolution)
plot(p2)
dev.off()
}
}

if (type == 3) {
y.lab = ifelse(is.null(y.lab) == F, y.lab, paste0("Weighted average of the absolute scores"))
x.lab = ifelse(is.null(x.lab) == F, x.lab, paste0("Grain yield"))
if (class == "waasb") {

if (is.null(x.lim) == F) {
x.lim <- x.lim
} else {
x.lim <- c(min(x$model$Y) - (min(x$model$Y) *
  axis.expand - min(x$model$Y)), max(x$model$Y) +
  (max(x$model$Y) * axis.expand - max(x$model$Y)))
}
if (is.null(y.lim) == F) {
y.lim <- y.lim
} else {
y.lim <- c(min(x$model$WAASB) - (min(x$model$WAASB) *
  axis.expand - min(x$model$WAASB)), max(x$model$WAASB) +
  (max(x$model$WAASB) * axis.expand - max(x$model$WAASB)))
}
m1 <- mean(x$model$Y)
m2 <- mean(x$model$WAASB)
I <- grid::grobTree(grid::textGrob("I", x = 0.02,
y = 0.98, hjust = 0))
II <- grid::grobTree(grid::textGrob("II", x = 0.97,
y = 0.97, hjust = 0))
III <- grid::grobTree(grid::textGrob("III", x = 0.01,
y = 0.03, hjust = 0))
IV <- grid::grobTree(grid::textGrob("IV", x = 0.96,
y = 0.03, hjust = 0))
p3 <- ggplot2::ggplot(x$model, aes(Y, WAASB, shape = type,
fill = type)) + geom_vline(xintercept = m1, linetype = line.type,
color = col.line, size = size.line, alpha = line.alpha) +
geom_hline(yintercept = m2, linetype = line.type,
  color = col.line, size = size.line, alpha = line.alpha) +
geom_point(size = size.shape, stroke = size.bor.tick,
  aes(fill = type), alpha = col.alpha) + scale_shape_manual(labels = leg.lab,
values = c(shape.env, shape.gen)) + scale_fill_manual(labels = leg.lab,
values = c(col.env, col.gen)) + ggrepel::geom_text_repel(aes(Y,
WAASB, label = (Code)), size = size.tex.pa, col = c(rep(col.gen,
ngen), rep(col.env, nenv)), force = repulsion) +
theme %+replace% theme(aspect.ratio = 1, axis.text = element_text(size = size.tex.lab,
  colour = "black"), axis.title = element_text(size = size.tex.lab,
  colour = "black"), legend.text = element_text(size = size.tex.leg),
  plot.title = element_text(size = size.tex.lab,
hjust = 0, vjust = 1)) + labs(x = paste(x.lab),
y = paste(y.lab)) + scale_x_continuous(limits = x.lim,
breaks = x.breaks) + scale_y_continuous(limits = y.lim,
breaks = y.breaks) + annotation_custom(I) + annotation_custom(II) +
annotation_custom(III) + annotation_custom(IV)
}

if (class == "waas") {
if (is.null(x.lim) == F) {
x.lim <- x.lim
} else {
x.lim <- c(min(x$model$Y) - (min(x$model$Y) *
  axis.expand - min(x$model$Y)), max(x$model$Y) +
  (max(x$model$Y) * axis.expand - max(x$model$Y)))
}

if (is.null(y.lim) == F) {
y.lim <- y.lim
} else {
y.lim <- c(min(x$model$WAAS) - (min(x$model$WAAS) *
  axis.expand - min(x$model$WAAS)), max(x$model$WAAS) +
  (max(x$model$WAAS) * axis.expand - max(x$model$WAAS)))
}
m1 <- mean(x$model$Y)
m2 <- mean(x$model$WAAS)
I <- grid::grobTree(grid::textGrob("I", x = 0.02,
y = 0.98, hjust = 0))
II <- grid::grobTree(grid::textGrob("II", x = 0.97,
y = 0.97, hjust = 0))
III <- grid::grobTree(grid::textGrob("III", x = 0.01,
y = 0.03, hjust = 0))
IV <- grid::grobTree(grid::textGrob("IV", x = 0.96,
y = 0.03, hjust = 0))
p3 <- ggplot2::ggplot(x$model, aes(Y, WAAS, shape = type,
fill = type)) + geom_vline(xintercept = m1, linetype = line.type,
color = col.line, size = size.line, alpha = line.alpha) +
geom_hline(yintercept = m2, linetype = line.type,
  color = col.line, size = size.line, alpha = line.alpha) +
geom_point(size = size.shape, stroke = size.bor.tick,
  aes(fill = type), alpha = col.alpha) + scale_shape_manual(labels = leg.lab,
values = c(shape.env, shape.gen)) + scale_fill_manual(labels = leg.lab,
values = c(col.env, col.gen)) + ggrepel::geom_text_repel(aes(Y,
WAAS, label = (Code)), size = size.tex.pa, col = c(rep(col.gen,
ngen), rep(col.env, nenv)), force = repulsion) +
theme %+replace% theme(aspect.ratio = 1, axis.text = element_text(size = size.tex.lab,
  colour = "black"), axis.title = element_text(size = size.tex.lab,
  colour = "black"), legend.text = element_text(size = size.tex.leg),
  plot.title = element_text(size = size.tex.lab,
hjust = 0, vjust = 1)) + labs(x = paste(x.lab),
y = paste(y.lab)) + scale_x_continuous(limits = x.lim,
breaks = x.breaks) + scale_y_continuous(limits = y.lim,
breaks = y.breaks) + annotation_custom(I) + annotation_custom(II) +
annotation_custom(III) + annotation_custom(IV)
}

if (export == F | FALSE) {
return(p3)
} else if (file.type == "pdf") {
if (is.null(file.name)) {
pdf("GY x WAASB.pdf", width = width, height = height)
} else pdf(paste0(file.name, ".pdf"), width = width,
height = height)
plot(p3)
dev.off()
}

if (file.type == "tiff") {
if (is.null(file.name)) {
tiff(filename = "GY x WAAS.tiff", width = width,
  height = height, units = "in", compression = "lzw",
  res = resolution)
} else tiff(filename = paste0(file.name, ".tiff"),
width = width, height = height, units = "in",
compression = "lzw", res = resolution)
plot(p3)
dev.off()
}
}

if (type == 4) {
data = as.data.frame(x[["MeansGxE"]])
minim <- min(data$nominal)
y.lab = ifelse(is.null(y.lab) == F, y.lab, paste0("Nominal Yield (Mg/ha)"))
x.lab = ifelse(is.null(x.lab) == F, x.lab, paste0("Environment PC1 [square root of  (Mg/ha)]"))

if (is.null(x.lim) == F) {
x.lim <- x.lim
} else {
x.lim <- c(min(data$envPC1), max(data$envPC1))
}

if (is.null(y.lim) == F) {
y.lim <- y.lim
} else {
y.lim <- c(min(data$nominal), max(data$nominal))

}
p4 <- ggplot2::ggplot(data, aes(x = envPC1, y = nominal,
group = GEN)) + geom_line(size = 1, aes(colour = GEN),
data = subset(data, envPC1 %in% c(max(envPC1), min(envPC1)))) +
geom_point(aes(x = envPC1, y = minim), data = subset(data,
GEN == data[1, 2])) + ggrepel::geom_label_repel(data = subset(data,
envPC1 == min(envPC1)), aes(label = GEN, fill = GEN),
size = size.tex.pa, color = "white", force = repulsion,
segment.color = "#bbbbbb") + ggrepel::geom_text_repel(aes(x = envPC1,
y = minim, label = ENV), size = size.tex.pa, force = repulsion,
data = subset(data, GEN == data[1, 2])) + theme %+replace%
theme(legend.position = "none", axis.text = element_text(size = size.tex.lab,
colour = "black"), axis.title = element_text(size = size.tex.lab,
colour = "black"), plot.title = element_text(size = size.tex.lab,
hjust = 0, vjust = 1)) + scale_x_continuous(limits = x.lim,
breaks = x.breaks) + scale_y_continuous(limits = y.lim,
breaks = y.breaks) + labs(x = paste(x.lab), y = y.lab)

if (export == F | FALSE) {
return(p4)
} else if (file.type == "pdf") {
if (is.null(file.name)) {
pdf("Adaptative reponses (AMMI1 model).pdf",
  width = width, height = height)
} else pdf(paste0(file.name, ".pdf"), width = width,
height = height)
plot(p4)
dev.off()

}

if (file.type == "tiff") {
if (is.null(file.name)) {
tiff(filename = "Adaptative reponses (AMMI1 model).tiff",
  width = width, height = height, units = "in",
  compression = "lzw", res = resolution)
} else tiff(filename = paste0(file.name, ".tiff"),
width = width, height = height, units = "in",
compression = "lzw", res = resolution)
plot(p4)
dev.off()
}
}
}

