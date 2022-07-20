#' Plot scores in different graphical interpretations
#'
#' @description
#' `r badge('stable')`
#'
#' Plot scores of genotypes and environments in different graphical
#' interpretations.
#'
#' Biplots type 1 and 2 are well known in AMMI analysis. In the plot type 3, the
#' scores of both genotypes and environments are plotted considering the
#' response variable and the WAASB, an stability index that considers all
#' significant principal component axis of traditional AMMI models or all
#' principal component axis estimated with BLUP-interaction effects (Olivoto et
#' al. 2019). Plot type 4 may be used to better understand the well known
#' 'which-won-where' pattern, facilitating the recommendation of appropriate
#' genotypes targeted for specific environments, thus allowing the exploitation
#' of narrow adaptations.
#'
#' @param x An object fitted with the functions [performs_ammi()],
#'   [waas()], [waas_means()], or [waasb()].
#' @param var The variable to plot. Defaults to `var = 1` the first
#'   variable of `x`.
#' @param type type of biplot to produce
#' * `type = 1` The default. Produces an AMMI1 biplot (Y x PC1) to make
#' inferences related to stability and productivity.
#' * `type = 2` Produces an AMMI2 biplot (PC1 x PC2) to make inferences
#' related to the interaction effects. Use the arguments `first` or
#' `second` to change the default IPCA shown in the plot.
#' * `type = 3` Valid for objects of class `waas` or `waasb`,
#' produces a biplot showing the GY x WAASB.
#' * `type = 4` Produces a plot with the Nominal yield x Environment PC.
#' @param first,second The IPCA to be shown in the first (x) and second (y)
#'   axis. By default, IPCA1 is shown in the `x` axis and IPCA2 in the
#'   `y` axis. For example, use `second = "PC3"` to shown the IPCA3 in
#'   the `y` axis.
#' @param repel If `TRUE` (default), the text labels repel away from each
#'   other and away from the data points.
#' @param repulsion Force of repulsion between overlapping text labels. Defaults
#'   to `1`.
#' @param max_overlaps Exclude text labels that overlap too many things. Defaults to 20.
#' @param polygon Logical argument. If `TRUE`, a polygon is drawn when
#'   `type = 2`.
#' @param title Logical values (Defaults to `TRUE`) to include
#'   automatically generated titles
#' @param plot_theme The graphical theme of the plot. Default is
#'   `plot_theme = theme_metan()`. For more details, see
#'   [ggplot2::theme()].
#' @param axis.expand Multiplication factor to expand the axis limits by to
#'   enable fitting of labels. Default is `1.1`.
#' @param x.lim,y.lim The range of x and y axes, respectively. Default is
#'   `NULL` (maximum and minimum values of the data set). New values can be
#'   inserted as `x.lim = c(x.min, x.max)` or `y.lim = c(y.min,
#'   y.max)`.
#' @param x.breaks,y.breaks The breaks to be plotted in the x and y axes,
#'   respectively. Defaults to `waiver()` (automatic breaks). New values
#'   can be inserted, for example, as `x.breaks = c(0.1, 0.2, 0.3)` or
#'   `x.breaks = seq(0, 1, by = 0.2)`
#' @param x.lab,y.lab The label of x and y axes, respectively. Defaults to
#'   `NULL`, i.e., each plot has a default axis label. New values can be
#'   inserted as `x.lab = 'my label'`.
#' @param shape.gen,shape.env The shape for genotypes and environments
#'   indication in the biplot. Default is `21` (circle) for genotypes and
#'   `23` (diamond) for environments. Values must be between `21-25`:
#'   `21` (circle), `22` (square), `23` (diamond), `24` (up
#'   triangle), and `25` (low triangle).
#' @param size.shape.gen,size.shape.env The size of the shapes for genotypes and
#'   environments respectively. Defaults to `2.2`.
#' @param size.bor.tick The size of tick of shape. Default is `0.1`. The
#'   size of the shape will be `max(size.shape.gen, size.shape.env) +
#'   size.bor.tick`
#' @param size.tex.lab,size.tex.gen,size.tex.env The size of the text for axis
#'   labels (Defaults to 12), genotypes labels, and environments labels
#'   (Defaults to 3.5).
#' @param size.line The size of the line that indicate the means in the biplot.
#'   Default is `0.5`.
#' @param size.segm.line The size of the segment that start in the origin of the
#'   biplot and end in the scores values. Default is `0.5`.
#' @param col.bor.gen,col.bor.env The color of the shape's border for genotypes
#'   and environments, respectively.
#' @param col.line The color of the line that indicate the means in the biplot.
#'   Default is `'gray'`
#' @param col.gen,col.env The shape color for genotypes (Defaults to
#'   `'blue'`) and environments (`'forestgreen'`). Must be length
#'   one or a vector of colors with the same length of the number of
#'   genotypes/environments.
#' @param col.alpha.gen,col.alpha.env The alpha value for the color for
#'   genotypes and environments, respectively. Defaults to `NA`. Values must be
#'   between `0` (full transparency) to `1` (full color).
#' @param col.segm.gen,col.segm.env The color of segment for genotypes (Defaults
#'   to `transparent_color()`) and environments (Defaults to 'forestgreen'),
#'   respectively. Valid arguments for plots with `type = 1` or `type
#'   = 2` graphics.
#'
#' @param highlight Genotypes/environments to be highlight in the plot. Defaults
#'   to `NULL`.
#' @param col.highlight The color for shape/labels when a value is provided in
#'   `highlight.` Defaults to `"red"`.
#' @param col.alpha.highlight The alpha value for the color of the highlighted
#'   genotypes. Defaults to `1`.
#' @param size.tex.highlight The size of the text for the highlighted genotypes.
#'   Defaults to `5.5`.
#' @param size.shape.highlight The size of the shape for the highlighted
#'   genotypes. Defaults to `3.2`.
#' @param leg.lab The labs of legend. Default is `Gen` and `Env`.
#' @param line.type The type of the line that indicate the means in the biplot.
#'   Default is `'solid'`. Other values that can be attributed are:
#'   `'blank'`, no lines in the biplot, `'dashed', 'dotted',
#'   'dotdash', 'longdash', and 'twodash'`.
#' @param line.alpha The alpha value that combine the line with the background
#'   to create the appearance of partial or full transparency. Default is
#'   `0.4`. Values must be between '0' (full transparency) to '1' (full
#'   color).
#' @param resolution deprecated
#' @param export Export (or not) the plot. Default is `FALSE`. If `TRUE`, calls the
#'   [ggplot2::ggsave()] function.
#' @param file.name The name of the file for exportation, default is
#'   `NULL`, i.e. the files are automatically named.
#' @param file.type The type of file to be exported. Currently recognises the extensions
#'   eps/ps, tex, pdf, jpeg, tiff, png (default), bmp, svg and wmf (windows only).
#' @param width The width 'inch' of the plot. Default is `8`.
#' @param height The height 'inch' of the plot. Default is `7`.
#' @param color Should type 4 plot have colors? Default to `TRUE`.
#' @param ... Currently not used.
#' @md
#' @references
#' Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, V.S. Marchioro, V.Q. de
#' Souza, and E. Jost. 2019. Mean performance and stability in multi-environment
#' trials I: Combining features of AMMI and BLUP techniques. Agron. J.
#' 111:2949-2960. \doi{10.2134/agronj2019.03.0220}
#'
#' @return An object of class `gg, ggplot`.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso [plot_eigen()]
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' # AMMI model
#'model <- waas(data_ge,
#'              env = ENV,
#'              gen = GEN,
#'              rep = REP,
#'              resp = everything())
#'
#' # GY x PC1 for variable GY (default plot)
#' plot_scores(model)
#'
#' # PC1 x PC2 (variable HM)
#' #
#' plot_scores(model,
#'             polygon = TRUE,            # Draw a convex hull polygon
#'             var = "HM",                # or var = 2 to select variable
#'             highlight = c("G1", "G2"), # Highlight genotypes 2 and 3
#'             type = 2)                  # type of biplot
#'
#' # PC3 x PC4 (variable HM)
#' # Change size of plot fonts and colors
#' # Minimal theme
#'plot_scores(model,
#'            var = "HM",
#'            type = 2,
#'            first = "PC3",
#'            second = "PC4",
#'            col.gen = "black",
#'            col.env = "gray",
#'            col.segm.env = "gray",
#'            size.tex.gen = 5,
#'            size.tex.env = 2,
#'            size.tex.lab = 16,
#'            plot_theme = theme_metan_minimal())
#'
#' # WAASB index
#' waasb_model <- waasb(data_ge, ENV, GEN, REP, GY)
#'
#' # GY x WAASB
#' # Highlight genotypes 2 and 8
#' plot_scores(waasb_model,
#'             type = 3,
#'             highlight = c("G2", "G8"))
#' }
plot_scores <- function(x,
                        var = 1,
                        type = 1,
                        first = "PC1",
                        second = "PC2",
                        repel = TRUE,
                        repulsion = 1,
                        max_overlaps = 20,
                        polygon = FALSE,
                        title = TRUE,
                        plot_theme = theme_metan(),
                        axis.expand = 1.1,
                        x.lim = NULL,
                        y.lim = NULL,
                        x.breaks = waiver(),
                        y.breaks = waiver(),
                        x.lab = NULL,
                        y.lab = NULL,
                        shape.gen = 21,
                        shape.env = 23,
                        size.shape.gen = 2.2,
                        size.shape.env = 2.2,
                        size.bor.tick = 0.1,
                        size.tex.lab = 12,
                        size.tex.gen = 3.5,
                        size.tex.env = 3.5,
                        size.line = 0.5,
                        size.segm.line = 0.5,
                        col.bor.gen = "black",
                        col.bor.env = "black",
                        col.line = "black",
                        col.gen = "blue",
                        col.env = "forestgreen",
                        col.alpha.gen = 1,
                        col.alpha.env = 1,
                        col.segm.gen = transparent_color(),
                        col.segm.env = "forestgreen",
                        highlight = NULL,
                        col.highlight = "red",
                        col.alpha.highlight = 1,
                        size.tex.highlight = 5.5,
                        size.shape.highlight = 3.2,
                        leg.lab = c("Env", "Gen"),
                        line.type = "solid",
                        line.alpha = 0.9,
                        resolution = deprecated(),
                        file.type = "png",
                        export = FALSE,
                        file.name = NULL,
                        width = 8,
                        height = 7,
                        color = TRUE,
                        ...) {
  if(is_present(resolution)){
    warning("`resolution` is deprecated as of metan 1.17.")
  }
  varname <- names(x)[var]
  x <- x[[var]]
  df <-
    x$model %>%
    mutate(Color = ifelse(type == "ENV", col.env, col.gen))
  test <- !missing(highlight)
  if(!all(highlight %in% df$Code)){
    not_in_code <- highlight[which(!highlight %in% df$Code)]
    stop(paste(not_in_code, collapse = ", "), " not present in the labels. Please, check and fix it.", call. = FALSE)
  }
  df %<>%
    mutate(type2 = ifelse(Code %in% highlight, "Selected", type),
           Color = case_when(type2 == "Selected" ~ col.highlight,
                             type2 == "GEN" ~ col.gen,
                             type2 == "ENV" ~ col.env),
           size.text = case_when(type2 == "Selected" ~ size.tex.highlight,
                                 type2 == "GEN" ~ size.tex.gen,
                                 type2 == "ENV" ~ size.tex.env),
           size.shape = case_when(type2 == "Selected" ~ size.shape.highlight,
                                  type2 == "GEN" ~ size.shape.gen,
                                  type2 == "ENV" ~ size.shape.env),
           alpha.col = case_when(type2 == "Selected" ~ col.alpha.highlight,
                                 type2 == "GEN" ~ col.alpha.gen,
                                 type2 == "ENV" ~ col.alpha.env),
           alpha.col.line = ifelse(type2 == "ENV", col.alpha.env, 0))
  if (polygon == TRUE & type != 2) {
    stop("The polygon can be drawn with type 2 graphic only.", call. = FALSE)
  }
  if (inherits(x, "performs_ammi") & type == 3) {
    stop("Biplot type invalid. Type 3 biplot can only be made with objects of class 'waas' or 'waasb'.", call. = FALSE)
  }
  size.tex.leg <- max(size.tex.env, size.tex.gen)/0.2917
  class <- class(x)
  nenv <- nrow(subset(df, type == "ENV"))
  ngen <- nrow(subset(df, type == "GEN"))
  color.bor <- c(rep(col.bor.gen, ngen), rep(col.bor.env, nenv))

  if (type == 1) {
    df <- df
    if(!is.null(y.lab)){
      y.lab <- y.lab
    } else{
      if(class %in% c("waas", "performs_ammi")){
        y.lab <- paste0("PC1 (", round(x$PCA[1, 7], 2), "%)")
      }
      if(class == "waasb"){
        y.lab <- paste0("PC1 (", round(x$PCA[1, 3], 2), "%)")
      }
      if(class == "waas_means"){
        y.lab <- paste0("PC1 (", round(x$proportion[1], 2), "%)")
      }
    }
    x.lab = ifelse(is.null(x.lab) == F, x.lab, paste0(varname))
    if (!is.null(x.lim)) {
      x.lim <- x.lim
    } else {
      x.lim <- c(min(df$Y) - (min(df$Y) * axis.expand - min(df$Y)),
                 max(df$Y) + (max(df$Y) * axis.expand - max(df$Y)))
    }

    if (is.null(y.lim) == FALSE) {
      y.lim <- y.lim
    } else {
      y.lim <- c(min(df$PC1 * axis.expand), max(df$PC1 * axis.expand))
    }
    mean <- mean(df$Y)

    p1 <- ggplot2::ggplot(df, aes(Y, PC1)) +
      geom_vline(xintercept = mean(df$Y),
                 linetype = line.type,
                 color = col.line,
                 size = size.line,
                 alpha = line.alpha) +
      geom_hline(yintercept = 0,
                 linetype = line.type,
                 size = size.line,
                 color = col.line,
                 alpha = line.alpha) +
      geom_segment(data = df,
                   show.legend = FALSE,
                   aes(x = mean,
                       y = 0,
                       xend = Y,
                       yend = PC1,
                       size = type,
                       color = type,
                       group = type),
                   alpha = df$alpha.col.line) +
      {if(!test)geom_point(aes(shape = type, fill = type),
                           alpha = df$alpha.col,
                           size = df$size.shape,
                           stroke = size.bor.tick,
                           color = color.bor)} +
      {if(test)geom_point(aes(shape = type),
                          alpha = df$alpha.col,
                          size = df$size.shape,
                          fill = df$Color,
                          stroke = size.bor.tick,
                          color = color.bor)}+
      plot_theme %+replace%
      theme(aspect.ratio = 1,
            axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            legend.text = element_text(size = size.tex.leg)) +
      labs(x = x.lab, y = y.lab) +
      scale_x_continuous(limits = x.lim, breaks = x.breaks) +
      scale_y_continuous(limits = y.lim, breaks = y.breaks) +
      scale_color_manual(name = "", values = c(col.segm.env, col.segm.gen), theme(legend.position = "none")) +
      scale_size_manual(name = "", values = c(size.segm.line, size.segm.line), theme(legend.position = "none"))+
      scale_shape_manual(labels = leg.lab, values = c(shape.env, shape.gen)) +
      {if(!test)scale_fill_manual(labels = leg.lab, values = c(col.env, col.gen))} +
      {if(test)scale_fill_identity()} +
      {if(repel)geom_text_repel(aes(Y, PC1, label = Code),
                                color = df$Color,
                                force = repulsion,
                                size =  df$size.text,
                                alpha = df$alpha.col,
                                max.overlaps = max_overlaps)} +
      {if(!repel)geom_text(aes(Y, PC1, label = Code),
                           alpha = df$alpha.col,
                           color = df$Color,
                           size = df$size.text,
                           hjust = "outward",
                           vjust = "outward")} +
      {if(title)ggtitle("AMMI1 Biplot")}
    if (export == FALSE) {
      return(p1)
    } else{
      if(is.null(file.name)) {
        ggsave("GY x PC1.png", p1, width = width, height = height)
      } else{
        ggsave(paste0(file.name, ".", file.type), p1, width = width, height = height)
      }
    }
  }

  if (type == 2) {
    first <- tidy_strings(first, sep = "")
    second <- tidy_strings(second, sep = "")
    if(class %in% c("waas", "performs_ammi", "waasb")){
      PCA <- x$PCA %>% select_cols(PC, Proportion, Accumulated)
      if(!extract_string(first) == "PC"){
        stop("Argument 'first' invalid. Please, use 'PC1', 'PC2', ..., 'PCn'.", call. = FALSE)
      }
      if(!extract_string(second) == "PC"){
        stop("Argument 'second' invalid. Please, use 'PC1', 'PC2', ..., 'PCn'.", call. = FALSE)
      }
      if(extract_number(first) > nrow(PCA)){
        stop("The number of principal components in 'first' is greater than those in model (", nrow(PCA), ").", call. = FALSE)
      }
      if(extract_number(second) > nrow(PCA)){
        stop("The number of principal components in 'second' is greater than those in model (", nrow(PCA), ").", call. = FALSE)
      }
      PCA <- subset(PCA, PC %in% c(first, second))
      if(!is.null(y.lab)){
        y.lab <- y.lab
      } else{
        if(class %in% c("waas", "performs_ammi", "waasb")){
          y.lab <- paste0(second, " (", round(PCA[which(PCA[, 1] == second), 2], 2), "%)")
        }
      }
      if(!is.null(x.lab)){
        x.lab <- x.lab
      } else{
        if(class %in% c("waas", "performs_ammi", "waasb")){
          x.lab <- paste0(first, " (", round(PCA[which(PCA[, 1] == first), 2], 2), "%)")
        }
      }
    } else{
      y.lab <- paste0(second, " (", round(x$proportion[as.numeric(substr(second, 3, nchar(second)))], 2), "%)")
      x.lab <- paste0(first, " (", round(x$proportion[as.numeric(substr(first, 3, nchar(first)))], 2), "%)")
    }
    df <- df %>% select_cols(type, Code, all_of(first), all_of(second), Color, size.text,
                             size.shape, alpha.col, alpha.col.line)
    if (!is.null(x.lim)) {
      x.lim <- x.lim
    } else {
      x.lim <- c(min(df[[first]] * axis.expand), max(df[[first]] * axis.expand))
    }
    if (!is.null(y.lim)) {
      y.lim <- y.lim
    } else {
      y.lim <- c(min(df[[second]] * axis.expand), max(df[[second]] * axis.expand))
    }
    p2 <-
      ggplot(df, aes(!!sym(first), !!sym(second), shape = type, fill = type)) +
      geom_vline(xintercept = 0,
                 linetype = line.type,
                 color = col.line,
                 size = size.line,
                 alpha = line.alpha) +
      geom_hline(yintercept = 0,
                 linetype = line.type,
                 color = col.line, size = size.line,
                 alpha = line.alpha) +
      geom_segment(data = df,
                   show.legend = FALSE,
                   aes(x = 0,
                       y = 0,
                       xend = !!sym(first),
                       yend = !!sym(second),
                       size = type,
                       color = type,
                       group = type),
                   alpha = df$alpha.col.line) +
      {if(!test)geom_point(aes(shape = type, fill = type),
                           alpha = df$alpha.col,
                           size = df$size.shape,
                           stroke = size.bor.tick,
                           color = color.bor)} +
      {if(test)geom_point(aes(shape = type),
                          alpha = df$alpha.col,
                          size = df$size.shape,
                          fill = df$Color,
                          stroke = size.bor.tick,
                          color = color.bor)}+
      plot_theme %+replace%
      theme(aspect.ratio = 1,
            axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            legend.text = element_text(size = size.tex.leg)) +
      labs(x = x.lab, y = y.lab) +
      scale_x_continuous(limits = x.lim, breaks = x.breaks) +
      scale_y_continuous(limits = y.lim, breaks = y.breaks) +
      scale_color_manual(name = "", values = c(col.segm.env, col.segm.gen), theme(legend.position = "none")) +
      scale_size_manual(name = "", values = c(size.segm.line, size.segm.line), theme(legend.position = "none"))+
      scale_shape_manual(labels = leg.lab, values = c(shape.env, shape.gen)) +
      {if(!test)scale_fill_manual(labels = leg.lab, values = c(col.env, col.gen))} +
      {if(test)scale_fill_identity()} +
      {if(repel)geom_text_repel(aes(!!sym(first), !!sym(second), label = Code),
                                force = repulsion,
                                color = df$Color,
                                size = df$size.text,
                                alpha = df$alpha.col,
                                max.overlaps = max_overlaps)} +
      {if(!repel)geom_text(aes(!!sym(first), !!sym(second), label = Code),
                           color = df$Color,
                           size =  df$size.text,
                           alpha = df$alpha.col,
                           hjust = "outward",
                           vjust = "outward")} +
      {if(title)ggtitle("AMMI2 Biplot")}
    if (polygon == TRUE) {
      gen <- data.frame(subset(df, type == "GEN"))
      coordgenotype <- data.frame(subset(df, type == "GEN"))[, 3:4]
      coordenviroment <- data.frame(subset(df, type == "ENV"))[, 3:4]
      hull <- chull(gen[, 3:4])
      indice <- c(hull, hull[1])
      segs <- NULL
      limx <- x.lim
      limy <- y.lim
      i <- 1
      while (is.na(indice[i + 1]) == FALSE) {
        m <- (coordgenotype[indice[i], 2] - coordgenotype[indice[i + 1], 2])/
          (coordgenotype[indice[i], 1] - coordgenotype[indice[i + 1], 1])
        mperp <- -1/m
        c2 <- coordgenotype[indice[i + 1], 2] - m * coordgenotype[indice[i + 1], 1]
        xint <- -c2/(m - mperp)
        xint <- ifelse(xint < 0,
                       min(coordenviroment[, 1], coordgenotype[, 1]),
                       max(coordenviroment[, 1], coordgenotype[, 1]))
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
      p2 <-
        p2 +
        geom_segment(aes(x = X1, y = X2),
                     xend = 0,
                     yend = 0,
                     linetype = 2,
                     size = size.segm.line,
                     color = col.gen,
                     data = segs,
                     show.legend = FALSE,
                     inherit.aes = FALSE,
                     alpha = line.alpha) +
        geom_polygon(data = gen[indice, ],
                     fill = NA,
                     col = col.gen,
                     linetype = 2)
    }
    if (export == FALSE) {
      return(p2)
    } else{
      if(is.null(file.name)) {
        ggsave("PC1 x PC2.png", p2, width = width, height = height)
      } else{
        ggsave(paste0(file.name, ".", file.type), p2, width = width, height = height)
      }
    }
  }


  if (type == 3) {
    y.lab = ifelse(!is.null(y.lab), y.lab, paste0("Weighted average of the absolute scores"))
    x.lab = ifelse(!is.null(x.lab), x.lab, paste0(varname))
    if (class == "waasb") {
      if (!is.null(x.lim)) {
        x.lim <- x.lim
      } else {
        x.lim <- c(min(df$Y) - (min(df$Y) * axis.expand - min(df$Y)),
                   max(df$Y) + (max(df$Y) * axis.expand - max(df$Y)))
      }
      if (!is.null(y.lim)) {
        y.lim <- y.lim
      } else {
        y.lim <- c(min(df$WAASB) - (min(df$WAASB) * axis.expand - min(df$WAASB)),
                   max(df$WAASB) + (max(df$WAASB) * axis.expand - max(df$WAASB)))
      }
      m1 <- mean(df$Y)
      m2 <- mean(df$WAASB)

      p3 <- ggplot(df, aes(Y, WAASB, shape = type, fill = type))+
        geom_vline(xintercept = m1,
                   linetype = line.type,
                   color = col.line,
                   size = size.line,
                   alpha = line.alpha) +
        geom_hline(yintercept = m2,
                   linetype = line.type,
                   color = col.line,
                   size = size.line,
                   alpha = line.alpha)+
        {if(!test)geom_point(aes(shape = type, fill = type),
                             alpha = df$alpha.col,
                             size = df$size.shape,
                             stroke = size.bor.tick,
                             color = color.bor)} +
        {if(test)geom_point(aes(shape = type),
                            alpha = df$alpha.col,
                            size = df$size.shape,
                            fill = df$Color,
                            stroke = size.bor.tick,
                            color = color.bor)}+

        {if(repel)geom_text_repel(aes(Y, WAASB, label = (Code)),
                                  color = df$Color,
                                  force = repulsion,
                                  size =  df$size.text,
                                  alpha = df$alpha.col,
                                  max.overlaps = max_overlaps)} +
        {if(!repel)geom_text(aes(Y, WAASB, label = (Code)),
                             color = df$Color,
                             size =  df$size.text,
                             alpha = df$alpha.col,
                             hjust = "outward",
                             vjust = "outward")} +
        {if(!test)scale_fill_manual(labels = leg.lab, values = c(col.env, col.gen))} +
        {if(test)scale_fill_identity()}
    }
    if (class %in% c("waas", "waas_means")) {
      if (!is.null(x.lim)) {
        x.lim <- x.lim
      } else {
        x.lim <- c(min(df$Y) - (min(df$Y) * axis.expand - min(df$Y)),
                   max(df$Y) + (max(df$Y) * axis.expand - max(df$Y)))
      }
      if (is.null(y.lim) == FALSE) {
        y.lim <- y.lim
      } else {
        y.lim <- c(min(df$WAAS) - (min(df$WAAS) * axis.expand - min(df$WAAS)),
                   max(df$WAAS) + (max(df$WAAS) * axis.expand - max(df$WAAS)))
      }
      m1 <- mean(df$Y)
      m2 <- mean(df$WAAS)
      p3 <-
        ggplot(df, aes(Y, WAAS, shape = type, fill = type))+
        geom_vline(xintercept = m1,
                   linetype = line.type,
                   color = col.line,
                   size = size.line,
                   alpha = line.alpha) +
        geom_hline(yintercept = m2,
                   linetype = line.type,
                   color = col.line,
                   size = size.line,
                   alpha = line.alpha)+
        {if(!test)geom_point(aes(shape = type, fill = type),
                             alpha = df$alpha.col,
                             size = df$size.shape,
                             stroke = size.bor.tick,
                             color = color.bor)} +
        {if(test)geom_point(aes(shape = type),
                            alpha = df$alpha.col,
                            size = df$size.shape,
                            fill = df$Color,
                            stroke = size.bor.tick,
                            color = color.bor)}+
        {if(repel)geom_text_repel(aes(Y, WAAS, label = (Code)),
                                  force = repulsion,
                                  color = df$Color,
                                  size =  df$size.text,
                                  alpha = df$alpha.col,
                                  max.overlaps = max_overlaps)} +
        {if(!repel)geom_text(aes(Y, WAAS, label = Code),
                             color = df$Color,
                             size =  df$size.text,
                             alpha = df$alpha.col,
                             hjust = "outward",
                             vjust = "outward")} +
        {if(!test)scale_fill_manual(labels = leg.lab, values = c(col.env, col.gen))} +
        {if(test)scale_fill_identity()}

    }
    p3 <-
      p3 +
      scale_shape_manual(labels = leg.lab, values = c(shape.env, shape.gen)) +
      plot_theme %+replace%
      theme(aspect.ratio = 1,
            axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            legend.text = element_text(size = size.tex.leg)) +
      labs(x = x.lab, y = y.lab) +
      scale_x_continuous(limits = x.lim, breaks = x.breaks) +
      scale_y_continuous(limits = y.lim, breaks = y.breaks) +
      annotate("text",
               x = x.lim[1],
               y = y.lim[2],
               hjust = 3,
               vjust = 0,
               label = "I") +
      annotate("text",
               x = x.lim[2],
               y = y.lim[2],
               hjust = -1,
               vjust = 0,
               label = "II") +
      annotate("text",
               x = x.lim[1],
               y = y.lim[1],
               hjust = 1,
               vjust = 1,
               label = "III") +
      annotate("text",
               x = x.lim[2],
               y = y.lim[1],
               hjust = 0,
               vjust = 1,
               label = "IV")
    if(title == TRUE){
      p3 <- p3 + ggtitle(ifelse(class == "waasb", "Y x WAASB biplot", "Y x WAAS biplot"))
    }
    if (export == FALSE) {
      return(p3)
    } else{
      if(is.null(file.name)) {
        ggsave("GY x WAASB.png", p3, width = width, height = height)
      } else{
        ggsave(paste0(file.name, ".", file.type), p3, width = width, height = height)
      }
    }
  }


  if (type == 4) {
    if(class == "waas_means"){
      EscENV <- subset(df, type ==  "ENV") %>%
        select_cols(Code, Y, PC1) %>%
        rename(ENV = Code)
      EscGEN <- subset(df, type ==  "GEN") %>%
        select_cols(Code, Y, PC1) %>%
        rename(GEN = Code)
      data <- x$ge_means
      data <- suppressMessages(
        suppressWarnings(
          mutate(data,
                 envPC1 = left_join(data, EscENV %>% select(ENV, PC1))$PC1,
                 genPC1 = left_join(data, EscGEN %>% select(GEN, PC1))$PC1,
                 nominal = left_join(data, EscGEN %>% select(GEN, Y))$Y + genPC1 * envPC1)
        )
      )%>%
        as.data.frame()
    } else{
      data <- as.data.frame(x[["MeansGxE"]])
    }
    minim <- min(data$nominal)
    y.lab = ifelse(is.null(y.lab) == F, y.lab, paste0("Nominal Yield (Mg/ha)"))
    x.lab = ifelse(is.null(x.lab) == F, x.lab, paste0("Environment PC1 [square root of  (Mg/ha)]"))

    if (!is.null(x.lim)) {
      x.lim <- x.lim
    } else {
      x.lim <- c(min(data$envPC1), max(data$envPC1))
    }

    if (is.null(y.lim) == FALSE) {
      y.lim <- y.lim
    } else {
      y.lim <- c(min(data$nominal), max(data$nominal))
    }
    p4 <- ggplot2::ggplot(data, aes(x = envPC1, y = nominal, group = GEN)) +
      {if(color)geom_line(data = subset(data, envPC1 %in% c(max(envPC1), min(envPC1))),
                          aes(colour = GEN),
                          size = 0.8)} +
      {if(!color)geom_line(data = subset(data, envPC1 %in% c(max(envPC1), min(envPC1))),
                           size = 0.8)} +
      {if(repel)geom_label_repel(data = subset(data, envPC1 == min(envPC1)),
                                 aes(label = GEN, color = GEN),
                                 size = size.tex.gen,
                                 force = repulsion,
                                 alpha = rep(col.alpha.gen, ngen))} +
      {if(repel)geom_text_repel(data = subset(data, GEN == data[1, 2]),
                                aes(x = envPC1, y = minim, label = ENV),
                                size = size.tex.env,
                                force = repulsion,
                                alpha = rep(col.alpha.env, nenv),
                                max.overlaps = max_overlaps)} +
      {if(!repel)geom_label(data = subset(data, envPC1 == min(envPC1)),
                            aes(label = GEN, color = GEN),
                            size = size.tex.gen,
                            alpha = rep(col.alpha.gen, ngen),
                            hjust = "center",
                            vjust = "top")} +
      {if(!repel)geom_text(data = subset(data, GEN == data[1, 2]),
                           aes(x = envPC1, y = minim, label = ENV),
                           size = size.tex.env,
                           alpha = rep(col.alpha.env, nenv),
                           vjust = -1)} +
      geom_point(data = subset(data, GEN == data[1, 2]),
                 aes(x = envPC1, y = minim),
                 shape = 17,
                 size = 2.8) +
      plot_theme %+replace%
      theme(legend.position = "none",
            axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black")) +
      scale_x_continuous(limits = x.lim, breaks = x.breaks) +
      scale_y_continuous(limits = y.lim, breaks = y.breaks,
                         expand = expansion(mult = c(0.003, 0.1))) +
      labs(x = x.lab, y = y.lab)
    if(title == TRUE){
      p4 <- p4 + ggtitle("Nominal yield plot")
    }
    if (export == FALSE) {
      return(p4)
    } else{
      if(is.null(file.name)) {
        ggsave("Adaptative reponses (AMMI1 model).png", p4, width = width, height = height)
      } else{
        ggsave(paste0(file.name, ".", file.type), p4, width = width, height = height)
      }
    }
  }
}

