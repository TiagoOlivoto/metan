#' Genotype plus genotype-by-environment model
#'
#' Produces genotype plus genotype-by-environment model based on a multi-environment
#' trial dataset containing at least the columns for genotypes, environments and one
#' response variable or a two-way table.
#'
#'
#' @param .data The dataset containing the columns related to Environments, Genotypes
#' and the response variable(s).
#' @param env The name of the column that contains the levels of the environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure a vector of variables may be used. For example \code{resp
#'   = c(var1, var2, var3)}. Select helpers are also supported.
#' @param centering The centering method. Must be one of the \code{'none | 0'}, for no
#'  centering; \code{'global | 1'}, for global centered (E+G+GE); \code{'environment | 2'} (default),
#'  for environment-centered (G+GE); or \code{'double | 3'}, for double centred (GE).
#'   A biplot cannot be produced with models produced without centering.
#' @param scaling The scaling method. Must be one of the \code{'none | 0'} (default), for no scaling;
#'  or \code{'sd | 1'}, where each value is divided by the standard deviation of its corresponding
#'   environment (column). This will put all environments roughly he same rang of values.
#'
#' @param svp The method for singular value partitioning. Must be one of the \code{'genotype | 1'},
#'  (The singular value is entirely partitioned into the genotype eigenvectors, also called row
#'  metric preserving); \code{'environment | 2'}, default, (The singular value is entirely partitioned into the
#'  environment eigenvectors, also called column metric preserving); or \code{'symmetrical | 3'}
#'  (The singular value is symmetrically partitioned into the genotype and the environment eigenvectors
#'  This SVP is most often used in AMMI analysis and other biplot analysis, but it is not ideal for
#'  visualizing either the relationship among genotypes or that among the environments).
#'
#' @param ... Arguments passed to the function
#'   \code{\link{impute_missing_val}()} for imputation of missing values in case
#'   of unbalanced data.
#'
#' @return The function returns a list of class \code{gge} containing the following objects
#'
#'  * \strong{coordgen} The coordinates for genotypes for all components.
#'
#'  * \strong{coordenv} The coordinates for environments for all components.
#'
#'  * \strong{eigenvalues} The vector of eigenvalues.
#'
#'  * \strong{totalvar} The overall variance.
#'
#'  * \strong{labelgen} The name of the genotypes.
#'
#'  * \strong{labelenv} The names of the environments.
#'
#'  * \strong{labelaxes} The axes labels.
#'
#'  * \strong{ge_mat} The data used to produce the model (scaled and centered).
#'
#'  * \strong{centering} The centering method.
#'
#'  * \strong{scaling} The scaling method.
#'
#'  * \strong{svp} The singular value partitioning method.
#'
#'  * \strong{d} The factor used to generate in which the ranges of genotypes and environments
#'  are comparable when singular value partitioning is set to 'genotype' or 'environment'.
#'  * \strong{grand_mean} The grand mean of the trial.
#'  * \strong{mean_gen} A vector with the means of the genotypes.
#'  * \strong{mean_env} A vector with the means of the environments.
#'  * \strong{scale_var} The scaling vector when the scaling method is \code{'sd'}.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Yan, W., and M.S. Kang. 2003. GGE biplot analysis: a graphical tool for breeders,
#'  geneticists, and agronomists. CRC Press.
#'
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' mod <- gge(data_ge, ENV, GEN, GY)
#' plot(mod)
#'
#' # GGE model for all numeric variables
#' mod2 <- gge(data_ge2, ENV, GEN, resp = everything())
#' plot(mod2, var = "ED")
#'
#' # If we have a two-way table with the mean values for
#' # genotypes and environments
#'
#' table <- make_mat(data_ge, GEN, ENV, GY) %>% round(2)
#' table
#' make_long(table) %>%
#' gge(ENV, GEN, Y) %>%
#' plot()
#'}
gge <- function(.data,
                env,
                gen,
                resp,
                centering = "environment",
                scaling = "none",
                svp = "environment",
                ...) {
  factors  <-
    .data %>%
    select({{env}}, {{gen}}) %>%
    mutate(across(everything(), as.factor))
  vars <- .data %>% select({{resp}}, -names(factors))
  vars %<>% select_numeric_cols()
  factors %<>% set_names("ENV", "GEN")
  listres <- list()
  nvar <- ncol(vars)
  for (var in 1:nvar) {
    ge_mat <-
      factors %>%
      mutate(Y = vars[[var]]) %>%
      make_mat(GEN, ENV, Y) %>%
      as.matrix()
    if(has_na(ge_mat)){
      ge_mat <- impute_missing_val(ge_mat, verbose = FALSE, ...)$.data
      warning("Data imputation used to fill the GxE matrix", call. = FALSE)
    }
    grand_mean <- mean(ge_mat)
    mean_env <- colMeans(ge_mat)
    mean_gen <- rowMeans(ge_mat)
    scale_val <- apply(ge_mat, 2, sd)
    labelgen <- rownames(ge_mat)
    labelenv <- colnames(ge_mat)
    if (any(is.na(ge_mat))) {
      stop("missing data in input data frame")
    }
    if (any(apply(ge_mat, 2, is.numeric) == FALSE)) {
      stop("not all columns are of class 'numeric'")
    }
    if (!(centering %in% c("none", "environment", "global", "double") |
          centering %in% 0:3)) {
      warning(paste("Centering method", centering, "not found; defaulting to environment centered"))
      centering <- "environment"
    }
    if (!(svp %in% c("genotype", "environment", "symmetrical") |
          svp %in% 1:3)) {
      warning(paste("svp method", svp, "not found; defaulting to column metric preserving"))
      svp <- "environment"
    }
    if (!(scaling %in% c("none", "sd") | scaling %in% 0:1)) {
      warning(paste("scaling method", scaling, "not found; defaulting to no scaling"))
      sd <- "none"
    }
    labelaxes <- paste("PC", 1:ncol(diag(svd(ge_mat)$d)), sep = "")
    # Centering
    if (centering == 1 | centering == "global") {
      ge_mat <- ge_mat - mean(ge_mat)
    }
    if (centering == 2 | centering == "environment") {
      ge_mat <- sweep(ge_mat, 2, colMeans(ge_mat))
    }
    if (centering == 3 | centering == "double") {
      grand_mean <- mean(ge_mat)
      mean_env <- colMeans(ge_mat)
      mean_gen <- rowMeans(ge_mat)
      for (i in 1:nrow(ge_mat)) {
        for (j in 1:ncol(ge_mat)) {
          ge_mat[i, j] <- ge_mat[i, j] + grand_mean - mean_env[j] -
            mean_gen[i]
        }
      }
    }
    # Scaling
    if (scaling == 1 | scaling == "sd") {
      ge_mat <- sweep(ge_mat, 2, apply(ge_mat, 2, sd), FUN = "/")
    }
    # Singular value partitioning
    if (svp == 1 | svp == "genotype") {
      coordgen <- svd(ge_mat)$u %*% diag(svd(ge_mat)$d)
      coordenv <- svd(ge_mat)$v
      d1 <- (max(coordenv[, 1]) - min(coordenv[, 1]))/(max(coordgen[,
                                                                    1]) - min(coordgen[, 1]))
      d2 <- (max(coordenv[, 2]) - min(coordenv[, 2]))/(max(coordgen[,
                                                                    2]) - min(coordgen[, 2]))
      coordenv <- coordenv/max(d1, d2)
    }
    if (svp == 2 | svp == "environment") {
      coordgen <- svd(ge_mat)$u
      coordenv <- svd(ge_mat)$v %*% diag(svd(ge_mat)$d)
      d1 <- (max(coordgen[, 1]) - min(coordgen[, 1]))/(max(coordenv[,
                                                                    1]) - min(coordenv[, 1]))
      d2 <- (max(coordgen[, 2]) - min(coordgen[, 2]))/(max(coordenv[,
                                                                    2]) - min(coordenv[, 2]))
      coordgen <- coordgen/max(d1, d2)
    }
    if (svp == 3 | svp == "symmetrical") {
      coordgen <- svd(ge_mat)$u %*% diag(sqrt(svd(ge_mat)$d))
      coordenv <- svd(ge_mat)$v %*% diag(sqrt(svd(ge_mat)$d))
    }
    eigenvalues <- svd(ge_mat)$d
    totalvar <- round(as.numeric(sum(eigenvalues^2)), 2)
    varexpl <- round(as.numeric((eigenvalues^2/totalvar) * 100),
                     2)
    if (svp == "genotype" | svp == "environment") {
      d <- max(d1, d2)
    } else {
      d <- NULL
    }
    tmp <- structure(
      list(coordgen = coordgen, coordenv = coordenv, eigenvalues = eigenvalues,
           totalvar = totalvar, varexpl = varexpl, labelgen = labelgen,
           labelenv = labelenv, labelaxes = labelaxes, ge_mat = ge_mat,
           centering = centering, scaling = scaling, svp = svp,
           d = d, grand_mean = grand_mean, mean_gen = mean_gen,
           mean_env = mean_env, scale_val = scale_val),
      class = "gge")
    if (nvar > 1) {
      listres[[paste(names(vars[var]))]] <-  tmp
    } else {
      listres[[paste(names(vars[var]))]] <- tmp
    }
  }
  return(structure(listres, class = "gge"))
}








#' Create GGE biplots
#'
#' Produces a ggplot2-based GGE biplot based on a model of class \code{gge}.
#' Since the output is an object of class \code{ggplot}, all stylistic attributes
#' of the output can be customized using the power of plot customization provided by
#' ggplot2.
#'
#'@param x An object of class \code{gge}
#'@param var The variable to plot. Defaults to \code{var = 1} the first variable
#'  of \code{x}.
#'@param type The type of biplot to produce.
#' \enumerate{
#' \item Basic biplot.
#' \item Mean performance vs. stability.
#' \item Which-won-where.
#' \item Discriminativeness vs. representativeness.
#' \item Examine an environment.
#' \item Ranking environments.
#' \item Examine a genotype.
#' \item Ranking genotypes.
#' \item Compare two genotypes.
#' \item Relationship among environments}
#' @param sel_env,sel_gen The name of the environment and genotype to examine
#'   when \code{type = 5} and   \code{type = 7}, respectively. Must be a string
#'   which matches a environment or genotype label.
#' @param sel_gen1,sel_gen2 The name of genotypes to compare between when
#'   \code{type = 9}. Must be a string present in the genotype's name.
#' @param shape.gen,shape.env The shape for genotype and environment indication
#'   in the biplot. Defaults to \code{shape.gen = 21} (circle) for genotypes and
#'   \code{shape.env = 23} (rhombus) for environments. Values must be between
#'   \code{21-25}: \code{21} (circle), \code{22} (square), \code{23} (rhombus),
#'   \code{24} (up triangle), and \code{25} (low triangle).
#' @param size.shape The size of the shape (both for genotypes and
#'   environments). Defaults to \code{2.2}.
#' @param size.shape.win The size of the shape for winners genotypes when
#'   \code{type = 3}. Defaults to \code{3.2}.
#' @param size.bor.tick The size of tick of shape. Default is \code{0.3}. The
#'   size of the shape will be \code{size.shape + size.bor.tick}
#' @param col.gen,col.env Color for genotype and environment attributes in the
#'   biplot. Defaults to \code{col.gen = 'blue'} and \code{col.env =
#'   'forestgreen'}
#' @param col.alpha The alpha value for the color. Defaults to \code{1}. Values
#'   must be between \code{0} (full transparency) to \code{1} (full color).
#' @param col.circle,col.alpha.circle The color and alpha values for the circle
#'   lines. Defaults to \code{'gray'} and \code{0.4}, respectively.
#' @param leg.lab The labs of legend. Default is \code{c('Env', 'Gen')}.
#' @param size.text.gen,size.text.env,size.text.lab The size of the text for
#'   genotypes, environments and labels, respectively.
#' @param size.line The size of the line in biplots (Both for segments and circles).
#' @param large_label The text size to use for larger labels where \code{type =
#'   3}, used for the outermost genotypes and where \code{type = 9}, used for
#'   the two selected genotypes. Defaults to 4.5
#' @param axis_expand multiplication factor to expand the axis limits by to
#'   enable fitting of labels. Defaults to 1.2
#' @param title Logical values (Defaults to \code{TRUE}) to include
#'   automatically generated information in the plot such as singular value
#'   partitioning, scaling and centering.
#' @param plot_theme The graphical theme of the plot. Default is
#'   \code{plot_theme = theme_metan()}. For more details, see
#'   \code{\link[ggplot2]{theme}}.
#' @param ... Currently not used.
#' @return A ggplot2-based biplot.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Yan, W., and M.S. Kang. 2003. GGE biplot analysis: a graphical
#'   tool for breeders, geneticists, and agronomists. CRC Press.
#' @method plot gge
#' @importFrom ggforce geom_arc
#' @export
#' @return An object of class \code{gg, ggplot}.
#' @examples
#' \donttest{
#' library(metan)
#' mod <- gge(data_ge, ENV, GEN, GY)
#' plot(mod)
#' plot(mod,
#'      type = 2,
#'      col.gen = 'blue',
#'      col.env = 'red',
#'      size.text.gen = 2)
#' }
plot.gge <- function(x,
                     var = 1,
                     type = 1,
                     sel_env = NA,
                     sel_gen = NA,
                     sel_gen1 = NA,
                     sel_gen2 = NA,
                     shape.gen = 21,
                     shape.env = 23,
                     size.shape = 2.2,
                     size.shape.win = 3.2,
                     size.bor.tick = 0.3,
                     col.gen = "blue",
                     col.env = "forestgreen",
                     col.alpha = 1,
                     col.circle = "gray",
                     col.alpha.circle = 0.5,
                     leg.lab = c("Env", "Gen"),
                     size.text.gen = 4,
                     size.text.env = 4,
                     size.text.lab = 12,
                     size.line = 0.8,
                     large_label = 4.5,
                     axis_expand = 1.2,
                     title = TRUE,
                     plot_theme = theme_metan(),
                     ...) {
  if(any(class(x) == "gtb")){
    if(all(leg.lab %in% c("Gen", "Env")))
      leg.lab <- c("Gen", "Trait")
  }
  model <- x[[var]]
  if (!class(model) == "gge") {
    stop("The model must be of class 'gge'")
  }
  coord_gen <- model$coordgen[, c(1, 2)]
  coord_env <- model$coordenv[, c(1, 2)]
  varexpl <- model$varexpl[c(1, 2)]
  labelgen <- model$labelgen
  labelenv <- model$labelenv
  labelaxes <- model$labelaxes[c(1, 2)]
  Data <- model$ge_mat
  centering <- model$centering
  svp <- model$svp
  scaling <- model$scaling
  ngen <- nrow(coord_gen)
  nenv <- nrow(coord_env)
  if (centering == "none") {
    stop("It is not possible to create a GGE biplot with a model produced without centering")
  }
  plotdata <- data.frame(rbind(data.frame(coord_gen,
                                          type = "genotype",
                                          label = labelgen),
                               data.frame(coord_env,
                                          type = "environment",
                                          label = labelenv)))
  colnames(plotdata)[1:2] <- c("d1", "d2")
  xlim <- c(min(plotdata$d1 * axis_expand),
            max(plotdata$d1 * axis_expand))
  ylim <- c(min(plotdata$d2 * axis_expand),
            max(plotdata$d2 * axis_expand))
  if (which(c(diff(xlim), diff(ylim)) == max(c(diff(xlim), diff(ylim)))) == 1) {
    xlim1 <- xlim
    ylim1 <- c(ylim[1] - (diff(xlim) - diff(ylim))/2,
               ylim[2] + (diff(xlim) - diff(ylim))/2)
  }
  if (which(c(diff(xlim), diff(ylim)) == max(c(diff(xlim), diff(ylim)))) == 2) {
    ylim1 <- ylim
    xlim1 <- c(xlim[1] - (diff(ylim) - diff(xlim))/2,
               xlim[2] + (diff(ylim) - diff(xlim))/2)
  }
  # Base plot
  P1 <-
    ggplot(data = plotdata, aes(x = d1, y = d2, group = "type")) +
    scale_color_manual(values = c(col.env, col.gen)) +
    scale_size_manual(values = c(size.text.gen, size.text.env)) +
    xlab(paste(labelaxes[1], " (", round(varexpl[1], 2), "%)", sep = "")) +
    ylab(paste(labelaxes[2]," (", round(varexpl[2], 2), "%)", sep = "")) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    scale_x_continuous(limits = xlim1, expand = c(0, 0)) +
    scale_y_continuous(limits = ylim1, expand = c(0, 0)) +
    coord_fixed() +
    plot_theme %+replace%
    theme(axis.text = element_text(size = size.text.lab, colour = "black"),
          axis.title = element_text(size = size.text.lab, colour = "black"))
  # Basic plot
  if (type == 1) {
    P2 <- P1 +
      geom_segment(xend = 0,
                   yend = 0,
                   size = size.line,
                   col = alpha(col.env, col.alpha),
                   data = subset(plotdata, type == "environment")) +
      geom_point(aes(d1, d2, fill = type, shape = type),
                 size = size.shape,
                 stroke = size.bor.tick,
                 alpha = col.alpha) +
      scale_shape_manual(labels = leg.lab,
                         values = c(shape.env, shape.gen)) +
      scale_fill_manual(labels = leg.lab,
                        values = c(col.env, col.gen)) +
      geom_text_repel(aes(col = type,
                          label = label,
                          size = type),
                      show.legend = FALSE,
                      col = c(rep(col.gen, ngen), rep(col.env, nenv))) +
      plot_theme %+replace%
      theme(axis.text = element_text(size = size.text.lab, colour = "black"),
            axis.title = element_text(size = size.text.lab, colour = "black"))
    if (title == TRUE) {
      if(any(class(x) == "gtb")){
        ggt <- ggtitle("GT Biplot")
      } else{
        ggt <- ggtitle("GGE Biplot")
      }
    }
  }
  # Mean vs. stability
  if (type == 2) {
    med1 <- mean(coord_env[, 1])
    med2 <- mean(coord_env[, 2])
    x1 <- NULL
    for (i in 1:nrow(Data)) {
      x <- solve(matrix(c(-med2, med1, med1, med2), nrow = 2),
                 matrix(c(0, med2 * coord_gen[i, 2] + med1 * coord_gen[i, 1]), ncol = 1))
      x1 <- rbind(x1, t(x))
    }
    plotdata$x1_x <- NA
    plotdata$x1_x[plotdata$type == "genotype"] <- x1[, 1]
    plotdata$x1_y <- NA
    plotdata$x1_y[plotdata$type == "genotype"] <- x1[, 2]
    P2 <- P1 +
      geom_segment(aes(xend = x1_x, yend = x1_y),
                   color = col.gen,
                   linetype = 2,
                   size = size.line,
                   data = subset(plotdata, type == "genotype")) +
      geom_abline(intercept = 0,
                  slope = med2/med1,
                  color = col.env,
                  size = size.line) +
      geom_abline(intercept = 0,
                  slope = -med1/med2,
                  color = col.env,
                  size = size.line) +
      geom_segment(x = 0,
                   y = 0,
                   xend = med1,
                   yend = med2,
                   arrow = arrow(length = unit(0.3, "cm")),
                   size = size.line,
                   color = col.env) +
      geom_text(aes(color = type,
                    label = label,
                    size = type),
                show.legend = FALSE)
    if (title == TRUE) {
      ggt <- ggtitle("Mean vs. Stability")
    }
  }
  # Which-won-where
  if (type == 3) {
    indice <- c(grDevices::chull(coord_gen[, 1], coord_gen[, 2]))
    polign <- data.frame(coord_gen[indice, ])
    indice <- c(indice, indice[1])
    segs <- NULL
    limx <- layer_scales(P1)$x$limits
    limy <- layer_scales(P1)$y$limits
    i <- 1
    while (is.na(indice[i + 1]) == FALSE) {
      m <- (coord_gen[indice[i], 2] - coord_gen[indice[i + 1], 2])/(coord_gen[indice[i], 1] - coord_gen[indice[i + 1], 1])
      mperp <- -1/m
      c2 <- coord_gen[indice[i + 1], 2] - m * coord_gen[indice[i + 1], 1]
      xint <- -c2/(m - mperp)
      xint <- ifelse(xint < 0, min(coord_env[, 1], coord_gen[, 1]), max(coord_env[, 1], coord_gen[, 1]))
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
    colnames(polign) <- c("X1", "X2")
    winners <- plotdata[plotdata$type == "genotype", ][indice[-1], ] %>% add_cols(win = "yes")
    others <- anti_join(plotdata, winners, by = "label") %>% add_cols(win = "no")
    df_winners <- rbind(winners, others)
    P2 <- P1 +
      geom_polygon(data = polign,
                   aes(x = X1, y = X2),
                   fill = NA,
                   col = col.gen,
                   size = size.line) +
      geom_segment(data = segs,
                   aes(x = X1, y = X2),
                   xend = 0,
                   yend = 0,
                   linetype = 2,
                   color = col.gen,
                   size = size.line) +
      geom_point(data = subset(df_winners, win == "no"),
                 aes(d1, d2, fill = type, shape = type),
                 size = size.shape,
                 stroke = size.bor.tick,
                 alpha = col.alpha,
                 show.legend = FALSE) +
      geom_text_repel(data = subset(df_winners, win == "no"),
                      aes(col = type, label = label, size = type),
                      show.legend = FALSE) +
      geom_point(data = subset(df_winners, win == "yes"),
                 aes(d1, d2, fill = type, shape = type),
                 size = size.shape.win,
                 stroke = size.bor.tick,
                 alpha = col.alpha) +
      geom_text_repel(data = subset(df_winners, win == "yes"),
                      aes(col = type, label = label),
                      show.legend = FALSE,
                      fontface = "bold",
                      size = large_label)
    if("genotype" %in% others$type){
      P2 <-
        P2 +
        scale_shape_manual(labels = leg.lab, values = c(shape.env, shape.gen)) +
        scale_fill_manual(labels = leg.lab, values = c(col.env, col.gen))
    } else{
      P2 <-
        suppressMessages(
          P2 +
            scale_shape_manual(labels = leg.lab[c(2,1)], values = c(shape.env, shape.gen)) +
            scale_fill_manual(labels = leg.lab[c(2,1)], values = c(col.env, col.gen)) +
            scale_color_manual(labels = leg.lab[c(2,1)], values = c(col.env, col.gen)) +
            scale_size_manual(labels = leg.lab[c(2,1)], values = c(size.text.env, size.text.gen))
        )
    }
    if (title == TRUE) {
      ggt <- ggtitle("Which-won-where pattern")
    }
  }
  # Discrimination vs. representativeness
  if (type == 4) {
    circles <- data.frame(x0 = 0,
                          y0 = 0,
                          start = 0,
                          end = pi * 2,
                          radio = 1:5 * max((max(coord_env[1, ]) - min(coord_env[1, ])),
                                            (max(coord_env[2, ]) - min(coord_env[2, ])))/10)
    P2 <- P1 +
      geom_arc(data = circles,
               aes(r = radio,
                   x0 = x0,
                   y0 = y0,
                   start = start,
                   end = end),
               color = col.circle,
               alpha = col.alpha.circle,
               size = size.line,
               inherit.aes = F,
               na.rm = TRUE) +
      geom_segment(xend = 0,
                   yend = 0,
                   x = mean(coord_env[, 1]),
                   y = mean(coord_env[, 2]),
                   arrow = arrow(ends = "first", length = unit(0.1, "inches")),
                   size = size.line,
                   color = col.env) +
      geom_abline(intercept = 0,
                  slope = mean(coord_env[, 2])/mean(coord_env[, 1]),
                  color = col.env,
                  size = size.line) +
      geom_point(aes(d1, d2, fill = type, shape = type),
                 size = size.shape,
                 stroke = size.bor.tick,
                 alpha = col.alpha) +
      scale_shape_manual(labels = leg.lab, values = c(shape.env, shape.gen)) +
      scale_fill_manual(labels = leg.lab,  values = c(col.env, col.gen)) +
      geom_segment(data = subset(plotdata, type == "environment"),
                   xend = 0,
                   yend = 0,
                   color = alpha(col.env, col.alpha),
                   size = size.line) +
      geom_text_repel(aes(col = type, label = label, size = type),
                      show.legend = FALSE,
                      color = c(rep(col.gen, ngen), rep(col.env, nenv)))
    if (title == TRUE) {
      ggt <- ggtitle("Discriminativeness vs. representativeness")
    }
  }
  # Examine an environment
  if (type == 5) {
    if (!sel_env %in% labelenv) {
      stop(paste("The environment", sel_env, "is not in the list of environment labels"))
    }
    venvironment <- labelenv == sel_env
    x1 <- NULL
    for (i in 1:nrow(Data)) {
      x <- solve(matrix(c(-coord_env[venvironment, 2],
                          coord_env[venvironment, 1],
                          coord_env[venvironment, 1],
                          coord_env[venvironment, 2]),
                        nrow = 2),
                 matrix(c(0, coord_env[venvironment, 1] * coord_gen[i, 1] +
                            coord_env[venvironment, 2] * coord_gen[i, 2]), ncol = 1))
      x1 <- rbind(x1, t(x))
    }
    plotdata$x1_x <- NA
    plotdata$x1_x[plotdata$type == "genotype"] <- x1[, 1]
    plotdata$x1_y <- NA
    plotdata$x1_y[plotdata$type == "genotype"] <- x1[, 2]
    P2 <- P1 + geom_segment(data = subset(plotdata, type == "genotype"),
                            aes(xend = x1_x, yend = x1_y),
                            col = col.gen,
                            linetype = "dotted") +
      geom_abline(slope = coord_env[venvironment, 2]/coord_env[venvironment, 1],
                  intercept = 0,
                  color = alpha(col.env, col.alpha),
                  size = size.line) +
      geom_abline(slope = -coord_env[venvironment, 1]/coord_env[venvironment, 2],
                  intercept = 0,
                  col = alpha(col.env, col.alpha),
                  size = size.line) +
      geom_segment(data = subset(plotdata, type == "environment" & label == sel_env),xend = 0,
                   yend = 0,
                   col = alpha(col.env, col.alpha),
                   size = size.line,
                   arrow = arrow(ends = "first", length = unit(0.5, "cm"))) +
      geom_text(data = subset(plotdata, type == "genotype"),
                aes(label = label),
                show.legend = FALSE,
                color = col.gen,
                size = size.text.gen) +
      geom_text(data = subset(plotdata, type == "environment" & label %in% sel_env),
                aes(label = label),
                show.legend = FALSE,
                col = col.env,
                size = size.text.env,
                fontface = "bold",
                hjust = "outward",
                vjust = "outward") +
      geom_point(data = subset(plotdata, type == "environment" & label == sel_env),
                 aes(d1, d2),
                 shape = 1,
                 size = 4)
    if (title == TRUE) {
      ggt <- ggtitle(paste("Environment:", sel_env))
    }
  }
  # Ranking environments
  if (type == 6) {
    med1 <- mean(coord_env[, 1])
    med2 <- mean(coord_env[, 2])
    mod <- max((coord_env[, 1]^2 + coord_env[, 2]^2)^0.5)
    xcoord <- sign(med1) * (mod^2/(1 + med2^2/med1^2))^0.5
    ycoord <- (med2/med1) * xcoord
    circles <- data.frame(x0 = xcoord,
                          y0 = ycoord,
                          start = 0,
                          end = pi * 2,
                          radio = 1:8 * ((xcoord - med1)^2 + (ycoord - med2)^2)^0.5/3)
    P2 <- P1 + geom_abline(intercept = 0,
                           slope = med2/med1,
                           col = col.env,
                           size = size.line) +
      geom_abline(intercept = 0,
                  slope = -med1/med2,
                  color = col.env,
                  size = size.line) +
      geom_arc(data = circles,
               aes(r = radio,
                   x0 = x0,
                   y0 = y0,
                   start = start,
                   end = end),
               col = col.circle,
               alpha = col.alpha.circle,
               inherit.aes = F,
               na.rm = TRUE) +
      geom_segment(x = 0,
                   y = 0,
                   xend = xcoord,
                   yend = ycoord,
                   arrow = arrow(length = unit(0.15, "inches")),
                   size = size.line,
                   color = col.env) +
      geom_point(fill = col.env,
                 shape = shape.env,
                 size = size.shape,
                 stroke = size.bor.tick,
                 alpha = col.alpha,
                 data = subset(plotdata, type == "environment")) +
      geom_text_repel(data = subset(plotdata, type == "environment"),
                      aes(label = label),
                      size = size.text.env,
                      show.legend = FALSE,
                      col = col.env) +
      geom_point(aes(xcoord, ycoord), shape = 1, size = 4)
    if (title == TRUE) {
      ggt <- ggtitle("Ranking Environments")
    }
  }
  # Examine a genotype
  if (type == 7) {
    if (!sel_gen %in% labelgen) {
      stop(paste("The genotype", sel_gen, "is not in the list of genotype labels"))
    }
    vgenotype <- labelgen == sel_gen
    x1 <- NULL
    for (i in 1:ncol(Data)) {
      x <- solve(matrix(c(-coord_gen[vgenotype, 2],
                          coord_gen[vgenotype, 1],
                          coord_gen[vgenotype, 1],
                          coord_gen[vgenotype, 2]),
                        nrow = 2),
                 matrix(c(0, coord_gen[vgenotype, 1] * coord_env[i, 1] +
                            coord_gen[vgenotype, 2] * coord_env[i, 2]),
                        ncol = 1))
      x1 <- rbind(x1, t(x))
    }
    plotdata$x1_x <- NA
    plotdata$x1_x[plotdata$type == "environment"] <- x1[, 1]
    plotdata$x1_y <- NA
    plotdata$x1_y[plotdata$type == "environment"] <- x1[, 2]
    P2 <- P1 + geom_segment(data = subset(plotdata, type == "environment"),
                            aes(xend = x1_x, yend = x1_y),
                            col = col.env,
                            linetype = "dotted") +
      geom_abline(slope = coord_gen[vgenotype, 2]/coord_gen[vgenotype, 1],
                  intercept = 0,
                  color = alpha(col.gen, col.alpha),
                  size = size.line) +
      geom_abline(slope = -coord_gen[vgenotype, 1]/coord_gen[vgenotype, 2],
                  intercept = 0,
                  col = alpha(col.gen, col.alpha),
                  size = size.line) +
      geom_segment(data = subset(plotdata, type == "genotype" & label == sel_gen),
                   xend = 0,
                   yend = 0,
                   col = alpha(col.gen, col.alpha),
                   arrow = arrow(ends = "first", length = unit(0.5, "cm")),
                   size = size.line) +
      geom_text(data = subset(plotdata, type == "environment"),
                aes(label = label),
                show.legend = FALSE,
                color = col.env,
                size = size.text.gen) +
      geom_text(data = subset(plotdata, type == "genotype" & label %in% sel_gen),
                aes(label = label),
                show.legend = FALSE,
                color = col.gen,
                size = size.text.gen,
                fontface = "bold",
                hjust = "outward",
                vjust = "outward") +
      geom_point(data = subset(plotdata, type == "genotype" & label == sel_gen),
                 aes(d1, d2),
                 shape = 1,
                 size = 4)
    if (title == TRUE) {
      ggt <- ggtitle(paste("Genotype:", sel_gen))
    }
  }
  # Ranking genotypes
  if (type == 8) {
    med1 <- mean(coord_env[, 1])
    med2 <- mean(coord_env[, 2])
    coordx <- 0
    coordy <- 0
    for (i in 1:nrow(Data)) {
      x <- solve(matrix(c(-med2, med1, med1, med2), nrow = 2),
                 matrix(c(0, med2 * coord_gen[i, 2] + med1 * coord_gen[i, 1]), ncol = 1))
      if (sign(x[1]) == sign(med1)) {
        if (abs(x[1]) > abs(coordx)) {
          coordx <- x[1]
          coordy <- x[2]
        }
      }
    }
    circles <- data.frame(x0 = coordx,
                          y0 = coordy,
                          radio = 1:10 * ((coordx - med1)^2 + (coordy - med2)^2)^0.5/3)
    P2 <- P1 +
      geom_abline(intercept = 0,
                  slope = med2/med1,
                  col = col.gen,
                  size = size.line) +
      geom_abline(intercept = 0,
                  slope = -med1/med2,
                  color = col.gen,
                  size = size.line) +
      geom_arc(data = circles,
               aes(r = radio,
                   x0 = x0,
                   y0 = y0,
                   start = 0,
                   end = pi * 2),
               color = col.circle,
               alpha = col.alpha.circle,
               size = size.line,
               inherit.aes = F) +
      geom_segment(x = 0,
                   y = 0,
                   xend = coordx,
                   yend = coordy,
                   arrow = arrow(length = unit(0.15, "inches")),
                   size = size.line,
                   color = col.gen) +
      geom_point(aes(d1, d2, fill = type, shape = type),
                 size = size.shape,
                 stroke = size.bor.tick,
                 alpha = col.alpha) +
      scale_shape_manual(labels = leg.lab, values = c(shape.env, shape.gen)) +
      scale_fill_manual(labels = leg.lab, values = c(col.env, col.gen)) +
      geom_text_repel(aes(col = type, label = label, size = type),
                      show.legend = FALSE,
                      color = c(rep(col.gen, ngen), rep(col.env, nenv))) +
      geom_point(aes(coordx, coordy), shape = 1, size = 3)
    if (title == TRUE) {
      ggt <- ggtitle("Ranking Genotypes")
    }
  }
  # Compare two genotypes
  if (type == 9) {
    if (!sel_gen1 %in% labelgen) {
      stop(paste("The genotype", sel_gen1, "is not in the list of genotype labels"))
    }
    if (!sel_gen2 %in% labelgen) {
      stop(paste("The genotype", sel_gen2, "is not in the list of genotype labels"))
    }
    if (sel_gen1 == sel_gen2) {
      stop(paste("It is not possible to compare a genotype to itself"))
    }
    vgenotype1 <- labelgen == sel_gen1
    vgenotype2 <- labelgen == sel_gen2
    P2 <- P1 +
      geom_segment(x = plotdata$d1[plotdata$label == sel_gen1 & plotdata$type == "genotype"],
                   xend = plotdata$d1[plotdata$label == sel_gen2 & plotdata$type == "genotype"],
                   y = plotdata$d2[plotdata$label == sel_gen1 & plotdata$type == "genotype"],
                   yend = plotdata$d2[plotdata$label == sel_gen2 & plotdata$type == "genotype"],
                   col = col.gen,
                   size = size.line) +
      geom_abline(intercept = 0,
                  slope = -(coord_gen[vgenotype1, 1] - coord_gen[vgenotype2, 1])/
                    (coord_gen[vgenotype1, 2] - coord_gen[vgenotype2, 2]),
                  color = col.gen,
                  size = size.line) +
      geom_text(aes(label = label),
                show.legend = FALSE,
                data = subset(plotdata, type == "environment"),
                col = col.env,
                size = size.text.env) +
      geom_text(data = subset(plotdata, type == "genotype" & !label %in% c(sel_gen1, sel_gen2)),
                aes(label = label),
                show.legend = FALSE,
                col = alpha(col.gen, alpha = col.alpha),
                size = size.text.gen) +
      geom_text(data = subset(plotdata, type == "genotype" & label %in% c(sel_gen1, sel_gen2)),
                aes(label = label),
                show.legend = FALSE,
                col = col.gen,
                size = large_label,
                hjust = "outward",
                vjust = "outward") +
      geom_point(data = subset(plotdata, type == "genotype" & label %in% c(sel_gen1, sel_gen2)),
                 shape = shape.gen,
                 fill = col.gen,
                 size = size.shape) +
      plot_theme
    if (title == TRUE) {
      ggt <- ggtitle(paste("Comparison of Genotype", sel_gen1, "with Genotype", sel_gen2))
    }
  }
  # Relationship among environments
  if (type == 10) {
    P2 <- P1 +
      geom_segment(data = subset(plotdata, type == "environment"),
                   xend = 0,
                   yend = 0,
                   col = alpha(col.env, col.alpha),
                   size = size.line) +
      geom_point(data = subset(plotdata, type == "environment"),
                 shape = shape.env,
                 fill = col.env,
                 size = size.shape,
                 stroke = size.bor.tick,
                 alpha = col.alpha) +
      geom_text_repel(data = subset(plotdata, type == "environment"),
                      aes(label = label),
                      show.legend = FALSE,
                      col = col.env,
                      size = size.text.env)
    if (title == TRUE) {
      if(any(class(x) == "gtb")){
        ggt <- ggtitle("Relationship Among Traits")
      } else{
        ggt <- ggtitle("Relationship Among Environments")
      }
    }
  }
  if (title == T) {
    scal_text <- ifelse(scaling == 1 | scaling == "sd", "Scaling = 1", "Scaling = 0")
    cent_text <-
      case_when(
        centering == 1 | centering == "global" ~ "Centering = 1",
        centering == 2 | centering == "environment" | centering == "trait" ~ "Centering = 2",
        centering == 3 | centering == "double" ~ "Centering = 3",
        FALSE ~ "No Centering"
      )
    svp_text <-
      case_when(
        svp == 1 | svp == "genotype" ~ "SVP = 1",
        svp == 2 | svp == "environment" | svp == "trait" ~ "SVP = 2",
        svp == 3 | svp == "symmetrical" ~ "SVP = 3"
      )
    annotationtxt <- paste(scal_text, ", ", cent_text, ", ", svp_text, sep = "")
    P2 <- P2 + ggtitle(label = ggt, subtitle = annotationtxt)
  }
  return(P2)
}







#' Predict a two-way table based on GGE model
#'
#' Predict the means for a genotype-vs-environment trial based on a Genotype
#' plus Genotype-vs-Environment interaction (GGE) model.
#'
#' This function is used to predict the response variable of a two-way table
#' (for examples the yielding of g genotypes in e environments) based on GGE
#' model. This prediction is based on the number of principal components used.
#' For more details see Yan and Kang (2007).
#'
#' @param object An object of class \code{gge}.
#' @param naxis The the number of principal components to be used in the
#'   prediction. Generally, two axis may be used. In this case, the estimated
#'   values will be those shown in the biplot.
#' @param output The type of output. It must be one of the \code{'long'}
#'   (default) returning a long-format table with the columns for environment
#'   (ENV), genotypes (GEN) and response variable (Y); or \code{'wide'} to
#'   return a two-way table with genotypes in the row, environments in the
#'   columns, filled by the estimated values.
#' @param ... Currently not used.
#' @return A two-way table with genotypes in rows and environments in columns if
#'   \code{output = "wide"} or a long format (columns ENV, GEN and Y) if
#'   \code{output = "long"} with the predicted values by the GGE model.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Yan, W., and M.S. Kang. 2003. GGE biplot analysis: a graphical
#'   tool for breeders, geneticists, and agronomists. CRC Press.
#' @method predict gge
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' mod <- gge(data_ge, GEN, ENV, c(GY, HM))
#' predict(mod)
#' }
#'
predict.gge <- function(object, naxis = 2, output = "wide", ...) {
  if (has_class(object, "gtb")) {
    stop("The object must be of class 'gge'.")
  }
  listres <- list()
  varin <- 1
  for (var in 1:length(object)) {
    objectin <- object[[var]]
    if (naxis > min(dim(objectin$coordenv))) {
      stop("The number of principal components cannot be greater than min(g, e), in this case ",
           min(dim(objectin$coordenv)))
    }
    # SVP
    if (objectin$svp == "environment" | objectin$svp == 2) {
      pred <- (objectin$coordgen[, 1:naxis] * (objectin$d)) %*%
        t(objectin$coordenv[, 1:naxis])
    }
    if (objectin$svp == "genotype" | objectin$svp == 1) {
      pred <- (objectin$coordgen[, 1:naxis] %*% t(objectin$coordenv[,
                                                                    1:naxis] * (objectin$d)))
    }
    if (objectin$svp == "symmetrical" | objectin$svp == 3) {
      pred <- (objectin$coordgen[, 1:naxis] %*% t(objectin$coordenv[,
                                                                    1:naxis]))
    }
    # Scaling
    if (objectin$scaling == "sd" | objectin$scaling == 1) {
      pred <- sweep(pred, 2, objectin$scale_val, FUN = "*")
    }
    # Centering
    if (objectin$centering == "global" | objectin$centering == 1) {
      pred <- pred + objectin$grand_mean
    }
    if (objectin$centering == "environment" | objectin$centering ==
        2) {
      pred <- sweep(pred, 2, objectin$mean_env, FUN = "+")
    }
    if (objectin$centering == "double" | objectin$centering == 3) {
      for (i in 1:nrow(pred)) {
        for (j in 1:ncol(pred)) {
          pred[i, j] <- pred[i, j] - objectin$grand_mean +
            objectin$mean_env[j] + objectin$mean_gen[i]
        }
      }
    }
    rownames(pred) <- objectin$labelgen
    colnames(pred) <- objectin$labelenv
    if (output == "wide") {
      temp <- as_tibble(pred, rownames = NA)
    }
    if (output == "long") {
      temp <-
        pred %>%
        make_long() %>%
        as_tibble()
    }
    listres[[paste(names(object[var]))]] <- temp
  }
  return(listres)
}
