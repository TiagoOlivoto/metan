#' Reorder a correlation matrix
#' @description
#' `r badge('stable')`
#'
#' Reorder the correlation matrix according to the correlation coefficient by
#' using  hclust for hierarchical clustering order. This is useful to identify
#' the hidden pattern in the matrix.
#'
#' @param x The correlation matrix
#'
#' @return The ordered correlation matrix
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(metan)
#' cor_mat <- corr_coef(data_ge2, PH, EH, CD, CL, ED, NKR)
#' cor_mat$cor
#' reorder_cormat(cor_mat$cor)
#' }
#'
reorder_cormat <- function(x){
  if(!is.matrix(x) | nrow(x) != ncol(x)){
    stop("object 'x' must be a square matrix.", call. = FALSE)
  }
  hc <- hclust(as.dist((1 - x) / 2))
  x <- x[hc$order, hc$order]
  return(x)
}
NULL



#' Linear and partial correlation coefficients
#'
#' Computes Pearson's linear correlation or partial correlation with p-values
#'
#' The partial correlation coefficient is a technique based on matrix operations
#' that allow us to identify the association between two variables by removing
#' the effects of the other set of variables present (Anderson 2003) A
#' generalized way to estimate the partial correlation coefficient between two
#' variables (i and j ) is through the simple correlation matrix that involves
#' these two variables and m other variables from which we want to remove the
#' effects. The estimate of the partial correlation coefficient between i and j
#' excluding the effect of m other variables is given by:
#' \loadmathjax
#' \mjsdeqn{r_{ij.m} = \frac{{- {a_{ij}}}}{{\sqrt {{a_{ii}}{a_{jj}}}}}}
#'
#' Where \mjseqn{r_{ij.m}} is the partial correlation coefficient between
#' variables i and j, without the effect of the other m variables;
#' \mjseqn{a_{ij}} is the ij-order element of the inverse of the linear
#' correlation matrix; \mjseqn{a_{ii}}, and \mjseqn{a_{jj}} are the elements of
#' orders ii and jj, respectively, of the inverse of the simple correlation
#' matrix.
#'
#' @param data The data set. It understand grouped data passed from
#'   [dplyr::group_by()].
#' @param ... Variables to use in the correlation. If no variable is informed
#'   all the numeric variables from `data` are used.
#' @param type The type of correlation to be computed. Defaults to `"linear"`.
#'   Use `type = "partial"` to compute partial correlation.
#' @param method a character string indicating which partial correlation
#'   coefficient is to be computed. One of "pearson" (default), "kendall", or
#'   "spearman"
#' @param use an optional character string giving a method for computing
#'   covariances in the presence of missing values. See [stats::cor] for more
#'   details
#'@param by One variable (factor) to compute the function by. It is a shortcut
#'  to [dplyr::group_by()].This is especially useful, for example,
#'  to compute correlation matrices by levels of a factor.
#' @param verbose Logical argument. If `verbose = FALSE` the code is run
#'  silently.
#' @return A list with the correlation coefficients and p-values
#' @importFrom dplyr between
#' @references
#' Anderson, T. W. 2003. An introduction to multivariate statistical analysis.
#' 3rd ed. Wiley-Interscience.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#' \donttest{
#' library(metan)
#'
#' # All numeric variables
#' all <- corr_coef(data_ge2)
#'
#'
#' # Select variable
#' sel <-
#'   corr_coef(data_ge2,
#'             EP, EL, CD, CL)
#' sel$cor
#'
#' # Select variables, partial correlation
#' sel <-
#'   corr_coef(data_ge2,
#'             EP, EL, CD, CL,
#'             type = "partial")
#' sel$cor
#'
#' }
#'
corr_coef <- function(data,
                      ...,
                      type = c("linear", "partial"),
                      method = c("pearson", "kendall", "spearman"),
                      use = c("pairwise.complete.obs", "everything", "complete.obs"),
                      by = NULL,
                      verbose = TRUE){


  if (!missing(by)){
    if(length(as.list(substitute(by))[-1L]) != 0){
      stop("Only one grouping variable can be used in the argument 'by'.\nUse 'group_by()' to pass '.data' grouped by more than one variable.", call. = FALSE)
    }
    data <- group_by(data, {{by}})
  }
  if(is_grouped_df(data)){
    corr_coef <-
      data %>%
        doo(corr_coef,
            ...,
            type = type,
            method = method,
            use = use,
            verbose = verbose)
      return(set_class(corr_coef, c("tbl_df",  "corr_coef_group", "tbl",  "data.frame")))
    } else{
      type <- rlang::arg_match(type)
      method <- rlang::arg_match(method)
      use <- rlang::arg_match(use)

      if(missing(...)){
        x <- select_numeric_cols(data)
      }
      if(!missing(...)){
        x <- select(data, ...) %>% select_numeric_cols()
      }
      if(has_na(x)){
        rlang::warn("Missing values in the data. Option 'pairwise.complete.obs' used to compute correlation with all complete pairs of observations.")
      }
      if(type == "linear"){
        apply_r <- function(A, FUN, ...) {
          mapply(function(a, B)
            lapply(B, function(x) FUN(a, x, ...)),
            a = A,
            MoreArgs = list(B = A))
        }
        pval <- suppressWarnings(apply(apply_r(x, cor.test, method = method), 1:2, function(x) x[[1]]$p.value))
        r <-  cor(x, method = method, use = use)

        # apply(apply_r(x, cor, method = method, use = use), 1:2, function(x) x[[1]])
      } else{
        # compute p-values
        get_pval <- function(r, n, nvar, df){
          pval_mat <- r
          for(i in 1:nrow(r)){
            for(j in 1:ncol(r)){
              if(i == j){
                pval_mat[i, j] <- 0
              }
              pval_mat[i, j] <- 2 * (1 - pt(abs(r[i, j]/(sqrt(1 - r[i, j]^2)) * sqrt(n - nvar)), df = df))
            }
          }
          return(pval_mat)
        }

        cor.x <- cor(x, method = method, use = use)
        n <- nrow(x)
        nvar <- ncol(cor.x)
        df <- n - nvar
        if (df < 0) {
          warning("The number of variables is higher than the number of individuals. Hypothesis testing will not be made.",
                  call. = FALSE)
        }
        m <- as.matrix(cor.x)
        X.resid <- -(solve_svd(m))
        diag(X.resid) <- 1/(1 - (1 - 1/diag(solve_svd(m))))
        r <- cov2cor(X.resid)
        pval <- get_pval(r, n, nvar, df)
      }
      return(structure(list(cor = r, pval = pval), class = "corr_coef"))
    }
  }



  #' Print an object of class corr_coef
  #'
  #' Print the `corr_coef` object in two ways. By default, the results
  #' are shown in the R console. The results can also be exported to the directory
  #' into a *.txt file.
  #'
  #'
  #' @param x An object of class `corr_coef`
  #' @param export A logical argument. If `TRUE`, a *.txt file is exported to
  #'   the working directory.
  #' @param file.name The name of the file if `export = TRUE`
  #' @param digits The significant digits to be shown.
  #' @param ... Options used by the tibble package to format the output. See
  #'   [tibble::formatting()] for more details.
  #' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
  #' @method print corr_coef
  #' @export
  #' @examples
  #' \donttest{
  #' library(metan)
  #' sel <- corr_coef(data_ge2, EP, EL, CD, CL)
  #' print(sel)
  #' }
  print.corr_coef <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
    if (export == TRUE) {
      file.name <- ifelse(is.null(file.name) == TRUE, "corr_coef print", file.name)
      sink(paste0(file.name, ".txt"))
    }
    cat("---------------------------------------------------------------------------\n")
    cat("Pearson's correlation coefficient\n")
    cat("---------------------------------------------------------------------------\n")
    print(x$cor, digits = digits)
    cat("---------------------------------------------------------------------------\n")
    cat("p-values for the correlation coefficients\n")
    cat("---------------------------------------------------------------------------\n")
    print(x$pval, digits = digits)
    cat("\n\n\n")
    if (export == TRUE) {
      sink()
    }
  }
  NULL


  #' Create a correlation heat map
  #'
  #' Create a correlation heat map for object of class `corr_coef`
  #'
  #' @param x The data set.
  #' @param type The type of heat map to produce. Either `lower` (default) to
  #'   produce a lower triangle heat map or `upper` to produce an upper
  #'   triangular heat map.
  #' @param diag Plot diagonal elements? Defaults to `FALSE`.
  #' @param reorder Reorder the correlation matrix to identify the hidden pattern?
  #'   Defaults to `TRUE`.
  #' @param signif How to show significant correlations. If `"stars"` is used
  #'   (default), stars are used showing the significance at 0.05 ("*"), 0.01
  #'   ("**") and 0.001 ("***") probability error. If `signif = "pval"`, then
  #'   the p-values are shown.
  #' @param show The correlations to show. Either `all` (default) or `signif`
  #'   (only significant correlations).
  #' @param p_val The p-value to the correlation significance.
  #' @param caption Logical. If `TRUE` (Default) includes a caption with the
  #'   significance meaning for stars.
  #' @param digits.cor,digits.pval The significant digits to show for correlations
  #'   and p-values, respectively.
  #' @param col.low,col.mid,col.high The color for the low (-1), mid(0) and high
  #'   (1) points in the color key. Defaults to `blue`, `white`, and
  #'   `red`, respectively.
  #' @param lab.x.position,lab.y.position The position of the x and y axis label.
  #'   Defaults to `"bottom"` and `"right"` if `type = "lower"` or
  #'   `"top"` and `"left"` if `type = "upper"`.
  #' @param legend.position The legend position in the plot.
  #' @param legend.title The title of the color key. Defaults to `"Pearson's
  #'   Correlation"`.
  #' @param size.text.cor The size of the text for correlation values. Defaults to 3.
  #' @param size.text.signif The size of the text for significance values (stars or p-values). Defaults to 3.
  #' @param size.text.lab The size of the text for labels. Defaults to 10.
  #' @param ... Currently not used.
  #' @method plot corr_coef
  #' @return An object of class `gg, ggplot`
  #' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
  #' @export
  #' @examples
  #' \donttest{
  #' library(metan)
  #' # All numeric variables
  #' x <- corr_coef(data_ge2)
  #' plot(x)
  #' plot(x, reorder = FALSE)
  #'
  #' # Select variables
  #' sel <- corr_coef(data_ge2, EP, EL, CD, CL)
  #' plot(sel,
  #'      type = "upper",
  #'      reorder = FALSE,
  #'      size.text.lab = 14,
  #'      size.text.plot = 5)
  #' }

  plot.corr_coef <- function(x,
                             type = "lower",
                             diag = FALSE,
                             reorder = TRUE,
                             signif = c("stars", "pval"),
                             show = c("all", "signif"),
                             p_val = 0.05,
                             caption = TRUE,
                             digits.cor = 2,
                             digits.pval = 3,
                             col.low = "red",
                             col.mid = "white",
                             col.high = "blue",
                             lab.x.position = NULL,
                             lab.y.position = NULL,
                             legend.position = NULL,
                             legend.title = "Pearson's\nCorrelation",
                             size.text.cor = 3,
                             size.text.signif = 3,
                             size.text.lab = 10,
                             ...){
    signif <- rlang::arg_match(signif)
    show <- rlang::arg_match(show)

    pval <- x$pval
    correl <- x$cor
    if(reorder == TRUE){
      correl <- reorder_cormat(correl)
      pval <- pval[rownames(correl), colnames(correl)]
    }

    nonsign <- which(pval >= p_val)
    if(show == "signif"){
      correl[nonsign] <- NA
    }

    if(type == "lower"){
      pval <- make_lower_tri(pval)
      correl <- make_lower_tri(correl)
    }
    if(type == "upper"){
      pval <- make_upper_tri(pval)
      correl <- make_upper_tri(correl)
    }
    bind_data <-
      expand.grid(dimnames(correl)) %>%
      mutate(cor = as.vector(correl),
             pval = as.vector(pval),
             pval_star = pval) %>%
      set_names("v1", "v2", "cor", "pval", "pval_star") %>%
      dplyr::filter(!is.na(pval))
    if(signif == "stars"){
      bind_data <-
        bind_data |>
        mutate(pval_star = case_when(pval < 0.001 ~ "***",
                                     between(pval, 0.001, 0.01) ~ "**",
                                     between(pval, 0.01, 0.05) ~ "*",
                                     pval >= 0.05 ~ "ns")) %>%
        set_names("v1", "v2", "cor", "pval", "pval_star") %>%
        dplyr::filter(!is.na(pval))
    } else{
      bind_data <-
        expand.grid(dimnames(correl)) %>%
        mutate(cor = as.vector(correl),
               pval = as.vector(pval),
               pval_star <- case_when(
                 pval < 0.001 ~ "***",
                 between(pval, 0.001, 0.01) ~ "**",
                 between(pval, 0.01, 0.05) ~ "*",
                 pval >= 0.05 ~ "ns"
               )) %>%
        set_names("v1", "v2", "cor", "pval", "pval_star") %>%
        dplyr::filter(!is.na(pval))
    }
    if(show == "signif"){
      bind_data <-
        bind_data |>
        mutate(pval_star = ifelse(pval >= p_val, NA, pval_star),
               pval = ifelse(pval >= p_val, NA, pval))
    }

    if(diag == FALSE){
      bind_data <- dplyr::filter(bind_data, !v1 == v2)
    }
    if(type == "lower"){
      lab.y.position <- ifelse(missing(lab.y.position), "right", lab.y.position)
      lab.x.position <- ifelse(missing(lab.x.position), "bottom", lab.x.position)
    }
    if(type == "upper"){
      lab.y.position <- ifelse(missing(lab.y.position), "left", lab.y.position)
      lab.x.position <- ifelse(missing(lab.x.position), "top", lab.x.position)
    }
    axis.text.x.hjust <- ifelse(lab.x.position == "bottom", 1, 0)
    if(missing(legend.position) & type == "lower"){
      legend.position <- c(0.2, 0.85)
    }
    if(missing(legend.position) & type == "upper"){
      legend.position <- c(0.8, 0.17)
    }
    if(!missing(legend.position)){
      legend.position <- legend.position
    }


    p <-
      ggplot(bind_data, aes(v1, v2, fill = cor)) +
      geom_tile(aes(fill = cor),
                colour = "gray") +
      geom_text(aes(label = ifelse(!is.na(cor), format(round(cor, digits.cor), nsmall = digits.cor), "")),
                vjust = 0,
                size = size.text.cor) +
      {if(signif == "stars")geom_text(aes(label = ifelse(!is.na(pval_star), pval_star, "")),
                                      vjust = 1.5,
                                      size = size.text.signif)} +
      {if(signif == "pval")geom_text(aes(label = ifelse(!is.na(pval), as.vector(format(signif(pval, digits = digits.pval), nsmall = digits.pval)), "")),
                                     vjust = 1.5,
                                     size = size.text.signif)} +
      scale_fill_gradient2(low = col.low,
                           high = col.high,
                           mid = col.mid,
                           midpoint = 0,
                           limit = c(-1, 1),
                           space = "Lab",
                           na.value = "transparent",
                           name = legend.title)+
      theme_bw() +
      xlab(NULL) +
      ylab(NULL) +
      scale_y_discrete(expand = expansion(mult = c(0,0)),
                       position = lab.y.position)+
      scale_x_discrete(expand = expansion(mult = c(0,0)),
                       position = lab.x.position) +
      coord_fixed() +
      theme(axis.ticks = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.background = element_blank(),
            legend.position = legend.position,
            legend.direction = "horizontal",
            axis.text = element_text(color = "black", size = size.text.lab),
            axis.text.x = element_text(angle = 45,
                                       vjust = 1,
                                       hjust = axis.text.x.hjust))+
      guides(fill = guide_colorbar(barwidth = 6,
                                   barheight = 1,
                                   title.position = "top",
                                   title.hjust = 0.5))
    if(signif == "stars" & caption == TRUE){
      p <-
        p + labs(caption = c("ns p >= 0.05; * p < 0.05; ** p < 0.01; and *** p < 0.001"))
    }
    suppressWarnings(return(p))
  }
  NULL



  #' Network plot of a correlation matrix
  #'
  #' Produces a network plot of a correlation matrix or an object computed with
  #' [corr_coef()]. Variables that are more highly correlated appear closer
  #' together and are joined by stronger (more opaque) and wider paths.  The proximity of the
  #' points is determined using multidimensional clustering, also known as
  #' principal coordinates analysis (Gower, 1966). The color of the paths also
  #' indicates the sign of the correlation (blue for positive and red for
  #' negative).
  #'
  #' @param model A model computed with [corr_coef()] or a symmetric matrix, often
  #'   produced with [stats::cor()].
  #' @param min_cor Number to indicate the minimum value of correlations to plot
  #'   (0-1 in absolute terms). By default, all the correlations are plotted when
  #'   `model` is a matrix, and significant correlations (p-value < 0.05) when
  #'   `model` is an object computed with [corr_coef()].
  #' @param show The correlations to be shown when `model` is an object computed
  #'   with [corr_coef()]. Either `"signif"` (default) to show only significant
  #'   correlations or `"all"` to show all the correlations.
  #' @param p_val The p-value to indicate significant correlations. Defaults to
  #'   `0.05`.
  #' @param legend The type of legend. Either `"full"` (ranges from -1 to +1) or
  #'   `"range"` (ranges according to the data range). Defaults to `"full"`.
  #' @param colours A vector of colors to use for n-color gradient.
  #' @param legend_width The width of the legend (considering `position = "right"`)
  #' @param legend_height The height of the legend (considering `position = "right"`)
  #' @param legend_position The legend position. Defaults to `"right"`.
  #' @param curved Shows curved paths? Defaults to `TRUE`.
  #' @param angle A numeric value between 0 and 180, giving an amount to skew the
  #'   control points of the curve. Values less than 90 skew the curve towards the
  #'   start point and values greater than 90 skew the curve towards the end
  #'   point.
  #' @param curvature A numeric value giving the amount of curvature. Negative
  #'   values produce left-hand curves, positive values produce right-hand curves,
  #'   and zero produces a straight line.
  #' @param expand_x,expand_y Vector of multiplicative range expansion factors. If
  #'   length 1, both the lower and upper limits of the scale are expanded
  #'   outwards by mult. If length 2, the lower limit is expanded by `mult[1]` and
  #'   the upper limit by `mult[2]`.
  #' @references
  #' Gower, J.C. 1966. Some Distance Properties of Latent Root and Vector Methods
  #' Used in Multivariate Analysis. Biometrika 53(3/4): 325–338.
  #' \doi{10.2307/2333639}
  #' @return A `ggplot` object
  #' @export
  #'
  #' @examples
  #' cor <- corr_coef(iris)
  #' network_plot(cor)
  #' network_plot(cor,
  #'              show = "all",
  #'              curved = FALSE,
  #'              legend_position = "bottom",
  #'              legend = "range")
  #'
  network_plot <- function(model,
                           min_cor = NULL,
                           show = c("signif", "all"),
                           p_val = 0.05,
                           legend = c("full", "range"),
                           colours = c("red", "white", "blue"),
                           legend_width = 1,
                           legend_height = 15,
                           legend_position = c("right", "left", "top", "bottom"),
                           curved = TRUE,
                           angle = 90,
                           curvature = 0.5,
                           expand_x = 0.25,
                           expand_y = 0.25){
    # Adapted from https://github.com/tidymodels/corrr/blob/683687e5ff061d6094ef09e14006b218a02c3f06/R/cor_df.R
    legend <- rlang::arg_match(legend)
    show <- rlang::arg_match(show)
    legend_position <- rlang::arg_match(legend_position)
    curvature <- ifelse(curved == TRUE, curvature, 0)

    # check the class of the model
    if(metan::has_class(model, "corr_coef")){
      correls <- model$cor
      pvals <- model$pval
      distance <- 1 - abs(correls)
    } else{
      if(isSymmetric(model)){
        dimmat <- dim(model)
        correls <- model
        pvals <- matrix(rep(0, prod(dimmat)), nrow = dimmat[[1]], ncol = dimmat[[2]])
        distance <- 1 - abs(correls)
        rlang::warn("It seems that a correlation matrix was used, so, all correlations are shown. Use `min_cor` to define a cut point.")
      }
    }
    if (!is.null(min_cor) ) {
      if(min_cor < 0 || min_cor > 1){
        rlang::abort("min_cor must be a value ranging from zero to one.")
      }
    }
    points <- if (ncol(correls) == 1) {
      # 1 var: a single central point
      matrix(c(0, 0), ncol = 2, dimnames = list(colnames(correls)))
    } else if (ncol(correls) == 2) {
      # 2 vars: 2 opposing points
      matrix(c(0, -0.1, 0, 0.1), ncol = 2, dimnames = list(colnames(correls)))
    } else {
      # More than 2 vars: uses stats::cmdscale()
      suppressWarnings(stats::cmdscale(distance, k = 2))
    }

    if (ncol(points) < 2) {
      cont_flag <- FALSE
      shift_matrix <- matrix(1,
                             nrow = nrow(correls),
                             ncol = ncol(correls)
      )
      diag(shift_matrix) <- 0

      for (shift in 10^(-6:-1)) {
        shifted_distance <- distance + shift * shift_matrix
        points <- suppressWarnings(stats::cmdscale(shifted_distance))
        if (ncol(points) > 1) {
          cont_flag <- TRUE
          break
        }
      }
      if (!cont_flag){
        rlang::abort("Can't generate network plot.\nAttempts to generate 2-d coordinates failed.")
      }
      rlang::warn("Plot coordinates derived from correlation matrix have dimension < 2.\nPairwise distances have been adjusted to facilitate plotting.")
    }
    points <- data.frame(points)
    colnames(points) <- c("x", "y")
    points$id <- rownames(points)

    # matrix of p-values
    pval_mat <- make_lower_tri(pvals)
    nonsign <- which(pval_mat >= p_val)

    # Create a proximity matrix of the paths to be plotted.
    proximity <- abs(correls) |> make_lower_tri()

    # show only significant correlations
    if(show == "signif"){
      proximity[nonsign] <- NA
    }
    # define a cut point for the correlations
    if(!is.null(min_cor)){
      proximity[proximity < min_cor] <- NA
    }

    # Produce a data frame of data needed for plotting the paths.
    n_paths <- sum(!is.na(proximity))
    paths <- data.frame(matrix(nrow = n_paths, ncol = 6))
    colnames(paths) <- c("x", "y", "xend", "yend", "proximity", "sign")
    path <- 1
    for (i in 1:nrow(proximity)) {
      for (j in 1:ncol(proximity)) {
        path_proximity <- proximity[i, j]
        if (!is.na(path_proximity)) {
          path_sign <- sign(correls[i, j])
          x <- points$x[i]
          y <- points$y[i]
          xend <- points$x[j]
          yend <- points$y[j]
          paths[path, ] <- c(x, y, xend, yend, path_proximity, path_sign)
          path <- path + 1
        }
      }
    }

    if(legend == "full"){
      leg_range = c(-1, 1)
    } else{
      leg_range = c(min(correls[row(correls)!=col(correls)]),
                    max(correls[row(correls)!=col(correls)]))
    }
    # return(paths)
    # create the plot
    ggplot() +
      geom_curve(data = paths,
                 aes(x = x,
                     y = y,
                     xend = xend,
                     yend = yend,
                     alpha = proximity * sign,
                     size = proximity,
                     colour = proximity * sign),
                 ncp = 200,
                 angle = angle,
                 curvature = curvature,
                 lineend = "round") +
      scale_alpha(limits = c(0, 1)) +
      scale_size(limits = c(0, 1)) +
      scale_colour_gradientn(limits = leg_range, colors = colours) +
      geom_point(data = points,
                 aes(x, y),
                 size = 3,
                 shape = 19,
                 colour = "white",
                 alpha = 0.3) +
      ggrepel::geom_text_repel(data = points,
                               aes(x, y, label = id),
                               fontface = "bold",
                               size = 5,
                               segment.size = 0.0,
                               segment.color = "transparent") +
      scale_x_continuous(expand = expansion(expand_x)) +
      scale_y_continuous(expand = expansion(expand_y)) +
      theme_void() +
      theme(plot.margin=unit(c(0, 15, 0, 0), 'pt'),
            legend.position = legend_position) +
      {if(legend_position %in% c("left", "right"))
        guides(size = "none",
               alpha = "none",
               color = guide_colourbar(draw.ulim = TRUE,
                                       draw.llim = TRUE,
                                       frame.colour = "black",
                                       ticks = TRUE,
                                       nbin = 10,
                                       barwidth = legend_width,
                                       barheight = legend_height,
                                       direction = 'vertical'))} +
      {if(legend_position %in% c("top", "bottom"))
        guides(size = "none",
               alpha = "none",
               color = guide_colourbar(draw.ulim = TRUE,
                                       draw.llim = TRUE,
                                       frame.colour = "black",
                                       ticks = TRUE,
                                       nbin = 10,
                                       barwidth = legend_height,
                                       barheight = legend_width,
                                       direction = 'horizontal'))} +
      {if (legend != "none") labs(colour = NULL)} +
      {if (legend == "none") theme(legend.position = "none")}
  }



  #' Focus on section of a correlation matrix
  #'
  #' Select a set of variables from a correlation matrix to keep as the columns,
  #' and exclude these or all other variables from the rows.
  #'
  #' @param model A model computed with [corr_coef()] or a symmetric matrix, often
  #'   produced with [stats::cor()].
  #' @param ... One or more unquoted variable name separated by commas. Variable
  #'   names can be used as if they were positions in the data frame, so
  #'   expressions like `x:y` can be used to select a range of variables.
  #'
  #' @return A tibble
  #' @export
  #'
  #' @examples
  #' corr_coef(data_ge2) |> corr_focus(PH)
  corr_focus <- function(model, ...) {
    if(metan::has_class(model, "corr_coef")){
      x <- model$cor |> as.data.frame()
    } else{
      if(is.data.frame(model)){
        rlang::abort("Invalid format. `model` must be a correlation matrix or an object computed with `corr_coef().`")
      }
      x <- as.data.frame(model)
    }
    x <- dplyr::select(x, ...)
    row_name <- rownames(x)
    x <- rownames_to_column(x, "var")
    vars <- colnames(x)[-1]
    dplyr::filter(x, !var  %in% vars) |>
      as_tibble()
  }
