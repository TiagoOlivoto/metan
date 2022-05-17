#' Multi-trait mean performance and stability index
#' @description
#' `r badge('experimental')`
#'
#' Computes the multi-trait stability index proposed by Olivoto et al. (2019)
#' considering different parametric and non-parametric stability indexes.
#'
#'
#' @param model An object of class `mtmps`
#' @param SI An integer (0-100). The selection intensity in percentage of the
#' total number of genotypes.
#' @param mineval The minimum value so that an eigenvector is retained in the
#' factor analysis.
#' @param verbose If `verbose = TRUE` (Default), some results are shown in the
#'   console.
#' @return An object of class `mtmps` with the following items:
#' * **data** The data used to compute the factor analysis.
#' * **cormat** The correlation matrix among the environments.
#' * **PCA** The eigenvalues and explained variance.
#' * **FA** The factor analysis.
#' * **KMO** The result for the Kaiser-Meyer-Olkin test.
#' * **MSA** The measure of sampling adequacy for individual variable.
#' * **communalities** The communalities.
#' * **communalities_mean** The communalities' mean.
#' * **initial_loadings** The initial loadings.
#' * **finish_loadings** The final loadings after varimax rotation.
#' * **canonical_loadings** The canonical loadings.
#' * **scores_gen** The scores for genotypes in all retained factors.
#' * **scores_ide** The scores for the ideotype in all retained factors.
#' * **MTSI** The multi-trait mean performance and stability index.
#' * **contri_fac** The relative contribution of each factor on the MTSI
#' value. The lower the contribution of a factor, the close of the ideotype the
#' variables in such factor are.
#' * **contri_fac_rank, contri_fac_rank_sel** The rank for the contribution
#' of each factor for all genotypes and selected genotypes, respectively.
#' * **sel_dif_trait, sel_dif_stab, sel_dif_mps** A data frame containing the
#' selection differential (gains) for the mean performance, stability index, and
#' mean performance and stability index, respectively. The following variables
#' are shown.
#'   - `VAR`: the trait's name.
#'   - `Factor`: The factor that traits where grouped into.
#'   - `Xo`: The original population mean.
#'   - `Xs`: The mean of selected genotypes.
#'   - `SD` and `SDperc`: The selection differential and selection differential in
#'   percentage, respectively.
#'   - `h2`: The broad-sense heritability.
#'   - `SG` and `SGperc`: The selection gains and selection gains in percentage,
#'   respectively.
#'   - `sense`: The desired selection sense.
#'   - `goal`: selection gains match desired sense? 100 for yes and 0 for no.
#' * **stat_dif_trait, stat_dif_stab, stat_dif_mps** A data frame with the
#' descriptive statistic for the selection gains for the mean performance,
#' stability index, and mean performance and stability index, respectively. The
#' following columns are shown by sense.
#'    - `sense`: The desired selection sense.
#'    - `variable`: the trait's name.
#'    - `min`: the minimum value for the selection gain.
#'    - `mean`: the mean value for the selection gain.
#'    - `ci`: the confidence interval for the selection gain.
#'    - `sd.amo`: the standard deviation for the selection gain.
#'    - `max`: the maximum value for the selection gain.
#'    - `sum`: the sum of the selection gain.
#' * **sel_gen** The selected genotypes.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @references Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, B.G. Sari, and
#'   M.I. Diel. 2019. Mean performance and stability in multi-environment trials
#'   II: Selection based on multiple traits. Agron. J. 111:2961-2969.
#' \doi{10.2134/agronj2019.03.0220}
#' @seealso [mgidi()], [mps()], [get_model_data()]
#' @examples
#' \donttest{
#' library(metan)
#' # The same approach as mtsi()
#' # mean performance and stability for GY and HM
#' # mean performance: The genotype's BLUP
#' # stability: the WAASB index (lower is better)
#' # weights: equal for mean performance and stability
#'
#' model <-
#' mps(data_ge,
#'     env = ENV,
#'     gen = GEN,
#'     rep = REP,
#'     resp = everything())
#' selection <- mtmps(model)
#'
#' # gains for stability
#' selection$sel_dif_stab
#'
#' # gains for mean performance
#' selection$sel_dif_trait
#'}
mtmps <- function(model,
                  SI = 15,
                  mineval = 1,
                  verbose = TRUE) {
  if(has_class(model, "mps_group")){
    bind <-
      model %>%
      mutate(data = map(data, ~.x %>%
                          mtmps(SI = SI,
                                mineval = mineval,
                                verbose = verbose)))
    return(set_class(bind, c("tbl_df",  "mtmps_group",  "tbl",  "data.frame")))
  } else{
    data <- model[["mps_ind"]] %>% column_to_rownames("GEN")
    if(has_na(data)){
      stop("Missing values for traits ")
    }
    rescaled <- model$sense_mper
    names(rescaled) <- names(data)
    ideotype.D <- rep(100, ncol(data))
    names(ideotype.D) <- names(data)
    df_ideotype <-
      data.frame(rescaled) %>%
      rownames_to_column("VAR") %>%
      set_names("VAR", "sense")
    rescaled_stab <- model$sense_stab
    names(rescaled_stab) <- names(data)
    df_ideotype_stab <-
      data.frame(rescaled_stab) %>%
      rownames_to_column("VAR") %>%
      set_names("VAR", "sense")
    if (is.null(SI)) {
      ngs <- NULL
    } else {
      ngs <- round(nrow(data) * (SI/100), 0)
    }
    observed <- model$observed
    means <- model$mps_ind %>% column_to_rownames("GEN")
    cor.means <- cor(means)
    eigen.decomposition <- eigen(cor.means)
    eigen.values <- eigen.decomposition$values
    eigen.vectors <- eigen.decomposition$vectors
    colnames(eigen.vectors) <- paste("PC", 1:ncol(cor.means), sep = "")
    rownames(eigen.vectors) <- colnames(means)
    if (length(eigen.values[eigen.values >= mineval]) == 1) {
      eigen.values.factors <- as.vector(c(as.matrix(sqrt(eigen.values[eigen.values >= mineval]))))
      initial_loadings <- cbind(eigen.vectors[, eigen.values >= mineval] * eigen.values.factors)
      A <- initial_loadings
    } else {
      eigen.values.factors <-
        t(replicate(ncol(cor.means), c(as.matrix(sqrt(eigen.values[eigen.values >= mineval])))))
      initial_loadings <- eigen.vectors[, eigen.values >= mineval] * eigen.values.factors
      A <- varimax(initial_loadings)[[1]][]
    }
    partial <- solve_svd(cor.means)
    k <- ncol(means)
    seq_k <- seq_len(ncol(means))
    for (j in seq_k) {
      for (i in seq_k) {
        if (i == j) {
          next
        } else {
          partial[i, j] <- -partial[i, j]/sqrt(partial[i, i] * partial[j, j])
        }
      }
    }
    KMO <- sum((cor.means[!diag(k)])^2)/(sum((cor.means[!diag(k)])^2) + sum((partial[!diag(k)])^2))
    MSA <- unlist(lapply(seq_k, function(i) {
      sum((cor.means[i, -i])^2)/(sum((cor.means[i, -i])^2) + sum((partial[i, -i])^2))
    }))
    names(MSA) <- colnames(means)
    colnames(A) <- paste("FA", 1:ncol(initial_loadings), sep = "")
    pca <- tibble(PC = paste("PC", 1:ncol(means), sep = ""),
                  Eigenvalues = eigen.values,
                  `Variance (%)` = (eigen.values/sum(eigen.values)) * 100,
                  `Cum. variance (%)` = cumsum(`Variance (%)`))
    Communality <- diag(A %*% t(A))
    Uniquenesses <- 1 - Communality
    fa <- cbind(A, Communality, Uniquenesses) %>% as_tibble(rownames = NA) %>%  rownames_to_column("VAR")
    z <- scale(means, center = FALSE, scale = apply(means, 2, sd))
    canonical_loadings <- t(t(A) %*% solve_svd(cor.means))
    scores <- z %*% canonical_loadings
    colnames(scores) <- paste("FA", 1:ncol(scores), sep = "")
    pos.var.factor <- which(abs(A) == apply(abs(A), 1, max), arr.ind = TRUE)
    var.factor <- lapply(1:ncol(A), function(i) {
      rownames(pos.var.factor)[pos.var.factor[, 2] == i]
    })
    names(var.factor) <- paste("FA", 1:ncol(A), sep = "")
    names.pos.var.factor <- rownames(pos.var.factor)
    ideotypes.matrix <- t(as.matrix(ideotype.D))/apply(means, 2, sd)
    rownames(ideotypes.matrix) <- "ID1"
    ideotypes.scores <- ideotypes.matrix %*% canonical_loadings
    gen_ide <- sweep(scores, 2, ideotypes.scores, "-")
    MTSI <- sort(apply(gen_ide, 1, function(x) sqrt(sum(x^2))), decreasing = FALSE)
    contr.factor <- data.frame((sqrt(gen_ide^2)/apply(gen_ide, 1, function(x) sum(sqrt(x^2)))) * 100) %>%
      rownames_to_column("GEN") %>%
      as_tibble()
    means.factor <- means[, names.pos.var.factor]
    observed <- observed[, names.pos.var.factor]
    contri_long <- pivot_longer(contr.factor, -GEN)
    if (!is.null(ngs)) {
      selected <- names(MTSI)[1:ngs]
      sel_dif <- tibble(VAR = names(pos.var.factor[, 2]),
                        Factor = paste("FA", as.numeric(pos.var.factor[, 2])),
                        Xo = colMeans(means.factor),
                        Xs = colMeans(means.factor[names(MTSI)[1:ngs], ]),
                        SD = Xs - Xo,
                        SDperc = (Xs - Xo) / Xo * 100)
      stat_dif_mps <-
        desc_stat(sel_dif, SDperc, stats = c("min, mean, ci, sd.amo, max, sum"))
      sel_dif_mean <-
        tibble(VAR = names(pos.var.factor[, 2]),
               Factor = paste("FA", as.numeric(pos.var.factor[, 2])),
               Xo = colMeans(observed),
               Xs = colMeans(observed[names(MTSI)[1:ngs], ]),
               SD = Xs - colMeans(observed),
               SDperc = (Xs - colMeans(observed)) / colMeans(observed) * 100) %>%
        left_join(df_ideotype, by = "VAR") %>%
        mutate(sense = case_when(sense == "l" ~ "decrease",
                                 sense == "a" ~ "average",
                                 sense == "h" ~ "increase"),
               goal = case_when(
                 sense == "decrease" & SDperc < 0 ~ 100,
                 sense == "increase" & SDperc > 0 ~ 100,
                 sense == "average" & SDperc == 0 ~ 100,
                 TRUE ~ 0
               )) %>%
        left_join(model$h2, by = "VAR") %>%
        add_cols(SG = SD * h2,
                 SGperc = SG / Xo * 100,
                 .after = "SDperc") %>%
        reorder_cols(h2, .after  = "SDperc")
      stat_gain <-
        desc_stat(sel_dif_mean,
                  by = sense,
                  any_of(c("SDperc", "SGperc")),
                  stats = c("min, mean, ci, sd.amo, max, sum"))
      waasb_index <- model$stability %>% rownames_to_column("GEN")
      waasb_selected <- colMeans(subset(waasb_index, GEN %in% selected) %>% select_numeric_cols())
      sel_dif_stab <-
        tibble(
          VAR = names(waasb_selected),
          Xo = colMeans(waasb_index %>% select_numeric_cols()),
          Xs = waasb_selected,
          SD = Xs - Xo,
          SDperc = (Xs - Xo) / Xo * 100) %>%
        left_join(df_ideotype_stab, by = "VAR") %>%
        mutate(sense = case_when(sense == "l" ~ "decrease",
                                 sense == "a" ~ "average",
                                 sense == "h" ~ "increase"),
               goal = case_when(
                 sense == "decrease" & SDperc < 0 ~ 100,
                 sense == "increase" & SDperc > 0 ~ 100,
                 sense == "average" & SDperc == 0 ~ 100,
                 TRUE ~ 0
               ))
      stat_dif_stab <-
        desc_stat(sel_dif_stab, SDperc,
                  stats = c("min, mean, ci, sd.amo, max, sum"))
      contri_fac_rank_sel <-
        contri_long %>%
        subset(GEN %in% selected) %>%
        ge_winners(name, GEN, value, type = "ranks", better = "l") %>%
        split_factors(ENV) %>%
        map_dfc(~.x %>% pull())
    }
    if (is.null(ngs)) {
      stat_dif_stab <- NULL
      stat_dif_mps <- NULL
      sel_dif <- NULL
      sel_dif_stab <- NULL
      sel_dif_mean <- NULL
      selected <- NULL
      contri_fac_rank_sel <- NULL
    }
    if (verbose) {
      cat("\n-------------------------------------------------------------------------------\n")
      cat("Principal Component Analysis\n")
      cat("-------------------------------------------------------------------------------\n")
      print(pca)
      cat("-------------------------------------------------------------------------------\n")
      cat("Factor Analysis - factorial loadings after rotation-\n")
      cat("-------------------------------------------------------------------------------\n")
      print(fa)
      cat("-------------------------------------------------------------------------------\n")
      cat("Comunalit Mean:", mean(Communality), "\n")
      cat("-------------------------------------------------------------------------------\n")
      if (!is.null(ngs)) {
        cat("Selection differential for the mean performance and stability index\n")
        cat("-------------------------------------------------------------------------------\n")
        print(sel_dif)
        cat("-------------------------------------------------------------------------------\n")
        cat("Selection differential for the mean of the variables\n")
        cat("-------------------------------------------------------------------------------\n")
        print(sel_dif_mean)
        cat("------------------------------------------------------------------------------\n")
        cat("Selected genotypes\n")
        cat("-------------------------------------------------------------------------------\n")
        cat(selected)
        cat("\n-------------------------------------------------------------------------------\n")
      }
    }
    contri_fac_rank <-
      contri_long %>%
      ge_winners(name, GEN, value, type = "ranks", better = "l") %>%
      split_factors(ENV) %>%
      map_dfc(~.x %>% pull())
    list(data = data,
         cormat = as.matrix(cor.means),
         PCA = pca,
         FA = fa,
         KMO = KMO,
         MSA = MSA,
         communalities = Communality,
         communalities_mean = mean(Communality),
         initial_loadings = data.frame(initial_loadings) %>% rownames_to_column("VAR") %>% as_tibble(),
         finish_loadings = data.frame(A) %>% rownames_to_column("VAR") %>% as_tibble(),
         canonical_loadings = data.frame(canonical_loadings) %>% rownames_to_column("VAR") %>% as_tibble(),
         scores_gen = data.frame(scores) %>% rownames_to_column("GEN") %>% as_tibble(),
         scores_ide = data.frame(ideotypes.scores) %>% rownames_to_column("GEN") %>% as_tibble(),
         MTSI = as_tibble(MTSI, rownames = NA) %>% rownames_to_column("Genotype") %>% rename(MTSI = value),
         contri_fac = contr.factor,
         contri_fac_rank = contri_fac_rank,
         contri_fac_rank_sel = contri_fac_rank_sel,
         sel_dif_trait = sel_dif_mean,
         stat_dif_trait = stat_gain,
         sel_dif_stab = sel_dif_stab,
         stat_dif_stab = stat_dif_stab,
         sel_dif_mps = sel_dif,
         stat_dif_mps = stat_dif_mps,
         sel_gen = selected) %>%
      set_class("mtmps") %>%
      return()
  }
}



#' Plot the multi-trait stability index
#'
#' Makes a radar plot showing the multitrait stability index proposed by Olivoto
#' et al. (2019)
#'
#'
#' @param x An object computed with [mps()].
#' @param SI An integer (0-100). The selection intensity in percentage of the
#'   total number of genotypes.
#' @param type The type of the plot. Defaults to `"index"`. Use `type
#'   = "contribution"` to show the contribution of each factor to the MTMPS
#'   index of the selected genotypes.
#' @param position The position adjustment when `type = "contribution"`.
#'   Defaults to `"fill"`, which shows relative proportions at each trait
#'   by stacking the bars and then standardizing each bar to have the same
#'   height. Use `position = "stack"` to plot the MGIDI index for each
#'   genotype.
#' @param genotypes When `type = "contribution"` defines the genotypes to
#'   be shown in the plot. By default (`genotypes = "selected"` only
#'   selected genotypes are shown. Use `genotypes = "all"` to plot the
#'   contribution for all genotypes.)
#' @param title Logical values (Defaults to `TRUE`) to include
#'   automatically generated titles.
#' @param radar Logical argument. If true (default) a radar plot is generated
#'   after using `coord_polar()`.
#' @param arrange.label Logical argument. If `TRUE`, the labels are
#'   arranged to avoid text overlapping. This becomes useful when the number of
#'   genotypes is large, say, more than 30.
#' @param x.lab,y.lab The labels for the axes x and y, respectively. x label is
#'   set to null when a radar plot is produced.
#' @param size.point The size of the point in graphic. Defaults to 2.5.
#' @param size.line The size of the line in graphic. Defaults to 0.7.
#' @param size.text The size for the text in the plot. Defaults to 10.
#' @param width.bar The width of the bars if `type = "contribution"`.
#'   Defaults to 0.75.
#' @param n.dodge The number of rows that should be used to render the x labels.
#'   This is useful for displaying labels that would otherwise overlap.
#' @param check.overlap Silently remove overlapping labels, (recursively)
#'   prioritizing the first, last, and middle labels.
#' @param invert Logical argument. If `TRUE`, rotate the plot.
#' @param col.sel The colour for selected genotypes. Defaults to `"red"`.
#' @param col.nonsel The colour for nonselected genotypes. Defaults to `"black"`.
#' @param legend.position The position of the legend.
#' @param ... Other arguments to be passed from  [ggplot2::theme()].
#' @return An object of class `gg, ggplot`.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method plot mtmps
#' @export
#' @references Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, B.G. Sari, and M.I. Diel. 2019. Mean performance and stability in multi-environment trials II: Selection based on multiple traits. Agron. J. (in press).
#' @examples
#' \donttest{
#' library(metan)
#' model <-
#' mps(data_ge,
#'     env = ENV,
#'     gen = GEN,
#'     rep = REP,
#'     resp = everything())
#' selection <- mtmps(model)
#' plot(selection)
#'}
#'
#'
plot.mtmps <- function(x,
                       SI = 15,
                       type = "index",
                       position = "fill",
                       genotypes = "selected",
                       title = TRUE,
                       radar = TRUE,
                       arrange.label = FALSE,
                       x.lab = NULL,
                       y.lab = NULL,
                       size.point = 2.5,
                       size.line = 0.7,
                       size.text = 10,
                       width.bar = 0.75,
                       n.dodge = 1,
                       check.overlap = FALSE,
                       invert = FALSE,
                       col.sel = "red",
                       col.nonsel = "black",
                       legend.position = "bottom",
                       ...) {
  if(!type %in% c("index", "contribution")){
    stop("The argument index must be one of the 'index' or 'contribution'", call. = FALSE)
  }
  if(!genotypes %in% c("selected", "all")){
    stop("The argument 'genotypes' must be one of the 'selected' or 'all'", call. = FALSE)
  }
  if(type == "index"){
    data <- x$MTSI %>%
      add_cols(sel = "Selected")
    data[["sel"]][(round(nrow(data) * (SI/100), 0) + 1):nrow(data)] <- "Nonselected"
    cutpoint <- max(subset(data, sel == "Selected")$MTSI)
    p <-
      ggplot(data = data, aes(x = reorder(Genotype, -MTSI), y = MTSI)) +
      geom_hline(yintercept = cutpoint, col = col.sel, size = size.line) +
      geom_path(colour = "black", group = 1, size = size.line) +
      geom_point(size = size.point,
                 stroke = size.point / 10,
                 aes(fill = sel),
                 shape = 21,
                 colour = "black",
      ) +
      scale_x_discrete() +
      scale_y_reverse() +
      theme_minimal()  +
      theme(legend.position = legend.position,
            legend.title = element_blank(),
            axis.title.x = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_line(size = size.line / 2),
            axis.text = element_text(colour = "black"),
            text = element_text(size = size.text),
            ...) +
      labs(y = "Multitrait stability index") +
      scale_fill_manual(values = c(col.nonsel, col.sel))
    if (radar == TRUE) {
      if(arrange.label == TRUE){
        tot_gen <- length(unique(data$Genotype))
        fseq <- c(1:(tot_gen/2))
        sseq <- c((tot_gen/2 + 1):tot_gen)
        fang <- c(90 - 180/length(fseq) * fseq)
        sang <- c(-90 - 180/length(sseq) * sseq)
        p <-
          p +
          coord_polar()  +
          theme(axis.text.x = suppressMessages(suppressWarnings(element_text(angle = c(fang, sang)))), ...)
      } else{
        p <- p + coord_polar()
      }
    }
  } else{
    x.lab <- ifelse(!missing(x.lab), x.lab, "Selected genotypes")
    y.lab <- ifelse(!missing(y.lab), y.lab, "Proportion")
    if(genotypes == "selected"){
      data <-
        x$contri_fac %>%
        subset(GEN %in% x$sel_gen) %>%
        droplevels()
    } else{
      data <- x$contri_fac
    }
    data %<>% pivot_longer(-GEN)
    if(radar == TRUE){
      p <-
        ggplot(data, aes(x = GEN, y = value)) +
        geom_polygon(aes(group = name, color = name),
                     fill = NA,
                     size = size.line) +
        geom_polygon(aes(group = 1, x = GEN, y = 100 / length(unique(name))),
                     fill = NA,
                     color = "black",
                     linetype = 2,
                     size = size.line,
                     show.legend = FALSE) +
        geom_line(aes(group = name, color = name), size = size.line) +
        theme_minimal()  +
        theme(strip.text.x = element_text(size = size.text),
              axis.text.x = element_text(color = "black", size = size.text),
              axis.ticks.y = element_blank(),
              panel.grid = element_line(size = size.line / 2),
              axis.text.y = element_text(size = size.text, color = "black"),
              legend.position = legend.position,
              legend.title = element_blank(),
              ...) +
        labs(x = NULL,
             y = "Contribution of each factor to the MTSI index") +
        {if(title)ggtitle("The strengths and weaknesses for genotypes")} +
        scale_y_reverse() +
        guides(color = guide_legend(nrow = 1)) +
        coord_radar()
      if(arrange.label == TRUE){
        tot_gen <- length(unique(data$GEN))
        fseq <- c(1:(tot_gen/2))
        sseq <- c((tot_gen/2 + 1):tot_gen)
        fang <- c(90 - 180/length(fseq) * fseq)
        sang <- c(-90 - 180/length(sseq) * sseq)
        p <- p  +
          theme(axis.text.x = suppressMessages(suppressWarnings(element_text(angle = c(fang, sang)))), ...)
      }
    } else{
      p <-
        ggplot(data, aes(GEN, value, fill = name))+
        geom_bar(stat = "identity",
                 position = position,
                 color = "black",
                 size = size.line,
                 width = width.bar) +
        scale_y_continuous(expand = expansion(0))+
        theme_metan()  +
        theme(legend.position = legend.position,
              axis.ticks = element_line(size = size.line),
              plot.margin = margin(0.5, 0.5, 0, 0, "cm"),
              panel.border = element_rect(size = size.line),
              ...)+
        scale_x_discrete(guide = guide_axis(n.dodge = n.dodge, check.overlap = check.overlap),
                         expand = expansion(0))+
        labs(x = x.lab,
             y = y.lab) +
        {if(title)ggtitle("The strengths and weaknesses for genotypes")} +
        guides(guide_legend(nrow = 1))
      if(invert == TRUE){
        p <- p + coord_flip()
      }
    }
  }
  return(p)
}



#' Print an object of class mtmps
#'
#' Print a `mtmps` object in two ways. By default, the results are shown in
#' the R console. The results can also be exported to the directory.
#'
#' @param x An object of class `mtmps`.
#' @param export A logical argument. If `TRUE|T`, a *.txt file is exported
#'   to the working directory
#' @param file.name The name of the file if `export = TRUE`
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   [`tibble::print()`][tibble::formatting] for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print mtmps
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' model <-
#' mps(data_ge,
#'     env = ENV,
#'     gen = GEN,
#'     rep = REP,
#'     resp = everything())
#' selection <- mtmps(model)
#' print(selection)
#' }
print.mtmps <- function(x,
                        export = FALSE,
                        file.name = NULL,
                        digits = 4, ...) {
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "mtmps print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  cat("-------------------- Correlation matrix used used in factor analysis -----------------\n")
  print(x$cormat)
  cat("\n")
  cat("---------------------------- Principal component analysis -----------------------------\n")
  print(x$PCA)
  cat("\n")
  cat("--------------------------------- Initial loadings -----------------------------------\n")
  print(x$initial_loadings)
  cat("\n")
  cat("-------------------------- Loadings after varimax rotation ---------------------------\n")
  print(x$finish_loadings)
  cat("\n")
  cat("--------------------------- Scores for genotypes-ideotype -----------------------------\n")
  print(rbind(x$scores_gen, x$scores_ide))
  cat("\n")
  cat("---------------------------- Multitrait stability index ------------------------------\n")
  print(x$MTSI)
  cat("\n")
  cat("------------------------- Selection differential (variables) --------------------------\n")
  print(x$sel_dif_trait)
  cat("\n")
  cat("-------------------------------- Selected genotypes -----------------------------------\n")
  cat(x$sel_gen)
  cat("\n")
  if (export == TRUE) {
    sink()
  }
}
