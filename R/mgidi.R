#' Multitrait Genotype-Ideotype Distance Index
#'
#' @description
#' `r badge('stable')`
#'
#' Computes the multi-trait genotype-ideotype distance index, MGIDI, (Olivoto
#' and Nardino, 2020), used to select genotypes in plant breeding programs based
#' on multiple traits.The MGIDI index is computed as follows:
#' \loadmathjax
#' \mjsdeqn{MGIDI_i = \sqrt{\sum\limits_{j = 1}^f(F_{ij} - {F_j})^2}}
#'
#' where \mjseqn{MGIDI_i} is the multi-trait genotype-ideotype distance index
#' for the *i*th genotype; \mjseqn{F_{ij}} is the score of the *i*th genotype in
#' the *j*th factor (*i = 1, 2, ..., g; j = 1, 2, ..., f*), being *g* and *f*
#' the number of genotypes and factors, respectively, and \mjseqn{F_j} is the
#' *j*th score of the ideotype. The genotype with the lowest MGIDI is then
#' closer to the ideotype and therefore should presents desired values for all
#' the analyzed traits.
#'
#' @param .data An object fitted with the function [gafem()],
#'   [gamem()] or a two-way table with BLUPs for genotypes in each
#'   trait (genotypes in rows and traits in columns). In the last case, row
#'   names must contain the genotypes names.
#' @param use_data Define which data to use if `.data` is an object of
#'   class `gamem`. Defaults to `"blup"` (the BLUPs for genotypes).
#'   Use `"pheno"` to use phenotypic means instead BLUPs for computing the
#'   index.
#' @param SI An integer (0-100). The selection intensity in percentage of the
#' total number of genotypes.
#' @param mineval The minimum value so that an eigenvector is retained in the
#' factor analysis.
#' @param ideotype A vector of length `nvar` where `nvar` is the
#'   number of variables used to plan the ideotype. Use `'h'` to indicate
#'   the traits in which higher values are desired or `'l'` to indicate the
#'   variables in which lower values are desired. For example, `ideotype =
#'   c("h, h, h, h, l")` will consider that the ideotype has higher values for
#'   the first four traits and lower values for the last trait. If `.data`
#'   is a model fitted with the functions [gafem()] or
#'   [gamem()], the order of the traits will be the declared in the
#'   argument `resp` in those functions.
#' @param weights Optional weights to assign for each trait in the selection
#'   process. It must be a numeric vector of length equal to the number of
#'   traits in `.data`. By default (`NULL`) a numeric vector of weights equal to
#'   1 is used, i.e., all traits have the same weight in the selection process.
#'   It is suggested weights ranging from 0 to 1. The weights will then shrink
#'   the ideotype vector toward 0. This is useful, for example, to prioritize
#'   grain yield rather than a plant-related trait in the selection process.
#' @param use The method for computing covariances in the presence of missing
#'   values. Defaults to `complete.obs`, i.e., missing values are handled
#'   by casewise deletion.
#' @param verbose If `verbose = TRUE` (Default) then some results are
#' shown in the console.
#' @return An object of class `mgidi` with the following items:
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
#' * **gen_ide** The distance between the scores of each genotype with the
#' ideotype.
#' * **MGIDI** The multi-trait genotype-ideotype distance index.
#' * **contri_fac** The relative contribution of each factor on the MGIDI
#' value. The lower the contribution of a factor, the close of the ideotype the
#' variables in such factor are.
#' * **contri_fac_rank, contri_fac_rank_sel** The rank for the contribution
#' of each factor for all genotypes and selected genotypes, respectively.
#' * **sel_dif** The selection differential for the variables.
#' * **stat_gain** A descriptive statistic for the selection gains. The
#' minimum, mean, confidence interval, standard deviation, maximum, and sum of
#' selection gain values are computed. If traits have negative and positive
#' desired gains, the statistics are computed for by strata.
#' * **sel_gen** The selected genotypes.
#' @md
#' @references Olivoto, T., and Nardino, M. (2020). MGIDI: toward an effective
#'   multivariate selection in biological experiments. Bioinformatics.
#' \doi{10.1093/bioinformatics/btaa981}
#' @importFrom tidyselect any_of all_of
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'\donttest{
#' library(metan)
#'
#'# simulate a data set
#'# 10 genotypes
#'# 5 replications
#'# 4 traits
#' df <-
#'  g_simula(ngen = 10,
#'           nrep = 5,
#'           nvars = 4,
#'           gen_eff = 35,
#'           seed = c(1, 2, 3, 4))
#'
#'# run a mixed-effect model (genotype as random effect)
#' mod <-
#'   gamem(df,
#'         gen = GEN,
#'         rep = REP,
#'         resp = everything())
#'# BLUPs for genotypes
#' gmd(mod, "blupg")
#'
#'# Compute the MGIDI index
#'# Default options (all traits with positive desired gains)
#'# Equal weights for all traits
#'mgidi_ind <- mgidi(mod)
#'gmd(mgidi_ind, "MGIDI")
#'
#'# Higher weight for traits V1 and V4
#'# This will increase the probability of selecting H7 and H9
#'# 30% selection pressure
#' mgidi_ind2 <-
#'    mgidi(mod,
#'          weights = c(1, .2, .2, 1),
#'          SI = 30)
#'gmd(mgidi_ind2, "MGIDI")
#'
#'# plot the contribution of each factor on the MGIDI index
#'p1 <- plot(mgidi_ind, type = "contribution")
#'p2 <- plot(mgidi_ind2, type = "contribution")
#'p1 + p2
#'
#'# Positive desired gains for V1, V2 and V3
#'# Negative desired gains for V4
#'mgidi_ind3 <-
#'   mgidi(mod,
#'        ideotype = c("h, h, h, l"))
#'
#'}
mgidi <- function(.data,
                  use_data = "blup",
                  SI = 15,
                  mineval = 1,
                  ideotype = NULL,
                  weights = NULL,
                  use = "complete.obs",
                  verbose = TRUE) {
  if(has_class(.data, c("gamem_group", "gafem_group", "waasb_group"))){
    bind <-
      .data %>%
      mutate(data = map(data, ~.x %>%
                          mgidi(use_data = use_data,
                                SI = SI,
                                mineval = mineval,
                                ideotype = ideotype,
                                use = use,
                                verbose = verbose,
                                weights = weights)))
    return(set_class(bind, c("tbl_df",  "mgidi_group", "mgidi", "tbl",  "data.frame")))
  } else{
    d <- match.call()
    if(!use_data %in% c("blup", "pheno")){
      stop("Argument 'use_data = ", d["use_data"], "'", "invalid. It must be either 'blup' or 'pheno'.")
    }
    if(has_class(.data, c("gamem", "waasb"))){
      data <-
        gmd(.data, ifelse(use_data == "blup", "blupg", "data"), verbose = FALSE) %>%
        mean_by(GEN) %>%
        column_to_rownames("GEN")
    } else if(has_class(.data, "gafem")){
      data <-
        gmd(.data, "Y", verbose = FALSE) %>%
        mean_by(GEN) %>%
        column_to_rownames("GEN")
    } else{
      if(has_class(.data, c("data.frame", "matrix")) & !has_rownames(.data)){
        stop("object '", d[[".data"]], "' must have rownames.", call. = FALSE)
      }
      if(any(sapply(.data, function(x){is.numeric(x)})== FALSE)){
        stop("All variables in '", d[[".data"]], "' must be numeric.",call. = FALSE)
      }
      data <- .data
    }
    if (length(data) == 1) {
      stop("The multi-trait stability index cannot be computed with one single variable.", call. = FALSE)
    }
    if(is.null(ideotype)){
      rescaled <- rep(100, length(data))
      ideotype.D <- rep(100, length(data))
      names(ideotype.D) <- names(data)
    } else{
      rescaled <- unlist(strsplit(ideotype, split="\\s*(\\s|,)\\s*")) %>%
        all_lower_case()
      if(length(rescaled) != length(data)){
        stop("Ideotype must have length ", ncol(data), ", the number of columns in data")
      }
      if(!all(rescaled %in% c("h", "l", "m"))){
        stop("argument 'ideotype' must have 'h', 'l', or 'm' only", call. = FALSE)
      }
      ideotype.D <- ifelse(rescaled == "m", 50, 100)
      names(ideotype.D) <- names(data)
      rescaled <- case_when(
        rescaled == "h" ~ 100,
        rescaled == "l" ~ 0,
        TRUE ~ 100)
    }
    if (is.null(SI)) {
      ngs <- NULL
    } else {
      ngs <- round(nrow(data) * (SI/100), 0)
    }
    means <- data.frame(matrix(ncol = ncol(data), nrow = nrow(data)))
    rownames(means) <- rownames(data)
    vars <- tibble(VAR = colnames(data),
                   sense = rescaled) %>%
      mutate(sense = ifelse(sense == 0, "decrease", "increase"))
    for (i in 1:ncol(data)) {
      means[i] <- resca(values = data[i], new_max = rescaled[i], new_min = 100 - rescaled[i])
      colnames(means) <- colnames(data)
    }
    if(has_na(means)){
      warning("Missing values observed in the table of means. Using complete observations to compute the correlation matrix.", call. = FALSE)
    }
    if(is.null(weights)){
      weights <- rep(1, ncol(data))
    }
    cor.means <- cor(means, use = use)
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
    z <- scale(means, center = FALSE, scale = apply(means, 2, sd, na.rm = TRUE))
    canonical_loadings <- t(t(A) %*% solve_svd(cor.means))
    scores <- z %*% canonical_loadings
    pos.var.factor <- which(abs(A) == apply(abs(A), 1, max), arr.ind = TRUE)
    var.factor <- lapply(1:ncol(A), function(i) {
      rownames(pos.var.factor)[pos.var.factor[, 2] == i]
    })
    names(var.factor) <- paste("FA", 1:ncol(A), sep = "")
    names.pos.var.factor <- rownames(pos.var.factor)
    ideotypes.matrix <- t(as.matrix(ideotype.D))/apply(means, 2, sd, na.rm = TRUE) * weights
    rownames(ideotypes.matrix) <- "ID1"
    ideotypes.scores <- ideotypes.matrix %*% canonical_loadings
    gen_ide <- sweep(scores, 2, ideotypes.scores, "-")
    MGIDI <- apply(gen_ide, 1, function(x){sqrt(sum(x^2))}) %>% sort(decreasing = FALSE)
    contr.factor <-
      data.frame((sqrt(gen_ide^2)/apply(gen_ide, 1, function(x) sum(sqrt(x^2)))) * 100) %>%
      rownames_to_column("GEN") %>%
      as_tibble()
    means.factor <- means[, names.pos.var.factor]
    observed <- means[, names.pos.var.factor]
    contri_long <- pivot_longer(contr.factor, -GEN)
    contri_fac_rank <-
      contri_long %>%
      ge_winners(name, GEN, value, type = "ranks", better = "l") %>%
      split_factors(ENV) %>%
      map_dfc(~.x %>% pull())
    if (!is.null(ngs)) {
      selected <- names(MGIDI)[1:ngs]
      data_order <- data[colnames(observed)]
      sel_dif_mean <-
        tibble(VAR = names(pos.var.factor[, 2]),
               Factor = paste("FA", as.numeric(pos.var.factor[, 2]), sep = ""),
               Xo = colMeans(data_order, na.rm = TRUE),
               Xs = colMeans(data_order[selected, ], na.rm = TRUE),
               SD = Xs - colMeans(data_order, na.rm = TRUE),
               SDperc = (Xs - colMeans(data_order, na.rm = TRUE)) / abs(colMeans(data_order, na.rm = TRUE)) * 100)

      if(has_class(.data, c("gamem", "gafem"))){
        h2 <- gmd(.data, "h2", verbose = FALSE)
        sel_dif_mean <-
          left_join(sel_dif_mean, h2, by = "VAR") %>%
          add_cols(SG = SD * h2,
                   SGperc = SG / Xo * 100)
      }
      sel_dif_mean <-
        sel_dif_mean %>%
        left_join(vars, by = "VAR") %>%
        mutate(goal = case_when(
          sense == "decrease" & SDperc < 0  |  sense == "increase" & SDperc > 0 ~ 100,
          TRUE ~ 0
        ))
      stat_gain <-
        sel_dif_mean %>%
        group_by(sense) %>%
        summarise(across(any_of(c("SDperc", "SGperc")),
                         list(n = ~n(),
                              min = min,
                              mean = mean,
                              max = max,
                              sum = sum,
                              sd = sd))) %>%
        pivot_longer(-sense) %>%
        separate(name, into = c("variable", "stat")) %>%
        pivot_wider(names_from = stat, values_from = value)

      contri_fac_rank_sel <-
        contri_long %>%
        subset(GEN %in% selected) %>%
        ge_winners(name, GEN, value, type = "ranks", better = "l") %>%
        split_factors(ENV) %>%
        map_dfc(~.x %>% pull())
    } else{
      sel_dif_mean <- NULL
      contri_fac_rank_sel <- NULL
    }
    if (verbose) {
      cat("\n-------------------------------------------------------------------------------\n")
      cat("Principal Component Analysis\n")
      cat("-------------------------------------------------------------------------------\n")
      print(round_cols(pca))
      cat("-------------------------------------------------------------------------------\n")
      cat("Factor Analysis - factorial loadings after rotation-\n")
      cat("-------------------------------------------------------------------------------\n")
      print(round_cols(fa))
      cat("-------------------------------------------------------------------------------\n")
      cat("Comunalit Mean:", mean(Communality), "\n")
      cat("-------------------------------------------------------------------------------\n")
      if (!is.null(ngs)) {
        cat("Selection differential \n")
        cat("-------------------------------------------------------------------------------\n")
        print(sel_dif_mean)
        cat("------------------------------------------------------------------------------\n")
        cat("Selected genotypes\n")
        cat("-------------------------------------------------------------------------------\n")
        cat(selected)
        cat("\n-------------------------------------------------------------------------------\n")
      }
    }

    return(structure(list(data = rownames_to_column(data, "GEN"),
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
                          gen_ide = as_tibble(gen_ide, rownames = NA) %>% rownames_to_column("GEN"),
                          MGIDI = as_tibble(MGIDI, rownames = NA) %>% rownames_to_column("Genotype") %>% rename(MGIDI = value),
                          contri_fac = contr.factor,
                          contri_fac_rank = contri_fac_rank,
                          contri_fac_rank_sel = contri_fac_rank_sel,
                          sel_dif = sel_dif_mean,
                          stat_gain = stat_gain,
                          sel_gen = selected),
                     class = "mgidi"))
  }
}




#' Plot the multi-trait genotype-ideotype distance index
#'
#' Makes a radar plot showing the multi-trait genotype-ideotype distance index
#'
#'
#' @param x An object of class `mgidi`
#' @param SI An integer (0-100). The selection intensity in percentage of the
#'   total number of genotypes.
#' @param radar Logical argument. If true (default) a radar plot is generated
#'   after using `coord_polar()`.
#' @param type The type of the plot. Defaults to `"index"`. Use `type
#'   = "contribution"` to show the contribution of each factor to the MGIDI
#'   index of the selected genotypes/treatments.
#' @param position The position adjustment when `type = "contribution"`.
#'   Defaults to `"fill"`, which shows relative proportions at each trait
#'   by stacking the bars and then standardizing each bar to have the same
#'   height. Use `position = "stack"` to plot the MGIDI index for each
#'   genotype/treatment.
#' @param rotate Logical argument. If `rotate = TRUE` the plot is rotated,
#'   i.e., traits in y axis and value in the x axis.
#' @param genotypes When `type = "contribution"` defines the genotypes to
#'   be shown in the plot. By default (`genotypes = "selected"` only
#'   selected genotypes are shown. Use `genotypes = "all"` to plot the
#'   contribution for all genotypes.)
#' @param n.dodge The number of rows that should be used to render the x labels.
#'   This is useful for displaying labels that would otherwise overlap.
#' @param check.overlap Silently remove overlapping labels, (recursively)
#'   prioritizing the first, last, and middle labels.
#' @param x.lab,y.lab The labels for the axes x and y, respectively. x label is
#'   set to null when a radar plot is produced.
#' @param title The plot title when `type = "contribution"`.
#' @param arrange.label Logical argument. If `TRUE`, the labels are
#'   arranged to avoid text overlapping. This becomes useful when the number of
#'   genotypes is large, say, more than 30.
#' @param size.point The size of the point in graphic. Defaults to 2.5.
#' @param size.line The size of the line in graphic. Defaults to 0.7.
#' @param size.text The size for the text in the plot. Defaults to 10.
#' @param width.bar The width of the bars if `type = "contribution"`.
#'   Defaults to 0.75.
#' @param col.sel The colour for selected genotypes. Defaults to `"red"`.
#' @param col.nonsel The colour for nonselected genotypes. Defaults to `"gray"`.
#' @param legend.position The position of the legend.
#' @param ... Other arguments to be passed from  [ggplot2::theme()].
#' @return An object of class `gg, ggplot`.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method plot mgidi
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' model <- gamem(data_g,
#'                gen = GEN,
#'                rep = REP,
#'                resp = c(KW, NR, NKE, NKR))
#' mgidi_index <- mgidi(model)
#' plot(mgidi_index)
#'}
#'
#'
plot.mgidi <- function(x,
                       SI = 15,
                       radar = TRUE,
                       type = "index",
                       position = "fill",
                       rotate = FALSE,
                       genotypes = "selected",
                       n.dodge = 1,
                       check.overlap = FALSE,
                       x.lab = NULL,
                       y.lab = NULL,
                       title = NULL,
                       arrange.label = FALSE,
                       size.point = 2.5,
                       size.line = 0.7,
                       size.text = 10,
                       width.bar = 0.75,
                       col.sel = "red",
                       col.nonsel = "gray",
                       legend.position = "bottom",
                       ...) {
  if(!type %in% c("index", "contribution")){
    stop("The argument index must be one of the 'index' or 'contribution'", call. = FALSE)
  }
  if(!genotypes %in% c("selected", "all")){
    stop("The argument 'genotypes' must be one of the 'selected' or 'all'", call. = FALSE)
  }
  if(type == "index"){
    x.lab <- ifelse(!missing(x.lab), x.lab, "Genotypes")
    y.lab <- ifelse(!missing(y.lab), y.lab, "Multi-trait genotype-ideotype distance index")
    data <- x$MGIDI %>% add_cols(sel = "Selected")
    data[["sel"]][(round(nrow(data) * (SI/100), 0) + 1):nrow(data)] <- "Nonselected"
    cutpoint <- max(subset(data, sel == "Selected")$MGIDI)
    p <-
      ggplot(data = data, aes(x = reorder(Genotype, -MGIDI), y = MGIDI)) +
      geom_hline(yintercept = cutpoint, col = col.sel, size = size.line) +
      geom_path(colour = "black", group = 1, size = size.line) +
      geom_point(size = size.point,
                 aes(fill = sel),
                 shape = 21,
                 colour = "black",
                 stroke  = size.point / 10) +
      scale_x_discrete() +
      scale_y_reverse() +
      theme_minimal() +
      theme(legend.position = legend.position,
            legend.title = element_blank(),
            panel.grid = element_line(size = size.line / 2),
            panel.border = element_blank(),
            axis.text = element_text(colour = "black"),
            text = element_text(size = size.text),
            ...) +
      labs(y = y.lab,
           x = x.lab) +
      scale_fill_manual(values = c(col.nonsel, col.sel))
    if (radar == TRUE) {
      p <-
        p +
        coord_polar() +
        theme(axis.title.x = element_blank(), ...)
      if(arrange.label == TRUE){
        tot_gen <- length(unique(data$Genotype))
        fseq <- c(1:(tot_gen/2))
        sseq <- c((tot_gen/2 + 1):tot_gen)
        fang <- c(90 - 180/length(fseq) * fseq)
        sang <- c(-90 - 180/length(sseq) * sseq)
        p <- p +
          theme(axis.text.x = suppressMessages(suppressWarnings(element_text(angle = c(fang, sang)))), ...)
      }
    }
  } else{
    if(genotypes == "selected"){
      data <-
        x$contri_fac %>%
        subset(GEN %in% x$sel_gen)
      data$GEN <-
        factor(data$GEN, levels = x$sel_gen)
    } else{
      data <- x$contri_fac
    }
    data %<>%
      pivot_longer(-GEN) %>%
      arrange(GEN)
    title <- ifelse(is.null(title), "Strengths and weaknesses view", title)
    y.lab <- ifelse(!missing(y.lab), y.lab, "Contribution to the MGIDI")
    if(radar == TRUE){
      p <-
        ggplot(data, aes(x = GEN, y = value)) +
        geom_polygon(aes(group = name, color = name), fill = NA, size = size.line) +
        geom_polygon(aes(group = 1, x = GEN, y = 100 / length(unique(name))),
                     fill = NA,
                     color = "black",
                     linetype = 2,
                     size = size.line,
                     show.legend = FALSE) +
        geom_line(aes(group = name, color = name), size = size.line) +
        theme_minimal() +
        theme(strip.text.x = element_text(size = size.text),
              axis.text.x = element_text(color = "black", size = size.text),
              axis.ticks.y = element_blank(),
              panel.grid = element_line(size = size.line / 2),
              axis.text.y = element_text(size = size.text, color = "black"),
              legend.position = legend.position,
              legend.title = element_blank(),
              ...) +
        labs(title = title,
             x = NULL,
             y = y.lab) +
        scale_y_reverse() +
        guides(color = guide_legend(nrow = 1)) +
        coord_radar()
      if(arrange.label == TRUE){
        tot_gen <- length(unique(data$GEN))
        fseq <- c(1:(tot_gen/2))
        sseq <- c((tot_gen/2 + 1):tot_gen)
        fang <- c(90 - 180/length(fseq) * fseq)
        sang <- c(-90 - 180/length(sseq) * sseq)
        p <- p +
          theme(axis.text.x = suppressMessages(suppressWarnings(element_text(angle = c(fang, sang)))), ...)
      }
    } else{
      x.lab <- ifelse(!missing(x.lab), x.lab, "Selected genotypes")
      y.lab <- ifelse(!missing(y.lab), y.lab, "Proportion")
      p <-
        ggplot(data, aes(GEN, value, fill = name))+
        geom_bar(stat = "identity",
                 position = position,
                 color = "black",
                 size = size.line,
                 width = width.bar) +
        scale_y_continuous(expand = expansion(c(0, ifelse(position == "fill", 0, 0.05))))+
        theme_metan() +
        theme(legend.position = legend.position,
              axis.ticks = element_line(size = size.line),
              plot.margin = margin(0.5, 0.5, 0, 0, "cm"),
              panel.border = element_rect(size = size.line),
              ...)+
        scale_x_discrete(guide = guide_axis(n.dodge = n.dodge, check.overlap = check.overlap),
                         expand = expansion(0))+
        labs(x = x.lab, y = y.lab)+
        guides(guide_legend(nrow = 1)) +
        ggtitle(title)
      if(rotate == TRUE){
        p <- p + coord_flip()
      }
    }
  }
  return(p)
}


#' Print an object of class mgidi
#' Print a `mgidi` object in two ways. By default, the results are shown in
#' the R console. The results can also be exported to the directory.
#'
#' @param x An object of class `mgidi`.
#' @param export A logical argument. If `TRUE|T`, a *.txt file is exported
#'   to the working directory
#' @param file.name The name of the file if `export = TRUE`
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   [`tibble::print()`][tibble::formatting] for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print mgidi
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' model <- gamem(data_g,
#'                gen = GEN,
#'                rep = REP,
#'                resp = c(KW, NR, NKE, NKR))
#' mgidi_index <- mgidi(model)
#' print(mgidi_index)
#' }
print.mgidi <- function(x,
                        export = FALSE,
                        file.name = NULL,
                        digits = 4,
                        ...) {
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "mgidi print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  cat("-------------------------------------------------------------------------------\n")
  cat("Correlation matrix used used in factor analysis \n")
  cat("-------------------------------------------------------------------------------\n")
  print(x$cormat, digits = 2)
  cat("-------------------------------------------------------------------------------\n")
  cat("Principal component analysis \n")
  cat("-------------------------------------------------------------------------------\n")
  print(x$PCA)
  cat("-------------------------------------------------------------------------------\n")
  cat("Initial loadings \n")
  cat("-------------------------------------------------------------------------------\n")
  print(x$initial_loadings)
  cat("-------------------------------------------------------------------------------\n")
  cat("Loadings after varimax rotation \n")
  cat("-------------------------------------------------------------------------------\n")
  print(x$finish_loadings)
  cat("-------------------------------------------------------------------------------\n")
  cat("Scores for genotypes-ideotype \n")
  cat("-------------------------------------------------------------------------------\n")
  print(rbind(x$scores_gen, x$scores_ide))
  cat("-------------------------------------------------------------------------------\n")
  cat("Multi-trait genotype-ideotype distance index \n")
  cat("-------------------------------------------------------------------------------\n")
  print(x$MGIDI)
  cat("-------------------------------------------------------------------------------\n")
  cat("Selection differential \n")
  cat("-------------------------------------------------------------------------------\n")
  print(x$sel_dif)
  cat("-------------------------------------------------------------------------------\n")
  cat("Selected genotypes \n")
  cat("-------------------------------------------------------------------------------\n")
  cat(x$sel_gen)
  if (export == TRUE) {
    sink()
  }
}
