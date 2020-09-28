#' Multi-trait stability index
#'
#' Computes the multi-trait stability index proposed by Olivoto et al. (2019)
#'
#'
#' @param .data An object of class \code{waasb} or \code{waas}.
#' @param index If \code{index = 'waasby'} (default) both stability and mean
#'   performance are considered. If \code{index = 'waasb'} the multi-trait index
#'   will be computed considering the stability of genotypes only.  More details
#'   can be seen in \code{\link{waasb}} and \code{\link{waas}} functions.
#' @param SI An integer (0-100). The selection intensity in percentage of the
#' total number of genotypes.
#' @param mineval The minimum value so that an eigenvector is retained in the
#' factor analysis.
#' @param verbose If \code{verbose = TRUE} (Default) then some results are
#' shown in the console.
#' @return An object of class \code{mtsi} with the following items:
#' * \strong{data} The data used to compute the factor analysis.
#' * \strong{cormat} The correlation matrix among the environments.
#' * \strong{PCA} The eigenvalues and explained variance.
#' * \strong{FA} The factor analysis.
#' * \strong{KMO} The result for the Kaiser-Meyer-Olkin test.
#' * \strong{MSA} The measure of sampling adequacy for individual variable.
#' * \strong{communalities} The communalities.
#' * \strong{communalities_mean} The communalities' mean.
#' * \strong{initial_loadings} The initial loadings.
#' * \strong{finish_loadings} The final loadings after varimax rotation.
#' * \strong{canonical_loadings} The canonical loadings.
#' * \strong{scores_gen} The scores for genotypes in all retained factors.
#' * \strong{scores_ide} The scores for the ideotype in all retained factors.
#' * \strong{MTSI} The multi-trait stability index.
#' * \strong{contri_fac} The relative contribution of each factor on the MTSI value.
#' The lower the contribution of a factor, the close of the ideotype the variables in such
#' factor are.
#' * \strong{contri_fac_rank, contri_fac_rank_sel} The rank for the contribution
#' of each factor for all genotypes and selected genotypes, respectively.
#' * \strong{sel_dif} The selection differential for the WAASBY index.
#' * \strong{mean_sd} The mean for the differential selection.
#' * \strong{sel_dif_var} The selection differential for the variables.
#' * \strong{total_sel_dif} The total selection differential.
#' * \strong{sel_dif_waasb} The selection differential for WAASB index.
#' * \strong{sel_gen} The selected genotypes.
#' @md
#' @importFrom purrr map_dfc
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @references Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, B.G. Sari, and
#'   M.I. Diel. 2019. Mean performance and stability in multi-environment trials
#'   II: Selection based on multiple traits. Agron. J. 111:2961-2969.
#' \doi{10.2134/agronj2019.03.0220}
#' @examples
#'\donttest{
#' library(metan)
#' # Based on stability only, for both GY and HM, higher is better
#' mtsi_model <- waasb(data_ge,
#'                     env = ENV,
#'                     gen = GEN,
#'                     rep = REP,
#'                     resp = c(GY, HM))
#' mtsi_index <- mtsi(mtsi_model, index = 'waasb')
#'
#'
#' # Based on mean performance and stability (using pipe operator %>%)
#' # GY: higher is better
#' # HM: lower is better
#'
#'mtsi_index2 <- data_ge %>%
#'  waasb(ENV, GEN, REP,
#'        resp = c(GY, HM),
#'        mresp = c(100, 0)) %>%
#'  mtsi()
#'}
mtsi <- function(.data,
                 index = "waasby",
                 SI = 15,
                 mineval = 1,
                 verbose = TRUE) {
  if (!index %in% c("waasb", "waasby")) {
    stop("The argument 'index' must be of of the 'waasb' or 'waasby'.", call. = FALSE)
  }
  if (length(.data) == 1) {
    stop("The multi-trait stability index cannot be computed with one single variable.", call. = FALSE)
  }
  if (index == "waasby") {
    ideotype.D <- rep(100, length(.data))
  }

  if (has_class(.data, c("waas", "waas_means"))){
    if (index == "waasb") {
      data <- gmd(.data, "WAAS", verbose = FALSE) %>% as.data.frame()
    }
    if (index == "waasby") {
      data <- gmd(.data, "WAASY", verbose = FALSE) %>% as.data.frame()
    }
  }
  if (class(.data) == "waasb") {
    if (index == "waasb") {
      data <- gmd(.data, "WAASB", verbose = FALSE) %>% as.data.frame()
    }
    if (index == "waasby") {
      data <- gmd(.data, "WAASBY", verbose = FALSE) %>% as.data.frame()
    }
  }
  if(has_na(data)){
    stop("Missing values for traits ")
  }
  if (is.null(SI)) {
    ngs <- NULL
  } else {
    ngs <- round(nrow(data) * (SI/100), 0)
  }
  observed <-
    gmd(.data, "Y", verbose = FALSE) %>%
    column_to_rownames("GEN")
  means <- data[, 2:ncol(data)]
  rownames(means) <- data[, 1]
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
  rownames(scores) <- data[, 1]
  pos.var.factor <- which(abs(A) == apply(abs(A), 1, max), arr.ind = TRUE)
  var.factor <- lapply(1:ncol(A), function(i) {
    rownames(pos.var.factor)[pos.var.factor[, 2] == i]
  })
  names(var.factor) <- paste("FA", 1:ncol(A), sep = "")
  names.pos.var.factor <- rownames(pos.var.factor)
  if (index == "waasb") {
    ideotype.D <- apply(means, 2, min)
  } else {
    names(ideotype.D) <- colnames(means)
  }
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
    mean_sd_ind <- apply(sel_dif[, 3:6], 2, mean)
    sel_dif_mean <-
      tibble(VAR = names(pos.var.factor[, 2]),
             Factor = paste("FA", as.numeric(pos.var.factor[, 2])),
             Xo = colMeans(observed),
             Xs = colMeans(observed[names(MTSI)[1:ngs], ]),
             SD = Xs - colMeans(observed),
             SDperc = (Xs - colMeans(observed)) / colMeans(observed) * 100) %>%
      left_join(
        gmd(.data, "details", verbose = FALSE) %>%
          pivot_longer(-Parameters) %>%
          subset(Parameters == "mresp") %>%
          remove_cols(Parameters) %>%
          set_names("VAR", "sense"),
        by = "VAR") %>%
      mutate(sense = ifelse(sense == 0, "decrease", "increase"),
             goal = case_when(
               sense == "decrease" & SDperc < 0  |  sense == "increase" & SDperc > 0 ~ 100,
               TRUE ~ 0
             ))
    if (class(.data) == "waasb") {
      h2 <- gmd(.data, "h2", verbose = FALSE)
      sel_dif_mean <-
        left_join(sel_dif_mean, h2, by = "VAR") %>%
        add_cols(SG = SD * h2,
                 SGperc = SG / Xo * 100,
                 .after = "SDperc") %>%
        reorder_cols(h2, .after  = "SDperc")
    }
    tota_gain <-
      sum_by(sel_dif_mean, sense) %>%
      select_cols(sense, one_of(c("SD", "SDperc", "SG", "SGperc")))
    what <- ifelse(has_class(.data, "WAASB"), "WAAS", "WAASB")
    waasb_index <- gmd(.data, what, verbose = FALSE)
    waasb_selected <- colMeans(subset(waasb_index, GEN %in% selected) %>% select_numeric_cols())
    sel_dif_waasb <-
      tibble(
        TRAIT = names(waasb_selected),
        Xo = colMeans(waasb_index %>% select_numeric_cols()),
        Xs = waasb_selected,
        SD = Xs - Xo,
        SDperc = (Xs - Xo) / Xo * 100)
    contri_fac_rank_sel <-
      contri_long %>%
      subset(GEN %in% selected) %>%
      ge_winners(name, GEN, value, type = "ranks", better = "l") %>%
      split_factors(ENV) %>%
      map_dfc(~.x %>% pull())
  }
  if (is.null(ngs)) {
    sel_dif <- NULL
    sel_dif_waasb <- NULL
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
      cat("Selection differential for the ", index, "index\n")
      cat("-------------------------------------------------------------------------------\n")
      print(sel_dif)
      cat("------------------------------------------------------------------------------\n")
      cat("Mean of selection differential\n")
      cat("-------------------------------------------------------------------------------\n")
      print(mean_sd_ind)
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

  return(structure(list(data = data,
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
                        sel_dif = sel_dif,
                        mean_sd = mean_sd_ind,
                        sel_dif_var = sel_dif_mean,
                        total_sel_dif = tota_gain,
                        sel_dif_waasb = sel_dif_waasb,
                        sel_gen = selected),
                   class = "mtsi"))
}







#' Plot the multi-trait stability index
#'
#' Makes a radar plot showing the multitrait stability index proposed by Olivoto
#' et al. (2019)
#'
#'
#' @param x An object of class \code{mtsi}
#' @param SI An integer [0-100]. The selection intensity in percentage of the
#'   total number of genotypes.
#' @param type The type of the plot. Defaults to \code{"index"}. Use \code{type
#'   = "contribution"} to show the contribution of each factor to the MGIDI
#'   index of the selected genotypes.
#' @param position The position adjustment when \code{type = "contribution"}.
#'   Defaults to \code{"fill"}, which shows relative proportions at each trait
#'   by stacking the bars and then standardizing each bar to have the same
#'   height. Use \code{position = "stack"} to plot the MGIDI index for each
#'   genotype.
#' @param genotypes When \code{type = "contribution"} defines the genotypes to
#'   be shown in the plot. By default (\code{genotypes = "selected"} only
#'   selected genotypes are shown. Use \code{genotypes = "all"} to plot the
#'   contribution for all genotypes.)
#' @param radar Logical argument. If true (default) a radar plot is generated
#'   after using \code{coord_polar()}.
#' @param arrange.label Logical argument. If \code{TRUE}, the labels are
#'   arranged to avoid text overlapping. This becomes useful when the number of
#'   genotypes is large, say, more than 30.
#' @param x.lab,y.lab The labels for the axes x and y, respectively. x label is
#'   set to null when a radar plot is produced.
#' @param size.point The size of the point in graphic. Defaults to 2.5.
#' @param size.line The size of the line in graphic. Defaults to 0.7.
#' @param size.text The size for the text in the plot. Defaults to 10.
#' @param width.bar The width of the bars if \code{type = "contribution"}.
#'   Defaults to 0.75.
#' @param n.dodge The number of rows that should be used to render the x labels.
#'   This is useful for displaying labels that would otherwise overlap.
#' @param check.overlap Silently remove overlapping labels, (recursively)
#'   prioritizing the first, last, and middle labels.
#' @param invert Logical argument. If \code{TRUE}, rotate the plot.
#' @param col.sel The colour for selected genotypes. Defaults to \code{"red"}.
#' @param col.nonsel The colour for nonselected genotypes. Defaults to \code{"black"}.
#' @param ... Other arguments to be passed from ggplot2::theme().
#' @return An object of class \code{gg, ggplot}.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method plot mtsi
#' @export
#' @references Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, B.G. Sari, and M.I. Diel. 2019. Mean performance and stability in multi-environment trials II: Selection based on multiple traits. Agron. J. (in press).
#' @examples
#' \donttest{
#' library(metan)
#' mtsi_model <- waasb(data_ge, ENV, GEN, REP, resp = c(GY, HM))
#' mtsi_index <- mtsi(mtsi_model)
#' plot(mtsi_index)
#'}
#'
#'
plot.mtsi <- function(x,
                      SI = 15,
                      type = "index",
                      position = "fill",
                      genotypes = "selected",
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
                      ...) {
  if (!class(x) == "mtsi") {
    stop("The object 'x' is not of class 'mtsi'")
  }
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
  cutpoint <- max(subset(data, sel == "sel_gen")$MTSI)
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
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_line(size = size.line),
          axis.text = element_text(colour = "black"),
          text = element_text(size = size.text)) +
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
        coord_polar() +
        theme(axis.text.x = element_text(angle = c(fang, sang)), legend.margin = margin(-120, 0, 0, 0), ...)
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
        geom_polygon(aes(group = name, color = name), fill = NA, size = size.line) +
        geom_line(aes(group = name, color = name), size = size.line) +
        geom_point(aes(group = name, color = name), size = size.line * 2) +
        theme_minimal() +
        theme(strip.text.x = element_text(size = size.text),
              axis.text.x = element_text(color = "black", size = size.text),
              axis.ticks.y = element_blank(),
              panel.grid = element_line(size = size.line),
              axis.text.y = element_text(size = size.text, color = "black"),
              legend.position = "bottom",
              legend.title = element_blank()) +
        labs(title = "The strengths and weaknesses for genotypes",
             x = NULL,
             y = "Contribution of each factor to the MTSI index") +
        scale_y_reverse() +
        guides(color = guide_legend(nrow = 1)) +
        coord_radar()

    } else{
    p <-
      ggplot(data, aes(GEN, value, fill = name))+
      geom_bar(stat = "identity",
               position = position,
               color = "black",
               size = size.line,
               width = width.bar) +
      scale_y_continuous(expand = expansion(0))+
      theme_metan()+
      theme(legend.position = "bottom",
            axis.ticks = element_line(size = size.line),
            plot.margin = margin(0.5, 0.5, 0, 0, "cm"),
            panel.border = element_rect(size = size.line))+
      scale_x_discrete(guide = guide_axis(n.dodge = n.dodge, check.overlap = check.overlap),
                       expand = expansion(0))+
      labs(title = "The strengths and weaknesses for genotypes",
           x = x.lab,
           y = y.lab)+
      guides(guide_legend(nrow = 1))
    if(invert == TRUE){
      p <- p + coord_flip()
    }
    }
  }
  return(p)
}





#' Print an object of class mtsi
#'
#' Print a \code{mtsi} object in two ways. By default, the results are shown in
#' the R console. The results can also be exported to the directory.
#'
#' @param x An object of class \code{mtsi}.
#' @param export A logical argument. If \code{TRUE|T}, a *.txt file is exported
#'   to the working directory
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print mtsi
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' # Based on stability only
#' MTSI_MODEL <- waasb(data_ge,
#'   resp = c(GY, HM),
#'   gen = GEN,
#'   env = ENV,
#'   rep = REP
#' )
#'
#' MTSI_index <- mtsi(MTSI_MODEL)
#' print(MTSI_index)
#' }
print.mtsi <- function(x, export = FALSE, file.name = NULL, digits = 4, ...) {
  if (!class(x) == "mtsi") {
    stop("The object must be of class 'mtsi'")
  }
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "mtsi print", file.name)
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
  cat("--------------------------- Selection differential (index) ----------------------------\n")
  print(x$sel_dif)
  cat("\n")
  cat("-------------------------- Mean of Selection differential -----------------------------\n")
  print(x$mean_sd)
  cat("\n")
  cat("------------------------- Selection differential (variables) --------------------------\n")
  print(x$sel_dif_var)
  cat("\n")
  cat("-------------------------------- Selected genotypes -----------------------------------\n")
  cat(x$sel_gen)
  cat("\n")
  if (export == TRUE) {
    sink()
  }
}
