#' Genotype-Ideotype Distance Index
#'
#' @description
#' Computes the multi-trait genotype-ideotype distance index (MGIDI). MGIDI can
#' be seen as the multi-trait stability index (Olivoto et al., 2019) computed
#' with weight for mean performance equals to 100. The MGIDI indes is computed
#' as follows:
#' \loadmathjax
#' \mjsdeqn{MGIDI_i = \sqrt{\sum\limits_{j = 1}^f(F_{ij} - {F_j})^2}}
#'
#' where \mjseqn{MGIDI_i} is the multi-trait genotype-ideotype distance index
#' for the \emph{i}th genotype; \mjseqn{F_{ij}} is the score of the \emph{i}th
#' genotype in the \emph{j}th factor (\emph{i = 1, 2, ..., g; j = 1, 2, ...,
#' f}), being \emph{g} and \emph{f} the number of genotypes and factors,
#' respectively, and \mjseqn{F_j} is the \emph{j}th score of the ideotype. The
#' genotype with the lowest MGIDI is then closer to the ideotype and therefore
#' presents desired values for all the analyzed traits.
#'
#' @param .data An object fitted with the function \code{\link{gafem}()},
#'   \code{\link{gamem}()} or a two-way table with BLUPs for genotypes in each
#'   trait (genotypes in rows and traits in columns). In the last case, row
#'   names must contain the genotypes names.
#' @param use_data Define which data to use if \code{.data} is an object of
#'   class \code{gamem}. Defaults to \code{"blup"} (the BLUPs for genotypes).
#'   Use \code{"pheno"} to use phenotypic means instead BLUPs for computing the
#'   index.
#' @param SI An integer (0-100). The selection intensity in percentage of the
#' total number of genotypes.
#' @param mineval The minimum value so that an eigenvector is retained in the
#' factor analysis.
#' @param ideotype A vector of length \code{nvar} where \code{nvar} is the
#'   number of variables used to plan the ideotype. Use \code{'h'} to indicate
#'   the traits in which higher values are desired or \code{'l'} to indicate the
#'   variables in which lower values are desired. For example, \code{ideotype =
#'   c("h, h, h, h, l")} will consider that the ideotype has higher values for
#'   the first four traits and lower values for the last trait. If \code{.data}
#'   is a model fitted with the functions \code{\link{gafem}()} or
#'   \code{\link{gamem}()}, the order of the traits will be the declared in the
#'   argument \code{resp} in those functions.
#' @param use The method for computing covariances in the presence of missing
#'   values. Defaults to \code{complete.obs}, i.e., missing values are handled
#'   by casewise deletion.
#' @param verbose If \code{verbose = TRUE} (Default) then some results are
#' shown in the console.
#' @return An object of class \code{mgidi} with the following items:
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
#' * \strong{gen_ide} The distance between the scores of each genotype with the
#' ideotype.|
#' * \strong{MGIDI} The multi-trait genotype-ideotype distance index.
#' * \strong{contri_fac} The relative contribution of each factor on the MGIDI
#' value.
#' The lower the contribution of a factor, the close of the ideotype the
#' variables in such factor are.
#' * \strong{sel_dif} The selection differential for the variables.
#' * \strong{total_gain} The selection differential for the variables.
#' * \strong{sel_gen} The selected genotypes.
#' @md
#' @references Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, B.G. Sari, and
#'   M.I. Diel. 2019. Mean performance and stability in multi-environment trials
#'   II: Selection based on multiple traits. Agron. J. 111:2961-2969.
#' \doi{10.2134/agronj2019.03.0220}
#' @importFrom tidyselect any_of all_of
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'\donttest{
#' library(metan)
#'
#' model <- gamem(data_g,
#'                gen = GEN,
#'                rep = REP,
#'                resp = c(NR, KW, CW, CL, NKE, TKW, PERK, PH))
#'
#' # Selection for increase all variables
#' mgidi_model <- mgidi(model)
#'
#'
#'
#' # plot the contribution of each factor on the MGIDI index
#' plot(mgidi_model, type = "contribution")
#'
#'}
mgidi <- function(.data,
                  use_data = "blup",
                  SI = 15,
                  mineval = 1,
                  ideotype = NULL,
                  use = "complete.obs",
                  verbose = TRUE) {
  d <- match.call()
  if(!use_data %in% c("blup", "pheno")){
    stop("Argument 'use_data = ", d["use_data"], "'", "invalid. It must be either 'blup' or 'pheno'.")
  }
  if(has_class(.data, c("gamem", "waasb"))){
    data <-
      gmd(.data, ifelse(use_data == "blup", "blupg", "data"), verbose = FALSE) %>%
      means_by(GEN) %>%
      column_to_rownames("GEN")
  } else if(has_class(.data, "gafem")){
    data <-
      gmd(.data, "Y", verbose = FALSE) %>%
      means_by(GEN) %>%
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
  colnames(scores) <- paste("FA", 1:ncol(scores), sep = "")
  rownames(scores) <- rownames(means)
  pos.var.factor <- which(abs(A) == apply(abs(A), 1, max), arr.ind = TRUE)
  var.factor <- lapply(1:ncol(A), function(i) {
    rownames(pos.var.factor)[pos.var.factor[, 2] == i]
  })
  names(var.factor) <- paste("FA", 1:ncol(A), sep = "")
  names.pos.var.factor <- rownames(pos.var.factor)
  ideotypes.matrix <- t(as.matrix(ideotype.D))/apply(means, 2, sd, na.rm = TRUE)
  rownames(ideotypes.matrix) <- "ID1"
  ideotypes.scores <- ideotypes.matrix %*% canonical_loadings
  gen_ide <- sweep(scores, 2, ideotypes.scores, "-")
  MGIDI <- sort(apply(gen_ide, 1, function(x) sqrt(sum(x^2))), decreasing = FALSE)
  contr.factor <- data.frame((sqrt(gen_ide^2)/apply(gen_ide, 1, function(x) sum(sqrt(x^2)))) * 100) %>%
    rownames_to_column("Gen") %>%
    as_tibble()
  means.factor <- means[, names.pos.var.factor]
  observed <- means[, names.pos.var.factor]
  if (!is.null(ngs)) {
    data_order <- data[colnames(observed)]
    sel_dif_mean <-
      tibble(VAR = names(pos.var.factor[, 2]),
             Factor = paste("FA", as.numeric(pos.var.factor[, 2])),
             Xo = colMeans(data_order, na.rm = TRUE),
             Xs = colMeans(data_order[names(MGIDI)[1:ngs], ], na.rm = TRUE),
             SD = Xs - colMeans(data_order, na.rm = TRUE),
             SDperc = (Xs - colMeans(data_order, na.rm = TRUE)) / colMeans(data_order, na.rm = TRUE) * 100)

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
    total_gain <-
      desc_stat(sel_dif_mean,
                by = sense,
                any_of(c("SDperc", "SGperc")),
                stats = c("min, mean, max, sum"))


  } else{
    sel_dif_mean <- NULL
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
      cat(names(MGIDI)[1:ngs])
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
                        sel_dif = sel_dif_mean,
                        total_gain = total_gain,
                        sel_gen = names(MGIDI)[1:ngs]),
                   class = "mgidi"))
}



#' Plot the multi-trait genotype-ideotype distance index
#'
#' Makes a radar plot showing the multi-trait genotype-ideotype distance index
#'
#'
#' @param x An object of class \code{mgidi}
#' @param SI An integer [0-100]. The selection intensity in percentage of the
#'   total number of genotypes.
#' @param radar Logical argument. If true (default) a radar plot is generated
#'   after using \code{coord_polar()}.
#' @param type The type of the plot. Defaults to \code{"index"}. Use \code{type
#'   = "contribution"} to show the contribution of each factor to the MGIDI
#'   index of the selected genotypes/treatments.
#' @param position The position adjustment when \code{type = "contribution"}.
#'   Defaults to \code{"fill"}, which shows relative proportions at each trait
#'   by stacking the bars and then standardizing each bar to have the same
#'   height. Use \code{position = "stack"} to plot the MGIDI index for each
#'   genotype/treatment.
#' @param rotate Logical argument. If \code{rotate = TRUE} the plot is rotated,
#'   i.e., traits in y axis and value in the x axis.
#' @param genotypes When \code{type = "contribution"} defines the genotypes to
#'   be shown in the plot. By default (\code{genotypes = "selected"} only
#'   selected genotypes are shown. Use \code{genotypes = "all"} to plot the
#'   contribution for all genotypes.)
#' @param n.dodge The number of rows that should be used to render the x labels.
#'   This is useful for displaying labels that would otherwise overlap.
#' @param check.overlap Silently remove overlapping labels, (recursively)
#'   prioritizing the first, last, and middle labels.
#' @param invert Deprecated argument as of 1.8.0. Use \code{rotate} instead.
#' @param x.lab,y.lab The labels for the axes x and y, respectively. x label is
#'   set to null when a radar plot is produced.
#' @param title The plot title when \code{type = "contribution"}.
#' @param arrange.label Logical argument. If \code{TRUE}, the labels are
#'   arranged to avoid text overlapping. This becomes useful when the number of
#'   genotypes is large, say, more than 30.
#' @param size.point The size of the point in graphic. Defaults to 2.5.
#' @param size.line The size of the line in graphic. Defaults to 0.7.
#' @param size.text The size for the text in the plot. Defaults to 10.
#' @param width.bar The width of the bars if \code{type = "contribution"}.
#'   Defaults to 0.75.
#' @param col.sel The colour for selected genotypes. Defaults to \code{"red"}.
#' @param col.nonsel The colour for nonselected genotypes. Defaults to \code{"black"}.
#' @param ... Other arguments to be passed from ggplot2::theme().
#' @return An object of class \code{gg, ggplot}.
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
                       invert = NULL,
                       x.lab = NULL,
                       y.lab = NULL,
                       title = NULL,
                       arrange.label = FALSE,
                       size.point = 2.5,
                       size.line = 0.7,
                       size.text = 10,
                       width.bar = 0.75,
                       col.sel = "red",
                       col.nonsel = "black",
                       ...) {
  if(!is.null(invert)){
    warning("Argument 'invert' is deprecated. Use 'replacement' instead.", call. = FALSE)
    rotate <- invert
  }
  if(!type %in% c("index", "contribution")){
    stop("The argument index must be one of the 'index' or 'contribution'", call. = FALSE)
  }
  if(!genotypes %in% c("selected", "all")){
    stop("The argument 'genotypes' must be one of the 'selected' or 'all'", call. = FALSE)
  }
  if(type == "index"){
  if (!class(x) == "mgidi") {
    stop("The object 'x' is not of class 'mgidi'")
  }
  x.lab <- ifelse(!missing(x.lab), x.lab, "Genotypes")
  y.lab <- ifelse(!missing(y.lab), y.lab, "Multi-trait genotype-ideotype distance index")
  data <- x$MGIDI %>% add_cols(sel = "Selected")
  data[["sel"]][(round(nrow(data) * (SI/100), 0) + 1):nrow(data)] <- "Nonselected"
  cutpoint <- max(subset(data, sel == "Selected")$MGIDI)

  p <-
    ggplot(data = data, aes(x = reorder(Genotype, -MGIDI), y = MGIDI)) +
    geom_hline(yintercept = cutpoint, col = col.sel, size = size.line) +
    geom_path(colour = "black", group = 1, size = size.line) +
    geom_point(size = size.point, aes(fill = sel), shape = 21, colour = "black", stroke  = size.point / 10) +
    scale_x_discrete() +
    scale_y_reverse() +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
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
      p <-
        p +
        theme(axis.text.x = element_text(angle = c(fang, sang)),
              legend.margin = margin(-120, 0, 0, 0), ...)
    }
  }
  } else{
    x.lab <- ifelse(!missing(x.lab), x.lab, "Selected genotypes")
    y.lab <- ifelse(!missing(y.lab), y.lab, "Proportion")
    if(genotypes == "selected"){
    data <-
      x$contri_fac %>%
      subset(Gen %in% x$sel_gen)
    data$Gen <-
      factor(data$Gen, levels = x$sel_gen)
    } else{
      data <- x$contri_fac
    }
    data %<>% pivot_longer(-Gen)
    title <- ifelse(is.null(title), "The strengths and weaknesses view of genotypes", title)
    p <-
      ggplot(data, aes(Gen, value, fill = name))+
      geom_bar(stat = "identity",
               position = position,
               color = "black",
               size = size.line,
               width = width.bar) +
        scale_y_continuous(expand = expansion(c(0, ifelse(position == "fill", 0, 0.05))))+
        theme_metan()+
        theme(legend.position = "bottom",
              axis.ticks = element_line(size = size.line),
              plot.margin = margin(0.5, 0.5, 0, 0, "cm"),
              panel.border = element_rect(size = size.line))+
        scale_x_discrete(guide = guide_axis(n.dodge = n.dodge, check.overlap = check.overlap),
                         expand = expansion(0))+
        labs(x = x.lab, y = y.lab)+
        guides(guide_legend(nrow = 1)) +
        ggtitle(title)
    if(rotate == TRUE){
      p <- p + coord_flip()
    }
  }
  return(p)
}


#' Print an object of class mgidi
#' Print a \code{mgidi} object in two ways. By default, the results are shown in
#' the R console. The results can also be exported to the directory.
#'
#' @param x An object of class \code{mgidi}.
#' @param export A logical argument. If \code{TRUE|T}, a *.txt file is exported
#'   to the working directory
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
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
  if (!class(x) == "mgidi") {
    stop("The object must be of class 'mgidi'")
  }
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
