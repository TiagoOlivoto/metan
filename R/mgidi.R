#' Genotype-Ideotype Distance Index
#'
#' Computes the multi-trait genotype-ideotype distance index (MGIDI). MGIDI can
#' be seen as the multi-trait stability index (Olivoto et al., 2019) computed
#' with weight for mean performance equals to 100.
#'
#'
#' @param .data An object fitted with \code{\link{gamem}}, \code{\link{gamem_met}}, or a two-way
#'   table with BLUPs for genotypes in each trait (genotypes in rows and traits
#'   in columns). In the last case, row names must contain the genotypes names.
#' @param ideotype A vector of length \code{nvar} where \code{nvar} is the
#'   number of variables used to plan the ideotype. Use \code{'h'} to indicate
#'   the variables with increase in selction or \code{'l'} to indicate the
#'   variables with reduction in selection. For example, \code{ideotype = c("h,
#'   h, l, h, l")}.
#' @param SI An integer (0-100). The selection intensity in percentage of the
#' total number of genotypes.
#' @param mineval The minimum value so that an eigenvector is retained in the
#' factor analysis.
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
#' * \strong{communalities.mean} The communalities' mean.
#' * \strong{initial.loadings} The initial loadings.
#' * \strong{finish.loadings} The final loadings after varimax rotation.
#' * \strong{canonical.loadings} The canonical loadings.
#' * \strong{scores.gen} The scores for genotypes in all retained factors.
#' * \strong{scores.ide} The scores for the ideotype in all retained factors.
#' * \strong{MGIDI} The multi-trait stability index.
#' * \strong{contri.fac} The relative contribution of each factor on the MGIDI value.
#' The lower the contribution of a factor, the close of the ideotype the variables in such
#' factor are.
#' * \strong{sel.dif} The selection differential for the variables.
#' * \strong{Selected} The selected genotypes.
#' @md
#' @references Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, B.G. Sari, and
#'   M.I. Diel. 2019. Mean performance and stability in multi-environment trials
#'   II: Selection based on multiple traits. Agron. J. 111:2961-2969.
#' \href{https://acsess.onlinelibrary.wiley.com/doi/full/10.2134/agronj2019.03.0221}{doi:10.2134/agronj2019.03.0220}
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
#' # Selection for increase all variables
#' mgidi_model <- mgidi(model)
#'
#'}
mgidi <- function(.data,
                  SI = 15,
                  mineval = 1,
                  ideotype = NULL,
                  verbose = TRUE) {
  d <- match.call()
  if(has_class(.data, c("gamem", "waasb"))){
    data <- gmd(.data, "blupg", verbose = FALSE) %>%
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
  ideotype.D <- rep(100, length(data))
  names(ideotype.D) <- names(data)
  if(is.null(ideotype)){
    rescaled <- rep(100, length(data))
  } else{
    rescaled <- unlist(strsplit(ideotype, split="\\s*(\\s|,)\\s*")) %>%
    all_lower_case()
    if(length(rescaled) != length(data)){
      stop("Ideotype must have length ", ncol(data), ", the number of columns in data")
    }
    if(!all(rescaled %in% c("h", "l"))){
      stop("argument 'ideotype' must have 'h' or 'l' only", call. = FALSE)
    }
    rescaled <- ifelse(rescaled == "h", 100, 0)
  }
  if (is.null(SI)) {
    ngs <- NULL
  } else {
    ngs <- round(nrow(data) * (SI/100), 0)
  }
  means <- data.frame(matrix(ncol = ncol(data), nrow = nrow(data)))
  rownames(means) <- rownames(data)
  for (i in 1:ncol(data)) {
    means[i] <- resca(values = data[i], new_max = rescaled[i], new_min = 100 - rescaled[i])
    colnames(means) <- colnames(data)
  }
  cor.means <- cor(means)
  eigen.decomposition <- eigen(cor.means)
  eigen.values <- eigen.decomposition$values
  eigen.vectors <- eigen.decomposition$vectors
  colnames(eigen.vectors) <- paste("PC", 1:ncol(cor.means), sep = "")
  rownames(eigen.vectors) <- colnames(means)
  if (length(eigen.values[eigen.values >= mineval]) == 1) {
    eigen.values.factors <- as.vector(c(as.matrix(sqrt(eigen.values[eigen.values >= mineval]))))
    initial.loadings <- cbind(eigen.vectors[, eigen.values >= mineval] * eigen.values.factors)
    A <- initial.loadings
  } else {
    eigen.values.factors <-
      t(replicate(ncol(cor.means), c(as.matrix(sqrt(eigen.values[eigen.values >= mineval])))))
    initial.loadings <- eigen.vectors[, eigen.values >= mineval] * eigen.values.factors
    A <- varimax(initial.loadings)[[1]][]
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
  colnames(A) <- paste("FA", 1:ncol(initial.loadings), sep = "")
  pca <- tibble(PC = paste("PC", 1:ncol(means), sep = ""),
                Eigenvalues = eigen.values,
                `Variance (%)` = (eigen.values/sum(eigen.values)) * 100,
                `Cum. variance (%)` = cumsum(`Variance (%)`))
  Communality <- diag(A %*% t(A))
  Uniquenesses <- 1 - Communality
  fa <- cbind(A, Communality, Uniquenesses) %>% as_tibble(rownames = NA) %>%  rownames_to_column("VAR")
  z <- scale(means, center = FALSE, scale = apply(means, 2, sd))
  canonical.loadings <- t(t(A) %*% solve_svd(cor.means))
  scores <- z %*% canonical.loadings
  colnames(scores) <- paste("FA", 1:ncol(scores), sep = "")
  rownames(scores) <- rownames(means)
  pos.var.factor <- which(abs(A) == apply(abs(A), 1, max), arr.ind = TRUE)
  var.factor <- lapply(1:ncol(A), function(i) {
    rownames(pos.var.factor)[pos.var.factor[, 2] == i]
  })
  names(var.factor) <- paste("FA", 1:ncol(A), sep = "")
  names.pos.var.factor <- rownames(pos.var.factor)
  ideotypes.matrix <- t(as.matrix(ideotype.D))/apply(means, 2, sd)
  rownames(ideotypes.matrix) <- "ID1"
  ideotypes.scores <- ideotypes.matrix %*% canonical.loadings
  gen_ide <- sweep(scores, 2, ideotypes.scores, "-")
  MGIDI <- sort(apply(gen_ide, 1, function(x) sqrt(sum(x^2))), decreasing = FALSE)
  contr.factor <- data.frame((sqrt(gen_ide^2)/apply(gen_ide, 1, function(x) sum(sqrt(x^2)))) * 100) %>%
    rownames_to_column("Gen") %>%
    as_tibble()
  means.factor <- means[, names.pos.var.factor]
  observed <- means[, names.pos.var.factor]
  if (!is.null(ngs)) {
    data_order <- data[colnames(observed)]
    sel.dif.mean <- tibble(VAR = names(pos.var.factor[, 2]),
                           Factor = paste("FA", as.numeric(pos.var.factor[, 2])),
                           xo = colMeans(data_order),
                           Xs = colMeans(data_order[names(MGIDI)[1:ngs], ]),
                           SD = Xs - colMeans(data_order),
                           SDperc = (Xs - colMeans(data_order)) / colMeans(data_order) * 100)
  }
  if (is.null(ngs)) {
    sel.dif <- NULL
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
      cat("Selection differential \n")
      cat("-------------------------------------------------------------------------------\n")
      print(sel.dif.mean)
      cat("------------------------------------------------------------------------------\n")
      cat("Selected genotypes\n")
      cat("-------------------------------------------------------------------------------\n")
      cat(names(MGIDI)[1:ngs])
      cat("\n-------------------------------------------------------------------------------\n")
    }
  }
  return(structure(list(data = as_tibble(data),
                        cormat = as.matrix(cor.means),
                        PCA = pca,
                        FA = fa,
                        KMO = KMO,
                        MSA = MSA,
                        communalities = Communality,
                        communalities.mean = mean(Communality),
                        initial.loadings = data.frame(initial.loadings) %>% rownames_to_column("VAR") %>% as_tibble(),
                        finish.loadings = data.frame(A) %>% rownames_to_column("VAR") %>% as_tibble(),
                        canonical.loadings = data.frame(canonical.loadings) %>% rownames_to_column("VAR") %>% as_tibble(),
                        scores.gen = data.frame(scores) %>% rownames_to_column("GEN") %>% as_tibble(),
                        scores.ide = data.frame(ideotypes.scores) %>% rownames_to_column("GEN") %>% as_tibble(),
                        MGIDI = as_tibble(MGIDI, rownames = NA) %>% rownames_to_column("Genotype") %>% rename(MGIDI = value),
                        contri.fac = contr.factor,
                        sel.dif = sel.dif.mean,
                        Selected = names(MGIDI)[1:ngs]),
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
#' @param arrange.label Logical argument. If \code{TRUE}, the labels are
#'   arranged to avoid text overlapping. This becomes useful when the number of
#'   genotypes is large, say, more than 30.
#' @param size.point The size of the point in graphic. Defaults to 2.5.
#' @param col.sel The colour for selected genotypes.
#' @param col.nonsel The colour for nonselected genotypes.
#' @param size.text The size for the text in the plot. Defaults to 10.
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
plot.mgidi <- function(x, SI = 15, radar = TRUE, arrange.label = FALSE, size.point = 2.5,
                      col.sel = "red", col.nonsel = "black", size.text = 10, ...) {
  if (!class(x) == "mgidi") {
    stop("The object 'x' is not of class 'mgidi'")
  }
  data <- x$MGIDI %>%
    add_cols(sel = "Selected")
  data[["sel"]][(round(nrow(data) * (SI/100), 0) + 1):nrow(data)] <- "Nonselected"
  cutpoint <- max(subset(data, sel == "Selected")$MGIDI)
  p <- ggplot(data = data, aes(x = reorder(Genotype, -MGIDI),
                               y = MGIDI)) + geom_hline(yintercept = cutpoint, col = col.sel) +
    geom_path(colour = "black", group = 1) + geom_point(size = size.point,
                                                        aes(fill = sel), shape = 21, colour = "black") + scale_x_discrete() +
    scale_y_reverse() + theme_minimal() + theme(legend.position = "bottom",
                                                legend.title = element_blank(), axis.title.x = element_blank(),
                                                panel.border = element_blank(), axis.text = element_text(colour = "black"),
                                                text = element_text(size = size.text)) + labs(y = "Multi-trait genotype-ideotype distance index") +
    scale_fill_manual(values = c(col.nonsel, col.sel))
  if (radar == TRUE) {
    if(arrange.label == TRUE){
      tot_gen <- length(unique(data$Genotype))
      fseq <- c(1:(tot_gen/2))
      sseq <- c((tot_gen/2 + 1):tot_gen)
      fang <- c(90 - 180/length(fseq) * fseq)
      sang <- c(-90 - 180/length(sseq) * sseq)
      p <- p + coord_polar() + theme(axis.text.x = element_text(angle = c(fang,
                                                                          sang)), legend.margin = margin(-120, 0, 0, 0), ...)
    } else{
      p <- p + coord_polar()
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
print.mgidi <- function(x, export = FALSE, file.name = NULL, digits = 4, ...) {
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
  print(x$cormat)
  cat("-------------------------------------------------------------------------------\n")
  cat("Principal component analysis \n")
  cat("-------------------------------------------------------------------------------\n")
  print(x$PCA)
  cat("-------------------------------------------------------------------------------\n")
  cat("Initial loadings \n")
  cat("-------------------------------------------------------------------------------\n")
  print(x$initial.loadings)
  cat("-------------------------------------------------------------------------------\n")
  cat("Loadings after varimax rotation \n")
  cat("-------------------------------------------------------------------------------\n")
  print(x$finish.loadings)
  cat("-------------------------------------------------------------------------------\n")
  cat("Scores for genotypes-ideotype \n")
  cat("-------------------------------------------------------------------------------\n")
  print(rbind(x$scores.gen, x$scores.ide))
  cat("-------------------------------------------------------------------------------\n")
  cat("Multi-trait genotype-ideotype distance index \n")
  cat("-------------------------------------------------------------------------------\n")
  print(x$MGIDI)
  cat("-------------------------------------------------------------------------------\n")
  cat("Selection differential \n")
  cat("-------------------------------------------------------------------------------\n")
  print(x$sel.dif)
  cat("-------------------------------------------------------------------------------\n")
  cat("Selected genotypes \n")
  cat("-------------------------------------------------------------------------------\n")
  cat(x$Selected)
  if (export == TRUE) {
    sink()
  }
}
