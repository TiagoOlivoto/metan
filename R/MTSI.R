#' Multi-trait stability index
#'
#' Computes the multitrait stability index proposed by Olivoto et al. (2019)
#'
#'
#' @param .data An object of class \code{waasb} or \code{waas}.
#' @param index If \code{index = "waasb"} (Default) the multitrait index will
#' be computed considering the stability of genotypes only. If \code{index =
#' "waasby"} both stability and mean performance are considered. More details
#' can be seen in \code{\link{waasb}} and \code{\link{waas}} functions.
#' @param SI An integer [0-100]. The selection intensity in percentage of the
#' total number of genotypes.
#' @param mineval The minimum value so that an eigenvector is retained in the
#' factor analysis.
#' @param verbose If \code{verbose = TRUE} (Default) then some results are
#' shown in the console.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @references Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, B.G. Sari, and M.I. Diel. 2019. Mean performance and stability in multi-environment trials II: Selection based on multiple traits. Agron. J.doi:10.2134/agronj2019.03.0221.
#' @examples
#'
#' library(metan)
#' library(dplyr)
#'
#' # Based on stability only, for both GY and HM, higher is better
#' mtsi_model = waasb(data_ge,
#'                    env = ENV,
#'                    gen = GEN,
#'                    rep = REP,
#'                    resp = c(GY, HM))
#' mtsi_index = mtsi(mtsi_model)
#'
#'
#' # Based on mean performance and stability (using pipe operator %>%)
#' # GY: higher is better
#' # HM: lower is better
#' mtsi_index2 = data_ge %>%
#'               waasb(ENV, GEN, REP, c(GY, HM), mresp = c(100, 0)) %>%
#'               mtsi(index = "waasby")
#'
mtsi <- function(.data, index = "waasb", SI = 15, mineval = 1, verbose = TRUE) {
  if(!index %in% c("waasb", "waasby")){
    stop("The argument 'index' must be of of the 'waasb' or 'waasby'.")
  }
    if (length(.data) == 1) {
        stop("The multitrait stability index cannot be computed with one single variable.")
    }
    if (index == "waasby") {
        ideotype.D <- rep(100, length(.data))
    }

    if (class(.data) == "waas") {
        if (index == "waasb") {
            bind <- data.frame(do.call(cbind, lapply(.data, function(x) {
                val <- x[["model"]][["WAAS"]]
            })))
        }
        if (index == "waasby") {
            bind <- data.frame(do.call(cbind, lapply(.data, function(x) {
                val <- x[["model"]][["WAASY"]]
            })))
        }
        bind$gen <- .data[[1]][["model"]][["Code"]]
        bind$type <- .data[[1]][["model"]][["type"]]
        data <- data.frame(subset(bind, type == "GEN") %>% select(-type) %>% select(gen,
            everything()))
    }
    if (class(.data) == "waasb") {
        if (index == "waasb") {
            bind <- data.frame(do.call(cbind, lapply(.data, function(x) {
                val <- x[["model"]][["WAASB"]]
            })))
        }
        if (index == "waasby") {
            bind <- data.frame(do.call(cbind, lapply(.data, function(x) {
                val <- x[["model"]][["WAASBY"]]
            })))
        }
        bind$gen <- .data[[1]][["model"]][["Code"]]
        bind$type <- .data[[1]][["model"]][["type"]]
        data <- data.frame(subset(bind, type == "GEN") %>% select(-type) %>% select(gen,
            everything()))
    }
    if (is.null(SI)) {
        ngs <- NULL
    } else {
        ngs <- round(nrow(data) * (SI/100), 0)
    }
  observed <- data.frame(do.call(cbind, lapply(.data, function(x) {
    val <- x[["model"]][["Y"]]
  })))
  observed$gen <- .data[[1]][["model"]][["Code"]]
  observed$type <- .data[[1]][["model"]][["type"]]
  observed %<>% dplyr::filter(type == "GEN") %>%
    select(-type) %>%
    column_to_rownames("gen")

    means <- data[, 2:ncol(data)]
    rownames(means) <- data[, 1]
    cor.means <- cor(means)
    eigen.decomposition <- eigen(cor.means)
    eigen.values <- eigen.decomposition$values
    eigen.vectors <- eigen.decomposition$vectors
    colnames(eigen.vectors) <- paste("PC", 1:ncol(cor.means), sep = "")
    rownames(eigen.vectors) <- colnames(means)
    if (length(eigen.values[eigen.values >= mineval]) == 1) {
        eigen.values.factors <- as.vector(c(as.matrix(sqrt(eigen.values[eigen.values >=
            mineval]))))
        initial.loadings <- cbind(eigen.vectors[, eigen.values >= mineval] * eigen.values.factors)
        A <- initial.loadings
    } else {
        eigen.values.factors <- t(replicate(ncol(cor.means), c(as.matrix(sqrt(eigen.values[eigen.values >=
            mineval])))))
        initial.loadings <- eigen.vectors[, eigen.values >= mineval] * eigen.values.factors
        A <- varimax(initial.loadings)[[1]][]
    }
    partial <- solve(cor.means)
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
    variance <- (eigen.values/sum(eigen.values)) * 100
    cumulative.var <- cumsum(eigen.values/sum(eigen.values)) * 100
    pca <- cbind(eigen.values, variance, cumulative.var)
    colnames(pca) <- c("Eigenvalues", "Variance (%)", "Cum. variance (%)")
    rownames(pca) <- paste("PC", 1:ncol(means), sep = "")
    Communality <- diag(A %*% t(A))
    Uniquenesses <- 1 - Communality
    fa <- cbind(A, Communality, Uniquenesses)
    z <- scale(means, center = F, scale = apply(means, 2, sd))
    canonical.loadings <- t(t(A) %*% solve(cor.means))
    scores <- z %*% canonical.loadings
    colnames(scores) <- paste("FA", 1:ncol(scores), sep = "")
    rownames(scores) <- data[, 1]
    pos.var.factor <- which(abs(A) == apply(abs(A), 1, max), arr.ind = T)
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
    ideotypes.scores <- ideotypes.matrix %*% canonical.loadings
    gen_ide <- sweep(scores, 2, ideotypes.scores, "-")
    MTSI <- sort(apply(gen_ide, 1, function(x) sqrt(sum(x^2))), decreasing = F)
    contr.factor <- (gen_ide^2/apply(gen_ide, 1, function(x) sum(x^2))) * 100
    means.factor <- means[, names.pos.var.factor]
    observed <- observed[, names.pos.var.factor]

    if (!is.null(ngs)) {
        selection.diferential <- data.frame(cbind(Factor = pos.var.factor[, 2],
                                                  Xo = colMeans(means.factor),
                                                  Xs = colMeans(means.factor[names(MTSI)[1:ngs], ]),
                                                  SD = colMeans(means.factor[names(MTSI)[1:ngs], ]) - colMeans(means.factor),
                                                  SDperc = (colMeans(means.factor[names(MTSI)[1:ngs], ]) - colMeans(means.factor))/colMeans(means.factor) * 100))
        selection.diferential[, 1] <- paste("FA", selection.diferential[, 1], sep = "")
        mean_sd_ind <- apply(selection.diferential[, 2:5], 2, mean)
        sel.dif.mean <- data.frame(cbind(Factor = pos.var.factor[, 2],
                                         Xo = colMeans(observed),
                                         Xs = colMeans(observed[names(MTSI)[1:ngs], ]),
                                         SD = colMeans(observed[names(MTSI)[1:ngs], ]) - colMeans(observed),
                                         SDperc = (colMeans(observed[names(MTSI)[1:ngs], ]) - colMeans(observed))/colMeans(observed) * 100))
        sel.dif.mean[, 1] <- paste("FA", sel.dif.mean[, 1], sep = "")

    }
    if (is.null(ngs)) {
        selection.diferential <- NULL
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
        cat("Multitrait stability index\n")
        cat("-------------------------------------------------------------------------------\n")
        print(round(MTSI, 4))
        cat("-------------------------------------------------------------------------------\n")
        if (!is.null(ngs)) {
            cat("Selection differential for the ", index, "index\n")
            cat("-------------------------------------------------------------------------------\n")
            print(selection.diferential)
            cat("------------------------------------------------------------------------------\n")
            cat("Mean of selection differential\n")
            cat("-------------------------------------------------------------------------------\n")
            print(mean_sd_ind)
            cat("-------------------------------------------------------------------------------\n")
            cat("Selection differential for the mean of the variables\n")
            cat("-------------------------------------------------------------------------------\n")
            print(sel.dif.mean)
            cat("------------------------------------------------------------------------------\n")
            cat("Selected genotypes\n")
            cat("-------------------------------------------------------------------------------\n")
            cat(names(MTSI)[1:ngs])
            cat("\n-------------------------------------------------------------------------------\n")
        }
    }
    return(
      structure(
        list(data = data,
             cormat = as.matrix(cor.means),
             PCA = data.frame(pca),
             FA = data.frame(fa),
             KMO = KMO,
             MSA = MSA,
             comunalits = Communality,
             comunalits.mean = mean(Communality),
             initial.loadings = initial.loadings,
             finish.loadings = A,
             canonical.loadings = canonical.loadings,
             scores.gen = scores,
             scores.ide = ideotypes.scores,
             MTSI = MTSI,
             contri.fac = data.frame(contr.factor),
             selection.diferential = selection.diferential,
             mean.sd = mean_sd_ind,
             sel.dif.obs = sel.dif.mean,
             Selected = names(MTSI)[1:ngs]),
        class = "mtsi"))
}
