#' Environment stratification
#' @description
#' `r badge('stable')`
#'
#' Computes environment stratification based on factor analysis.
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s)
#' @param env The name of the column that contains the levels of the
#' environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param resp The response variable(s). To analyze multiple variables in a
#' single procedure use, for example, `resp = c(var1, var2, var3)`.
#' @param use The method for computing covariances in the presence of missing
#'   values. Defaults to `complete.obs`, i.e., missing values are handled
#'   by casewise deletion.
#' @param mineval The minimum value so that an eigenvector is retained in the
#' factor analysis.
#' @param verbose Logical argument. If `verbose = FALSE` the code will run
#'   silently.
#' @return An object of class `env_stratification` which is a list with one
#'   element per analyzed trait. For each trait, the following values are given.
#' * `data` The genotype-environment means.
#' * `cormat`: The correlation matrix among the environments.
#' * `PCA`: The eigenvalues and explained variance.
#' * `FA`: The factor analysis.
#' * `env_strat`: The environmental stratification.
#' * `mega_env_code`: The environments within each mega-environment.
#' * `mega_env_stat`: The statistics for each mega-environment.
#' * `KMO`: The result for the Kaiser-Meyer-Olkin test.
#' * `MSA`: The measure of sampling adequacy for individual variable.
#' * `communalities_mean`: The communalities' mean.
#' * `initial_loadings`: The initial loadings.
#' @references Murakami, D.M.D., and C.D.C. Cruz. 2004. Proposal of
#'   methodologies for environment stratification and analysis of genotype
#'   adaptability. Crop Breed. Appl. Biotechnol. 4:7-11.
#' @md
#' @author Tiago Olivoto, \email{tiagoolivoto@@gmail.com}
#' @importFrom tidyr chop
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' model <-
#' env_stratification(data_ge,
#'                    env = ENV,
#'                    gen = GEN,
#'                    resp = everything())
#' gmd(model)
#'
#'}
#' @seealso [env_dissimilarity()]
#'
#'
env_stratification <- function(.data,
                               env,
                               gen,
                               resp,
                               use = "complete.obs",
                               mineval = 1,
                               verbose = TRUE) {
    factors  <-
        .data %>%
        select({{env}}, {{gen}}) %>%
        as_factor(everything())
    vars <- .data %>% select({{resp}}, -names(factors))
    vars %<>% select_numeric_cols()
    factors %<>% set_names("ENV", "GEN")
    listres <- list()
    nvar <- ncol(vars)
    for (var in 1:nvar) {
        data <- factors %>%
            mutate(Y = vars[[var]])
        if(has_na(data)){
            data <- remove_rows_na(data, verbose = verbose)
            has_text_in_num(data)
        }
        means <- make_mat(data, GEN, ENV, Y)
        cor.means <- cor(means, use = use)
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
            eigen.values.factors <- t(replicate(ncol(cor.means), c(as.matrix(sqrt(eigen.values[eigen.values >= mineval])))))
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
                    partial[i, j] <- -partial[i, j]/sqrt(partial[i,
                                                                 i] * partial[j, j])
                }
            }
        }
        KMO <- sum((cor.means[!diag(k)])^2)/(sum((cor.means[!diag(k)])^2) +
                                                 sum((partial[!diag(k)])^2))
        MSA <- unlist(lapply(seq_k, function(i) {
            sum((cor.means[i, -i])^2)/(sum((cor.means[i, -i])^2) +
                                           sum((partial[i, -i])^2))
        }))
        names(MSA) <- colnames(means)
        colnames(A) <- paste("FA", 1:ncol(initial.loadings), sep = "")
        pca <- tibble(PCA = paste("PC", 1:ncol(means), sep = ""),
                      Eigenvalues = eigen.values,
                      Variance = (eigen.values/sum(eigen.values)) * 100,
                      Cumul_var = cumsum(Variance))
        Communality <- diag(A %*% t(A))
        Uniquenesses <- 1 - Communality
        fa <- data.frame(Env = names(means), A, Communality, Uniquenesses)
        z <- scale(means, center = FALSE, scale = apply(means, 2, sd))
        canonical.loadings <- t(t(A) %*% solve_svd(cor.means))
        pos.var.factor <- which(abs(A) == apply(abs(A), 1, max),arr.ind = TRUE)
        var.factor <- lapply(1:ncol(A), function(i) {
            rownames(pos.var.factor)[pos.var.factor[, 2] == i]
        })
        names(var.factor) <- paste("ME", 1:ncol(A), sep = "")
        names.pos.var.factor <- rownames(pos.var.factor)
        means.factor <- means[, names.pos.var.factor]
        genv <- tibble(ENV = names(means.factor),
                       MEGA_ENV = paste("ME", pos.var.factor[, 2], sep = ""),
                       MEAN = colMeans(means.factor, na.rm = TRUE),
                       MIN = apply(means.factor, 2, min, na.rm = TRUE),
                       MAX = apply(means.factor, 2, max, na.rm = TRUE),
                       CV = (apply(means.factor, 2, sd, na.rm = TRUE)/apply(means.factor, 2, mean, na.rm = TRUE)) * 100)
        colnames(initial.loadings) <- paste("ME", 1:ncol(initial.loadings), sep = "")
        temp <- list(data = means,
                     cormat = as.matrix(cor.means),
                     PCA = as_tibble(pca) %>% colnames_to_upper(),
                     FA = as_tibble(fa) %>% colnames_to_upper(),
                     env_strat = as_tibble(genv),
                     mega_env_code = genv %>% select_cols(1:2) %>%  chop(ENV) %>% as.data.frame(),
                     mega_env_stat = genv %>% means_by(MEGA_ENV, verbose = verbose, na.rm = TRUE) %>% remove_cols(verbose),
                     KMO = KMO,
                     MSA = MSA,
                     communalities_mean = mean(Communality),
                     initial_loadings = as_tibble(cbind(Env = names(means), initial.loadings))) %>%
            add_class("env_stratification")
        listres[[paste(names(vars[var]))]] <- temp
    }
    return(structure(listres, class = "env_stratification"))
}


#' Plot the env_stratification model
#'
#' This function plots the correlation between environments generated with
#' [env_stratification()]
#' @param x An object of class `env_stratification`
#' @param var The variable to plot. Defaults to `var = 1` the first
#'   variable of `x`.
#' @param ... Further arguments passed to [plot.corr_coef()]
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso [env_dissimilarity()]
#' @method plot env_stratification
#' @return An object of class `gg, ggplot`.
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' model <-
#' env_stratification(data_ge,
#'                    env = ENV,
#'                    gen = GEN,
#'                    resp = GY)
#' plot(model)
#'}
plot.env_stratification <- function(x,
                                    var = 1,
                                    ...) {
    if (!has_class(x, "env_stratification")) {
        stop("The object 'x' is not of class 'ge_factanal'")
    }
    x <- x[[var]]

    cormat <- x[["data"]] %>% corr_coef()
    p <- plot(cormat, ...)
    return(p)
}


#' Print the env_stratification model
#'
#' Print an object of class `ge_factanal` in two ways. By default, the results are
#' shown in the R console. The results can also be exported to the directory.
#'
#'
#' @param x An object of class `env_stratification`.
#' @param export A logical argument. If `TRUE`, a *.txt file is exported to
#'   the working directory
#' @param file.name The name of the file if `export = TRUE`
#' @param digits The significant digits to be shown.
#' @param ... Currently not used.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print env_stratification
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' model <-
#' env_stratification(data_ge,
#'                    env = ENV,
#'                    gen = GEN,
#'                    resp = GY)
#' print(model)
#' }
print.env_stratification <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
    if (!class(x) == "env_stratification") {
        stop("The object must be of class 'env_stratification'")
    }
    opar <- options(pillar.sigfig = digits)
    on.exit(options(opar))
    if (export == TRUE) {
        file.name <- ifelse(is.null(file.name) == TRUE, "env_stratification print", file.name)
        sink(paste0(file.name, ".txt"))
    }
    for (i in 1:length(x)) {
        var <- x[[i]]
        cat("Variable", names(x)[i], "\n")
        cat("------------------------------------------------------------------------------------\n")
        cat("Environment stratification based on factor analysis\n")
        cat("------------------------------------------------------------------------------------\n")
        a <- var$mega_env_code
        print.data.frame(a, row.names = FALSE)
        cat("------------------------------------------------------------------------------------\n")
        cat("Statistic by environment \n")
        cat("------------------------------------------------------------------------------------\n")
        print.data.frame(var$env_strat %>% round_cols(digits = digits), row.names = FALSE)
        cat("------------------------------------------------------------------------------------\n")
        cat("Statistic by mega-environment (mean values) \n")
        cat("------------------------------------------------------------------------------------\n")
        print.data.frame(var$mega_env_stat %>% round_cols(digits = digits), row.names = FALSE)
        cat("------------------------------------------------------------------------------------\n")
        cat("Mean = mean; Min = minimum; Max = maximum; CV = coefficient of variation (%)\n")
        cat("------------------------------------------------------------------------------------\n")
        cat("\n\n")
    }
    if (export == TRUE) {
        sink()
    }
}
