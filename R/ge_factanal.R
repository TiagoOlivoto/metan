#' Stability analysis and environment stratification
#'
#' This function computes the stability analysis and environmental stratification
#' using factor analysis as proposed by Murakami and Cruz (2004).
#'
#' @param .data The dataset containing the columns related to Environments, Genotypes,
#'              replication/block and response variable(s)
#' @param env The name of the column that contains the levels of the
#' environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#' replications/blocks
#' @param resp The response variable(s). To analyze multiple variables in a
#' single procedure use, for example, \code{resp = c(var1, var2, var3)}.
#' @param mineval The minimum value so that an eigenvector is retained in the
#' factor analysis.
#' @param verbose Logical argument. If \code{verbose = FALSE} the code will run silently.
#' @return An object of class \code{ge_factanal} with the following items:
#' \item{data}{The data used to compute the factor analysis.}
#' \item{cormat}{The correlation matrix among the environments.}
#' \item{PCA}{The eigenvalues and explained variance.}
#' \item{FA}{The factor analysis.}
#' \item{env_strat}{The environmental stratification.}
#' \item{KMO}{The result for the Kaiser-Meyer-Olkin test.}
#' \item{MSA}{The measure of sampling adequacy for individual variable.}
#' \item{communalities}{The communalities.}
#' \item{communalities.mean}{The communalities' mean.}
#' \item{initial.loadings}{The initial loadings.}
#' \item{finish.loadings}{The final loadings after varimax rotation.}
#' \item{canonical.loadings}{The canonical loadings.}
#' \item{scores.gen}{The scores for genotypes for the first and second factors.}
#' @references Murakami, D.M.D., and C.D.C. Cruz. 2004. Proposal of methodologies for
#' environment stratification and analysis of genotype adaptability.
#' Crop Breed. Appl. Biotechnol. 4:7-11.
#'
#' @author Tiago Olivoto, \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' model <- ge_factanal(data_ge2,
#'                      env = ENV,
#'                      gen = GEN,
#'                      rep = REP,
#'                      resp = PH)
#'}
#' @seealso \code{\link{superiority}, \link{ecovalence}, \link{ge_stats}, \link{ge_reg}}
#'
#'
ge_factanal <- function(.data, env, gen, rep, resp, mineval = 1,
                        verbose = TRUE) {
    factors  <-
        .data %>%
        select({{env}}, {{gen}}, {{rep}}) %>%
        mutate_all(as.factor)
    vars <- .data %>% select({{resp}}, -names(factors))
    vars %<>% select_numeric_cols()
    factors %<>% set_names("ENV", "GEN", "REP")
    listres <- list()
    nvar <- ncol(vars)
    for (var in 1:nvar) {
        data <- factors %>%
            mutate(mean = vars[[var]])
        if(has_na(data)){
            data <- remove_rows_na(data)
            has_text_in_num(data)
        }
        means <- make_mat(data, GEN, ENV, mean)
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
        variance <- (eigen.values/sum(eigen.values)) * 100
        cumulative.var <- cumsum(eigen.values/sum(eigen.values)) *
            100
        pca <- data.frame(PCA = paste("PC", 1:ncol(means), sep = ""),
                          Eigenvalues = eigen.values, Variance = variance,
                          Cumul_var = cumulative.var)
        Communality <- diag(A %*% t(A))
        Uniquenesses <- 1 - Communality
        fa <- data.frame(Env = names(means), A, Communality,
                         Uniquenesses)
        z <- scale(means, center = FALSE, scale = apply(means, 2,
                                                        sd))
        canonical.loadings <- t(t(A) %*% solve_svd(cor.means))
        scores <- z %*% canonical.loadings
        colnames(scores) <- paste("FA", 1:ncol(scores), sep = "")
        rownames(scores) <- rownames(means)
        pos.var.factor <- which(abs(A) == apply(abs(A), 1, max),
                                arr.ind = TRUE)
        var.factor <- lapply(1:ncol(A), function(i) {
            rownames(pos.var.factor)[pos.var.factor[, 2] == i]
        })
        names(var.factor) <- paste("FA", 1:ncol(A), sep = "")
        names.pos.var.factor <- rownames(pos.var.factor)
        means.factor <- means[, names.pos.var.factor]
        genv <- data.frame(Env = names(means.factor),
                           Factor = paste("FA", pos.var.factor[, 2], sep = ""),
                           Mean = colMeans(means.factor),
                           Min = apply(means.factor, 2, min),
                           Max = apply(means.factor, 2, max),
                           CV = (apply(means.factor, 2, sd)/apply(means.factor, 2, mean)) * 100)
        colnames(initial.loadings) <- paste("FA", 1:ncol(initial.loadings), sep = "")
        if(ncol(scores) < 2){
            warning("The number of retained factors is ",ncol(scores),
                    ".\nA plot with the scores cannot be obtained.\nUse 'mineval' to increase the number of factors retained", call. = FALSE)
        }
        temp <- (structure(list(data = as_tibble(data),
                                cormat = as.matrix(cor.means),
                                PCA = as_tibble(pca),
                                FA = as_tibble(fa),
                                env_strat = as_tibble(genv),
                                KMO = KMO,
                                MSA = MSA,
                                communalities = Communality,
                                communalities.mean = mean(Communality),
                                initial.loadings = as_tibble(cbind(Env = names(means), as_tibble(initial.loadings))),
                                finish.loadings = as_tibble(cbind(Env = names(means), as_tibble(A))),
                                canonical.loadings = as_tibble(cbind(Env = names(means), as_tibble(canonical.loadings))),
                                scores.gen = as_tibble(cbind(Gen = rownames(means), as_tibble(scores)))), class = "ge_factanal"))
        listres[[paste(names(vars[var]))]] <- temp
    }
    return(structure(listres, class = "ge_factanal"))
}







#' Plot the ge_factanal model
#'
#' This function plot the scores for genotypes obtained in the factor analysis
#' to interpret the stability
#'
#' @param x An object of class \code{ge_factanal}
#' @param var The variable to plot. Defaults to \code{var = 1} the first
#'   variable of \code{x}.
#' @param plot_theme The graphical theme of the plot. Default is
#'   \code{plot_theme = theme_metan()}. For more details, see
#'   \code{\link[ggplot2]{theme}}.
#' @param x.lim The range of x-axis. Default is \code{NULL} (maximum and minimum
#'   values of the data set). New arguments can be inserted as \code{x.lim =
#'   c(x.min, x.max)}.
#' @param x.breaks The breaks to be plotted in the x-axis. Default is
#'   \code{authomatic breaks}. New arguments can be inserted as \code{x.breaks =
#'   c(breaks)}
#' @param x.lab The label of x-axis. Each plot has a default value. New
#'   arguments can be inserted as \code{x.lab = "my label"}.
#' @param y.lim The range of x-axis. Default is \code{NULL}. The same arguments
#'   than \code{x.lim} can be used.
#' @param y.breaks The breaks to be plotted in the x-axis. Default is
#'   \code{authomatic breaks}. The same arguments than \code{x.breaks} can be
#'   used.
#' @param y.lab The label of y-axis. Each plot has a default value. New
#'   arguments can be inserted as \code{y.lab = "my label"}.
#' @param shape The shape for genotype indication in the plot. Default is
#'   \code{1} (circle). Values between  \code{21-25}: \code{21} (circle),
#'   \code{22} (square), \code{23} (diamond), \code{24} (up triangle), and
#'   \code{25} (low triangle) allows a color for fill the shape.
#' @param col.shape The shape color for genotypes. Must be one value or a vector
#'   of colors with the same length of the number of genotypes. Default is
#'   \code{"gray30"}. Other values can be attributed. For example,
#'   \code{transparent_color()}, will make a plot with only an outline around the
#'   shape area.
#' @param col.alpha The alpha value for the color. Default is \code{1}. Values
#'   must be between \code{0} (full transparency) to \code{1} (full color).
#' @param size.shape The size of the shape (both for genotypes and
#'   environments). Default is \code{2.2}.
#' @param size.bor.tick The size of tick of shape. Default is \code{0.3}. The
#'   size of the shape will be \code{size.shape + size.bor.tick}
#' @param size.tex.lab The size of the text in the axes text and labels. Default
#'   is \code{12}.
#' @param size.tex.pa The size of the text of the plot area. Default is
#'   \code{3.5}.
#' @param force.repel Force of repulsion between overlapping text labels.
#'   Defaults to 1.
#' @param line.type The type of the line that indicate the means in the biplot.
#'   Default is \code{"solid"}. Other values that can be attributed are:
#'   \code{"blank"}, no lines in the biplot, \code{"dashed", "dotted",
#'   "dotdash", "longdash", and "twodash"}.
#' @param line.alpha The alpha value that combine the line with the background
#'   to create the appearance of partial or full transparency. Default is
#'   \code{0.4}. Values must be between "0" (full transparency) to "1" (full
#'   color).
#' @param col.line The color of the line that indicate the means in the biplot.
#'   Default is \code{"gray"}
#' @param size.line The size of the line that indicate the means in the biplot.
#'   Default is \code{0.5}.
#' @param ... Currently not used..
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso \code{\link{ge_factanal}}
#' @method plot ge_factanal
#' @return An object of class \code{gg, ggplot}.
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' library(ggplot2)
#' model = ge_factanal(data_ge2,
#'                     env = ENV,
#'                     gen = GEN,
#'                     rep = REP,
#'                     resp = PH)
#' plot(model)
#'
#' plot(model,
#'      size.shape = 3,
#'      force.repel = 10,
#'      col.shape = "orange",
#'      col.line = "red")
#'}
plot.ge_factanal <- function(x, var = 1, plot_theme = theme_metan(), x.lim = NULL, x.breaks = waiver(),
                             x.lab = NULL, y.lim = NULL, y.breaks = waiver(), y.lab = NULL,
                             shape = 21, col.shape = "gray30", col.alpha = 1, size.shape = 2.2,
                             size.bor.tick = 0.3, size.tex.lab = 12, size.tex.pa = 3.5,
                             force.repel = 1, line.type = "dashed", line.alpha = 1,
                             col.line = "black", size.line = 0.5,  ...) {
    x <- x[[var]]
    if (!class(x) == "ge_factanal") {
        stop("The object 'x' is not of class 'ge_factanal'")
    }
    data <- data.frame(scores.gen)
    if(ncol(data) == 2){
        stop("A plot cannot be generated with only one factor. \nUse 'mineval' argument in 'ge_factanal()' to increase the number of factors retained.", call. = FALSE)
    }
    if (is.null(y.lab) == FALSE) {
        y.lab <- y.lab
    } else {
        y.lab <- paste("Factor 2 (",round(x$PCA$Variance[[2]],2), "%)", sep = "")
    }
    if (is.null(x.lab) == FALSE) {
        x.lab <- x.lab
    } else {
        x.lab <- paste("Factor 1 (",round(x$PCA$Variance[[1]],2), "%)", sep = "")
    }

    p <- ggplot(data = data, aes(x = FA1, y = FA2)) +
        geom_hline(yintercept = mean(data[,3]), linetype = line.type, color = col.line, size = size.line, alpha = line.alpha)+
        geom_vline(xintercept = mean(data[,2]), linetype = line.type, color = col.line, size = size.line, alpha = line.alpha)+
        geom_point(shape = shape, size = size.shape, fill = col.shape, stroke = size.bor.tick, alpha = col.alpha)+
        labs(x = x.lab, y = y.lab)+
        geom_text_repel(aes(label = Gen), size = size.tex.pa, force = force.repel)+
        scale_x_continuous(limits = x.lim, breaks = x.breaks) +
        scale_y_continuous(limits = y.lim, breaks = y.breaks) +
        plot_theme %+replace%
        theme(aspect.ratio = 1,
              axis.text = element_text(size = size.tex.lab, color = "black"),
              axis.title = element_text(size = size.tex.lab, color = "black"),
              axis.ticks = element_line(color = "black"))
    return(p)
}







#' Print an object of class ge_factanal
#'
#' Print the \code{ge_factanal} object in two ways. By default, the results are
#' shown in the R console. The results can also be exported to the directory.
#'
#'
#' @param x An object of class \code{ge_factanal}.
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported to
#'   the working directory
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print ge_factanal
#' @export
#' @examples
#' \donttest{
#' model <- ge_factanal(data_ge2,
#'   env = ENV,
#'   gen = GEN,
#'   rep = REP,
#'   resp = PH
#' )
#' print(model)
#' }
print.ge_factanal <- function(x, export = FALSE, file.name = NULL, digits = 4, ...) {
    if (!class(x) == "ge_factanal") {
        stop("The object must be of class 'ge_factanal'")
    }
    opar <- options(pillar.sigfig = digits)
    on.exit(options(opar))
    if (export == TRUE) {
        file.name <- ifelse(is.null(file.name) == TRUE, "ge_factanal print", file.name)
        sink(paste0(file.name, ".txt"))
    }
    for (i in 1:length(x)) {
        var <- x[[i]]
        cat("Variable", names(x)[i], "\n")
        cat("------------------------------------------------------------------------------------\n")
        cat("Correlation matrix among environments\n")
        cat("------------------------------------------------------------------------------------\n")
        print(as_tibble(var$cormat, rownames = "ENV"))
        cat("------------------------------------------------------------------------------------\n")
        cat("Eigenvalues and explained variance\n")
        cat("------------------------------------------------------------------------------------\n")
        print(var$PCA)
        cat("------------------------------------------------------------------------------------\n")
        cat("Initial loadings\n")
        cat("------------------------------------------------------------------------------------\n")
        print(var$initial.loadings)
        cat("------------------------------------------------------------------------------------\n")
        cat("Loadings after varimax rotation and commonalities\n")
        cat("------------------------------------------------------------------------------------\n")
        print(var$FA)
        cat("------------------------------------------------------------------------------------\n")
        cat("Environmental stratification based on factor analysis\n")
        cat("------------------------------------------------------------------------------------\n")
        print(var$env_strat)
        cat("------------------------------------------------------------------------------------\n")
        cat("Mean = mean; Min = minimum; Max = maximum; CV = coefficient of variation (%)\n")
        cat("------------------------------------------------------------------------------------\n")
        cat("\n\n\n")
    }
    if (export == TRUE) {
        sink()
    }
}
