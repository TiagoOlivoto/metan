#' Clustering analysis
#' @description
#' `r badge('stable')`
#' \loadmathjax
#' Performs clustering analysis with selection of variables.
#' @details
#' When `selvar = TRUE` a variable selection algorithm is executed. The
#' objective is to select a group of variables that most contribute to explain
#' the variability of the original data. The selection of the variables is based
#' on eigenvalue/eigenvectors solution based on the following steps.
#' 1. compute the distance matrix and the cophenetic correlation with the original
#' variables (all numeric variables in dataset);
#' 2. compute the eigenvalues and eigenvectors of the correlation matrix between
#' the variables;
#' 3. Delete the variable with the largest weight (highest eigenvector in
#' the lowest eigenvalue);
#' 4. Compute the distance matrix and cophenetic correlation with the remaining
#' variables;
#' 5. Compute the Mantel's correlation between the obtained distances matrix and
#' the original distance matrix;
#' 6. Iterate steps 2 to 5 *p* - 2 times, where *p* is the number of original
#' variables.
#'
#' At the end of the *p* - 2 iterations, a summary of the models is returned.
#' The distance is calculated with the variables that generated the model with
#' the largest cophenetic correlation. I suggest a careful evaluation aiming at
#' choosing a parsimonious model, i.e., the one with the fewer number of
#' variables, that presents acceptable cophenetic correlation and high
#' similarity with the original distances.
#'
#'@param .data The data to be analyzed. It can be a data frame, possible with
#'  grouped data passed from [dplyr::group_by()].
#' @param ... The variables in `.data` to compute the distances. Set to
#'   `NULL`, i.e., all the numeric variables in `.data` are used.
#'@param by One variable (factor) to compute the function by. It is a shortcut
#'  to [dplyr::group_by()]. To compute the statistics by more than
#'  one grouping variable use that function.
#' @param scale Should the data be scaled before computing the distances? Set to
#'   FALSE. If TRUE, then, each observation will be divided by the standard
#'   deviation of the variable \mjseqn{Z_{ij} = X_{ij} / sd_j}
#' @param selvar Logical argument, set to `FALSE`. If `TRUE`, then an
#'   algorithm for selecting variables is implemented. See the section
#'   **Details** for additional information.
#' @param verbose Logical argument. If `TRUE` (default) then the results
#'   for variable selection are shown in the console.
#' @param distmethod The distance measure to be used. This must be one of
#'   `'euclidean'`, `'maximum'`, `'manhattan'`,
#'   `'canberra'`, `'binary'`, `'minkowski'`, `'pearson'`,
#'   `'spearman'`, or `'kendall'`. The last three are
#'   correlation-based distance.
#' @param clustmethod The agglomeration method to be used. This should be one of
#'   `'ward.D'`, `'ward.D2'`, `'single'`, `'complete'`,
#'   `'average'` (= UPGMA), `'mcquitty'` (= WPGMA), `'median'` (=
#'   WPGMC) or `'centroid'` (= UPGMC).
#' @param nclust The number of clusters to be formed. Set to `NA`
#' @return
#' * **data** The data that was used to compute the distances.
#'
#' * **cutpoint** The cutpoint of the dendrogram according to Mojena (1977).
#'
#' * **distance** The matrix with the distances.
#'
#' * **de** The distances in an object of class `dist`.
#'
#' * **hc** The hierarchical clustering.
#'
#' * **Sqt** The total sum of squares.
#'
#' * **tab** A table with the clusters and similarity.
#'
#' * **clusters** The sum of square and the mean of the clusters for each
#'   variable.
#'
#' * **cofgrap** If `selectvar = TRUE`, then, `cofpgrap` is a
#'   ggplot2-based graphic showing the cophenetic correlation for each model
#'   (with different number of variables). Else, will be a `NULL` object.
#'
#' * **statistics** If `selectvar = TRUE`, then, `statistics` shows
#'   the summary of the models fitted with different number of variables,
#'   including cophenetic correlation, Mantel's correlation with the original
#'   distances (all variables) and the p-value associated with the Mantel's
#'   test. Else, will be a `NULL` object.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Mojena, R. 2015. Hierarchical grouping methods and stopping
#'   rules: an evaluation. Comput. J. 20:359-363. \doi{10.1093/comjnl/20.4.359}
#'
#' @export
#' @examples
#' \donttest{
#' library(metan)
#'
#' # All rows and all numeric variables from data
#' d1 <- clustering(data_ge2)
#'
#' # Based on the mean for each genotype
#' mean_gen <-
#'  data_ge2 %>%
#'  mean_by(GEN) %>%
#'  column_to_rownames("GEN")
#'
#' d2 <- clustering(mean_gen)
#'
#'
#' # Select variables for compute the distances
#' d3 <- clustering(mean_gen, selvar = TRUE)
#'
#' # Compute the distances with standardized data
#' # Define 4 clusters
#' d4 <- clustering(data_ge,
#'                  by = ENV,
#'                  scale = TRUE,
#'                  nclust = 4)
#'
#'}
clustering <- function(.data,
                       ...,
                       by = NULL,
                       scale = FALSE,
                       selvar = FALSE,
                       verbose = TRUE,
                       distmethod = "euclidean",
                       clustmethod = "average",
                       nclust = NA) {
  if (!missing(by)){
    if(length(as.list(substitute(by))[-1L]) != 0){
      stop("Only one grouping variable can be used in the argument 'by'.\nUse 'group_by()' to pass '.data' grouped by more than one variable.", call. = FALSE)
    }
    .data <- group_by(.data, {{by}})
  }
  if(is_grouped_df(.data)){
    results <- .data %>%
      doo(clustering,
          ...,
          scale = scale,
          selvar = selvar,
          verbose = verbose,
          distmethod = distmethod,
          clustmethod = clustmethod,
          nclust = nclust)
    return(add_class(results, "group_clustering"))
  }
  if (scale == TRUE && selvar == TRUE) {
    stop("It is not possible to execute the algorithm for variable selection when 'scale = TRUE'. Please, verify.")
  }
  if (!distmethod %in% c("euclidean", "maximum", "manhattan",
                         "canberra", "binary", "minkowski", "pearson", "spearman",
                         "kendall")) {
    stop("The argument 'distmethod' is incorrect. It should be one of the 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski', 'pearson', 'spearman', or 'kendall'.")
  }
  if (!clustmethod %in% c("complete", "ward.D", "ward.D2",
                          "single", "average", "mcquitty", "median", "centroid")) {
    stop("The argument 'distmethod' is incorrect. It should be one of the 'ward.D', 'ward.D2', 'single', 'average', 'mcquitty', 'median' or 'centroid'.")
  }
  if (missing(...)){
    data <- select_numeric_cols(.data)
    if(has_na(data)){
      data <- remove_rows_na(data)
      has_text_in_num(data)
    }
  } else{
    data <- select(.data, ...) %>%
      select_numeric_cols()
    if(has_na(data)){
      data <- remove_rows_na(data)
      has_text_in_num(data)
    }
  }
  if (scale == TRUE) {
    data <- data.frame(scale(data, center = FALSE, scale = apply(data, 2, sd, na.rm = TRUE)))
  } else {
    data <- data
  }
  if (selvar == TRUE) {
    n <- (ncol(data) - 1)
    statistics <- data.frame(matrix(nrow = n, ncol = 6))
    ModelEstimates <- list()
    modelcode <- 1
    namesv <- "-"
    original <- data
    if (distmethod %in% c("pearson", "spearman", "kendall")) {
      dein <- as.dist(cor(t(data), method = distmethod))
    } else {
      dein <- dist(data, method = distmethod, diag = T,
                   upper = T)
    }
    pb <- progress(max = n)
    for (i in 1:n) {
      if (distmethod %in% c("pearson", "spearman", "kendall")) {
        de <- as.dist(cor(t(data), method = distmethod))
      } else {
        de <- dist(data, method = distmethod, diag = T,
                   upper = T)
      }
      hc <- hclust(de, method = clustmethod)
      d2 <- cophenetic(hc)
      cof <- cor(d2, de)
      mant <- mantel_test(de, dein, nboot = 1000)
      mantc <- mant[[1]]
      mantp <- mant[[3]]
      evect <- data.frame(t(prcomp(data)$rotation))
      var <- abs(evect)[nrow(evect), ]
      names <- apply(var, 1, function(x) which(x == max(x)))
      npred <- ncol(data)
      statistics[i, 1] <- paste("Model", modelcode)
      statistics[i, 2] <- namesv
      statistics[i, 3] <- cof
      statistics[i, 4] <- npred
      statistics[i, 5] <- mantc
      statistics[i, 6] <- mantp
      mat <- as.matrix(de)
      mat <- as.data.frame(mat)
      Results <- list(nvars = npred, excluded = namesv,
                      namevars = names(data), distance = mat, cormantel = mantc,
                      pvmant = mantp)
      namesv <- names(data[names])
      data2 <- data.frame(data[-(match(c(namesv), names(data)))])
      data <- data2
      ModelEstimates[[paste("Model", modelcode)]] <- Results
      names(statistics) <- c("Model", "excluded", "cophenetic",
                             "remaining", "cormantel", "pvmantel")
      if (verbose == TRUE) {
        run_progress(pb, actual = i,
                     sleep = 0.1,
                     text = paste0(namesv, " excluded in this step"))
      }
      modelcode <- modelcode + 1
    }
    cat("--------------------------------------------------------------------------\n")
    cat("\nSummary of the adjusted models", "\n")
    cat("--------------------------------------------------------------------------\n")
    print(statistics, row.names = F)
    cat("--------------------------------------------------------------------------\n")
    model <- statistics$Model[which.max(statistics$cophenetic)]
    predvar <- ModelEstimates[[model]]$namevars
    data <- data.frame(original[(match(c(predvar), names(original)))])
    if (verbose == TRUE) {
      cat("Suggested variables to be used in the analysis\n")
      cat("--------------------------------------------------------------------------\n")
      cat("The clustering was calculated with the ",
          model, "\nThe variables included in this model were...\n",
          predvar, "\n")
      cat("--------------------------------------------------------------------------\n\n")
    }
  } else {
    data <- data
    cofgrap <- NULL
    statistics <- NULL
  }
  if (distmethod %in% c("pearson", "spearman", "kendall")) {
    de <- as.dist(1 - cor(t(data), method = distmethod))
  } else {
    de <- dist(data, method = distmethod, diag = T, upper = T)
  }
  mat <- as.matrix(de)
  hc <- hclust(de, method = clustmethod)
  d2 <- cophenetic(hc)
  cof <- cor(d2, de)
  k <- 1.25
  pcorte <- mean(hc$height) + k * sd(hc$height)
  if (!is.na(nclust)) {
    groups <- cutree(hc, k = nclust)
    Mgroups <- cbind(data, groups)
    distance <- hc$height[(length(hc$height) - nclust):length(hc$height)]
    Sim <- (1 - distance/max(de))
    Passos <- 1:length(Sim)
    Simgroups <- length(Sim):1
    similarity <- Sim * 100
    Tab <- cbind(Passos, Simgroups, round(similarity,
                                          3), round(distance, 2))
    colnames(Tab) <- c("Steps", "Groups", "Similarity",
                       "Distance")
    TabResgroups <- NULL
    MGr <- cbind(data, groups)
    for (i in 1:nclust) {
      NewGroups <- subset(MGr, groups == i)
      GrupCalc <- NewGroups[, 1:(ncol(NewGroups) -
                                   1)]
      Qtd.Elementos <- nrow(NewGroups)
      if (Qtd.Elementos == 1)
        Media <- GrupCalc else Media <- apply(GrupCalc, 2, mean)
      if (Qtd.Elementos == 1)
        SqG <- 0 else SqG <- sum(sweep(GrupCalc, 2, Media)^2)
      TabResgroups <- rbind(TabResgroups, c(i, Qtd.Elementos,
                                            SqG, Media))
    }
    colnames(TabResgroups) <- c("Cluster", "Number of Elements",
                                "Sum_sq", paste(colnames(TabResgroups[, 4:(ncol(TabResgroups))])))
  } else {
    TabResgroups <- NULL
    Tab <- NULL
    Mgroups <- NULL
  }
  Sqt <- sum(sweep(data, 2, apply(data, 2, mean))^2)
  return(list(data = Mgroups,
              cutpoint = pcorte,
              distance = mat,
              de = de,
              hc = as.dendrogram(hc),
              cophenetic = cof,
              Sqt = Sqt,
              tab = as.data.frame(Tab),
              clusters = as.data.frame(TabResgroups),
              statistics = statistics) %>%
           set_class("clustering"))
}



#' Plot an object of class clustering
#'
#' Plot an object of class clustering
#'
#'
#' @param x An object of class `clustering`
#' @param horiz Logical indicating if the dendrogram should be drawn
#'   horizontally or not.
#' @param type The type of plot. Must be one of the 'dendrogram' or
#'   'cophenetic'.
#' @param ... Other arguments passed from the function `plot.dendrogram` or
#'   `abline`.
#' @return An object of class `gg, ggplot` if `type == "cophenetic"`.
#' @method plot clustering
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' mean_gen <-
#'  data_ge2 %>%
#'  mean_by(GEN) %>%
#'  column_to_rownames("GEN")
#'
#' d <- clustering(mean_gen)
#' plot(d, xlab = "Euclidean Distance")
#'}
plot.clustering <- function(x, horiz = TRUE, type = "dendrogram", ...){
  if (type == "dendrogram"){
    plot(x$hc, horiz = horiz, ...)
    if(horiz == TRUE){
      abline(v = x$cutpoint, ...)
    }
    if(horiz == FALSE){
      abline(h = x$cutpoint, ...)
    }
  }
  if (type == "cophenetic"){
    ggplot2::ggplot(x$statistics, ggplot2::aes(x = remaining, y = cophenetic))+
      ggplot2::geom_point(size = 3)+
      ggplot2::theme_bw()+
      ggplot2::geom_line(size = 1)+
      ggplot2::theme(axis.ticks.length = unit(.2, "cm"),
                     axis.text = ggplot2::element_text(size = 12, colour = "black"),
                     axis.title = ggplot2::element_text(size = 12, colour = "black"),
                     axis.ticks = ggplot2::element_line(colour = "black"),
                     plot.margin = margin(0.5, 0.5, 0.2, 0.6, "cm"),
                     axis.title.y = ggplot2::element_text(margin = margin(r=16)),
                     legend.title = ggplot2::element_blank(),
                     legend.text = ggplot2::element_text(size=12),
                     panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1),
                     panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank(),
                     panel.grid.minor.y = ggplot2::element_blank())+
      ggplot2::labs(x = "Number of variables", y = "Cophenetic correlation")
  }
}
