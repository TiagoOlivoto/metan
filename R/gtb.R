#' Genotype by trait biplot
#' @description
#' `r badge('stable')`
#'
#' Produces a genotype-by-trait biplot model. From a genotype by environment by
#' trait three-way table, genotype-by-trait tables in any single environment,
#' across all environments, or across a subset of the environments can be
#' generated and visually studied using biplots. The model for biplot analysis
#' of genotype by trait data is the singular value decomposition of
#' trait-standardized two-way table.
#'
#'
#' @param .data The dataset containing the columns related to Genotypes and the
#'   response variable(s).
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param resp The response variables, i.e., `resp = c(var1, var2, var3)`.
#'   Select helpers can also be used.
#' @param centering The centering method. Must be one of the `'none | 0'`,
#'   for no centering; `'global | 1'`, for global centered (T+G+GT);
#'   `'trait | 2'` (default), for trait-centered (G+GT); or `'double |
#'   3'`, for double centred (GT). A biplot cannot be produced with models
#'   produced without centering.
#' @param scaling The scaling method. Must be one of the `'none | 0'`, for
#'   no scaling; or `'sd | 1'` (default), where each value is divided by
#'   the standard deviation of its corresponding trait (column). This will put
#'   all traits roughly he same rang of values.
#' @param svp The method for singular value partitioning. Must be one of the
#'   `'genotype | 1'`, (The singular value is entirely partitioned into the
#'   genotype eigenvectors, also called row metric preserving); `'trait |
#'   2'`, default, (The singular value is entirely partitioned into the trait
#'   eigenvectors, also called column metric preserving); or `'symmetrical
#'   | 3'` (The singular value is symmetrically partitioned into the genotype
#'   and the trait eigenvectors This SVP is most often used in AMMI analysis and
#'   other biplot analysis, but it is not ideal for visualizing either the
#'   relationship among genotypes or that among the traits).
#'
#'
#' @return The function returns a list of class `gge` that is compatible with the function `plot()` used in [gge()].
#'
#'  * **coordgen** The coordinates for genotypes for all components.
#'
#'  * **coordenv** The coordinates for traits for all components.
#'
#'  * **eigenvalues** The vector of eigenvalues.
#'
#'  * **totalvar** The overall variance.
#'
#'  * **labelgen** The name of the genotypes.
#'
#'  * **labelenv** The names of the traits.
#'
#'  * **labelaxes** The axes labels.
#'
#'  * **gt_mat** The data used to produce the model (scaled and centered).
#'
#'  * **centering** The centering method.
#'
#'  * **scaling** The scaling method.
#'
#'  * **svp** The singular value partitioning method.
#'
#'  * **d** The factor used to generate in which the ranges of genotypes
#'  and traits are comparable when singular value partitioning is set to
#'  'genotype' or 'trait'.
#'  * **grand_mean** The grand mean of the trial.
#'  * **mean_gen** A vector with the means of the genotypes.
#'  * **mean_env** A vector with the means of the traits.
#'  * **scale_var** The scaling vector when the scaling method is `'sd'`.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Yan, W., and M.S. Kang. 2003. GGE biplot analysis: a graphical tool for breeders,
#'  geneticists, and agronomists. CRC Press.
#'
#' @export
#' @examples
#' \donttest{
#' library(metan)

#' # GT biplot for all numeric variables
#' mod <- gtb(data_ge2, GEN, resp = contains("E"))
#' plot(mod)
#'
#'}
gtb <- function(.data,
                gen,
                resp,
                centering = "trait",
                scaling = "sd",
                svp = "trait") {
  factors  <-
    .data %>%
    select({{gen}})
  vars <- .data %>%
    select({{resp}}, -names(factors)) %>%
    select_numeric_cols()
  factors %<>% set_names("GEN")
  data <-
    cbind(factors, vars)
  if(has_na(data)){
    data <- remove_rows_na(data)
    has_text_in_num(data)
  }
    gt_mat <-
    means_by(data, GEN) %>%
    column_to_rownames("GEN") %>%
    as.matrix()
    grand_mean <- mean(gt_mat)
    mean_trait <- colMeans(gt_mat)
    mean_gen <- rowMeans(gt_mat)
    scale_val <- apply(gt_mat, 2, sd)
    labelgen <- rownames(gt_mat)
    labelenv <- colnames(gt_mat)
    if (any(is.na(gt_mat))) {
      stop("missing data in input data frame")
    }
    if (any(apply(gt_mat, 2, is.numeric) == FALSE)) {
      stop("not all columns are of class 'numeric'")
    }
    if (!(centering %in% c("none", "trait", "global", "double") |
          centering %in% 0:3)) {
      warning(paste("Centering method", centering, "not found; defaulting to environment centered"))
      centering <- "trait"
    }
    if (!(svp %in% c("genotype", "trait", "symmetrical") |
          svp %in% 1:3)) {
      warning(paste("svp method", svp, "not found; defaulting to column metric preserving"))
      svp <- "trait"
    }
    if (!(scaling %in% c("none", "sd") | scaling %in% 0:1)) {
      warning(paste("scaling method", scaling, "not found; defaulting to no scaling"))
      sd <- "none"
    }
    labelaxes <- paste("PC", 1:ncol(diag(svd(gt_mat)$d)), sep = "")
    # Centering
    if (centering == 1 | centering == "global") {
      gt_mat <- gt_mat - mean(gt_mat)
    }
    if (centering == 2 | centering == "trait") {
      gt_mat <- sweep(gt_mat, 2, colMeans(gt_mat))
    }
    if (centering == 3 | centering == "double") {
      grand_mean <- mean(gt_mat)
      mean_trait <- colMeans(gt_mat)
      mean_gen <- rowMeans(gt_mat)
      for (i in 1:nrow(gt_mat)) {
        for (j in 1:ncol(gt_mat)) {
          gt_mat[i, j] <- gt_mat[i, j] + grand_mean - mean_trait[j] - mean_gen[i]
        }
      }
    }
    # Scaling
    if (scaling == 1 | scaling == "sd") {
      gt_mat <- sweep(gt_mat, 2, apply(gt_mat, 2, sd), FUN = "/")
    }
    # Singular value partitioning
    if (svp == 1 | svp == "genotype") {
      coordgen <- svd(gt_mat)$u %*% diag(svd(gt_mat)$d)
      coordenv <- svd(gt_mat)$v
      d1 <- (max(coordenv[, 1]) - min(coordenv[, 1]))/(max(coordgen[, 1]) - min(coordgen[, 1]))
      d2 <- (max(coordenv[, 2]) - min(coordenv[, 2]))/(max(coordgen[, 2]) - min(coordgen[, 2]))
      coordenv <- coordenv/max(d1, d2)
    }
    if (svp == 2 | svp == "trait") {
      coordgen <- svd(gt_mat)$u
      coordenv <- svd(gt_mat)$v %*% diag(svd(gt_mat)$d)
      d1 <- (max(coordgen[, 1]) - min(coordgen[, 1]))/(max(coordenv[, 1]) - min(coordenv[, 1]))
      d2 <- (max(coordgen[, 2]) - min(coordgen[, 2]))/(max(coordenv[, 2]) - min(coordenv[, 2]))
      coordgen <- coordgen/max(d1, d2)
    }
    if (svp == 3 | svp == "symmetrical") {
      coordgen <- svd(gt_mat)$u %*% diag(sqrt(svd(gt_mat)$d))
      coordenv <- svd(gt_mat)$v %*% diag(sqrt(svd(gt_mat)$d))
    }
    eigenvalues <- svd(gt_mat)$d
    totalvar <- round(as.numeric(sum(eigenvalues^2)), 2)
    varexpl <- round(as.numeric((eigenvalues^2/totalvar) * 100),
                     2)
    if (svp == "genotype" | svp == "trait") {
      d <- max(d1, d2)
    } else {
      d <- NULL
    }
    tmp <-
      list(coordgen = coordgen, coordenv = coordenv, eigenvalues = eigenvalues,
           totalvar = totalvar, varexpl = varexpl, labelgen = labelgen,
           labelenv = labelenv, labelaxes = labelaxes, ge_mat = gt_mat,
           centering = centering, scaling = scaling, svp = svp,
           d = d, grand_mean = grand_mean, mean_gen = mean_gen,
           mean_env = mean_trait, scale_val = scale_val) %>%
      set_class(c("gge", "gtb"))
  return(set_class(list(mod = tmp), c("gge", "gtb")))
}
