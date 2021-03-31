#' Genotype by yield*trait biplot
#' @description
#' `r badge('stable')`
#'
#' Produces a Genotype by Yield*Trait biplot (GTY) proposed by Yan and
#' Fregeau-Reid (2018).
#'
#'
#' @param .data The dataset containing the columns related to Genotypes, Yield,
#'   and Traits.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param yield The column containing the yield values.
#' @param traits The column(s) with the *traits* values. Defaults to
#'   *NULL*. In this case, all numeric traits in `.data`, except that
#'   in `yield` are selected. To select specific traits from `.data`,
#'   use a list of unquoted comma-separated variable names (e.g. *traits =
#'   c(var1, var2, var3)*), an specific range of variables, (e.g. *traits =
#'   c(var1:var3)*), or even a select helper like `starts_with("N")`.
#' @param ideotype A vector of `"h"` or `"l"` with the same length of
#'   `traits` to define which trait is desired to increase or decrease. By
#'   default (`ideotype = NULL`) for all numeric traits in `traits`
#'   are assumed that high values is desirable. Following the order of the
#'   traits selected in `traits`, use `"h"` to indicate the traits in
#'   which higher values are desired or `"l"` to indicate the variables in
#'   which lower values are desired. Then, `yield ` will be multiplied by
#'   traits with `"h"` and divided by traits with `"l"` to generate
#'   the Genotype by yield*trait table. For example, `ideotype = c("h, h,
#'   l")` will assume that the ideotype has higher values for the first two
#'   traits and lower values for the last trait.
#' @param weight The weight assumed for each trait. Similar to `ideotype`
#'   argument, provide a numeric vector of the same length of `traits`.
#'   Suggested values are between 0 and 2.
#' @param prefix The prefix used in the biplot for the yield*trait combinations.
#'   Defaults to `"Y"`.
#' @param centering The centering method. Must be one of the `'none | 0'`,
#'   for no centering; `'global | 1'`, for global centered (T+G+GYT);
#'   `'trait | 2'` (default), for trait-centered (G+GYT); or `'double |
#'   3'`, for double centered (GYT). A biplot cannot be produced with models
#'   produced without centering.
#' @param scaling The scaling method. Must be one of the `'none | 0'`, for
#'   no scaling; or `'sd | 1'` (default), so that the mean for each trait
#'   or yield-trait combination becomes 0 and the variance becomes unit.
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
#'  * **data** The Genotype by yield\*trait (GYT) data.
#'
#'  * **ge_mat** The Genotype by yield\*trait (GYT) data  (scaled and centered).
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
#' @references
#'  Yan, W., & Fregeau-Reid, J. (2018). Genotype by Yield\*Trait (GYT) Biplot: a
#'  Novel Approach for Genotype Selection based on Multiple Traits. Scientific
#'  Reports, 8(1), 1-10. \doi{10.1038/s41598-018-26688-8}
#'
#' @export
#' @examples
#' \donttest{
#' library(metan)

#' # GYT biplot for all numeric traits of 'data_g'
#' # KW (kernel weight) considered as 'yield',
#' mod <- gytb(data_g, GEN, KW)
#' plot(mod)
#'
#'}
gytb <- function(.data,
                 gen,
                 yield,
                 traits = everything(),
                 ideotype = NULL,
                 weight = NULL,
                 prefix = "Y",
                 centering = "trait",
                 scaling = "sd",
                 svp = "trait") {

  factors  <-
    .data %>%
    select({{gen}}) %>%
    set_names("GEN")
  yld <-
    .data %>%
    select({{yield}}, -{{gen}}) %>%
    select_numeric_cols() %>%
    pull()
  others <-
    .data %>%
    select({{traits}}, -{{gen}}, -{{yield}}) %>%
    select_numeric_cols()
  if(is.null(weight)){
    weights <- rep(1, length(others))
  } else{
    if(length(weight) != length(others)){
      stop("weight must have length ", ncol(others), ", the number of traits in 'traits' argument.")
    }
    weights <-  weight
  }
  if(is.null(ideotype)){
    ideotype <- rep("h", length(others))
  } else{
    ideotype <- unlist(strsplit(ideotype, split="\\s*(\\s|,)\\s*")) %>%
      all_lower_case()
    if(length(ideotype) != length(others)){
      stop("Ideotype must have length ", ncol(others), ", the number of traits in 'traits' argument.")
    }
    if(!all(ideotype %in% c("h", "l"))){
      stop("argument 'ideotype' must have 'h' or 'l' only", call. = FALSE)
    }
  }
  varinc <- others[, which(ideotype == "h")]
  if(ncol(varinc) > 0 ){
    varinc <-
      varinc %>%
      apply(MARGIN = 2, function(x){
        yld * x
      })
  } else{
    varinc <- varinc
  }
  colnames(varinc) <- paste(prefix, "*", colnames(varinc), sep = "")

  vardec <- others[, which(ideotype == "l")]
  if(ncol(vardec) > 0 ){
    vardec <-
      vardec %>%
      apply(MARGIN = 2, function(x){
        yld / x
      })
  } else{
    vardec <- vardec
  }
  colnames(vardec) <- paste(prefix, "/", colnames(vardec), sep = "")
  data <- cbind(factors, varinc, vardec)
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
    scaling <- "none"
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
  gt_mat <- sweep(gt_mat, 2, STATS = weights, FUN = "*")

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
  if (svp %in% c("genotype", "trait", 1, 2)) {
    d <- max(d1, d2)
  } else {
    d <- NULL
  }
  tmp <-
    list(data = data,
         ge_mat = gt_mat,
         coordgen = coordgen,
         coordenv = coordenv,
         eigenvalues = eigenvalues,
         totalvar = totalvar,
         varexpl = varexpl,
         labelgen = labelgen,
         labelenv = labelenv,
         labelaxes = labelaxes,
         centering = centering,
         scaling = scaling,
         svp = svp,
         d = d,
         grand_mean = grand_mean,
         mean_gen = mean_gen,
         mean_env = mean_trait,
         scale_val = scale_val) %>%
    set_class(c("gge", "gytb"))
  return(set_class(list(mod = tmp), c("gge", "gytb")))
}
