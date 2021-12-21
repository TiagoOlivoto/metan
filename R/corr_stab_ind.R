#' Correlation between stability indexes
#' @description
#' `r badge('stable')`
#'
#' Computes the Spearman's rank correlation between the parametric and
#' nonparametric stability indexes computed with the function
#' [ge_stats()].
#'
#' @param x An object of class `ge_stats`.
#' @param stats The statistics to compute the correlation. See the section
#'   **Details** for more information.
#' @param plot Plot the heat map with the correlations? Defaults to `TRUE`.
#' @param ... Other arguments to be passed to the function
#'   [plot.corr_coef()].
#' @details The argument `stats` is used to chose the statistics to show the
#'   ranks. Allowed values are `"all"` (All statistics, default), `"par"`
#'   (Parametric statistics), `"nonpar"` (Non-parametric statistics), `"ammi"`
#'   (AMMI-based stability statistics), or the following values that can be
#'   combined into comma-separated character vector. `"Y"` (Response variable),
#'   `"Var"` (Genotype's variance), `"Shukla"` (Shukla's variance), `"Wi_g",
#'   "Wi_f", "Wi_u"` (Annichiarrico's genotypic confidence index for all,
#'   favorable and unfavorable environments, respectively), `"Ecoval"` (Wricke's
#'   ecovalence), `"Sij"` (Deviations from the joint-regression analysis),
#'   `"R2"` (R-squared from the joint-regression analysis), `"ASTAB"` (AMMI
#'   Based Stability Parameter), `"ASI"` (AMMI Stability Index), `"ASV"`
#'   (AMMI-stability value), `"AVAMGE"` (Sum Across Environments of Absolute
#'   Value of GEI Modelled by AMMI ), `"Da"` (Annicchiarico's D Parameter
#'   values), `"Dz"` (Zhang's D Parameter), `"EV"` (Sums of the Averages of the
#'   Squared Eigenvector Values), `"FA"` (Stability Measure Based on Fitted AMMI
#'   Model), `"MASV"` (Modified AMMI Stability Value), `"SIPC"` (Sums of the
#'   Absolute Value of the IPC Scores), `"Za"` (Absolute Value of the Relative
#'   Contribution of IPCs to the Interaction), `"WAAS"` (Weighted average of
#'   absolute scores), `"HMGV"` (Harmonic mean of the genotypic value), `"RPGV"`
#'   (Relative performance of the genotypic values), `"HMRPGV"` (Harmonic mean
#'   of the relative performance of the genotypic values), `"Pi_a", "Pi_f",
#'   "Pi_u"` (Superiority indexes for all, favorable and unfavorable
#'   environments, respectively), `"Gai"` (Geometric adaptability index), `"S1"`
#'   (mean of the absolute rank differences of a genotype over the n
#'   environments), `"S2"` (variance among the ranks over the k environments),
#'   `"S3"` (sum of the absolute deviations), `"S6"` (relative sum of squares of
#'   rank for each genotype), `"N1", "N2", "N3", "N4"` (Thennarasu"s
#'   statistics)).
#' @return A list with the data (ranks) correlation, p-values and a heat map showing the
#'   correlation coefficients.
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'
#' @examples
#' \donttest{
#' library(metan)
#' model <- ge_stats(data_ge, ENV, GEN, REP, GY)
#' a <- corr_stab_ind(model)
#' }
#'
corr_stab_ind <- function(x, stats = "all", plot = TRUE, ...){
  all_s <- c("Y", "Var", "Shukla", "Wi_g", "Wi_f", "Wi_u", "Ecoval", "Sij", "R2","ASI", "ASV", "AVAMGE", "DA","DZ","EV","FA","MASI","MASV","SIPC","ZA","WAAS","HMGV", "RPGV", "HMRPGV", "Pi_a", "Pi_f", "Pi_u", "Gai", "S1", "S2", "S3", "S6", "N1", "N2", "N3", "N4")
  par_s <- c("Y", "Var", "Shukla", "Wi_g", "Wi_f", "Wi_u", "Ecoval", "Sij", "R2","ASI", "ASV", "AVAMGE", "DA","DZ","EV","FA","MASI","MASV","SIPC","ZA","WAAS","HMGV", "RPGV", "HMRPGV")
  nonpar_s <- c("Y", "Pi_a", "Pi_f", "Pi_u", "Gai", "S1", "S2", "S3", "S6", "N1", "N2", "N3", "N4" )
  ammi_s <- c("Y", "ASI", "ASV", "AVAMGE", "DA","DZ","EV","FA","MASI","MASV","SIPC","ZA","WAAS")
  if(!stats %in% c("all", "par", "nonpar", "ammi")){
    stats = unlist(strsplit(stats, split=", "))
  } else {
    if(any(stats == "all")){
      stats = all_s
    }
    if(any(stats == "par")){
      stats = par_s
    }
    if(any(stats == "nonpar")){
      stats = nonpar_s
    }
    if(any(stats == "ammi")){
      stats = ammi_s
    }
  }
  if(any(!stats %in% c("all", "par", "nonpar", "ammi", all_s)) == TRUE){
    stop("Argument 'stats' with invalid values. See ?corr_stab_ind for more details.", call. = FALSE)
  }
  bind <- do.call(
    cbind,
    lapply(x, function(x) {
      x %>% select(contains("_R"))
    })) %>%
    as_tibble() %>%
    mutate(gen = x[[1]][["GEN"]]) %>%
    pivot_longer(cols = contains(".")) %>%
    separate(name, into = c("var", "stat"), sep = "(\\.)") %>%
    separate(stat, into = "stat", sep = "_(?=[^_]*$)", extra = "drop") %>%
    pivot_wider(values_from = value, names_from = stat) %>%
    select(-c(var, gen))
  bind <- select(bind, stats)
  corr <- corr_coef(bind)
  p <- plot(corr, legend.title = "Speraman's\nCorrelation", ...)
  if(plot == TRUE){
    suppressWarnings(print(p))
  }
  invisible(structure(list(data = bind,
                        corr = corr$cor,
                        pval = corr$pval,
                        plot = p),
                   class = "corr_stab"))
}
