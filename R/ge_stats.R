#' Statistics for genotype-vs-environment interaction
#'
#' Computes \strong{(i)} within-environment analysis of variance, GEI effect,
#' GEI means, and genotype plus GEI effects; \strong{(ii)} parametric statistics
#' including Annicchiarico's genotypic confidence index (1992), Ecovalence
#' (Wricke, 1965), regression-based stability (Eberhart and Russell., 1966),
#' Shukla's stability variance parameter (1972); and \strong{(iii)}
#' nonparametric Fox's stability function (Fox et al. 1990), and  superiority
#' index (Lin and Binns, 1988)
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#'   replications/blocks.
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure use, for example, \code{resp = c(var1, var2, var3)}.
#' @param verbose Logical argument. If \code{verbose = FALSE} the code will run
#'   silently.
#' @param prob The probability error assumed.
#' @return An object of class \code{ge_stats} with the following items for each
#'   variable:
#' * \strong{individual} The individual analysis of variance. Call the function
#' \code{\link{anova_ind}} internally.
#' * \strong{ge_mean}, ge_effect The genotype-vs-environment means and effects,
#' respectively.
#' * \strong{gge_effect} The genotype plus genotype-vs-environment effects.
#' * \strong{ge_stats} The variance, sum of squares and mean squares for each
#' genotype.
#' * \strong{ER} The Eberhart and Russell regression model. Call the function
#' \code{\link{ge_reg}} internally.
#' * \strong{ANN} The Annicchiarico's genotypic confidence index. Call the
#' function \code{\link{Annicchiarico}} internally.
#' * \strong{Ecoval} The Wrike's ecovalence. Call the function
#' \code{\link{ecovalence}} internally.
#' * \strong{Shukla} The Shukla's stability variance. Call the function
#' \code{\link{Shukla}} internally.
#' * \strong{Fox} The Fox's stability function. Call the function
#' \code{\link{Fox}} internally.
#' * \strong{Superiority} The Lin and Binns' superiority index. Call the
#' function \code{\link{superiority}} internally.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references
#'   Annicchiarico, P. 1992. Cultivar adaptation and recommendation from alfalfa
#'   trials in Northern Italy. Journal of Genetic \& Breeding, 46:269-278
#'
#'   Eberhart, S.A., and W.A. Russell. 1966. Stability parameters for
#'   comparing Varieties. Crop Sci. 6:36-40.
#'   \href{https://www.crops.org/publications/cs/abstracts/6/1/CS0060010036}{doi:10.2135/cropsci1966.0011183X000600010011x}
#'
#'   Fox, P.N., B. Skovmand, B.K. Thompson, H.J. Braun, and R.
#'   Cormier. 1990. Yield and adaptation of hexaploid spring triticale.
#'   Euphytica 47:57-64.
#'   \href{https://link.springer.com/article/10.1007/BF00040364}{doi:10.1007/BF00040364}.
#'
#'   Kang, M.S., and H.N. Pham. 1991. Simultaneous Selection for High
#'   Yielding and Stable Crop Genotypes. Agron. J. 83:161.
#'   \href{https://dl.sciencesocieties.org/publications/aj/abstracts/83/1/AJ0830010161}{doi:10.2134/agronj1991.00021962008300010037x}.
#'
#'   Lin, C.S., and M.R. Binns. 1988. A superiority measure of cultivar
#'   performance for cultivar x location data. Can. J. Plant Sci. 68:193-198.
#'   \href{http://pubs.aic.ca/doi/abs/10.4141/cjps88-018}{doi:10.4141/cjps88-018}
#'
#'   Shukla, G.K. 1972. Some statistical aspects of partitioning
#'   genotype-environmental components of variability. Heredity. 29:238-245.
#'   \href{http://www.nature.com/articles/hdy197287}{doi:10.1038/hdy.1972.87}.
#'
#'   Wricke, G. 1965. Zur berechnung der okovalenz bei sommerweizen und hafer.
#'   Z. Pflanzenzuchtg 52:127-138.
#'
#'
#' @export
#' @examples
#'
#' library(metan)
#'
#' model = ge_stats(.data = data_ge2,
#'                  env = ENV,
#'                  gen = GEN,
#'                  rep = REP,
#'                  resp = PH)
#'
#' #Alternatively, using the forward-pipe operator %>%
#' #More than one trait
#'
#' model2 = data_ge2 %>%
#'          ge_stats(ENV, GEN, REP, c(PH, EH, EL))
#'
#'
ge_stats = function(.data,
                    env,
                    gen,
                    rep,
                    resp,
                    verbose = TRUE,
                    prob = 0.05){
  datain <- .data
  GEN <- factor(eval(substitute(gen), eval(datain)))
  ENV <- factor(eval(substitute(env), eval(datain)))
  REP <- factor(eval(substitute(rep), eval(datain)))
  listres <- list()
  d <- match.call()
  nvar <- as.numeric(ifelse(length(d$resp) > 1, length(d$resp) - 1, length(d$resp)))
  for (var in 2:length(d$resp)) {
    if (length(d$resp) > 1) {
      Y <- eval(substitute(resp)[[var]], eval(datain))
      varnam = paste(d$resp[var])
    } else {
      Y <- eval(substitute(resp), eval(datain))
      varnam = paste(d$resp)
    }
    data <- data.frame(ENV, GEN, REP, Y)
    names(data) = c("ENV", "GEN", "REP", "mean")
individual <- data %>% anova_ind(ENV, GEN, REP, mean)
ge_mean = make_mat(data, GEN, ENV, mean)
ge_effect = ge_effects(data, ENV, GEN, REP, mean)[[1]]
gge_effect = ge_effects(data, ENV, GEN, REP, mean, type = "gge")[[1]]
Mean = apply(ge_mean, 1, mean)
Variance = rowSums(apply(ge_mean, 2, function(x) (x - Mean)^2))
GENSS = (rowSums(ge_mean^2) - (rowSums(ge_mean)^2)/length(unique(ENV)))*length(unique(REP))
mod = anova(lm(mean ~ REP/ENV + ENV * GEN, data = data))
GENMS = GENSS / mod[2, 1]
Fcal =  GENMS / mod[6, 3]
pval = pf(Fcal, mod[2, 1], mod[6, 1], lower.tail = FALSE)
ec_mod = ecovalence(data, ENV, GEN, REP, mean, verbose = FALSE)[[1]]
lb_mod = superiority(data, ENV, GEN, REP, mean, verbose = FALSE)[[1]]
an_mod = Annicchiarico(data, ENV, GEN, REP, mean, verbose = FALSE, prob = prob)
shu_mod <- Shukla(data, ENV, GEN, REP, mean, verbose = FALSE)[[1]]
fox_mod <- Fox(data, ENV, GEN, REP, mean, verbose = FALSE)[[1]]
ANN = list(environments = an_mod[[1]]$environments,
           general = an_mod[[1]]$general,
           favorable = an_mod[[1]]$favorable,
           unfavorable = an_mod[[1]]$unfavorable)

er_mod = ge_reg(data, ENV, GEN, REP, mean, verbose = FALSE)
ER = list(anova = er_mod[[1]]$anova,
          regression = er_mod[[1]]$regression,
          plot = er_mod[[1]]$plot)

ge_stat = data.frame(Mean = Mean,
                     Variance = Variance,
                     SQE_Gi = GENSS,
                     MSE_Gi = GENMS,
                     Fcal = Fcal,
                     P_val = pval) %>%
  rownames_to_column("GEN")


temp = list(individual = individual[[1]][[1]],
            ge_mean = as_tibble(ge_mean, rownames = NA),
            ge_effect = ge_effect,
            gge_effect = gge_effect,
            ge_stats = as_tibble(ge_stat, rownames = NA),
            ER = ER,
            ANN = ANN,
            Ecoval = ec_mod,
            Shukla = shu_mod,
            Fox = fox_mod,
            Superiority = lb_mod)

if (length(d$resp) > 1) {
  listres[[paste(d$resp[var])]] <- temp
  if (verbose == TRUE) {
    cat("Evaluating variable", paste(d$resp[var]), round((var - 1)/(length(d$resp) -
                                                                      1) * 100, 1), "%", "\n")
  }
} else {
  listres[[paste(d$resp)]] <- temp
}
  }
  return(structure(listres, class = "ge_stats"))
}
