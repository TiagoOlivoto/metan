#' Statistics for genotype-vs-environment interaction
#'
#' Compute parametric and nonparametric statistics for genotype-vs-environment
#' interaction (GEI), including within-environment analysis of
#' variance, GEI effect, GEI means, genotype plus GEI effects, Ecovalence (Wricke, 1965)
#' regression-based stability (Eberhart and Russell., 1966), and superiority
#' index (Lin and Binns, 1988)
#'
#' @param .data The dataset containing the columns related to Environments,
#' Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#' environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#' replications/blocks.
#' @param resp The response variable(s). To analyze multiple variables in a
#' single procedure use, for example, \code{resp = c(var1, var2, var3)}.
#' @param verbose Logical argument. If \code{verbose = FALSE} the code will run
#' silently.
#' @param prob The probability error assumed.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Eberhart, S.A., and W.A. Russell. 1966. Stability parameters for comparing Varieties.
#' Crop Sci. 6:36-40. \href{https://www.crops.org/publications/cs/abstracts/6/1/CS0060010036}{doi:10.2135/cropsci1966.0011183X000600010011x}
#'
#' Lin, C.S., and M.R. Binns. 1988. A superiority measure of
#' cultivar performance for cultivar x location data. Can. J. Plant Sci.
#' 68:193-198.
#' \href{http://pubs.aic.ca/doi/abs/10.4141/cjps88-018}{doi:10.4141/cjps88-018}
#'
#' Wricke, G. 1965. Zur berechnung der okovalenz bei sommerweizen und hafer. Z.
#' Pflanzenzuchtg 52:127-138.
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
data2 =  data  %>%
  dplyr::group_by(ENV, GEN) %>%
  dplyr::summarise(mean = mean(mean)) %>%
  as.data.frame()
data3 = mutate(data2,
               ge = residuals(lm(mean ~ ENV + GEN, data = data2)),
               gge = residuals(lm(mean ~ ENV, data = data2)))
ge_mean = make_mat(data3, GEN, ENV, mean)
ge_effect = make_mat(data3, GEN, ENV, ge)
gge_effect = make_mat(data3, GEN, ENV, gge)
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
            ge_effect = as_tibble(ge_effect, rownames = NA),
            gge_effect = as_tibble(gge_effect, rownames = NA),
            ge_stats = as_tibble(ge_stat, rownames = NA),
            ER = ER,
            ANN = ANN,
            Ecoval = ec_mod,
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
