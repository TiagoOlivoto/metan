#' Annicchiarico's genotypic confidence index
#'
#' Stability analysis using the known genotypic confidence index
#' (Annicchiarico, 1992).
#'
#'
#' @param .data The dataset containing the columns related to Environments,
#' Genotypes, replication/block and response variable(s)
#' @param env The name of the column that contains the levels of the
#' environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#' replications/blocks
#' @param resp The response variable(s). To analyze multiple variables in a
#' single procedure use, for example, \code{resp = c(var1, var2, var3)}.
#' @param prob The probability of error assumed.
#' @param verbose Logical argument. If \code{verbose = FALSE} the code will run silently.
#' @author Tiago Olivoto, \email{tiagoolivoto@@gmail.com}
#' @seealso \code{\link{superiority}, \link{ecovalence}, \link{ge_stats}}
#' @references Annicchiarico, P. 1992. Cultivar adaptation and recommendation
#' from alfalfa trials in Northern Italy. J. Genet. Breed. 46:269-278.
#' @export
#' @examples
#'
#' library(metan)
#' Ann = Annicchiarico(data_ge2,
#'                    env = ENV,
#'                    gen = GEN,
#'                    rep = REP,
#'                    resp = PH)
#' summary(Ann)
#'
#'
Annicchiarico = function(.data,
                         env,
                         gen,
                         rep,
                         resp,
                         prob = 0.05,
                         verbose = TRUE){
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
    ge_mean =  data  %>%
      dplyr::group_by(ENV, GEN) %>%
      dplyr::summarise(mean = mean(mean))
    environments = data  %>%
      dplyr::group_by(ENV) %>%
      dplyr::summarise(Mean = mean(mean))
    environments = mutate(environments,
                          index = Mean - mean(environments$Mean),
                          class = ifelse(index < 0, "unfavorable", "favorable")) %>%
      as.data.frame()
    data = suppressMessages(left_join(data, environments %>% select(ENV, class)))
    mat_g  = make_mat(data, row = GEN, col = ENV, value = mean)
    rp_g = sweep(mat_g, 2, colMeans(mat_g), "/")*100
    Wi_g = rowMeans(rp_g)  - qnorm(1- prob) * apply(rp_g, 1, sd)
    general = data.frame(Genotype = rownames(mat_g),
                         Mean = rowMeans(mat_g),
                         Mean_rp = rowMeans(rp_g),
                         Sd_rp = apply(rp_g, 1, sd),
                         Wi = Wi_g,
                         order_Wi = order(-Wi_g))
    ge_mf = subset(data, class == "favorable")
    mat_f  = dplyr::select_if(make_mat(ge_mf, row = GEN, col = ENV, value = mean), function(x) !any(is.na(x)))
    rp_f = sweep(mat_f, 2, colMeans(mat_f), "/")*100
    Wi_f = rowMeans(rp_f)  - qnorm(1- prob) * apply(rp_f, 1, sd)
    favorable = data.frame(Genotype = rownames(mat_f),
                           Mean = rowMeans(mat_f),
                           Mean_rp = rowMeans(rp_f),
                           Sd_rp = apply(rp_f, 1, sd),
                           Wi = Wi_f,
                           order_Wi = order(-Wi_f))
    ge_mu = subset(data, class == "unfavorable")
    mat_u  = dplyr::select_if(make_mat(ge_mu, row = GEN, col = ENV, value = mean), function(x) !any(is.na(x)))
    rp_u = sweep(mat_u, 2, colMeans(mat_u), "/")*100
    Wi_u = rowMeans(rp_u)  - qnorm(1- prob) * apply(rp_u, 1, sd)
    unfavorable = data.frame(Genotype = rownames(mat_u),
                             Mean = rowMeans(mat_u),
                             Mean_rp = rowMeans(rp_u),
                             Sd_rp = apply(rp_u, 1, sd),
                             Wi = Wi_u,
                             order_Wi = order(-Wi_u))
    temp = list(environments = environments,
                general = general,
                favorable = favorable,
                unfavorable = unfavorable)
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
  return(structure(listres, class = "Annicchiarico"))
}
