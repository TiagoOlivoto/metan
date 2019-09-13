#' Stability analysis based on Wricke's model
#'
#' The function computes the ecovalence (Wricke, 1965) for stability analysis.
#'
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
#' @return An object of class \code{ecovalence} containing the results for each
#' variable used in the argument \code{resp}.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Wricke, G. 1965. Zur berechnung der okovalenz bei sommerweizen
#' und hafer. Z. Pflanzenzuchtg 52:127-138.
#' @export
#' @examples
#'
#' \dontrun{
#' library(metan)
#' out = ecovalence(data_ge2,
#'                  env = ENV,
#'                  gen = GEN,
#'                  rep = REP,
#'                  resp = PH)
#' }
#'
ecovalence <- function(.data, env, gen, rep, resp, verbose = TRUE) {
  datain <- .data
  GEN <- factor(eval(substitute(gen), eval(datain)))
  ENV <- factor(eval(substitute(env), eval(datain)))
  REP <- factor(eval(substitute(rep), eval(datain)))
  listres <- list()
  d <- match.call()
  nvar <- as.numeric(ifelse(length(d$resp) > 1, length(d$resp) -
                              1, length(d$resp)))
  for (var in 2:length(d$resp)) {
    if (length(d$resp) > 1) {
      Y <- eval(substitute(resp)[[var]], eval(datain))
    } else {
      Y <- eval(substitute(resp), eval(datain))
    }
    data <- data.frame(ENV, GEN, REP, Y)
    names(data) <- c("ENV", "GEN", "REP", "mean")
    data2 <- data %>% dplyr::group_by(ENV, GEN) %>% dplyr::summarise(mean = mean(mean)) %>%
      as.data.frame()
    data3 <- mutate(data2, ge = residuals(lm(mean ~ ENV +
                                               GEN, data = data2)))
    ge_effect <- make_mat(data3, GEN, ENV, ge)
    Ecoval <- rowSums(ge_effect^2 * length(unique(REP)))
    Ecov_perc <- (Ecoval/sum(Ecoval)) * 100
    rank <- rank(Ecoval)
    temp <- as_tibble(cbind(ge_effect, Ecoval, Ecov_perc,
                            rank), rownames = NA)
    if (length(d$resp) > 1) {
      listres[[paste(d$resp[var])]] <- temp
      if (verbose == TRUE) {
        cat("Evaluating variable", paste(d$resp[var]),
            round((var - 1)/(length(d$resp) - 1) * 100,
                  1), "%", "\n")
      }
    } else {
      listres[[paste(d$resp)]] <- temp
    }
  }
  return(structure(listres, class = "ecovalence"))
}
