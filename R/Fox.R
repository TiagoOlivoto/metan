#' Fox's stability function
#'
#' Performs a stability analysis based on the criteria of Fox et al. (1990),
#' using the statistical "TOP third" only. A stratified ranking of the genotypes
#' at each environment is done. The proportion of locations at which the
#' genotype occurred in the top third are expressed in the output.
#'
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
#' @return An object of class \code{Fox}, which is a list containing the results
#'   for each variable used in the argument \code{resp}. For each variable, a
#'   tibble with the following columns is returned.
#' * \strong{GEN} the genotype's code.
#' * \strong{mean} the mean for the response variable.
#' * \strong{TOP} The proportion of locations at which the
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Fox, P.N., B. Skovmand, B.K. Thompson, H.J. Braun, and R.
#'   Cormier. 1990. Yield and adaptation of hexaploid spring triticale.
#'   Euphytica 47:57-64.
#'   \href{https://link.springer.com/article/10.1007/BF00040364}{doi:10.1007/BF00040364}.
#'
#' @export
#' @examples
#'
#' library(metan)
#' out = Fox(data_ge2,
#'           env = ENV,
#'           gen = GEN,
#'           rep = REP,
#'           resp = PH)
#'
Fox <- function(.data, env, gen, rep, resp, verbose = TRUE) {
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
    } else {
      Y <- eval(substitute(resp), eval(datain))
    }
    data <- data.frame(ENV, GEN, REP, Y)
    temp <- data %>%
      group_by(ENV) %>%
      mutate(grank = rank(-Y)) %>%
      group_by(GEN) %>%
      summarise(mean = mean(Y),
                TOP = sum(grank <= 3)) %>%
      as_tibble()

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
  return(structure(listres, class = "Fox"))
}
