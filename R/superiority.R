#' Lin e Binns' superiority index
#'
#' Nonparametric stability analysis using the superiority index proposed by Lin
#' & Binns, (1992).
#'
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s)
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#'   replications/blocks
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure use, for example, \code{resp = c(var1, var2, var3)}.
#' @param verbose Logical argument. If \code{verbose = FALSE} the code will run
#'   silently.
#' @return An object of class \code{superiority} where each element is the
#'   result of one variable and contains the following items:
#'
#' * \strong{environments} The mean for each environment, the environment index
#' and classification as favorable and unfavorable environments.
#' * \strong{index} The superiority index computed for all, favorable and
#' unfavorable environments.
#'
#' @md
#' @author Tiago Olivoto, \email{tiagoolivoto@@gmail.com}
#' @seealso \code{\link{Annicchiarico}, \link{ecovalence}, \link{ge_stats}}
#' @references Lin, C.S., and M.R. Binns. 1988. A superiority measure of
#'   cultivar performance for cultivar x location data. Can. J. Plant Sci.
#'   68:193-198.
#'   \href{http://pubs.aic.ca/doi/abs/10.4141/cjps88-018}{doi:10.4141/cjps88-018}
#'
#' @export
#' @examples
#'
#' library(metan)
#' out = superiority(data_ge2,
#'                  env = ENV,
#'                  gen = GEN,
#'                  rep = REP,
#'                  resp = PH)
#'
#'
superiority <- function(.data, env, gen, rep, resp, verbose = TRUE) {
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
      varnam <- paste(d$resp[var])
    } else {
      Y <- eval(substitute(resp), eval(datain))
      varnam <- paste(d$resp)
    }
    data <- data.frame(ENV, GEN, REP, Y)
    names(data) <- c("ENV", "GEN", "REP", "mean")
    ge_mean <- data %>% dplyr::group_by(ENV, GEN) %>% dplyr::summarise(mean = mean(mean))
    environments <- data %>% dplyr::group_by(ENV) %>% dplyr::summarise(Mean = mean(mean))
    environments <- mutate(environments, index = Mean - mean(environments$Mean),
                           class = ifelse(index < 0, "unfavorable", "favorable")) %>%
      as.data.frame()
    data <- suppressMessages(left_join(data, environments %>%
                                         select(ENV, class)))
    lin_fun <- function(mat) {
      P <- apply(mat, 1, function(x) {
        sum((x - apply(mat, 2, max))^2)/(2 * length(x))
      })
      return(P)
    }
    mat_g <- make_mat(data, row = GEN, col = ENV, value = mean)
    ge_mf <- subset(data, class == "favorable")
    mat_f <- dplyr::select_if(make_mat(ge_mf, row = GEN,
                                       col = ENV, value = mean), function(x) !any(is.na(x)))
    ge_mu <- subset(data, class == "unfavorable")
    mat_u <- dplyr::select_if(make_mat(ge_mu, row = GEN,
                                       col = ENV, value = mean), function(x) !any(is.na(x)))
    temp <- list(environments = tibble(environments), index = tibble(Genotypes = rownames(mat_g),
                                                                     Pi_all = lin_fun(mat_g), Or_a = rank(lin_fun(mat_g)),
                                                                     Pi_favorable = lin_fun(mat_f), Or_f = rank(lin_fun(mat_f)),
                                                                     Pi_unfavorable = lin_fun(mat_u), Or_u = rank(lin_fun(mat_u))))
    rownames(temp) <- NULL
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
  return(structure(listres, class = "superiority"))
}
