#' Adjusted Coefficient of Variation as yield stability index
#'
#' Performs a stability analysis based on the scale-adjusted coefficient of
#' variation (Doring and Reckling, 2018). For more details see
#' \code{\link{acv}()}
#'
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure use, for example, \code{resp = c(var1, var2, var3)}.
#' @param verbose Logical argument. If \code{verbose = FALSE} the code will run
#'   silently.
#' @return An object of class \code{ge_acv}, which is a list containing the
#'   results for each variable used in the argument \code{resp}. For each
#'   variable, a tibble with the following columns is returned.
#' * \strong{GEN} the genotype's code.
#' * \strong{ACV} The adjusted coefficient of variation
#' * \strong{ACV_R} The rank for the ACV value.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Doring, T.F., and M. Reckling. 2018. Detecting global trends of
#'   cereal yield stability by adjusting the coefficient of variation. Eur. J.
#'   Agron. 99: 30-36. \doi{10.1016/j.eja.2018.06.007}
#'
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' out <- ge_acv(data_ge2, ENV, GEN, c(EH, PH, EL, CD, ED, NKE))
#' gmd(out)
#' }
#'
ge_acv <- function(.data, env, gen, resp, verbose = TRUE) {
  factors  <-
    .data %>%
    select({{env}}, {{gen}}) %>%
    mutate(across(everything(), as.factor))
  vars <-
    .data %>%
    select({{resp}}, -names(factors)) %>%
    select_numeric_cols()
  factors %<>% set_names("ENV", "GEN")
  listres <- list()
  nvar <- ncol(vars)
  if (verbose == TRUE) {
    pb <- progress_bar$new(
      format = "Evaluating the variable :what [:bar]:percent",
      clear = FALSE, total = nvar, width = 90)
  }
  for (var in 1:nvar) {
    data <-
      factors %>%
      mutate(Y = vars[[var]])
    if(has_na(data)){
      data <- remove_rows_na(data)
      has_text_in_num(data)
    }
    temp <- make_mat(data, ENV, GEN, Y)
    varamo <- sapply(temp, var)
    means <- sapply(temp, mean)
    results <-
      acv(means, varamo) %>%
      select_cols(acv) %>%
      colnames_to_upper() %>%
      add_cols(GEN = names(means),
               ACV_R = rank(ACV)) %>%
      reorder_cols(GEN, .before = ACV)
    if (verbose == TRUE) {
      pb$tick(tokens = list(what = names(vars[var])))
    }
    listres[[paste(names(vars[var]))]] <- results
  }
  return(structure(listres, class = "ge_acv"))
}
