#' Geometric adaptability index
#'
#' Performs a stability analysis based on the geometric mean (GAI), according to
#' the following model (Mohammadi and Amri, 2008):
#' \loadmathjax
#'  \mjsdeqn{GAI = \sqrt[E]{{\mathop {\bar Y}\nolimits_1  + \mathop {\bar Y}\nolimits_2  + ... + \mathop {\bar Y}\nolimits_i }}}
#'  where \mjseqn{\bar Y_1}, \mjseqn{\bar Y_2}, and \mjseqn{\bar Y_i} are
#'  the mean yields of the first, second and \emph{i}-th genotypes across
#'  environments, and E is the number of environments
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
#' @return An object of class \code{gai}, which is a list containing the results
#'   for each variable used in the argument \code{resp}. For each variable, a
#'   tibble with the following columns is returned.
#' * \strong{GEN} the genotype's code.
#' * \strong{GAI} Geometric adaptability index
#' * \strong{GAI_R} The rank for the GAI value.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references
#' Mohammadi, R., & Amri, A. (2008). Comparison of parametric and non-parametric
#' methods for selecting stable and adapted durum wheat genotypes in variable
#' environments. Euphytica, 159(3), 419-432.
#' \doi{10.1007/s10681-007-9600-6}.
#'
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' out <- gai(data_ge2, ENV, GEN, REP, c(EH, PH, EL, CD, ED, NKE))
#' }
#'
gai <- function(.data, env, gen, rep, resp, verbose = TRUE) {
  factors  <-
    .data %>%
    select({{env}}, {{gen}}, {{rep}}) %>%
    mutate(across(everything(), as.factor))
  vars <-
    .data %>%
    select({{resp}}, -names(factors)) %>%
    select_numeric_cols()
  factors %<>% set_names("ENV", "GEN", "REP")
  listres <- list()
  nvar <- ncol(vars)
  if (verbose == TRUE) {
    pb <- progress_bar$new(
      format = "Evaluating the variable :what [:bar]:percent",
      clear = FALSE, total = nvar, width = 90)
  }
  for (var in 1:nvar) {
    data <- factors %>%
      mutate(Y = vars[[var]])
    if(has_na(data)){
      data <- remove_rows_na(data)
      has_text_in_num(data)
    }
    temp <-
      make_mat(data, ENV, GEN, Y) %>%
      gmean() %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("GEN") %>%
      mutate(rank = rank(-V1)) %>%
      set_names("GEN", "GAI", "GAI_R")
    if (verbose == TRUE) {
      pb$tick(tokens = list(what = names(vars[var])))
    }
    listres[[paste(names(vars[var]))]] <- temp
  }
  return(structure(listres, class = "gai"))
}
