#' @title Power Law Residuals as yield stability index
#' \loadmathjax
#' @description
#' `r badge('stable')`
#'
#' Performs a stability analysis based on the Power Law Residuals (POLAR)
#' statistics (Doring et al., 2015). POLAR is the residuals from the linear
#' regression of \mjseqn{log(\sigma^2}) against \mjseqn{log(\mu}) and can be
#' used as a measure of crop stability with lower stability (relative to all
#' samples with that mean yield) indicated by more positive POLAR values, and
#' higher stability (relative to all samples with that mean yield) indicated by
#' more negative POLAR values.
#'
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure use, for example, `resp = c(var1, var2, var3)`.
#' @param base The base with respect to which logarithms are computed. Defaults
#'   to `10`.
#' @param verbose Logical argument. If `verbose = FALSE` the code will run
#'   silently.
#' @return An object of class `ge_acv`, which is a list containing the
#'   results for each variable used in the argument `resp`. For each
#'   variable, a tibble with the following columns is returned.
#' * **GEN** the genotype's code.
#' * **POLAR** The Power Law Residuals
#' * **POLAR_R** The rank for the ACV value.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Doring, T.F., S. Knapp, and J.E. Cohen. 2015. Taylor's power law
#' and the stability of crop yields. F. Crop. Res. 183: 294-302.
#' \doi{10.1016/j.fcr.2015.08.005}
#'
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' out <- ge_polar(data_ge2, ENV, GEN, c(EH, PH, EL, CD, ED, NKE))
#' gmd(out)
#' }
#'
ge_polar <- function(.data, env, gen, resp, base = 10, verbose = TRUE) {
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
    pb <- progress(max = nvar, style = 4)
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
    if(has_na(temp)){
      warning("Missing values in the GxE matrix\nIndex was computed after removing them.", call. = FALSE)
    }
    varamo <- sapply(temp, var, na.rm = TRUE)
    means <- sapply(temp, mean, na.rm = TRUE)
    mi <- log(means, base = base)
    vi <- log(varamo, base = 10)
    mod <- lm(vi ~ mi)
    ui <- residuals(mod)
    results <- data.frame(
      GEN = names(means),
      POLAR = ui,
      POLAR_R = rank(ui)
    )
    if (verbose == TRUE) {
      run_progress(pb,
                   actual = var,
                   text = paste("Evaluating trait", names(vars[var])))
    }
    listres[[paste(names(vars[var]))]] <- results
  }
  return(structure(listres, class = "ge_polar"))
}
