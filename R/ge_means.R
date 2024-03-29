#' Genotype-environment means
#' @description
#' `r badge('stable')`
#'
#' Computes genotype-environment interaction means
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, and the response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param resp The response variable(s). To analyze multiple variables at once,
#'   a vector of variables may be used. For example `resp = c(var1, var2,
#'   var3)`. Select helpers are also allowed.
#' @return A list where each element is the result for one variable containing:
#' * **ge_means**: A two-way table with the means for genotypes (rows) and
#' environments (columns).
#' * **gen_means**: A tibble with the means for genotypes.
#' * **env_means**: A tibble with the means for environments.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' means_ge <- ge_means(data_ge, ENV, GEN, resp = everything())
#'
#' # Genotype-environment interaction means
#' get_model_data(means_ge)
#'
#' # Environment means
#' get_model_data(means_ge, what = "env_means")
#'
#' # Genotype means
#' get_model_data(means_ge, what = "gen_means")
#'
#' }
#'
ge_means <- function(.data, env, gen, resp) {
  factors  <-
    .data %>%
    select({{env}}, {{gen}}) %>%
    mutate(across(everything(), as.factor))
  vars <- .data %>% select({{resp}}, -names(factors))
  vars %<>% select_numeric_cols()
  factors %<>% set_names("ENV", "GEN")
  listres <- list()
  nvar <- ncol(vars)
  for (var in 1:nvar) {
    data <- factors %>%
      add_cols(Mean = vars[[var]])
    if(has_na(data)){
      data <- remove_rows_na(data)
      has_text_in_num(data)
    }
    ge_m <- data %>%
      mean_by(ENV, GEN)
    g_m <- mean_by(data, GEN)
    e_m <- mean_by(data, ENV)
    listres[[paste(names(vars[var]))]] <-
      structure(list(ge_means_mat =  make_mat(ge_m, GEN, ENV, Mean),
                     ge_means_long = ge_m,
                     gen_means = g_m,
                     env_means = e_m),
                class = "ge_means")
  }
  return(structure(listres, class = "ge_means"))
}
