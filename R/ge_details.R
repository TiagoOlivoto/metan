#' Details for genotype-environment trials
#' @description
#' `r badge('stable')`
#'
#' Provide details for genotype-environment trials
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure a vector of variables may be used. For example `resp
#'   = c(var1, var2, var3)`. Select helpers are also allowed.
#'
#' @return A tibble with the following results for each variable:
#' * `Mean`: The grand mean.
#' * `SE`: The standard error of the mean.
#' * `SD`: The standard deviation.
#' * `CV`: The coefficient of variation.
#' * `Min,Max`: The minimum and maximum value, indicating the genotype and environment of occurence.
#' * `MinENV, MinGEN`: The environment and genotype with the lower mean.
#' * `MaxENV, MaxGEN`: The environment and genotype with the higher mean.
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(metan)
#' details <- ge_details(data_ge2, ENV, GEN, everything())
#' print(details)
#' }
ge_details <- function(.data, env, gen, resp){
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
      mutate(Y = vars[[var]])
    if(has_na(data)){
      data <- remove_rows_na(data)
      has_text_in_num(data)
    }
    env_data <- mean_by(data, ENV, na.rm = TRUE) %>%
      add_cols(TYPE = "Env") %>%
      rename(CODE = "ENV")
    gen_data <- mean_by(data, GEN, na.rm = TRUE) %>%
      add_cols(TYPE = "Gen")%>%
      rename(CODE = "GEN")
    df_bind <-
      rbind(env_data, gen_data) %>%
      column_to_first(TYPE, CODE)
    min_group <-
      df_bind %>%
      group_by(TYPE) %>%
      top_n(1, -Y) %>%
      select(TYPE, CODE, Y) %>%
      slice(1) %>%
      as.data.frame()
    max_group <-
      df_bind %>%
      group_by(TYPE) %>%
      top_n(1, Y) %>%
      select(TYPE, CODE, Y) %>%
      slice(1) %>%
      as.data.frame()
    min <- data %>% top_n(1, -Y) %>% select(ENV, GEN, Y) %>% slice(1)
    max <- data %>% top_n(1, Y) %>% select(ENV, GEN, Y) %>% slice(1)
    desc_st <- desc_stat(data, stats = c("mean, se, sd.pop, cv"), verbose = FALSE, na.rm = TRUE)
    temp <- tibble(Parameters = c("Mean", "SE", "SD", "CV", "Min", "Max", "MinENV", "MaxENV", "MinGEN", "MaxGEN"),
                   Values = c(round(desc_st[1, 2], 2),
                              round(desc_st[1, 3], 2),
                              round(desc_st[1, 4], 2),
                              round(desc_st[1, 5], 2),
                              paste0(round(min[3], 2), " (", min$GEN, " in ", min$ENV,")"),
                              paste0(round(max$Y, 2), " (", max$GEN, " in ", max$ENV,")"),
                              paste0(min_group[1,2], " (", round(min_group[1,3], 2),")"),
                              paste0(max_group[1,2], " (", round(max_group[1,3], 2),")"),
                              paste0(min_group[2,2], " (", round(min_group[2,3], 2), ") "),
                              paste0(max_group[2,2], " (", round(max_group[2,3], 2), ") ")))
    listres[[paste(names(vars[var]))]] <- temp
  }
  bind <-
    do.call(cbind, lapply(listres, function(x) {
    val <- x[["Values"]] %>% as.character()
  })) %>%
    as_tibble() %>%
    mutate(Parameters = listres[[1]][["Parameters"]]) %>%
    select(Parameters, everything())
  return(bind)
}
