#' @title Utilities for data organization
#'
#' @description
#' `r badge('experimental')`
#'
#' Useful function for data organization before statistical analysis
#' * `add_seq_block()`: Add a column with sequential block numeration in
#' multi-environment data sets.
#' * `recode_factor()`: Recode a factor column. A sequential numbering (with
#' possible prefix) is used to identify each level.
#' * `df_to_selegen_54()`: Given a multi-environment data with environment,
#' genotype, and replication, format the data to be used in the Selegen software
#' (model 54).
#' @name utils_data_org
#' @param data A data frame.
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes,
#'   that will be treated as random effect.
#' @param rep The name of the column that contains the levels of the
#'   replications/blocks.
#' @param new_factor The name of the new column created.
#' @param prefix An optional prefix to bind with the new factor.
#' @param factor A column to recode.
#' @param verbose Logical argument. If `verbose = FALSE` the code will run
#'   silently.
#' @md
#' @importFrom dplyr distinct
#' @references
#' Resende, M.D. V. 2016. Software Selegen-REML/BLUP: a useful tool for plant
#' breeding. Crop Breed. Appl. Biotechnol. 16(4): 330–339.
#' \doi{10.1590/1984-70332016v16n4a49}.
#'
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(metan)
#' df_ge <- ge_simula(ngen = 2,
#'                    nenv = 3,
#'                    nrep = 2) %>%
#'          add_cols(ENV = c(rep("CACIQUE", 4),
#'                           rep("FREDERICO", 4),
#'                           rep("SANTA_MARIA", 4)))
#' df_ge
#'
#' # Add sequential block numbering over environments
#' add_seq_block(df_ge, ENV, REP, prefix = "B")
#'
#' # Recode the 'ENV' column to "ENV1", "ENV2", and so on.
#' recode_factor(df_ge,
#'               factor = ENV,
#'               prefix = "ENV",
#'               new_factor = ENV_CODE)
#'
#' # Format the data to be used in the Selegen software (model 54)
#' df <- df_to_selegen_54(df_ge, ENV, GEN, REP) %>%
#' recode_factor(ENV, prefix = "E", new_factor = ENV)
#' }
#'
add_seq_block <- function(data,
                          env,
                          rep,
                          new_factor = BLOCK,
                          prefix = "",
                          verbose = TRUE){
  call <- match.call()
  df <- data %>% arrange({{env}}, {{rep}})
  temp <-
    df %>%
    arrange({{env}}, {{rep}}) %>%
    group_by({{env}}) %>%
    count({{rep}}) %>%
    ungroup()
  temp <-
    mutate(temp, temp = 1:nrow(temp)) %>%
    select(n, temp) %>%
    as.matrix()
  block_vector <- c()
  for(i in 1:nrow(temp)){
    block_vector <- append(block_vector, replicate(temp[i, 1], temp[i, 2]))
  }
  df %<>%
    mutate({{new_factor}} := paste(prefix, block_vector, sep = ""),
           .after = {{rep}})
  if(verbose == TRUE){
    message("The data `", call[["data"]], "` has been arranged according to the `",
            call[["env"]], "` and `", call[["rep"]], "` columns.")
  }
  return(df)
}

#' @name utils_data_org
#' @export
recode_factor <- function(data,
                          factor,
                          new_factor = CODE,
                          prefix = "",
                          verbose = TRUE){
  call <- match.call()
  a <- distinct(data, {{factor}}) %>% pull() %>% as.character()
  temp <-
    data %>%
    group_by({{factor}}) %>%
    count({{factor}}) %>%
    ungroup() %>%
    as.data.frame()
  temp <-
    temp %>%
    mutate(code = 1:nrow(temp))
  row.names(temp) <- temp[,1]
  temp[,1] <- NULL
  block_vector <- c()
  for(i in 1:nrow(temp)){
    block_vector <- append(block_vector, replicate(temp[i, 1], temp[i, 2]))
  }
  data %<>%
    arrange({{factor}}) %>%
    mutate({{new_factor}} := paste(prefix, block_vector, sep = ""),
           .after = {{factor}})
  if(verbose == TRUE){
    message("The data `", call[["data"]], "` has been arranged according to the `",
            call[["factor"]], "` column.")
  }
  return(data)
}

#' @name utils_data_org
#' @export
df_to_selegen_54 <- function(data,
                             env,
                             gen,
                             rep,
                             verbose = TRUE){
  call <- match.call()
  if(verbose == TRUE){
    message("The data `", call[["data"]], "` has been arranged according to the `",
            call[["env"]], "` and `", call[["rep"]], "` columns.")
  }
  add_seq_block(data,
                env = {{env}},
                rep = {{rep}},
                new_factor = REP,
                verbose = FALSE) %>%
    add_cols(Parcela = 1:nrow(data),
             repeticao = concatenate(., {{env}}, {{gen}}, pull = TRUE),
             obs = 1) %>%
    select({{env}}, Parcela, {{gen}}, REP, repeticao, obs, everything()) %>%
    colnames_to_upper()
}

