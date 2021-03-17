#' @name data_simula
#' @title Simulate genotype and genotype-environment data
#' @description
#' `r badge('experimental')`
#'
#' * `g_simula()` simulate replicated genotype data.
#' * `ge_simula()` simulate replicated genotype-environment data.
#'
#' @details The functions simulate genotype or genotype-environment data given a
#'   desired number of genotypes, environments and effects. All effects are
#'   sampled from an uniform distribution. For example, given 10 genotypes, and
#'   `gen_eff = 30`, the genotype effects will be sampled as `runif(10, min =
#'   -30, max = 30)`. Use the argument `seed` to ensure reproducibility. If more
#'   than one trait is used (`nvars > 1`), the effects and seed can be passed as
#'   a numeric vector. Single numeric values will be recycled with a warning
#'   when more than one trait is used.
#'
#' @param ngen The number of genotypes.
#' @param nenv The number of environments.
#' @param nrep The number of replications.
#' @param nvars The number of traits.
#' @param gen_eff The genotype effect.
#' @param env_eff The environment effect
#' @param rep_eff The replication effect
#' @param ge_eff The genotype-environment interaction effect.
#' @param res_eff The residual effect. The effect is sampled from a normal
#'   distribution with zero mean and standard deviation equal to `res_eff`.
#'   Be sure to change `res_eff` when changin the `intercept` scale.
#' @param intercept The intercept.
#' @param seed The seed.
#' @md
#' @importFrom tidyr expand_grid
#' @return A data frame with the simulated traits
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(metan)
#' # Genotype data (5 genotypes and 3 replicates)
#' gen_data <-
#'    g_simula(ngen = 5,
#'             nrep = 3,
#'             seed = 1)
#' gen_data
#' inspect(gen_data, plot = TRUE)
#'
#' aov(V1 ~ GEN + REP, data = gen_data) %>% anova()
#'
#' # Genotype-environment data
#' # 5 genotypes, 3 environments, 4 replicates and 2 traits
#' df <-
#' ge_simula(ngen = 5,
#'           nenv = 3,
#'           nrep = 4,
#'           nvars = 2,
#'           seed = 1)
#' ge_plot(df, ENV, GEN, V1)
#' aov(V1 ~ ENV*GEN + ENV/REP, data = df) %>% anova()
#'
#' # Change genotype effect (trait 1 with fewer differences among genotypes)
#' # Define different intercepts for the two traits
#' df2 <-
#' ge_simula(ngen = 10,
#'           nenv = 3,
#'           nrep = 4,
#'           nvars = 2,
#'           gen_eff = c(1, 50),
#'           intercept = c(80, 1500),
#'           seed = 1)
#' ge_plot(df2, ENV, GEN, V2)
#' }
ge_simula <- function(ngen,
                      nenv,
                      nrep,
                      nvars = 1,
                      gen_eff = 20,
                      env_eff = 15,
                      rep_eff = 5,
                      ge_eff = 10,
                      res_eff = 5,
                      intercept = 100,
                      seed = NULL){
  if(length(gen_eff) != 1 & length(gen_eff) != nvars){
    stop("Argument 'gen_eff' must have length 1 or the same length of 'nvars'.", call. = FALSE)
  }
  if(length(gen_eff) == 1 & length(gen_eff) != nvars){
    warning("'gen_eff = ", gen_eff, "' recycled for all the ", nvars, " traits.", call. = FALSE)
  }

  if(length(env_eff) != 1 & length(env_eff) != nvars ){
    stop("Argument 'env_eff' must have length 1 or the same length of 'nvars'.", call. = FALSE)
  }
  if(length(env_eff) == 1 & length(env_eff) != nvars){
    warning("'env_eff = ", env_eff, "' recycled for all the ", nvars, " traits.", call. = FALSE)
  }

  if(length(rep_eff) != 1 & length(rep_eff) != nvars){
    stop("Argument 'rep_eff' must have length 1 or the same length of 'nvars'.", call. = FALSE)
  }
  if(length(rep_eff) == 1 & length(rep_eff) != nvars){
    warning("'rep_eff = ", rep_eff, "' recycled for all the ", nvars, " traits.", call. = FALSE)
  }

  if(length(ge_eff) != 1 & length(ge_eff) != nvars){
    stop("Argument 'ge_eff' must have length 1 or the same length of 'nvars'.", call. = FALSE)
  }
  if(length(ge_eff) == 1 & length(ge_eff) != nvars){
    warning("'ge_eff = ", ge_eff, "' recycled for all the ", nvars, " traits.", call. = FALSE)
  }
  if(length(res_eff) != 1 & length(res_eff) != nvars){
    stop("Argument 'res_eff' must have length 1 or the same length of 'nvars'.", call. = FALSE)
  }
  if(length(res_eff) == 1 & length(res_eff) != nvars){
    warning("'res_eff = ", res_eff, "' recycled for all the ", nvars, " traits.", call. = FALSE)
  }
  if(length(intercept) != 1 & length(intercept) != nvars){
    stop("Argument 'intercept' must have length 1 or the same length of 'nvars'.", call. = FALSE)
  }
  if(length(intercept) == 1 & length(intercept) != nvars){
    warning("'intercept = ", intercept, "' recycled for all the ", nvars, " traits.", call. = FALSE)
  }
  if(!missing(seed) & length(seed) != 1 & length(seed) != nvars){
    stop("Argument 'seed' must have length 1 or the same length of 'nvars'.", call. = FALSE)
  }
  if(length(seed) == 1 & length(seed) != nvars){
    warning("'seed = ", seed, "' recycled for all the ", nvars, " traits.", call. = FALSE)
  }
  compute_one_var <-
    function(ngen = ngen,
             nenv = nenv,
             nrep = nrep,
             gen_eff = gen_eff,
             env_eff = env_eff,
             rep_eff = rep_eff,
             ge_eff = ge_eff,
             res_eff = res_eff,
             intercept = intercept,
             seed = seed){
      if(missing(seed)){
        seed <- .Random.seed[3]
      } else{
        seed <- seed
      }
      df <-
        expand_grid(ENV = paste("E", 1:nenv, sep = ""),
                    GEN = paste("H", 1:ngen, sep = ""),
                    REP = paste("B", 1:nrep, sep = ""))
      set.seed(seed)
      env_eff <-
        data.frame(ENV = paste("E", 1:nenv, sep = ""),
                   env_eff = runif(nenv, -env_eff, env_eff))
      gen_eff <-
        data.frame(GEN = paste("H", 1:ngen, sep = ""),
                   gen_eff = runif(ngen, -gen_eff, gen_eff))
      block_eff <-
        data.frame(REP = paste("B", 1:nrep, sep = ""),
                   block_eff = runif(nrep, -rep_eff, rep_eff))
      df2 <-
        df %>%
        left_join(env_eff, by = "ENV") %>%
        left_join(gen_eff, by = "GEN") %>%
        left_join(block_eff, by = "REP") %>%
        add_cols(ge_eff = rep(runif(nenv * ngen, -ge_eff, ge_eff), each = nrep)) %>%
        add_cols(resid = rnorm(nenv * ngen * nrep, 0, res_eff)) %>%
        add_cols(V1 = intercept) %>%
        mutate(across(V1, ~.x + gen_eff + env_eff + block_eff + ge_eff + resid)) %>%
        remove_cols(env_eff:resid)
      return(df2)
    }
  dfs <- list()
  for(i in 1:nvars){
    temp <- compute_one_var(ngen = ngen,
                            nenv = nenv,
                            nrep = nrep,
                            gen_eff = gen_eff[ifelse(length(gen_eff) != 1, i, 1)],
                            env_eff = env_eff[ifelse(length(env_eff) != 1, i, 1)],
                            rep_eff = rep_eff[ifelse(length(rep_eff) != 1, i, 1)],
                            ge_eff = ge_eff[ifelse(length(ge_eff) != 1, i, 1)],
                            res_eff = res_eff[ifelse(length(res_eff) != 1, i, 1)],
                            intercept = intercept[ifelse(length(intercept) != 1, i, 1)],
                            seed = seed[ifelse(length(seed) != 1, i, 1)])
    dfs[[paste("V", i, sep = "")]] <- temp
  }
  dfs_bind <-
    dfs %>%
    reduce(left_join, by = c("ENV", "GEN", "REP")) %>%
    set_names("ENV", "GEN", "REP", paste("V", 1:nvars, sep = "")) %>%
    as_factor(1:3)
  return(dfs_bind)
}
#' @name data_simula
#' @export
g_simula <- function(ngen,
                     nrep,
                     nvars = 1,
                     gen_eff = 20,
                     rep_eff = 5,
                     res_eff = 5,
                     intercept = 100,
                     seed = NULL){
  if(length(gen_eff) != 1 & length(gen_eff) != nvars){
    stop("Argument 'gen_eff' must have length 1 or the same length of 'nvars'.", call. = FALSE)
  }
  if(length(gen_eff) == 1 & length(gen_eff) != nvars){
    warning("'gen_eff = ", gen_eff, "' recycled for all the ", nvars, " traits.", call. = FALSE)
  }
  if(length(rep_eff) != 1 & length(rep_eff) != nvars){
    stop("Argument 'rep_eff' must have length 1 or the same length of 'nvars'.", call. = FALSE)
  }
  if(length(rep_eff) == 1 & length(rep_eff) != nvars){
    warning("'rep_eff = ", rep_eff, "' recycled for all the ", nvars, " traits.", call. = FALSE)
  }
  if(length(res_eff) != 1 & length(res_eff) != nvars){
    stop("Argument 'res_eff' must have length 1 or the same length of 'nvars'.", call. = FALSE)
  }
  if(length(res_eff) == 1 & length(res_eff) != nvars){
    warning("'res_eff = ", res_eff, "' recycled for all the ", nvars, " traits.", call. = FALSE)
  }
  if(length(intercept) != 1 & length(intercept) != nvars){
    stop("Argument 'intercept' must have length 1 or the same length of 'nvars'.", call. = FALSE)
  }
  if(length(intercept) == 1 & length(intercept) != nvars){
    warning("'intercept = ", intercept, "' recycled for all the ", nvars, " traits.", call. = FALSE)
  }
  if(!missing(seed) & length(seed) != 1 & length(seed) != nvars){
    stop("Argument 'seed' must have length 1 or the same length of 'nvars'.", call. = FALSE)
  }
  if(length(seed) == 1 & length(seed) != nvars){
    warning("'seed = ", seed, "' recycled for all the ", nvars, " traits.", call. = FALSE)
  }
  compute_one_var <-
    function(ngen = ngen,
             nrep = nrep,
             gen_eff = gen_eff,
             rep_eff = rep_eff,
             res_eff = res_eff,
             intercept = intercept,
             seed = seed){
      if(missing(seed)){
        seed <- .Random.seed[3]
      } else{
        seed <- seed
      }
      df <-
        expand_grid(GEN = paste("H", 1:ngen, sep = ""),
                    REP = paste("B", 1:nrep, sep = ""))
      set.seed(seed)
      gen_eff <-
        data.frame(GEN = paste("H", 1:ngen, sep = ""),
                   gen_eff = runif(ngen, -gen_eff, gen_eff))
      block_eff <-
        data.frame(REP = paste("B", 1:nrep, sep = ""),
                   block_eff = runif(nrep, -rep_eff, rep_eff))
      df2 <-
        df %>%
        left_join(gen_eff, by = "GEN") %>%
        left_join(block_eff, by = "REP") %>%
        add_cols(resid = rnorm(ngen * nrep, 0, res_eff)) %>%
        add_cols(V1 = intercept) %>%
        mutate(across(V1, ~.x + gen_eff + block_eff + resid)) %>%
        remove_cols(gen_eff:resid)
      return(df2)
    }
  dfs <- list()
  for(i in 1:nvars){
    temp <- compute_one_var(ngen = ngen,
                            nrep = nrep,
                            gen_eff = gen_eff[ifelse(length(gen_eff) != 1, i, 1)],
                            rep_eff = rep_eff[ifelse(length(rep_eff) != 1, i, 1)],
                            res_eff = res_eff[ifelse(length(res_eff) != 1, i, 1)],
                            intercept = intercept[ifelse(length(intercept) != 1, i, 1)],
                            seed = seed[ifelse(length(seed) != 1, i, 1)])
    dfs[[paste("V", i, sep = "")]] <- temp
  }
  dfs_bind <-
    dfs %>%
    reduce(left_join, by = c("GEN", "REP")) %>%
    set_names("GEN", "REP", paste("V", 1:nvars, sep = "")) %>%
    as_factor(1:2)
  return(dfs_bind)
}
