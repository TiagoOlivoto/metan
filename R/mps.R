#' Mean performance and stability in multi-environment trials
#' @description
#' `r badge('experimental')`
#'
#' This function implements the weighting method between mean performance and
#' stability (Olivoto et al., 2019) considering different parametric and
#' non-parametric stability indexes.
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#'   replications/blocks.
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure a vector of variables may be used. For example `resp
#'   = c(var1, var2, var3)`.
#' @param block Defaults to `NULL`. In this case, a randomized complete
#'   block design is considered. If block is informed, then an alpha-lattice
#'   design is employed considering block as random to make use of inter-block
#'   information, whereas the complete replicate effect is always taken as
#'   fixed, as no inter-replicate information was to be recovered (Mohring et
#'   al., 2015).
#' @param by One variable (factor) to compute the function by. It is a shortcut
#'   to [dplyr::group_by()].This is especially useful, for example, when the
#'   researcher want to analyze environments within mega-environments. In this
#'   case, an object of class mps_grouped is returned.
#' @param random The effects of the model assumed to be random. Defaults to
#'   `random = "gen"`. See [gamem_met()] to see the random effects assumed
#'   depending on the experimental design of the trials.
#' @param verbose Logical argument. If `verbose = FALSE` the code will run
#'   silently.
#' @param performance Wich considers as mean performance. Either `blupg` (for
#'   Best Linear Unbiased Prediction) or `blueg` (for Best Linear Unbiased
#'   Estimation)
#' @param stability The stability method. One of the following:
#' * `"waasb"` The weighted average of absolute scores (Olivoto et al. 2019).
#' * `"ecovalence"` The Wricke's ecovalence (Wricke, 1965).
#' * `"Shukla"` The Shukla's stability variance parameter (Shukla, 1972).
#' * `"hmgv"` The harmonic mean of genotypic values (Resende, 2007).
#' * `"s2di"` The deviations from the Eberhart and Russell regression (Eberhart
#' and Russell, 1966).
#' * `"r2"` The determination coefficient of the Eberhart and Russell regression
#' (Eberhart and Russell, 1966)..
#' * `"rmse"` The root mean squared error of the Eberhart and Russell regression
#' (Eberhart and Russell, 1966).
#' * `"wi"` Annicchiarico's genotypic confidence index (Annicchiarico, 1992).
#' * `"polar"` Power Law Residuals as yield stability index (Doring et al.,
#' 2015).
#' * `"acv"` Adjusted Coefficient of Variation (Doring and Reckling, 2018)
#' * `"pi"` Lin e Binns' superiority index (Lin and Binns, 1988).
#' * `"gai"` Geometric adaptability index (Mohammadi and Amri, 2008).
#' * `"s1", "s2", "s3", and "s6"` Huehn's stability statistics (Huehn, 1979).
#' * `"n1", "n2", "n3", and "n4"` Thennarasu's stability statistics (Thennarasu,
#' 1995).
#' * `"asv", "ev", "za", and "waas"` AMMI-based stability indexes (see
#' [AMMI_indexes()]).
#'
#' @param ideotype_mper,ideotype_stab The new maximum value after rescaling the
#'   response variable/stability index. By default, all variables in `resp` are
#'   rescaled so that de maximum value is 100 and the minimum value is 0 (i.e.,
#'   `ideotype_mper = NULL` and `ideotype_stab = NULL`). It must be a character
#'   vector of the same length of `resp` if rescaling is assumed to be different
#'   across variables, e.g., if for the first variable smaller values are better
#'   and for the second one, higher values are better, then `ideotype_mper =
#'   c("l, h")` must be used. For stability index in which lower values are
#'   better, use `ideotype_stab = "l"`. Character value of length 1 will be
#'   recycled with a warning message.
#' @param wmper The weight for the mean performance. By default, all variables
#'   in `resp` have equal weights for mean performance and stability (i.e.,
#'   `wmper = 50`). It must be a numeric vector of the same length of `resp` to
#'   assign different weights across variables, e.g., if for the first variable
#'   equal weights for mean performance and stability are assumed and for the
#'   second one, a higher weight for mean performance (e.g. 65) is assumed, then
#'   `wmper = c(50, 65)` must be used. Numeric value of length 1 will be
#'   recycled with a warning message.
#' @param verbose Logical argument. If `verbose = FALSE` the code will run
#'   silently.
#' @references
#' Annicchiarico, P. 1992. Cultivar adaptation and recommendation from alfalfa
#' trials in Northern Italy. J. Genet. Breed. 46:269-278.
#'
#' Doring, T.F., S. Knapp, and J.E. Cohen. 2015. Taylor's power law and the
#' stability of crop yields. F. Crop. Res. 183: 294-302.
#' \doi{10.1016/j.fcr.2015.08.005}
#'
#' Doring, T.F., and M. Reckling. 2018. Detecting global trends of cereal yield
#' stability by adjusting the coefficient of variation. Eur. J. Agron. 99:
#' 30-36. \doi{10.1016/j.eja.2018.06.007}
#'
#' Eberhart, S.A., and W.A. Russell. 1966. Stability parameters for comparing Varieties.
#' Crop Sci. 6:36-40. \doi{10.2135/cropsci1966.0011183X000600010011x}
#'
#' Huehn, V.M. 1979. Beitrage zur erfassung der phanotypischen stabilitat. EDV
#' Med. Biol. 10:112.
#'
#' Lin, C.S., and M.R. Binns. 1988. A superiority measure of cultivar
#' performance for cultivar x location data. Can. J. Plant Sci. 68:193-198.
#' \doi{10.4141/cjps88-018}
#'
#' Mohammadi, R., & Amri, A. (2008). Comparison of parametric and non-parametric
#' methods for selecting stable and adapted durum wheat genotypes in variable
#' environments. Euphytica, 159(3), 419-432. \doi{10.1007/s10681-007-9600-6}
#'
#' Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, V.S. Marchioro, V.Q. de
#' Souza, and E. Jost. 2019. Mean performance and stability in multi-environment
#' trials I: Combining features of AMMI and BLUP techniques. Agron. J.
#' \doi{10.2134/agronj2019.03.0220}
#'
#' Resende MDV (2007) Matematica e estatistica na analise de experimentos e no
#' melhoramento genetico. Embrapa Florestas, Colombo
#'
#' Shukla, G.K. 1972. Some statistical aspects of partitioning
#' genotype-environmental components of variability. Heredity. 29:238-245.
#' \doi{10.1038/hdy.1972.87}
#'
#' Thennarasu, K. 1995. On certain nonparametric procedures for studying
#' genotype x environment interactions and yield stability. Ph.D. thesis. P.J.
#' School, IARI, New Delhi, India.
#'
#' Wricke, G. 1965. Zur berechnung der okovalenz bei sommerweizen und hafer. Z.
#' Pflanzenzuchtg 52:127-138.
#'
#' @return An object of class `mps` with the following items.
#' * `observed`: The observed value on a genotype-mean basis.
#' * `performance`: The performance for genotypes (BLUPs or BLUEs)
#' * `performance_res`: The rescaled values of genotype's performance,
#' considering `ideotype_mper`.
#' * `stability`: The stability for genotypes, chosen with argument `stability`.
#' * `stability_res`: The rescaled values of genotype's stability, considering
#' `ideotype_stab`.
#' * `mps_ind`: The mean performance and stability for the traits.
#' * `h2`: The broad-sense heritability for the traits.
#' * `perf_method`: The method for measuring genotype's performance.
#' * `wmper`: The weight for the mean performance.
#' * `sense_mper`: The goal for genotype's performance (`l` = lower, `h` = higher).
#' * `stab_method`: The method for measuring genotype's stability.
#' * `wstab`: The weight for the mean stability.
#' * `sense_stab`: The goal for genotype's stability (`l` = lower, `h` = higher).
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso [mtsi()], [mtmps()], [mgidi()]
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' # The same approach as mtsi()
#' # mean performance and stability for GY and HM
#' # mean performance: The genotype's BLUP
#' # stability: the WAASB index (lower is better)
#' # weights: equal for mean performance and stability
#'
#' model <-
#' mps(data_ge,
#'     env = ENV,
#'     gen = GEN,
#'     rep = REP,
#'     resp = everything())
#'
#' # The mean performance and stability after rescaling
#' model$mps_ind
#'}
mps <- function(.data,
                env,
                gen,
                rep,
                resp,
                block = NULL,
                by = NULL,
                random = "gen",
                performance = "blupg",
                stability = "waasb",
                ideotype_mper = NULL,
                ideotype_stab = NULL,
                wmper = NULL,
                verbose = TRUE) {

  if (!missing(by)){
    if(length(as.list(substitute(by))[-1L]) != 0){
      stop("Only one grouping variable can be used in the argument 'by'.\nUse 'group_by()' to pass '.data' grouped by more than one variable.", call. = FALSE)
    }
    .data <- group_by(.data, {{by}})
  }
  if(is_grouped_df(.data)){
    if(!missing(block)){
      results <-
        .data %>%
        doo(mps,
            env = {{env}},
            gen = {{gen}},
            rep = {{rep}},
            resp = {{resp}},
            block = {{block}},
            random = random,
            performance = performance,
            stability = stability,
            ideotype_mper = ideotype_mper,
            ideotype_stab = ideotype_stab,
            wmper = wmper,
            verbose = verbose)
    } else{
      results <-
        .data %>%
        doo(mps,
            env = {{env}},
            gen = {{gen}},
            rep = {{rep}},
            resp = {{resp}},
            random = random,
            performance = performance,
            stability = stability,
            ideotype_mper = ideotype_mper,
            ideotype_stab = ideotype_stab,
            wmper = wmper,
            verbose = verbose)
    }
    return(set_class(results, c("tbl_df",  "mps_group", "mps", "tbl",  "data.frame")))
  }
  stab <- c("waasb", "ecovalence", "hmgv", "s2di", "r2", "rmse", "wi", "polar", "acv", "pi", "gai", "s1", "s2", "s3", "s6", "n1", "n2", "n3", "n4", "asv", "ev", "za", "waas", "shukla")
  if(!stability %in%  stab){
    stop("`stability` must be one of the ", paste(stab, collapse = ", "), call. = FALSE)
  }
  if(is.numeric(ideotype_mper)){
    stop("Using a numeric vector in `ideotype_mper` is deprecated as of metan 1.9.0. use 'h' or 'l' instead.\nOld code: 'ideotype_mper = c(100, 100, 0)'.\nNew code: 'ideotype_mper = c(\"h, h, l\")", call. = FALSE)
  }

  if(is.numeric(ideotype_stab)){
    stop("Using a numeric vector in `ideotype_stab` is deprecated as of metan 1.9.0. use 'h' or 'l' instead.\nOld code: 'ideotype_stab = c(100, 100, 0)'.\nNew code: 'ideotype_stab = c(\"h, h, l\")", call. = FALSE)
  }
  if(missing(block)){
    mod <-
      gamem_met(.data = .data,
                env = {{env}},
                gen = {{gen}},
                rep = {{rep}},
                resp = {{resp}},
                random = random,
                verbose = verbose)
  } else{
    mod <-
      gamem_met(.data = .data,
                env = {{env}},
                gen = {{gen}},
                rep = {{rep}},
                resp = {{resp}},
                block = {{block}},
                random = random,
                verbose = verbose)
  }
  # mean performance
  observed <-
    gmd(mod, "data", verbose = FALSE) %>%
    means_by(GEN) %>%
    remove_rownames() %>%
    column_to_rownames("GEN")
  mperf <-
    switch(performance,
           blupg =  gmd(mod, "blupg", verbose = FALSE) %>% column_to_rownames("GEN"),
           blueg =  gmd(mod, "blueg", verbose = FALSE) %>% column_to_rownames("GEN")
    )
  nvar <- ncol(mperf)
  if(is.null(ideotype_mper)){
    rescaled <- replicate(ncol(mperf), 100)
    ideotype.D <- replicate(ncol(mperf), 100)
    names(ideotype.D) <- names(data)
  } else{
    rescaled <-
      unlist(strsplit(ideotype_mper, split="\\s*(\\s|,)\\s*")) %>%
      all_lower_case()
    if (length(rescaled) != ncol(mperf)) {
      warning("Invalid length in `ideotype_mper`. Setting `ideotype_mper = '", ideotype_mper[[1]],
              "'` to all the ", nvar, " variables.", call. = FALSE)
      rescaled <- replicate(nvar, ideotype_mper[[1]])
    }
    if(!all(rescaled %in% c("h", "l", "m"))){
      stop("argument `ideotype_mper` must have 'h', 'l', or 'm' only", call. = FALSE)
    }
    ideotype.D <- ifelse(rescaled == "m", 50, 100)
    names(ideotype.D) <- colnames(mperf)
    rescaled <- case_when(
      rescaled == "h" ~ 100,
      rescaled == "l" ~ 0,
      TRUE ~ 100)
  }
  mperf_res <- data.frame(matrix(ncol = ncol(mperf), nrow = nrow(mperf)))
  rownames(mperf_res) <- rownames(mperf)
  colnames(mperf_res) <- colnames(mperf)
  for (i in 1:ncol(mperf)) {
    mperf_res[i] <- resca(values = mperf[i], new_max = rescaled[i], new_min = 100 - rescaled[i])
  }
  sense_mper <- rescaled

  # stability
  stab <-
    switch(stability,
           # Shukla's stability variance parameter
           shukla = Shukla(.data,
                           env = {{env}},
                           gen = {{gen}},
                           rep = {{rep}},
                           resp = {{resp}},
                           verbose = FALSE) %>%
             gmd("ShuklaVar", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           # AMMI-based indexes
           waas = performs_ammi(.data,
                              env = {{env}},
                              gen = {{gen}},
                              rep = {{rep}},
                              resp = {{resp}},
                              verbose = FALSE) %>%
             AMMI_indexes() %>%
             gmd("WAAS", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           za = performs_ammi(.data,
                              env = {{env}},
                              gen = {{gen}},
                              rep = {{rep}},
                              resp = {{resp}},
                              verbose = FALSE) %>%
             AMMI_indexes() %>%
             gmd("ZA", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           ev = performs_ammi(.data,
                                env = {{env}},
                                gen = {{gen}},
                                rep = {{rep}},
                                resp = {{resp}},
                                verbose = FALSE) %>%
             AMMI_indexes() %>%
             gmd("EV", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           sipc = performs_ammi(.data,
                               env = {{env}},
                               gen = {{gen}},
                               rep = {{rep}},
                               resp = {{resp}},
                               verbose = FALSE) %>%
             AMMI_indexes() %>%
             gmd("SIPC", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           asv = performs_ammi(.data,
                               env = {{env}},
                               gen = {{gen}},
                               rep = {{rep}},
                               resp = {{resp}},
                               verbose = FALSE) %>%
             AMMI_indexes() %>%
             gmd("ASV", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           # Thennarasu's stability statistics
           n1 = Thennarasu(.data,
                           env = {{env}},
                           gen = {{gen}},
                           resp = {{resp}},
                           verbose = FALSE) %>%
             gmd("N1", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           n2 = Thennarasu(.data,
                           env = {{env}},
                           gen = {{gen}},
                           resp = {{resp}},
                           verbose = FALSE) %>%
             gmd("N2", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           n3 = Thennarasu(.data,
                           env = {{env}},
                           gen = {{gen}},
                           resp = {{resp}},
                           verbose = FALSE) %>%
             gmd("N3", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           n4 = Thennarasu(.data,
                           env = {{env}},
                           gen = {{gen}},
                           resp = {{resp}},
                           verbose = FALSE) %>%
             gmd("N4", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           # Huehn's stability statistics
           s1 = Huehn(.data,
                      env = {{env}},
                      gen = {{gen}},
                      resp = {{resp}},
                      verbose = FALSE) %>%
             gmd("S1", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           s2 = Huehn(.data,
                      env = {{env}},
                      gen = {{gen}},
                      resp = {{resp}},
                      verbose = FALSE) %>%
             gmd("S2", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           s3 = Huehn(.data,
                      env = {{env}},
                      gen = {{gen}},
                      resp = {{resp}},
                      verbose = FALSE) %>%
             gmd("S3", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           s6 = Huehn(.data,
                      env = {{env}},
                      gen = {{gen}},
                      resp = {{resp}},
                      verbose = FALSE) %>%
             gmd("S6", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           # Geometric adaptability index
           gai = gai(.data,
                     env = {{env}},
                     gen = {{gen}},
                     resp = {{resp}},
                     verbose = FALSE) %>%
             gmd(verbose = FALSE) %>%
             column_to_rownames("GEN"),
           # Lin e Binns' superiority index
           pi = superiority(.data,
                              env = {{env}},
                              gen = {{gen}},
                              resp = {{resp}},
                              verbose = FALSE) %>%
             gmd(verbose = FALSE) %>%
             column_to_rownames("GEN"),
           # Annicchiarico's genotypic confidence index
           wi = Annicchiarico(.data,
                            env = {{env}},
                            gen = {{gen}},
                            rep = {{rep}},
                            resp = {{resp}},
                            verbose = FALSE) %>%
             gmd(verbose = FALSE) %>%
             column_to_rownames("GEN"),
           # Power Law Residuals
           polar = ge_polar(.data,
                        env = {{env}},
                        gen = {{gen}},
                        resp = {{resp}},
                        verbose = FALSE) %>%
             gmd(verbose = FALSE) %>%
             column_to_rownames("GEN"),
           # adjusted coefficient of variation
           acv = ge_acv(.data,
                        env = {{env}},
                        gen = {{gen}},
                        resp = {{resp}},
                        verbose = FALSE) %>%
             gmd(verbose = FALSE) %>%
             column_to_rownames("GEN"),
           # WAASB
           waasb = waasb(.data = .data,
                         env = {{env}},
                         gen = {{gen}},
                         rep = {{rep}},
                         resp = {{resp}},
                         random = random,
                         verbose = FALSE) %>%
             gmd("WAASB", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           #Ecovalence
           ecovalence = ecovalence(.data,
                                   env = {{env}},
                                   gen = {{gen}},
                                   rep = {{rep}},
                                   resp = {{resp}},
                                   verbose = FALSE) %>%
             gmd(verbose = FALSE) %>%
             column_to_rownames("GEN"),
           # hmgv
           hmgv = blup_indexes(mod) %>%
             suppressWarnings() %>%
             gmd("HMGV", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           # deviations from the regression
           s2di = ge_reg(.data,
                         env = {{env}},
                         gen = {{gen}},
                         rep = {{rep}},
                         resp = {{resp}},
                         verbose = FALSE) %>%
             gmd("s2di", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           # RMSE
           rmse = ge_reg(.data,
                         env = {{env}},
                         gen = {{gen}},
                         rep = {{rep}},
                         resp = {{resp}},
                         verbose = FALSE) %>%
             gmd("RMSE", verbose = FALSE) %>%
             column_to_rownames("GEN"),
           # R2
           r2 = ge_reg(.data,
                       env = {{env}},
                       gen = {{gen}},
                       rep = {{rep}},
                       resp = {{resp}},
                       verbose = FALSE) %>%
             gmd("R2", verbose = FALSE) %>%
             column_to_rownames("GEN")
    )

  if(is.null(ideotype_stab)){
    rescaled <- replicate(ncol(stab), 0)
    ideotype.D <- replicate(ncol(stab), 0)
    names(ideotype.D) <- names(data)
  } else{
    rescaled <-
      unlist(strsplit(ideotype_stab, split="\\s*(\\s|,)\\s*")) %>%
      all_lower_case()
    if (length(rescaled) != ncol(mperf)) {
      warning("Invalid length in `ideotype_stab`. Setting `ideotype_stab = '", ideotype_stab[[1]],
              "'` to all the ", nvar, " variables.", call. = FALSE)
      rescaled <- replicate(nvar, ideotype_stab[[1]])
    }
    if(!all(rescaled %in% c("h", "l", "m"))){
      stop("`ideotype_stab` must have 'h', 'l', or 'm' only", call. = FALSE)
    }
    ideotype.D <- ifelse(rescaled == "m", 50, 100)
    names(ideotype.D) <- colnames(stab)
    rescaled <- case_when(
      rescaled == "h" ~ 100,
      rescaled == "l" ~ 0,
      TRUE ~ 100)
  }
  stab_res <- data.frame(matrix(ncol = ncol(stab), nrow = nrow(stab)))
  rownames(stab_res) <- rownames(stab)
  colnames(stab_res) <- colnames(stab)
  for (i in 1:ncol(stab)) {
    stab_res[i] <- resca(values = stab[i], new_max = rescaled[i], new_min = 100 - rescaled[i])
  }

  # weight mean performance and stability

  if (is.null(wmper)) {
    peso_perf <- replicate(nvar, 50)
    peso_stab <- 100 - peso_perf
  } else {
    peso_perf <- wmper
    peso_stab <- 100 - peso_perf
    if (length(wmper) != nvar) {
      warning("Invalid length in `wmper`. Setting `wmper = ", wmper[[1]],
              "` to all the ", nvar, " variables.", call. = FALSE)
      peso_perf <- replicate(nvar, wmper[[1]])
      peso_stab <- 100 - peso_perf
    }
    if (min(wmper) < 0 | max(wmper) > 100) {
      stop("The range of the numeric vector `wmper` must be between 0 and 100.")
    }
  }

  listres <- list()
  for (i in 1:nvar) {
    weghting <-
      tibble(gen = rownames(mperf_res),
             mper = as.numeric(mperf_res[[i]]),
             stab = as.numeric(stab_res[[i]]),
             mps_ind = ((mper * peso_perf[i]) + (stab * peso_stab[i]))/(peso_perf[i] + peso_stab[i])) %>%
      select_cols(gen, mps_ind) %>%
      set_names("GEN", names(mperf_res)[i])
    listres[[paste(names(mperf_res)[i])]] <- weghting
  }
  mps_ind <- listres %>% reduce(left_join, by = "GEN")
  if(verbose == TRUE){
    message("Mean performance: ", performance)
    message("Stability: ", stability)
  }
  return(
    list(
      observed = observed,
      performance = mperf,
      performance_res = mperf_res,
      stability = stab,
      stability_res = stab_res,
      mps_ind = mps_ind,
      h2 = gmd(mod, "h2", verbose = FALSE),
      perf_method = performance,
      wmper = peso_perf,
      sense_mper = ifelse(sense_mper == 0, "l", "h"),
      stab_method = stability,
      wstab = peso_stab,
      sense_stab = ifelse(rescaled == 0, "l", "h")
    ) %>% set_class("mps")
  )
}

