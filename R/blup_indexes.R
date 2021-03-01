#' Stability indexes based on a mixed-effect model
#' @description
#' `r badge('stable')`
#'
#' * [hmgv()] Computes the harmonic mean of genotypic values (HMGV).
#' * [rpgv()] Computes the relative performance of the genotypic values (RPGV).
#' * [hmrpgv()] Computes the harmonic mean of the relative performance of
#' genotypic values (HMRPGV).
#' * [blup_indexes()] Is a wrapper around the above functions that also computes
#' the WAASB index (Olivoto et al. 2019) if an object computed with [waasb()] is
#' used as input data.
#' @name blup_indexes
#' @details
#' The indexes computed with this function have been used to select genotypes
#' with stability performance in a mixed-effect model framework. Some examples
#' are in Alves et al (2018), Azevedo Peixoto et al. (2018), Dias et al. (2018)
#' and Colombari Filho et al. (2013).
#'
#'\loadmathjax
#' The HMGV index is computed as
#'\mjsdeqn{HMG{V_i} = \frac{E}{{\sum\limits_{j = 1}^E {\frac{1}{{G{v_{ij}}}}}}}}
#'
#' where \mjseqn{E} is the number of environments included in the analysis,
#' \mjseqn{Gv_{ij}} is the genotypic value (BLUP) for the ith genotype in the
#' jth environment.
#'
#' The RPGV index is computed as
#' \mjsdeqn{RPGV_i = \frac{1}{E}{\sum\limits_{j = 1}^E {Gv_{ij}} /\mathop \mu
#' \nolimits_j }}
#'
#' The HMRPGV index is computed as
#' \mjsdeqn{HMRPG{V_i} = \frac{E}{{\sum\limits_{j = 1}^E
#' {\frac{1}{{G{v_{ij}}/{\mu _j}}}} }}}
#'
#' @param model An object of class `waasb` computed with [waasb()] or
#'   [gamem_met()].
#' @return
#' A tibble containing the indexes.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references
#'   Alves, R.S., L. de Azevedo Peixoto, P.E. Teodoro, L.A. Silva, E.V.
#'   Rodrigues, M.D.V. de Resende, B.G. Laviola, and L.L. Bhering. 2018.
#'   Selection of Jatropha curcas families based on temporal stability and
#'   adaptability of genetic values. Ind. Crops Prod. 119:290-293.
#'   \doi{10.1016/J.INDCROP.2018.04.029}
#'
#'   Azevedo Peixoto, L. de, P.E. Teodoro, L.A. Silva, E.V. Rodrigues, B.G.
#'   Laviola, and L.L. Bhering. 2018. Jatropha half-sib family selection with
#'   high adaptability and genotypic stability. PLoS One 13:e0199880.
#'   \doi{10.1371/journal.pone.0199880}
#'
#'   Colombari Filho, J.M., M.D.V. de Resende, O.P. de Morais, A.P. de Castro,
#'   E.P. Guimaraes, J.A. Pereira, M.M. Utumi, and F. Breseghello. 2013. Upland
#'   rice breeding in Brazil: a simultaneous genotypic evaluation of stability,
#'   adaptability and grain yield. Euphytica 192:117-129.
#'   \doi{10.1007/s10681-013-0922-2}
#'
#'   Dias, P.C., A. Xavier, M.D.V. de Resende, M.H.P. Barbosa, F.A. Biernaski,
#'   R.A. Estopa. 2018. Genetic evaluation of Pinus taeda clones from somatic
#'   embryogenesis and their genotype x environment interaction. Crop Breed.
#'   Appl. Biotechnol. 18:55-64.
#'   \doi{10.1590/1984-70332018v18n1a8}
#'
#'   Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, V.S. Marchioro, V.Q. de
#'   Souza, and E. Jost. 2019. Mean performance and stability in
#'   multi-environment trials I: Combining features of AMMI and BLUP techniques.
#'   Agron. J. 111:2949-2960. \doi{10.2134/agronj2019.03.0220}
#'
#'   Resende MDV (2007) Matematica e estatistica na analise de experimentos e no
#'   melhoramento genetico. Embrapa Florestas, Colombo
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' res_ind <- waasb(data_ge,
#'                  env = ENV,
#'                  gen = GEN,
#'                  rep = REP,
#'                  resp = c(GY, HM))
#' model_indexes <- blup_indexes(res_ind)
#'
#' # Alternatively using the pipe operator %>%
#' res_ind <- data_ge %>%
#'            waasb(ENV, GEN, REP, c(GY, HM)) %>%
#'            blup_indexes()
#'}
#'
#'
# Harmonic mean genotypic values
hmgv <- function(model){
    if (!has_class(model, "waasb")) {
        stop("The object 'model' must be an object of class \"waasb\"")
    }
    hmgv <- function(model){
        GEPRED <-
            model[["BLUPint"]] %>%
            make_mat(GEN, ENV, Predicted)
        HMGV <- tibble(GEN = rownames(GEPRED),
                       Y = model[["MeansGxE"]] %>% means_by(GEN) %>% pull(Y),
                       HMGV = apply(GEPRED, 1, FUN = hmean, na.rm = TRUE),
                       HMGV_R = rank(-HMGV))
        return(HMGV)
    }
    map(model, hmgv) %>%
        set_class("blup_ind") %>%
        return()
}
# Relative performance of genotypic values
#' @name blup_indexes
#' @export
rpgv <- function(model){
    if (!has_class(model, "waasb")) {
        stop("The object 'model' must be an object of class \"waasb\"")
    }
    rpgv <- function(model){
        GEPRED <-
            model[["BLUPint"]] %>%
            make_mat(GEN, ENV, Predicted)
        Y <- model[["MeansGxE"]]
        GEMEAN <- make_mat(Y, GEN, ENV, Y)
        ovmean <- mean(Y$Y)
        mean_env <- apply(GEMEAN, 2, FUN = mean)
        RPGV <-
            tibble(GEN = rownames(GEPRED),
                   Y = model[["MeansGxE"]] %>% means_by(GEN) %>% pull(Y),
                   RPGV = apply(t(t(GEPRED)/mean_env), 1, mean, na.rm = TRUE)) %>%
            add_cols(RPGV_Y = RPGV * ovmean,
                     RPGV_R = rank(-RPGV_Y))
        return(RPGV)
    }
    map(model, rpgv) %>%
        set_class("blup_ind") %>%
        return()
}
# Harmonic mean of the relative performance of genotypic values
#' @name blup_indexes
#' @export
hmrpgv <- function(model){
    if (!has_class(model, "waasb")) {
        stop("The object 'model' must be an object of class \"waasb\"")
    }
    hmrpgv <- function(model){
        GEPRED <-
            model[["BLUPint"]] %>%
            make_mat(GEN, ENV, Predicted)
        Y <- model[["MeansGxE"]]
        GEMEAN <- make_mat(Y, GEN, ENV, Y)
        ovmean <- mean(Y$Y)
        mean_env <- apply(GEMEAN, 2, FUN = mean)
        HMRPGV <-
            tibble(GEN = rownames(GEPRED),
                   Y = Y %>% means_by(GEN) %>% pull(Y),
                   HMRPGV = apply(t(t(GEPRED)/mean_env), 1, hmean, na.rm = TRUE)) %>%
            mutate(HMRPGV_Y = HMRPGV * ovmean,
                   HMRPGV_R = rank(-HMRPGV_Y))
        return(HMRPGV)
    }
    map(model, hmrpgv) %>%
        set_class("blup_ind") %>%
        return()
}
#' @name blup_indexes
#' @export
Resende_indexes <- function(model) {
    if (!is(model, "waasb")) {
        stop("The object 'model' must be an object of class \"waasb\"")
    }
    deprecated_error("1.13.0", "metan::Resende_indexes()", "metan::blup_indexes()")
}
#' @name blup_indexes
#' @export
blup_indexes <- function(model) {
    if (!is(model, "waasb")) {
        stop("The object 'model' must be an object of class \"waasb\"")
    }
    test <- "model" %in% names(model[[1]])
    hmgv <- hmgv(model) %>% rbind_fill_id(.id = "TRAIT")
    rpgv <- rpgv(model) %>%  rbind_fill_id(.id = "TRAIT")
    hmrpgv <- hmrpgv(model) %>%  rbind_fill_id(.id = "TRAIT")
    bind <-
        hmgv %>%
        left_join(rpgv, by = c("TRAIT", "GEN", "Y")) %>%
        left_join(hmrpgv, by = c("TRAIT", "GEN", "Y"))
    if(test){
        waasb_ind <-
            gmd(model, "WAASB", verbose = FALSE) %>%
            pivot_longer(-GEN, names_to = "TRAIT") %>%
            rename(WAASB = value)
        waasb_ran <-
            gmd(model, "OrWAASB", verbose = FALSE) %>%
            pivot_longer(-GEN, names_to = "TRAIT") %>%
            rename(WAASB_R = value)
        bind <-
            bind %>%
            left_join(waasb_ind, by = c("TRAIT", "GEN")) %>%
            left_join(waasb_ran, by = c("TRAIT", "GEN"))
    }
    if(!test){
        warning("The WAASB index was not computed.\nUse an object computed with `waasb()` to get this index.", call. = FALSE)
    }
    bind <-
        split_factors(bind, TRAIT) %>%
        set_class("blup_ind")
    return(bind)
}
