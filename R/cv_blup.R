#' Cross-validation procedure
#'
#' Cross-validation for blup prediction.
#'
#' This function provides a cross-validation procedure for mixed models using
#' replicate-based data. By default, complete blocks are randomly selected
#' within each evironment. In each iteraction, the original dataset is split up
#' into two datasets: training and validation data. The "training" set has all
#' combinations (genotype x environment) with the number of replications
#' informed in \code{nrepval}. The "validation" set has the remaining
#' replication. The estimated values are compared with the "validation" data
#' and the Root Means Square Prediction Difference is computed. At the end of
#' boots, a list is returned.
#'
#' @param .data The dataset containing the columns related to Environments,
#' Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#' environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#' replications/blocks.
#' @param resp The response variable.
#' @param nboot The number of resamples to be used in the cross-validation. Defaults to 100.
#' @param nrepval The number of replicates (r) from total number of replicates
#' (R) to be used in the modeling dataset. Only one replicate is used as
#' validating data each step, so, \code{Nrepval} must be equal \code{R-1}
#' @param random The effects of the model assumed to be random. Default is
#' \code{random = "gen"} (genotype and genotype-vs-environment as random
#' effects. Other values allowed are \code{random = "env"} (environment,
#' genotype-vs-environment and block-within-environment random effects) or
#' \code{random = "all"} all effects except the intercept are assumed to be
#' random effects.
#' @param verbose A logical argument to define if a progress bar is shown.
#' Default is \code{TRUE}.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso \code{\link{cv_ammi}}, \code{\link{cv_ammif}}
#' @export
#' @examples
#'
#' \dontrun{
#' library(metan)
#' model = cv_blup(data_ge,
#'                         env = ENV,
#'                         gen = GEN,
#'                         rep = REP,
#'                         resp = GY,
#'                         nboot = 100,
#'                         nrepval = 2)
#'
#' # Alternatively using the pipe operator %>%
#' library(dplyr)
#' model = data_ge %>%
#'         cv_blup(ENV, GEN, REP, GY, 100, 2)
#'
#' }
#'
cv_blup <- function(.data, env, gen, rep, resp, nboot = 100, nrepval, random = "gen", verbose = TRUE) {
    data  = .data %>% select(ENV = !!enquo(env),
                             GEN = !!enquo(gen),
                             REP = !!enquo(rep),
                             Y = !!enquo(resp))
    data = mutate(data, ID = rownames(data))
    Nbloc <- length(unique(data$REP))
    Nenv <- length(unique(data$ENV))

    if (nrepval != Nbloc - 1) {
        stop("The number replications used for validation must be equal to total number of replications -1 (In this case ",
             (Nbloc - 1), ").")
    }

    if (verbose == TRUE) {
        pb <- winProgressBar(title = "the model is being built, please, wait.",
                             min = 1, max = nboot, width = 570)
    }

    RMSPDres <- data.frame(RMSPD = matrix(NA, nboot, 1))
    for (b in 1:nboot) {
        tmp = split_factors(data, ENV, keep_factors = TRUE, verbose = FALSE)
        modeling = do.call(rbind,
                           lapply(tmp, function(x){
                               X2 <- sample(unique(data$REP), nrepval, replace = F)
                               x %>%
                                   dplyr::group_by(!!enquo(gen)) %>%
                                   dplyr::filter(unique(data$REP) %in% c(X2))
                           })
        ) %>% as.data.frame()
        rownames(modeling) <- modeling$ID
        testing <- suppressWarnings(dplyr::anti_join(data, modeling, by = c("ENV",
                                                                            "GEN", "REP", "Y", "ID")))%>%
            arrange(ENV, GEN, REP) %>% as.data.frame()
        MEDIAS <- data.frame(modeling %>% dplyr::group_by(ENV, GEN) %>% dplyr::summarise(Y = mean(Y)))


        if(random == "gen"){
            model <- suppressWarnings(suppressMessages(lme4::lmer(Y ~ REP %in% ENV +
                                                                      ENV +
                                                                      (1 | GEN) +
                                                                      (1 | GEN:ENV), data = modeling)))
        }
        if(random == "env"){
            model <- suppressWarnings(suppressMessages(lme4::lmer(Y ~ GEN +
                                                                      (1 | ENV/REP) +
                                                                      (1 | GEN:ENV), data = modeling)))
        }
        if(random == "all"){
            model <- suppressWarnings(suppressMessages(lme4::lmer(Y ~ 1 +
                                                                      (1 | GEN) +
                                                                      (1 | ENV/REP) +
                                                                      (1 | GEN:ENV), data = modeling)))
        }

        validation <- data.frame(mutate(modeling, pred = predict(model)) %>% dplyr::group_by(ENV,
                                                                                             GEN) %>% dplyr::summarise(pred = mean(pred)))
        validation <- mutate(validation, error = pred - testing$Y)

        RMSPD <- sqrt(sum(validation$error^2)/length(validation$error))
        RMSPDres[, 1][b] <- RMSPD
        if (verbose == TRUE) {
            ProcdAtua <- b
            setWinProgressBar(pb, b, title = paste("Estimating BLUPs for ", ProcdAtua,
                                                   " of ", nboot, " validation datasets", "-", round(b/nboot *
                                                                                                         100, 1), "% Concluded -"))
        }
    }
    if(random == "gen"){
        MOD = "BLUP_g"
    }
    if(random == "env"){
        MOD = "BLUP_e"
    }
    if(random == "all"){
        MOD = "BLUP_ge"
    }
    RMSPDres <- dplyr::mutate(RMSPDres, MODEL = MOD) %>%
        dplyr::select(MODEL, everything())
    RMSPDmean <- RMSPDres %>%
        dplyr::group_by(MODEL) %>%
        dplyr::summarise(mean = mean(RMSPD),
                         sd = sd(RMSPD),
                         se = sd(RMSPD)/sqrt(n()),
                         Q2.5 = quantile(RMSPD, 0.025),
                         Q97.5 = quantile(RMSPD, 0.975))
    if (verbose == TRUE) {
        close(pb)
    }
    return(structure(list(RMSPD = RMSPDres, RMSPDmean = RMSPDmean), class = "cv_blup"))
}
