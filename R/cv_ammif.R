#' Cross-validation for estimation of all AMMI-family models
#'
#' Cross-validation for estimation of all AMMI-family models
#'
#' This function provides a complete cross-validation of replicate-based data
#' using AMMI-family models. By default, the first validation is carried out
#' considering the AMMIF (all possible axis used). Considering this model, the
#' original dataset is split up into two datasets: training set and validation
#' set. The "training" set has all combinations (genotype x environment) with
#' the number of replications informed in \code{nrepval}. The dataset
#' "validation" set has the remaining replication. The splitting of the dataset
#' into modeling and validating data depends on the design informed. For
#' Completely Randomized Block Design (default), completely blocks are selected
#' within environments. The remained block serves validation data. If
#' \code{design = "RCD"} is informed, completely randomly samples are made for
#' each genotype-by-environment combination. The estimated values (depending on
#' the \code{naxis} informed) are compared with the "validation" data. the Root
#' Mean Square Prediction Difference (RMSPD) is computed. At the end of boots,
#' a list is returned.
#'
#' @param .data The dataset containing the columns related to Environments,
#' Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#' environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#' replications/blocks.
#' @param resp The response variable.
#' @param nboot The number of resamples to be used in the cross-validation
#' @param design The experimental desig to be considered. Default is
#' \code{RCBD} (Randomized complete Block Design). For Completely Randomized
#' Designs inform \code{design = "CRD"}.
#' @param nrepval The number of replicates (r) from total number of replicates
#' (R) to be used in the modeling dataset. Only one replicate is used as
#' validating data each step, so, \code{Nrepval} must be equal \code{R-1}
#' @param verbose A logical argument to define if a progress bar is shown.
#' Default is \code{TRUE}.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso \code{\link{cv_ammi}, \link{cv_blup}}
#' @export
#' @examples
#'
#' \dontrun{
#' library(metan)
#' model = cv_ammif(data_ge,
#'                          env = ENV,
#'                          gen = GEN,
#'                          rep = REP,
#'                          resp = GY,
#'                          nboot = 100,
#'                          nrepval = 2)
#'
#' # Alternatively (and more intuitively) using the pipe operator %>%
#' library(dplyr)
#' model = data_ge %>%
#'         cv_ammif(ENV, GEN, REP, GY, 100, 2)
#' }
#'
cv_ammif <- function(.data, env, gen, rep, resp, nboot = 50,
                     design = "RCBD", nrepval, verbose = TRUE) {
    if (!design %in% c("RCBD", "CRD")) {
        stop("Incorrect experimental design informed! Plesease inform RCBD for randomized complete block or CRD for completely randomized design.")
    }
    data  = .data %>% select(ENV = !!enquo(env),
                             GEN = !!enquo(gen),
                             REP = !!enquo(rep),
                             Y = !!enquo(resp))
    data = mutate(data, ID = rownames(data))
    Nenv <- length(unique(data$ENV))
    Ngen <- length(unique(data$GEN))
    Nbloc <- length(unique(data$REP))
    minimo <- min(Nenv, Ngen) - 1
    naxisvalidation <- minimo + 1
    totalboot <- naxisvalidation * nboot
    initial <- 0
    NAXIS <- minimo
    AMMIval <- data.frame(matrix(NA, nboot, 1))


    if (nrepval != Nbloc - 1) {
        stop("The number replications used for validation must be equal to total number of replications -1 (In this case ",
             (Nbloc - 1), ").")
    }

    if (verbose == TRUE) {
        pb <- winProgressBar(title = "the model is being built, please, wait.",
                             min = 1, max = totalboot, width = 570)
    }
    condition = (design == "CRD")
    condition2 = (design == "RCBD")
    for (y in 1:naxisvalidation) {
        RMSPDres <- data.frame(matrix(NA, nboot, 1))

        for (b in 1:nboot) {
            if (condition) {
                X <- sample(1:10000, 1)
                set.seed(X)
                modeling <- data %>% dplyr::group_by(ENV, GEN) %>% dplyr::sample_n(nrepval,
                                                                                   replace = F)
                modeling <- as.data.frame(modeling[order(modeling$ID), ])
                rownames(modeling) <- modeling$ID
            }
            if (condition2) {
                tmp = split_factors(data, !!enquo(env), keep_factors = TRUE, verbose = FALSE)
                modeling = do.call(rbind,
                                   lapply(tmp, function(x){
                                       X2 <- sample(unique(data$REP), nrepval, replace = F)
                                       x %>%
                                           dplyr::group_by(!!enquo(gen)) %>%
                                           dplyr::filter(unique(data$REP) %in% c(X2))
                                   })
                ) %>% as.data.frame()
                rownames(modeling) <- modeling$ID
            }

            testing <- suppressWarnings(dplyr::anti_join(data, modeling, by = c("ENV", "GEN", "REP", "Y", "ID"))) %>%
                arrange(ENV, GEN, REP) %>% as.data.frame()
            MEDIAS <- data.frame(modeling %>% dplyr::group_by(ENV, GEN) %>%
                                     dplyr::summarise(Y = mean(Y)))
            residual <- residuals(lm(Y ~ ENV + GEN, data = MEDIAS))
            s <- svd(t(matrix(residual, Nenv, byrow = T)))
            MGEN = model.matrix(~factor(testing$GEN) - 1)
            MENV = model.matrix(~factor(testing$ENV) - 1)
            MEDIAS %<>% mutate(Ypred = Y - residual,
                             ResAMMI =  ((MGEN %*% s$u[, 1:NAXIS]) *
                                             (MENV %*% s$v[, 1:NAXIS])) %*%
                                 s$d[1:NAXIS],
                             YpredAMMI = Ypred + ResAMMI,
                             testing = testing$Y,
                             error = YpredAMMI - testing,
                             errrorAMMI0 = Ypred - testing)

            if (NAXIS == 0) {
                RMSPDres[, 1][b] <- sqrt(sum(MEDIAS$errrorAMMI0^2)/length(MEDIAS$errrorAMMI0))
            } else {
                RMSPDres[, 1][b] <- sqrt(sum(MEDIAS$error^2)/length(MEDIAS$error))
            }

            if (NAXIS == minimo) {
                ACTUAL <- "AMMIF"
            } else ACTUAL <- sprintf("AMMI%.0f", NAXIS)
            initial <- initial + 1
            if (verbose == TRUE) {
                ProcdAtua <- b
                setWinProgressBar(pb, initial, title = paste("|Family = ", ACTUAL,
                                                             "| Validating ", ProcdAtua, " of ", nboot, "validation datasets, considering",
                                                             NAXIS, "axes", "-", round(initial/totalboot * 100, 1), "% Concluded -"))
            }
        }

        if (NAXIS == minimo) {
            AMMIval[["AMMIF"]] <- RMSPDres[, 1]
        } else AMMIval[[sprintf("AMMI%.0f", NAXIS)]] <- RMSPDres[, 1]
        NAXIS <- NAXIS - 1
        initial <- initial
    }
    if (verbose == TRUE) {
        close(pb)
    }

    RMSPD <- stack(AMMIval[, -c(1)]) %>% dplyr::select(ind, everything())
    names(RMSPD) = c("MODEL", "RMSPD")
    RMSPDmean <- RMSPD %>% dplyr::group_by(MODEL) %>% dplyr::summarise(mean = mean(RMSPD))
    RMSPDmean <- RMSPDmean[order(RMSPDmean$mean), ]
    return(structure(list(RMSPD = RMSPD, RMSPDmean = RMSPDmean, Estimated = MEDIAS,
                          Modeling = modeling, Testing = testing), class = "cv_ammif"))
}
