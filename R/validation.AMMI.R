validation.AMMI <- function(.data, env, gen, rep, resp, design = "RCBD", nboot, nrepval,
                            naxis, verbose = TRUE) {

    if (!design %in% c("RCBD", "CRD")) {
        stop("Incorrect experimental design informed! Plesease inform RCBD for randomized complete block or CRD for completely randomized design.")
    }

    RMSPDres <- data.frame(RMSPD = matrix(0, nboot, 1))
    Y <- eval(substitute(resp), eval(.data))
    GEN <- factor(eval(substitute(gen), eval(.data)))
    ENV <- factor(eval(substitute(env), eval(.data)))
    REP <- factor(eval(substitute(rep), eval(.data)))
    REPS <- eval(substitute(rep), eval(.data))
    data <- data.frame(ENV, GEN, REP, Y)
    data <- mutate(data, ID = rownames(data))
    Nenv <- length(unique(ENV))
    Ngen <- length(unique(GEN))
    Nbloc <- length(unique(REP))
    minimo <- min(Nenv, Ngen) - 1

    if (naxis > minimo) {
        stop("The number of axis to be used must be lesser than or equal to ",
             minimo, " [min(GEN-1;ENV-1)]")
    }

    if (nrepval != Nbloc - 1) {
        stop("The number replications used for validation must be equal to total number of replications -1 (In this case ",
             (Nbloc - 1), ").")
    }

    if (verbose == TRUE) {
        pb <- winProgressBar(title = "the model is being built, please, wait.",
                             min = 1, max = nboot, width = 570)
    }
    condition = (design == "CRD")
    condition2 = (design == "RCBD")

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
            tmp = group_factors(data, !!enquo(env), keep_factors = TRUE, verbose = FALSE)
            modeling = do.call(rbind,
                               lapply(tmp, function(x){
                                   X2 <- sample(unique(REPS), nrepval, replace = F)
                                   x %>%
                                       dplyr::group_by(!!enquo(gen)) %>%
                                       dplyr::filter(REP %in% c(X2))
                               })
            ) %>% as.data.frame()
            rownames(modeling) <- modeling$ID
        }
        testing <- dplyr::anti_join(data, modeling, by = c("ENV", "GEN",
                                                           "REP", "Y", "ID"))
        testing <- testing[order(testing[, 1], testing[, 2], testing[, 3]),
                           ]
        MEDIAS <- data.frame(modeling %>% dplyr::group_by(ENV, GEN) %>%
                                 dplyr::summarise(Y = mean(Y)))
        residual <- residuals(lm(Y ~ ENV + GEN, data = MEDIAS))
        s <- svd(t(matrix(residual, Nenv, byrow = T)))
        MEDIAS <- mutate(MEDIAS,
                         Ypred = Y - residual,
                         ResAMMI =  ((model.matrix(~factor(testing$GEN) - 1) %*% s$u[, 1:naxis]) *
                                         (model.matrix(~factor(testing$ENV) - 1) %*% s$v[, 1:naxis])) %*%
                             s$d[1:naxis],
                         YpredAMMI = Ypred + ResAMMI,
                         testing = testing$Y,
                         error = YpredAMMI - testing,
                         errrorAMMI0 = Ypred - testing)
        if (naxis == 0) {
            RMSPDres[, 1][b] <- sqrt(sum(MEDIAS$errrorAMMI0^2)/length(MEDIAS$errrorAMMI0))
        } else {
            RMSPDres[, 1][b] <- sqrt(sum(MEDIAS$error^2)/length(MEDIAS$error))
        }

        if (verbose == TRUE) {
            ProcdAtua <- b
            setWinProgressBar(pb, b, title = paste("Validating ", ProcdAtua,
                                                   " of ", nboot, "validation datasets, considering", naxis, "axes",
                                                   "-", round(b/nboot * 100, 1), "% Concluded -"))
        }

    }
    if (verbose == TRUE) {
        close(pb)
    }

    RSMEmean <- mean(RMSPDres$RMSPD)
    return(structure(list(RMSPD = RMSPDres, RSMEmean = RSMEmean, Estimated = MEDIAS,
                          Modeling = modeling, Testing = testing), class = "validation.AMMI"))
}
