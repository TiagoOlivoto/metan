validation.blup <- function(.data, env, gen, rep, resp, nboot, nrepval, verbose = TRUE) {


    RMSPDres <- data.frame(RMSPD = matrix(".", nboot, 1))
    for (n in c(1, 1:ncol(RMSPDres))) {
        RMSPDres[, n] <- as.numeric(RMSPDres[, n])
    }
    Y <- eval(substitute(resp), eval(.data))
    GEN <- factor(eval(substitute(gen), eval(.data)))
    ENV <- factor(eval(substitute(env), eval(.data)))
    REP <- factor(eval(substitute(rep), eval(.data)))
    REPS <- eval(substitute(rep), eval(.data))
    data <- data.frame(ENV, GEN, REP, Y)
    data <- mutate(data, ID = rownames(data))
    Nbloc <- length(unique(REP))
    Nenv <- length(unique(ENV))

    if (nrepval != Nbloc - 1) {
        stop("The number replications used for validation must be equal to total number of replications -1 (In this case ",
            (Nbloc - 1), ").")
    } else {

        if (verbose == TRUE) {
            pb <- winProgressBar(title = "the model is being built, please, wait.",
                min = 1, max = nboot, width = 570)
        }
        for (b in 1:nboot) {
            temp <- data.frame(matrix(".", 0, ncol(data)))
            for (n in c(4:5)) {
                temp[, n] <- as.numeric(temp[, n])
            }
            names(temp) <- names(data)
            actualenv <- 0
            for (K in 1:Nenv) {
                X <- sample(1:10000, 1)
                set.seed(X)
                X2 <- sample(unique(REPS), nrepval, replace = F)
                names <- factor(data$ENV)
                names <- levels(names)[actualenv + 1]
                actualenv <- actualenv + 1
                temp2 <- dplyr::filter(data, ENV == names)
                modeling <- temp2 %>% dplyr::group_by(GEN) %>% dplyr::filter(REP %in%
                  c(X2))
                modeling <- as.data.frame(modeling)
                modeling <- rbind(temp, modeling)
                temp <- modeling
            }
            rownames(modeling) <- modeling$ID
            testing <- suppressWarnings(dplyr::anti_join(data, modeling, by = c("ENV",
                "GEN", "REP", "Y", "ID")))
            testing <- testing[order(testing[, 1], testing[, 2], testing[, 3]), ]
            MEDIAS <- data.frame(modeling %>% dplyr::group_by(ENV, GEN) %>% dplyr::summarise(Y = mean(Y)))

            model <- suppressWarnings(suppressMessages(lme4::lmer(Y ~ REP %in% ENV +
                (1 | GEN) + ENV + (1 | GEN:ENV), data = modeling)))
            validation <- data.frame(mutate(modeling, pred = predict(model)) %>% dplyr::group_by(ENV,
                GEN) %>% dplyr::summarise(pred = mean(pred)))
            validation <- mutate(validation, error = pred - testing$Y)

            RMSPD <- sqrt(sum(validation$error^2)/length(validation$error))
            RMSPDres[, 1][b] <- RMSPD
            if (verbose == TRUE) {
                ProcdAtua <- b
                setWinProgressBar(pb, b, title = paste("Estimating BLUPs for ", ProcdAtua,
                  " of ", nboot, " total validation datasets", "-", round(b/nboot *
                    100, 1), "% Concluded -"))
            }
        }
        RMSPDres <- dplyr::mutate(RMSPDres, MODEL = "BLUP") %>% dplyr::select(MODEL,
            everything())
        RMSPDmean <- RMSPDres %>% dplyr::group_by(MODEL) %>% dplyr::summarise(mean = mean(RMSPD))
        if (verbose == TRUE) {
            close(pb)
        }
        return(structure(list(RMSPD = RMSPDres, RMSPDmean = RMSPDmean), class = "validation.blup"))
    }
}

