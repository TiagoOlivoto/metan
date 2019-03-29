AMMI_indexes <- function(.data, order.y = NULL) {

    if (is.null(order.y)) {
        order.y <- rep("h", length(.data))
    }
    if (!is(.data, "WAAS.AMMI")) {
        stop("The object 'x' must be an object of class \"WAAS.AMMI\"")
    }
    if (any(!order.y %in% c("h", "l")) == TRUE) {
        stop("The arguments in 'order.y' must be one of the 'h' or 'l'")
    }
    if (length(order.y) != length(.data)) {
        stop("The lenght of argument 'order.y' must be equal the length of 'x'")
    }
    listres <- list()
    varin <- 1
    for (var in 1:length(.data)) {

        model <- .data[[var]]

        n <- sum(model$PCA$`Pr(>F)` <= 0.05, na.rm = TRUE)
        meange <- model$MeansGxE
        effects <- residuals(lm(Y ~ ENV + GEN, data = meange))
        meange$residual <- effects
        ge <- by(meange[, 7], meange[, c(2, 1)], function(x) sum(x, na.rm = TRUE))
        ge <- array(ge, dim(ge), dimnames(ge))
        svdge <- svd(ge)
        gamma.n <- svdge$u[, 1:n]
        theta.n <- model$PCA$Percent[1:n]/100

        # sum of absolute scores
        PCA <- data.frame(model$model)
        SCOR <- PCA[PCA[, 1] == "GEN", c(seq(4, n + 3))]
        if (n == 1) {
            SIPC <- abs(SCOR)
        } else {
            SIPC <- unname(rowSums(apply(SCOR, 2, FUN = abs)))
        }
        mean <- PCA[PCA[, 1] == "GEN", c(2:3)]
        rS <- rank(SIPC)

        if (length(order.y) == 1) {
            if (order.y == "h") {
                rY <- rank(-mean[2])
            }
            if (order.y == "l") {
                rY <- rank(mean[2])
            }
        } else {
            if (order.y[[varin]] == "h") {
                rY <- rank(-mean[2])
            }
            if (order.y[[varin]] == "l") {
                rY <- rank(mean[2])
            }
        }
        varin <- varin + 1
        ssiSIPC <- rS + rY


        # Za
        if (n == 1) {
            Za <- abs(gamma.n * theta.n)
        } else {
            Za <- rowSums(abs(gamma.n %*% diag(theta.n)))
        }
        rZA <- rank(Za)
        ssiZA <- rZA + rY


        # averages of the squared eigenvector
        if (n == 1) {
            EV <- gamma.n^2/n
        } else {
            EV <- rowSums(gamma.n^2/n)
        }
        rEV <- rank(EV)
        ssiEV <- rEV + rY

        # AMMI stability values
        pc <- model$anova$`Sum Sq`[5]/model$anova$`Sum Sq`[6]
        SCOR2 <- PCA[PCA[, 1] == "GEN", c(seq(4, 2 + 3))]
        ASV <- sqrt((pc * SCOR2[, 1])^2 + SCOR2[, 2]^2)
        rASV <- rank(ASV)
        ssiASV <- rASV + rY

        temp <- data.frame(GEN = mean[1], Y = mean[2], rY = rY, ASV = ASV, rASV = rASV,
            ssiASV = ssiASV, SIPC = SIPC, rSIPC = rS, ssiSIPC = ssiSIPC, EV = EV,
            rEV = rEV, ssiEV = ssiEV)

        listres[[paste(names(.data[var]))]] <- temp

    }

    invisible(listres)
}

