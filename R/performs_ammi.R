#' Additive Main effects and Multiplicative Interaction
#'
#' Compute the Additive Main effects and Multiplicative interaction. This is a
#' helper function for other procedures performed in the \pkg{metan} package
#' such as \code{\link{waas}}, \code{\link{cv_ammi}}, and
#' \code{\link{cv_ammif}}.
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments
#' @param gen The name of the column that contains the levels of the genotypes
#' @param rep The name of the column that contains the levels of the
#'   replications/blocks
#' @param resp The response variable
#' @return
#' * \strong{ANOVA} The analysis of variance for the AMMI model.
#'
#' * \strong{analysis} The principal component analysis
#'
#' * \strong{means} means of genotype vs environment
#'
#' *\strong{biplot} scores for genotypes and environments in all the possible
#' axes.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'
#' library(metan)
#' ammi_model = performs_ammi(data_ge, ENV, GEN, REP, GY)
#'
#'
performs_ammi <- function(.data, env, gen, rep, resp) {
    data <- .data %>% select(!!enquo(env), !!enquo(gen), !!enquo(rep),
                             !!enquo(resp)) %>% as.data.frame()
    names(data) <- c("ENV", "GEN", "REP", "Y")
    GEN <- factor(data[, 2])
    ENV <- factor(data[, 1])
    REP <- factor(data[, 3])
    Y <- eval(data[, 4])
    nenv <- length(unique(ENV))
    ngen <- length(unique(GEN))
    minimo <- min(ngen, nenv) - 1
    nrep <- length(unique(REP))
    model <- aov(Y ~ ENV + REP %in% ENV + GEN + ENV:GEN)
    df <- fortify(model)
    datares <- model$model
    datares$factors <- paste(datares$ENV, datares$GEN)
    residuals <- cbind(datares, df %>% select(fitted = .fitted,
                                              resid = .resid, stdres = .stdresid))
    if (minimo < 2) {
        stop("The analysis AMMI is not possible. Both genotypes and environments must have more than two levels.")
    }
    mm <- anova(model)
    nn <- mm[2, ]
    mm[2, ] <- mm[3, ]
    mm[3, ] <- nn
    row.names(mm)[2] <- "REP(ENV)"
    row.names(mm)[3] <- "GEN     "
    mm[1, 4] <- mm[1, 3]/mm[2, 3]
    mm[1, 5] <- 1 - pf(mm[1, 4], mm[1, 1], mm[2, 1])
    anova <- mm
    probint <- anova[4, 5]
    DFE <- df.residual(model)
    MSE <- deviance(model)/DFE
    MEANS <- data %>% group_by(ENV, GEN) %>% summarise(Y = mean(Y)) %>%
        ungroup()
    residual <- residuals(lm(Y ~ ENV + GEN, data = MEANS))
    MEANS %<>% mutate(RESIDUAL = residual)
    s <- svd(t(matrix(residual, nenv, byrow = T)))
    U <- s$u[, 1:minimo]
    LL <- diag(s$d[1:minimo])
    V <- s$v[, 1:minimo]
    SS <- (s$d[1:minimo]^2) * nrep
    SUMA <- sum(SS)
    percent <- (1/SUMA) * SS * 100
    DFAMMI <- replicate(minimo, 0)
    acum <- DFAMMI
    MSAMMI <- DFAMMI
    F.AMMI <- DFAMMI
    PROBF <- DFAMMI
    acumula <- 0
    for (i in 1:(minimo)) {
        DF <- (ngen - 1) + (nenv - 1) - (2 * i - 1)
        if (DF <= 0)
            break
        DFAMMI[i] <- DF
        acumula <- acumula + percent[i]
        acum[i] <- acum[i] + acumula
        MSAMMI[i] <- SS[i]/DFAMMI[i]
        if (MSE > 0)
            F.AMMI[i] <- round(MSAMMI[i]/MSE, 2) else F.AMMI[i] <- NA
        if (DFE > 0)
            PROBF[i] <- round(1 - pf(F.AMMI[i], DFAMMI[i], DFE),
                              4) else PROBF[i] <- NA
    }
    percent <- round(percent, 1)
    acum <- round(acum, 1)
    SS <- round(SS, 5)
    MSAMMI <- round(MSAMMI, 5)
    SSAMMI <- data.frame(percent, acum, Df = DFAMMI, `Sum Sq` = SS,
                         `Mean Sq` = MSAMMI, `F value` = F.AMMI, Pr.F = PROBF)
    nssammi <- nrow(SSAMMI)
    SSAMMI <- SSAMMI[SSAMMI$Df > 0, ]
    nss <- nrow(SSAMMI)
    row.names(SSAMMI) <- paste("PC", 1:nss, sep = "")
    SCOREG <- U %*% LL^0.5
    SCOREE <- V %*% LL^0.5
    colnames(SCOREG) <- colnames(SCOREE) <- paste("PC", 1:minimo,
                                                  sep = "")
    bplot <- MEANS %>% group_by(GEN) %>% summarise(Y = mean(Y)) %>%
        mutate(type = "GEN") %>% rename(Code = GEN) %>% cbind(.,
                                                              SCOREG) %>% rbind(., MEANS %>% group_by(ENV) %>% summarise(Y = mean(Y)) %>%
                                                                                    mutate(type = "ENV") %>% rename(Code = ENV) %>% cbind(.,
                                                                                                                                          SCOREE)) %>% select(type, Code, everything())
    PC <- SSAMMI %>% dplyr::select(-percent, -acum, everything())
    resid <- as.data.frame(anova[nrow(anova), ])
    rownames(resid) <- "Residuals"
    sum <- as.data.frame(anova[nrow(anova), ])
    sum$Df <- sum(anova$Df)
    sum$`Sum Sq` <- sum(anova$`Sum Sq`)
    sum$`Mean Sq` <- sum$`Sum Sq`/sum$Df
    rownames(sum) <- "Total"
    ERRO <- rbind(resid, sum)
    names(PC) <- paste(c("Df", "Sum Sq", "Mean Sq", "F value",
                         "Pr(>F)", "Percent", "Accumul"))
    anova <- rbind_fill(mm[-nrow(mm), ], PC, ERRO)
    invisible(structure(list(ANOVA = anova, analysis = PC, means = MEANS,
                             biplot = bplot, residuals = residuals, probint = probint),
                        class = "WAASB"))
}
