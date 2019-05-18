## WAASratio.AMMI.R
#' Different scenarios of stability and mean performance
#'
#' This function computes the WAASBY index in AMMI analysis in
#' different combinations of weights for stability and productivity.
#'
#' This function is very similar to the \code{WAASBYratio}. The main difference
#' is that here, the WAASBY values are computed considering a traditional AMMI
#' model.
#'
#' @param .data The dataset containing the columns related to Environments,
#' Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#' environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#' replications/blocks.
#' @param resp The response variable.
#' @param mresp A numeric value that will be the new maximum value after reescaling.
#' By default, the variable in \code{resp} is rescaled so that the original maximum
#'  and minimum values are 100 and 0, respectively.
#' @param p.valuePC The p-value for considering the PC significant. Default is
#' 0.05. The number of significant Principal Components to be used for
#' calculating the WAASB will be chosen based on this probability.
#' @param increment The range of the increment for WAAS/GY ratio. Default is 5.
#' The function compute the WAASY values starting with a weight o 100 for
#' stability and 0 for response variable. With the default, the first scenario
#' will be a WAAS/GY ratio = 100/0. In the next scenario, the WAASY values are
#' computed based on a WAAS/GY ratio = 95/5.
#' @param saveWAASY Automatically save the WAASY values when the wheight for
#' WAAS (stability) in the WAAS/GY ratio is "saveWAASY". Default is 100. The
#' value of "saveWAASY" must be multiple of "Increment". If this assumption is
#' not valid, an error will be occour.
#' @param progbar A logical argument to define if a progress bar is shown.
#' Default is \code{TRUE}.
#' @return \item{anova}{Joint analysis of variance for the main effects and
#' Principal Component analysis of the interaction effect.}
#'
#' \item{PC}{Principal Component Analysis.}
#'
#' \item{MeansGxE}{The means of genotypes in the environments, with observed,
#' predicted and residual values.}
#'
#' \item{WAAS}{A data frame with the response variable, the scores of all
#' Principal Components, the estimates of Weighted Average of Absolute Scores,
#' and WAASY (the index that consider the weights for stability and
#' productivity in the genotype ranking.}
#'
#' \item{WAASY}{The values of the WAASY estimated when the wheight for the
#' stability in the loop match with argument "saveWAASY".}
#'
#' \item{WAASY.values}{All the values of WAASY estimated in the different
#' scenarios of WAAS/GY weighting ratio.}
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @templateVar fun WAASBYratio
#' @template template-depr_fun
NULL

#' @templateVar old WAASratio.AMMI
#' @templateVar new wsmp
#' @template template-depr_pkg
#'
#' @export
#'
WAASratio.AMMI <- function(.data, env, gen, rep, resp, mresp = 100, p.valuePC = 0.05, increment = 5,
    saveWAASY = 50, progbar = TRUE) {
    PesoWAAS <- 100
    PesoResp <- 0
    minresp <- 100 - mresp
    test <- PesoWAAS%%increment == 0
    test2 <- saveWAASY%%increment == 0
    if(!mresp %in% c(100, 0)){
        stop("The value 'mresp' must be 0 or 100.")
    }
    if (test == FALSE) {
        stop("The argument 'increment = ", increment, "' is invalid. Please, note that this value must result in an integer in the expression '100 / increment'. Please, consider changing the values.")
    }
    if (test2 == FALSE) {
        stop("The argument 'saveWAASY = ", saveWAASY, "' must be divisible by 'increment' (",
             increment, "). Please, consider changing the values.")
    }
    data = .data
    Y <- eval(substitute(resp), eval(data))
    ENV <- factor(eval(substitute(env), eval(data)))
    GEN <- factor(eval(substitute(gen), eval(data)))
    REP <- factor(eval(substitute(rep), eval(data)))
    Nenv <- length(unique(ENV))
    Ngen <- length(unique(GEN))
    ncomb <- (100/increment) + 1
    CombWAASY <- data.frame(type = matrix(".", (Ngen + Nenv), 1))
    WAASY.Values <- list()
    model <- performs_ammi(ENV, GEN, REP, Y)
    anova <- model$ANOVA
    PC <- model$analysis
    MeansGxE <- model$means
    totalcomb <- ncomb * nrow(PC)
    initial <- 0
    if (progbar == TRUE) {
        pb <- winProgressBar(title = "the model is being built, please, wait.",
                             min = 1, max = totalcomb, width = 570)
    }
    for (k in 1:ncomb) {
        Escores <- model$biplot
        Escores <- Escores %>% mutate(Code = row.names(Escores)) %>%
                   dplyr::select(type, Code, everything())
        SigPC1 <- nrow(PC[which(PC[, 5] < p.valuePC), ])
        Pesos <- as.data.frame( model$analysis[6][c(1:SigPC1), ])
        colnames(Pesos) <- "Percent"
        for (i in 1:ncomb) {
            WAAS <- Escores %>%
                select(contains("PC")) %>%
                abs() %>%
                t() %>%
                as.data.frame() %>%
                slice(1:SigPC1) %>%
                mutate(Percent = Pesos$Percent)

            WAASAbsInicial = Escores %>% mutate(WAAS = sapply(WAAS[, -ncol(WAAS)], weighted.mean, w = WAAS$Percent)) %>%
                group_by(type) %>%
                mutate(PctResp = (mresp - minresp)/(max(Y) - min(Y)) * (Y - max(Y)) + mresp,
                       PctWAAS = (minresp - mresp)/(max(WAAS) - min(WAAS)) * (WAAS - max(WAAS)) + minresp,
                       wRes = PesoResp,
                       wWAAS = PesoWAAS,
                       OrResp = rank(-Y),
                       OrWAAS = rank(WAAS),
                       OrPC1 = rank(abs(PC1)),
                       WAASY = ((PctResp * wRes) + (PctWAAS * wWAAS))/(wRes + wWAAS),
                       OrWAASY = rank(-WAASY)) %>%
                dplyr::ungroup()
            inicial <- as.data.frame(WAASAbsInicial$OrWAAS)
            colnames(inicial) <- paste0(SigPC1, "PCA")
        SigPC2 <- 1
        for (j in 1:nrow(PC)) {
            WAAS <- Escores %>%
                select(contains("PC")) %>%
                abs() %>%
                t() %>%
                as.data.frame() %>%
                slice(1:SigPC2) %>%
                mutate(Percent = Pesos$Percent[1:SigPC2])

            WAASAbs = Escores %>% mutate(WAAS = sapply(WAAS[, -ncol(WAAS)], weighted.mean, w = WAAS$Percent)) %>%
                group_by(type) %>%
                mutate(PctResp = (mresp - minresp)/(max(Y) - min(Y)) * (Y - max(Y)) + mresp,
                       PctWAAS = (minresp - mresp)/(max(WAAS) - min(WAAS)) * (WAAS - max(WAAS)) + minresp,
                       wRes = PesoResp,
                       wWAAS = PesoWAAS,
                       OrResp = rank(-Y),
                       OrWAAS = rank(WAAS),
                       OrPC1 = rank(abs(PC1)),
                       WAASY = ((PctResp * wRes) + (PctWAAS * wWAAS))/(wRes + wWAAS),
                       OrWAASY = rank(-WAASY)) %>%
                dplyr::ungroup()
            results <- as.data.frame(WAASAbs$OrWAAS)
            names(results) <- paste0(SigPC2, "PCA")
            final <- cbind(results, inicial)
            inicial <- final
            SigPC2 %<>% + 1
            ProcdAtua <- j
            initial <- initial + 1
            if (progbar == TRUE) {
                setWinProgressBar(pb, initial, title = paste("Obtaining the Ranks considering",
                                                             ProcdAtua, " of ", nrow(PC), "Principal components:", "|WAAS:",
                                                             PesoWAAS, "% ", "GY:", PesoResp, "%|", round(initial/totalcomb *
                                                                                                                   100, 2), "% Concluded -"))
            }
        }
        initial <- initial
        WAAS <- WAASAbsInicial
        WAAS$type <- ifelse(WAAS$type == "GEN", "Genotype", "Environment")
        CombWAASY[[sprintf("%.0f/%.0f", PesoWAAS, PesoResp)]] <- WAASAbsInicial$OrWAASY
        WAASY.Values[[paste("WAAS/GY", PesoWAAS, "/", PesoResp)]] <- data.frame(WAAS)
        PesoResp <- PesoResp + increment
        PesoWAAS <- PesoWAAS - increment

        if (PesoWAAS + increment == saveWAASY) {
            genotypes = WAAS %>%
                dplyr::filter(type == "Genotype") %>%
                dplyr::select(Code, wRes, wWAAS, WAASY) %>%
                dplyr::arrange(WAASY) %>%
                mutate(Mean = ifelse(WAASY < mean(WAASY), "below", "above"))
        }
    }
        if (progbar == TRUE) {
            close(pb)
        }
        Rank <- final[, -(SigPC2)]
        Names <- WAASAbsInicial %>% select(type, Code, OrResp, OrPC1, OrWAAS)
        Rank <- cbind(Names, Rank)
        hetdata <- Rank %>% dplyr::filter(type == "GEN")
        rownames(hetdata) <- hetdata$Code
        hetdata %<>% select(contains("PCA")) %>% as.matrix()
        CombWAASY %<>% select(-type) %>%
            mutate(type = Names$type, Code = Names$Code) %>%
            select(type, Code, everything()) %>%
            dplyr::filter(type == "GEN")
        hetcomb <- CombWAASY
        rownames(hetcomb) <- CombWAASY$Code
        hetcomb %<>% select(contains("/")) %>% as.matrix()
        CorrRank <- Rank %>% select(type, Code, OrResp, OrPC1, OrWAAS) %>% dplyr::filter(type == "GEN")
        CorcombWAASY <- as.data.frame(cbind(CorrRank, hetcomb))
        rownames(CorcombWAASY) <- CorcombWAASY$Code
        CorcombWAASY %<>% select(-type, -Code) %>%
            rename(Y = OrResp, PCA1 = OrPC1, WAASB = OrWAAS)%>%
            as.matrix()
    PC1 <- Pesos[1, 1]
    PC2 <- Pesos[2, 1]
    mean <- mean(WAAS$Y)
    return(structure(list(anova = anova, PC = PC, MeansGxE = MeansGxE, WAAS = WAAS,
                          WAASY = genotypes, WAASxGY = WAASY.Values,hetcomb = hetcomb, hetdata = hetdata,
                          Ranks = Rank), class = "WAASratio.AMMI"))
}
}
