#' Weighting stability and mean performance in different scenarios
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
#' @export
#' @examples
#'
#' \dontrun{
#' # Default, with increment of 5 and saving the WAASY values when weight is 50
#' wratio = WAASratio.AMMI(data_ge,
#'                         resp = GY,
#'                         gen = GEN,
#'                         env = ENV,
#'                         rep = REP)
#'
#' # Incrementing 2-by-2
#' wratio2 = WAASratio.AMMI(data_ge,
#'                          resp = GY,
#'                          gen = GEN,
#'                          env = ENV,
#'                          rep = REP,
#'                          increment = 50)
#' }
#'
WAASratio.AMMI <- function(.data, env, gen, rep, resp, p.valuePC = 0.05, increment = 5,
    saveWAASY = 50, progbar = TRUE) {
    data = .data
    PesoWAAS <- 100
    PesoResp <- 0
    Y <- eval(substitute(resp), eval(data))
    ENV <- factor(eval(substitute(env), eval(data)))
    GEN <- factor(eval(substitute(gen), eval(data)))
    REP <- factor(eval(substitute(rep), eval(data)))
    Nenv <- length(unique(ENV))
    Ngen <- length(unique(GEN))
    ncomb <- (100/increment) + 1

    test <- PesoWAAS%%increment == 0
    test2 <- saveWAASY%%increment == 0

    if (test == FALSE) {
        stop("The argument 'increment = ", increment, "' is invalid. Please, note that this value must result in an integer in the expression '100 / increment'. Please, consider changing the values.")
    } else {

        if (test2 == FALSE) {
            stop("The argument 'saveWAASY = ", saveWAASY, "' must be divisible by 'increment' (",
                increment, "). Please, consider changing the values.")
        } else {

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
                Escores <- cbind(Code = row.names(Escores), Escores)
                Escores <- Escores %>% dplyr::select(type, everything())
                SigPC1 <- nrow(PC[which(PC[, 5] < p.valuePC), ])
                Pesos <- model$analysis[1]
                Pesos <- as.data.frame(Pesos[c(1:SigPC1), ])
                colnames(Pesos) <- "Percent"
                WAAS <- Escores
                WAASAbs <- Escores
                for (i in 4:ncol(WAAS)) {
                  WAAS[, i] <- abs(WAAS[i])
                }
                t_WAAS <- data.frame(t(WAAS))
                colnames(t_WAAS) <- rownames(WAAS)
                rownames(t_WAAS) <- colnames(WAAS)
                t_WAAS <- t_WAAS[-c(1, 2, 3), ]
                t_WAAS <- t_WAAS[c(1:SigPC1), ]
                t_WAAS <- cbind(t_WAAS, Pesos)
                for (i in 1:ncol(t_WAAS)) {
                  t_WAAS[, i] <- as.numeric(as.character(t_WAAS[, i]))
                }
                Ponderado <- t(as.data.frame(sapply(t_WAAS[, -ncol(t_WAAS)], weighted.mean,
                  w = t_WAAS$Percent)))
                rownames(Ponderado) <- c("WAAS")
                t_WAAS <- subset(t_WAAS, select = -Percent)
                colnames(Ponderado) <- colnames(t_WAAS)
                t_WAAS <- rbind(t_WAAS, Ponderado)
                t_WAAS2 <- data.frame(t(t_WAAS))
                colnames(t_WAAS2) <- rownames(t_WAAS)
                rownames(t_WAAS2) <- colnames(t_WAAS)
                WAASAbs <- cbind(WAASAbs, subset(t_WAAS2, select = WAAS))
                WAASAbs2 <- subset(WAASAbs, type == "ENV")
                WAASAbs2$PctResp <- resca(WAASAbs2$Y, 0, 100)
                WAASAbs2$PctWAAS <- resca(WAASAbs2$WAAS, 100, 0)
                WAASAbs3 <- subset(WAASAbs, type == "GEN")
                WAASAbs3$PctResp <- resca(WAASAbs3$Y, 0, 100)
                WAASAbs3$PctWAAS <- resca(WAASAbs3$WAAS, 100, 0)
                WAASAbs <- rbind(WAASAbs3, WAASAbs2) %>% dplyr::group_by(type) %>%
                  dplyr::mutate(OrResp = rank(-Y), OrWAAS = rank(WAAS), OrPC1 = rank(abs(PC1)),
                    PesRes = as.vector(PesoResp), PesWAAS = as.vector(PesoWAAS), WAASY = (PctResp *
                      PesRes + PctWAAS * PesWAAS)/(PesRes + PesWAAS))
                WAASAbsInicial <- WAASAbs %>% dplyr::group_by(type) %>% dplyr::mutate(OrWAASY = rank(-WAASY)) %>%
                  dplyr::ungroup()
                inicial <- as.data.frame(WAASAbsInicial$OrWAAS)
                colnames(inicial) <- paste0(SigPC1, "PC")
                SigPC2 <- 1
                for (j in 1:nrow(PC)) {
                  Escores <- model$biplot
                  Escores <- cbind(Code = row.names(Escores), Escores)
                  Escores <- Escores %>% dplyr::select(type, everything())
                  Pesos <- model$analysis[1]
                  Pesos <- as.data.frame(Pesos[c(1:SigPC2), ])
                  colnames(Pesos) <- "Percent"
                  WAAS <- Escores
                  WAASAbs <- Escores
                  for (i in 4:ncol(WAAS)) {
                    WAAS[, i] <- abs(WAAS[i])
                  }
                  t_WAAS <- data.frame(t(WAAS))
                  colnames(t_WAAS) <- rownames(WAAS)
                  rownames(t_WAAS) <- colnames(WAAS)
                  t_WAAS <- t_WAAS[-c(1, 2, 3), ]
                  t_WAAS <- t_WAAS[c(1:SigPC2), ]
                  t_WAAS <- cbind(t_WAAS, Pesos)
                  for (i in 1:ncol(t_WAAS)) {
                    t_WAAS[, i] <- as.numeric(as.character(t_WAAS[, i]))
                  }
                  Ponderado <- t(as.data.frame(sapply(t_WAAS[, -ncol(t_WAAS)], weighted.mean,
                    w = t_WAAS$Percent)))
                  rownames(Ponderado) <- c("WAAS")
                  t_WAAS <- subset(t_WAAS, select = -Percent)
                  colnames(Ponderado) <- colnames(t_WAAS)
                  t_WAAS <- rbind(t_WAAS, Ponderado)
                  t_WAAS2 <- data.frame(t(t_WAAS))
                  colnames(t_WAAS2) <- rownames(t_WAAS)
                  rownames(t_WAAS2) <- colnames(t_WAAS)
                  WAASAbs <- cbind(WAASAbs, subset(t_WAAS2, select = WAAS))
                  WAASAbs2 <- subset(WAASAbs, type == "ENV")
                  WAASAbs2$PctResp <- resca(WAASAbs2$Y, 0, 100)
                  WAASAbs2$PctWAAS <- resca(WAASAbs2$WAAS, 100, 0)
                  WAASAbs3 <- subset(WAASAbs, type == "GEN")
                  WAASAbs3$PctResp <- resca(WAASAbs3$Y, 0, 100)
                  WAASAbs3$PctWAAS <- resca(WAASAbs3$WAAS, 100, 0)
                  WAASAbs <- rbind(WAASAbs3, WAASAbs2) %>% dplyr::mutate(PesRes = as.vector(PesoResp),
                    PesWAAS = as.vector(PesoWAAS), WAASY = (PctResp * PesRes + PctWAAS *
                      PesWAAS)/(PesRes + PesWAAS)) %>% dplyr::group_by(type) %>% dplyr::mutate(OrResp = rank(-Y),
                    OrWAAS = rank(WAAS), OrPC1 = rank(abs(PC1)), OrWAASY = rank(-WAASY)) %>%
                    dplyr::ungroup()
                  results <- as.data.frame(WAASAbs$OrWAAS)
                  names(results) <- paste0(SigPC2, "PC")
                  final <- cbind(results, inicial)
                  inicial <- final
                  SigPC2 <- SigPC2 + 1
                  ProcdAtua <- j
                  initial <- initial + 1
                  Sys.sleep(0.1)
                  if (progbar == TRUE) {
                    setWinProgressBar(pb, initial, title = paste("Obtaining the Ranks considering",
                      ProcdAtua, " of ", nrow(PC), "Principal components:", "|WAAS:",
                      PesoWAAS, "% ", "GY:", PesoResp, "%|", "-", round(initial/totalcomb *
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
                  genotypes <- subset(WAAS, type == "Genotype")
                  genotypes <- genotypes[, c("Code", "PesRes", "PesWAAS", "WAASY")]
                  genotypes <- genotypes[order(genotypes$WAASY), ]
                  genotypes$Code <- factor(genotypes$Code, levels = genotypes$Code)
                  genotypes$Mean <- ifelse(genotypes$WAASY < mean(genotypes$WAASY),
                    "below", "above")
                }
            }
            Rank <- final[, -(SigPC2)]
            Names <- WAASAbsInicial[, c("type", "Code", "OrResp", "OrPC1", "OrWAAS")]
            Rank <- cbind(Names, Rank)
            hetdata <- as.data.frame(subset(Rank, type == "GEN"))
            rownames(hetdata) <- hetdata$Code
            hetdata <- hetdata[, -c(1, 2, 3, 4)]
            colnames(hetdata)[1:1] <- c("WAAS")
            hetdata <- as.matrix(hetdata)
            CombWAASY <- CombWAASY[, -1]
            CombWAASY <- cbind(Names[, c(1, 2)], CombWAASY)
            hetcomb <- as.data.frame(subset(CombWAASY, type == "GEN"))
            rownames(hetcomb) <- hetcomb$Code
            hetcomb <- hetcomb[, -c(1, 2)]
            hetcomb <- as.matrix(hetcomb)
            CorrRank <- Rank[, c("type", "Code", "OrResp", "OrWAAS")]
            CorrCWA <- CombWAASY[, -c(1, 2)]
            CorcombWAASY <- as.data.frame(cbind(CorrRank, CorrCWA))
            CorcombWAASY <- subset(CorcombWAASY, type == "GEN")
            rownames(CorcombWAASY) <- CorcombWAASY$Code
            CorcombWAASY <- CorcombWAASY[, -c(1)]
            PC1 <- Pesos[1, 1]
            PC2 <- Pesos[2, 1]
            mean <- mean(WAAS$Y)
            if (progbar == TRUE) {
                close(pb)
            }
            return(structure(list(anova = anova, PC = PC, MeansGxE = MeansGxE, WAAS = WAAS,
                WAASxGY = WAASY.Values, WAASY = genotypes, hetcomb = hetcomb, hetdata = hetdata,
                Ranks = Rank), class = "WAASratio.AMMI"))
        }
    }
}
