#' Weighting stability and mean performance in different scenarios
#'
#' This function computes the WAASBY index in mixed-effec model analysis in
#' different combinations of weights for stability and productivity.
#'
#' This function considers both stability (weighted average of absolute scores
#' based on SVD of BLUP-interaction effects matrix) and productivitye for
#' genotype ranking. This function provide the option of attributing weights
#' for stability and productive in genotype ranking. This is important
#' depending on the goal of a selection strategy. For example, if a a goal of a
#' breeding program is to select a genotype whith high yielding (independeltly
#' on the stability performance), that genotype with the first rank in an
#' WAASB/GY = 0/100 ratio should be selected. The reciprocal is true. Aiming at
#' selecting a high-stable genotype (independentely on the productivity), that
#' genotype with the first rank in an WAASB/GY = 100/0 ratio should be
#' selected. By defalut, the increment on the WAASB/GY ratio is equal to 5. In
#' other words, twenty one different combinations are computed. Each
#' combination, the genotypes are ranked regarding the WAASY value.
#'
#' @param .data The dataset containing the columns related to Environments,
#' Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#' environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#' replications/blocks.
#' @param resp The response variable.
#' @param increment The range of the increment for WAASB/GY ratio. Default is
#' 5. The function compute the WAASBY values starting with a weight o 100 for
#' stability and 0 for response variable. With the default, the first scenario
#' will be a WAASB/GY ratio = 100/0. In the next scenario, the WAASBY values
#' are computed based on a WAASB/GY ratio = 95/5.
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
#' library(metan)
#'
#' # Default, with increment of 5 and saving the WAASY values when weight is 50
#' wratio = WAASBYratio(data_ge,
#'                      resp = GY,
#'                      gen = GEN,
#'                      env = ENV,
#'                      rep = REP)
#'
#' # Incrementing 2-by-2
#' wratio2 = WAASBYratio(data_ge,
#'                       resp = GY,
#'                       gen = GEN,
#'                       env = ENV,
#'                       rep = REP,
#'                       increment = 50)
#' }
#'
WAASBYratio <- function(.data, env, gen, rep, resp, increment = 5, saveWAASY = 50,
    progbar = TRUE) {
  data = .data
    PesoWAAS <- 100
    PesoResp <- 0
    Y <- eval(substitute(resp), eval(data))
    GEN <- factor(eval(substitute(gen), eval(data)))
    ENV <- factor(eval(substitute(env), eval(data)))
    REP <- factor(eval(substitute(rep), eval(data)))
    Nenv <- length(unique(ENV))
    Ngen <- length(unique(GEN))
    minimo <- min(Nenv, Ngen) - 1
    ncomb <- (100/increment) + 1
    totalcomb <- ((100/increment) + 1) * minimo
    CombWAASY <- data.frame(type = matrix(".", (Ngen + Nenv), 1))
    ovmean <- mean(Y)

    test <- PesoWAAS%%increment == 0
    test2 <- saveWAASY%%increment == 0

    if (test == FALSE) {
        stop("The argument 'increment = ", increment, "' is invalid. Please, note that this value must result in an integer in the expression '100 / increment'. Please, consider changing the values.")
    } else {

        if (test2 == FALSE) {
            stop("The argument 'saveWAASY = ", saveWAASY, "' must be divisible by 'increment' (",
                increment, "). Please, consider changing the values.")
        } else {

            WAASY.Values <- list()
            initial <- 0
            model <- suppressWarnings(suppressMessages(lme4::lmer(Y ~ REP %in% ENV +
                ENV + (1 | GEN) + (1 | GEN:ENV))))
            summ <- summary(model)
            bups <- lme4::ranef(model)
            blupGEN <- bups$GEN
            blupINT <- bups$`GEN:ENV`
            blups <- data.frame(Names = rownames(blupINT))
            blups <- data.frame(do.call("rbind", strsplit(as.character(blups$Names),
                ":", fixed = TRUE)))
            blups <- blups %>% select(-X1, everything())
            blups <- cbind(blups, blupINT)
            names(blups) <- c("Code", "GEN", "BLUPge")
            blups <- dplyr::arrange(blups, Code)
            intmatrix <- by(blups[, 3], blups[, c(2, 1)], function(x) sum(x, na.rm = TRUE))
            s <- svd(intmatrix)
            U <- s$u[, 1:minimo]
            LL <- diag(s$d[1:minimo])
            V <- s$v[, 1:minimo]
            Eigenvalue <- data.frame(Eigenvalue = s$d[1:minimo]^2)
            Eigenvalue <- mutate(Eigenvalue, Proportion = s$d[1:minimo]^2/sum(s$d[1:minimo]^2) *
                100)
            Eigenvalue <- mutate(group_by(Eigenvalue, Proportion), Accumulated = cumsum(Proportion))
            Eigenvalue$PC <- rownames(Eigenvalue)
            Eigenvalue <- Eigenvalue %>% select(PC, everything())
            SCOREG <- U %*% LL^0.5
            SCOREE <- V %*% LL^0.5
            Escores <- rbind(SCOREG, SCOREE)
            colnames(Escores) <- paste("PC", 1:minimo, sep = "")
            raw <- data.frame(ENV, GEN, Y)
            MEDIAS <- tapply(raw[, 3], raw[, c(1, 2)], mean)
            xx <- rownames(MEDIAS)
            yy <- colnames(MEDIAS)
            fila <- length(xx)
            col <- length(yy)
            total <- fila * col
            x <- character(length = total)
            y <- character(length = total)
            z <- numeric(length = total)
            k <- 0
            for (i in 1:fila) {
                for (j in 1:col) {
                  k <- k + 1
                  x[k] <- xx[i]
                  y[k] <- yy[j]
                  z[k] <- MEDIAS[i, j]
                }
            }
            MEDIAS <- data.frame(ENV = x, GEN = y, Y = z)
            OUTMED <- by(MEDIAS[, 3], MEDIAS[, c(2, 1)], function(x) sum(x, na.rm = TRUE))
            MEscores <- Escores[1:Ngen, ]
            NEscores <- Escores[(Ngen + 1):(Ngen + Nenv), ]
            MGEN <- data.frame(type = "GEN", Code = rownames(OUTMED), Y = apply(OUTMED,
                1, mean), MEscores)
            MENV <- data.frame(type = "ENV", Code = colnames(OUTMED), Y = apply(OUTMED,
                2, mean), NEscores)
            Escores <- rbind(MGEN, MENV)
            Pesos <- data.frame(Percent = Eigenvalue$Proportion)
            if (progbar == TRUE) {
                pb <- winProgressBar(title = "the model is being built, please, wait.",
                  min = 1, max = totalcomb, width = 570)
            }
            for (k in 1:ncomb) {
                WAAS <- Escores
                WAASAbs <- Escores
                for (i in 4:ncol(WAAS)) {
                  WAAS[, i] <- abs(WAAS[i])
                }
                t_WAAS <- data.frame(t(WAAS))
                colnames(t_WAAS) <- rownames(WAAS)
                rownames(t_WAAS) <- colnames(WAAS)
                t_WAAS <- t_WAAS[-c(1, 2, 3), ]
                t_WAAS <- t_WAAS[c(1:minimo), ]
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
                colnames(inicial) <- paste0(minimo, "PC")
                SigPC2 <- 1
                for (j in 1:minimo) {
                  Escores <- rbind(MGEN, MENV)
                  Pesos <- data.frame(Percent = Eigenvalue$Proportion)
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
                  names(results) <- paste0(SigPC2, "PCA")
                  final <- cbind(results, inicial)
                  inicial <- final
                  SigPC2 <- SigPC2 + 1
                  ProcdAtua <- j
                  Sys.sleep(0.1)
                  initial <- initial + 1
                  if (progbar == TRUE) {
                    setWinProgressBar(pb, initial, title = paste("Obtaining the Ranks considering",
                      ProcdAtua, " of ", minimo, "Principal components:", "|WAAS:",
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
            if (progbar == TRUE) {
                close(pb)
            }
            Rank <- final[, -(SigPC2)]
            Names <- WAASAbsInicial[, c("type", "Code", "OrResp", "OrPC1", "OrWAAS")]
            Rank <- cbind(Names, Rank)
            hetdata <- as.data.frame(subset(Rank, type == "GEN"))
            rownames(hetdata) <- hetdata$Code
            hetdata <- hetdata[, -c(1, 2, 3, 4, 5)]
            hetdata <- as.matrix(hetdata)
            CombWAASY <- CombWAASY[, -1]
            CombWAASY <- cbind(Names[, c(1, 2)], CombWAASY)
            hetcomb <- as.data.frame(subset(CombWAASY, type == "GEN"))
            rownames(hetcomb) <- hetcomb$Code
            hetcomb <- hetcomb[, -c(1, 2)]
            hetcomb <- as.matrix(hetcomb)
            CorrRank <- Rank[, c("type", "Code", "OrResp", "OrPC1", "OrWAAS")]
            CorrRank <- subset(CorrRank, type == "GEN")
            CorcombWAASY <- as.data.frame(cbind(CorrRank, hetcomb))
            rownames(CorcombWAASY) <- CorcombWAASY$Code
            CorcombWAASY <- CorcombWAASY[, -c(1, 2)]
            names(CorcombWAASY)[1] <- "Y"
            names(CorcombWAASY)[2] <- "PCA1"
            names(CorcombWAASY)[3] <- "WAASB"
            CorcombWAASY <- as.matrix(CorcombWAASY)
            PC1 <- Pesos[1, 1]
            PC2 <- Pesos[2, 1]
            mean <- mean(WAAS$Y)
            return(structure(list(MeansGxE = MEDIAS, WAASxGY = WAASY.Values, WAAS = WAAS,
                WAASY = genotypes, hetcomb = hetcomb, hetdata = hetdata, CorcombWAASY = CorcombWAASY,
                Ranks = Rank), class = "WAASBYratio"))
            if (progbar == TRUE) {
            }
        }
    }
}
