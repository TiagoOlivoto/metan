#' Different scenarios of stability and mean performance
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
#' @param mresp A numeric value that will be the new maximum value after reescaling.
#' By default, the variable in \code{resp} is rescaled so that the original maximum
#'  and minimum values are 100 and 0, respectively.
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
WAASBYratio <- function(.data, env, gen, rep, resp,  mresp = 100,
                         increment = 10, saveWAASY = 50, progbar = TRUE) {
  if(!mresp %in% c(100, 0)){
    stop("The value 'mresp' must be 0 or 100.")
  }
  PesoWAAS <- 100
  PesoResp <- 0
  minresp <- 100 - mresp
  Y <- eval(substitute(resp), eval(.data))
  GEN <- factor(eval(substitute(gen), eval(.data)))
  ENV <- factor(eval(substitute(env), eval(.data)))
  REP <- factor(eval(substitute(rep), eval(.data)))
  data <- data.frame(ENV, GEN, REP, Y)
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
  }
  if (test2 == FALSE) {
    stop("The argument 'saveWAASY = ", saveWAASY, "' must be divisible by 'increment' (",
         increment, "). Please, consider changing the values.")
  }
  WAASY.Values <- list()
  initial <- 0
  model <- suppressWarnings(suppressMessages(lme4::lmer(Y ~ REP %in% ENV +
                                                          ENV + (1 | GEN) + (1 | GEN:ENV))))
  summ <- summary(model)
  bups <- lme4::ranef(model)
  blups <- data.frame(Names = rownames(bups$`GEN:ENV`))
  blups = blups %>%
    data.frame(do.call("rbind",
                       strsplit(as.character(blups$Names),
                                ":", fixed = TRUE))) %>%
    dplyr::select(-Names) %>%
    dplyr::select(-X1, everything()) %>%
    dplyr::mutate(BLUPge = bups[[1]]$`(Intercept)`) %>%
    dplyr::rename(Code = X2, GEN = X1) %>%
    dplyr::arrange(Code)
  intmatrix <- by(blups[, 3], blups[, c(2, 1)], function(x) sum(x, na.rm = TRUE))
  s <- svd(intmatrix)
  U <- s$u[, 1:minimo]
  LL <- diag(s$d[1:minimo])
  V <- s$v[, 1:minimo]
  Eigenvalue <- data.frame(Eigenvalue = s$d[1:minimo]^2) %>%
    dplyr::mutate(Proportion = s$d[1:minimo]^2/sum(s$d[1:minimo]^2) * 100,
                  Accumulated = cumsum(Proportion),
                  PC = paste("PC", 1:minimo, sep = "")) %>%
    dplyr::select(PC, everything())
  SCOREG <- U %*% LL^0.5
  SCOREE <- V %*% LL^0.5
  colnames(SCOREG) <- colnames(SCOREE) <- paste("PC", 1:minimo, sep = "")
  MEDIAS <- data %>% select(ENV, GEN, Y) %>%
    dplyr::group_by(ENV, GEN) %>%
    dplyr::summarise(Y = mean(Y)) %>%
    ungroup()

  MGEN = MEDIAS %>% group_by(GEN) %>% summarise(Y = mean(Y)) %>% mutate(type = "GEN")
  MGEN = cbind(MGEN, SCOREG) %>% rename(Code = GEN)
  MENV = MEDIAS %>% group_by(ENV) %>% summarise(Y = mean(Y)) %>% mutate(type = "ENV")
  MENV = cbind(MENV, SCOREE) %>% rename(Code = ENV)
  Escores <- rbind(MGEN, MENV) %>% select(type, everything())
  Pesos <- data.frame(Percent = Eigenvalue$Proportion)
  if (progbar == TRUE) {
    pb <- winProgressBar(title = "the model is being built, please, wait.",
                         min = 1, max = totalcomb, width = 570)
  }
  for (k in 1:ncomb) {
    WAASB <- Escores %>%
      select(contains("PC")) %>%
      abs() %>%
      t() %>%
      as.data.frame() %>%
      mutate(Percent = Pesos$Percent)
    WAASAbsInicial = Escores %>% mutate(WAASB = sapply(WAASB[, -ncol(WAASB)], weighted.mean, w = WAASB$Percent)) %>%
      group_by(type) %>%
      mutate(PctResp = (mresp - minresp)/(max(Y) - min(Y)) * (Y - max(Y)) + mresp,
             PctWAASB = (minresp - mresp)/(max(WAASB) - min(WAASB)) * (WAASB - max(WAASB)) + minresp,
             wRes = PesoResp,
             wWAASB = PesoWAAS,
             OrResp = rank(-Y),
             OrWAASB = rank(WAASB),
             OrPC1 = rank(abs(PC1)),
             WAASBY = ((PctResp * wRes) + (PctWAASB * wWAASB))/(wRes + wWAASB),
             OrWAASBY = rank(-WAASBY)) %>%
      dplyr::ungroup()
    inicial <- as.data.frame(WAASAbsInicial$OrWAASB)
    colnames(inicial) <- paste0(minimo, "PC")

    SigPC2 <- 1
    for (j in 1:minimo) {
      WAASB <- Escores %>%
        select(contains("PC")) %>%
        abs() %>%
        t() %>%
        as.data.frame() %>%
        slice(1:SigPC2) %>%
        mutate(Percent = Pesos$Percent[1:SigPC2])


      WAASAbs = Escores %>% mutate(WAASB = sapply(WAASB[, -ncol(WAASB)], weighted.mean, w = WAASB$Percent)) %>%
        group_by(type) %>%
        mutate(PctResp = (mresp - minresp)/(max(Y) - min(Y)) * (Y - max(Y)) + mresp,
               PctWAASB = (minresp - mresp)/(max(WAASB) - min(WAASB)) * (WAASB - max(WAASB)) + minresp,
               wRes = PesoResp,
               wWAASB = PesoWAAS,
               OrResp = rank(-Y),
               OrWAASB = rank(WAASB),
               OrPC1 = rank(abs(PC1)),
               WAASBY = ((PctResp * wRes) + (PctWAASB * wWAASB))/(wRes + wWAASB),
               OrWAASBY = rank(-WAASBY)) %>%
        dplyr::ungroup()
      results <- as.data.frame(WAASAbs$OrWAASB)
      names(results) <- paste0(SigPC2, "PCA")
      final <- cbind(results, inicial)
      inicial <- final
      SigPC2 <- SigPC2 + 1
      ProcdAtua <- j
      initial <- initial + 1
      if (progbar == TRUE) {
        setWinProgressBar(pb, initial, title = paste("Obtaining the Ranks considering",
                                                     ProcdAtua, " of ", minimo, "Principal components:", "|WAAS:",
                                                     PesoWAAS, "% ", "GY:", PesoResp, "%|", round(initial/totalcomb *
                                                                                                    100, 2), "% Concluded -"))
      }
    }
    initial <- initial
    WAAS <- WAASAbsInicial
    WAAS$type <- ifelse(WAAS$type == "GEN", "Genotype", "Environment")
    CombWAASY[[sprintf("%.0f/%.0f", PesoWAAS, PesoResp)]] <- WAASAbsInicial$OrWAASBY
    WAASY.Values[[paste("WAAS/GY", PesoWAAS, "/", PesoResp)]] <- data.frame(WAAS)
    PesoResp <- PesoResp + increment
    PesoWAAS <- PesoWAAS - increment
    if (PesoWAAS + increment == saveWAASY) {
      genotypes = WAAS %>%
        dplyr::filter(type == "Genotype") %>%
        dplyr::select(Code, wRes, wWAASB, WAASBY) %>%
        dplyr::arrange(WAASBY) %>%
        mutate(Mean = ifelse(WAASBY < mean(WAASBY), "below", "above"))
    }
  }
  if (progbar == TRUE) {
    close(pb)
  }
  Rank <- final[, -(SigPC2)]
  Names <- WAASAbsInicial %>% select(type, Code, OrResp, OrPC1, OrWAASB)
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
  CorrRank <- Rank %>% select(type, Code, OrResp, OrPC1, OrWAASB) %>% dplyr::filter(type == "GEN")
  CorcombWAASY <- as.data.frame(cbind(CorrRank, hetcomb))
  rownames(CorcombWAASY) <- CorcombWAASY$Code
  CorcombWAASY %<>% select(-type, -Code) %>%
    rename(Y = OrResp, PCA1 = OrPC1, WAASB = OrWAASB)%>%
    as.matrix()
  PC1 <- Pesos[1, 1]
  PC2 <- Pesos[2, 1]
  mean <- mean(WAAS$Y)
  return(structure(list(MeansGxE = MEDIAS, WAASxGY = WAASY.Values, WAAS = WAAS,
                        WAASY = genotypes, hetcomb = hetcomb, hetdata = hetdata, CorcombWAASY = CorcombWAASY,
                        Ranks = Rank), class = "WAASBYratio"))

}

