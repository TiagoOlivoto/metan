#' Cross-validation procedure
#'
#' Cross-validation for estimation of AMMI models
#'
#' For each iteration, the original dataset is split into two datasets:
#' modeling and validating data. The dataset 'modeling' has all combinations
#' (genotype x environment) with the number of replications informed in
#' \code{nrepval}. The dataset 'validating' has one replication. The splitting
#' of the dataset into modeling and validating data depends on the design
#' informed. For Completely Randomized Block Design (default), completely blocks
#' are selected within environments. The remained block serves validation data.
#' If \code{design = 'RCD'} is informed, completely randomly samples are made
#' for each genotype-by-environment combination. The estimated values
#' (depending on NAXIS informed) are compared with the 'validating' data. the
#' Root Means Square error is computed. At the end of boots, a list is returned
#' with the following values.
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
#' @param design The experimental desig to be considered. Default is
#' \code{RCBD} (Randomized complete Block Design). For Completely Randomized
#' Designs inform \code{design = 'CRD'}.
#' @param nrepval The number of replicates (r) from total number of replicates
#' (R) to be used in the modeling dataset. Only one replicate is used as
#' validating data each step, so, \code{Nrepval} must be equal \code{R-1}
#' @param naxis The number of axis to be considered for estimation of GE
#' effects.
#' @param verbose A logical argument to define if a progress bar is shown.
#' Default is \code{TRUE}.
#' @return \item{RMSE}{A vector with Nboot-estimates of the root mean squared
#' error estimated with the difference between predicted and validating data.}
#'
#' \item{RSMEmean}{The mean of RMSE estimates.}
#'
#' \item{Estimated}{A data frame that contain the values (predicted, observed,
#' validation) of the last loop.}
#'
#' \item{Modeling}{The dataset used as modeling data in the last loop.}
#'
#' \item{Testing}{The dataset used as testing data in the last loop.}
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso \code{\link{cv_ammif}, \link{cv_blup}}
#' @importFrom progress progress_bar
#' @export
#' @examples
#'
#' \dontrun{
#' library(metan)
#' model = cv_ammi(data_ge,
#'                 env = ENV,
#'                 gen = GEN,
#'                 rep = REP,
#'                 resp = GY,
#'                 nboot = 100,
#'                 nrepval = 2,
#'                 naxis = 2)
#'
#' # Alternatively using the pipe operator %>%
#' model = data_ge %>%
#'         cv_ammiF(ENV, GEN, REP, GY, 100, 2, 2)
#'
#' }
#'
#'
cv_ammi <- function(.data, env, gen, rep, resp, nboot = 100,
                    design = "RCBD", nrepval, naxis, verbose = TRUE) {
  if (!design %in% c("RCBD", "CRD")) {
    stop("Incorrect experimental design informed! Plesease inform RCBD for randomized complete block or CRD for completely randomized design.")
  }
  data <- .data %>% select(ENV = !!enquo(env), GEN = !!enquo(gen),
                           REP = !!enquo(rep), Y = !!enquo(resp))
  RMSPDres <- data.frame(RMSPD = matrix(0, nboot, 1))
  data <- mutate(data, ID = rownames(data))
  Nenv <- length(unique(data$ENV))
  Ngen <- length(unique(data$GEN))
  Nbloc <- length(unique(data$REP))
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
    pb <- progress_bar$new(
      format = "Validating :current of :total sets [:bar] :percent (:elapsedfull -:eta left)",
      clear = FALSE, total = nboot, width = 80)
  }
  condition <- (design == "CRD")
  condition2 <- (design == "RCBD")
  for (b in 1:nboot) {
    if (condition) {
      X <- sample(1:10000, 1)
      set.seed(X)
      modeling <- data %>% dplyr::group_by(ENV, GEN) %>%
        dplyr::sample_n(nrepval, replace = FALSE) %>% arrange(ID) %>%
        as.data.frame()
      rownames(modeling) <- modeling$ID
    }
    if (condition2) {
      tmp <- split_factors(data, ENV, keep_factors = TRUE,
                           verbose = FALSE)
      modeling <- do.call(rbind, lapply(tmp, function(x) {
        X2 <- sample(unique(data$REP), nrepval, replace = FALSE)
        x %>% dplyr::group_by(GEN) %>% dplyr::filter(unique(data$REP) %in%
                                                       c(X2))
      })) %>% as.data.frame()
      rownames(modeling) <- modeling$ID
    }
    testing <- dplyr::anti_join(data, modeling, by = c("ENV",
                                                       "GEN", "REP", "Y", "ID")) %>% arrange(ENV, GEN, REP) %>%
      as.data.frame()
    MEDIAS <- data.frame(modeling %>% dplyr::group_by(ENV,
                                                      GEN) %>% dplyr::summarise(Y = mean(Y)))
    residual <- residuals(lm(Y ~ ENV + GEN, data = MEDIAS))
    s <- svd(t(matrix(residual, Nenv, byrow = T)))
    MEDIAS %<>% mutate(Ypred = Y - residual, ResAMMI = ((model.matrix(~factor(testing$GEN) -
                                                                        1) %*% s$u[, 1:naxis]) * (model.matrix(~factor(testing$ENV) -
                                                                                                                 1) %*% s$v[, 1:naxis])) %*% s$d[1:naxis], YpredAMMI = Ypred +
                         ResAMMI, testing = testing$Y, error = YpredAMMI -
                         testing, errrorAMMI0 = Ypred - testing)
    RMSPDres[, 1][b] <- ifelse(naxis == 0, sqrt(sum(MEDIAS$errrorAMMI0^2)/length(MEDIAS$errrorAMMI0)),
                               sqrt(sum(MEDIAS$error^2)/length(MEDIAS$error)))
    if (verbose == TRUE) {
      pb$tick()
    }
  }
  RMSPDmean <- summarise(MODEL = paste("AMMI", naxis, sep = ""),
                         RMSPDres, mean = mean(RMSPD), sd = sd(RMSPD), se = sd(RMSPD)/sqrt(n()),
                         Q2.5 = quantile(RMSPD, 0.025), Q97.5 = quantile(RMSPD,
                                                                         0.975))
  RMSPDres <- RMSPDres %>% mutate(MODEL = paste("AMMI", naxis, sep = "")) %>%
    dplyr::select(MODEL, everything())
  return(structure(list(RMSPD = RMSPDres, RMSPDmean = RMSPDmean,
                        Estimated = MEDIAS, Modeling = modeling, Testing = testing),
                   class = "cv_ammi"))
}
