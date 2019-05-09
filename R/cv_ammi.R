#' Cross-validation for estimation of AMMI model
#'
#' Cross-validation for estimation of AMMI models
#'
#' For each iteration, the original dataset is split into two datasets:
#' modeling and validating data. The dataset "modeling" has all combinations
#' (genotype x environment) with the number of replications informed in
#' \code{nrepval}. The dataset "validating" has one replication. The splitting
#' of the dataset into modeling and validating data depends on the design
#' informed. For Completely Randomized Block Design (default), compltely blocks
#' are selected within environments. The remained block serves validation data.
#' If \code{design = "RCD"} is informed, completely randomly samples are made
#' for each genotype-by-environment combination. The estimated values
#' (depending on NAXIS informed) are compared with the "validating" data. the
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
#' @param design The experimental desig to be considered. Default is
#' \code{RCBD} (Randomized complete Block Design). For Completely Randomized
#' Designs inform \code{design = "CRD"}.
#' @param nboot The number of resamples to be used in the cross-validation
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
#' @export
#' @examples
#'
#' \dontrun{
#' library(METAAB)
#' model = cv_ammi(data_ge,
#'                         env = ENV,
#'                         gen = GEN,
#'                         rep = REP,
#'                         resp = GY,
#'                         nboot = 100,
#'                         nrepval = 2,
#'                         naxis = 2)
#'
#' # Alternatively using the pipe operator %>%
#' library(dplyr)
#' model = data_ge %>%
#'         cv_ammiF(ENV, GEN, REP, GY, 100, 2, 2)
#'
#' }
#'
cv_ammi <- function(.data, env, gen, rep, resp, design = "RCBD", nboot, nrepval,
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
    RMSPDres = RMSPDres %>% mutate(MODEL = paste("AMMI", naxis, sep = "")) %>%
        dplyr::select(MODEL, everything())

    return(structure(list(RMSPD = RMSPDres, RSMEmean = RSMEmean, Estimated = MEDIAS,
                          Modeling = modeling, Testing = testing), class = "cv_ammi"))
}
