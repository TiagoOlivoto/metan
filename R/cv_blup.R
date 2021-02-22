#' Cross-validation procedure
#' @description
#' `r badge('stable')`
#'
#' Cross-validation for blup prediction.
#'
#' This function provides a cross-validation procedure for mixed models using
#' replicate-based data. By default, complete blocks are randomly selected
#' within each environment. In each iteration, the original dataset is split up
#' into two datasets: training and validation data. The 'training' set has all
#' combinations (genotype x environment) with R - 1 replications. The
#' 'validation' set has the remaining replication. The estimated values are
#' compared with the 'validation' data and the Root Means Square Prediction
#' Difference (Olivoto et al. 2019) is computed. At the end of boots, a list is
#' returned.
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#'   replications/blocks. **AT LEAST THREE REPLICATES ARE REQUIRED TO
#'   PERFORM THE CROSS-VALIDATION**.
#' @param resp The response variable.
#' @param block Defaults to `NULL`. In this case, a randomized complete
#'   block design is considered. If block is informed, then a resolvable
#'   alpha-lattice design (Patterson and Williams, 1976) is employed. See how
#'   fixed and random effects are considered, see the section **Details**.
#' @param nboot The number of resamples to be used in the cross-validation.
#'   Defaults to 200
#' @param random The effects of the model assumed to be random. See
#'   **Details** for more information.
#' @param verbose A logical argument to define if a progress bar is shown.
#'   Default is `TRUE`.
#' @details Six models may be fitted depending upon the values in `block`
#'   and `random` arguments.
#'   *  **Model 1:** `block = NULL` and `random = "gen"` (The
#'   default option). This model considers a Randomized Complete Block Design in
#'   each environment assuming genotype and genotype-environment interaction as
#'   random effects. Environments and blocks nested within environments are
#'   assumed to fixed factors.
#'
#'   *  **Model 2:** `block = NULL` and `random = "env"`. This
#'   model considers a Randomized Complete Block Design in each environment
#'   treating environment, genotype-environment interaction, and blocks nested
#'   within environments as random factors. Genotypes are assumed to be fixed
#'   factors.
#'
#'   *  **Model 3:** `block = NULL` and `random = "all"`. This
#'   model considers a Randomized Complete Block Design in each environment
#'   assuming a random-effect model, i.e., all effects (genotypes, environments,
#'   genotype-vs-environment interaction and blocks nested within environments)
#'   are assumed to be random factors.
#'
#'   *  **Model 4:** `block` is not `NULL` and `random =
#'   "gen"`. This model considers an alpha-lattice design in each environment
#'   assuming genotype, genotype-environment interaction, and incomplete blocks
#'   nested within complete replicates as random to make use of inter-block
#'   information (Mohring et al., 2015). Complete replicates nested within
#'   environments and environments are assumed to be fixed factors.
#'
#'   *  **Model 5:** `block` is not `NULL` and `random =
#'   "env"`. This model considers an alpha-lattice design in each environment
#'   assuming genotype as fixed. All other sources of variation (environment,
#'   genotype-environment interaction, complete replicates nested within
#'   environments, and incomplete blocks nested within replicates) are assumed
#'   to be random factors.
#'
#'   *  **Model 6:** `block` is not `NULL` and `random =
#'   "all"`. This model considers an alpha-lattice design in each environment
#'   assuming all effects, except the intercept, as random factors.
#'
#' **IMPORTANT:**  An error is returned if any combination of
#' genotype-environment has a different number of replications than observed in
#' the trial.
#'
#' @return An object of class `cv_blup` with the following items: *
#' **RMSPD**: A vector with nboot-estimates of the root mean squared
#' prediction difference between predicted and validating data. *
#' **RMSPDmean** The mean of RMSPDmean estimates.
#' @references Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, V.S. Marchioro,
#'   V.Q. de Souza, and E. Jost. 2019. Mean performance and stability in
#'   multi-environment trials I: Combining features of AMMI and BLUP techniques.
#'   Agron. J. 111:2949-2960.
#' \doi{10.2134/agronj2019.03.0220}
#'
#' @references Patterson, H.D., and E.R. Williams. 1976. A new class of
#' resolvable incomplete block designs. Biometrika 63:83-92.
#'
#' @references Mohring, J., E. Williams, and H.-P. Piepho. 2015. Inter-block
#'   information: to recover or not to recover it? TAG. Theor. Appl. Genet.
#'   128:1541-54.
#'   \doi{10.1007/s00122-015-2530-0}
#'
#'
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso [cv_ammi()], [cv_ammif()]
#' @export
#' @importFrom utils capture.output
#' @examples
#'
#' \donttest{
#' library(metan)
#' model <- cv_blup(data_ge,
#'                  env = ENV,
#'                  gen = GEN,
#'                  rep = REP,
#'                  resp = GY,
#'                  nboot = 5)
#'
#' }
#'
cv_blup <- function(.data,
                    env,
                    gen,
                    rep,
                    resp,
                    block = NULL,
                    nboot = 200,
                    random = "gen",
                    verbose = TRUE) {
    if(missing(block)){
        data <- .data %>%
            dplyr::select(ENV = {{env}},
                          GEN = {{gen}},
                          REP = {{rep}},
                          Y = {{resp}}) %>%
            mutate(across(1:3, as.factor))
        data <- add_row_id(data)
        Nbloc <- nlevels(data$REP)
        nrepval <- Nbloc - 1
        if (Nbloc <= 2) {
            stop("At least three replicates are required to perform the cross-validation.", call. = FALSE)
        }
        test <-
            data %>%
            n_by(ENV, GEN) %>%
        mutate(test =
            case_when(Y == 0  ~ FALSE,
                      Y > 0 & Y != Nbloc ~ TRUE,
                      TRUE ~ FALSE)
        )
        if(any(test$test == TRUE)){
            df_test <-
                data.frame(
                    test[which(test$test == TRUE),] %>%
                        remove_cols(row_id,REP, test) %>%
                        rename(n = Y)
                )
            message(paste0(capture.output(df_test), collapse = "\n"))
            stop("Combinations of genotype and environment with different number of replication than observed in the trial (", Nbloc, ")", call. = FALSE)
        }
        data <- remove_rows_na(data, verbose = FALSE)
        if (verbose == TRUE) {
            pb <- progress(max = nboot, style = 4)
        }
        RMSPDres <- base::data.frame(RMSPD = matrix(NA, nboot, 1))
        model_formula <- dplyr::case_when(
            random == "gen" ~ paste("Y ~ ENV/REP + (1 | GEN) + (1 | GEN:ENV)"),
            random == "env" ~ paste("Y ~ GEN + (1 | ENV/REP) + (1 | GEN:ENV)"),
            random == "all" ~ paste("Y ~ 1 + (1 | GEN) + (1 | ENV/REP) + (1 | GEN:ENV)")
        )
        model_formula = stats::as.formula(model_formula)
        MOD <- dplyr::case_when(
            random == "gen" ~ "BLUP_g_RCBD",
            random == "env" ~ "BLUP_e_RCBD",
            random == "all" ~ "BLUP_ge_RCBD"
        )
        for (b in 1:nboot) {
            tmp <- split_factors(data, ENV, keep_factors = TRUE, verbose = FALSE)
            modeling <- do.call(rbind, lapply(tmp, function(x) {
                X2 <- sample(unique(data$REP), nrepval, replace = FALSE)
                x %>%
                    group_by(GEN) %>%
                    dplyr::filter(REP %in% X2)
            })) %>%
                base::as.data.frame()
            rownames(modeling) <- modeling$row_id
            testing <- suppressWarnings(anti_join(data, modeling, by = c("ENV", "GEN", "REP", "Y", "row_id"))) %>%
                arrange(ENV, GEN, REP) %>%
                base::as.data.frame()
            MEDIAS <-
                modeling %>%
                means_by(ENV, GEN) %>%
                as.data.frame()

            model <- suppressWarnings(suppressMessages(lme4::lmer(model_formula, data = modeling)))

            validation <-
                modeling %>%
                mutate(pred = predict(model)) %>%
                means_by(ENV, GEN) %>%
                mutate(error = pred - testing$Y)
            RMSPD <- sqrt(sum(validation$error^2)/length(validation$error))
            RMSPDres[, 1][b] <- RMSPD
            if (verbose == TRUE) {
                run_progress(pb,
                             text = paste("Validating", b, "of", nboot, "sets"),
                             actual = b)
            }
        }
    }
    if(!missing(block)){
        data <- .data %>% select(ENV = {{env}},
                                 GEN = {{gen}},
                                 REP = {{rep}},
                                 BLOCK = {{block}},
                                 Y = {{resp}})
        data <- add_row_id(data)
        Nbloc <- nlevels(data$REP)
        nrepval <- Nbloc - 1
        if (Nbloc <= 2) {
            stop("At least three replicates are required to perform the cross-validation.", call. = FALSE)
        }
        test <-
            data %>%
            n_by(ENV, GEN) %>%
            mutate(test =
                       case_when(Y == 0  ~ FALSE,
                                 Y > 0 & Y != Nbloc ~ TRUE,
                                 TRUE ~ FALSE)
            )
        if(any(test$test == TRUE)){
            df_test <-
                data.frame(
                    test[which(test$test == TRUE),] %>%
                        remove_cols(row_id,REP, test) %>%
                        rename(n = Y)
                )
            message(paste0(capture.output(df_test), collapse = "\n"))
            stop("Combinations of genotype and environment with different number of replication than observed in the trial (", Nbloc, ")", call. = FALSE)
        }
        data <- remove_rows_na(data, verbose = FALSE)
        if (verbose == TRUE) {
            pb <- progress(max = nboot, style = 4)
        }
        RMSPDres <- data.frame(RMSPD = matrix(NA, nboot, 1))
        model_formula <- case_when(
            random == "gen" ~ paste("Y ~  (1 | GEN) + ENV / REP + (1|BLOCK:(REP:ENV))  + (1 | GEN:ENV)"),
            random == "env" ~ paste("Y ~ GEN + (1|ENV/REP/BLOCK) + (1 | GEN:ENV)"),
            random == "all" ~ paste("Y ~  (1 | GEN) + (1|ENV/REP/BLOCK) + (1 | GEN:ENV)")
        )
        model_formula = as.formula(model_formula)
        MOD <- case_when(
            random == "gen" ~ "BLUP_g_Alpha",
            random == "env" ~ "BLUP_e_Alpha",
            random == "all" ~ "BLUP_ge_Alpha"
        )
        for (b in 1:nboot) {
            tmp <- split_factors(data, ENV, keep_factors = TRUE, verbose = FALSE)
            modeling <- do.call(rbind, lapply(tmp, function(x) {
                X2 <- sample(unique(data$REP), nrepval, replace = FALSE)
                x %>%
                    group_by(GEN) %>%
                    dplyr::filter(REP %in% X2)
            })) %>%
                as.data.frame()
            rownames(modeling) <- modeling$row_id
            testing <- suppressWarnings(anti_join(data, modeling, by = c("ENV", "GEN", "REP", "BLOCK", "Y", "row_id"))) %>%
                arrange(ENV, GEN, REP, BLOCK) %>%
                as.data.frame()


            model <- suppressWarnings(suppressMessages(lme4::lmer(model_formula, data = modeling)))

            validation <-
                modeling %>%
                mutate(pred = predict(model)) %>%
                means_by(ENV, GEN) %>%
                mutate(error = pred - testing$Y,
                       code = testing$GEN)
            RMSPD <- sqrt(sum(validation$error^2)/length(validation$error))
            RMSPDres[, 1][b] <- RMSPD
            if (verbose == TRUE) {
                run_progress(pb,
                             text = paste("Validating", b, "of", nboot, "sets"),
                             actual = b)
            }
        }
    }

    RMSPDres <- RMSPDres %>%
        mutate(MODEL = MOD) %>%
        select(MODEL, everything())
    RMSPDmean <- RMSPDres %>%
        group_by(MODEL) %>%
        summarise(mean = mean(RMSPD),
                  sd = sd(RMSPD),
                  se = sd(RMSPD)/sqrt(n()),
                  Q2.5 = quantile(RMSPD, 0.025),
                  Q97.5 = quantile(RMSPD, 0.975),
                  .groups = "drop")
    return(structure(list(RMSPD = RMSPDres, RMSPDmean = RMSPDmean),
                     class = "cvalidation"))
}
