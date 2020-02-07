#' Additive Main effects and Multiplicative Interaction
#'
#' Compute the Additive Main effects and Multiplicative interaction. This
#' function also serves as a helper function for other procedures performed in
#' the \pkg{metan} package such as \code{\link{waas}} and \code{\link{wsmp}}
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments
#' @param gen The name of the column that contains the levels of the genotypes
#' @param rep The name of the column that contains the levels of the
#'   replications/blocks
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure a vector of variables may be used. For example \code{resp
#'   = c(var1, var2, var3)}.
#' @param block Defaults to \code{NULL}. In this case, a randomized complete
#'   block design is considered. If block is informed, then a resolvable
#'   alpha-lattice design (Patterson and Williams, 1976) is employed.
#'   \strong{All effects, except the error, are assumed to be fixed.}
#' @param verbose Logical argument. If \code{verbose = FALSE} the code will run
#'   silently.
#' @return
#' * \strong{ANOVA}: The analysis of variance for the AMMI model.
#'
#' * \strong{PCA}: The principal component analysis
#'
#' * \strong{MeansGxE}: The means of genotypes in the environments
#'
#' * \strong{model}: scores for genotypes and environments in all the possible
#' axes.
#'  * \strong{augment:} Information about each observation in the dataset. This
#'  includes predicted values in the \code{fitted} column, residuals in the
#'  \code{resid} column, standardized residuals in the \code{stdres} column,
#'  the diagonal of the 'hat' matrix in the \code{hat}, and standard errors for
#'  the fitted values in the \code{se.fit} column.
#' @md
#' @references Patterson, H.D., and E.R. Williams. 1976. A new class of
#' resolvable incomplete block designs. Biometrika 63:83-92.
#' \href{https://doi.org/10.1093/biomet/63.1.83}{doi:10.1093/biomet/63.1.83}
#' @seealso \code{\link{waas}} \code{\link{waas_means}} \code{\link{waasb}} \code{\link{get_model_data}}
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' model <- performs_ammi(data_ge, ENV, GEN, REP, resp = c(GY, HM))
#'
#' # GY x PC1 (variable GY)
#' plot_scores(model)
#'
#' # PC1 x PC2 (variable HM)
#' plot_scores(model,
#'             var = 2, # or "HM"
#'             type = 2)
#'
#' # Nominal yield plot (variable GY)
#' # Draw a convex hull polygon
#' plot_scores(model, type = 4)
#'
#'}
performs_ammi <- function(.data,
                          env,
                          gen,
                          rep,
                          resp,
                          block = NULL,
                          verbose = TRUE) {
    if(!missing(block)){
        factors  <- .data %>%
            select(ENV = {{env}},
                   GEN = {{gen}},
                   REP = {{rep}},
                   BLOCK = {{block}}) %>%
            mutate_all(as.factor)
    } else{
        factors  <- .data %>%
            select(ENV = {{env}},
                   GEN = {{gen}},
                   REP = {{rep}}) %>%
            mutate_all(as.factor)
    }
    vars <- .data %>% select({{resp}}, -names(factors))
    has_text_in_num(vars)
    vars %<>% select_numeric_cols()
    listres <- list()
    nvar <- ncol(vars)
    for (var in 1:nvar) {
        data <- factors %>%
            mutate(mean = vars[[var]])
        nenv <- nlevels(data$ENV)
        ngen <- nlevels(data$GEN)
        nrep <- nlevels(data$REP)
        minimo <- min(ngen, nenv) - 1
        if(missing(block)){
            model <- aov(mean ~ GEN + ENV + GEN:ENV + ENV/REP, data = data)
            influence <- lm.influence(model)
            augment <- model$model %>%
                add_cols(hat = influence$hat,
                         sigma = influence$sigma,
                         fitted = predict(model),
                         resid = residuals(model),
                         stdres = rstandard(model),
                         se.fit = predict(model, se.fit = TRUE)$se.fit) %>%
                add_cols(factors = concatenate(., GEN, REP, pull = TRUE)) %>%
                column_to_first(ENV, GEN, REP)
        } else{
            model <- aov(mean ~ GEN + ENV + GEN:ENV + ENV/REP/BLOCK, data = data)
            influence <- lm.influence(model)
            augment <- model$model %>%
                add_cols(hat = influence$hat,
                         sigma = influence$sigma,
                         fitted = predict(model),
                         resid = residuals(model),
                         stdres = rstandard(model),
                         se.fit = predict(model, se.fit = TRUE)$se.fit) %>%
                add_cols(factors = concatenate(., GEN, REP, BLOCK, pull = TRUE)) %>%
                column_to_first(ENV, GEN, REP, BLOCK)
        }
        if (minimo < 2) {
            stop("The analysis AMMI is not possible. Both genotypes and environments must have more than two levels.")
        }
        if(missing(block)){
            anova <- anova(model) %>%
                as.data.frame() %>%
                rownames_to_column("Source") %>%
                select_rows(2, 4, 1, 3, 5)
            anova[2, 1] <- "REP(ENV)"
            anova[1, 5] <- anova[1, 4] / anova[2, 4]
            anova[1, 6] <- 1 - pf(anova[1, 5], anova[1, 2], anova[2, 2])
            probint <- anova[4, 6]
        } else{
            anova <- anova(model) %>%
                as.data.frame() %>%
                rownames_to_column("Source") %>%
                select_rows(2, 4, 5, 1, 3, 6)
            anova[2, 1] <- "REP(ENV)"
            anova[3, 1] <- "BLOCK(REP*ENV)"
            probint <- anova[5, 6]
        }
        DFE <- df.residual(model)
        MSE <- deviance(model)/DFE
        MEANS <- data %>%
            group_by(ENV, GEN) %>%
            summarise(Y = mean(mean)) %>%
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
            if (DF <= 0){
                break
            }
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
        SSAMMI <- data.frame(percent, acum,
                             Df = DFAMMI,
                             `Sum Sq` = SS,
                             `Mean Sq` = MSAMMI,
                             `F value` = F.AMMI,
                             Pr.F = PROBF)
        nssammi <- nrow(SSAMMI)
        SSAMMI <- SSAMMI[SSAMMI$Df > 0, ]
        nss <- nrow(SSAMMI)
        row.names(SSAMMI) <- paste("PC", 1:nss, sep = "")
        SCOREG <- U %*% LL^0.5
        SCOREE <- V %*% LL^0.5
        colnames(SCOREG) <- colnames(SCOREE) <- paste("PC", 1:minimo, sep = "")
        bplot <- MEANS %>%
            group_by(GEN) %>%
            summarise(Y = mean(Y)) %>%
            mutate(type = "GEN") %>%
            rename(Code = GEN) %>%
            cbind(., SCOREG) %>%
            rbind(., MEANS %>%
                      group_by(ENV) %>%
                      summarise(Y = mean(Y)) %>%
                      mutate(type = "ENV") %>%
                      rename(Code = ENV) %>%
                      cbind(., SCOREE)) %>%
            column_to_first(type, Code)
        PC <- SSAMMI %>%
            column_to_last(percent, acum) %>%
            rownames_to_column("Source") %>%
            rename(`Sum Sq` = Sum.Sq,
                   `Mean Sq` = Mean.Sq,
                   `F value` = F.value,
                   `Pr(>F)` = Pr.F,
                   Percent = percent,
                   Accumul = acum)
        resid <- anova[nrow(anova), ]
        anova <- rbind_fill(anova[-nrow(anova), ], PC, resid) %>%
            add_rows(Source = "Total",
                     Df = sum(anova$Df),
                     `Sum Sq` = sum(anova$`Sum Sq`),
                     `Mean Sq` = `Sum Sq` / Df,
                     `F value` = NA,
                     `Pr(>F)` = NA,
                     Percent = NA,
                     Accumul = NA) %>%
            as_tibble()
        MeansGxE <- MEANS[, 1:3]
        EscGEN <- subset(bplot, type == "GEN")
        names(EscGEN)[2] <- "GEN"
        names(EscGEN)[3] <- "y"
        EscENV <- subset(bplot, type == "ENV")
        names(EscENV)[2] <- "ENV"
        MeansGxE <- suppressMessages(
            suppressWarnings(
                mutate(MeansGxE,
                       envPC1 = left_join(MeansGxE, EscENV %>% select(ENV, PC1))$PC1,
                       genPC1 = left_join(MeansGxE, EscGEN %>% select(GEN, PC1))$PC1,
                       nominal = left_join(MeansGxE, EscGEN %>% select(GEN, y))$y + genPC1 * envPC1)
            )
        )
        listres[[paste(names(vars[var]))]] <-
            structure(list(ANOVA = anova,
                           PCA = PC %>% rename(PC = Source),
                           MeansGxE = MeansGxE,
                           model = bplot %>% as_tibble(),
                           augment = augment,
                           probint = probint),
                      class = "performs_ammi")
        if (verbose == TRUE) {
            cat("variable", paste(names(vars[var])),"\n")
            cat("---------------------------------------------------------------------------\n")
            cat("AMMI analysis table\n")
            cat("---------------------------------------------------------------------------\n")
            print(as.data.frame(anova), digits = 3, row.names = FALSE)
            cat("---------------------------------------------------------------------------\n\n")
        }
    }
    if (verbose == TRUE) {
        if (length(which(unlist(lapply(listres, function(x) {
            x[["probint"]]
        })) > 0.05)) > 0) {
            cat("------------------------------------------------------------\n")
            cat("Variables with nonsignificant GxE interaction\n")
            cat(names(which(unlist(lapply(listres, function(x) {
                x[["probint"]]
            })) > 0.05)), "\n")
            cat("------------------------------------------------------------\n")
        } else {
            cat("All variables with significant (p < 0.05) genotype-vs-environment interaction\n")
        }
        cat("Done!\n")
    }
    invisible(structure(listres, class = "performs_ammi"))
}


#' Several types of residual plots
#'
#' Residual plots for a output model of class \code{performs_ammi}. Seven types
#' of plots are produced: (1) Residuals vs fitted, (2) normal Q-Q plot for the
#' residuals, (3) scale-location plot (standardized residuals vs Fitted Values),
#' (4) standardized residuals vs Factor-levels, (5) Histogram of raw residuals
#' and (6) standardized residuals vs observation order, and (7) 1:1 line plot.
#'
#'
#' @param x An object of class \code{performs_ammi}.
#' @param ... Additional arguments passed on to the function
#'   \code{\link{residual_plots}}
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method plot performs_ammi
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' model <- performs_ammi(data_ge, ENV, GEN, REP, GY)
#' plot(model)
#' plot(model,
#'      which = c(3, 5),
#'      nrow = 2,
#'      labels = TRUE,
#'      size.lab.out = 4,
#'      align = "v")
#'}
#'
plot.performs_ammi <- function(x, ...) {
    residual_plots(x,  ...)
}








#' Predict the means of a performs_ammi object
#'
#' Predict the means of a performs_ammi object considering a specific number of axis.
#'
#' This function is used to predict the response variable of a two-way table
#' (for examples the yielding of the i-th genotype in the j-th environment)
#' based on AMMI model. This prediction is based on the number of multiplicative
#' terms used. If \code{naxis = 0}, only the main effects (AMMI0) are used. In
#' this case, the predicted mean will be the predicted value from OLS
#' estimation. If \code{naxis = 1} the AMMI1 (with one multiplicative term) is
#' used for predicting the response variable. If \code{naxis =
#' min(gen-1;env-1)}, the AMMIF is fitted and the predicted value will be the
#' cell mean, i.e. the mean of R-replicates of the i-th genotype in the j-th
#' environment. The number of axis to be used must be carefully chosen.
#' Procedures based on Postdictive success (such as Gollobs's d.f.) or
#' Predictive success (such as cross-validation) should be used to do this. This
#' package provide both. \code{\link{performs_ammi}} function compute
#' traditional AMMI analysis showing the number of significant axis. On the
#' other hand, \code{\link{cv_ammif}} function provide a cross-validation,
#' estimating the RMSPD of all AMMI-family models, based on resampling
#' procedures.
#'
#' @param object An object of class performs_ammi
#' @param naxis The the number of axis to be use in the prediction. If
#'   \code{object} has more than one variable, then \code{naxis} must be a
#'   vector.
#' @param ... Additional parameter for the function
#' @return A list where each element is the predicted values by the AMMI model
#'   for each variable.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method predict performs_ammi
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' model <- performs_ammi(data_ge, ENV, GEN, REP,
#'                        resp = c(GY, HM))
#' # Predict GY with 3 IPCA and HM with 1 IPCA
#' predict <- predict(model, naxis = c(3, 1))
#' }
#'
predict.performs_ammi <- function(object, naxis = 2, ...) {
    cal <- match.call()
    if (class(object) != "performs_ammi") {
        stop("The objectin must be an objectin of the class 'performs_ammi'")
    }
    if (length(object) != length(naxis)) {
        stop("The argument 'naxix = ", cal[3], "' must length of ",
             length(object), ", the same number of variables in object '",
             cal[2], "'.")
    }
    listres <- list()
    varin <- 1
    for (var in 1:length(object)) {
        objectin <- object[[var]]
        MEDIAS <- objectin$MeansGxE %>% select(ENV, GEN, Y)
        Nenv <- length(unique(MEDIAS$ENV))
        Ngen <- length(unique(MEDIAS$GEN))
        minimo <- min(Nenv, Ngen) - 1
        if (naxis[var] > minimo) {
            stop("The number of axis to be used must be lesser than or equal to min(GEN-1;ENV-1), in this case, ",
                 minimo, ".")
        } else {
            if (naxis[var] == 0) {
                stop("Invalid argument. The AMMI0 model is calculated automatically. Please, inform naxis > 0")
            } else {
                ovmean <- mean(MEDIAS$Y)
                x1 <- model.matrix(~factor(MEDIAS$ENV) - 1)
                z1 <- model.matrix(~factor(MEDIAS$GEN) - 1)
                modelo1 <- lm(Y ~ ENV + GEN, data = MEDIAS)
                MEDIAS <- mutate(MEDIAS, resOLS = residuals(modelo1))
                intmatrix <- t(matrix(MEDIAS$resOLS, Nenv, byrow = T))
                s <- svd(intmatrix)
                if (length(object) > 1) {
                    U <- s$u[, 1:naxis[var]]
                    LL <- s$d[1:naxis[var]]
                    V <- s$v[, 1:naxis[var]]
                } else {
                    U <- s$u[, 1:naxis]
                    LL <- s$d[1:naxis]
                    V <- s$v[, 1:naxis]
                }
                temp <- mutate(MEDIAS,
                               Ypred = Y - resOLS,
                               ResAMMI = ((z1 %*% U) * (x1 %*% V)) %*% LL,
                               YpredAMMI = Ypred + ResAMMI,
                               AMMI0 = Ypred) %>%
                    as_tibble()
                listres[[paste(names(object[var]))]] <- temp
            }
        }
    }
    invisible(structure(listres, class = "performs_ammi"))
}






#' Print an object of class performs_ammi
#'
#' Print the \code{performs_ammi} object in two ways. By default, the results are shown
#' in the R console. The results can also be exported to the directory.
#'
#'
#' @param x An object of class \code{performs_ammi}.
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported to
#'   the working directory
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print performs_ammi
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' model <- performs_ammi(data_ge, ENV, GEN, REP,
#'                        resp = c(GY, HM))
#' print(model)
#' }
print.performs_ammi <- function(x, export = FALSE, file.name = NULL, digits = 4, ...) {
    if (!class(x) == "performs_ammi") {
        stop("The object must be of class 'performs_ammi'")
    }
    if (export == TRUE) {
        file.name <- ifelse(is.null(file.name) == TRUE, "performs_ammi print", file.name)
        sink(paste0(file.name, ".txt"))
    }
    opar <- options(pillar.sigfig = digits)
    on.exit(options(opar))
    for (i in 1:length(x)) {
        var <- x[[i]]
        cat("Variable", names(x)[i], "\n")
        cat("---------------------------------------------------------------------------\n")
        cat("AMMI analysis table\n")
        cat("---------------------------------------------------------------------------\n")
        print(var$ANOVA)
        cat("---------------------------------------------------------------------------\n")
        cat("Scores for genotypes and environments\n")
        cat("---------------------------------------------------------------------------\n")
        print(var$model)
        cat("\n\n\n")
    }
    if (export == TRUE) {
        sink()
    }
}
