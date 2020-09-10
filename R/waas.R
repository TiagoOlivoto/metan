#' Weighted Average of Absolute Scores
#'
#' Compute the Weighted Average of Absolute Scores for AMMI analysis (Olivoto et
#' al., 2019).
#'
#' This function compute the weighted average of absolute scores, estimated as
#' follows:
#'\loadmathjax
#'\mjsdeqn{WAAS_i = \sum_{k = 1}^{p} |IPCA_{ik} \times EP_k|/ \sum_{k =
#'1}^{p}EP_k}
#'
#' where \mjseqn{WAAS_i} is the weighted average of absolute scores of the
#' \emph{i}th genotype; \mjseqn{IPCA_{ik}} is the score of the \emph{i}th genotype
#' in the \emph{k}th IPCA; and \mjseqn{EP_k} is the explained variance of the *k*th
#' IPCA for \emph{k = 1,2,..,p}, considering \emph{p} the number of significant
#' PCAs, or a declared number of PCAs. For example if \code{prob = 0.05}, all
#' axis that are significant considering this probability level are used. The
#' number of axis can be also informed by declaring \code{naxis = x}. This will
#' override the number of significant axes according to the argument code{prob}.
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#'   replications/blocks.
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure a vector of variables may be used. For example \code{resp
#'   = c(var1, var2, var3)}.
#' @param block Defaults to \code{NULL}. In this case, a randomized complete
#'   block design is considered. If block is informed, then a resolvable
#'   alpha-lattice design (Patterson and Williams, 1976) is employed.
#'   \strong{All effects, except the error, are assumed to be fixed.}
#' @param mresp  The new maximum value after rescaling the response variable. By
#'   default, all variables in \code{resp} are rescaled so that de maximum value
#'   is 100 and the minimum value is 0 (i.e., \code{mresp = 100}). It must be a
#'   numeric vector of the same length of \code{resp} if rescaling is assumed to
#'   be different across variables, e.g., if for the first variable smaller
#'   values are better and for the second one, higher values are better, then
#'   \code{mresp = c(0, 100)} must be used. Numeric value of length 1 will be
#'   recycled with a warning message.
#' @param wresp The weight for the response variable(s) for computing the WAASBY
#'   index. By default, all variables in \code{resp} have equal weights for mean
#'   performance and stability (i.e., \code{wresp = 50}). It must be a numeric
#'   vector of the same length of \code{resp} to assign different weights across
#'   variables, e.g., if for the first variable equal weights for mean
#'   performance and stability are assumed and for the second one, a higher
#'   weight for mean performance (e.g. 65) is assumed, then \code{wresp = c(50,
#'   65)} must be used. Numeric value of length 1 will be recycled with a
#'   warning message.
#' @param prob The p-value for considering an interaction principal component
#'   axis significant.
#' @param naxis The number of IPCAs to be used for computing the WAAS index.
#'   Default is \code{NULL} (Significant IPCAs are used). If values are
#'   informed, the number of IPCAS will be used independently on its
#'   significance. Note that if two or more variables are included in
#'   \code{resp}, then \code{naxis} must be a vector.
#' @param ind_anova Logical argument set to \code{FALSE}. If \code{TRUE} an
#'   within-environment ANOVA is performed.
#' @param verbose Logical argument. If \code{verbose = FALSE} the code is run
#'   silently.
#' @references Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, V.S. Marchioro,
#'   V.Q. de Souza, and E. Jost. 2019a. Mean performance and stability in
#'   multi-environment trials I: Combining features of AMMI and BLUP techniques.
#'   Agron. J. 111:2949-2960. \doi{10.2134/agronj2019.03.0220}
#'
#' @return An object of class \code{waas} with the following items for each
#'   variable:
#'
#' * \strong{individual} A within-environments ANOVA considering a fixed-effect
#' model.
#' * \strong{model} A data frame with the response variable, the scores of all
#' Principal Components, the estimates of Weighted Average of Absolute Scores,
#' and WAASY (the index that consider the weights for stability and productivity
#' in the genotype ranking.
#'
#' * \strong{MeansGxE} The means of genotypes in the environments
#'
#' * \strong{PCA} Principal Component Analysis.
#'
#' * \strong{anova} Joint analysis of variance for the main effects and
#' Principal Component analysis of the interaction effect.
#'
#' * \strong{Details} A list summarizing the results. The following information
#' are showed. \code{WgtResponse}, the weight for the response variable in
#' estimating WAASB, \code{WgtWAAS} the weight for stability, \code{Ngen} the
#' number of genotypes, \code{Nenv} the number of environments, \code{OVmean}
#' the overall mean, \code{Min} the minimum observed (returning the genotype and
#' environment), \code{Max} the maximum observed, \code{Max} the maximum
#' observed, \code{MinENV} the environment with the lower mean, \code{MaxENV}
#' the environment with the larger mean observed, \code{MinGEN} the genotype
#' with the lower mean, \code{MaxGEN} the genotype with the larger.
#'  * \strong{augment:} Information about each observation in the dataset. This
#'  includes predicted values in the \code{fitted} column, residuals in the
#'  \code{resid} column, standardized residuals in the \code{stdres} column,
#'  the diagonal of the 'hat' matrix in the \code{hat}, and standard errors for
#'  the fitted values in the \code{se.fit} column.
#' * \strong{probint} The p-value for the genotype-vs-environment interaction.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso \code{\link{waas_means}} \code{\link{waasb}} \code{\link{get_model_data}}
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' #===============================================================#
#' # Example 1: Analyzing all numeric variables considering p-value#
#' # <= 0.05 to compute the WAAS.                                  #
#' #===============================================================#
#'model <- waas(data_ge,
#'              env = ENV,
#'              gen = GEN,
#'              rep = REP,
#'              resp = everything())
#' # Residual plot (first variable)
#' plot(model)
#'
#' # Get the WAAS index
#' get_model_data(model, "WAAS")
#'
#' # Plot WAAS and response variable
#' plot_scores(model, type = 3)
#'
#'
#' #===============================================================#
#' # Example 2: Declaring the number of axis to be used for        #
#' # computing WAAS and assigning a larger weight for the response #
#' # variable when computing the WAASBY index.                     #
#' #===============================================================#
#'
#' model2 <- waas(data_ge,
#'                env = ENV,
#'                gen = GEN,
#'                rep = REP,
#'                resp = everything(),
#'                naxis = 1, # Only to compare with PC1
#'                wresp = 60)
#' # Get the WAAS index (it will be |PC1|)
#' get_model_data(model2)
#'
#' # Get values for IPCA1
#' get_model_data(model2, "PC1")
#'
#'
#' #===============================================================#
#' # Example 3: Analyzing GY and HM assuming a random-effect model.#
#' # Smaller values for HM and higher values for GY are better.    #
#' # To estimate WAASBY, higher weight for the GY (60%) and lower  #
#' # weight for HM (40%) are considered for mean performance.      #
#' #===============================================================#
#'
#' model3 <- waas(data_ge,
#'                env = ENV,
#'                gen = GEN,
#'                rep = REP,
#'                resp = c(GY, HM),
#'                mresp = c(100, 0),
#'                wresp = c(60, 40))
#'
#'
#' # Get the ranks for the WAASY index
#' get_model_data(model3, what = "OrWAASY")
#'}
#'
waas <- function(.data,
                 env,
                 gen,
                 rep,
                 resp,
                 block = NULL,
                 mresp = NULL,
                 wresp = NULL,
                 prob = 0.05,
                 naxis = NULL,
                 ind_anova = FALSE,
                 verbose = TRUE) {
    if(!missing(block)){
        factors  <- .data %>%
            select({{env}},
                   {{gen}},
                   {{rep}},
                   {{block}}) %>%
            mutate(across(everything(), as.factor))
    } else{
        factors  <- .data %>%
            select({{env}},
                   {{gen}},
                   {{rep}}) %>%
            mutate(across(everything(), as.factor))
    }
    vars <- .data %>% select({{resp}}, -names(factors))
    vars %<>% select_numeric_cols()
    if(!missing(block)){
        factors %<>% set_names("ENV", "GEN", "REP", "BLOCK")
    } else{
        factors %<>% set_names("ENV", "GEN", "REP")
    }
    nvar <- ncol(vars)
    if (!is.null(naxis)) {
        if (length(naxis) != nvar) {
            warning("Invalid length in 'naxis'. Setting naxis = ", naxis[[1]],
                    " to all the ", nvar, " variables.", call. = FALSE)
            naxis <- replicate(nvar, naxis[[1]])
        }
    }
    if (is.null(mresp)) {
        mresp <- replicate(nvar, 100)
        minresp <- 100 - mresp
    } else {
        mresp <- mresp
        minresp <- 100 - mresp
        if (length(mresp) != nvar) {
            warning("Invalid length in 'mresp'. Setting mresp = ", mresp[[1]],
                    " to all the ", nvar, " variables.", call. = FALSE)
            mresp <- replicate(nvar, mresp[[1]])
            minresp <- 100 - mresp
        }
        if (sum(mresp == 100) + sum(mresp == 0) != nvar) {
            stop("The values of the numeric vector 'mresp' must be 0 or 100.")
        }
    }
    if (is.null(wresp)) {
        PesoResp <- replicate(nvar, 50)
        PesoWAASB <- 100 - PesoResp
    } else {
        PesoResp <- wresp
        PesoWAASB <- 100 - PesoResp
        if (length(wresp) != nvar) {

            warning("Invalid length in 'wresp'. Setting wresp = ", wresp[[1]],
                    " to all the ", nvar, " variables.", call. = FALSE)
            PesoResp <- replicate(nvar, wresp[[1]])
            PesoWAASB <- 100 - PesoResp
        }
        if (min(wresp) < 0 | max(wresp) > 100) {
            stop("The range of the numeric vector 'wresp' must be equal between 0 and 100.")
        }
    }




    listres <- list()
    vin <- 0
    for (var in 1:nvar) {
        data <- factors %>%
            mutate(Y = vars[[var]])
        check_labels(data)
        if(has_na(data)){
            data <- remove_rows_na(data)
            has_text_in_num(data)
        }
        Nenv <- length(unique(data$ENV))
        Ngen <- length(unique(data$GEN))
        minimo <- min(Nenv, Ngen) - 1
        vin <- vin + 1
        if(ind_anova == TRUE){
            individual <- data %>% anova_ind(ENV, GEN, REP, Y)
        } else{
            individual = NULL
        }
        if(missing(block)){
            model <- performs_ammi(data, ENV, GEN, REP, Y, verbose = FALSE)[[1]]
        } else{
            model <- performs_ammi(data, ENV, GEN, REP, Y, block = BLOCK, verbose = FALSE)[[1]]
        }
        PC <- model$PCA
        Escores <- model$model
        MeansGxE <- model$MeansGxE
        if (is.null(naxis)) {
            SigPC1 <- nrow(PC[which(PC[, 6] < prob), ])
        } else {
            SigPC1 <- naxis[vin]
        }
        if (SigPC1 > minimo) {
            stop("The number of axis to be used must be lesser than or equal to ",
                 minimo, " [min(GEN-1;ENV-1)]")
        } else {
            Pesos <- slice(model$PCA[7], 1:SigPC1)
            WAAS <-
                Escores %>%
                select(contains("PC")) %>%
                abs() %>%
                t() %>%
                as.data.frame() %>%
                slice(1:SigPC1) %>%
                mutate(Percent = Pesos$Proportion)
            WAASAbs <- mutate(Escores, WAAS = sapply(WAAS[, -ncol(WAAS)], weighted.mean, w = WAAS$Percent))
            if (nvar > 1) {
                WAASAbs %<>%
                    group_by(type) %>%
                    mutate(PctResp = (mresp[vin] - minresp[vin])/(max(Y) - min(Y)) * (Y - max(Y)) + mresp[vin],
                           PctWAAS = (0 - 100)/(max(WAAS) - min(WAAS)) * (WAAS - max(WAAS)) + 0,
                           wRes = PesoResp[vin],
                           wWAAS = PesoWAASB[vin],
                           OrResp = rank(-Y),
                           OrWAAS = rank(WAAS),
                           OrPC1 = rank(abs(PC1)),
                           WAASY = ((PctResp * wRes) + (PctWAAS * wWAAS))/(wRes + wWAAS),
                           OrWAASY = rank(-WAASY)) %>%
                    ungroup()
            } else {
                WAASAbs %<>%
                    group_by(type) %>%
                    mutate(PctResp = (mresp - minresp)/(max(Y) - min(Y)) * (Y - max(Y)) + mresp,
                           PctWAAS = (0 - 100)/(max(WAAS) - min(WAAS)) * (WAAS - max(WAAS)) + 0,
                           wRes = PesoResp,
                           wWAAS = PesoWAASB,
                           OrResp = rank(-Y),
                           OrWAAS = rank(WAAS),
                           OrPC1 = rank(abs(PC1)),
                           WAASY = ((PctResp * wRes) + (PctWAAS * wWAAS))/(wRes + wWAAS),
                           OrWAASY = rank(-WAASY)) %>%
                    ungroup()
            }
            min_group <- Escores %>% group_by(type) %>% top_n(1, -Y) %>% select(type, Code, Y) %>% slice(1) %>% as.data.frame()
            max_group <- Escores %>% group_by(type) %>% top_n(1, Y) %>% select(type, Code, Y) %>% slice(1) %>% as.data.frame()
            min <- MeansGxE %>% top_n(1, -Y) %>% select(ENV, GEN, Y) %>% slice(1)
            max <- MeansGxE %>% top_n(1, Y) %>% select(ENV, GEN, Y) %>% slice(1)
            Details <- tibble(Parameters = c("Ngen", "Nenv", "OVmean","Min", "Max", "MinENV", "MaxENV", "MinGEN", "MaxGEN", "SigPC"),
                              Values = c(Ngen, Nenv, round(mean(MeansGxE$Y), 4),
                                         paste0(round(min[3], 4), " (", min$GEN, " in ", min$ENV,")"),
                                         paste0(round(max$Y, 4), " (", max$GEN, " in ", max$ENV,")"),
                                         paste0(min_group[1,2], " (", round(min_group[1,3], 3),")"),
                                         paste0(max_group[1,2], " (", round(max_group[1,3], 3),")"),
                                         paste0(min_group[2,2], " (", round(min_group[2,3], 3), ") "),
                                         paste0(max_group[2,2], " (", round(max_group[2,3], 3), ") "),
                                         SigPC1))
            listres[[paste(names(vars[var]))]] <-
                structure(list(individual = individual[[1]],
                               model = WAASAbs,
                               MeansGxE = MeansGxE,
                               PCA = as_tibble(PC, rownames = NA),
                               anova = model$ANOVA,
                               Details = Details,
                               augment = as_tibble(model$augment, rownames = NA),
                               probint = model$probint),
                          class = "waas")

            if (verbose == TRUE) {
                cat("variable", paste(names(vars[var])),"\n")
                cat("---------------------------------------------------------------------------\n")
                cat("AMMI analysis table\n")
                cat("---------------------------------------------------------------------------\n")
                print(as.data.frame(model$ANOVA), digits = 3, row.names = FALSE)
                cat("---------------------------------------------------------------------------\n\n")
            }
        }
    }
    if (verbose == TRUE) {
        if (length(which(unlist(lapply(listres, function(x) {
            x[["probint"]]
        })) > prob)) > 0) {
            cat("------------------------------------------------------------\n")
            cat("Variables with nonsignificant GxE interaction\n")
            cat(names(which(unlist(lapply(listres, function(x) {
                x[["probint"]]
            })) > prob)), "\n")
            cat("------------------------------------------------------------\n")
        } else {
            cat("All variables with significant (p < 0.05) genotype-vs-environment interaction\n")
        }
        cat("Done!\n")
    }
    invisible(structure(listres, class = "waas"))
}



#' Several types of residual plots
#'
#' Residual plots for a output model of class \code{waas}. Seven types
#' of plots are produced: (1) Residuals vs fitted, (2) normal Q-Q plot for the
#' residuals, (3) scale-location plot (standardized residuals vs Fitted Values),
#' (4) standardized residuals vs Factor-levels, (5) Histogram of raw residuals
#' and (6) standardized residuals vs observation order, and (7) 1:1 line plot.
#'
#'
#' @param x An object of class \code{waas}.
#' @param ... Additional arguments passed on to the function
#'   \code{\link{residual_plots}}
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method plot waas
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' model <- waas(data_ge, ENV, GEN, REP, GY)
#' plot(model)
#' plot(model,
#'      which = c(3, 5),
#'      nrow = 2,
#'      labels = TRUE,
#'      size.lab.out = 4)
#' }
#'
plot.waas <- function(x, ...) {
    residual_plots(x,  ...)
}


#' Print an object of class waas
#'
#' Print the \code{waas} object in two ways. By default, the results are shown
#' in the R console. The results can also be exported to the directory.
#'
#'
#' @param x An object of class \code{waas}.
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported to
#'   the working directory
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print waas
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' model <- waas(data_ge,
#'   resp = c(GY, HM),
#'   gen = GEN,
#'   env = ENV,
#'   rep = REP
#' )
#' print(model)
#' }
print.waas <- function(x, export = FALSE, file.name = NULL, digits = 4, ...) {
    if (!class(x) == "waas") {
        stop("The object must be of class 'waas'")
    }
    if (export == TRUE) {
        file.name <- ifelse(is.null(file.name) == TRUE, "waas print", file.name)
        sink(paste0(file.name, ".txt"))
    }
    opar <- options(pillar.sigfig = digits)
    on.exit(options(opar))
    for (i in 1:length(x)) {
        var <- x[[i]]
        cat("Variable", names(x)[i], "\n")
        cat("---------------------------------------------------------------------------\n")
        cat("Individual analysis of variance\n")
        cat("---------------------------------------------------------------------------\n")
        print(var$individual$individual, ...)
        cat("---------------------------------------------------------------------------\n")
        cat("AMMI analysis table\n")
        cat("---------------------------------------------------------------------------\n")
        print(var$anova, ...)
        cat("---------------------------------------------------------------------------\n")
        cat("Weighted average of the absolute scores\n")
        cat("---------------------------------------------------------------------------\n")
        print(var$model, ...)
        cat("---------------------------------------------------------------------------\n")
        cat("Some information regarding the analysis\n")
        cat("---------------------------------------------------------------------------\n")
        print(var$Details, ...)
        cat("\n\n\n")
    }
    if (export == TRUE) {
        sink()
    }
}







#' Predict the means of a waas object
#'
#' Predict the means of a waas object considering a specific number of axis.
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
#' Predictive sucess (such as cross-validation) should be used to do this. This
#' package provide both. \code{\link{waas}} function compute traditional AMMI
#' analysis showing the number of significant axis. On the other hand,
#' \code{\link{cv_ammif}} function provide a cross-validation, estimating the
#' RMSPD of all AMMI-family models, based on resampling procedures.
#'
#' @param object An object of class waas
#' @param naxis The the number of axis to be use in the prediction. If
#'   \code{object} has more than one variable, then \code{naxis} must be a
#'   vector.
#' @param ... Additional parameter for the function
#' @return A list where each element is the predicted values by the AMMI model
#'   for each variable.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method predict waas
#' @export
#' @examples
#' \donttest{
#' library(metan)
#'model <- waas(data_ge,
#'              env = ENV,
#'              gen = GEN,
#'              rep = REP,
#'              resp = c(GY, HM))
#' # Predict GY with 3 IPCA and HM with 1 IPCA
#' predict <- predict(model, naxis = c(3, 1))
#' predict
#' }
#'
predict.waas <- function(object, naxis = 2, ...) {
    cal <- match.call()
    if (class(object) != "waas") {
        stop("The objectin must be an objectin of the class 'waas'")
    }
    if (length(object) != length(naxis)) {
        warning("Invalid length in 'naxis'. Setting 'mresp = ", naxis[[1]],
                "' to all the ", length(object), " variables.", call. = FALSE)
        naxis <- replicate(length(object), naxis[[1]])
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
                 minimo, ".", call. = FALSE)
        } else {
            if (naxis[var] == 0) {
                warning("Predicted values of AMMI0 model are in colum 'pred_ols'.", call. = FALSE)
            }
            ovmean <- mean(MEDIAS$Y)
            x1 <- model.matrix(~factor(MEDIAS$ENV) - 1)
            z1 <- model.matrix(~factor(MEDIAS$GEN) - 1)
            modelo1 <- lm(Y ~ ENV + GEN, data = MEDIAS)
            MEDIAS <- mutate(MEDIAS, res_ols = residuals(modelo1))
            intmatrix <- t(matrix(MEDIAS$res_ols, Nenv, byrow = T))
            s <- svd(intmatrix)
            U <- s$u[, 1:naxis[var]]
            LL <- s$d[1:naxis[var]]
            V <- s$v[, 1:naxis[var]]
            temp <-
                MEDIAS %>%
                add_cols(pred_ols = Y - res_ols,
                         res_ammi = ((z1 %*% U) * (x1 %*% V)) %*% LL,
                         pred_ammi = pred_ols + res_ammi)
            listres[[paste(names(object[var]))]] <- temp

        }
    }
    invisible(listres)
}
