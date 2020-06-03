#' Cross-validation procedure
#'
#' Cross-validation for estimation of AMMI models
#'
#' THe original dataset is split into two datasets: training set and validation
#' set. The 'training' set has all combinations (genotype x environment) with
#' N-1 replications. The 'validation' set has the remaining replication. The
#' splitting of the dataset into modeling and validation sets depends on the
#' design informed. For Completely Randomized Block Design (default), and
#' alpha-lattice design (declaring \code{block} arguments), complete replicates
#' are selected within environments. The remained replicate serves as validation
#' data. If \code{design = 'RCD'} is informed, completely randomly samples are
#' made for each genotype-by-environment combination (Olivoto et al. 2019). The
#' estimated values considering \code{naxis}-Interaction Principal Component
#' Axis are compared with the 'validation' data. The Root Mean Square Prediction
#' Difference (RMSPD) is computed. At the end of boots, a list is returned.
#'
#' \strong{IMPORTANT:} If the data set is unbalanced (i.e., any genotype missing
#' in any environment) the function will return an error. An error is also
#' observed if any combination of genotype-environment has a different number of
#' replications than observed in the trial.
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#'   replications/blocks. \strong{AT LEAST THREE REPLICATES ARE REQUIRED TO
#'   PERFORM THE CROSS-VALIDATION}.
#' @param resp The response variable.
#' @param block Defaults to \code{NULL}. In this case, a randomized complete
#'   block design is considered. If block is informed, then a resolvable
#'   alpha-lattice design (Patterson and Williams, 1976) is employed.
#'   \strong{All effects, except the error, are assumed to be fixed.}
#' @param nboot The number of resamples to be used in the cross-validation.
#'   Defaults to 200.
#' @param naxis The number of axis to be considered for estimation of GE
#'   effects.
#' @param design The experimental design. Defaults to \code{RCBD} (Randomized
#'   complete Block Design). For Completely Randomized Designs inform
#'   \code{design = 'CRD'}.
#' @param verbose A logical argument to define if a progress bar is shown.
#'   Default is \code{TRUE}.
#' @references Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, V.S. Marchioro,
#'   V.Q. de Souza, and E. Jost. 2019. Mean performance and stability in
#'   multi-environment trials I: Combining features of AMMI and BLUP techniques.
#'   Agron. J. 111:2949-2960.
#' \href{https://acsess.onlinelibrary.wiley.com/doi/abs/10.2134/agronj2019.03.0220}{doi:10.2134/agronj2019.03.0220}
#'
#' @references Patterson, H.D., and E.R. Williams. 1976. A new class of
#' resolvable incomplete block designs. Biometrika 63:83-92.
#' \href{https://doi.org/10.1093/biomet/63.1.83}{doi:10.1093/biomet/63.1.83}
#'
#' @return An object of class \code{cv_ammi} with the following items: *
#' \strong{RMSPD}: A vector with nboot-estimates of the Root Mean Squared
#' Prediction Difference between predicted and validating data.
#'
#' * \strong{RMSPDmean}: The mean of RMSPDmean estimates.
#'
#' * \strong{Estimated}: A data frame that contain the values (predicted,
#' observed, validation) of the last loop.
#'
#' * \strong{Modeling}: The dataset used as modeling data in the last loop
#'
#' * \strong{Testing}: The dataset used as testing data in the last loop.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso \code{\link{cv_ammif}, \link{cv_blup}}
#' @importFrom progress progress_bar
#' @export
#' @examples
#'
#' \donttest{
#' library(metan)
#' model <- cv_ammi(data_ge,
#'                 env = ENV,
#'                 gen = GEN,
#'                 rep = REP,
#'                 resp = GY,
#'                 nboot = 5,
#'                 naxis = 2)
#' }
#'
#'
cv_ammi <- function(.data, env, gen, rep, resp, block = NULL, naxis = 2, nboot = 200, design = "RCBD", verbose = TRUE) {
  if(missing(block)){
    if (!design %in% c("RCBD", "CRD")) {
      stop("Incorrect experimental design informed! Plesease inform RCBD for randomized complete block or CRD for completely randomized design.")
    }
    data <- .data %>%
      dplyr::select(ENV = {{env}},
                    GEN = {{gen}},
                    REP = {{rep}},
                    Y = {{resp}})%>%
      mutate_at(1:3, as.factor)
    RMSPDres <- data.frame(RMSPD = matrix(0, nboot, 1))
    data <- tibble::rowid_to_column(data)
    Nenv <- length(unique(data$ENV))
    Ngen <- length(unique(data$GEN))
    Nbloc <- length(unique(data$REP))
    if (Nbloc <= 2) {
      stop("At least three replicates are required to perform the cross-validation.")
    }
    nrepval <- Nbloc - 1
    minimo <- min(Nenv, Ngen) - 1
    if (naxis > minimo) {
      stop("The number of axis to be used must be lesser than or equal to ",
           minimo, " [min(GEN-1;ENV-1)]")
    }
    if (verbose == TRUE) {
      pb <- progress_bar$new(
        format = "Validating :current of :total sets [:bar] :percent (:elapsedfull -:eta left)",
        clear = FALSE, total = nboot, width = 80)
    }
    test <-
      data %>%
      n_by(ENV, GEN) %>%
      mutate(test =
               case_when(Y > 0 & Y != Nbloc ~ TRUE,
                         TRUE ~ FALSE)
      )
    if(any(test$test == TRUE)){
      df_test <-
        data.frame(
          test[which(test$test == TRUE),] %>%
            remove_cols(rowid,REP, test) %>%
            rename(n = Y)
        )
      message(paste0(capture.output(df_test), collapse = "\n"))
      stop("Combinations of genotype and environment with different number of replication than observed in the trial (", Nbloc, ")", call. = FALSE)
    }
    if(!is_balanced_trial(data, ENV, GEN, Y)){
      stop("AMMI analysis cannot be computed with unbalanced data.", call. = FALSE)
    }
    for (b in 1:nboot) {
      if (design == "CRD") {
        X <- sample(1:10000, 1)
        set.seed(X)
        modeling <- data %>%
          dplyr::group_by(ENV, GEN) %>%
          dplyr::sample_n(nrepval, replace = FALSE) %>%
          arrange(rowid) %>%
          as.data.frame()
        rownames(modeling) <- modeling$rowid
      }
      if (design == "RCBD") {
        tmp <- split_factors(data, ENV, keep_factors = TRUE, verbose = FALSE)
        modeling <- do.call(rbind, lapply(tmp, function(x) {
          X2 <- sample(unique(data$REP), nrepval, replace = FALSE)
          x %>% dplyr::group_by(GEN) %>% dplyr::filter(REP %in% X2)
        })) %>% as.data.frame()
        rownames(modeling) <- modeling$rowid
      }
      testing <- anti_join(data, modeling, by = c("ENV", "GEN", "REP", "Y", "rowid")) %>%
        arrange(ENV, GEN) %>%
        as.data.frame()
      MEDIAS <-
        modeling %>%
        means_by(ENV, GEN) %>%
        as.data.frame()
      residual <-
        modeling %>%
        means_by(ENV, GEN) %>%
        ungroup() %>%
        mutate(residuals = residuals(lm(Y ~ ENV + GEN, data = .))) %>%
        pull(residuals)
      s <- svd(t(matrix(residual, Nenv, byrow = T)))
      MGEN <- model.matrix(~factor(testing$GEN) - 1)
      MENV <- model.matrix(~factor(testing$ENV) - 1)
      MEDIAS %<>% mutate(Ypred = Y - residual,
                         ResAMMI = ((MGEN %*% s$u[, 1:naxis]) * (MENV %*% s$v[, 1:naxis])) %*% s$d[1:naxis],
                         YpredAMMI = Ypred + ResAMMI,
                         testing = testing$Y,
                         error = YpredAMMI - testing,
                         errrorAMMI0 = Ypred - testing)
      RMSPDres[, 1][b] <- ifelse(naxis == 0, sqrt(sum(MEDIAS$errrorAMMI0^2)/length(MEDIAS$errrorAMMI0)),
                                 sqrt(sum(MEDIAS$error^2)/length(MEDIAS$error)))
      if (verbose == TRUE) {
        pb$tick()
      }
    }
    RMSPDmean <- summarise(RMSPDres,
                           MODEL = paste("AMMI", naxis, sep = ""),
                           mean = mean(RMSPD),
                           sd = sd(RMSPD),
                           se = sd(RMSPD)/sqrt(n()),
                           Q2.5 = quantile(RMSPD, 0.025),
                           Q97.5 = quantile(RMSPD, 0.975))
    RMSPDres <- RMSPDres %>%
      mutate(MODEL = paste("AMMI", naxis, sep = "")) %>%
      dplyr::select(MODEL, everything())
  }

  if(!missing(block)){
    data <- .data %>%
      dplyr::select(ENV = {{env}},
                    GEN = {{gen}},
                    REP = {{rep}},
                    BLOCK = {{block}},
                    Y = {{resp}})
    RMSPDres <- data.frame(RMSPD = matrix(0, nboot, 1))
    data <- tibble::rowid_to_column(data)
    Nenv <- nlevels(data$ENV)
    Ngen <- nlevels(data$GEN)
    Nbloc <- nlevels(data$REP)
    if (Nbloc <= 2) {
      stop("At least three replicates are required to perform the cross-validation.")
    }
    nrepval <- Nbloc - 1
    minimo <- min(Nenv, Ngen) - 1
    if (naxis > minimo) {
      stop("The number of axis to be used must be lesser than or equal to ",
           minimo, " [min(GEN-1;ENV-1)]")
    }
    test <-
      data %>%
      n_by(ENV, GEN) %>%
      mutate(test =
               case_when(Y > 0 & Y != Nbloc ~ TRUE,
                         TRUE ~ FALSE)
      )
    if(any(test$test == TRUE)){
      df_test <-
        data.frame(
          test[which(test$test == TRUE),] %>%
            remove_cols(rowid,REP, test) %>%
            rename(n = Y)
        )
      message(paste0(capture.output(df_test), collapse = "\n"))
      stop("Combinations of genotype and environment with different number of replication than observed in the trial (", Nbloc, ")", call. = FALSE)
    }
    if(!is_balanced_trial(data, ENV, GEN, Y)){
      stop("AMMI analysis cannot be computed with unbalanced data.", call. = FALSE)
    }
    if (verbose == TRUE) {
      pb <- progress_bar$new(
        format = "Validating :current of :total sets [:bar] :percent (:elapsedfull -:eta left)",
        clear = FALSE, total = nboot, width = 80)
    }
    for (b in 1:nboot) {
      tmp <- split_factors(data, ENV, keep_factors = TRUE, verbose = FALSE)
      modeling <- do.call(rbind, lapply(tmp, function(x) {
        X2 <- sample(unique(data$REP), nrepval, replace = FALSE)
        x %>% dplyr::group_by(GEN) %>% dplyr::filter(REP %in% X2)
      })) %>% as.data.frame()
      rownames(modeling) <- modeling$rowid
      testing <- anti_join(data, modeling, by = c("ENV", "GEN", "REP", "BLOCK", "Y", "rowid")) %>%
        arrange(ENV, GEN) %>%
        as.data.frame()
      MEDIAS <- modeling %>%
        group_by(ENV, GEN) %>%
        summarise(Y = mean(Y)) %>%
        as.data.frame()
      residual <- modeling %>%
        mutate(residuals = residuals(lm(Y ~ ENV/REP + REP:BLOCK + ENV +  GEN, data = .))) %>%
        group_by(ENV, GEN) %>%
        summarise(residuals = mean(residuals)) %>%
        ungroup() %>%
        pull(residuals)
      s <- svd(t(matrix(residual, Nenv, byrow = T)))
      MGEN <- model.matrix(~factor(testing$GEN) - 1)
      MENV <- model.matrix(~factor(testing$ENV) - 1)
      MEDIAS %<>% mutate(Ypred = Y - residual,
                         ResAMMI = ((MGEN %*% s$u[, 1:naxis]) * (MENV %*% s$v[, 1:naxis])) %*% s$d[1:naxis],
                         YpredAMMI = Ypred + ResAMMI,
                         testing = testing$Y,
                         error = YpredAMMI - testing,
                         errrorAMMI0 = Ypred - testing)
      RMSPDres[, 1][b] <- ifelse(naxis == 0, sqrt(sum(MEDIAS$errrorAMMI0^2)/length(MEDIAS$errrorAMMI0)),
                                 sqrt(sum(MEDIAS$error^2)/length(MEDIAS$error)))
      if (verbose == TRUE) {
        pb$tick()
      }
    }
    RMSPDmean <- summarise(RMSPDres,
                           MODEL = paste("a_AMMI", naxis, sep = ""),
                           mean = mean(RMSPD),
                           sd = sd(RMSPD),
                           se = sd(RMSPD)/sqrt(n()),
                           Q2.5 = quantile(RMSPD, 0.025),
                           Q97.5 = quantile(RMSPD, 0.975))
    RMSPDres <- RMSPDres %>%
      mutate(MODEL = paste("a_AMMI", naxis, sep = "")) %>%
      dplyr::select(MODEL, everything())
  }

  return(structure(list(RMSPD = RMSPDres, RMSPDmean = RMSPDmean,
                        Estimated = MEDIAS, Modeling = modeling, Testing = testing),
                   class = "cvalidation"))
}








#' Plot the RMSPD of a cross-validation procedure
#'
#' Boxplot showing the Root Means Square Prediction Difference of of a cross
#' validation procedure.
#'
#' Five statistics are shown in this type of plot. The lower and upper hinges
#' correspond to the first and third quartiles (the 25th and 75th percentiles).
#' The upper whisker extends from the hinge to the largest value no further than
#' 1.5 * IQR from the hinge (where IQR is the inter-quartile range). The lower
#' whisker extends from the hinge to the smallest value at most 1.5 * IQR of the
#' hinge. Data beyond the end of the whiskers are considered outlying points.
#'
#' @param x An object of class \code{cvalidation} fitted with the functions
#'   \code{\link{cv_ammi}}, \code{\link{cv_ammif}}, \code{\link{cv_blup}}, or a
#'   bound object fitted with \code{\link{bind_cv}}.
#' @param violin Define if a violin plot is used with boxplot. Default is 'TRUE'
#' @param export Export (or not) the plot. Default is \code{T}.
#' @param order_box Logical argument. If \code{TRUE} then the boxplots will be
#'   ordered according to the values of the RMSPD.
#' @param x.lab The label of x-axis. New arguments can be inserted as
#'   \code{x.lab = 'my x label'}.
#' @param y.lab The label of y-axis. New arguments can be inserted as
#'   \code{y.lab = 'my y label'}.
#' @param size.tex.lab The size of the text in axis text and labels.
#' @param file.type The type of file to be exported. Default is \code{pdf},
#'   Graphic can also be exported in \code{*.tiff} format by declaring
#'   \code{file.type = 'tiff'}.
#' @param file.name The name of the file for exportation, default is
#'   \code{NULL}, i.e. the files are automatically named.
#' @param plot_theme The graphical theme of the plot. Default is
#'   \code{plot_theme = theme_metan()}. For more details,see
#'   \code{\link[ggplot2]{theme}}.
#' @param width The width 'inch' of the plot. Default is \code{6}.
#' @param height The height 'inch' of the plot. Default is \code{6}.
#' @param resolution The resolution of the plot. Parameter valid if
#'   \code{file.type = 'tiff'} is used. Default is \code{300} (300 dpi)
#' @param col.violin Parameter valid if \code{violin = T}. Define the color of
#'   the violin plot. Default is 'gray90.
#' @param col.boxplot Define the color for boxplot. Default is 'gray70'.
#' @param col.boxplot.win Define the color for boxplot of the best model.
#'   Default is 'cyan'.
#' @param width.boxplot The width of boxplots. Default is \code{0.2}.
#' @param x.lim The range of x-axis. Default is \code{NULL} (maximum and minimum
#'   values of the data set). New arguments can be inserted as \code{x.lim =
#'   c(x.min, x.max)}.
#' @param x.breaks The breaks to be plotted in the x-axis. Default is
#'   \code{authomatic breaks}. New arguments can be inserted as \code{x.breaks =
#'   c(breaks)}
#' @param ... Currently not used.
#' @return An object of class \code{gg, ggplot}.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method plot cvalidation
#' @export
#' @examples
#'
#' \donttest{
#' validation <- cv_ammif(data_ge,
#'                        resp = GY,
#'                        gen = GEN,
#'                        env = ENV,
#'                        rep = REP,
#'                        nboot = 5)
#' plot(validation)
#' }
#'
#'
plot.cvalidation <- function(x, violin = FALSE, export = FALSE, order_box =  FALSE,
                             x.lab = NULL, y.lab = NULL, size.tex.lab = 12, file.type = "pdf",
                             file.name = NULL, plot_theme = theme_metan(), width = 6, height = 6,
                             resolution = 300, col.violin = "gray90", col.boxplot = "gray70",
                             col.boxplot.win = "cyan", width.boxplot = 0.6, x.lim = NULL,
                             x.breaks = waiver(), ...) {
  if (!class(x) == "cvalidation") {
    stop("The object 'x' must be of class 'cvalidation'.")
  }
  y.lab <- ifelse(missing(y.lab),
                  expression(paste("Root mean square prediction difference (Mg ha"^-1, ")")),
                  y.lab)
  x.lab <- ifelse(missing(x.lab), "Model", x.lab)
  a <- suppressMessages(x$RMSPD %>%
                          group_by(MODEL) %>%
                          summarise(RMSPD = mean(RMSPD)) %>%
                          top_n(-1))
  mod <- paste(a$MODEL[1])
  if (violin == TRUE) {
    dodge <- position_dodge(width = 1)
    if(order_box == TRUE){
      p1 <- ggplot(x$RMSPD, aes(x = reorder(MODEL, RMSPD), y = RMSPD))
    } else {
      p1 <- ggplot(x$RMSPD, aes(x = MODEL, y = RMSPD))
    }
    p1 = p1 +
      geom_violin(position = dodge, fill = col.violin) +
      stat_boxplot(geom = "errorbar", width=width.boxplot/2, na.rm = TRUE)+
      geom_boxplot(width = width.boxplot,
                   position = dodge,
                   color = "black",
                   outlier.alpha = 0.75,
                   outlier.color = "black",
                   fill = col.boxplot) +
      geom_boxplot(data = x$RMSPD[x$RMSPD$MODEL == mod, ],
                   aes(x = MODEL, y = RMSPD),
                   width = width.boxplot,
                   position = dodge,
                   color = "black",
                   outlier.alpha = 0.75,
                   outlier.color = "black",
                   fill = col.boxplot.win) +
      stat_summary(fun = mean,
                   geom = "point",
                   shape = 23,
                   fill = "black") +
      plot_theme %+replace%
      theme(axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black")) +
      coord_flip() +
      scale_y_continuous(limits = x.lim, breaks = x.breaks) +
      labs(x = x.lab, y = y.lab)
  } else {
    dodge <- position_dodge(width = 1)
    if(order_box == TRUE){
      p1 <- ggplot(x$RMSPD, aes(x = reorder(MODEL, RMSPD), y = RMSPD))
    } else {
      p1 <- ggplot(x$RMSPD, aes(x = MODEL, y = RMSPD))
    }
    p1 = p1 +
      stat_boxplot(geom = "errorbar", width=width.boxplot/2, na.rm = TRUE)+
      geom_boxplot(width = width.boxplot,
                   position = dodge,
                   color = "black",
                   outlier.alpha = 0.75,
                   outlier.color = "black",
                   fill = col.boxplot) +
      geom_boxplot(data = x$RMSPD[x$RMSPD$MODEL == mod, ],
                   aes(x = MODEL, y = RMSPD),
                   width = width.boxplot,
                   position = dodge,
                   outlier.alpha = 0.75,
                   color = "black",
                   outlier.color = "black",
                   fill = col.boxplot.win) +
      stat_summary(fun = mean,
                   geom = "point",
                   shape = 23,
                   fill = "black") +
      plot_theme %+replace%
      theme(axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black")) +
      coord_flip() +
      scale_y_continuous(limits = x.lim, breaks = x.breaks) +
      labs(x = x.lab, y = y.lab)
  }
  if (export == FALSE) {
    return(p1)
  } else if (file.type == "pdf") {
    if (is.null(file.name)) {
      pdf("RMSPD validation.pdf", width = width, height = height)
    } else pdf(paste0(file.name, ".pdf"), width = width, height = height)
    plot(p1)
    dev.off()
  }
  if (file.type == "tiff") {
    if (is.null(file.name)) {
      tiff(filename = "RMSPD validation.tiff", width = width,
           height = height, units = "in", compression = "lzw",
           res = resolution)
    } else tiff(filename = paste0(file.name, ".tiff"), width = width,
                height = height, units = "in", compression = "lzw",
                res = resolution)
    plot(p1)
    dev.off()
  }
}
