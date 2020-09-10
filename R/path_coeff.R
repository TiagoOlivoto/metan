#' Path coefficients with minimal multicollinearity
#'
#' Computes direct and indirect effects in path analysis. An algorithm to select
#' a set of predictors with minimal multicollinearity and high explanatory power
#' is implemented.
#'
#' When \code{brutstep = TRUE}, first, the algorithm will select a set of
#' predictors with minimal multicollinearity. The selection is based on the
#' variance inflation factor (VIF). An iterative process is performed until the
#' maximum VIF observed is less than \code{maxvif}. The variables selected in
#' this iterative process are then used in a series of stepwise-based
#' regressions. The first model is fitted and p-1 predictor variables are
#' retained (p is the number of variables selected in the iterative process. The
#' second model adjusts a regression considering p-2 selected variables, and so
#' on until the last model, which considers only two variables. Three objects
#' are created. \code{Summary}, with the process summary, \code{Models},
#' containing the aforementioned values for all the adjusted models; and
#' \code{Selectedpred}, a vector with the name of the selected variables in the
#' iterative process.
#'
#' @param .data The data. Must be a data frame or a grouped data passed from
#'   \code{\link[dplyr]{group_by}()}
#' @param resp The dependent variable.
#' @param by One variable (factor) to compute the function by. It is a shortcut
#'   to \code{\link[dplyr]{group_by}()}. To compute the statistics by more than
#'   one grouping variable use that function.
#' @param pred The predictor variables, set to \code{everything()}, i.e., the
#'   predictor variables are all the numeric variables in the data except that
#'   in \code{resp}.
#' @param exclude Logical argument, set to false. If \code{exclude = TRUE}, then
#'   the variables in \code{pred} are deleted from the data, and the analysis
#'   will use as predictor those that remained, except that in \code{resp}.
#' @param correction Set to \code{NULL}. A correction value (k) that will be
#'   added into the diagonal elements of the \bold{X'X} matrix aiming at
#'   reducing the harmful problems of the multicollinearity in path analysis
#'   (Olivoto et al., 2017)
#' @param knumber When \code{correction = NULL}, a plot showing the values of
#'   direct effects in a set of different k values (0-1) is produced.
#'   \code{knumber} is the number of k values used in the range of 0 to 1.
#' @param brutstep Logical argument, set to \code{FALSE}. If true, then an
#'   algorithm will select a subset of variables with minimal multicollinearity
#'   and fit a set of possible models. See the \bold{Details} section for more
#'   information.
#' @param maxvif The maximum value for the Variance Inflation Factor (cut point)
#'   that will be accepted. See the \bold{Details} section for more information.
#' @param missingval How to deal with missing values. For more information,
#'   please see \code{\link[stats]{cor}()}.
#' @param plot_res If \code{TRUE}, create a scatter plot of residual against
#'   predicted value and a normal Q-Q plot.
#' @param verbose If \code{verbose = TRUE} then some results are shown in the
#'   console.
#' @param ... Additional arguments passed on to \code{\link[stats]{plot.lm}}
#' @return An object of class \code{path_coeff, group_path, or brute_path} with
#'   the following items:
#' * \strong{Corr.x} A correlation matrix between the predictor variables.
#' * \strong{Corr.y} A vector of correlations between each predictor variable with
#' the dependent variable.
#' * \strong{Coefficients} The path coefficients. Direct effects are the diagonal
#' elements, and the indirect effects those in the off-diagonal elements
#' (column)
#' * \strong{Eigen} Eigenvectors and eigenvalues of the \code{Corr.x.}
#' * \strong{VIF} The Variance Inflation Factors.
#' * \strong{plot} A ggplot2-based graphic showing the direct effects in 21
#' different k values.
#' * \strong{Predictors} The predictor variables used in the model.
#' * \strong{CN} The Condition Number, i.e., the ratio between the highest and
#' lowest eigenvalue.
#' * \strong{Det} The matrix determinant of the \code{Corr.x.}.
#' * \strong{R2} The coefficient of determination of the model.
#' * \strong{Residual} The residual effect of the model.
#' * \strong{Response} The response variable.
#' * \strong{weightvar} The order of the predictor variables with the highest weight
#' (highest eigenvector) in the lowest eigenvalue.
#'
#' If \code{.data} is a grouped data passed from \code{\link[dplyr]{group_by}()}
#' then the results will be returned into a list-column of data frames,
#' containing:
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references
#' Olivoto, T., V.Q. Souza, M. Nardino, I.R. Carvalho, M. Ferrari, A.J.
#' Pelegrin, V.J. Szareski, and D. Schmidt. 2017. Multicollinearity in path
#' analysis: a simple method to reduce its effects. Agron. J. 109:131-142.
#' \doi{10.2134/agronj2016.04.0196}
#'
#' @export
#' @examples
#' \donttest{
#' library(metan)
#'
#' # Using KW as the response variable and all other ones as predictors
#' pcoeff <- path_coeff(data_ge2, resp = KW)
#'
#'
#' # Declaring the predictors
#' # Create a residual plot with 'plot_res = TRUE'
#' pcoeff2 <- path_coeff(data_ge2,
#'                       resp = KW,
#'                       pred = c(PH, EH, NKE, TKW),
#'                       plot_res = TRUE)
#'
#'
#' # Selecting variables to be excluded from the analysis
#'pcoeff3 <- path_coeff(data_ge2,
#'                      resp = KW,
#'                      pred = c(NKR, PERK, KW, NKE),
#'                      exclude = TRUE)
#'
#'
#' # Selecting a set of predictors with minimal multicollinearity
#' # Maximum variance Inflation factor of 5
#'pcoeff4 <- path_coeff(data_ge2,
#'                      resp = KW,
#'                      brutstep = TRUE,
#'                      maxvif = 5)
#'
#'
#' # When one analysis should be carried out for each environment
#' # Using the forward-pipe operator %>%
#'pcoeff5 <- path_coeff(data_ge2, resp = KW, by = ENV)
#'}
#'
#'
path_coeff <- function(.data,
                       resp,
                       by = NULL,
                       pred = everything(),
                       exclude = FALSE,
                       correction = NULL,
                       knumber = 50,
                       brutstep = FALSE,
                       maxvif = 10,
                       missingval = "pairwise.complete.obs",
                       plot_res = FALSE,
                       verbose = TRUE,
                       ...) {
  if (missing(resp) == TRUE) {
    stop("A response variable (dependent) must be declared.")
  }
  kincrement <- 1/knumber
  if (!missing(by)){
    if(length(as.list(substitute(by))[-1L]) != 0){
      stop("Only one grouping variable can be used in the argument 'by'.\nUse 'group_by()' to pass '.data' grouped by more than one variable.", call. = FALSE)
    }
    .data <- group_by(.data, {{by}})
  }
  if(is_grouped_df(.data)){
    results <- .data %>%
      doo(path_coeff,
          resp = {{resp}},
          pred = {{pred}},
          exclude = exclude,
          correction = correction,
          knumber = knumber,
          brutstep = brutstep,
          maxvif = maxvif,
          missingval = missingval,
          verbose = verbose)
    return(set_class(results, c("group_path", "tbl_df", "tbl",  "data.frame")))
  }
  data <- select_numeric_cols(.data)
  nam <- names(.data)
  if (brutstep == FALSE) {
    if (exclude == TRUE) {
      pr <- data %>% select(-{{pred}}, -{{resp}})
    } else {
      pr <- data %>% select({{pred}}, -{{resp}})
    }
    names <- colnames(pr)
    y <- data %>% select({{resp}})
    if(plot_res == TRUE){
    dfs <- cbind(y, pr)
    form <- as.formula(paste(names(y), "~ ."))
    mod <- lm(form, data = dfs)
    opar <- par(mfrow = c(1, 2))
    on.exit(par(opar))
    plot(mod, which = c(1, 2), ...)
    }
    nam_resp <- names(y)
    cor.y <- cor(pr, y, use = missingval)
    cor.x <- cor(pr, use = missingval)
    if (is.null(correction) == FALSE) {
      diag(cor.x) <- diag(cor.x) + correction
    } else {
      cor.x <- cor(pr, use = missingval)
    }
    if (is.null(correction) == TRUE) {
      betas <- data.frame(matrix(nrow = knumber, ncol = length(pr) + 1))
      cc <- 0
      nvar <- length(pr) + 1
      for (i in 1:knumber) {
        cor.x2 <- cor.x
        diag(cor.x2) <- diag(cor.x2) + cc
        betas[i, 1] <- cc
        betas[i, 2:nvar] <- t(solve_svd(cor.x2, cor.y))
        cc <- cc + kincrement
      }
      names(betas) <- paste0(c("K", names(pr)))
      xx <- betas$K
      yy <- colnames(betas)
      fila <- knumber
      col <- length(yy)
      total <- fila * (col - 1)
      x <- character(length = total)
      y <- character(length = total)
      z <- numeric(length = total)
      k <- 0
      for (i in 1:fila) {
        for (j in 2:col) {
          k <- k + 1
          x[k] <- xx[i]
          y[k] <- yy[j]
          z[k] <- betas[i, j]
        }
      }
      x <- as.numeric(x)
      betas <- data.frame(K = x, VAR = y, direct = z)
      p1 <- ggplot2::ggplot(betas, ggplot2::aes(K,
                                                direct, col = VAR)) + ggplot2::geom_line(size = 1) +
        ggplot2::theme_bw() + ggplot2::theme(axis.ticks.length = unit(0.2,
                                                                      "cm"), axis.text = element_text(size = 12,
                                                                                                      colour = "black"), axis.title = element_text(size = 12,
                                                                                                                                                   colour = "black"), axis.ticks = element_line(colour = "black"),
                                             legend.position = "bottom", plot.margin = margin(0.1,
                                                                                              0.1, 0.1, 0.1, "cm"), legend.title = element_blank(),
                                             panel.border = element_rect(colour = "black",
                                                                         fill = NA, size = 1), panel.grid.major.x = element_blank(),
                                             panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(),
                                             panel.grid.minor.y = element_blank()) + ggplot2::labs(x = "k values",
                                                                                                   y = "Beta values") + ggplot2::geom_abline(intercept = 0,
                                                                                                                                             slope = 0) + ggplot2::scale_x_continuous(breaks = seq(0,
                                                                                                                                                                                                   1, by = 0.1))
    } else {
      p1 <- "No graphic generated due to correction value"
    }
    eigen <- eigen(cor.x)
    Det <- det(cor.x)
    NC <- max(eigen$values)/min(eigen$values)
    Aval <- data.frame(eigen$values)
    names(Aval) <- "Eigenvalues"
    Avet <- data.frame(t(eigen$vectors))
    names(Avet) <- names
    AvAvet <- cbind(Aval, Avet)
    Direct <- solve(cor.x, cor.y)
    n <- ncol(cor.x)
    Coeff <- data.frame(cor.x)
    for (i in 1:n) {
      for (j in 1:n) {
        Coeff[i, j] <- Direct[j] * cor.x[i, j]
      }
    }
    Residual <- 1 - t(Direct) %*% cor.y
    R2 <- t(Direct) %*% cor.y
    VIF <- data.frame(diag(solve_svd(cor.x)))
    names(VIF) <- "VIF"
    if (verbose == TRUE) {
      if (NC > 1000) {
        cat(paste0("Severe multicollinearity. \n",
                   "Condition Number = ", round(NC, 3), "\n",
                   "Please, consider using a correction factor, or use 'brutstep = TRUE'. \n"))
      }
      if (NC < 100) {
        cat(paste0("Weak multicollinearity. \n", "Condition Number = ",
                   round(NC, 3), "\n", "You will probably have path coefficients close to being unbiased. \n"))
      }
      if (NC > 100 & NC < 1000) {
        cat(paste0("Moderate multicollinearity! \n",
                   "Condition Number = ", round(NC, 3), "\n",
                   "Please, cautiosely evaluate the VIF and matrix determinant.",
                   "\n"))
      }
    }
    last <- data.frame(weight = t(AvAvet[c(nrow(AvAvet)),
                                         ])[-c(1), ])
    abs <- data.frame(weight = abs(last[, "weight"]))
    rownames(abs) <- rownames(last)
    last <- abs[order(abs[, "weight"], decreasing = T),
                , drop = FALSE]
    weightvarname <- paste(rownames(last), collapse = " > ")
    temp <- structure(list(Corr.x = as_tibble(cor.x, rownames = NA),
                           Corr.y = as_tibble(cor.y, rownames = NA),
                           Coefficients = cbind(as_tibble(t(Coeff), rownames = NA), as_tibble(cor.y, rownames = NA) %>% set_names("linear")),
                           Eigen = as_tibble(AvAvet),
                           VIF =  rownames_to_column(VIF, "VAR") %>% as_tibble(),
                           plot = p1,
                           Predictors = names(pr),
                           CN = NC,
                           Det = Det,
                           R2 = R2,
                           Residual = Residual,
                           Response = nam_resp,
                           weightvar = weightvarname),
                      class = "path_coeff")
    return(temp)
  }
  if (brutstep == TRUE) {
    yyy <- data %>% dplyr::select({{resp}}) %>% as.data.frame()
    nam_resp <- names(yyy)
    xxx <- data %>% dplyr::select(-{{resp}}) %>% as.data.frame()
    cor.xx <- cor(xxx, use = missingval)
    VIF <- data.frame(VIF = diag(solve_svd(cor.xx))) %>%
      rownames_to_column("VAR") %>%
      arrange(VIF)
    repeat {
      VIF2 <- slice(VIF, -n())
      xxx2 <- data[VIF2$VAR]
      VIF3 <- data.frame(VIF = diag(solve_svd(cor(xxx2, use = missingval))))%>%
        rownames_to_column("VAR") %>%
        arrange(VIF)
      if (max(VIF3$VIF) <= maxvif)
        break
      VIF <- VIF3
    }
    xxx <- data[VIF3$VAR] %>% as.data.frame()
    selectedpred <- VIF3$VAR
    npred <- ncol(xxx) - 1
    statistics <- data.frame(matrix(nrow = npred - 1,
                                    ncol = 8))
    ModelEstimates <- list()
    modelcode <- 1
    nproced <- npred - 1
    if (verbose == TRUE) {
      cat("--------------------------------------------------------------------------\n")
      cat(paste("The algorithm has selected a set of ",
                nrow(VIF3), " predictors with largest VIF = ",
                round(max(VIF3$VIF), 3), ".", sep = ""), "\n")
      cat("Selected predictors:", paste0(selectedpred),"\n")
      cat(paste("A forward stepwise-based selection procedure will fit", nproced, "models.\n"))
      cat("--------------------------------------------------------------------------\n")
    }
    for (i in 1:nproced) {
      FDSel <- FWDselect::selection(x = xxx, y = unlist(yyy),
                                    q = npred, method = "lm", criterion = "aic",
                                    cluster = FALSE)
      predstw <- FDSel$Variable_names
      x <- data[, c(predstw)]
      names <- colnames(x)
      y <- data %>% dplyr::select({{resp}})
      nam_resp <- names(y)
      cor.y <- cor(x, y, use = missingval)
      cor.x <- cor(x, use = missingval)
      if (is.null(correction) == FALSE) {
        diag(cor.x) <- diag(cor.x) + correction
      } else {
        cor.x <- cor(x, use = missingval)
      }
      if (is.null(correction) == T) {
        betas <- data.frame(matrix(nrow = knumber,
                                   ncol = length(predstw) + 1))
        cc <- 0
        nvar <- length(predstw) + 1
        for (i in 1:knumber) {
          cor.x2 <- cor.x
          diag(cor.x2) <- diag(cor.x2) + cc
          betas[i, 1] <- cc
          betas[i, 2:nvar] <- t(solve_svd(cor.x2, cor.y))
          cc <- cc + kincrement
        }
        names(betas) <- paste0(c("K", predstw))
        xx <- betas$K
        yy <- colnames(betas)
        fila <- knumber
        col <- length(yy)
        total <- fila * (col - 1)
        x <- character(length = total)
        y <- character(length = total)
        z <- numeric(length = total)
        k <- 0
        for (i in 1:fila) {
          for (j in 2:col) {
            k <- k + 1
            x[k] <- xx[i]
            y[k] <- yy[j]
            z[k] <- betas[i, j]
          }
        }
        x <- as.numeric(x)
        betas <- data.frame(K = x, VAR = y, direct = z)
        p1 <- ggplot2::ggplot(betas, ggplot2::aes(K,
                                                  direct, col = VAR)) + ggplot2::geom_line(size = 1) +
          ggplot2::theme_bw() + ggplot2::theme(axis.ticks.length = unit(0.2,
                                                                        "cm"), axis.text = element_text(size = 12,
                                                                                                        colour = "black"), axis.title = element_text(size = 12,
                                                                                                                                                     colour = "black"), axis.ticks = element_line(colour = "black"),
                                               legend.position = "bottom", plot.margin = margin(0.1,
                                                                                                0.1, 0.1, 0.1, "cm"), legend.title = element_blank(),
                                               panel.border = element_rect(colour = "black",
                                                                           fill = NA, size = 1), panel.grid.major.x = element_blank(),
                                               panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(),
                                               panel.grid.minor.y = element_blank()) + ggplot2::labs(x = "k values",
                                                                                                     y = expression(paste(beta, " values"))) +
          ggplot2::geom_abline(intercept = 0, slope = 0) +
          ggplot2::scale_x_continuous(breaks = seq(0,
                                                   1, by = 0.1))
      } else {
        p1 <- "No graphic generated due to correction value"
      }
      eigen <- eigen(cor.x)
      Det <- det(cor.x)
      NC <- max(eigen$values)/min(eigen$values)
      Aval <- as.data.frame(eigen$values)
      names(Aval) <- "Eigenvalues"
      Avet <- as.data.frame(eigen$vectors)
      names(Avet) <- names
      AvAvet <- cbind(Aval, Avet)
      Direct <- solve(cor.x, cor.y)
      n <- ncol(cor.x)
      Coeff <- data.frame(cor.x)
      for (i in 1:n) {
        for (j in 1:n) {
          Coeff[i, j] <- Direct[j] * cor.x[i, j]
        }
      }
      Residual <- 1 - t(Direct) %*% cor.y
      R2 <- t(Direct) %*% cor.y
      VIF <- data.frame(diag(solve_svd(cor.x)))
      names(VIF) <- "VIF"
      last <- data.frame(weight = t(AvAvet[c(nrow(AvAvet)),
                                           ])[-c(1), ])
      abs <- data.frame(weight = abs(last[, "weight"]))
      rownames(abs) <- rownames(last)
      last <- abs[order(abs[, "weight"], decreasing = T),
                  , drop = FALSE]
      weightvarname <- paste(rownames(last), collapse = " > ")
      cor.y <- data.frame(cor.y)
      Results <- list(Corr.x = data.frame(cor.x),
                      Corr.y = data.frame(cor.y),
                      Coefficients = cbind(data.frame(t(Coeff)), data.frame(cor.y) %>% set_names("linear")),
                      Eigen = AvAvet,
                      VIF = VIF,
                      plot = p1,
                      Predictors = predstw,
                      CN = NC,
                      Det = Det,
                      R2 = R2,
                      Residual = Residual,
                      Response = nam_resp,
                      weightvar = weightvarname)
      ModelEstimates[[paste("Model_", modelcode, sep = "")]] <- set_class(Results, "path_coeff")
      statistics[i, 1] <- paste("Model", modelcode, sep = "")
      statistics[i, 2] <- FDSel$Information_Criterion
      statistics[i, 3] <- npred
      statistics[i, 4] <- NC
      statistics[i, 5] <- Det
      statistics[i, 6] <- R2
      statistics[i, 7] <- Residual
      statistics[i, 8] <- max(VIF)
      if (verbose == TRUE) {
        cat(paste("Adjusting the model ", modelcode,
                  " with ", npred, " predictors (", round(modelcode/nproced *
                                                            100, 2), "% concluded)", "\n", sep = ""))
      }
      npred <- npred - 1
      modelcode <- modelcode + 1
    }
    statistics %<>% remove_rows(1) %>%
      set_names("Model", "AIC", "Numpred", "CN", "Determinant", "R2", "Residual", "maxVIF") %>%
      arrange(Model) %>%
      tidy_strings()
    # arrange(statistics, Model) %>% tidy_strings()
    # # statistics <- statistics[-c(1), ]
    # # names(statistics) <- c("Model", "AIC", "Numpred", "CN", "Determinant", "R2", "Residual", "maxVIF")
    if (verbose == TRUE) {
      cat("Done!\n")
      cat("--------------------------------------------------------------------------\n")
      cat("Summary of the adjusted models", "\n")
      cat("--------------------------------------------------------------------------\n")
      print(statistics, digits = 3, row.names = FALSE)
      cat("--------------------------------------------------------------------------\n\n")
    }
    temp <- list(Models = set_class(ModelEstimates, "path_coeff"),
                 Summary = statistics,
                 Selectedpred = selectedpred)
    return(set_class(temp, c("brute_path", "path_coeff")))
  }
}








#' Print an object of class path_coeff
#'
#' Print an object generated by the function 'path_coeff()'. By default, the
#' results are shown in the R console. The results can also be exported to the
#' directory.
#'
#'
#' @param x An object of class \code{path_coeff} or \code{group_path}.
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported to
#'   the working directory
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print path_coeff
#' @export
#' @examples
#' \donttest{
#' library(metan)
#'
#' # KW as dependent trait and all others as predictors
#' pcoeff <- path_coeff(data_ge2, resp = KW)
#' print(pcoeff)
#'
#' # Call the algorithm for selecting a set of predictors
#' # With minimal multicollinearity (no VIF larger than 5)
#' pcoeff2 <- path_coeff(data_ge2,
#'                       resp = KW,
#'                       brutstep = TRUE,
#'                       maxvif = 5)
#' print(pcoeff2)
#' }
print.path_coeff <- function(x, export = FALSE, file.name = NULL, digits = 4, ...) {
  args <- match.call()
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "path_coeff print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  if (!has_class(x, c("path_coeff", "brute_path"))) {
    stop("The object 'x' must be of class 'path_coeff' or 'group_path'.")
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  if (has_class(x, "brute_path")) {
    df <- as_tibble(x[["Summary"]])
    print(df)
    message("Go to '", args[["x"]], " > s' to select a specific model", sep = "")

  } else{
    cat("----------------------------------------------------------------------------------------------\n")
    cat("Correlation matrix between the predictor traits\n")
    cat("----------------------------------------------------------------------------------------------\n")
    print(x$Corr.x)
    cat("----------------------------------------------------------------------------------------------\n")
    cat("Vector of correlations between dependent and each predictor\n")
    cat("----------------------------------------------------------------------------------------------\n")
    print(t(x$Corr.y))
    cat("----------------------------------------------------------------------------------------------\n")
    cat("Multicollinearity diagnosis and goodness-of-fit\n")
    cat("----------------------------------------------------------------------------------------------\n")
    cat("Condition number: ", round(x$CN, 4), "\n")
    cat("Determinant:      ", round(x$Det, 8), "\n")
    cat("R-square:         ", round(x$R2, 4), "\n")
    cat("Residual:         ", round(x$Residual, 4), "\n")
    cat("Response:         ", paste(x$Response), "\n")
    cat("Predictors:       ", paste(x$Predictors), "\n")
    cat("----------------------------------------------------------------------------------------------\n")
    cat("Variance inflation factors\n")
    cat("----------------------------------------------------------------------------------------------\n")
    print(x$VIF)
    cat("----------------------------------------------------------------------------------------------\n")
    cat("Eigenvalues and eigenvectors\n")
    cat("----------------------------------------------------------------------------------------------\n")
    print(x$Eigen)
    cat("----------------------------------------------------------------------------------------------\n")
    cat("Variables with the largest weight in the eigenvalue of smallest magnitude\n")
    cat("----------------------------------------------------------------------------------------------\n")
    cat(x$weightvar, "\n")
    cat("----------------------------------------------------------------------------------------------\n")
    cat("Direct (diagonal) and indirect (off-diagonal) effects\n")
    cat("----------------------------------------------------------------------------------------------\n")
    print(x$Coefficients)
    cat("----------------------------------------------------------------------------------------------\n")
  }
  if (export == TRUE) {
    sink()
  }
}
