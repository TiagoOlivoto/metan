#' @title Path coefficients with minimal multicollinearity
#' @name path_coeff
#' @description
#' `r badge('stable')`
#'
#' * `path_coeff()` computes a path analysis using a data frame as input data.
#' * `path_coeff_seq()` computes a sequential path analysis using primary and secondary traits.
#' * `path_coeff_mat()` computes a path analysis using correlation matrices as
#' input data.
#'
#' @details In `path_coeff()`, when `brutstep = TRUE`, an algorithm to
#'   select a set of predictors with minimal multicollinearity and high
#'   explanatory power is implemented. first, the algorithm will select a set of
#'   predictors with minimal multicollinearity. The selection is based on the
#'   variance inflation factor (VIF). An iterative process is performed until
#'   the maximum VIF observed is less than `maxvif`. The variables selected
#'   in this iterative process are then used in a series of stepwise-based
#'   regressions. The first model is fitted and p-1 predictor variables are
#'   retained (p is the number of variables selected in the iterative process.
#'   The second model adjusts a regression considering p-2 selected variables,
#'   and so on until the last model, which considers only two variables. Three
#'   objects are created. `Summary`, with the process summary,
#'   `Models`, containing the aforementioned values for all the adjusted
#'   models; and `Selectedpred`, a vector with the name of the selected
#'   variables in the iterative process.
#'
#' @param .data The data. Must be a data frame or a grouped data passed from
#'   [dplyr::group_by()]
#' @param resp <[`tidy-select`][dplyr_tidy_select]> The dependent trait.
#' @param cor_mat Matrix of correlations containing both dependent and
#'   independent traits.
#' @param pred <[`tidy-select`][dplyr_tidy_select]> The predictor traits. set to
#'   `everything()`, i.e., the predictor traits are all the numeric traits in
#'   the data except that in `resp`. To select multiple traits, use a
#'   comma-separated vector of names, (e.g., `pred = c(V1, V2, V2)`), an
#'   interval of trait names, (e.g., `pred = c(V1:V3)`), or even a select helper
#'   (e.g., `pred = starts_with("V")`).
#' @param chain_1,chain_2 <[`tidy-select`][dplyr_tidy_select]> The traits used
#'   in the first (primary) and second (secondary) chain.
#' @param by One variable (factor) to compute the function by. It is a shortcut
#'   to [dplyr::group_by()]. To compute the statistics by more than
#'   one grouping variable use that function.
#' @param exclude Logical argument, set to false. If `exclude = TRUE`, then
#'   the traits in `pred` are deleted from the data, and the analysis
#'   will use as predictor those that remained, except that in `resp`.
#' @param correction Set to `NULL`. A correction value (k) that will be
#'   added into the diagonal elements of the **X'X** matrix aiming at
#'   reducing the harmful problems of the multicollinearity in path analysis
#'   (Olivoto et al., 2017)
#' @param knumber When `correction = NULL`, a plot showing the values of
#'   direct effects in a set of different k values (0-1) is produced.
#'   `knumber` is the number of k values used in the range of 0 to 1.
#' @param brutstep Logical argument, set to `FALSE`. If true, then an
#'   algorithm will select a subset of variables with minimal multicollinearity
#'   and fit a set of possible models. See the **Details** section for more
#'   information.
#' @param maxvif The maximum value for the Variance Inflation Factor (cut point)
#'   that will be accepted. See the **Details** section for more information.
#' @param missingval How to deal with missing values. For more information,
#'   please see [stats::cor()].
#' @param plot_res If `TRUE`, create a scatter plot of residual against
#'   predicted value and a normal Q-Q plot.
#' @param verbose If `verbose = TRUE` then some results are shown in the
#'   console.
#' @param ... Depends on the function used:
#' * For `path_coeff()` additional arguments passed on to [stats::plot.lm()].
#' * For `path_coeff_seq()` additional arguments passed on to [path_coeff].
#' @return Depends on the function used:
#' * `path_coeff()`, returns a list with the following items:
#'    - **Corr.x** A correlation matrix between the predictor variables.
#'    - **Corr.y** A vector of correlations between each predictor variable
#' with the dependent variable.
#'    - **Coefficients** The path coefficients. Direct effects are the
#' diagonal elements, and the indirect effects those in the off-diagonal
#' elements (lines).
#'    - **Eigen** Eigenvectors and eigenvalues of the `Corr.x.`
#'    - **VIF** The Variance Inflation Factors.
#'    - **plot** A ggplot2-based graphic showing the direct effects in 21
#' different k values.
#'    - **Predictors** The predictor variables used in the model.
#'    - **CN** The Condition Number, i.e., the ratio between the highest and
#' lowest eigenvalue.
#'    - **Det** The matrix determinant of the `Corr.x.`.
#'    - **R2** The coefficient of determination of the model.
#'    - **Residual** The residual effect of the model.
#'    - **Response** The response variable.
#'    - **weightvar** The order of the predictor variables with the highest
#' weight (highest eigenvector) in the lowest eigenvalue.
#' * `path_coeff_seq()` returns a list with the following objects
#'    - **resp_fc** an object of class `path_coeff` with the results for the
#'    analysis with dependent trait and first chain predictors.
#'    - **resp_sc** an object of class `path_coeff` with the results for the
#'    analysis with dependent trait and second chain predictors.
#'    - **resp_sc2** The path coefficients of second chain predictors and the
#'    dependent trait through the first chain predictors
#'    - **fc_sc_list** A list of objects with the path analysis using each trait
#'    in the first chain as dependent and second chain as predictors.
#'    - **fc_sc_coef** The coefficients between first- and second-chain traits.
#'    - **cor_mat** A correlation matrix between the analyzed traits.
#' If `.data` is a grouped data passed from [dplyr::group_by()]
#' then the results will be returned into a list-column of data frames.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references
#' Olivoto, T., V.Q. Souza, M. Nardino, I.R. Carvalho, M. Ferrari, A.J.
#' Pelegrin, V.J. Szareski, and D. Schmidt. 2017. Multicollinearity in path
#' analysis: a simple method to reduce its effects. Agron. J. 109:131-142.
#' \doi{10.2134/agronj2016.04.0196}
#'
#' Olivoto, T., M. Nardino, I.R. Carvalho, D.N. Follmann, M. Ferrari, et al.
#' 2017. REML/BLUP and sequential path analysis in estimating genotypic values
#' and interrelationships among simple maize grain yield-related traits. Genet.
#' Mol. Res. 16(1): gmr16019525. \doi{10.4238/gmr16019525}
#'
#' @export
#' @examples
#' library(metan)
#'
#' # Using KW as the response variable and all other ones as predictors
#' pcoeff <- path_coeff(data_ge2, resp = KW)
#'
#' # The same as above, but using the correlation matrix
#' cor_mat <- cor(data_ge2 %>% select_numeric_cols())
#' pcoeff2 <- path_coeff_mat(cor_mat, resp = KW)
#'
#' # Declaring the predictors
#' # Create a residual plot with 'plot_res = TRUE'
#' pcoeff3<- path_coeff(data_ge2,
#'                       resp = KW,
#'                       pred = c(PH, EH, NKE, TKW),
#'                       plot_res = TRUE)
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
#' pcoeff5 <- path_coeff(data_ge2, resp = KW, by = ENV)
#'
#'
#'# sequential path analysis
#'# KW as dependent trait
#'# NKE and TKW as primary predictors
#'# PH, EH, EP, and EL as secondary traits
#' pcoeff6 <-
#'  path_coeff_seq(data_ge2,
#'                resp = KW,
#'                chain_1 = c(NKE, TKW),
#'                chain_2 = c(PH, EH, EP, EL))
#' pcoeff6$resp_sc$Coefficients
#' pcoeff6$resp_sc2
#'
#'
path_coeff <- function(.data,
                       resp,
                       pred = everything(),
                       by = NULL,
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
        betas[i, 2:nvar] <- t(solve_svd(cor.x2) %*% cor.y)
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
      p1 <-
        ggplot(betas, aes(K, direct, col = VAR)) +
        geom_line(size = 1) +
        theme_bw() +
        theme(axis.ticks.length = unit(0.2, "cm"),
              axis.text = element_text(size = 12, colour = "black"),
              axis.title = element_text(size = 12, colour = "black"),
              axis.ticks = element_line(colour = "black"),
              legend.position = "bottom", plot.margin = margin(0.1,  0.1, 0.1, 0.1, "cm"),
              legend.title = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA, size = 1),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank()) +
        labs(x = "k values", y = expression(beta~values)) +
        geom_abline(intercept = 0,  slope = 0) +
        scale_x_continuous(breaks = seq(0, 1, by = 0.1))
    } else {
      p1 <- "No graphic generated due to correction value"
    }
    if(det(cor.x) < 1e-15){
      warning("System is computationally singular. Check the collinearity of predictors with `metan::colindiag()`.", call. = FALSE)
    }
    eigen <- eigen(cor.x)
    Det <- det(cor.x)
    NC <- max(eigen$values)/min(eigen$values)
    Aval <- data.frame(eigen$values)
    names(Aval) <- "Eigenvalues"
    Avet <- data.frame(t(eigen$vectors))
    names(Avet) <- names
    AvAvet <- cbind(Aval, Avet)
    Direct <- solve_svd(cor.x %*% t(cor.x)) %*% cor.x %*% cor.y
    n <- ncol(cor.x)
    Coeff <- data.frame(cor.x)
    for (i in 1:n) {
      for (j in 1:n) {
        Coeff[i, j] <- Direct[j] * cor.x[i, j]
      }
    }
    Residual <- sqrt(1 - t(Direct) %*% cor.y)
    R2 <- t(Direct) %*% cor.y
    VIF <- data.frame(diag(solve_svd(cor.x)))
    names(VIF) <- "VIF"
    if (verbose == TRUE) {
      if (NC > 1000 | NC < 0) {
        cat(paste0("Severe multicollinearity. \n",
                   "Condition Number: ", round(NC, 3), "\n",
                   "Consider using a correction factor with 'correction' argument.\n",
                   "Consider identifying collinear traits with `non_collinear_vars()`\n"))
      }
      if (NC >= 0 & NC <= 100) {
        cat(paste0("Weak multicollinearity. \n", "Condition Number: ",
                   round(NC, 3), "\n", "You will probably have path coefficients close to being unbiased. \n"))
      }
      if (NC > 100 & NC < 1000) {
        cat(paste0("Moderate multicollinearity! \n",
                   "Condition Number: ", round(NC, 3), "\n",
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
                           Coefficients = cbind(as_tibble(Coeff, rownames = NA), as_tibble(cor.y, rownames = NA) %>% set_names("linear")),
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
      data_sel <- cbind(yyy, xxx)
      model_sel <- select_pred(data_sel,
                               resp = names(yyy),
                               npred = npred)
      predstw <- model_sel$predictors
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
          betas[i, 2:nvar] <- t(solve_svd(cor.x2) %*% cor.y)
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
        p1 <-
          ggplot(betas, aes(K, direct, col = VAR)) +
          geom_line(size = 1) +
          theme_bw() +
          theme(axis.ticks.length = unit(0.2, "cm"),
                axis.text = element_text(size = 12, colour = "black"),
                axis.title = element_text(size = 12, colour = "black"),
                axis.ticks = element_line(colour = "black"),
                legend.position = "bottom", plot.margin = margin(0.1,  0.1, 0.1, 0.1, "cm"),
                legend.title = element_blank(),
                panel.border = element_rect(colour = "black", fill = NA, size = 1),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.minor.y = element_blank()) +
          labs(x = "k values", y = expression(beta~values)) +
          geom_abline(intercept = 0,  slope = 0) +
          scale_x_continuous(breaks = seq(0, 1, by = 0.1))
      } else {
        p1 <- "No graphic generated due to correction value"
      }
      if(det(cor.x) < 1e-15){
        warning("System is computationally singular. Check the collinearity of predictors with `metan::colindiag()`.", call. = FALSE)
      }
      eigen <- eigen(cor.x)
      Det <- det(cor.x)
      NC <- max(eigen$values)/min(eigen$values)
      Aval <- as.data.frame(eigen$values)
      names(Aval) <- "Eigenvalues"
      Avet <- as.data.frame(eigen$vectors)
      names(Avet) <- names
      AvAvet <- cbind(Aval, Avet)
      Direct <- solve_svd(cor.x %*% t(cor.x)) %*% cor.x %*% cor.y
      n <- ncol(cor.x)
      Coeff <- data.frame(cor.x)
      for (i in 1:n) {
        for (j in 1:n) {
          Coeff[i, j] <- Direct[j] * cor.x[i, j]
        }
      }
      Residual <- sqrt(1 - t(Direct) %*% cor.y)
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
                      Coefficients = cbind(data.frame(Coeff), data.frame(cor.y) %>% set_names("linear")),
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
      statistics[i, 2] <- model_sel$AIC
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



#' @name path_coeff
#' @export
path_coeff_mat <- function(cor_mat,
                           resp,
                           correction = NULL,
                           knumber = 50,
                           verbose = TRUE){
  cor_mat <- as.matrix(cor_mat)
  if(!isSymmetric(cor_mat)){
    stop("Object '", test["cor_mat"], "' must be a symmetric.", call. = FALSE)
  }
  cor.y <-
    cor_mat %>%
    as.data.frame() %>%
    select({{resp}}) %>%
    subset(rownames(.) != colnames(.)) %>%
    as.matrix()
  cor.x <-
    cor_mat %>%
    as.data.frame() %>%
    remove_cols({{resp}}) %>%
    subset(rownames(.) != colnames(cor.y)) %>%
    as.matrix()
  kincrement <- 1/knumber
  if (is.null(correction) == FALSE) {
    diag(cor.x) <- diag(cor.x) + correction
  } else {
    cor.x <- cor.x
  }
  names <- colnames(cor.x)
  if (is.null(correction) == TRUE) {
    betas <- data.frame(matrix(nrow = knumber, ncol = ncol(cor.x) + 1))
    cc <- 0
    nvar <- ncol(cor.x) + 1
    for (i in 1:knumber) {
      cor.x2 <- cor.x
      diag(cor.x2) <- diag(cor.x2) + cc
      betas[i, 1] <- cc
      betas[i, 2:nvar] <- t(solve_svd(cor.x2) %*% cor.y)
      cc <- cc + kincrement
    }
    names(betas) <- paste0(c("K", colnames(cor.x)))
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
    p1 <-
      ggplot(betas, aes(K, direct, col = VAR)) +
      geom_line(size = 1) +
      theme_bw() +
      theme(axis.ticks.length = unit(0.2, "cm"),
            axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.ticks = element_line(colour = "black"),
            legend.position = "bottom", plot.margin = margin(0.1,  0.1, 0.1, 0.1, "cm"),
            legend.title = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, size = 1),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank()) +
      labs(x = "k values", y = "Beta values") + geom_abline(intercept = 0,  slope = 0) +
      scale_x_continuous(breaks = seq(0, 1, by = 0.1))
  } else {
    p1 <- "No graphic generated due to correction value"
  }
  if(det(cor.x) < 1e-15){
    warning("System is computationally singular. Check the collinearity of predictors with `metan::colindiag()`.", call. = FALSE)
  }
  eigen <- eigen(cor.x)
  Det <- det(cor.x)
  NC <- max(eigen$values)/min(eigen$values)
  Aval <- data.frame(eigen$values)
  names(Aval) <- "Eigenvalues"
  Avet <- data.frame(t(eigen$vectors))
  names(Avet) <- names
  AvAvet <- cbind(Aval, Avet)
  Direct <- solve_svd(cor.x %*% t(cor.x)) %*% cor.x %*% cor.y
  n <- ncol(cor.x)
  Coeff <- data.frame(cor.x)
  for (i in 1:n) {
    for (j in 1:n) {
      Coeff[i, j] <- Direct[j] * cor.x[i, j]
    }
  }
  Residual <- sqrt(1 - t(Direct) %*% cor.y)
  R2 <- t(Direct) %*% cor.y
  VIF <- data.frame(diag(solve_svd(cor.x)))
  names(VIF) <- "VIF"
  if (verbose == TRUE) {
    if (NC > 1000 | NC <= 0) {
      cat(paste0("Severe multicollinearity. \n",
                 "Condition Number = ", round(NC, 3), "\n",
                 "Please, consider using a correction factor, or use 'brutstep = TRUE'. \n"))
    }
    if (NC >= 0 & NC <= 100) {
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
                         Coefficients = cbind(as_tibble(Coeff, rownames = NA), as_tibble(cor.y, rownames = NA) %>% set_names("linear")),
                         Eigen = as_tibble(AvAvet),
                         VIF =  rownames_to_column(VIF, "VAR") %>% as_tibble(),
                         plot = p1,
                         Response = colnames(cor.y),
                         Predictors = colnames(cor.x),
                         CN = NC,
                         Det = Det,
                         R2 = R2,
                         Residual = Residual,
                         weightvar = weightvarname),
                    class = "path_coeff")
  return(temp)
}


#' @name path_coeff
#' @export
path_coeff_seq <- function(.data,
                           resp,
                           chain_1,
                           chain_2,
                           by = NULL,
                           verbose = TRUE,
                           ...){
  if (!missing(by)) {
    if (length(as.list(substitute(by))[-1L]) != 0) {
      stop("Only one grouping variable can be used in the argument 'by'.\nUse 'group_by()' to pass '.data' grouped by more than one variable.", call. = FALSE)
    }
    .data <- group_by(.data, {{by}})
  }
  if (is_grouped_df(.data)){
    results <- .data %>%
      doo(path_coeff_seq,
          resp = {{resp}},
          chain_1 = {{chain_1}},
          chain_2 = {{chain_2}},
          verbose = verbose,
          ...)
    return(set_class(results, c("group_path_seq", "tbl_df", "tbl",  "data.frame")))
  }
  make_names <- function(vars_c2){
    mat_names <- matrix(rep(paste0("indirect_", vars_c2), length(vars_c2)),
                        ncol = length(vars_c2),
                        nrow = length(vars_c2))
    diag(mat_names) <- "direct"
    mat_names <- rbind(mat_names, paste0("linear"))
    return(as.vector(mat_names))
  }
  cor_coef <- corr_coef(.data, {{resp}}, {{chain_1}}, {{chain_2}})$cor
  # main and first chain
  if(isTRUE(verbose)){
    cat("========================================================\n")
    cat("Collinearity diagnosis of first chain predictors\n")
    cat("========================================================\n")
  }
  resp_fc <- path_coeff(.data,
                        resp = {{resp}},
                        pred = {{chain_1}},
                        verbose = verbose,
                        ...)
  # main and second chain
  if(isTRUE(verbose)){
    cat("========================================================\n")
    cat("Collinearity diagnosis of second chain predictors\n")
    cat("========================================================\n")
  }
  resp_sc <- path_coeff(.data,
                        resp = {{resp}},
                        pred = {{chain_2}},
                        verbose = verbose,
                        ...)
  ncoef <- length(as.matrix(resp_sc$Coefficients))
  vars_c1 <- resp_fc$Predictors
  vars_c2 <- resp_sc$Predictors
  nfc <- length(vars_c1)
  fc_sc <- list()
  res <- as.data.frame(matrix(nrow = ncoef, ncol = nfc))
  names_effects <- make_names(vars_c2)
  for (i in 1:nfc) {
    tmp <- path_coeff(.data,
                      resp =  vars_c1[[i]],
                      pred = {{chain_2}},
                      verbose = FALSE,
                      ...)
    fc_sc[[vars_c1[[i]]]] <- tmp
    coefs <- as.matrix(tmp$Coefficients)
    n <- 0
    for (j in 1:nrow(coefs)) {
      for (k in 1:ncol(coefs)) {
        n <- n + 1
        res[n, i] <- coefs[j, k]
      }
    }
  }
  stats <- rbind(
    sapply(fc_sc, function(x){
      x$R2
    }),
    sapply(fc_sc, function(x){
      x$Residual
    })
  )
  names(stats) <- vars_c1
  names(res) <- vars_c1
  res <- rbind(res, stats)
  res$trait <- c(rep(vars_c2, each = length(vars_c2) + 1), "R2", "Residual")
  res$effects <- c(names_effects, "R2", "Residual")
  res <- res[c("trait", "effects", setdiff(colnames(res), c("effects", "trait")))]

  # coefficients of second chain
  coef_chain_2 <- resp_sc$Coefficients %>% remove_cols(linear) %>% as.matrix()
  # direct effects first chain and main trait
  fcmv <- resp_fc$Coefficients
  dir_fcmv <- diag(as.matrix(fcmv[, 1:(ncol(fcmv) - 1)]))
  # direct effects first chain and second chain
  fcsc <- res[1:(nrow(res) - 2), ]
  dir_fcsc <-
    fcsc[which(fcsc$effects == "direct"), -2] %>%
    remove_rownames() %>%
    column_to_rownames("trait")

  effects <- list()
  for (i in 1:length(vars_c2)) {
    sec_trait <- vars_c2[i]
    npt <- length(vars_c1)
    nst <- length(vars_c2)
    if (nst < 2) {
      var_dif <- setdiff(vars_c2, sec_trait)
      rnam <- c("direct", "total")
    } else{
      var_dif <- setdiff(vars_c2, sec_trait)
      rnam <- c("direct", paste0("indirect_",var_dif), "total")
    }
    tmp <- matrix(nrow = nst, ncol = npt)
    dir_ef <- dir_fcsc[i, ]
    for (j in 1:nst) {
      if (j == 1) {
        nv <- 0
      } else{
        nv <- nv + 1
      }
      for (k in 1:npt) {
        if (j == 1) {
          tmp[j, k] <- dir_fcmv[k] * dir_ef[[k]]
        } else{
          var_ind <- var_dif[nv]
          cor_i <- which(rownames(cor_coef) == sec_trait)
          cor_j <- which(colnames(cor_coef) == var_ind)
          tmp[j, k] <- dir_fcmv[k] * dir_fcsc[which(rownames(dir_fcsc) == var_ind), k] * cor_coef[cor_i, cor_j]
        }
      }
    }
    colnames(tmp) <- vars_c1
    tmp_coef <- coef_chain_2[i, ][c(i, setdiff(1:length(coef_chain_2[i, ]), i))]
    tmp <-
      as.data.frame(tmp) %>%
      mutate(total =  as.vector(tmp_coef),
             residual = rowSums(tmp) - total) %>%
      reorder_cols(residual, .before = total)
    tmp <-
      rbind(tmp, colSums(tmp)) %>%
      mutate(trait = sec_trait,
             effect = rnam,
             .before = 1)
    effects[[sec_trait]] <- tmp
  }
  return(list(resp_fc = resp_fc,
              resp_sc = resp_sc,
              resp_sc2 =  rbind_fill_id(effects),
              fc_sc_list = fc_sc,
              fc_sc_coef = res,
              cor_mat = cor_coef) %>%
           set_class("path_coeff_seq"))
}

#' Print an object of class path_coeff
#'
#' Print an object generated by the function 'path_coeff()'. By default, the
#' results are shown in the R console. The results can also be exported to the
#' directory.
#'
#'
#' @param x An object of class `path_coeff` or `group_path`.
#' @param export A logical argument. If `TRUE`, a *.txt file is exported to
#'   the working directory
#' @param file.name The name of the file if `export = TRUE`
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   [`tibble::print()`][tibble::formatting] for more details.
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
    print.data.frame(x$Corr.x, digits = digits)
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



#' Plots an object of class `path_coeff`
#'
#' Plots an object generated by `path_coeff()`. Options includes the path
#' coefficients, variance inflaction factor and the beta values with different
#' values of 'k' values added to the diagonal of the correlation matrix of
#' explanatory traits. See more on `Details` section.
#'
#' The plot `which = "coef"` (default) is interpreted as follows:
#' * The direct effects are shown in the diagonal (highlighted with a thicker
#' line). In the example, the direct effect of NKE on KW is 0.718.
#' * The indirect effects are shown in the line. In the example, the indirect effect of EH on KW through TKW is 0.396.
#' * The linear correlation (direct + indirect) is shown in the last column.
#'
#' @param x An object of class `path_coeff` or `group_path`.
#' @param which Which to plot: one of `'coef'`, `'vif'`, or `'betas'`.
#' @param size.text.plot,size.text.labels The size of the text for plot area and
#'   labels, respectively.
#' @param digits The significant digits to be shown.
#' @param ... Further arguments passed on to [ggplot2::theme()].
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method plot path_coeff
#' @export
#' @examples
#' \donttest{
#' library(metan)
#'
#' # KW as dependent trait and all others as predictors
#' # PH, EH, NKE, and TKW as predictors
#'
#' pcoeff <-
#'   path_coeff(data_ge2,
#'              resp = KW,
#'              pred = c(PH, EH, NKE, TKW))
#' plot(pcoeff)
#' plot(pcoeff, which = "vif")
#' plot(pcoeff, which = "betas")
#' }
plot.path_coeff <- function(x,
                            which = "coef",
                            size.text.plot = 4,
                            size.text.labels = 10,
                            digits = 3,
                            ...){
  if(!which %in% c("coef", "vif", "betas")){
    warning("'which' must be one of 'coef', 'vif', or 'betas'. Plotting the coefficients.", call. = FALSE)
    which <- "coef"
  }
  if(which == "coef"){
    mat <- x$Coefficients
    r2 <- format(round(x$R2, digits = digits), nsmall = digits)
    res <- format(round(x$Residual, digits = digits), nsmall = digits)

    helper_plot <- function(mat) {
      mat <- make_long(mat)
      fcts <- as.character(unique(factor(mat$GEN)))
      mat <-
        mat %>%
        mutate(lwid = ifelse((ENV == GEN), 1.5, 0.5)) %>%
        mutate(ENV =  factor(as.factor(ENV), levels = c(fcts, "linear"))) %>%
        mutate(GEN = factor(as.factor(GEN), levels = rev(fcts))) %>%
        arrange(lwid)
      p <-
        ggplot(mat, aes(x = ENV, y = GEN, fill = Y)) +
        geom_tile(colour = "black",
                  size = mat$lwid) +
        scale_x_discrete(position = "top") +
        scale_fill_gradient2(low = "red",
                             mid = "white",
                             high = "blue",
                             limits = c(-1, 1)) +
        geom_text(aes(label = format(round(Y, digits = digits), nsmall = digits)),
                  size = size.text.plot) +

        labs(x = NULL,
             y = NULL,
             caption = expr(paste("Response: ", !!x$Response, ~~ "|",
                                  ~~R^2, ": ", !!r2, ~~ "|",
                                  ~~Residual, ": ", !!res, sep = ""))) +
        theme_minimal() +
        theme(legend.position = 'bottom',
              legend.title = element_blank(),
              legend.key.width = unit(60, "pt"),
              legend.key.height = unit(10, "pt"),
              axis.text = element_text(color = "black", size = size.text.labels),
              ...) +
        coord_fixed()
      return(p)
    }
    return(helper_plot(mat))
  }
  if(which == "vif"){
    vifs <- x$VIF
    p <-
      ggplot(vifs, aes(x = VAR, y = VIF))+
      geom_col(color = "black", fill = "gray", width = 0.7) +
      geom_text(aes(label = round(VIF, 2)),
                angle = 90, hjust = -0.2,
                size = size.text.plot) +
      geom_hline(yintercept = 10, linetype = 2, col = 'red') +
      labs(x = "Traits",
           y = "Variance Inflation Factor",
           caption = "The red dashed line shows the value of VIF = 10.") +
      scale_y_continuous(expand = expansion(c(0, 0.15))) +
      theme(axis.text = element_text(color = "black",
                                     size = size.text.labels),
            panel.grid.minor = element_blank()) +
      theme_bw() +
      theme(axis.ticks.length = unit(0.2, "cm"),
            axis.text = element_text(size = 12, colour = "black"),
            axis.title = element_text(size = 12, colour = "black"),
            axis.ticks = element_line(colour = "black"),
            panel.border = element_rect(colour = "black", fill = NA, size = 1),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            ...)
    return(p)
  }
  if(which == "betas"){
    x$plot
  }
}
