#' Eberhart and Russell's regression model
#' @description
#' `r badge('stable')`
#'
#' Regression-based stability analysis using the Eberhart and Russell (1966) model.
#'
#' @param .data The dataset containing the columns related to Environments, Genotypes,
#'              replication/block and response variable(s)
#' @param env The name of the column that contains the levels of the
#' environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#' replications/blocks
#' @param resp The response variable(s). To analyze multiple variables in a
#' single procedure use, for example, `resp = c(var1, var2, var3)`.
#' @param verbose Logical argument. If `verbose = FALSE` the code will run
#'   silently.
#' @return An object of class `ge_reg` with the folloing items for each
#'   variable:
#' * data: The data with means for genotype and environment combinations and the
#' environment index
#' * anova: The analysis of variance for the regression model.
#' * regression: A data frame with the following columns: `GEN`, the genotypes;
#' `b0` and `b1` the intercept and slope of the regression, respectively;
#' `t(b1=1)` the calculated t-value; `pval_t` the p-value for the t test; `s2di`
#' the deviations from the regression (stability parameter); `F(s2di=0)` the
#' F-test for the deviations; `pval_f` the p-value for the F test; `RMSE` the
#' root-mean-square error; `R2` the determination coefficient of the regression.
#' * b0_variance: The variance of b0.
#' * b1_variance: The variance of b1.
#' @md
#' @seealso [superiority()], [ecovalence()], [ge_stats()]
#' @author Tiago Olivoto, \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#' \donttest{
#' library(metan)
#'reg <- ge_reg(data_ge2,
#'              env = ENV,
#'              gen = GEN,
#'              rep = REP,
#'              resp = PH)
#'plot(reg)
#'
#'}
#' @references Eberhart, S.A., and W.A. Russell. 1966. Stability parameters for comparing Varieties.
#' Crop Sci. 6:36-40. \doi{10.2135/cropsci1966.0011183X000600010011x}

ge_reg = function(.data,
                  env,
                  gen,
                  rep,
                  resp,
                  verbose = TRUE){
  factors  <-
    .data %>%
    select({{env}}, {{gen}}, {{rep}}) %>%
    mutate(across(everything(), as.factor))
  vars <-
    .data %>%
    select({{resp}}, -names(factors)) %>%
    select_numeric_cols()
  factors %<>% set_names("ENV", "GEN", "REP")
  listres <- list()
  nvar <- ncol(vars)
  if (verbose == TRUE) {
    pb <- progress(max = nvar, style = 4)
  }
  for (var in 1:nvar) {
    data <-
      factors %>%
      mutate(Y = vars[[var]])
    if(has_na(data)){
      data <- remove_rows_na(data)
      has_text_in_num(data)
    }
    data2 <-
      data %>%
      means_by(ENV, GEN) %>%
      as.data.frame()
    model1 <- lm(Y ~ GEN + ENV + ENV/REP + ENV * GEN, data = data)
    modav <- anova(model1)
    mydf <-
      data %>%
      means_by(GEN, ENV)
    iamb <-
      data %>%
      means_by(ENV) %>%
      add_cols(IndAmb = Y - mean(Y)) %>%
      select_cols(ENV, IndAmb)

    iamb2 <-
      data %>%
      means_by(ENV, GEN) %>%
      left_join(iamb, by = "ENV")
    meandf <- make_mat(mydf, GEN, ENV, Y) %>% rownames_to_column("GEN")
    matx <- make_mat(mydf, GEN, ENV, Y) %>% as.matrix()
    iij <- apply(matx, 2, mean) - mean(matx)
    sumij2 <- sum((iij)^2)
    YiIj <- matx %*% iij
    bij <- YiIj/sumij2
    svar <- (apply(matx^2, 1, sum)) - (((apply(matx, 1, sum))^2)/ncol(matx))
    bYijIj <- bij * YiIj
    dij <- svar - bYijIj
    pred <- apply(matx, 1, mean) + bij %*% iij
    gof <- function(x, y){
      R2 = NULL
      RMSE = NULL
      for (i in 1:nrow(x)){
        R2[i] =  cor(x[i, ], y[i, ])^2
        RMSE[i] = sqrt(sum((x[i, ] - y[i, ])^2)/ncol(x))
      }
      return(list(R2 = R2, RMSE = RMSE))
    }
    S2e <- modav$"Mean Sq"[5]
    nrep <- length(levels(data$REP))
    en <- length(levels(data$ENV))
    ge <- length(levels(mydf$GEN))
    S2di <- (dij/(en - 2)) - (S2e/nrep)
    amod2 <- anova(lm(Y ~ GEN + ENV, data = data2))
    # amod2 <- anova(model2)
    SSL <- amod2$"Sum Sq"[2]
    SSGxL <- amod2$"Sum Sq"[3]
    SS.L.GxL <- SSL + SSGxL
    SSL.Linear <- (1/length(levels(data$GEN))) * (colSums(matx) %*% iij)^2/sum(iij^2)
    SS.L.GxL.linear <- sum(bYijIj) - SSL.Linear
    Df <- c(en * ge - 1,
            ge - 1,
            ge * (en - 1),
            1,
            ge - 1,
            ge * (en - 2),
            replicate(length(dij), en - 2),
            en*(nrep - 1) * (ge - 1))
    poolerr <- modav$"Sum Sq"[5]/nrep
    sigma2 <- modav$"Mean Sq"[5]
    dferr <- modav$"Df"[5]
    vbo <- sigma2 / (en * nrep)
    vb1 <- sigma2 / (nrep * sumij2)
    tcal <- (bij - 1) / sqrt(vb1)
    ptcal <- 2 * pt(abs(tcal), dferr, lower.tail = FALSE)
    SSS <- c(sum(amod2$"Sum Sq"),
             amod2$"Sum Sq"[1],
             SSL + SSGxL,
             SSL.Linear,
             SS.L.GxL.linear,
             sum(dij),
             dij,
             poolerr) * nrep
    MSSS <- (SSS/Df)
    FVAL <- c(NA,
              MSSS[2]/MSSS[6],
              NA,
              NA,
              MSSS[5]/MSSS[6],
              NA,
              MSSS[7:(length(MSSS) - 1)]/MSSS[length(MSSS)],
              NA)
    PLINES <- 1 - pf(FVAL[7:(length(MSSS) - 1)], Df[7], Df[length(Df)])
    pval <- c(NA,
              1 - pf(FVAL[2], Df[2], Df[6]),
              NA,
              NA,
              1 - pf(FVAL[5], Df[5], Df[6]),
              NA,
              PLINES,
              NA)
    anovadf <- data.frame(Df,
                          `Sum Sq` = SSS,
                          `Mean Sq` = MSSS,
                          `F value` = FVAL,
                          `Pr(>F)` = pval,
                          check.names = FALSE)
    rownames(anovadf) <- c("Total", "GEN", "ENV + (GEN x ENV)", "ENV (linear)",
                           " GEN x ENV (linear)", "Pooled deviation",
                           levels(data$GEN), "Pooled error")
    temp <- structure(list(data = iamb2,
                           anova = as_tibble(rownames_to_column(anovadf, "SV")),
                           regression = tibble(GEN = levels(mydf$GEN),
                                               b0 = apply(matx, 1, mean),
                                               b1 = as.numeric(bij),
                                               `t(b1=1)` = tcal,
                                               pval_t = ptcal,
                                               s2di = as.numeric(S2di),
                                               `F(s2di=0)` = FVAL[7:(length(FVAL) - 1)],
                                               pval_f = PLINES,
                                               RMSE = gof(pred, matx)$RMSE,
                                               R2 = gof(pred, matx)$R2),
                           bo_variance = vbo,
                           b1_variance = vb1),
                      class = "ge_reg")
    if (verbose == TRUE) {
      run_progress(pb,
                   actual = var,
                   text = paste("Evaluating trait", names(vars[var])))
    }
    listres[[paste(names(vars[var]))]] <- temp
  }
  return(structure(listres, class = "ge_reg"))
}



#' Plot an object of class ge_reg
#'
#' Plot the regression model generated by the function `ge_reg`.
#'
#'
#' @param x An object of class `ge_factanal`
#' @param var The variable to plot. Defaults to `var = 1` the first
#'   variable of `x`.
#' @param type The type of plot to show. `type = 1` produces a plot with
#'   the environmental index in the x axis and the genotype mean yield in the y
#'   axis. `type = 2` produces a plot with the response variable in the x
#'   axis and the slope/deviations of the regression in the y axis.
#' @param plot_theme The graphical theme of the plot. Default is
#'   `plot_theme = theme_metan()`. For more details, see
#'   [ggplot2::theme()].
#' @param x.lim The range of x-axis. Default is `NULL` (maximum and minimum
#'   values of the data set). New arguments can be inserted as `x.lim =
#'   c(x.min, x.max)`.
#' @param x.breaks The breaks to be plotted in the x-axis. Default is
#'   `authomatic breaks`. New arguments can be inserted as `x.breaks =
#'   c(breaks)`
#' @param x.lab The label of x-axis. Each plot has a default value. New
#'   arguments can be inserted as `x.lab = "my label"`.
#' @param y.lim The range of x-axis. Default is `NULL`. The same arguments
#'   than `x.lim` can be used.
#' @param y.breaks The breaks to be plotted in the x-axis. Default is
#'   `authomatic breaks`. The same arguments than `x.breaks` can be
#'   used.
#' @param y.lab The label of y-axis. Each plot has a default value. New
#'   arguments can be inserted as `y.lab = "my label"`.
#' @param leg.position The position of the legend.
#' @param size.tex.lab The size of the text in the axes text and labels. Default
#'   is `12`.
#' @param ... Currently not used..
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso [ge_factanal()]
#' @method plot ge_reg
#' @return An object of class `gg, ggplot`.
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' model <- ge_reg(data_ge2, ENV, GEN, REP, PH)
#' plot(model)
#' }
#'
plot.ge_reg <- function(x,
                        var = 1,
                        type = 1,
                        plot_theme = theme_metan(),
                        x.lim = NULL,
                        x.breaks = waiver(),
                        x.lab = NULL,
                        y.lim = NULL,
                        y.breaks = waiver(),
                        y.lab = NULL,
                        leg.position = "right",
                        size.tex.lab = 12,
                        ...){
  x <- x[[var]]
  if (!class(x) == "ge_reg") {
    stop("The object 'x' is not of class 'ge_reg'", call. = FALSE)
  }
  if(!type  %in% c(1, 2)){
    stop("Argument 'type' must be either 1 or 2", call. = FALSE)
  }
  if(type == 1){
    y.lab <- ifelse(missing(y.lab), "Response variable", y.lab)
    x.lab <- ifelse(missing(x.lab), "Environmental index", x.lab)
    p <-
      ggplot(x$data, aes(x = IndAmb, y = Y))+
      geom_point(aes(colour = GEN), size = 1.5)+
      geom_smooth(aes(colour = GEN), method = "lm", formula = y ~ x, se = FALSE)+
      theme_bw()+
      labs(x = x.lab, y = y.lab)+
      plot_theme %+replace%
      theme(axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            axis.ticks = element_line(color = "black"),
            axis.ticks.length = unit(.2, "cm"),
            legend.position = leg.position)
    return(p)
  }
  if(type == 2){
    y.lab <- ifelse(missing(y.lab), "Slope of the regression", y.lab)
    x.lab <- ifelse(missing(x.lab), "Response variable", x.lab)
    p <-
      ggplot(x$regression, aes(x = b0, y = b1))+
      geom_point(size = 1.5)+
      geom_hline(yintercept = mean(x$regression$b1))+
      geom_text_repel(aes(label = GEN))+
      labs(x = x.lab, y = y.lab) +
      plot_theme

    p2 <-
      ggplot(x$regression, aes(x = b0, y = s2di))+
      geom_point(size = 1.5)+
      geom_hline(yintercept = mean(x$regression$s2di))+
      geom_text_repel(aes(label = GEN))+
      labs(x = x.lab, y = "Deviations from the regression") +
      plot_theme

    p + p2

  }
}

#' Print an object of class ge_reg
#'
#' Print the `ge_reg` object in two ways. By default, the results are shown
#' in the R console. The results can also be exported to the directory into a
#' *.txt file.
#'
#' @param x An object of class `ge_reg`.
#' @param export A logical argument. If `TRUE`, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if `export = TRUE`
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   [`tibble::print()`][tibble::formatting] for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print ge_reg
#' @export
#' @examples
#' \donttest{
#'
#' library(metan)
#' model <- ge_reg(data_ge2, ENV, GEN, REP, PH)
#' print(model)
#' }
print.ge_reg <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (!class(x) == "ge_reg") {
    stop("The object must be of class 'ge_reg'")
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "ge_reg print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    cat("Joint-regression Analysis of variance\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$anova)
    cat("---------------------------------------------------------------------------\n")
    cat("Regression parameters\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$regression)
    cat("---------------------------------------------------------------------------\n")
    cat("Variance of b0:", var[["bo_variance"]], "\n")
    cat("Variance of b1:", var[["b1_variance"]], "\n")
    cat("\n\n\n")
  }
  if (export == TRUE) {
    sink()
  }
}
