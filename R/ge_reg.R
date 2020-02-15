#' Eberhart and Russell's regression model
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
#' single procedure use, for example, \code{resp = c(var1, var2, var3)}.
#' @param verbose Logical argument. If \code{verbose = FALSE} the code will run silently.
#' @return An object of class \code{ge_reg} with the folloing items for each variable:
#' \item{data}{The data with means for genotype and environment combinations and the
#' environment index}
#' \item{anova}{The analysis of variance for the regression model.}
#' \item{regression}{The estimated coefficients of the regression model.}
#' @seealso \code{\link{superiority}, \link{ecovalence}, \link{ge_stats}}
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
#' Crop Sci. 6:36-40. \href{https://www.crops.org/publications/cs/abstracts/6/1/CS0060010036}{doi:10.2135/cropsci1966.0011183X000600010011x}.

ge_reg = function(.data,
                  env,
                  gen,
                  rep,
                  resp,
                  verbose = TRUE){
  factors  <-
    .data %>%
    select({{env}}, {{gen}}, {{rep}}) %>%
    mutate_all(as.factor)
  vars <- .data %>% select({{resp}}, -names(factors))
  has_text_in_num(vars)
  vars %<>% select_numeric_cols()
  factors %<>% set_names("ENV", "GEN", "REP")
  listres <- list()
  nvar <- ncol(vars)
  for (var in 1:nvar) {
    data <- factors %>%
      mutate(mean = vars[[var]])
    data2 =  data  %>%
      dplyr::group_by(ENV, GEN) %>%
      dplyr::summarise(mean = mean(mean)) %>%
      as.data.frame()
    model1 <- lm(mean ~ GEN + ENV + ENV/REP + ENV * GEN, data = data)
    modav <- anova(model1)
    mydf = data.frame(aggregate(mean ~ GEN + ENV, data = data, mean))
    myAgg = aggregate(mean ~ GEN, mydf, "c")
    iamb = data.frame(aggregate(mean ~ ENV, data = data, mean))
    iamb = dplyr::mutate(iamb, IndAmb = mean - mean(mean))
    iamb2 = data.frame(aggregate(mean ~ ENV + GEN, data = data, mean))
    iamb2 = suppressMessages(dplyr::mutate(iamb2,
                                           IndAmb = dplyr::left_join(iamb2, iamb %>% select(ENV, IndAmb))$IndAmb))
    matx <- myAgg$mean
    meandf = data.frame(GEN = myAgg$GEN, myAgg$mean)
    names(meandf) = c("GEN", levels(mydf$ENV))
    iij = apply(matx, 2, mean) - mean(matx)
    YiIj = matx %*% iij
    bij = YiIj/sum((iij)^2)
    svar = (apply(matx^2, 1, sum)) - (((apply(matx, 1, sum))^2)/ncol(matx))
    bYijIj = bij * YiIj
    dij = svar - bYijIj
    pred = apply(matx, 1, mean) + bij %*% iij
    gof = function(x, y){
      R2 = NULL
      RMSE = NULL
      for (i in 1:nrow(x)){
        R2[i] =  cor(x[i, ], y[i, ])^2
        RMSE[i] = sqrt(sum((x[i, ] - y[i, ])^2)/ncol(x))
      }
      return(list(R2 = R2, RMSE = RMSE))
    }
    S2e = modav$"Mean Sq"[5]
    rps = length(levels(data$REP))
    en = length(levels(data$ENV))
    S2di = (dij/(en - 2)) - (S2e/rps)
    data2 = data2
    model2 <- lm(mean ~ GEN + ENV, data = data2)
    amod2 <- anova(model2)
    SSL = amod2$"Sum Sq"[2]
    SSGxL = amod2$"Sum Sq"[3]
    SS.L.GxL = SSL + SSGxL
    SSL.Linear = (1/length(levels(data$GEN))) * (colSums(matx) %*% iij)^2/sum(iij^2)
    SS.L.GxL.linear = sum(bYijIj) - SSL.Linear
    ge = length(levels(mydf$GEN))
    Df <- c(en * ge - 1, ge - 1, ge * (en - 1), 1, ge - 1, ge * (en - 2),
            replicate(length(dij), en - 2), en * ge * (rps - 1))
    poolerr = modav$"Sum Sq"[5]/rps
    SSS <- c(sum(amod2$"Sum Sq"), amod2$"Sum Sq"[1], SSL + SSGxL,
             SSL.Linear, SS.L.GxL.linear, sum(dij), dij, poolerr) * rps
    MSSS = (SSS/Df)
    FVAL = c(NA, MSSS[2]/MSSS[6], NA, NA, MSSS[5]/MSSS[6], NA,
             MSSS[7:(length(MSSS) - 1)]/MSSS[length(MSSS)], NA)
    PLINES = 1 - pf(FVAL[7:(length(MSSS) - 1)], Df[7], Df[length(Df)])
    pval = c(NA, 1 - pf(FVAL[2], Df[2], Df[6]), NA, NA, 1 -
               pf(FVAL[5], Df[5], Df[6]), NA, PLINES, NA)
    anovadf <- data.frame(Df, `Sum Sq` = SSS, `Mean Sq` = MSSS,
                          `F value` = FVAL, `Pr(>F)` = pval, check.names = FALSE)
    rownames(anovadf) <- c("Total", "GEN", "ENV + (GEN x ENV)", "ENV (linear)",
                           " GEN x ENV (linear)", "Pooled deviation",
                           levels(data$GEN), "Pooled error")
    temp = structure(list(data = iamb2,
                          anova = as_tibble(rownames_to_column(anovadf, "SV")),
                          regression = tibble(GEN = levels(mydf$GEN),
                                              Y = apply(matx, 1, mean),
                                              slope = as.numeric(bij),
                                              deviations = as.numeric(S2di),
                                              RMSE = gof(pred, matx)$RMSE,
                                              R2 = gof(pred, matx)$R2)),
                     class = "ge_reg")
    if (nvar > 1) {
      listres[[paste(names(vars[var]))]] <- temp
      if (verbose == TRUE) {
        cat("Evaluating variable", paste(names(vars[var])),
            round((var - 1)/(length(vars) - 1) * 100, 1), "%", "\n")
      }
    } else {
      listres[[paste(names(vars[var]))]] <- temp
    }
  }
  return(structure(listres, class = "ge_reg"))
}








#' Plot an object of class ge_reg
#'
#' Plot the regression model generated by the function \code{ge_reg}.
#'
#'
#' @param x An object of class \code{ge_factanal}
#' @param var The variable to plot. Defaults to \code{var = 1} the first
#'   variable of \code{x}.
#' @param type The type of plot to show. \code{type = 1} produces a plot with
#'   the environmental index in the x axis and the genotype mean yield in the y
#'   axis. \code{type = 2} produces a plot with the response variable in the x
#'   axis and the slope of the regression in the y axis.
#' @param plot_theme The graphical theme of the plot. Default is
#'   \code{plot_theme = theme_metan()}. For more details, see
#'   \code{\link[ggplot2]{theme}}.
#' @param x.lim The range of x-axis. Default is \code{NULL} (maximum and minimum
#'   values of the data set). New arguments can be inserted as \code{x.lim =
#'   c(x.min, x.max)}.
#' @param x.breaks The breaks to be plotted in the x-axis. Default is
#'   \code{authomatic breaks}. New arguments can be inserted as \code{x.breaks =
#'   c(breaks)}
#' @param x.lab The label of x-axis. Each plot has a default value. New
#'   arguments can be inserted as \code{x.lab = "my label"}.
#' @param y.lim The range of x-axis. Default is \code{NULL}. The same arguments
#'   than \code{x.lim} can be used.
#' @param y.breaks The breaks to be plotted in the x-axis. Default is
#'   \code{authomatic breaks}. The same arguments than \code{x.breaks} can be
#'   used.
#' @param y.lab The label of y-axis. Each plot has a default value. New
#'   arguments can be inserted as \code{y.lab = "my label"}.
#' @param leg.position The position of the legend.
#' @param size.tex.lab The size of the text in the axes text and labels. Default
#'   is \code{12}.
#' @param ... Currently not used..
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso \code{\link{ge_factanal}}
#' @method plot ge_reg
#' @return An object of class \code{gg, ggplot}.
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
      ggplot(x$data, aes(x = IndAmb, y = mean))+
      geom_point(aes(colour = GEN), size = 1.5)+
      geom_smooth(aes(colour = GEN), method = "lm", se = FALSE)+
      theme_bw()+
      labs(x = x.lab, y = y.lab)+
      plot_theme %+replace%
      theme(axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            axis.ticks = element_line(color = "black"),
            axis.ticks.length = unit(.2, "cm"),
            legend.position = leg.position)
  }
  if(type == 2){
    y.lab <- ifelse(missing(y.lab), "Slope of the regression", y.lab)
    x.lab <- ifelse(missing(x.lab), "Response variable", x.lab)
    p <-
      ggplot(x$regression, aes(x = Y, y = slope))+
      geom_point(size = 1.5)+
      geom_hline(yintercept = mean(x$regression$slope))+
      theme_bw()+
      geom_text_repel(aes(label = GEN))+
      labs(x = x.lab, y = y.lab)+
      plot_theme %+replace%
      theme(axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            axis.ticks = element_line(color = "black"),
            axis.ticks.length = unit(.2, "cm"),
            legend.position = leg.position)
  }
  return(p)
}







#' Print an object of class ge_reg
#'
#' Print the \code{ge_reg} object in two ways. By default, the results are shown
#' in the R console. The results can also be exported to the directory into a
#' *.txt file.
#'
#' @param x An object of class \code{ge_reg}.
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
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
    cat("\n\n\n")
  }
  if (export == TRUE) {
    sink()
  }
}
