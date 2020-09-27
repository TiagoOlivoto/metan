#' Genotype analysis by mixed-effect models
#'
#' Analysis of genotypes in single experiments using mixed-effect models with
#' estimation of genetic parameters.
#'
#'
#' @param .data The dataset containing the columns related to, Genotypes,
#' replication/block and response variable(s).
#' @param gen The name of the column that contains the levels of the genotypes,
#'   that will be treated as random effect.
#' @param rep The name of the column that contains the levels of the
#'   replications (assumed to be fixed).
#' @param resp The response variable(s). To analyze multiple variables in a
#' single procedure a vector of variables may be used. For example \code{resp =
#' c(var1, var2, var3)}. Select helpers are also allowed.
#' @param block Defaults to \code{NULL}. In this case, a randomized complete
#'   block design is considered. If block is informed, then an alpha-lattice
#'   design is employed considering block as random to make use of inter-block
#'   information, whereas the complete replicate effect is always taken as
#'   fixed, as no inter-replicate information was to be recovered (Mohring et
#'   al., 2015).
#'@param by One variable (factor) to compute the function by. It is a shortcut
#'  to \code{\link[dplyr]{group_by}()}.This is especially useful, for example,
#'  when the researcher what to fit a mixed-effect model for each environment.
#'  In this case, an object of class gamem_grouped is returned.
#'  \code{\link{mgidi}} can then be used to compute the mgidi index within each
#'  environment.
#' @param prob The probability for estimating confidence interval for BLUP's
#'   prediction.
#' @param verbose Logical argument. If \code{verbose = FALSE} the code are run
#' silently.
#' @references Mohring, J., E. Williams, and H.-P. Piepho. 2015. Inter-block information:
#' to recover or not to recover it? TAG. Theor. Appl. Genet. 128:1541-54.
#'  \doi{10.1007/s00122-015-2530-0}

#' @return An object of class \code{gamem} or \code{gamem_grouped}, which is a
#'   list with the following items for each element (variable):
#'  * \strong{fixed:} Test for fixed effects.
#'
#'  * \strong{random:} Variance components for random effects.
#'
#'  * \strong{LRT:} The Likelihood Ratio Test for the random effects.
#'
#'  * \strong{BLUPgen:} The estimated BLUPS for genotypes
#'
#'  * \strong{ranef:} The random effects of the model
#'
#'  * \strong{Details:} A tibble with the following data: \code{Ngen}, the
#'  number of genotypes; \code{OVmean}, the grand mean; \code{Min}, the minimum
#'  observed (returning the genotype and replication/block); \code{Max} the
#'  maximum observed, \code{MinGEN} the winner genotype, \code{MaxGEN}, the
#'  loser genotype.
#'
#' * \strong{ESTIMATES:} A tibble with the values for the genotypic variance,
#' block-within-replicate variance (if an alpha-lattice design is used by
#' informing the block in \code{block}), the residual variance and their
#' respective contribution to the phenotypic variance; broad-sence heritability,
#' heritability on the entry-mean basis, genotypic coefficient of variation
#' residual coefficient of variation and ratio between genotypic and residual
#' coefficient of variation.
#'
#'  * \strong{residuals:} The residuals of the model.
#'
#'  * \strong{formula} The formula used to fit the model.
#'
#' @details \code{gamem} analyses data from a one-way genotype testing experiment.
#' By default, a randomized complete block design is used according to the following model:
#' \loadmathjax
#' \mjsdeqn{Y_{ij} = m + g_i + r_j + e_{ij}}
#' where \mjseqn{Y_{ij}} is the response variable of the ith genotype in the \emph{j}th block;
#'  \emph{m} is the grand mean (fixed); \mjseqn{g_i} is the effect of the \emph{i}th genotype
#'  (assumed to be random); \mjseqn{r_j} is the effect of the \emph{j}th replicate (assumed to be fixed);
#'  and \mjseqn{e_{ij}} is the random error.
#'
#' When \code{block} is informed, then a resolvable alpha design is implemented, according to the following model:
#'
#' \mjsdeqn{Y_{ijk} = m + g_i + r_j + b_{jk} + e_{ijk}}
#' where where \mjseqn{y_{ijk}} is the response variable of the \emph{i}th genotype in the
#' \emph{k}th block of the \emph{j}th replicate; \emph{m} is the intercept, \mjseqn{t_i} is
#'  the effect for the \emph{i}th genotype \mjseqn{r_j} is the effect of the \emph{j}th
#'  replicate, \mjseqn{b_{jk}} is the effect of the \emph{k}th incomplete block of
#'  the \emph{j}th replicate, and \mjseqn{e_{ijk}} is the plot error effect
#'  corresponding to \mjseqn{y_{ijk}}.
#'
#' @md
#' @importFrom cowplot draw_label ggdraw
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso \code{\link{get_model_data}} \code{\link{waasb}}
#' @export
#' @examples
#'\donttest{
#' library(metan)
#'
#' # fitting the model considering an RCBD
#' # Genotype as random effects
#'
#'rcbd <- gamem(data_g,
#'              gen = GEN,
#'              rep = REP,
#'              resp = c(PH, ED, EL, CL, CW, KW, NR, TKW, NKE))
#'
#' # Likelihood ratio test for random effects
#' get_model_data(rcbd, "lrt")
#'
#'
#' # Variance components
#' get_model_data(rcbd, "vcomp")
#'
#' # Genetic parameters
#' get_model_data(rcbd, "genpar")
#'
#' # random effects
#' get_model_data(rcbd, "ranef")
#'
#' # Predicted values
#' predict(rcbd)
#'
#' # fitting the model considering an alpha-lattice design
#' # Genotype and block-within-replicate as random effects
#' # Note that block effect was now informed.
#'
#' alpha <- gamem(data_alpha,
#'                gen = GEN,
#'                rep = REP,
#'                block = BLOCK,
#'                resp = YIELD)
#' # Genetic parameters
#' get_model_data(alpha, "genpar")
#'
#' # Random effects
#' get_model_data(alpha, "ranef")
#'}
#'
gamem <- function(.data,
                  gen,
                  rep,
                  resp,
                  block = NULL,
                  by = NULL,
                  prob = 0.05,
                  verbose = TRUE) {
  if (!missing(by)){
    if(length(as.list(substitute(by))[-1L]) != 0){
      stop("Only one grouping variable can be used in the argument 'by'.\nUse 'group_by()' to pass '.data' grouped by more than one variable.", call. = FALSE)
    }
    .data <- group_by(.data, {{by}})
  }
if(is_grouped_df(.data)){
  if(!missing(block)){
    results <-
      .data %>%
      doo(gamem,
          gen = {{gen}},
          rep = {{rep}},
          resp = {{resp}},
          block = {{block}},
          prob = prob,
          verbose = verbose)
  } else{
    results <-
      .data %>%
      doo(gamem,
          gen = {{gen}},
          rep = {{rep}},
          resp = {{resp}},
          prob = prob,
          verbose = verbose)
  }
  return(set_class(results, c("tbl_df",  "gamem_group", "tbl",  "data.frame")))
}
  # RCBD
  if (missing(block) == TRUE) {
    factors  <- .data %>%
      select({{gen}}, {{rep}}) %>%
      mutate(across(everything(), as.factor))
    vars <- .data %>% select({{resp}}, -names(factors))
    vars %<>% select_numeric_cols()
    factors %<>% set_names("GEN", "REP")
    listres <- list()
    nvar <- ncol(vars)
    if (verbose == TRUE) {
      pb <- progress_bar$new(
        format = "Evaluating the variable :what [:bar]:percent (:eta left )",
        clear = FALSE, total = nvar, width = 90)
    }
    model_formula <- "Y ~ REP + (1 | GEN)"
    ran_ef <- c("GEN")
    fix_ef <- c("REP")
    for (var in 1:nvar) {
      data <- factors %>%
        mutate(Y = vars[[var]])
      check_labels(data)
      if(has_na(data)){
        data <- remove_rows_na(data, verbose = verbose) %>% droplevels()
      }
      Ngen <- nlevels(data$GEN)
      Nbloc <- nlevels(data$REP)
      ovmean <- mean(data$Y)
      Complete <- suppressWarnings(suppressMessages(lmerTest::lmer(model_formula, data = data)))
      LRT <- lmerTest::ranova(Complete, reduce.terms = FALSE) %>%
        mutate(model = c("Complete", "Genotype")) %>%
        select(model, everything())
      fixed <- anova(Complete)
      random <- lme4::VarCorr(Complete) %>%
        as.data.frame() %>%
        select(1, 4) %>%
        arrange(grp) %>%
        rename(Group = grp, Variance = vcov)
      regen <- ranef(Complete, condVar = TRUE)
      GV <- as.numeric(random[1, 2])
      RV <- as.numeric(random[2, 2])
      FV <- GV + RV
      h2g <- GV / FV
      h2mg <- GV / (GV + RV / (Nbloc))
      AccuGen <- sqrt(h2mg)
      CVg <- (sqrt(GV) / ovmean) * 100
      CVr <- (sqrt(RV) / ovmean) * 100
      CVratio <- CVg / CVr
      PROB <- ((1 - (1 - prob)) / 2) + (1 - prob)
      t <- qt(PROB, Nbloc)
      Limits <- t * sqrt(((1 - AccuGen) * GV))
      GVper <- (GV / FV) * 100
      RVper <- (RV / FV) * 100
      ESTIMATES <- tibble(
        Parameters = c(
          "Gen_var", "Gen (%)", "Res_var",
          "Res (%)", "Phen_var", "H2", "h2mg",
          "Accuracy", "CVg", "CVr", "CV ratio"
        ),
        Values = c(GV, GVper, RV, RVper, FV, h2g, h2mg, AccuGen, CVg, CVr, CVratio)
      )
      data_factors <- data %>% select_non_numeric_cols()
      BLUPgen <-
        data.frame(GEN = data %>% get_levels(GEN),
                   BLUPg = regen$GEN$`(Intercept)`) %>%
        add_cols(Predicted = BLUPg + ovmean) %>%
        arrange(-Predicted) %>%
        add_cols(Rank = rank(-Predicted),
                 LL = Predicted - Limits,
                 UL = Predicted + Limits) %>%
        column_to_first(Rank)
      ranef <-
        suppressWarnings(
          left_join(data_factors, BLUPgen, by = "GEN") %>%
            select_cols(GEN, REP, BLUPg) %>%
            add_cols(Predicted = BLUPg + left_join(data_factors, means_by(data, REP), by = "REP")$Y)
        )
      min_gen <- data %>%
        group_by(GEN) %>%
        summarise(Y = mean(Y)) %>%
        top_n(1, -Y) %>%
        select(GEN, Y) %>%
        slice(1)
      max_gen <- data %>%
        group_by(GEN) %>%
        summarise(Y = mean(Y)) %>%
        top_n(1, Y) %>%
        select(GEN, Y) %>%
        slice(1)
      max <- data %>%
        top_n(1, Y) %>%
        slice(1)
      min <- data %>%
        top_n(1, -Y) %>%
        slice(1)
      Details <- tibble(
        Parameters = c("Ngen", "OVmean", "Min", "Max", "MinGEN", "MaxGEN"),
        Values = c(
          Ngen,
          round(mean(data$Y), 4),
          paste0(round(min$Y, 4), " (", min$GEN, " in ", min$REP, ")"),
          paste0(round(max$Y, 4), " (", max$GEN, " in ", max$REP, ")"),
          paste0(round(min_gen[1, 2], 4), " (", min_gen$GEN, ")"),
          paste0(round(max_gen[1, 2], 4), " (", max_gen$GEN, ")")
        )
      )
      residuals <- data.frame(fortify.merMod(Complete))
      temp <- structure(list(
        fixed = fixed %>% rownames_to_column("SOURCE") %>% as_tibble(),
        random = as_tibble(random),
        LRT = as_tibble(LRT),
        BLUPgen = BLUPgen,
        ranef = ranef,
        Details = as_tibble(Details),
        ESTIMATES = as_tibble(ESTIMATES),
        residuals = as_tibble(residuals),
        formula = model_formula
      ),
      class = "gamem"
      )
      if (verbose == TRUE) {
        pb$tick(tokens = list(what = names(vars[var])))
      }
      listres[[paste(names(vars[var]))]] <- temp
    }
  }
  # ALPHA-LATTICE
  if (missing(block) == FALSE) {
    factors  <- .data %>%
      select({{gen}}, {{rep}}, {{block}}) %>%
      mutate(across(everything(), as.factor))
    vars <- .data %>% select({{resp}}, -names(factors))
    vars %<>% select_numeric_cols()
    factors %<>% set_names("GEN", "REP", "BLOCK")
    listres <- list()
    nvar <- ncol(vars)
    if (verbose == TRUE) {
      pb <- progress_bar$new(
        format = "Evaluating the variable :what [:bar]:percent",
        clear = FALSE, total = nvar, width = 90)
    }
    model_formula <- "Y ~ (1 | GEN) + REP + (1 | REP:BLOCK)"
    ran_ef <- c("GEN, BLOCK(REP)")
    fix_ef <- c("REP")
    for (var in 1:nvar) {
      data <- factors %>%
        mutate(Y = vars[[var]])
      check_labels(data)
      if(has_na(data)){
        data <- remove_rows_na(data, verbose = verbose) %>% droplevels()
        has_text_in_num(data)
      }
      Ngen <- nlevels(data$GEN)
      Nbloc <- nlevels(data$REP)
      ovmean <- mean(data$Y)
      Complete <- suppressWarnings(suppressMessages(lmerTest::lmer(model_formula, data = data)))
      LRT <- lmerTest::ranova(Complete, reduce.terms = FALSE) %>%
        mutate(model = c("Complete", "Genotype", "rep:block")) %>%
        select(model, everything())
      fixed <- anova(Complete)
      random <- lme4::VarCorr(Complete) %>%
        as.data.frame() %>%
        select(1, 4) %>%
        arrange(grp) %>%
        rename(Group = grp, Variance = vcov)
      regen <- ranef(Complete, condVar = TRUE)
      GV <- as.numeric(random[1, 2])
      BRV <- as.numeric(random[2, 2])
      RV <- as.numeric(random[3, 2])
      FV <- GV + RV + BRV
      h2g <- GV / FV
      vv <- attr(regen$GEN, "postVar")
      vblup <- 2 * mean(vv)
      sg2 <- c(lme4::VarCorr(Complete)[["GEN"]])
      # H^2 measure proposed by Cullis
      h2mg <- 1 - (vblup / 2 / sg2)
      AccuGen <- sqrt(h2mg)
      CVg <- (sqrt(GV) / ovmean) * 100
      CVr <- (sqrt(RV) / ovmean) * 100
      CVratio <- CVg / CVr
      PROB <- ((1 - (1 - prob)) / 2) + (1 - prob)
      t <- qt(PROB, Nbloc)
      Limits <- t * sqrt(((1 - AccuGen) * GV))
      GVper <- (GV / FV) * 100
      BRper <- (BRV / FV) * 100
      RVper <- (RV / FV) * 100
      ESTIMATES <- tibble(
        Parameters = c(
          "Gen_var", "Gen (%)", "rep:block_var", "rep:block (%)", "Res_var",
          "Res (%)", "Phen_var", "H2", "h2mg", "Accuracy", "CVg", "CVr", "CV ratio"
        ),
        Values = c(GV, GVper, BRV, BRper, RV, RVper, FV, h2g, h2mg, AccuGen, CVg, CVr, CVratio)
      )
      data_factors <- data %>% select_non_numeric_cols()
      BLUPgen <-
        data.frame(GEN = data %>% get_levels(GEN),
                   BLUPg = regen$GEN$`(Intercept)`) %>%
        add_cols(Predicted = BLUPg + ovmean) %>%
        arrange(-Predicted) %>%
        add_cols(Rank = rank(-Predicted),
                 LL = Predicted - Limits,
                 UL = Predicted + Limits) %>%
        column_to_first(Rank)
      blupBWR <- data.frame(Names = rownames(regen$`REP:BLOCK`)) %>%
        separate(Names, into = c("REP", "BLOCK"), sep = ":") %>%
        add_cols(BLUPbre = regen$`REP:BLOCK`[[1]]) %>%
        to_factor(1:2)
      ranef <-
        suppressWarnings(
          left_join(data_factors, BLUPgen, by = "GEN") %>%
            left_join(blupBWR, by = c("REP", "BLOCK")) %>%
            select_cols(GEN, REP, BLOCK, BLUPg, BLUPbre) %>%
            add_cols(`BLUPg+bre` =  BLUPg + BLUPbre,
                     Predicted = `BLUPg+bre` + left_join(data_factors, means_by(data, REP), by = "REP")$Y)
        )
      min_gen <- data %>%
        group_by(GEN) %>%
        summarise(Y = mean(Y)) %>%
        top_n(1, -Y) %>%
        select(GEN, Y) %>%
        slice(1)
      max_gen <- data %>%
        group_by(GEN) %>%
        summarise(Y = mean(Y)) %>%
        top_n(1, Y) %>%
        select(GEN, Y) %>%
        slice(1)
      max <- data %>%
        top_n(1, Y) %>%
        slice(1)
      min <- data %>%
        top_n(1, -Y) %>%
        slice(1)
      Details <- tibble(
        Parameters = c("Ngen", "OVmean", "Min", "Max", "MinGEN", "MaxGEN"),
        Values = c(
          Ngen,
          round(mean(data$Y), 4),
          paste0(round(min$Y, 4), " (", min$GEN, " in ", min$BLOCK, " of ", min$REP, ")"),
          paste0(round(max$Y, 4), " (", max$GEN, " in ", max$BLOCK, " of ", max$REP, ")"),
          paste0(round(min_gen[1, 2], 4), " (", min_gen$GEN, ")"),
          paste0(round(max_gen[1, 2], 4), " (", max_gen$GEN, ")")
        )
      )
      residuals <- as_tibble(fortify.merMod(Complete))
      temp <- structure(list(
        fixed = fixed %>% rownames_to_column("SOURCE") %>% as_tibble(),
        random = as_tibble(random),
        LRT = as_tibble(LRT),
        BLUPgen = BLUPgen,
        ranef = ranef,
        Details = as_tibble(Details),
        ESTIMATES = as_tibble(ESTIMATES),
        residuals = as_tibble(residuals),
        formula = model_formula
      ),
      class = "gamem"
      )
      if (verbose == TRUE) {
        pb$tick(tokens = list(what = names(vars[var])))
      }
      listres[[paste(names(vars[var]))]] <- temp

    }
  }
  if (verbose == TRUE) {
    message("Method: REML/BLUP\n", appendLF = FALSE)
    message("Random effects: ", ran_ef, "\n", appendLF = FALSE)
    message("Fixed effects: ", fix_ef, "\n", appendLF = FALSE)
    message("Denominador DF: Satterthwaite's method\n", appendLF = FALSE)
    cat("---------------------------------------------------------------------------\n")
    cat("P-values for Likelihood Ratio Test of the analyzed traits\n")
    cat("---------------------------------------------------------------------------\n")
    print.data.frame(sapply(listres, function(x){
      x$LRT[["Pr(>Chisq)"]]
    }) %>%
      as.data.frame() %>%
      add_cols(model = listres[[1]][["LRT"]][["model"]]) %>%
      column_to_first(model), row.names = FALSE, digits = 3)
    cat("---------------------------------------------------------------------------\n")
    if (length(which(unlist(lapply(listres, function(x) {
      x[["LRT"]] %>% dplyr::filter(model == "Genotype") %>% pull(`Pr(>Chisq)`)
    })) > prob)) > 0) {
      cat("Variables with nonsignificant Genotype effect\n")
      cat(names(which(unlist(lapply(listres, function(x) {
        x[["LRT"]][which(x[["LRT"]][[1]] == "Genotype"), 7] %>% pull()
      })) > prob)), "\n")
      cat("---------------------------------------------------------------------------\n")
    } else {
      cat("All variables with significant (p < 0.05) genotype effect\n")
    }
  }
  invisible(structure(listres, class = "gamem"))
}








#' Print an object of class gamem
#'
#' Print the \code{gamem} object in two ways. By default, the results are shown
#' in the R console. The results can also be exported to the directory.
#'
#'
#' @param x An object fitted with the function \code{\link{gamem}} .
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported to
#'   the working directory
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print gamem
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' alpha <- gamem(data_alpha,
#'   gen = GEN,
#'   rep = REP,
#'   block = BLOCK,
#'   resp = YIELD
#' )
#'
#' print(alpha)
#' }
print.gamem <- function(x, export = FALSE, file.name = NULL, digits = 4, ...) {
  if (!class(x) == "gamem") {
    stop("The object must be of class 'gamem'")
  }
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "gamem print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    cat("Fixed-effect anova table\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$fixed, ...)
    cat("---------------------------------------------------------------------------\n")
    cat("Variance components for random effects\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$random, ...)
    cat("---------------------------------------------------------------------------\n")
    cat("Likelihood ratio test for random effects\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$LRT, ...)
    cat("---------------------------------------------------------------------------\n")
    cat("Details of the analysis\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$Details, ...)
    cat("---------------------------------------------------------------------------\n")
    cat("Genetic parameters\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$ESTIMATES, ...)
    cat("\n\n\n")
  }
  if (export == TRUE) {
    sink()
  }
}





#' Predict method for gamem fits
#'
#' Obtains predictions from an object fitted with \code{\link{gamem}}.
#'
#'
#' @param object An object of class \code{gamem}
#' @param ... Currently not used
#' @return A tibble with the predicted values for each variable in the model
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method predict gamem
#' @export
#' @examples
#'\donttest{
#' library(metan)
#'model <- gamem(data_g,
#'               gen = GEN,
#'               rep = REP,
#'               resp = everything())
#' predict(model)
#' }
#'
predict.gamem <- function(object, ...) {
  factors <- object[[1]][["ranef"]] %>% select_non_numeric_cols()
  numeric <- sapply(object, function(x){
    x[["ranef"]][["Predicted"]]
  })
  return(cbind(factors, numeric) %>% as_tibble())
}



#' Several types of residual plots
#'
#' Residual plots for a output model of class \code{gamem}. Six types of plots
#' are produced: (1) Residuals vs fitted, (2) normal Q-Q plot for the residuals,
#' (3) scale-location plot (standardized residuals vs Fitted Values), (4)
#' standardized residuals vs Factor-levels, (5) Histogram of raw residuals and
#' (6) standardized residuals vs observation order. For a \code{waasb} object,
#' normal Q-Q plot for random effects may also be obtained declaring \code{type
#' = 're'}
#'
#'
#' @param x An object of class \code{gamem}.
#' @param var The variable to plot. Defaults to \code{var = 1} the first
#'   variable of \code{x}.
#' @param type One of the \code{"res"} to plot the model residuals (default),
#'   \code{type = 're'} to plot normal Q-Q plots for the random effects, or
#'   \code{"vcomp"} to create a bar plot with the variance components.
#' @param position The position adjustment when \code{type = "vcomp"}. Defaults
#'   to \code{"fill"}, which shows relative proportions at each trait by
#'   stacking the bars and then standardizing each bar to have the same height.
#'   Use \code{position = "stack"} to plot the phenotypic variance for each
#'   trait.
#' @param rotate Logical argument. If \code{rotate = TRUE} the plot is rotated,
#'   i.e., traits in y axis and value in the x axis.
#' @param conf Level of confidence interval to use in the Q-Q plot (0.95 by
#' default).
#' @param out How the output is returned. Must be one of the 'print' (default)
#' or 'return'.
#' @param n.dodge The number of rows that should be used to render the x labels.
#'   This is useful for displaying labels that would otherwise overlap.
#' @param check.overlap Silently remove overlapping labels, (recursively)
#'   prioritizing the first, last, and middle labels.
#' @param labels Logical argument. If \code{TRUE} labels the points outside
#' confidence interval limits.
#' @param plot_theme The graphical theme of the plot. Default is
#'   \code{plot_theme = theme_metan()}. For more details, see
#'   \code{\link[ggplot2]{theme}}.
#' @param alpha The transparency of confidence band in the Q-Q plot. Must be a
#' number between 0 (opaque) and 1 (full transparency).
#' @param fill.hist The color to fill the histogram. Default is 'gray'.
#' @param col.hist The color of the border of the the histogram. Default is
#' 'black'.
#' @param col.point The color of the points in the graphic. Default is 'black'.
#' @param col.line The color of the lines in the graphic. Default is 'red'.
#' @param col.lab.out The color of the labels for the 'outlying' points.
#' @param size.line The size of the line in graphic. Defaults to 0.7.
#' @param size.text The size for the text in the plot. Defaults to 10.
#' @param width.bar The width of the bars if \code{type = "contribution"}.
#' @param size.lab.out The size of the labels for the 'outlying' points.
#' @param size.tex.lab The size of the text in axis text and labels.
#' @param size.shape The size of the shape in the plots.
#' @param bins The number of bins to use in the histogram. Default is 30.
#' @param which Which graphics should be plotted. Default is \code{which =
#' c(1:4)} that means that the first four graphics will be plotted.
#' @param ncol,nrow The number of columns and rows of the plot pannel. Defaults
#'   to \code{NULL}
#' @param align Specifies whether graphs in the grid should be horizontally
#'   (\code{"h"}) or vertically (\code{"v"}) aligned. \code{"hv"} (default)
#'   align in both directions, \code{"none"} do not align the plot.
#' @param ... Additional arguments passed on to the function
#'   \code{\link[cowplot]{plot_grid}}
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method plot gamem
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' model <- gamem(data_g,
#'                gen = GEN,
#'                rep = REP,
#'                resp = PH)
#' plot(model)
#'}
#'
plot.gamem <- function(x,
                       var = 1,
                       type = "res",
                       position = "fill",
                       rotate = FALSE,
                       conf = 0.95,
                       out = "print",
                       n.dodge = 1,
                       check.overlap = FALSE,
                       labels = FALSE,
                       plot_theme = theme_metan(),
                       alpha = 0.2,
                       fill.hist = "gray",
                       col.hist = "black",
                       col.point = "black",
                       col.line = "red",
                       col.lab.out = "red",
                       size.line = 0.7,
                       size.text = 10,
                       width.bar = 0.75,
                       size.lab.out = 2.5,
                       size.tex.lab = 10,
                       size.shape = 1.5,
                       bins = 30,
                       which = c(1:4),
                       ncol = NULL,
                       nrow = NULL,
                       align = "hv",
                       ...) {
  if(!type  %in% c("res", 're', "vcomp")){
    stop("Argument type = '", match.call()[["type"]], "' invalid. Use one of 'res', 're', or 'vcomp'", call. = FALSE)
  }
  if(type %in% c("vcomp", "re") && !class(x)  %in% c("waasb", "gamem")){
    stop("Arguments 're' and 'vcomp' valid for objects of class 'waasb' or 'gamem'. ")
  }
  if(is.numeric(var)){
    var_name <- names(x)[var]
  } else{
    var_name <- var
  }
  if(!var_name %in% names(x)){
    stop("Variable not found in ", match.call()[["x"]] , call. = FALSE)
  }
  if (type == "re" & max(which) >= 5) {
    stop("When type =\"re\", 'which' must be a value between 1 and 4")
  }
  if(type == "vcomp"){
    list <- lapply(x, function(x){
      x[["random"]] %>% select(Group, Variance)
    })
    vcomp <- suppressWarnings(
      lapply(seq_along(list),
             function(i){
               set_names(list[[i]], "Group", names(list)[i])
             }) %>%
        reduce(full_join, by = "Group") %>%
        pivot_longer(-Group))
    p1 <-
      ggplot(vcomp, aes(x = name, y = value, fill = Group)) +
      geom_bar(stat = "identity",
               position = position,
               color = "black",
               size = size.line,
               width = width.bar) +
      scale_y_continuous(expand = expansion(c(0, ifelse(position == "fill", 0, 0.05)))) +
      scale_x_discrete(guide = guide_axis(n.dodge = n.dodge, check.overlap = check.overlap)) +
      theme_bw()+
      theme(legend.position = "bottom",
            axis.ticks = element_line(size = size.line),
            axis.ticks.length = unit(0.2, "cm"),
            panel.grid = element_blank(),
            legend.title = element_blank(),
            strip.background = element_rect(fill = NA),
            text = element_text(size = size.text, colour = "black"),
            axis.text = element_text(size = size.text, colour = "black")) +
      theme(legend.position = "bottom") +
      labs(x = "Traits",
      y = ifelse(position == "fill", "Proportion of phenotypic variance", "Phenotypic variance"))
    if(rotate == TRUE){
      p1 <- p1 + coord_flip()
    }
    return(p1)
  }
  if (type == "res") {
    x <- x[[var]]
    df <- data.frame(x$residuals)
    df$id <- rownames(df)
    df <- data.frame(df[order(df$.scresid), ])
    P <- ppoints(nrow(df))
    df$z <- qnorm(P)
    n <- nrow(df)
    Q.x <- quantile(df$.scresid, c(0.25, 0.75))
    Q.z <- qnorm(c(0.25, 0.75))
    b <- diff(Q.x)/diff(Q.z)
    coef <- c(Q.x[1] - b * Q.z[1], b)
    zz <- qnorm(1 - (1 - conf)/2)
    SE <- (coef[2]/dnorm(df$z)) * sqrt(P * (1 - P)/n)
    fit.value <- coef[1] + coef[2] * df$z
    df$upper <- fit.value + zz * SE
    df$lower <- fit.value - zz * SE
    df$label <- ifelse(df$.scresid > df$.scresid | df$.scresid <
                         df$lower, rownames(df), "")
    df$factors <- paste(df$ENV, df$GEN)
    # Residuals vs .fitted
    p1 <- ggplot(df, aes(.fitted, .resid)) +
      geom_point(col = col.point, size = size.shape) +
      geom_smooth(se = F, method = "loess", col = col.line) +
      geom_hline(yintercept = 0, linetype = 2, col = "gray") +
      labs(x = "Fitted values", y = "Residual") +
      ggtitle("Residual vs fitted") + plot_theme %+replace%
      theme(axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            plot.title = element_text(size = size.tex.lab, hjust = 0, vjust = 1))
    if (labels != FALSE) {
      p1 <- p1 +
        ggrepel::geom_text_repel(aes(.fitted, .resid, label = (label)),
                                 color = col.lab.out,
                                 size = size.lab.out)
    } else {
      p1 <- p1
    }
    # normal qq
    p2 <- ggplot(df, aes(z, .scresid)) +
      geom_point(col = col.point, size = size.shape) +
      geom_abline(intercept = coef[1],
                  slope = coef[2],
                  size = 1,
                  col = col.line) +
      geom_ribbon(aes_(ymin = ~lower, ymax = ~upper),
                  alpha = 0.2) +
      labs(x = "Theoretical quantiles", y = "Sample quantiles") +
      ggtitle("Normal Q-Q") +
      plot_theme %+replace%
      theme(axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            plot.title = element_text(size = size.tex.lab, hjust = 0, vjust = 1))
    if (labels != FALSE) {
      p2 <- p2 + ggrepel::geom_text_repel(aes(z, .scresid, label = (label)),
                                          color = col.lab.out,
                                          size = size.lab.out)
    } else {
      p2 <- p2
    }
    # scale-location
    p3 <- ggplot(df, aes(.fitted, sqrt(abs(.resid)))) +
      geom_point(col = col.point, size = size.shape) +
      geom_smooth(se = F, method = "loess", col = col.line) +
      labs(x = "Fitted Values", y = expression(sqrt("|Standardized residuals|"))) +
      ggtitle("Scale-location") +
      plot_theme %+replace%
      theme(axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            plot.title = element_text(size = size.tex.lab, hjust = 0, vjust = 1))
    if (labels != FALSE) {
      p3 <- p3 + ggrepel::geom_text_repel(aes(.fitted, sqrt(abs(.resid)),
                                              label = (label)),
                                          color = col.lab.out,
                                          size = size.lab.out)
    } else {
      p3 <- p3
    }
    # Residuals vs Factor-levels
    p4 <- ggplot(df, aes(factors, .scresid)) +
      geom_point(col = col.point, size = size.shape) +
      geom_hline(yintercept = 0, linetype = 2, col = "gray") +
      labs(x = "Factor levels", y = "Standardized residuals") +
      ggtitle("Residuals vs factor-levels") +
      plot_theme %+replace%
      theme(axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            panel.grid.major.x = element_blank(),
            axis.text.x = element_text(color = "white"),
            plot.title = element_text(size = size.tex.lab, hjust = 0, vjust = 1))
    if (labels != FALSE) {
      p4 <- p4 + ggrepel::geom_text_repel(aes(factors,
                                              .scresid, label = (label)),
                                          color = col.lab.out,
                                          size = size.lab.out)
    } else {
      p4 <- p4
    }
    # Histogram of residuals
    p5 <- ggplot(df, aes(x = .resid)) +
      geom_histogram(bins = bins,
                     colour = col.hist,
                     fill = fill.hist,
                     aes(y = ..density..)) +
      stat_function(fun = dnorm,
                    color = col.line,
                    size = 1,
                    args = list(mean = mean(df$.resid),
                                sd = sd(df$.resid))) +
      labs(x = "Raw residuals", y = "Density") +
      ggtitle("Histogram of residuals") +
      plot_theme %+replace%
      theme(axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            plot.title = element_text(size = size.tex.lab, hjust = 0, vjust = 1))
    # Residuals vs order
    p6 <- ggplot(df, aes(as.numeric(id), .scresid, group = 1)) +
      geom_point(col = col.point, size = size.shape) +
      geom_line(col = col.line) +
      geom_hline(yintercept = 0,
                 linetype = 2,
                 col = col.line) +
      labs(x = "Observation order", y = "Standardized residuals") +
      ggtitle("Residuals vs observation order") +
      plot_theme %+replace%
      theme(axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            plot.title = element_text(size = size.tex.lab, hjust = 0, vjust = 1))
    p7 <- ggplot(df, aes(.fitted, Y)) +
      geom_point(col = col.point, size = size.shape) +
      facet_wrap(~GEN) + geom_abline(intercept = 0, slope = 1, col = col.line) +
      labs(x = "Fitted values", y = "Observed values") +
      ggtitle("1:1 line plot") +
      plot_theme %+replace%
      theme(axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            plot.title = element_text(size = size.tex.lab, hjust = 0, vjust = 1),
            panel.spacing = unit(0, "cm"))
    plots <- list(p1, p2, p3, p4, p5, p6, p7)
    p1 <-
      plot_grid(plotlist = plots[c(which)],
                ncol = ncol,
                nrow = nrow,
                align = align,
                ...)
    title <- ggdraw() +
      draw_label(var_name,
                 fontface = 'bold',
                 x = 0,
                 hjust = 0) +
      theme(plot.margin = margin(0, 0, 0, 7))
    p1 <-
      plot_grid(title, p1,
                ncol = 1,
                rel_heights = c(0.05, 1))
    return(p1)
  }
  if (type == "re") {
    x <- x[[var]]
    blups <-
      x$ranef %>%
      select_cols(contains("BLUP"))
    fact <-x$ranef %>% select_non_numeric_cols()
    qlist <- list()
    for (i in 1:ncol(blups)) {
      df <-
        data.frame(blups[i]) %>%
        distinct_all() %>%
        rowid_to_column(var = "id") %>%
        arrange(across(2))
      P <- ppoints(nrow(df))
      df$z <- qnorm(P)
      n <- nrow(df)
      Q.x <- quantile(df[[2]], c(0.25, 0.75))
      Q.z <- qnorm(c(0.25, 0.75))
      b <- diff(Q.x)/diff(Q.z)
      coef <- c(Q.x[1] - b * Q.z[1], b)
      zz <- qnorm(1 - (1 - conf)/2)
      SE <- (coef[2]/dnorm(df$z)) * sqrt(P * (1 - P)/n)
      fit.value <- coef[1] + coef[2] * df$z
      df %<>% add_cols(upper = fit.value + zz * SE,
                       lower = fit.value - zz * SE,
                       label = ifelse(df[[2]] > upper | df[[2]] < lower, id, ""),
                       intercept = coef[1],
                       slope = coef[2],
                       var = paste(names(blups[i]))
      ) %>%
        set_names("id",    "blup", "z",     "upper", "lower", "label", "intercept", "slope", "var")
      qlist[[paste(names(blups[i]))]] <- df
    }

    df <- do.call(rbind, qlist)
    # normal qq GEI effects
    p1 <- ggplot(df, aes(z, blup)) +
      geom_point(col = col.point, size = size.shape) +
      geom_abline(aes(intercept = intercept,
                      slope = slope),
                  size = 1, col = col.line) +
      geom_ribbon(aes_(ymin = ~lower, ymax = ~upper),
                  alpha = 0.2) +
      labs(x = "Theoretical quantiles", y = "Sample quantiles")+
      facet_wrap( ~var, scales = "free") +
      plot_theme %+replace%
      theme(axis.text = element_text(size = size.tex.lab, colour = "black"),
            axis.title = element_text(size = size.tex.lab, colour = "black"),
            plot.title.position = "plot",
            plot.title = element_text(size = size.tex.lab + 1, hjust = 0, vjust = 1, face = "bold"))
    if (labels != FALSE) {
      p1 <- p1 + ggrepel::geom_text_repel(aes(z, blup, label = (label)),
                                          color = col.lab.out,
                                          size = size.lab.out)+
        labs(title = var_name)
    } else {
      p1 <- p1 + labs(title = var_name)
    }
    return(p1)
  }
}
