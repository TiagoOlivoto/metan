#' Genotype analysis by fixed-effect models
#' @description
#' `r badge('stable')`
#'
#' One-way analysis of variance of genotypes conducted in both randomized
#' complete block and alpha-lattice designs.
#'
#'
#' @param .data The dataset containing the columns related to, Genotypes,
#' replication/block and response variable(s).
#' @param gen The name of the column that contains the levels of the genotypes, that will
#' be treated as random effect.
#' @param rep The name of the column that contains the levels of the replications
#' (assumed to be fixed).
#' @param resp The response variable(s). To analyze multiple variables in a
#' single procedure a vector of variables may be used. For example `resp =
#' c(var1, var2, var3)`. Select helpers are also allowed.
#' @param block Defaults to `NULL`. In this case, a randomized complete
#'   block design is considered. If block is informed, then a resolvable
#'   alpha-lattice design (Patterson and Williams, 1976) is employed.
#'   **All effects, except the error, are assumed to be fixed.** Use the
#'   function [gamem()] to analyze a one-way trial with mixed-effect
#'   models.
#'@param by One variable (factor) to compute the function by. It is a shortcut
#'  to [dplyr::group_by()].This is especially useful, for example,
#'  when the researcher want to fit a fixed-effect model for each environment.
#'  In this case, an object of class gafem_grouped is returned.
#'  [mgidi()] can then be used to compute the mgidi index within each
#'  environment.
#' @param prob The error probability. Defaults to 0.05.
#' @param verbose Logical argument. If `verbose = FALSE` the code are run
#' silently.
#' @references Patterson, H.D., and E.R. Williams. 1976. A new class of
#' resolvable incomplete block designs. Biometrika 63:83-92.

#' @return A list where each element is the result for one variable containing
#'   the following objects:
#'
#'  * **anova:** The one-way ANOVA table.
#'  * **model:** The model with of `lm`.
#'  * **augment:** Information about each observation in the dataset. This
#'  includes predicted values in the `fitted` column, residuals in the
#'  `resid` column, standardized residuals in the `stdres` column,
#'  the diagonal of the 'hat' matrix in the `hat`, and standard errors for
#'  the fitted values in the `se.fit` column.
#'  * **hsd:** The Tukey's 'Honest Significant Difference' for genotype
#'  effect.
#'  * **details:** A tibble with the following data: `Ngen`, the
#'  number of genotypes; `OVmean`, the grand mean; `Min`, the minimum
#'  observed (returning the genotype and replication/block); `Max` the
#'  maximum observed, `MinGEN` the loser winner genotype, `MaxGEN`,
#'  the winner genotype.
#'
#' @details `gafem` analyses data from a one-way genotype testing
#'   experiment. By default, a randomized complete block design is used
#'   according to the following model:
#' \loadmathjax
#' \mjsdeqn{Y_{ij} = m + g_i + r_j + e_{ij}}
#' where \mjseqn{Y_{ij}} is the response variable of the ith genotype in the
#' *j*th block; *m* is the grand mean (fixed); \mjseqn{g_i} is the effect
#' of the *i*th genotype; \mjseqn{r_j} is the effect of the *j*th
#' replicate; and \mjseqn{e_{ij}} is the random error.
#'
#' When `block` is informed, then a resolvable alpha design is implemented,
#' according to the following model:
#'
#' \mjsdeqn{Y_{ijk} = m + g_i + r_j + b_{jk} + e_{ijk}}
#'  where where \mjseqn{y_{ijk}}
#' is the response variable of the *i*th genotype in the *k*th block
#' of the *j*th replicate; *m* is the intercept, \mjseqn{t_i} is the
#' effect for the *i*th genotype \mjseqn{r_j} is the effect of the *j*th
#' replicate, \mjseqn{b_{jk}} is the effect of the *k*th incomplete block of
#' the *j*th replicate, and \mjseqn{e_{ijk}} is the plot error effect
#' corresponding to \mjseqn{y_{ijk}}. All effects, except the random error are
#' assumed to be fixed.
#'
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso [get_model_data()] [gamem()]
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' # RCBD
#'rcbd <- gafem(data_g,
#'              gen = GEN,
#'              rep = REP,
#'              resp = c(PH, ED, EL, CL, CW))
#'
#' # Fitted values
#' get_model_data(rcbd)
#'
#' # ALPHA-LATTICE DESIGN
#' alpha <- gafem(data_alpha,
#'               gen = GEN,
#'               rep = REP,
#'               block = BLOCK,
#'               resp = YIELD)
#'
#' # Fitted values
#' get_model_data(alpha)
#'
#'}
#'
gafem <- function(.data,
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
        doo(gafem,
            gen = {{gen}},
            rep = {{rep}},
            resp = {{resp}},
            block = {{block}},
            prob = prob,
            verbose = verbose)
    } else{
      results <-
        .data %>%
        doo(gafem,
            gen = {{gen}},
            rep = {{rep}},
            resp = {{resp}},
            prob = prob,
            verbose = verbose)
    }
    return(set_class(results, c("tbl_df",  "gafem_group", "tbl",  "data.frame")))
  }
  if (missing(block) == TRUE) {
    factors  <-
      .data %>%
      select({{gen}}, {{rep}}) %>%
      mutate(across(everything(), as.factor))
    vars <- .data %>% select({{resp}}, -names(factors))
    vars %<>% select_numeric_cols()
    factors %<>% set_names("GEN", "REP")
    listres <- list()
    nvar <- ncol(vars)
    if (verbose == TRUE) {
      pb <- progress(max = nvar, style = 4)
    }
    for (var in 1:nvar) {
      data <- factors %>%
        mutate(Y = vars[[var]])
      if(has_na(data)){
        data <- remove_rows_na(data)
        has_text_in_num(data)
      }
      Ngen <- nlevels(data$GEN)
      Nbloc <- nlevels(data$REP)
      ovmean <- mean(data$Y)
      model <- aov(Y ~ REP + GEN, data = data)
      anova <- anova(model) %>%
        rownames_to_column("Source") %>%
        as_tibble()
      min_gen <-
        data %>%
        mean_by(GEN) %>%
        top_n(1, -Y) %>%
        select(GEN, Y) %>%
        slice(1)
      max_gen <-
        data %>%
        mean_by(GEN) %>%
        top_n(1, Y) %>%
        select(GEN, Y) %>%
        slice(1)
      max <-
        data %>%
        top_n(1, Y) %>%
        slice(1)
      min <-
        data %>%
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
      influence <- lm.influence(model)
      augment <- model$model %>%
        add_cols(hat = influence$hat,
                 sigma = influence$sigma,
                 fitted = predict(model),
                 resid = residuals(model),
                 stdres = rstandard(model),
                 se.fit = predict(model, se.fit = TRUE)$se.fit) %>%
        add_cols(factors = concatenate(., GEN, REP, pull = TRUE)) %>%
        column_to_first(GEN, REP)
      temp <- structure(list(
        anova = anova,
        model = model,
        augment = augment,
        hsd = tukey_hsd(model, which = "GEN"),
        details = Details
      ),
      class = "gafem"
      )
      listres[[paste(names(vars[var]))]] <- temp
      if (verbose == TRUE) {
        run_progress(pb,
                     actual = var,
                     text = paste("Evaluating trait", names(vars[var])))

      }
    }
  }
  # ALPHA-LATTICE
  if (missing(block) == FALSE) {
    factors  <-
      .data %>%
      select({{gen}}, {{rep}}, {{block}}) %>%
      mutate(across(everything(), as.factor))
    vars <- .data %>% select({{resp}}, -names(factors))
    vars %<>% select_numeric_cols()
    factors %<>% set_names("GEN", "REP", "BLOCK")
    listres <- list()
    nvar <- ncol(vars)
    if (verbose == TRUE) {
      pb <- progress(max = nvar, style = 4)
    }
    for (var in 1:nvar) {
      data <- factors %>%
        mutate(Y = vars[[var]])
      if(has_na(data)){
        data <- remove_rows_na(data)
        has_text_in_num(data)
      }
      Ngen <- nlevels(data$GEN)
      Nbloc <- nlevels(data$REP)
      ovmean <- mean(data$Y)
      model <- aov(Y ~ REP + GEN + REP:BLOCK, data = data)
      anova <- anova(model) %>%
        rownames_to_column("Source") %>%
        as_tibble()
      anova[3, 1] <- "BLOCK(REP)"
      min_gen <-
        data %>%
        mean_by(GEN) %>%
        top_n(1, -Y) %>%
        select(GEN, Y) %>%
        slice(1)
      max_gen <-
        data %>%
        mean_by(GEN) %>%
        top_n(1, Y) %>%
        select(GEN, Y) %>%
        slice(1)
      max <-
        data %>%
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
      influence <- lm.influence(model)
      augment <- model$model %>%
        add_cols(hat = influence$hat,
                 sigma = influence$sigma,
                 fitted = predict(model),
                 resid = residuals(model),
                 stdres = rstandard(model),
                 se.fit = predict(model, se.fit = TRUE)$se.fit) %>%
        add_cols(factors = concatenate(., GEN, REP, BLOCK, pull = TRUE)) %>%
        column_to_first(GEN, REP, BLOCK)
      temp <- structure(list(
        anova = anova,
        model = model,
        augment = augment,
        hsd = tukey_hsd(model, which = "GEN"),
        details = Details
      ),
      class = "gafem"
      )
      listres[[paste(names(vars[var]))]] <- temp
      if (verbose == TRUE) {
        run_progress(pb,
                     actual = var,
                     text = paste("Evaluating trait", names(vars[var])))
      }
    }
  }
  if (verbose == TRUE) {
    design <- ifelse(missing(block), "rcbd", "alpha")
    cat("---------------------------------------------------------------------------\n")
    switch(design,
           rcbd =  cat("One-way ANOVA table (Randomized complete block design)\n"),
           alpha = cat("One-way ANOVA table (Alpha-lattice design)\n"))
    cat("---------------------------------------------------------------------------\n")
    print.data.frame(sapply(listres, function(x){
      x$anova[["Pr(>F)"]]
    }) %>%
      as.data.frame() %>%
      add_cols(model = listres[[1]][["anova"]][["Source"]]) %>%
      column_to_first(model), row.names = FALSE, digits = 3)
    if (length(which(unlist(lapply(listres, function(x) {
      x[["anova"]][1, 6]
    })) > prob)) > 0) {
      cat("---------------------------------------------------------------------------\n")
      cat("Variables with nonsignificant genotype effect\n")
      cat(names(which(unlist(lapply(listres, function(x) {
        as.numeric(x[["anova"]][2, 6])
      })) > prob)), "\n")
      cat("---------------------------------------------------------------------------\n\n")
    }
  }
  invisible(structure(listres, class = "gafem"))
}

#' Several types of residual plots
#'
#' Residual plots for a output model of class `gafem`. Seven types
#' of plots are produced: (1) Residuals vs fitted, (2) normal Q-Q plot for the
#' residuals, (3) scale-location plot (standardized residuals vs Fitted Values),
#' (4) standardized residuals vs Factor-levels, (5) Histogram of raw residuals
#' and (6) standardized residuals vs observation order, and (7) 1:1 line plot.
#'
#'
#' @param x An object of class `gafem`.
#' @param ... Additional arguments passed on to the function
#'   [residual_plots()]
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method plot gafem
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' model <- gafem(data_g, GEN, REP, PH)
#'
#' plot(model)
#' plot(model,
#'      which = c(3, 5),
#'      nrow = 2,
#'      labels = TRUE,
#'      size.lab.out = 4)
#' }
#'
plot.gafem <- function(x, ...) {
  residual_plots(x,  ...)
}
