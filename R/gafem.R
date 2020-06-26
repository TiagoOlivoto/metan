#' Genotype analysis by fixed-effect models
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
#' single procedure a vector of variables may be used. For example \code{resp =
#' c(var1, var2, var3)}. Select helpers are also allowed.
#' @param prob The error probability. Defaults to 0.05.
#' @param block Defaults to \code{NULL}. In this case, a randomized complete
#'   block design is considered. If block is informed, then a resolvable
#'   alpha-lattice design (Patterson and Williams, 1976) is employed.
#'   \strong{All effects, except the error, are assumed to be fixed.} Use the
#'   function \code{\link{gamem}} to analyze a one-way trial with mixed-effect
#'   models.
#' @param verbose Logical argument. If \code{verbose = FALSE} the code are run
#' silently.
#' @references Patterson, H.D., and E.R. Williams. 1976. A new class of
#' resolvable incomplete block designs. Biometrika 63:83-92.

#' @return A list where each element is the result for one variable containing
#'   the following objects:
#'
#'  * \strong{anova:} The one-way ANOVA table.
#'  * \strong{model:} The model with of \code{lm}.
#'  * \strong{augment:} Information about each observation in the dataset. This
#'  includes predicted values in the \code{fitted} column, residuals in the
#'  \code{resid} column, standardized residuals in the \code{stdres} column,
#'  the diagonal of the 'hat' matrix in the \code{hat}, and standard errors for
#'  the fitted values in the \code{se.fit} column.
#'  * \strong{hsd:} The Tukey's 'Honest Significant Difference' for genotype
#'  effect.
#'  * \strong{details:} A tibble with the following data: \code{Ngen}, the
#'  number of genotypes; \code{OVmean}, the grand mean; \code{Min}, the minimum
#'  observed (returning the genotype and replication/block); \code{Max} the
#'  maximum observed, \code{MinGEN} the loser winner genotype, \code{MaxGEN},
#'  the winner genotype.
#'
#' @details \code{gafem} analyses data from a one-way genotype testing
#'   experiment. By default, a randomized complete block design is used
#'   according to the following model:
#' \deqn{Y_{ij} = m + g_i + r_j + e_{ij}}
#' where \eqn{Y_{ij}} is the response variable of the ith genotype in the
#' \emph{j}th block; \emph{m} is the grand mean (fixed); \eqn{g_i} is the effect
#' of the \emph{i}th genotype; \eqn{r_j} is the effect of the \emph{j}th
#' replicate; and \eqn{e_{ij}} is the random error.
#'
#' When \code{block} is informed, then a resolvable alpha design is implemented,
#' according to the following model:
#'
#' \deqn{Y_{ijk} = m + g_i + r_j + b_{jk} + e_{ijk}} where where \eqn{y_{ijk}}
#' is the response variable of the \emph{i}th genotype in the \emph{k}th block
#' of the \emph{j}th replicate; \emph{m} is the intercept, \eqn{t_i} is the
#' effect for the \emph{i}th genotype \eqn{r_j} is the effect of the \emph{j}th
#' replicate, \eqn{b_{jk}} is the effect of the \emph{k}th incomplete block of
#' the \emph{j}th replicate, and \eqn{e_{ijk}} is the plot error effect
#' corresponding to \eqn{y_{ijk}}. All effects, except the random error are
#' assumed to be fixed.
#'
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso \code{\link{get_model_data}} \code{\link{gamem}}
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
                  prob = 0.05,
                  block = NULL,
                  verbose = TRUE) {
  # RCBD
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
        means_by(GEN) %>%
        top_n(1, -Y) %>%
        select(GEN, Y) %>%
        slice(1)
      max_gen <-
        data %>%
        means_by(GEN) %>%
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
        cat("variable", paste(names(vars[var])),"\n")
        cat("---------------------------------------------------------------------------\n")
        cat("One-way ANOVA table (randomized complete block design)\n")
        cat("---------------------------------------------------------------------------\n")
        print(as.data.frame(anova), digits = 3, row.names = FALSE)
        cat("---------------------------------------------------------------------------\n\n")
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
        means_by(GEN) %>%
        top_n(1, -Y) %>%
        select(GEN, Y) %>%
        slice(1)
      max_gen <-
        data %>%
        means_by(GEN) %>%
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
        cat("variable", paste(names(vars[var])),"\n")
        cat("---------------------------------------------------------------------------\n")
        cat("One-way ANOVA table (alpha-lattice design)\n")
        cat("---------------------------------------------------------------------------\n")
        print(as.data.frame(anova), digits = 3, row.names = FALSE)
        cat("---------------------------------------------------------------------------\n\n")
      }
    }
  }
  if (verbose == TRUE) {
    if (length(which(unlist(lapply(listres, function(x) {
      x[["anova"]][1, 6]
    })) > prob)) > 0) {
      cat("---------------------------------------------------------------------------\n")
      cat("Variables with nonsignificant genotype effect\n")
      cat(names(which(unlist(lapply(listres, function(x) {
        as.numeric(x[["anova"]][2, 6])
      })) > prob)), "\n")
      cat("---------------------------------------------------------------------------\n")
    }
    cat("Done!\n")
  }
  invisible(structure(listres, class = "gafem"))
}

#' Several types of residual plots
#'
#' Residual plots for a output model of class \code{gafem}. Seven types
#' of plots are produced: (1) Residuals vs fitted, (2) normal Q-Q plot for the
#' residuals, (3) scale-location plot (standardized residuals vs Fitted Values),
#' (4) standardized residuals vs Factor-levels, (5) Histogram of raw residuals
#' and (6) standardized residuals vs observation order, and (7) 1:1 line plot.
#'
#'
#' @param x An object of class \code{gafem}.
#' @param ... Additional arguments passed on to the function
#'   \code{\link{residual_plots}}
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
