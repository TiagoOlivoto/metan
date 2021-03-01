#' Joint analysis of variance
#' @description
#' `r badge('stable')`
#'
#' Performs a joint analysis of variance to check for the presence of
#' genotype-vs-environment interactions using both randomized complete block and
#' alpha-lattice designs.
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments. The analysis of variance is computed for each level of this
#'   factor.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#'   replications/blocks.
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure a vector of variables may be used. For example `resp
#'   = c(var1, var2, var3)`.
#' @param block Defaults to `NULL`. In this case, a randomized complete
#'   block design is considered. If block is informed, then a resolvable
#'   alpha-lattice design (Patterson and Williams, 1976) is employed.
#'   **All effects, except the error, are assumed to be fixed.**
#' @param verbose Logical argument. If `verbose = FALSE` the code will run
#'   silently.
#' @references Patterson, H.D., and E.R. Williams. 1976. A new class of
#' resolvable incomplete block designs. Biometrika 63:83-92.
#' @return A list where each element is the result for one variable containing
#'   the following objects:
#'
#'  * **anova:** The two-way ANOVA table
#'  * **model:** The model of class `lm`.
#'  * **augment:** Information about each observation in the dataset. This
#'  includes predicted values in the `fitted` column, residuals in the
#'  `resid` column, standardized residuals in the `stdres` column,
#'  the diagonal of the 'hat' matrix in the `hat`, and standard errors for
#'  the fitted values in the `se.fit` column.
#'  * **details:** A tibble with the following data: `Ngen`, the
#'  number of genotypes; `OVmean`, the grand mean; `Min`, the minimum
#'  observed (returning the genotype and replication/block); `Max` the
#'  maximum observed, `MinGEN` the loser winner genotype, `MaxGEN`,
#'  the winner genotype.
#'
#' @md
#' @seealso [get_model_data()] [anova_ind()]
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' # traditional usage approach
#' j_an <- anova_joint(data_ge,
#'                     env = ENV,
#'                     gen = GEN,
#'                     rep = REP,
#'                     resp = everything())
#' # Predicted values
#' get_model_data(j_an)
#'
#' # Details
#' get_model_data(j_an, "details")
#'}
#'
anova_joint <- function(.data,
                        env,
                        gen,
                        rep,
                        resp,
                        block = NULL,
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
  listres <- list()
  nvar <- ncol(vars)
  for (var in 1:nvar) {
    data <- factors %>%
      mutate(mean = vars[[var]])
    if(has_na(data)){
      data <- remove_rows_na(data)
      has_text_in_num(data)
    }
    msr <- data %>%
      split_factors(ENV, keep_factors = T)
    if(missing(block)){
      msr <- do.call(rbind,
                     lapply(msr, function(x){
                       summary(aov(mean ~ GEN + REP, data = x))[[1]][3,3]
                     }))
      model <- aov(mean ~ GEN + ENV + GEN:ENV + ENV/REP, data = data)
      anova <-
        anova(model) %>%
        rownames_to_column("Source") %>%
        select_rows(2, 4, 1, 3, 5) %>%
        as.data.frame()
      anova[2, 1] <- "REP(ENV)"
      anova[1, 5] <- anova[1, 4] / anova[2, 4]
      anova[1, 6] <- 1 - pf(anova[1, 5], anova[1, 2], anova[2, 2])
      CV <- tibble(Source = "CV(%)", Df = as.numeric(sqrt(anova[5, 4]) / mean(data$mean) * 100))
      msr <- tibble(Source = "MSR+/MSR-", Df = max(msr) / min(msr))
      ovmean <- tibble(Source = "OVmean", Df = mean(data$mean))
      temp <- rbind_fill_id(anova, CV, msr, ovmean)
    } else{
      msr <- do.call(rbind,
                     lapply(msr, function(x){
                       summary(aov(mean ~ GEN + REP + REP:BLOCK, data = x))[[1]][4,3]
                     }))
      model <- aov(mean ~ GEN + ENV + GEN:ENV + ENV/REP/BLOCK, data = data)
      anova <- anova(model) %>%
        as.data.frame() %>%
        rownames_to_column("Source") %>%
        select_rows(2, 4, 5, 1, 3, 6)
      anova[2, 1] <- "REP(ENV)"
      anova[3, 1] <- "BLOCK(REP*ENV)"
      CV <- tibble(Source = "CV(%)", Df = as.numeric(sqrt(anova[6, 4]) / mean(data$mean) * 100))
      msr <- tibble(Source = "MSR+/MSR-", Df = max(msr) / min(msr))
      ovmean <- tibble(Source = "OVmean", Df = mean(data$mean))
      temp <- rbind_fill_id(anova, CV, msr, ovmean)
    }
    influence <- lm.influence(model)
    augment <- model$model %>%
      add_cols(hat = influence$hat,
               sigma = influence$sigma,
               fitted = predict(model),
               resid = residuals(model),
               stdres = rstandard(model),
               se.fit = predict(model, se.fit = TRUE)$se.fit) %>%
      add_cols(factors = concatenate(., GEN, REP, pull = TRUE))
    if(missing(block)){
      augment <- column_to_first(augment, ENV, GEN, REP)
    } else{
      augment <- column_to_first(augment, ENV, GEN, REP, BLOCK)
    }
    listres[[paste(names(vars[var]))]] <- list(anova = temp,
                                               model = model,
                                               augment = augment,
                                               details = ge_details(data, ENV, GEN, mean))
    if (verbose == TRUE) {
      cat("variable", paste(names(vars[var])),"\n")
      cat("---------------------------------------------------------------------------\n")
      cat("Joint ANOVA table\n")
      cat("---------------------------------------------------------------------------\n")
      print(as.data.frame(temp), digits = 3, row.names = FALSE)
      cat("---------------------------------------------------------------------------\n\n")
    }
  }
  if (verbose == TRUE) {
    index <- ifelse(missing(block), 4, 5)
    if (length(which(unlist(lapply(listres, function(x) {
      x[["anova"]][index, 6]
    })) > 0.05)) > 0) {
      cat("---------------------------------------------------------------------------\n")
      cat("Variables with nonsignificant GxE interaction\n")
      cat(names(which(unlist(lapply(listres, function(x) {
        x[["anova"]][["Pr(>F)"]][index]
      })) > 0.05)), "\n")
      cat("---------------------------------------------------------------------------\n")
    } else {
      cat("All variables with significant (p < 0.05) genotype-vs-environment interaction\n")
    }
    cat("Done!\n")
  }
  return(structure(listres, class = "anova_joint"))
}



#' Several types of residual plots
#'
#' Residual plots for a output model of class `anova_joint`. Seven types
#' of plots are produced: (1) Residuals vs fitted, (2) normal Q-Q plot for the
#' residuals, (3) scale-location plot (standardized residuals vs Fitted Values),
#' (4) standardized residuals vs Factor-levels, (5) Histogram of raw residuals
#' and (6) standardized residuals vs observation order, and (7) 1:1 line plot.
#'
#'
#' @param x An object of class `anova_joint`.
#' @param ... Additional arguments passed on to the function
#'   [residual_plots()]
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method plot anova_joint
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' model <- anova_joint(data_ge, ENV, GEN, REP, GY)
#' plot(model)
#' plot(model,
#'      which = c(3, 5),
#'      nrow = 2,
#'      labels = TRUE,
#'      size.lab.out = 4)
#' }
#'
plot.anova_joint <- function(x, ...) {
  residual_plots(x,  ...)
}






#' Print an object of class anova_joint
#'
#' Print the `anova_joint` object in two ways. By default, the results are shown
#' in the R console. The results can also be exported to the directory into a
#' *.txt file.
#'
#' @param x An object of class `anova_joint`.
#' @param export A logical argument. If `TRUE`, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if `export = TRUE`
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   [`tibble::print()`][tibble::formatting] for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print anova_joint
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' model <- data_ge %>% anova_joint(ENV, GEN, REP, c(GY, HM))
#' print(model)
#' }
print.anova_joint <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (!class(x) == "anova_joint") {
    stop("The object must be of class 'anova_joint'")
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "anova_joint print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    print(var)
    cat("---------------------------------------------------------------------------\n")
    cat("\n\n\n")
  }
  if (export == TRUE) {
    sink()
  }
}
