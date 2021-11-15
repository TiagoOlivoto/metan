#' Genotype-environment analysis by mixed-effect models
#' @description
#' `r badge('stable')`
#'
#' Genotype analysis in multi-environment trials using mixed-effect or
#' random-effect models.
#'
#'
#' The nature of the effects in the model is chosen with the argument
#' `random`. By default, the experimental design considered in each
#' environment is a randomized complete block design. If `block` is
#' informed, a resolvable alpha-lattice design (Patterson and Williams, 1976) is
#' implemented. The following six models can be fitted depending on the values
#' of `random` and `block` arguments.
#'   *  **Model 1:** `block = NULL` and `random = "gen"` (The
#'   default option). This model considers a Randomized Complete Block Design in
#'   each environment assuming genotype and genotype-environment interaction as
#'   random effects. Environments and blocks nested within environments are
#'   assumed to fixed factors.
#'
#'   *  **Model 2:** `block = NULL` and `random = "env"`. This
#'   model considers a Randomized Complete Block Design in each environment
#'   treating environment, genotype-environment interaction, and blocks nested
#'   within environments as random factors. Genotypes are assumed to be fixed
#'   factors.
#'
#'   *  **Model 3:** `block = NULL` and `random = "all"`. This
#'   model considers a Randomized Complete Block Design in each environment
#'   assuming a random-effect model, i.e., all effects (genotypes, environments,
#'   genotype-vs-environment interaction and blocks nested within environments)
#'   are assumed to be random factors.
#'
#'   *  **Model 4:** `block` is not `NULL` and `random =
#'   "gen"`. This model considers an alpha-lattice design in each environment
#'   assuming genotype, genotype-environment interaction, and incomplete blocks
#'   nested within complete replicates as random to make use of inter-block
#'   information (Mohring et al., 2015). Complete replicates nested within
#'   environments and environments are assumed to be fixed factors.
#'
#'   *  **Model 5:** `block` is not `NULL` and `random =
#'   "env"`. This model considers an alpha-lattice design in each environment
#'   assuming genotype as fixed. All other sources of variation (environment,
#'   genotype-environment interaction, complete replicates nested within
#'   environments, and incomplete blocks nested within replicates) are assumed
#'   to be random factors.
#'
#'   *  **Model 6:** `block` is not `NULL` and `random =
#'   "all"`. This model considers an alpha-lattice design in each environment
#'   assuming all effects, except the intercept, as random factors.
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#'   replications/blocks.
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure a vector of variables may be used. For example `resp
#'   = c(var1, var2, var3)`.
#' @param block Defaults to `NULL`. In this case, a randomized complete
#'   block design is considered. If block is informed, then an alpha-lattice
#'   design is employed considering block as random to make use of inter-block
#'   information, whereas the complete replicate effect is always taken as
#'   fixed, as no inter-replicate information was to be recovered (Mohring et
#'   al., 2015).
#'@param by One variable (factor) to compute the function by. It is a shortcut
#'  to [dplyr::group_by()].This is especially useful, for example,
#'  when the researcher want to analyze environments within mega-environments.
#'  In this case, an object of class waasb_grouped is returned.
#' @param random The effects of the model assumed to be random. Defaults to
#'   `random = "gen"`. See **Details** to see the random effects
#'   assumed depending on the experimental design of the trials.
#' @param prob The probability for estimating confidence interval for BLUP's
#'   prediction.
#' @param verbose Logical argument. If `verbose = FALSE` the code will run
#'   silently.
#' @references
#' Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, V.S. Marchioro, V.Q. de
#' Souza, and E. Jost. 2019. Mean performance and stability in multi-environment
#' trials I: Combining features of AMMI and BLUP techniques. Agron. J.
#' 111:2949-2960.
#' \doi{10.2134/agronj2019.03.0220}
#'
#' Mohring, J., E. Williams, and H.-P. Piepho. 2015. Inter-block information: to
#' recover or not to recover it? TAG. Theor. Appl. Genet. 128:1541-54.
#' \doi{10.1007/s00122-015-2530-0}
#'
#' Patterson, H.D., and E.R. Williams. 1976. A new class of resolvable
#' incomplete block designs. Biometrika 63:83-92.
#'
#'
#' @return An object of class `waasb` with the following items for each
#'   variable:
#'
#'
#' * **fixed** Test for fixed effects.
#'
#' * **random** Variance components for random effects.
#'
#' * **LRT** The Likelihood Ratio Test for the random effects.
#'
#'
#' * **BLUPgen** The random effects and estimated BLUPS for genotypes (If
#' `random = "gen"` or `random = "all"`)
#'
#' * **BLUPenv** The random effects and estimated BLUPS for environments,
#' (If `random = "env"` or `random = "all"`).
#'
#' * **BLUPint** The random effects and estimated BLUPS of all genotypes in
#' all environments.
#'
#'
#' * **MeansGxE** The phenotypic means of genotypes in the environments.
#'
#' * **modellme** The mixed-effect model of class `lmerMod`.
#'
#' * **residuals** The residuals of the mixed-effect model.
#'
#' * **model_lm** The fixed-effect model of class `lm`.
#'
#' * **residuals_lm** The residuals of the fixed-effect model.
#'
#' * **Details** A list summarizing the results. The following information
#' are shown: `Nenv`, the number of environments in the analysis;
#' `Ngen` the number of genotypes in the analysis; `Mean` the grand
#' mean; `SE` the standard error of the mean; `SD` the standard
#' deviation. `CV` the coefficient of variation of the phenotypic means,
#' estimating WAASB, `Min` the minimum value observed (returning the
#' genotype and environment), `Max` the maximum value observed (returning
#' the genotype and environment); `MinENV` the environment with the lower
#' mean, `MaxENV` the environment with the larger mean observed,
#' `MinGEN` the genotype with the lower mean, `MaxGEN` the genotype
#' with the larger.
#'
#' * **ESTIMATES** A tibble with the genetic parameters (if `random =
#' "gen"` or `random = "all"`) with the following columns: `Phenotypic
#' variance` the phenotypic variance; `Heritability` the broad-sense
#' heritability; `GEr2` the coefficient of determination of the interaction
#' effects; `h2mg` the heritability on the mean basis;
#' `Accuracy` the selective accuracy; `rge` the genotype-environment
#' correlation; `CVg` the genotypic coefficient of variation; `CVr`
#' the residual coefficient of variation; `CV ratio` the ratio between
#' genotypic and residual coefficient of variation.
#'
#'  * **formula** The formula used to fit the mixed-model.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso [mtsi()] [waas()]
#'   [get_model_data()] [plot_scores()]
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' #===============================================================#
#' # Example 1: Analyzing all numeric variables assuming genotypes #
#' # as random effects                                             #
#' #===============================================================#
#'model <- gamem_met(data_ge,
#'                   env = ENV,
#'                   gen = GEN,
#'                   rep = REP,
#'                   resp = everything())
#' # Distribution of random effects (first variable)
#' plot(model, type = "re")
#'
#' # Genetic parameters
#' get_model_data(model, "genpar")
#'
#'
#'
#' #===============================================================#
#' # Example 2: Unbalanced trials                                  #
#' # assuming all factors as random effects                        #
#' #===============================================================#
#' un_data <- data_ge %>%
#'              remove_rows(1:3) %>%
#'              droplevels()
#'
#'model2 <- gamem_met(un_data,
#'                    env = ENV,
#'                    gen = GEN,
#'                    rep = REP,
#'                    random = "all",
#'                    resp = GY)
#'get_model_data(model2)
#' }
#'
gamem_met <- function(.data,
                      env,
                      gen,
                      rep,
                      resp,
                      block = NULL,
                      by = NULL,
                      random = "gen",
                      prob = 0.05,
                      verbose = TRUE) {
  if (!random %in% c("env", "gen", "all")) {
    stop("The argument 'random' must be one of the 'gen', 'env', or 'all'.")
  }
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
        doo(gamem_met,
            env = {{env}},
            gen = {{gen}},
            rep = {{rep}},
            resp = {{resp}},
            block = {{block}},
            random = random,
            prob = prob,
            verbose = verbose)
    } else{
      results <-
        .data %>%
        doo(gamem_met,
            env = {{env}},
            gen = {{gen}},
            rep = {{rep}},
            resp = {{resp}},
            random = random,
            prob = prob,
            verbose = verbose)
    }
    return(set_class(results, c("tbl_df",  "waasb_group", "tbl",  "data.frame")))
  }
  block_test <- missing(block)
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
  model_fixed <-
    case_when(
      block_test ~ paste("Y ~ GEN + ENV/REP + GEN:ENV"),
     !block_test ~ paste("Y ~ GEN + ENV/REP/BLOCK + GEN:ENV")
    )
  model_formula <-
    case_when(
      random == "gen" & block_test ~ paste("Y ~ ENV/REP + (1 | GEN) + (1 | GEN:ENV)"),
      random == "env" & block_test ~ paste("Y ~ GEN + (1 | ENV/REP) + (1 | GEN:ENV)"),
      random == "all" & block_test ~ paste("Y ~ (1 | GEN) + (1 | ENV/REP) + (1 | GEN:ENV)"),
      random == "gen" & !block_test ~ paste("Y ~  (1 | GEN) + ENV / REP + (1|BLOCK:(REP:ENV))  + (1 | GEN:ENV)"),
      random == "env" & !block_test ~ paste("Y ~ 0 + GEN + (1| ENV/REP/BLOCK)  + (1 | GEN:ENV)"),
      random == "all" & !block_test ~ paste("Y ~  (1 | GEN) + (1|ENV/REP/BLOCK) + (1 | GEN:ENV)")
    )
  lrt_groups <-
    strsplit(
      case_when(
        random == "gen" & block_test ~ c("COMPLETE GEN GEN:ENV"),
        random == "env" & block_test ~ c("COMPLETE REP(ENV) ENV GEN:ENV"),
        random == "all" & block_test ~ c("COMPLETE GEN REP(ENV) ENV GEN:ENV"),
        random == "gen" & !block_test ~ c("COMPLETE GEN BLOCK(ENV:REP) GEN:ENV"),
        random == "env" & !block_test ~ c("COMPLETE BLOCK(ENV:REP) REP(ENV) ENV GEN:ENV"),
        random == "all" & !block_test ~ c("COMPLETE GEN BLOCK(ENV:REP) REP(ENV) ENV GEN:ENV")
      ), " ")[[1]]
  mod1 <- random == "gen" & block_test
  mod2 <- random == "gen" & !block_test
  mod3 <- random == "env" & block_test
  mod4 <- random == "env" & !block_test
  mod5 <- random == "all" & block_test
  mod6 <- random == "all" & !block_test
  nvar <- ncol(vars)
  listres <- list()
  vin <- 0
  if (verbose == TRUE) {
    pb <- progress(max = nvar, style = 4)
  }
  for (var in 1:nvar) {
    data <- factors %>%
      mutate(Y = vars[[var]])
    check_labels(data)
    if(has_na(data)){
      data <- remove_rows_na(data)
      has_text_in_num(data)
    }
    if(!is_balanced_trial(data, ENV, GEN, Y) && random == "env"){
      warning("Fitting a model with unbalanced data considering genotype as fixed effect is not suggested.", call. = FALSE)
    }
    Nenv <- nlevels(data$ENV)
    Ngen <- nlevels(data$GEN)
    Nrep <- nlevels(data$REP)
    vin <- vin + 1
    ovmean <- mean(data$Y)
    fixed_mod <- lm(model_fixed, data = data)
    Complete <-
      lmerTest::lmer(model_formula, data = data) %>%
      suppressWarnings() %>%
      suppressMessages()
    LRT <-
      lmerTest::ranova(Complete, reduce.terms = FALSE) %>%
      mutate(model = lrt_groups) %>%
      column_to_first(model) %>%
      suppressWarnings() %>%
      suppressMessages()
    fixed <- anova(Complete)
    var_eff <-
      lme4::VarCorr(Complete) %>%
      as.data.frame() %>%
      select_cols(1, 4) %>%
      arrange(grp) %>%
      rename(Group = grp, Variance = vcov) %>%
      add_cols(Percent = (Variance / sum(Variance)) * 100)
    if(random %in% c("gen", "all")){
      GV <- as.numeric(var_eff[which(var_eff[1] == "GEN"), 2])
      IV <- as.numeric(var_eff[which(var_eff[1] == "GEN:ENV"), 2])
      RV <- as.numeric(var_eff[which(var_eff[1] == "Residual"), 2])
      FV <- sum(var_eff$Variance)
      h2g <- GV/FV
      h2mg <- GV/(GV + IV/Nenv + RV/(Nenv * Nrep))
      GEr2 <- IV/(GV + IV + RV)
      AccuGen <- sqrt(h2mg)
      rge <- IV/(IV + RV)
      CVg <- (sqrt(GV)/ovmean) * 100
      CVr <- (sqrt(RV)/ovmean) * 100
      CVratio <- CVg/CVr
      PROB <- ((1 - (1 - prob))/2) + (1 - prob)
      t <- qt(PROB, 100)
      Limits <- t * sqrt(((1 - AccuGen) * GV))
      genpar <- tibble(Parameters = c("Phenotypic variance", "Heritability", "GEIr2", "h2mg",
                                      "Accuracy", "rge", "CVg", "CVr", "CV ratio"),
                       Values = c(FV, h2g, GEr2, h2mg, AccuGen, rge, CVg, CVr, CVratio))
    } else{
      genpar <- NULL
    }
    bups <- lme4::ranef(Complete)
    bINT <-
      data.frame(Names = rownames(bups[["GEN:ENV"]])) %>%
      separate(Names, into = c("GEN", "ENV"), sep = ":") %>%
      add_cols(BLUPge = bups[["GEN:ENV"]][[1]]) %>%
      as_factor(1:2)
    Details <-
      rbind(ge_details(data, ENV, GEN, Y),
            tribble(~Parameters,  ~Y,
                    "Ngen", Ngen,
                    "Nenv", Nenv)) %>%
      rename(Values = Y)
    if(mod1){
      ran_ef <- c("GEN, GEN:ENV")
      fix_ef <- c("ENV, REP(ENV)")
      rand_ef <-
      data_factors <- data %>% select_non_numeric_cols()
      BLUPgen <-
        means_by(data, GEN) %>%
        add_cols(BLUPg = bups$GEN$`(Intercept)`,
                 Predicted = BLUPg + ovmean,
                 Rank = rank(-Predicted),
                 LL = Predicted - Limits,
                 UL = Predicted + Limits) %>%
        arrange(-Predicted) %>%
        column_to_first(Rank)
      BLUPint <-
        suppressWarnings(
          left_join(data_factors, bINT, by = c("ENV", "GEN")) %>%
            left_join(BLUPgen, by = "GEN") %>%
            select(ENV, GEN, REP, BLUPg, BLUPge) %>%
            add_cols(`BLUPg+ge` = BLUPge + BLUPg,
                     Predicted = predict(Complete))
        )
      BLUPenv <- NULL
    } else if(mod2){
      ran_ef <- c("GEN, BLOCK(ENV:REP), GEN:ENV")
      fix_ef <- c("ENV, REP(ENV)")
      data_factors <- data %>% select_non_numeric_cols()
      BLUPgen <-
        means_by(data, GEN) %>%
        add_cols(BLUPg = bups$GEN$`(Intercept)`,
                 Predicted = BLUPg + ovmean,
                 Rank = rank(-Predicted),
                 LL = Predicted - Limits,
                 UL = Predicted + Limits) %>%
        arrange(-Predicted) %>%
        column_to_first(Rank)
      blupBRE <-
        data.frame(Names = rownames(bups$`BLOCK:(REP:ENV)`)) %>%
        separate(Names, into = c("BLOCK", "REP", "ENV"), sep = ":") %>%
        add_cols(BLUPbre = bups$`BLOCK:(REP:ENV)`[[1]]) %>%
        as_factor(1:3)
      BLUPint <-
        suppressWarnings(
          left_join(data_factors, bINT, by = c("ENV", "GEN")) %>%
            left_join(BLUPgen, by = "GEN") %>%
            left_join(blupBRE, by = c("ENV", "REP", "BLOCK")) %>%
            select(ENV, REP, BLOCK, GEN, BLUPg, BLUPge, BLUPbre) %>%
            add_cols(`BLUPg+ge+bre` = BLUPge + BLUPg + BLUPbre,
                     Predicted = `BLUPg+ge+bre` + left_join(data_factors, data %>% means_by(ENV, REP), by = c("ENV", "REP"))$Y)
        )
      BLUPenv <- NULL
    } else if (mod3){
      ran_ef <- c("REP(ENV), ENV, GEN:ENV")
      fix_ef <- c("GEN")
      data_factors <- data %>% select_non_numeric_cols()
      BLUPgen <- NULL
      BLUPenv <-
        means_by(data, ENV) %>%
        add_cols(BLUPe =  bups$ENV$`(Intercept)`,
                 Predicted = BLUPe + ovmean) %>%
        arrange(-Predicted) %>%
        add_cols(Rank = rank(-Predicted)) %>%
        column_to_first(Rank)
      blupRWE <-
        data.frame(Names = rownames(bups$`REP:ENV`)) %>%
        separate(Names, into = c("REP", "ENV"), sep = ":") %>%
        add_cols(BLUPre = bups$`REP:ENV`[[1]]) %>%
        as_factor(1:2)
      BLUPint <-
        suppressWarnings(
          left_join(data_factors, bINT, by = c("ENV", "GEN")) %>%
            left_join(BLUPenv, by = "ENV") %>%
            left_join(blupRWE, by = c("ENV", "REP")) %>%
            select(ENV, GEN, REP, BLUPe, BLUPge, BLUPre) %>%
            add_cols(`BLUPge+e+re` = BLUPge + BLUPe + BLUPre,
                     Predicted = `BLUPge+e+re` + left_join(data_factors, means_by(data, GEN), by = c("GEN"))$Y)
        )
    } else if (mod4){
      ran_ef <- c("BLOCK(ENV:REP), REP(ENV), ENV, GEN:ENV")
      fix_ef <- c("GEN")
      data_factors <- data %>% select_non_numeric_cols()
      BLUPgen <- NULL
      BLUPenv <-
        means_by(data, ENV) %>%
        add_cols(BLUPe =  bups$ENV$`(Intercept)`,
                 Predicted = BLUPe + ovmean) %>%
        arrange(-Predicted) %>%
        add_cols(Rank = rank(-Predicted)) %>%
        column_to_first(Rank)
      blupRWE <-
        data.frame(Names = rownames(bups$`REP:ENV`)) %>%
        separate(Names, into = c("REP", "ENV"), sep = ":") %>%
        add_cols(BLUPre = bups$`REP:ENV`[[1]]) %>%
        as_factor(1:2)
      blupBRE <-
        data.frame(Names = rownames(bups$`BLOCK:(REP:ENV)`)) %>%
        separate(Names, into = c("BLOCK", "REP", "ENV"), sep = ":") %>%
        add_cols(BLUPbre = bups$`BLOCK:(REP:ENV)`[[1]]) %>%
        as_factor(1:3)
      genCOEF <-
        summary(Complete)[["coefficients"]] %>%
        as.data.frame() %>%
        rownames_to_column("GEN") %>%
        replace_string(GEN, pattern = "GEN") %>%
        rename(Y = Estimate) %>%
        as_factor(1)
      BLUPint <-
        suppressWarnings(
          left_join(data_factors, bINT, by = c("ENV", "GEN")) %>%
            left_join(BLUPenv, by = "ENV") %>%
            left_join(blupRWE, by = c("ENV", "REP")) %>%
            left_join(blupBRE, by = c("ENV", "REP", "BLOCK")) %>%
            select(ENV, REP, BLOCK, GEN, BLUPe, BLUPge, BLUPre, BLUPbre) %>%
            add_cols(`BLUPe+ge+re+bre` = BLUPge + BLUPe + BLUPre + BLUPbre,
                     Predicted = `BLUPe+ge+re+bre` + left_join(data_factors, genCOEF, by = "GEN")$Y)
        )
    } else if (mod5){
      ran_ef <- c("GEN, REP(ENV), ENV, GEN:ENV")
      fix_ef <- c("-")
      data_factors <- data %>% select_non_numeric_cols()
      BLUPgen <-
        means_by(data, GEN) %>%
        add_cols(BLUPg = bups$GEN$`(Intercept)`,
                 Predicted = BLUPg + ovmean,
                 Rank = rank(-Predicted),
                 LL = Predicted - Limits,
                 UL = Predicted + Limits) %>%
        arrange(-Predicted) %>%
        column_to_first(Rank)
      BLUPenv <-
        means_by(data, ENV) %>%
        add_cols(BLUPe =  bups$ENV$`(Intercept)`,
                 Predicted = BLUPe + ovmean,
                 Rank = rank(-Predicted)) %>%
        arrange(-Predicted) %>%
        column_to_first(Rank)
      blupRWE <- data.frame(Names = rownames(bups$`REP:ENV`)) %>%
        separate(Names, into = c("REP", "ENV"), sep = ":") %>%
        add_cols(BLUPre = bups$`REP:ENV`[[1]]) %>%
        arrange(ENV) %>%
        as_factor(1:2)
      BLUPint <-
        suppressWarnings(
          left_join(data_factors, bINT, by = c("ENV", "GEN")) %>%
            left_join(BLUPgen, by = "GEN") %>%
            left_join(BLUPenv, by = "ENV") %>%
            left_join(blupRWE, by = c("ENV", "REP")) %>%
            select(GEN, ENV, REP, BLUPe, BLUPg, BLUPge, BLUPre) %>%
            add_cols(`BLUPg+e+ge+re` = BLUPge + BLUPe + BLUPg + BLUPre,
                     Predicted = `BLUPg+e+ge+re` + ovmean)
        )
    } else if (mod6){
      ran_ef <- c("GEN, BLOCK(ENV:REP), REP(ENV), ENV, GEN:ENV")
      fix_ef <- c("-")
      data_factors <- data %>% select_non_numeric_cols()
      BLUPgen <-
        means_by(data, GEN) %>%
        add_cols(BLUPg = bups$GEN$`(Intercept)`,
                 Predicted = BLUPg + ovmean,
                 Rank = rank(-Predicted),
                 LL = Predicted - Limits,
                 UL = Predicted + Limits) %>%
        arrange(-Predicted) %>%
        column_to_first(Rank)
      BLUPenv <-
        means_by(data, ENV) %>%
        add_cols(BLUPe =  bups$ENV$`(Intercept)`,
                 Predicted = BLUPe + ovmean,
                 Rank = rank(-Predicted)) %>%
        arrange(-Predicted) %>%
        column_to_first(Rank)
      blupRWE <- data.frame(Names = rownames(bups$`REP:ENV`)) %>%
        separate(Names, into = c("REP", "ENV"), sep = ":") %>%
        add_cols(BLUPre = bups$`REP:ENV`[[1]]) %>%
        arrange(ENV) %>%
        as_factor(1:2)
      blupBRE <-
        data.frame(Names = rownames(bups$`BLOCK:(REP:ENV)`)) %>%
        separate(Names, into = c("BLOCK", "REP", "ENV"), sep = ":") %>%
        add_cols(BLUPbre = bups$`BLOCK:(REP:ENV)`[[1]]) %>%
        as_factor(1:3)
      BLUPint <-
        suppressWarnings(
          left_join(data_factors, bINT, by = c("ENV", "GEN")) %>%
            left_join(BLUPgen, by = "GEN") %>%
            left_join(BLUPenv, by = "ENV") %>%
            left_join(blupRWE, by = c("ENV", "REP")) %>%
            left_join(blupBRE, by = c("ENV", "REP", "BLOCK")) %>%
            select(GEN, ENV, REP, BLOCK, BLUPg, BLUPe, BLUPge, BLUPre, BLUPbre) %>%
            add_cols(`BLUPg+e+ge+re+bre` = BLUPg + BLUPge + BLUPe + BLUPre + BLUPbre,
                     Predicted = `BLUPg+e+ge+re+bre` + ovmean)
        )
    }
    residuals <- data.frame(fortify.merMod(Complete))
    residuals$reff <- BLUPint$BLUPge
    temp <- structure(list(fixed = fixed %>% rownames_to_column("SOURCE") %>% as_tibble(),
                           random = var_eff,
                           LRT = LRT,
                           BLUPgen = BLUPgen,
                           BLUPenv = BLUPenv,
                           BLUPint = BLUPint,
                           MeansGxE = means_by(data, ENV, GEN),
                           modellme = Complete,
                           residuals = as_tibble(residuals),
                           model_lm = fixed_mod,
                           residuals_lm = tibble(fortify(fixed_mod)),
                           Details = as_tibble(Details),
                           ESTIMATES = genpar,
                           formula = model_formula), class = "waasb")
    if (verbose == TRUE) {
      run_progress(pb,
                   actual = var,
                   text = paste("Evaluating trait", names(vars[var])))
    }
    listres[[paste(names(vars[var]))]] <- temp
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
      x[["LRT"]] %>% dplyr::filter(model == "GEN:ENV") %>% pull(`Pr(>Chisq)`)
    })) > prob)) > 0) {
      cat("Variables with nonsignificant GxE interaction\n")
      cat(names(which(unlist(lapply(listres, function(x) {
        x[["LRT"]][which(x[["LRT"]][[1]] == "GEN:ENV"), 7]
      })) > prob)), "\n")
      cat("---------------------------------------------------------------------------\n")
    } else {
      cat("All variables with significant (p < 0.05) genotype-vs-environment interaction\n")
    }
  }
  invisible(set_class(listres, "waasb"))
}
