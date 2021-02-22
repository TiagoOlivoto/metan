#' Get data from a model easily
#'
#' `r badge('experimental')`
#'
#' * `get_model_data()` Easily get data from some objects generated in the
#' **metan** package such as the WAASB and WAASBY indexes  (Olivoto et al.,
#' 2019a, 2019b) BLUPs, variance components, details of AMMI models and
#' AMMI-based stability statistics.
#' * `gmd()` Is a shortcut to `get_model_data`.
#' * `sel_gen()` Extracts the selected genotypes by a given index.
#' @param x An object created with the functions [AMMI_indexes()],
#'   [anova_ind()], [anova_joint()], [can_corr()] [ecovalence()],  [Fox()],
#'   [gai()], [gamem()],[gafem()], [ge_acv()], [ge_means()], [ge_reg()],
#'   [gytb()], [mgidi()], [performs_ammi()], [Resende_indexes()], [Shukla()],
#'   [superiority()], [waas()] or [waasb()].
#' @param what What should be captured from the model. See more in section
#'   **Details**.
#' @param type Chose if the statistics must be show by genotype (`type =
#'   "GEN"`, default) or environment (`TYPE = "ENV"`), when possible.
#' @param verbose Logical argument. If `verbose = FALSE` the code will run
#'   silently.
#' @return A tibble showing the values of the variable chosen in argument
#'   `what`.
#' @name get_model_data
#' @details
#' Bellow are listed the options allowed in the argument `what` depending
#' on the class of the object
#'
#'  **Objects of class `AMMI_indexes`:**
#' * `"ASV"` AMMI stability value.
#' * `"EV"` Averages of the squared eigenvector values.
#' * `"SIPC"` Sums of the absolute value of the IPCA scores.
#' * `"WAAS"` Weighted average of absolute scores (default).
#' * `"ZA"` Absolute value of the relative contribution of IPCAs to the
#' interaction.
#'
#'  **Objects of class `anova_ind`:**
#' * `"MEAN"`The mean value of the variable
#' * `"MSG", "FCG", "PFG"` The mean square, F-calculated and P-values for
#' genotype effect, respectively.
#' * `"MSB", "FCB", "PFB"` The mean square, F-calculated and P-values for
#' block effect in randomized complete block design.
#' * `"MSCR", "FCR", "PFCR"` The mean square, F-calculated and P-values for
#' complete replicates in alpha lattice design.
#' * `"MSIB_R", "FCIB_R", "PFIB_R"` The mean square, F-calculated and
#' P-values for incomplete blocks within complete replicates, respectively (for
#' alpha lattice design only).
#' * `"MSE"` The mean square of error.
#' * `"CV"` The coefficient of variation.
#' * `"h2"` The broad-sence heritability
#' * `"MSE"` The accucary of selection (square root of h2).
#'
#'
#'  **Objects of class `anova_joint` or `gafem`:**
#' * `"Y"` The observed values.
#' * `"h2"` The broad-sense heritability.
#' * `"Sum Sq"` Sum of squares.
#' * `"Mean Sq"` Mean Squares.
#' * `"F value"` F-values.
#' * `"Pr(>F)"` P-values.
#' * `".fitted"` Fitted values (default).
#' * `".resid"` Residuals.
#' * `".stdresid"` Standardized residuals.
#' * `".se.fit"` Standard errors of the fitted values.
#' * `"details"` Details.
#'
#'  **Objects of class `Annicchiarico` and `Schmildt`:**
#' * `"Sem_rp"` The standard error of the relative mean performance (Schmildt).
#' * `"Mean_rp"` The relative performance of the mean.
#' * `"rank"` The rank for genotypic confidence index.
#' * `"Wi"` The genotypic confidence index.
#'
#'  **Objects of class `can_corr`:**
#' * `"coefs"` The canonical coefficients (default).
#' * `"loads"` The canonical loadings.
#' * `"crossloads"` The canonical cross-loadings.
#' * `"canonical"` The canonical correlations and hypothesis testing.
#'
#'  **Objects of class `ecovalence`:**
#' * `"Ecoval"` Ecovalence value (default).
#' * `"Ecov_perc"` Ecovalence in percentage value.
#' * `"rank"` Rank for ecovalence.
#'
#'  **Objects of class `fai_blup`:** See the **Value** section of
#'  [fai_blup()] to see valid options for `what` argument.
#'
#'  **Objects of class `ge_acv`:**
#' * `"ACV"` The adjusted coefficient of variation (default).
#' * `"ACV_R"` The rank for adjusted coefficient of variation.
#'
#'  **Objects of class `ge_polar`:**
#' * `"POLAR"` The Power Law Residuals (default).
#' * `"POLAR_R"` The rank for Power Law Residuals.
#'
#'  **Objects of class `ge_reg`:**
#' * `GEN`: the genotypes.
#' * `b0` and `b1` (default): the intercept and slope of the regression,
#' respectively.
#' * `t(b1=1)`: the calculated t-value
#' * `pval_t`: the p-value for the t test.
#' *  `s2di` the deviations from the regression (stability parameter).
#' * `F(s2di=0)`: the F-test for the deviations.
#' *  `pval_f`: the p-value for the F test;
#' *  `RMSE` the root-mean-square error.
#' *  `R2` the determination coefficient of the regression.
#'
#'
#'  **Objects of class `ge_effects`:**
#' * For objects of class `ge_effects` no argument `what` is required.
#'
#'  **Objects of class `ge_means`:**
#' * `"ge_means"` Genotype-environment interaction means (default).
#' * `"env_means"` Environment means.
#' * `"gen_means"` Genotype means.
#'
#'  **Objects of class `gge`:**
#' * `"scores"` The scores for genotypes and environments for all the
#' analyzed traits (default).
#' * `"exp_var"` The eigenvalues and explained variance.
#'
#'  **Objects of class `gytb`:**
#' * `"gyt"` Genotype by yield*trait table (Default).
#' * `"stand_gyt"` The standardized (zero mean and unit variance) Genotype by yield*trait table.
#' * `"si"` The superiority index (sum standardized value across all yield*trait combinations).
#'
#'  **Objects of class `mgidi`:** See the **Value** section of
#'  [mgidi()] to see valid options for `what` argument.
#'
#'  **Objects of class `mtsi`:** See the **Value** section of
#'  [mtsi()] to see valid options for `what` argument.
#'
#'  **Objects of class `Shukla`:**
#' * `"rMean"` Rank for the mean.
#' * `"ShuklaVar"` Shukla's stablity variance (default).
#' * `"rShukaVar"` Rank for Shukla's stablity variance.
#' * `"ssiShukaVar"` Simultaneous selection index.
#'
#'  **Objects of class `sh`:** See the **Value** section of
#'  [Smith_Hazel()] to see valid options for `what` argument.
#'
#'  **Objects of class `Fox`:**
#' * `"TOP"` The proportion of locations at which the genotype occurred in
#' the top third (default).
#'
#'  **Objects of class `gai`:**
#' * `"GAI"` The geometric adaptability index (default).
#' * `"GAI_R"` The rank for the GAI values.
#'
#'  **Objects of class `superiority`:**
#' * `"Pi_a"` The superiority measure for all environments (default).
#' * `"R_a"` The rank for Pi_a.
#' * `"Pi_f"` The superiority measure for favorable environments.
#' * `"R_f"` The rank for Pi_f.
#' * `"Pi_u"` The superiority measure for unfavorable environments.
#' * `"R_u"` The rank for Pi_u.
#'
#'  **Objects of class `Huehn`:**
#' * `"S1"` Mean of the absolute rank differences of a genotype over the n
#' environments (default).
#' * `"S2"` variance among the ranks over the k environments.
#' * `"S3"` Sum of the absolute deviations.
#' * `"S6"` Relative sum of squares of rank for each genotype.
#' * `"S1_R"`, `"S2_R"`, `"S3_R"`, and  `"S6_R"`, the ranks
#' for S1, S2, S3, and S6, respectively.
#'
#'  **Objects of class `Thennarasu`:**
#' * `"N1"` First statistic (default).
#' * `"N2"` Second statistic.
#' * `"N3"` Third statistic.
#' * `"N4"` Fourth statistic.
#' * `"N1_R"`, `"N2_R"`, `"N3_R"`, and `"N4_R"`, The ranks
#' for the statistics.
#'
#'
#'  **Objects of class `performs_ammi`:**
#' * `"PC1", "PC2", ..., "PCn"` The values for the nth interaction
#' principal component axis.
#' * `"ipca_ss"` Sum of square for each IPCA.
#' * `"ipca_ms"` Mean square for each IPCA.
#' * `"ipca_fval"` F value for each IPCA.
#' * `"ipca_pval"` P-value for for each IPCA.
#' * `"ipca_expl"`  Explained sum of square for each IPCA (default).
#' * `"ipca_accum"` Accumulated explained sum of square.
#'
#'
#' **Objects of class `waas`, `waas_means`, and `waasb`:**
#' * `"PC1", "PC2", ..., "PCn"` The values for the nth interaction
#' principal component axis.
#' * `"WAASB"`  The weighted average of the absolute scores (default for
#' objects of class `waas`).
#' * `"PctResp"` The rescaled values of the response variable.
#' * `"PctWAASB"` The rescaled values of the WAASB.
#' * `"wResp"` The weight for the response variable.
#' * `"wWAASB"` The weight for the stability.
#' * `"OrResp"` The ranking regarding the response variable.
#' * `"OrWAASB"` The ranking regarding the WAASB.
#' * `"OrPC1"` The ranking regarding the first principal component axix.
#' * `"WAASBY"` The superiority index WAASBY.
#' * `"OrWAASBY"` The ranking regarding the superiority index.
#'
#'  **Objects of class `gamem` and `waasb`:**
#' * `"blupge"` for genotype-vs-environment's predicted mean (class waasb).
#' * `"blupg"` For genotype's predicted mean.
#' * `"data"` The data used.
#' * `"details"` The details of the trial.
#' * `"genpar"` Genetic parameters (default).
#' * `"gcov"` The genotypic variance-covariance matrix.
#' * `"pcov"` The phenotypic variance-covariance matrix.
#' * `"gcor"` The genotypic correlation matrix.
#' * `"pcor"` The phenotypic correlation matrix.
#' * `"h2"` The broad-sense heritability.
#' * `"lrt"` The likelihood-ratio test for random effects.
#' * `"vcomp"` The variance components for random effects.
#' * `"ranef"` Random effects.
#'
#'  **Objects of class `Res_ind`**
#' * `"HMGV"` For harmonic mean of genotypic values.
#' * `"RPGV or RPGV_Y"` For relative performance of genotypic values
#' * `"HMRPGV"` For harmonic mean of relative performance of genotypic values
#'
#'
#' @md
#' @importFrom dplyr starts_with matches case_when full_join arrange_if
#' @importFrom purrr reduce
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @references
#'
#' Annicchiarico, P. 1992. Cultivar adaptation and recommendation from alfalfa
#' trials in Northern Italy. J. Genet. Breed. 46:269-278.
#'
#' Dias, P.C., A. Xavier, M.D.V. de Resende, M.H.P. Barbosa, F.A. Biernaski,
#' R.A. Estopa. 2018. Genetic evaluation of Pinus taeda clones from somatic
#' embryogenesis and their genotype x environment interaction. Crop Breed. Appl.
#' Biotechnol. 18:55-64.
#' \doi{10.1590/1984-70332018v18n1a8}
#'
#' Azevedo Peixoto, L. de, P.E. Teodoro, L.A. Silva, E.V. Rodrigues, B.G.
#' Laviola, and L.L. Bhering. 2018. Jatropha half-sib family selection with high
#' adaptability and genotypic stability. PLoS One 13:e0199880.
#' \doi{10.1371/journal.pone.0199880}
#'
#' Eberhart, S.A., and W.A. Russell. 1966. Stability parameters for comparing
#' Varieties. Crop Sci. 6:36-40.
#' \doi{10.2135/cropsci1966.0011183X000600010011x}
#'
#' Fox, P.N., B. Skovmand, B.K. Thompson, H.J. Braun, and R. Cormier. 1990.
#' Yield and adaptation of hexaploid spring triticale. Euphytica 47:57-64.
#' \doi{10.1007/BF00040364}
#'
#' Huehn, V.M. 1979. Beitrage zur erfassung der phanotypischen stabilitat. EDV
#' Med. Biol. 10:112.
#'
#' Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, V.S. Marchioro, V.Q. de
#' Souza, and E. Jost. 2019a. Mean performance and stability in
#' multi-environment trials I: Combining features of AMMI and BLUP techniques.
#' Agron. J. 111:2949-2960.
#' \doi{10.2134/agronj2019.03.0220}
#'
#' Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, B.G. Sari, and M.I. Diel.
#' 2019b. Mean performance and stability in multi-environment trials II:
#' Selection based on multiple traits. Agron. J. 111:2961-2969.
#' \doi{10.2134/agronj2019.03.0221}
#'
#' Purchase, J.L., H. Hatting, and C.S. van Deventer. 2000.
#' Genotype vs environment interaction of winter wheat (Triticum aestivum L.)
#' in South Africa: II. Stability analysis of yield performance. South African
#' J. Plant Soil 17:101-107.
#' \doi{10.1080/02571862.2000.10634878}
#'
#' Resende MDV (2007) Matematica e estatistica na analise de experimentos e no
#' melhoramento genetico. Embrapa Florestas, Colombo
#'
#' Sneller, C.H., L. Kilgore-Norquest, and D. Dombek. 1997. Repeatability of
#' Yield Stability Statistics in Soybean. Crop Sci. 37:383-390.
#' \doi{10.2135/cropsci1997.0011183X003700020013x}
#'
#' Mohammadi, R., & Amri, A. (2008). Comparison of parametric and non-parametric
#' methods for selecting stable and adapted durum wheat genotypes in variable
#' environments. Euphytica, 159(3), 419-432.
#' \doi{10.1007/s10681-007-9600-6}
#'
#' Wricke, G. 1965. Zur berechnung der okovalenz bei sommerweizen und hafer. Z.
#' Pflanzenzuchtg 52:127-138.
#'
#' Zali, H., E. Farshadfar, S.H. Sabaghpour, and R. Karimizadeh. 2012.
#' Evaluation of genotype vs environment interaction in chickpea using measures
#' of stability from AMMI model. Ann. Biol. Res. 3:3126-3136.
#'
#' @seealso [AMMI_indexes()], [anova_ind()], [anova_joint()], [ecovalence()],
#'   [Fox()], [gai()], [gamem()], [gafem()], [ge_acv()], [ge_polar()]
#'   [ge_means()], [ge_reg()], [mgidi()], [mtsi()], [performs_ammi()],
#'   [Resende_indexes()], [Shukla()], [superiority()], [waas()], [waasb()]
#' @importFrom dplyr bind_rows
#' @importFrom purrr map_dfr
#' @examples
#' \donttest{
#' library(metan)
#'
#' #################### joint-regression analysis #####################
#' ge_r <- ge_reg(data_ge2,
#'                env = ENV,
#'                gen = GEN,
#'                rep =  REP,
#'                resp = c(PH, EH, CD, CL, ED))
#' get_model_data(ge_r)
#' # Significance of deviations from the regression
#' # Use gmd(), a shortcut for get_model_data
#' gmd(ge_r, "pval_f")
#'
#'
#' #################### WAASB index #####################
#' # Fitting the WAAS index
#' AMMI <- waasb(data_ge2,
#'               env = ENV,
#'               gen = GEN,
#'               rep = REP,
#'               resp = c(PH, ED, TKW, NKR))
#'
#' # Getting the weighted average of absolute scores
#' gmd(AMMI, what = "WAASB")
#'
#' # And the rank for the WAASB index.
#' gmd(AMMI, what = "OrWAASB")
#'
#'
#' #################### BLUP model #####################
#' # Fitting a mixed-effect model
#' # Genotype and interaction as random
#' blup <- gamem_met(data_ge2,
#'                   env = ENV,
#'                   gen = GEN,
#'                   rep = REP,
#'                   resp = c(PH, ED, TKW, NKR))
#'
#' # Getting p-values for likelihood-ratio test
#' gmd(blup, what = "lrt")
#'
#' # Getting the variance components
#' gmd(blup, what = "vcomp")
#'
#' # Getting the genetic parameters
#' gmd(blup)
#'
#' ### BLUP-based stability indexes ###
#' blup %>%
#' Resende_indexes() %>%
#' gmd("HMRPGV_R")
#'
#'}
#'
get_model_data <- function(x,
                           what = NULL,
                           type = "GEN",
                           verbose = TRUE) {
  call_f <- match.call()
  if (!has_class(x, c("waasb", "waasb_group", "waas","waas_means", "gamem", "performs_ammi",
                      "Res_ind", "AMMI_indexes", "ecovalence", "ge_reg", "Fox", "Shukla",
                      "superiority", "ge_effects", "gai", "Huehn", "Thennarasu",
                      "ge_stats", "Annicchiarico", "Schmildt", "ge_means", "anova_joint",
                      "gafem", "gafem_group", "gamem_group", "anova_ind", "gge", "can_cor",
                      "can_cor_group", "gytb", "ge_acv", "ge_polar", "mgidi", "mtsi",
                      "env_stratification", "fai_blup", "sh"))) {
    stop("Invalid class in object ", call_f[["x"]], ". See ?get_model_data for more information.")
  }
  if (!is.null(what) && substr(what, 1, 2) == "PC") {
    npc <- ncol(x[[1]][["model"]] %>%
                  select(starts_with("PC")) %>%
                  select(matches("PC\\d+")))
    npcwhat <- as.numeric(substr(what, 3, nchar(what)))
    if (npcwhat > npc) {
      stop("The number of principal components informed is greater than those in model (", npc, ").", call. = FALSE)
    }
  }
  check <- c(
    "blupg", "blupge", "Y", "WAASB", "PctResp", "PctWAASB", "wRes", "wWAASB", "OrResp", "OrWAASB",
    "OrPC1", "WAASBY", "OrWAASBY", "vcomp", "lrt", "details", "genpar", "ranef", "data", "gcov",
    "gcor", "pcov", "pcor", "fixed", "h2")
  check1 <- c("Y", "WAAS", "PctResp", "PctWAAS", "wRes", "wWAAS", "OrResp", "OrWAAS", "OrPC1", "WAASY", "OrWAASY")
  check2 <- paste("PC", 1:200, sep = "")
  check3 <- c("blupg", "blupge", "vcomp", "lrt", "genpar", "details", "ranef", "data", "gcov", "gcor", "pcov", "pcor", "fixed")
  check3.1 <- c("h2", "blupg", "vcomp", "lrt", "genpar", "details", "ranef", "data", "gcov", "gcor", "pcov", "pcor", "fixed")
  check4 <- c("Y", "WAASB", "PctResp", "PctWAASB", "wRes", "wWAASB",
              "OrResp", "OrWAASB", "OrPC1", "WAASBY", "OrWAASBY")
  check5 <- c("ipca_ss", "ipca_ms", "ipca_fval", "ipca_pval", "ipca_expl", "ipca_accum")
  check6 <- c("HMGV", "HMGV_R", "RPGV", "RPGV_Y", "RPGV_R", "HMRPGV", "HMRPGV_Y", "HMRPGV_R")
  check7 <- c("ASV", "SIPC", "EV", "ZA", "WAAS", "ASV_R", "SIPC_R", "EV_R", "ZA_R", "WAAS_R", "ASV_SSI", "SIPC_SSI", "EV_SSI", "ZA_SSI", "WAAS_SSI")
  check8 <- c("Ecoval", "Ecov_perc", "rank")
  check9 <- c("GEN", "b0", "b1", "t(b1=1)", "pval_t", "s2di", "F(s2di=0)", "pval_f", "RMSE", "R2")
  check10 <- c("TOP")
  check11 <- c("ShuklaVar", "rMean", "rShukaVar", "ssiShukaVar")
  check12 <- c("Pi_a", "R_a", "Pi_f", "R_f", "Pi_u", "R_u")
  check13 <- c("GAI", "GAI_R")
  check14 <- c("S1","S1_R", "S2", "S2_R", "S3", "S3_R", "S6", "S6_R")
  check15 <- c("N1", "N1_R", "N2", "N2_R", "N3", "N3_R", "N4", "N4_R")
  check16 <- c("stats", "ranks")
  check17 <- c("Mean_rp", "Sd_rp", "Wi", "rank")
  check18 <- c("Mean_rp", "Sem_rp", "Wi", "rank")
  check19 <- c("ge_means", "env_means", "gen_means")
  check20 <- c("Y", "h2", "Sum Sq", "Mean Sq", "F value", "Pr(>F)", "fitted", "resid", "stdres", "se.fit", "details")
  check21 <- c("MEAN", "MSG", "FCG", "PFG", "MSB", "FCB", "PFB", "MSCR", "FCR", "PFCR", "MSIB_R", "FCIB_R", "PFIB_R", "MSE", "CV", "h2", "AS")
  check22 <- c("scores", "exp_var")
  check23 <- c("coefs", "loads", "crossloads", "canonical")
  check24 <- c("gyt", "stand_gyt", "si")
  check25 <- c("ACV", "ACV_R")
  check26 <- c("POLAR", "POLAR_R")
  check27 <- c("data", "cormat", "PCA", "FA", "KMO", "MSA", "communalities",
               "communalities_mean", "initial_loadings", "finish_loadings",
               "canonical_loadings", "scores_gen", "scores_ide", "gen_ide",
               "MGIDI", "contri_fac", "contri_fac_rank", "contri_fac_rank_sel",
               "sel_dif", "stat_gain", "sel_gen")
  check28 <- c("data", "cormat", "PCA", "FA", "KMO", "MSA", "communalities",
               "communalities_mean", "initial_loadings", "finish_loadings",
               "canonical_loadings", "scores_gen", "scores_ide", "gen_ide",
               "MTSI", "contri_fac", "contri_fac_rank", "contri_fac_rank_sel",
               "sel_dif_trait", "stat_dif_trait", "sel_dif_waasb", "stat_dif_waasb",
               "sel_dif_waasby", "stat_dif_waasby", "sel_gen")
  check29 <- c("FA", "env_strat", "mega_env_stat")
  check30 <- c("data", "eigen", "FA", "canonical_loadings", "FAI", "sel_dif_trait",
               "sel_gen", "construction_ideotypes")
  check31 <- c("b", "index", "sel_dif_trait", "total_gain", "sel_gen", "gcov", "pcov")
  if (!is.null(what) && what %in% check3 && !has_class(x, c("waasb", "waasb_group", "gamem", "gamem_group", "gafem", "anova_joint"))) {
    stop("Invalid argument 'what'. It can only be used with an oject of class 'waasb' or 'gamem', 'gafem, or 'anova_joint'. Please, check and fix.")
  }
  if (!type %in% c("GEN", "ENV")) {
    stop("Argument 'type' invalid. It must be either 'GEN' or 'ENV'.")
  }

  if(has_class(x, "sh")){
    if (is.null(what)){
      what <- "sel_dif_trait"
    }
    if (!what %in% check31) {
      stop("Invalid value in 'what' for object of class 'sh'. Allowed are ", paste(check31, collapse = ", "), call. = FALSE)
    }
    bind <- x[[what]]
  }

  if(has_class(x, "fai_blup")){
    if (is.null(what)){
      what <- "sel_dif_trait"
    }
    if (!what %in% check30) {
      stop("Invalid value in 'what' for object of class 'fai_blup'. Allowed are ", paste(check30, collapse = ", "), call. = FALSE)
    }
    if(what == "sel_dif_trait"){
      bind <- x[[what]][[1]]
    } else{
      bind <- x[[what]]
    }
  }

  if(has_class(x, "env_stratification")){
    if (is.null(what)){
      what <- "env_strat"
    }
    if (!what %in% check29) {
      stop("Invalid value in 'what' for object of class 'env_stratification'. Allowed are ", paste(check29, collapse = ", "), call. = FALSE)
    }
    bind <- x %>% map_dfr(~ bind_rows(!!! .x %>% .[[what]]), .id = 'TRAIT')
  }
  if(has_class(x, "mtsi")){
    if (is.null(what)){
      what <- "sel_dif_trait"
    }
    if (!what %in% check28) {
      stop("Invalid value in 'what' for object of class 'mtsi'. Allowed are ", paste(check28, collapse = ", "), call. = FALSE)
    }
    bind <- x[[what]]
  }
  if(has_class(x, "mgidi")){
    if (is.null(what)){
      what <- "sel_dif"
    }
    if (!what %in% check27) {
      stop("Invalid value in 'what' for object of class 'mgidi'. Allowed are ", paste(check27, collapse = ", "), call. = FALSE)
    }
    if(has_class(x, "mgidi_group")){
      bind <-
        x %>%
        mutate(data = map(data, ~.x %>% .[[what]])) %>%
        unnest(data)
    } else{
      bind <- x[[what]]
    }
  }
  if (has_class(x,  "ge_polar")) {
    if (is.null(what)){
      what <- "POLAR"
    }
    if (!what %in% check26) {
      stop("Invalid value in 'what' for object of class 'ge_acv'. Allowed are ", paste(check26, collapse = ", "), call. = FALSE)
    }
    bind <- sapply(x, function(x) {
      x[[what]]
    }) %>%
      as_tibble() %>%
      mutate(GEN = x[[1]][["GEN"]]) %>%
      column_to_first(GEN)
  }

  if (has_class(x,  "ge_acv")) {
    if (is.null(what)){
      what <- "ACV"
    }
    if (!what %in% check25) {
      stop("Invalid value in 'what' for object of class 'ge_acv'. Allowed are ", paste(check25, collapse = ", "), call. = FALSE)
    }
    bind <- sapply(x, function(x) {
      x[[what]]
    }) %>%
      as_tibble() %>%
      mutate(GEN = x[[1]][["GEN"]]) %>%
      column_to_first(GEN)
  }

  if(has_class(x,  "gge") & length(class(x)) == 1){
    if (is.null(what)){
      what <- "scores"
    }
    if(has_class(x,  "gge") && !what %in% check22){
      stop("Invalid value in 'what' for an object of class '", class(x), "'. Allowed are ", paste(check22, collapse = ", "), call. = FALSE)
    }
    if(what == "scores"){
      npc <- length(x[[1]]$varexpl)
      bind <- lapply(x, function(x) {
        rbind(x$ coordgen %>%
                as.data.frame() %>%
                set_names(paste("PC", 1:npc, sep = "")) %>%
                add_cols(TYPE = "GEN",
                         CODE = x$labelgen,
                         .before = 1),
              x$coordenv %>%
                as.data.frame() %>%
                set_names(paste("PC", 1:npc, sep = "")) %>%
                add_cols(TYPE = "ENV",
                         CODE = x$labelenv,
                         .before = 1))
      })
    }
    if(what == "exp_var"){
      bind <- lapply(x, function(x) {
        tibble(PC = x$labelaxes,
               Eigenvalue = x$eigenvalues,
               Variance = x$varexpl,
               Accumulated = cumsum(Variance))
      })
    }
  }

  if(has_class(x, "gytb")){
    if (is.null(what)){
      what <- "gyt"
    }
    if(has_class(x, "gytb") && !what %in% check24){
      stop("Invalid value in 'what' for an object of class '", class(x), "'. Allowed are ", paste(check24, collapse = ", "), call. = FALSE)
    }

    if(what == "gyt"){
      bind <- x[["mod"]][["data"]]
    }
    if(what == "stand_gyt"){
      bind <- x[["mod"]][["ge_mat"]]
    }
    if(what == "si"){
      bind <-
        x[["mod"]][["ge_mat"]] %>%
        as.data.frame() %>%
        add_cols(SI = rowSums(.)) %>%
        rownames_to_column("GEN") %>%
        select(GEN, SI) %>%
        arrange(-SI)
    }
  }
  if(has_class(x, c("can_cor", "can_cor_group"))){
    if (is.null(what)){
      what <- "coefs"
    }
    if(has_class(x, c("can_cor", "can_cor_group")) && !what %in% check23){
      stop("Invalid value in 'what' for an object of class '", class(x), "'. Allowed are ", paste(check23, collapse = ", "), call. = FALSE)
    }
    fg_what <- case_when(
      what == "coefs" ~ "Coef_FG",
      what == "loads" ~ "Loads_FG",
      what == "crossloads" ~ "Crossload_FG"
    )
    sg_what <- case_when(
      what == "coefs" ~ "Coef_SG",
      what == "loads" ~ "Loads_SG",
      what == "crossloads" ~ "Crossload_SG"
    )
    if(has_class(x, "can_cor_group")){
      npairs <- ncol(x[["data"]][[1]][["Coef_FG"]])
      if(what == "canonical"){
        bind <-
          x %>%
          mutate(test = map(data, ~.x %>% .[["Sigtest"]])) %>%
          remove_cols(data) %>%
          unnest(test)
      } else{
        bind <-
          rbind(
            x %>%
              mutate(FG = map(data, ~.x %>% .[[fg_what]] %>%
                                as_tibble(rownames = NA) %>%
                                set_names(paste("CP", 1:npairs, sep = "")) %>%
                                rownames_to_column("VAR"))) %>%
              remove_cols(data) %>%
              unnest(FG) %>%
              add_cols(GROUP = "FG", .before = VAR),

            x %>%
              mutate(SG = map(data, ~.x %>% .[[sg_what]] %>%
                                as_tibble(rownames = NA) %>%
                                set_names(paste("CP", 1:npairs, sep = "")) %>%
                                rownames_to_column("VAR"))) %>%
              remove_cols(data) %>%
              unnest(SG) %>%
              add_cols(GROUP = "SG", .before = VAR)
          )
      }
    } else{
      npairs <- ncol(x[["Coef_FG"]])
      if(what == "canonical"){
        bind <-
          x[["Sigtest"]] %>%
          as_tibble(rownames = NA) %>%
          rownames_to_column("GROUP")
      } else{
        bind <-
          rbind(x[[fg_what]] %>%
                  as_tibble(rownames = NA) %>%
                  set_names(paste("CP", 1:npairs, sep = "")) %>%
                  rownames_to_column("VAR") %>%
                  add_cols(GROUP = "FG", .before = VAR),
                x[[sg_what]] %>%
                  as_tibble(rownames = NA) %>%
                  set_names(paste("CP", 1:npairs, sep = "")) %>%
                  rownames_to_column("VAR") %>%
                  add_cols(GROUP = "SG", .before = VAR)
          )
      }
    }
  }

  if (has_class(x, c("waasb", "waasb_group", "gamem", "gamem_group"))) {
    if (is.null(what)){
      what <- "genpar"
    }
    if(has_class(x, c("gamem_group", "waasb_group"))){
      bind <-
        x %>%
        mutate(bind = map(data, ~.x %>% gmd(what = what, verbose = verbose))) %>%
        unnest(bind) %>%
        remove_cols(data)
    } else{
      if(is.null(x[[1]][["ESTIMATES"]]) == TRUE && what == "genpar"){
        warning("Using what = 'genpar' is only possible for models fitted with random = 'gen' or random = 'all'\nSetting what to 'vcomp'.", call. = FALSE)
        what <- "vcomp"
      }
      if(has_class(x,  "gamem") && !what %in% check3.1){
        stop("Invalid value in 'what' for an object of class '", class(x), "'. Allowed are ", paste(check3.1, collapse = ", "), call. = FALSE)
      }
      if(has_class(x,  "waasb") && !what %in% check){
        stop("Invalid value in 'what' for an object of class '", class(x), "'. Allowed are ", paste(check, collapse = ", "), call. = FALSE)
      }
      if (has_class(x,  "waasb") & what %in% check4) {
        bind <- sapply(x, function(x) {
          x$model[[what]]
        }) %>%
          as_tibble() %>%
          mutate(GEN = x[[1]][["model"]][["Code"]],
                 TYPE = x[[1]][["model"]][["type"]]) %>%
          dplyr::filter(TYPE == {{type}}) %>%
          remove_cols(TYPE) %>%
          column_to_first(GEN)
      }
      if(what == "h2"){
        bind <-
          gmd(x, verbose = FALSE) %>%
          subset(Parameters == "h2mg") %>%
          remove_cols(1) %>%
          t() %>%
          as.data.frame() %>%
          rownames_to_column("VAR") %>%
          set_names("VAR", "h2")
      }
      if (what == "data") {
        factors <- x[[1]][["residuals"]] %>% select_non_numeric_cols()
        vars <- sapply(x, function(x) {
          val <- x[["residuals"]][["Y"]]
        }) %>%
          as_tibble()
        bind <- as_tibble(cbind(factors, vars))
      }
      if (what == "gcov") {
        data <- gmd(x, "data", verbose = FALSE)
        if(ncol(select_numeric_cols(data)) < 2){
          stop("Only one numeric variable. No matrix generated.", call. = FALSE)
        }
        fctrs <- names(select_non_numeric_cols(data))
        formula <-
          x[[1]][["formula"]] %>%
          replace_string(pattern = "Y", replacement = "value") %>%
          as.formula()
        gvar <-
          data %>%
          pivot_longer(-all_of(fctrs)) %>%
          group_by(name) %>%
          doo(~lmer(formula, data = .) %>% VarCorr()) %>%
          mutate(data = as.numeric(map(data, ~ .[["GEN"]])))
        factors <- select_non_numeric_cols(data)
        combined_vars <- comb_vars(data, verbose = FALSE)
        gcov <-
          cbind(factors, combined_vars) %>%
          pivot_longer(-all_of(fctrs)) %>%
          group_by(name) %>%
          doo(~lmer(formula, data = .) %>% VarCorr()) %>%
          mutate(data = as.numeric(map(data, ~ .[["GEN"]]))) %>%
          separate(name, into = c("v1", "v2"), sep = "x") %>%
          left_join(gvar, by = c("v1" = "name")) %>%
          left_join(gvar, by = c("v2" = "name")) %>%
          mutate(gcov = (data.x - data.y - data) / 2)
        gcov_mat <- diag(gvar$data, nrow = length(gvar$data), ncol = length(gvar$data))
        colnames(gcov_mat) <- rownames(gcov_mat) <- gvar$name
        for (i in 1:nrow(gcov)){
          gcov_mat[which(rownames(gcov_mat) == as.character(gcov[i, 1])),
                   which(colnames(gcov_mat) == as.character(gcov[i, 2]))] <- pull(gcov[i, 6])
        }
        for(i in 1:nrow(gcov_mat)){
          for(j in 1:ncol(gcov_mat)){
            if(gcov_mat[i, j] == 0){
              gcov_mat[i, j] <- gcov_mat[j, i]
            } else{
              gcov_mat[i, j] <- gcov_mat[i, j]
            }
          }
        }
        bind <- make_sym(gcov_mat, diag = diag(gcov_mat), make = "lower")
        bind <- bind[names(x), names(x)]
      }
      if (what == "gcor") {
        gcov <- gmd(x, "gcov", verbose = FALSE)
        bind <- matrix(NA, nrow = nrow(gcov), ncol = ncol(gcov))
        for(i in 1:nrow(gcov)){
          for(j in 1:ncol(gcov)){
            if(i == j){
              next
            } else{
              bind[i, j] <- gcov[i, j] / sqrt(gcov[i, i] * gcov[j, j])
            }
          }
        }
        diag(bind) <- 1
        rownames(bind) <- colnames(bind) <- rownames(gcov)
      }
      if (what == "pcov") {
        data <- gmd(x, "data", verbose = FALSE)
        if(ncol(select_numeric_cols(data)) < 2){
          stop("nly one numeric variable. No matrix generated.", call. = FALSE)
        }
        bind <-
          data %>%
          means_by(GEN) %>%
          remove_cols(GEN) %>%
          cov()
      }
      if (what == "pcor") {
        pcov <- gmd(x, "pcov", verbose = FALSE)
        bind <- matrix(NA, nrow = nrow(pcov), ncol = ncol(pcov))
        for(i in 1:nrow(pcov)){
          for(j in 1:ncol(pcov)){
            if(i == j){
              next
            } else{
              bind[i, j] <- pcov[i, j] / sqrt(pcov[i, i] * pcov[j, j])
            }
          }
        }
        diag(bind) <- 1
        rownames(bind) <- colnames(bind) <- rownames(pcov)
      }
      if (what == "fixed"){
        temps <- lapply(seq_along(x), function(i) {
          x[[i]][["fixed"]] %>%
            add_cols(VAR = names(x)[i]) %>%
            column_to_first(VAR)
        })
        names(temps) <- names(x)
        bind <- temps %>% reduce(full_join, by = names(temps[[1]]))
      }
      if (what == "vcomp") {
        bind <- sapply(x, function(x) {
          val <- x[["random"]][["Variance"]]
        }) %>%
          as_tibble() %>%
          mutate(Group = x[[1]][["random"]][["Group"]]) %>%
          column_to_first(Group)
      }
      if (what == "genpar") {
        bind <- sapply(x, function(x) {
          val <- x[["ESTIMATES"]][["Values"]]
        }) %>%
          as_tibble() %>%
          mutate(Parameters = x[[1]][["ESTIMATES"]][["Parameters"]]) %>%
          column_to_first(Parameters)
      }
      if (what == "details") {
        bind <- sapply(x, function(x) {
          val <- x[["Details"]][["Values"]] %>% as.character()
        }) %>%
          as_tibble() %>%
          mutate(Parameters = x[[1]][["Details"]][["Parameters"]]) %>%
          column_to_first(Parameters)
      }
      if (what == "lrt") {
        temps <- lapply(seq_along(x), function(i) {
          x[[i]][["LRT"]] %>%
            remove_rows_na(verbose = FALSE) %>%
            add_cols(VAR = names(x)[i]) %>%
            column_to_first(VAR)
        })
        names(temps) <- names(x)
        bind <- temps %>% reduce(full_join, by = names(temps[[1]]))
      }
      if (what %in% c("blupg", "blupge")) {
        if (what == "blupg") {
          list <- lapply(x, function(x){
            x[["BLUPgen"]] %>% select(GEN, Predicted)
          })
          bind <- suppressWarnings(
            lapply(seq_along(list),
                   function(i){
                     set_names(list[[i]], "GEN", names(list)[i])
                   }) %>%
              reduce(full_join, by = "GEN") %>%
              arrange(GEN)
          )
        }
        if (what == "blupge") {
          list <- lapply(x, function(x){
            x[["BLUPint"]] %>% select(ENV, GEN, Predicted)
          })
          bind <-  suppressWarnings(
            lapply(seq_along(list),
                   function(i){
                     set_names(list[[i]], "ENV", "GEN", names(list)[i])
                   }) %>%
              reduce(full_join, by = c("ENV", "GEN")) %>%
              arrange(ENV, GEN) %>%
              means_by(ENV, GEN)
          )
        }
      }
      if (what == "ranef") {
        dfs<-
          lapply(x, function(x){
            if(has_class(x,  "waasb")){
              int <- x[["BLUPint"]]
            } else{
              int <- x[["ranef"]]
            }
            factors <- int %>% select_non_numeric_cols()
            numeric <- int %>% select_cols(contains("BLUP"))
            df_list2 <- list()
            for (i in 1:ncol(numeric)){
              temp <-
                cbind(factors, numeric[i])
              var_names <- strsplit(case_when(names(temp[ncol(factors)+1])== "BLUPg" ~ "GEN",
                                              names(temp[ncol(factors)+1])== "BLUPe" ~ "ENV",
                                              names(temp[ncol(factors)+1])== "BLUPge" ~ c("ENV GEN"),
                                              names(temp[ncol(factors)+1])== "BLUPre" ~ c("ENV REP"),
                                              names(temp[ncol(factors)+1])== "BLUPg+ge" ~ c("ENV GEN"),
                                              names(temp[ncol(factors)+1])== "BLUPbre" ~ c("REP BLOCK"),
                                              names(temp[ncol(factors)+1])== "BLUPg+bre" ~ c("GEN REP BLOCK"),
                                              names(temp[ncol(factors)+1])== "BLUPg+ge+bre" ~ c("ENV REP BLOCK GEN"),
                                              names(temp[ncol(factors)+1])== "BLUPe+ge+re+bre" ~ c("ENV REP BLOCK GEN"),
                                              names(temp[ncol(factors)+1])== "BLUPg+e+ge+re+bre" ~ c("ENV REP BLOCK GEN"),
                                              names(temp[ncol(factors)+1])== "BLUPg+e+ge+re" ~ c("ENV GEN REP"),
                                              names(temp[ncol(factors)+1])== "BLUPge+e+re" ~ c("ENV GEN REP")),
                                    " ")[[1]]
              temp <-
                temp %>%
                select(all_of(var_names), last_col()) %>%
                distinct_all(.keep_all = TRUE)
              fact_nam <- sapply(colnames(temp %>% select_non_numeric_cols()), paste) %>%
                paste(., collapse = '_')
              df_list2[[paste(fact_nam)]] <- temp
            }
            return(df_list2)
          })
        nvcomp <- length(dfs[[1]])
        bind <- list()

        for(i in 1:nvcomp){
          var_names <- names(dfs[[1]][[i]] %>%  select_non_numeric_cols())
          index <- length(var_names)
          num <-
            lapply(seq_along(dfs),
                   function(j){
                     set_names(dfs[[j]][[i]], var_names, names(dfs)[j])
                   }) %>%
            reduce(full_join, by = var_names) %>%
            arrange(across(where(~!is.numeric(.x))))
          bind[[names(dfs[[1]])[i]]] <- num
        }
      }
    }
  }

  if (has_class(x, "anova_ind")) {
    if (is.null(what)){
      what <- "MEAN"
    }
    if (!what %in% c(check21)) {
      stop("Invalid value in 'what' for object of class, ", class(x), ". Allowed are ", paste(check21, collapse = ", "), call. = FALSE)
    }
    bind <- sapply(x, function(x) {
      x[["individual"]][[what]]
    }) %>%
      as_tibble() %>%
      mutate(ENV = x[[1]][["individual"]][["ENV"]]) %>%
      column_to_first(ENV)
  }
  if (has_class(x, c("anova_joint", "gafem", "gafem_group"))) {
    if(has_class(x, c("gafem_group"))){
      bind <-
        x %>%
        mutate(bind = map(data, ~.x %>% gmd(what = what, verbose = verbose))) %>%
        unnest(bind) %>%
        remove_cols(data)
    } else{
    if (is.null(what)){
      what <- "fitted"
    }
    if (!what %in% c(check20)) {
      stop("Invalid value in 'what' for object of class, ", class(x), ". Allowed are ", paste(check20, collapse = ", "), call. = FALSE)
    }
    if(what %in% c("Sum Sq", "Mean Sq", "F value", "Pr(>F)")){
      bind <- sapply(x, function(x) {
        x[["anova"]][[what]]
      }) %>%
        as_tibble()
      bind <- cbind(x[[1]][["anova"]] %>% select_non_numeric_cols(), bind) %>%
        remove_rows_na(verbose = FALSE)
    }
    if(what  == "h2"){
      bind <- sapply(x, function(x){
        MSG <- as.numeric(x[["anova"]][which(x[["anova"]][["Source"]] == "GEN"), 4])
        MSE <- as.numeric(x[["anova"]][which(x[["anova"]][["Source"]] == "Residuals"), 4])
        (MSG - MSE) / MSG
      }) %>%
        as.data.frame() %>%
        rownames_to_column("VAR") %>%
        set_names("VAR", "h2")
    }
    if(what %in% c("Y", "fitted", "resid", "stdres", "se.fit")){
      bind <- sapply(x, function(x){
        x[["augment"]][[what]]
      }) %>%  as_tibble()
      bind <- cbind(x[[1]][["augment"]] %>% select_non_numeric_cols(), bind) %>%
        as_tibble()
    }
    if(what == "details"){
      bind <- sapply(x, function(x){
        x[["details"]][[2]]
      }) %>%
        as_tibble() %>%
        mutate(Parameters = x[[1]][["details"]][["Parameters"]]) %>%
        column_to_first(Parameters)
    }
    }
  }

  if(has_class(x,  "ge_means")){
    if (is.null(what)){
      what <- "ge_means"
    }
    if (!what %in% c(check19)) {
      stop("Invalid value in 'what' for object of class 'ge_means'. Allowed are ", paste(check19, collapse = ", "), call. = FALSE)
    }
    if(what == "ge_means"){
      bind <- sapply(x, function(x) {
        x[["ge_means_long"]][["Mean"]]
      }) %>%
        as_tibble() %>%
        add_cols(ENV = x[[1]][["ge_means_long"]][["ENV"]],
                 GEN = x[[1]][["ge_means_long"]][["GEN"]]) %>%
        column_to_first(ENV, GEN)
    }
    if(what == "env_means"){
      bind <- sapply(x, function(x) {
        x[["env_means"]][["Mean"]]
      }) %>%
        as_tibble() %>%
        add_cols(ENV = x[[1]][["env_means"]][["ENV"]]) %>%
        column_to_first(ENV)
    }
    if(what == "gen_means"){
      bind <- sapply(x, function(x) {
        x[["gen_means"]][["Mean"]]
      }) %>%
        as_tibble() %>%
        add_cols(GEN = x[[1]][["gen_means"]][["GEN"]]) %>%
        column_to_first(GEN)
    }
  }

  if (has_class(x,  "Annicchiarico")) {
    if (is.null(what)){
      what <- "Wi"
    }
    if (!what %in% c(check17)) {
      stop("Invalid value in 'what' for object of class 'Annicchiarico'. Allowed are ", paste(check17, collapse = ", "), call. = FALSE)
    }
    bind <- sapply(x, function(x) {
      x[["general"]][[what]]
    }) %>%
      as_tibble() %>%
      mutate(GEN = x[[1]][["general"]][["GEN"]]) %>%
      column_to_first(GEN)
  }

  if (has_class(x,  "Schmildt")) {
    if (is.null(what)){
      what <- "Wi"
    }
    if (!what %in% c(check18)) {
      stop("Invalid value in 'what' for object of class 'Schmildt'. Allowed are ", paste(check18, collapse = ", "), call. = FALSE)
    }
    bind <- sapply(x, function(x) {
      x[["general"]][[what]]
    }) %>%
      as_tibble() %>%
      mutate(GEN = x[[1]][["general"]][["GEN"]]) %>%
      column_to_first(GEN)
  }

  if (has_class(x,  "ge_stats")) {
    if (is.null(what)){
      what <- "stats"
    }
    if (!what  %in%  check16) {
      stop("Invalid value in 'what' for object of class 'ge_stats'. Allowed are ", paste(check16, collapse = ", "), call. = FALSE)
    }
    bind <- do.call(cbind, lapply(x, function(x) {
      if(what == "stats"){
        x %>% select(-contains("_R"), -contains("GEN"))
      } else{
        x %>% select(contains("_R"))
      }
    })) %>%
      as_tibble() %>%
      mutate(GEN = x[[1]][["GEN"]]) %>%
      pivot_longer(cols = contains(".")) %>%
      separate(name, into = c("var", "stat"), sep = "\\.") %>%
      pivot_wider(values_from = value, names_from = stat) %>%
      column_to_first(var) %>%
      arrange(var)
  }

  if (has_class(x,  "Thennarasu")) {
    if (is.null(what)){
      what <- "N1"
    }
    if (!what %in% c(check15)) {
      stop("Invalid value in 'what' for object of class 'Thennarasu'. Allowed are ", paste(check15, collapse = ", "), call. = FALSE)
    }
    bind <- sapply(x, function(x) {
      x[[what]]
    }) %>%
      as_tibble() %>%
      mutate(GEN = x[[1]][["GEN"]]) %>%
      column_to_first(GEN)
  }

  if (has_class(x,  "Huehn")) {
    if (is.null(what)){
      what <- "S1"
    }
    if (!what %in% c(check14)) {
      stop("Invalid value in 'what' for object of class 'Huehn'. Allowed are ", paste(check14, collapse = ", "), call. = FALSE)
    }
    bind <- sapply(x, function(x) {
      x[[what]]
    }) %>%
      as_tibble() %>%
      mutate(GEN = x[[1]][["GEN"]]) %>%
      column_to_first(GEN)
  }

  if (has_class(x,  "gai")) {
    if (is.null(what)){
      what <- "GAI"
    }
    if (!what %in% check13) {
      stop("Invalid value in 'what' for object of class 'gai'. Allowed are ", paste(check13, collapse = ", "), call. = FALSE)
    }
    bind <- sapply(x, function(x) {
      x[[what]]
    }) %>%
      as_tibble() %>%
      mutate(GEN = x[[1]][["GEN"]]) %>%
      column_to_first(GEN)
  }

  if (has_class(x,  "ge_effects")) {
    bind <- sapply(x, function(x) {
      make_long(x)[[3]]
    }) %>%
      as_tibble()
    factors <- x[[1]] %>%
      make_long() %>%
      select(1:2)
    bind <- cbind(factors, bind)
  }

  if (has_class(x,  "superiority")) {
    if (is.null(what)){
      what <- "Pi_a"
    }
    if (!what %in% c(check12)) {
      stop("Invalid value in 'what' for object of class 'superiority'. Allowed are ", paste(check12, collapse = ", "), call. = FALSE)
    }
    bind <- sapply(x, function(x) {
      x[["index"]][[what]]
    }) %>%
      as_tibble() %>%
      mutate(GEN = x[[1]][["index"]][["GEN"]]) %>%
      column_to_first(GEN)
  }

  if (has_class(x,  "Shukla")) {
    if (is.null(what)){
      what <- "ShuklaVar"
    }
    if (!what %in% c(check11)) {
      stop("Invalid value in 'what' for object of class 'Shukla'. Allowed are ", paste(check11, collapse = ", "), call. = FALSE)
    }
    bind <- sapply(x, function(x) {
      x[[what]]
    }) %>%
      as_tibble() %>%
      mutate(GEN = x[[1]][["GEN"]]) %>%
      column_to_first(GEN)
  }

  if (has_class(x,  "Fox")) {
    if (is.null(what)){
      what <- "TOP"
    }
    if (!what %in% c(check10)) {
      stop("Invalid value in 'what' for object of class 'Fox'. Allowed are ", paste(check10, collapse = ", "), call. = FALSE)
    }
    bind <- sapply(x, function(x) {
      x[[what]]
    }) %>%
      as_tibble() %>%
      mutate(GEN = x[[1]][["GEN"]]) %>%
      column_to_first(GEN)
  }

  if (has_class(x,  "ge_reg")) {
    if (is.null(what)){
      what <- "b1"
    }
    if (!what %in% c(check9)) {
      stop("Invalid value in 'what' for object of class 'ge_reg'. Allowed are ", paste(check9, collapse = ", "), call. = FALSE)
    }
    bind <- sapply(x, function(x) {
      x[["regression"]][[what]]
    }) %>%
      as_tibble() %>%
      mutate(GEN = x[[1]][["regression"]][["GEN"]]) %>%
      column_to_first(GEN)
  }

  if (has_class(x,  "ecovalence")) {
    if (is.null(what)){
      what <- "Ecoval"
    }
    if (!what %in% check8) {
      stop("Invalid value in 'what' for object of class 'ecovalence'. Allowed are ", paste(check8, collapse = ", "), call. = FALSE)
    }
    bind <- sapply(x, function(x) {
      x[[what]]
    }) %>%
      as_tibble() %>%
      mutate(GEN = x[[1]][["GEN"]]) %>%
      column_to_first(GEN)
  }

  if (has_class(x,  "AMMI_indexes")) {
    if (is.null(what)){
      what <- "WAAS"
    }
    if (!what %in% c(check7)) {
      stop("Invalid value in 'what' for object of class 'AMMI_indexes'. Allowed are ", paste(check7, collapse = ", "), call. = FALSE)
    }
    bind <- sapply(x, function(x) {
      x[[what]]
    }) %>%
      as_tibble() %>%
      mutate(GEN = x[[1]][["GEN"]]) %>%
      column_to_first(GEN)
  }

  if (has_class(x,  "Res_ind")) {
    if (is.null(what)){
      what <- "HMRPGV"
    }
    if (!what %in% c(check6)) {
      stop("Invalid value in 'what' for object of class 'Res_ind'. Allowed are ", paste(check6, collapse = ", "), call. = FALSE)
    }
    bind <- sapply(x, function(x) {
      x[[what]]
    }) %>%
      as_tibble() %>%
      mutate(GEN = x[[1]][["GEN"]]) %>%
      column_to_first(GEN)
  }

  if (has_class(x,  "performs_ammi")) {
    if (is.null(what)){
      what <- "ipca_expl"
    }
    if (!what %in% c("Y", check2, check5)) {
      stop("Invalid value in 'what' for object of class 'performs_ammi'.")
    }
    if (what == "Y" | what %in% check2) {
      bind <- sapply(x, function(x) {
        x$model[[what]]
      }) %>%
        as_tibble() %>%
        mutate(GEN = x[[1]][["model"]][["Code"]],
               TYPE = x[[1]][["model"]][["type"]]) %>%
        dplyr::filter(TYPE == {{type}}) %>%
        remove_cols(TYPE) %>%
        column_to_first(GEN)
    }
    if (what  %in% check5) {
      what <- case_when(
        what == "ipca_ss" ~ "Sum Sq",
        what == "ipca_ms" ~ "Mean Sq",
        what == "ipca_fval" ~ "F value",
        what == "ipca_pval" ~ "Pr(>F)",
        what == "ipca_expl" ~ "Proportion",
        what == "ipca_accum" ~ "Accumulated"
      )
      bind <- sapply(x, function(x) {
        val <- x[["PCA"]][[what]]
      }) %>%
        as_tibble() %>%
        mutate(PC = x[[1]][["PCA"]][["PC"]],
               DF = x[[1]][["PCA"]][["Df"]]) %>%
        column_to_first(PC, DF)
    }
  }

  if (has_class(x, c("waas", "waas_means"))){
    if (is.null(what)){
      what <- "WAAS"
    }
    if (!what %in% c(check1, check2)) {
      stop("Invalid value in 'what' for object of class '", class(x), "'. Allowed are ", paste(check1, collapse = ", "), call. = FALSE)
    }
    if (what %in% check1 | what  %in% check2) {
      bind <- sapply(x, function(x) {
        x$model[[what]]
      }) %>%
        as_tibble() %>%
        mutate(GEN = x[[1]][["model"]][["Code"]],
               TYPE = x[[1]][["model"]][["type"]]) %>%
        dplyr::filter(TYPE == {{type}}) %>%
        remove_cols(TYPE) %>%
        column_to_first(GEN)
    }
    if (what == "details") {
      bind <- sapply(x, function(x) {
        val <- x[["Details"]][["Values"]] %>% as.character()
      }) %>%
        as_tibble() %>%
        mutate(Parameters = x[[1]][["Details"]][["Parameters"]]) %>%
        column_to_first(Parameters)
    }
  }

  if(verbose == TRUE){
    message("Class of the model: ", paste(class(x), collapse = ", "))
    message("Variable extracted: ", what)
  }
  return(bind)
}

#' @name get_model_data
#' @export
gmd <- function(x,
                what = NULL,
                type = "GEN",
                verbose = TRUE){
  get_model_data(x, what, type, verbose)
}

#' @name get_model_data
#' @export
sel_gen <- function(x){
  gmd(x, "sel_gen", verbose = FALSE)
}
