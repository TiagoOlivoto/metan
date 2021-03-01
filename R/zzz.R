#' @import ggplot2
#' @import stats
#' @import mathjaxr
#' @importFrom dplyr  select  everything  mutate  group_by group_keys
#'   group_split  left_join  sample_n  anti_join pull select_if arrange slice
#'   contains top_n summarise ungroup rename is_grouped_df desc n funs
#' @importFrom tibble tibble as_tibble
#' @importFrom ggrepel  geom_text_repel geom_label_repel
#' @importFrom magrittr  %>% %<>%
#' @importFrom GGally  ggally_text  eval_data_col  ggpairs
#' @importFrom lme4  ranef  VarCorr  fortify.merMod
#' @importFrom lmerTest  ranova  lmer
#' @importFrom grDevices  colorRampPalette  dev.off  pdf
#'             chull  tiff  boxplot.stats
#' @importFrom methods  is  as
#' @importFrom graphics  plot  boxplot  hist  par  mtext abline
#' @importFrom utils  head combn  stack data
#' @importFrom methods setClass setGeneric setMethod setRefClass
#'
NULL

#' @title Data for examples
#'
#' @description This dataset contains the means for grain yield of 10 genotypes cultivated
#' in 5 environments. The interaction effects for this data is found in
#' [int.effects()]
#'
#' @name meansGxE
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords data
NULL

#' Data for examples
#'
#' @name int.effects
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords data
NULL


#' @title Multi-environment trial of oat
#' @description This dataset contain data on two variables assessed in 10 genotypes growing
#' in in 11 environments. The experimental design was a RCBD with 3
#' replicates(blocks). This data provide examples for several functions of
#' \pkg{metan} package.
#' @format A tibble with 420 observations on the following 5 variables.
#' * **ENV** A factor with 14 levels; each level represents one cultivation
#' environment.
#' * **GEN** A factor with 10 levels; each level represents one genotype.
#' * **REP** A factor with 3 levels; each level represents one
#' replication/block.
#' * **GY** A continuous variable (grain yield) observed in each plot.
#' * **HM** A continuous variable (hectoliter mass) observed in each plot.
#' @md
#' @source Personal data
#' @name data_ge
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords data
NULL

#' @title Multi-environment trial of maize
#' @description This dataset contain data on 15 traits assessed in 13 maize hybrids growing
#' in 4 environments. The experimental design was a RCBD with 3 blocks and 1
#' replications per block. It may be used as example in several functions of
#' \pkg{metan} package.
#' @format A tibble with 156 observations on the following 18 variables.
#' * **ENV** A factor with 4 levels; each level represents one cultivation
#' environment.
#' * **GEN** A factor with 13 levels; each level represents one maize hybrid.
#' * **REP** A factor with 3 levels; each level represents one
#' replication/block.
#' * **PH** Plant height, in cm.
#' * **EH** Ear height, in cm.
#' * **EP** Ear position, i.e., the ratio EH/PH.
#' * **EL** Ear length, in cm.
#' * **ED** Ear diameter, in mm.
#' * **CL** Cob length, in cm.
#' * **CD** Cob diameter, in mm.
#' * **CW** Cob weight, in g.
#' * **KW** Kernel weight, in cm.
#' * **NR** Number of rows.
#' * **NKR** Number of kernels per row.
#' * **CDED** Cob diameter / Ear diameter ratio.
#' * **PERK** Percentage of kernels.
#' * **TKW** Thousand-kernel weight
#' * **NKE** Number of kernels per row.
#' @md
#' @name data_ge2
#' @source Personal data
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords data
NULL



#' @title Single maize trial
#' @description This dataset contain data on 15 traits assessed in 13 maize hybrids.
#' The experimental design was a RCBD with 3 blocks and 1
#' replications per block. It is used as an example in the function [gamem()]
#'  of the \pkg{metan} package.
#' @format A tibble with 39 observations on the following 17 variables.
#' * **GEN** A factor with 13 levels; each level represents one maize hybrid.
#' * **REP** A factor with 3 levels; each level represents one replication/block.
#' * **PH** Plant height, in cm.
#' * **EH** Ear height, in cm.
#' * **EP** Ear position, i.e., the ratio EH/PH.
#' * **EL** Ear length, in cm.
#' * **ED** Ear diameter, in mm.
#' * **CL** Cob length, in cm.
#' * **CD** Cob diameter, in mm.
#' * **CW** Cob weight, in g.
#' * **KW** Kernel weight, in cm.
#' * **NR** Number of rows.
#' * **NKR** Number of kernels per row.
#' * **CDED** Cob diameter / Ear diameter ratio.
#' * **PERK** Percentage of kernels.
#' * **TKW** Thousand-kernel weight
#' * **NKE** Number of kernels per row.
#' @md
#' @name data_g
#' @source Personal data
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords data
NULL



#' @title Data from an alpha lattice design
#' @description Alpha lattice design of spring oats
#' @format A tibble with 72 observations on the following 5 variables.
#' * **PLOT** Plot number
#' * **REP** Replicate code
#' * **BLOCK** Incomplete block code
#' * **GEN** Genotype code
#' * **YIELD** Observed dry matter yield (tonnes/ha)
#'
#' @details A spring oats trial grown in Craibstone. There were 24 varieties
#'  in 3 replicates, each consisting of 6 incomplete blocks of 4 plots.
#'  Planted in a resolvable alpha design. The plots were laid out in a single line.
#' @source J. A. John & E. R. Williams (1995). Cyclic and computer generated designs,
#'  Chapman and Hall, London. Page 146.
#' @name data_alpha
#' @md
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords data
NULL

.onAttach <- function(libname, pkgname) {
  vers <-  "v1.12.0.9000"
  packageStartupMessage("|========================================================|")
  packageStartupMessage("| Multi-Environment Trial Analysis (metan) ", vers, "  |")
  packageStartupMessage("| Author: Tiago Olivoto                                  |")
  packageStartupMessage("| Type 'citation('metan')' to know how to cite metan     |")
  packageStartupMessage("| Type 'vignette('metan_start')' for a short tutorial    |")
  packageStartupMessage("| Visit 'https://bit.ly/2TIq6JE' for a complete tutorial |")
  packageStartupMessage("|========================================================|")
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("..density..", ".fitted", ".resid", ".scresid", "Accumul",
      "Accumulated", "BLUPe", "BLUPg", "BLUPge", "CI", "Code", "ENV", "Eigenvalue",
      "FAI", "GEN", "GY_HMRPGV", "GY_RPGV", "Genotype", "MODEL", "Mean", "PC", "PC1",
      "PC2", "Parameters", "PctResp", "PctWAAS", "PctWAASB", "Percent", "PesRes",
      "PesWAAS", "PesWAASB", "Predicted", "Proportion", "RMSPD", "Rank", "ResAMMI",
      "WAAS", "WAASBY", "WAASY", "X1", "X2", "Y", "Ypred", "YpredAMMI", "envPC1",
      "factors", "gen", "genPC1", "ggee", "id", "label", "nominal", "pred", "resOLS",
      "sel", "stdres", "type", "y", "z", "Corr", "combn", "cov2cor", "partial",
      "pt", ".", "linear", "my_custom_cor", "my_custom_smooth", "K", "direct", "VAR",
      "eq", "IndAmb", "REP", "gge", "ind", "cophenetic", "remaining", "index", "ge",
      "FA1", "FA2", "Gen", "wRes", "wWAASB", "OrResp", "OrPC1", "OrWAASB", "wWAAS",
      "OrWAAS", ".stdresid", "WAASB", "grp", "Names", "ID", "MTSI", "Pair", "LL", "UL",
      "ci", "mean_var", "se", "x", "d1", "d2", "radio", "x0", "x1_x", "x1_y", "y0",
      "val", "Statistic", "Cluster", "CCP", "U1", "U2", "V1", "V2", "Var", "LEVEL",
      "GROUP", "Group", "model", "Variance (%)", "Xs", "Xo", "Pr(>Chisq)", "dat",
      "BLOCK", "rowid", "grank", "rMean", "rShukaVar", "value", "DF", "without", "RPGV_Y",
      "HMRPGV_Y", "Pairs", "name", "v1", "v2", "slope", "new_var", "TYPE", "CODE",
      "Sum.Sq", "Mean.Sq", "F.value", "Pr.F", "Sum Sq", "Df", "Source", "BLUPbre",
      "BLUPe+ge+re+bre", "BLUPg+e+ge+re", "BLUPg+e+ge+re+bre", "BLUPg+bre", "BLUPg+ge+bre",
      "BLUPge+e+re", "BLUPre", "Estimate", "HMRPGV", "RPGV", "Variance", "blup",
      "intercept", "lower", "upper", "pred_ols", "res_ammi", "res_ols", "pattern",
      "replacement", "comparison", "group1", "group2", "p.adj", "term", "rel_freq",
      "variable", "Model", "level", "RESIDUAL", "MGIDI", "sd.amo", "SD", "SG", "data.x",
      "data.y", "sense", "win", "where", "test", "FG", "SI", "ACV", "SDperc", "MEGA_ENV",
      "group", "hjust", "theta", "vjust", "x_raw", "y_raw", "row_id", "b0", "b1", "s2di",
      "TRAIT"))
  }
