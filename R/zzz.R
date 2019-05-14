#' @import ggplot2
#' @import stats
#' @importFrom grid grobTree textGrob grid.newpage pushViewport viewport grid.layout
#' @importFrom ggrepel  geom_text_repel
#' @importFrom magrittr  %>%
#' @importFrom GGally  ggally_text  eval_data_col  ggpairs
#' @importFrom ade4  mantel.randtest
#' @importFrom dendextend  set  rotate_DendSer
#' @importFrom lme4  ranef  VarCorr  fortify.merMod
#' @importFrom lmerTest  ranova  lmer
#' @importFrom dplyr  select  everything  mutate  group_by  group_by_if  group_keys
#'             group_split  left_join  sample_n  anti_join  enquo  pull  summarise_all
#'             select_if arrange
#' @importFrom gplots  heatmap.2
#' @importFrom rlang  eval_bare  expr
#' @importFrom FWDselect  selection
#' @importFrom grDevices  colorRampPalette  dev.off  pdf
#'             chull  tiff  boxplot.stats
#' @importFrom methods  is  as
#' @importFrom graphics  plot  boxplot  hist  par  mtext
#' @importFrom utils  head  setWinProgressBar  winProgressBar  combn  stack
#' @importFrom methods setClass setGeneric setMethod setRefClass
#'
NULL

#' A dataset with means of 10 genotypes cultivated in 5 environments
#'
#' This dataset contains the means for grain yield of 10 genotypes cultivated
#' in 5 environments. The interaction effects for this data is found in
#' \link{int.effects}
#'
#' @name meansGxE
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords data
NULL

#' Dataset with genotype-by-environment interaction effects
#'
#' @name int.effects
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords data
NULL


#' A replicate-based data of 10 genotypes cultivated in 14 environments
#'
#' This dataset contain data on two variables assessed in 10 genotypes growing
#' in in 11 environments. The experimental design was a RCBD with 3
#' replicates(blocks). This data provide examples for several functions of
#' \pkg{metan} package.
#'
#' @param ENV A factor with 14 levels; each level represents one cultivation
#' environment.
#' @param GEN A factor with 10 levels; each level represents one genotype.
#' @param REP A factor with 3 levels; each level represents one
#' replication/block.
#' @param GY A continuous variable (grain yield) observed in each plot.
#' @param HM A continuous variable (hectoliter mass) observed in each plot.
#'
#' @name data_ge
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords data
NULL

#' A replicate-based data of 10 genotypes cultivated in 14 environments
#'
#' This dataset contain data on two variables assessed in 10 genotypes growing
#' at 14 environments. The experimental design in each environment was a RCBD with 3
#' replicates(blocks). This data provide examples for several functions of
#' \pkg{metan} package.
#'
#' @param ENV A factor with 14 levels; each level represents one cultivation
#' environment.
#' @param GEN A factor with 10 levels; each level represents one genotype.
#' @param REP A factor with 3 levels; each level represents one
#' replication/block.
#' @param GY A continuous variable (grain yield) observed in each plot.
#' @param HM A continuous variable (hectoliter mass) observed in each plot.
#'
#' @name data_ge
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords data
NULL


#' A replicate-based data of 13 maize hybrids cultivated in 4 environments
#'
#' This dataset contain data on 15 traits assessed in 13 maize hybrids growing
#' in 4 environments. The experimental design was a RCBD with 3 blocks and 1
#' replications per block. It may be used as example in several functions of
#' \pkg{metan} package.
#'
#' @param ENV A factor with 4 levels; each level represents one cultivation
#' environment.
#' @param GEN A factor with 13 levels; each level represents one maize hybrid.
#' @param REP A factor with 3 levels; each level represents one
#' replication/block.
#' @param PH Plant height, in cm.
#' @param EH Ear height, in cm.
#' @param EP Ear position, i.e., the ratio EH/PH.
#' @param EL Ear length, in cm.
#' @param ED Ear diameter, in mm.
#' @param CL Cob length, in cm.
#' @param CD Cob diameter, in mm.
#' @param CW Cob weight, in g.
#' @param KW Kernel weight, in cm.
#' @param NR Number of rows.
#' @param NKR Number of kernels per row.
#' @param CDED Cob diameter / Ear diameter ratio.
#' @param PERK Percentage of kernels.
#' @param TKW Thousand-kernel weight
#' @param NKE Number of kernels per row.
#'
#' @name data_ge2
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords data
NULL

if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("..density..", ".fitted", ".resid", ".scresid", "Accumul",
        "Accumulated", "BLUPe", "BLUPg", "BLUPge", "CI", "Code", "ENV", "Eigenvalue",
        "FAI", "GEN", "GY_HMRPGV", "GY_RPGV", "Genotype", "MODEL", "Mean", "PC", "PC1",
        "PC2", "Parameters", "PctResp", "PctWAAS", "PctWAASB", "Percent", "PesRes",
        "PesWAAS", "PesWAASB", "Predicted", "Proportion", "RMSPD", "Rank", "ResAMMI",
        "WAAS", "WAASBY", "WAASY", "X1", "X2", "Y", "Ypred", "YpredAMMI", "envPC1",
        "factors", "gen", "genPC1", "ggee", "id", "label", "nominal", "pred", "resOLS",
        "sel", "stdres", "type", "y", "z", "Corr", "combn", "cov2cor", "partial",
        "pt", ".", "linear", "my_custom_cor", "my_custom_smooth", "K", "direct", "VAR",
        "eq", "IndAmb", "REP", "gge", "ind", "cophenetic", "remaining", "index", "ge",
        "FA1", "FA2", "Gen"))
}
