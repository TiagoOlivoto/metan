#' @import ggplot2
#' @import stats
#' @importFrom dplyr  select  everything  mutate  group_by  group_by_if  group_keys
#'             group_split  left_join  sample_n  anti_join  enquo  pull  summarise_all
#'             select_if arrange slice contains top_n summarise ungroup rename
#'             is_grouped_df mutate_at desc tibble as_tibble n summarise_if funs
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom grid grobTree textGrob grid.newpage pushViewport viewport grid.layout
#' @importFrom ggrepel  geom_text_repel
#' @importFrom magrittr  %>% %<>%
#' @importFrom GGally  ggally_text  eval_data_col  ggpairs
#' @importFrom ade4  mantel.randtest
#' @importFrom dendextend  set  rotate_DendSer
#' @importFrom lme4  ranef  VarCorr  fortify.merMod
#' @importFrom lmerTest  ranova  lmer
#' @importFrom gplots  heatmap.2
#' @importFrom FWDselect  selection
#' @importFrom grDevices  colorRampPalette  dev.off  pdf
#'             chull  tiff  boxplot.stats
#' @importFrom methods  is  as
#' @importFrom graphics  plot  boxplot  hist  par  mtext abline
#' @importFrom utils  head combn  stack
#' @importFrom methods setClass setGeneric setMethod setRefClass
#'
NULL

#' @title Data for examples
#'
#' @description This dataset contains the means for grain yield of 10 genotypes cultivated
#' in 5 environments. The interaction effects for this data is found in
#' \link{int.effects}
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


#' @title Data for examples
#' @description This dataset contain data on two variables assessed in 10 genotypes growing
#' in in 11 environments. The experimental design was a RCBD with 3
#' replicates(blocks). This data provide examples for several functions of
#' \pkg{metan} package.
#' @param ENV A factor with 14 levels; each level represents one cultivation
#' environment.
#' @param GEN A factor with 10 levels; each level represents one genotype.
#' @param REP A factor with 3 levels; each level represents one
#' replication/block.
#' @param GY A continuous variable (grain yield) observed in each plot.
#' @param HM A continuous variable (hectoliter mass) observed in each plot.
#' @name data_ge
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords data
NULL

#' @title Data for examples
#' @description This dataset contain data on 15 traits assessed in 13 maize hybrids growing
#' in 4 environments. The experimental design was a RCBD with 3 blocks and 1
#' replications per block. It may be used as example in several functions of
#' \pkg{metan} package.
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
#' @name data_ge2
#' @docType data
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @keywords data
NULL

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
      "val", "Statistic", "Cluster", "CCP", "U1", "U2", "V1", "V2", "Var", "LEVEL"))
}
