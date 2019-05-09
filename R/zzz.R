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
#' @importFrom gtools  mixedorder
#' @importFrom dplyr  select  everything  mutate  group_by  group_by_if  group_keys
#'             group_split  left_join  sample_n  anti_join  enquo  pull  summarise_all
#'             select_if
#' @importFrom gplots  heatmap.2
#' @importFrom rlang  eval_bare  expr
#' @importFrom FWDselect  selection
#' @importFrom grDevices  colorRampPalette  dev.off  pdf
#'             chull  tiff  boxplot.stats
#' @importFrom methods  is  as
#' @importFrom graphics  plot  boxplot  hist  par  mtext
#' @importFrom utils  head  setWinProgressBar  winProgressBar  combn  stack
#' @importFrom methods setClass setGeneric setMethod setRefClass

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
        "eq", "IndAmb", "REP", "gge", "ind", "cophenetic", "remaining", "index", "ge"))
}
