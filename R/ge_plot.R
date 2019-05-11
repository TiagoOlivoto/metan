#' Graphical analysis of genotype-vs-environment interaction
#'
#' This function produces a line plot for a graphical interpretation of the
#' genotype-vs-environment interaction. By default, environments are in the
#' x axis whereas the genotypes are depicted by different lines. The y axis
#' contains the value of the selected variable.
#'
#' @param .data The dataset containing the columns related to Environments,
#' Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#' environments
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param resp The response variable.
#' @param theme The graphical theme of the plot. Default is \code{theme =
#' theme_waasb()}. Please, see \code{?WAASB::theme_waasb}. An own theme can be
#' applied using the arguments: \code{theme = theme(some stuff
#' here)}. For more details, please, see \code{?ggplot2::theme}.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#' library(METAAB)
#' library(ggplot2)
#' ge_plot(data_ge2, ENV, GEN, PH)
#' ge_plot(data_ge2, ENV, GEN, PH)+
#' theme_gray()

ge_plot = function(.data, env, gen, resp, theme = theme_waasb()){
ggplot(.data, aes(x = !!enquo(env), y = !!enquo(resp))) +
      geom_point(aes(colour = !!enquo(gen), group = !!enquo(gen)), stat='summary', fun.y=mean) +
      stat_summary(aes(colour = !!enquo(gen), group = !!enquo(gen)), fun.y=mean, geom="line")+
      geom_point(stat = "summary", fun.y = mean, size = 3, shape = 18) +
      theme %+replace%
      theme(legend.positio = "right")
}
