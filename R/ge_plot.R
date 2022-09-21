#' Graphical analysis of genotype-vs-environment interaction
#' @description
#' `r badge('stable')`
#'
#' This function produces a line plot for a graphical interpretation of the
#' genotype-vs-environment interaction. By default, environments are in the x
#' axis whereas the genotypes are depicted by different lines. The y axis
#' contains the value of the selected variable. A heatmap can also be created.
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param resp The response variable.
#' @param type The type of plot `type = 1` for a heatmap or `type = 2`
#'   for a line plot.
#' @param values Show the values in the plot? Defaults to `TRUE`.
#' @param text_row_pos,text_col_pos The position of the text in the
#'   rows and columns. The defaults show the text at left and top.
#' @param average Show the average values for environments and genotypes?
#'   Defaults to `TRUE`.
#' @param order_g,order_e A charactere vector indicating the order of the levels
#'   for genotypes and environments, respectively. This can be used to change
#'   the default ordering of rows and columns.
#' @param xlab,ylab The labels for x and y axis, respectively.
#' @param width_bar,heigth_bar The width and heigth of the legend bar,
#'   respectively.
#' @param plot_theme The graphical theme of the plot. Default is
#'   `plot_theme = theme_metan()`. For more details,see
#'   [ggplot2::theme()].
#' @param colour Logical argument. If `FALSE` then the plot will not be
#'   colored.
#' @return An object of class `gg, ggplot`.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' ge_plot(data_ge2, ENV, GEN, PH)
#' ge_plot(data_ge, ENV, GEN, GY, type = 2)
#'}
ge_plot <- function(.data,
                    env,
                    gen,
                    resp,
                    type = 1,
                    values = TRUE,
                    text_col_pos = c("top", "bottom"),
                    text_row_pos = c("left", "right"),
                    average = TRUE,
                    order_g = NULL,
                    order_e = NULL,
                    xlab = NULL,
                    ylab = NULL,
                    width_bar = 1.5,
                    heigth_bar = 15,
                    plot_theme = theme_metan(),
                    colour = TRUE) {
  text_col_pos <- rlang::arg_match(text_col_pos)
  text_row_pos <- rlang::arg_match(text_row_pos)
  if(type == 1){
    if(isTRUE(average)){
      mat <-
        select_cols(.data,
                    ENV = {{env}},
                    GEN = {{gen}},
                    Y = {{resp}}) |>
        make_mat(GEN, ENV, Y) |>
        row_col_mean()

      colnames(mat)[ncol(mat)] <- "Average"
      rownames(mat)[nrow(mat)] <- "Average"
      if(is.null(order_e)){
        order_e <- colnames(mat)[-ncol(mat)]
      } else{
        order_e <- order_e
      }
      if(is.null(order_g)){
        order_g <- rownames(mat)[-nrow(mat)]
      } else{
        order_g <- order_g
      }
      df_long <-
        make_long(mat) |>
        as_factor(1:2) |>
        mutate(ENV = factor(ENV, levels = c(order_e, "Average")),
               GEN = factor(GEN, levels = c("Average", order_g)))
    } else{
      df_long <-
        select_cols(.data,
                    ENV = {{env}},
                    GEN = {{gen}},
                    Y = {{resp}}) |>
        means_by(ENV, GEN)
      if(is.null(order_e)){
        order_e <- levels(df_long$ENV)
      } else{
        order_e <- order_e
      }
      if(is.null(order_g)){
        order_g <- levels(df_long$GEN)
      } else{
        order_g <- order_g
      }
      df_long <-
        df_long |>
        mutate(ENV = factor(ENV, levels = order_e),
               GEN = factor(GEN, levels = order_g))
    }
    p <-
      ggplot(df_long, aes(ENV, GEN, fill = Y)) +
      geom_tile(color = "black")+
      {if(text_row_pos == "left")
        scale_y_discrete(expand = expansion(mult = c(0,0)))}+
      {if(text_row_pos != "left")
        scale_y_discrete(expand = expansion(mult = c(0,0)),
                         position = "right")}+
      {if(text_col_pos != "top")
        scale_x_discrete(expand = expansion(mult = c(0,0)))} +
      {if(text_col_pos == "top")
        scale_x_discrete(position = "top",
                         expand = expansion(0))} +
      scale_fill_viridis_c() +
      {if(isTRUE(values))geom_text(aes(label = round(Y, 1)),
                                   color = "black",
                                   size = 3)} +
      guides(fill = guide_colourbar(label = TRUE,
                                    draw.ulim = TRUE,
                                    draw.llim = TRUE,
                                    frame.colour = "black",
                                    ticks = TRUE,
                                    nbin = 10,
                                    label.position = "right",
                                    barwidth = width_bar,
                                    barheight = heigth_bar,
                                    direction = 'vertical'))+
      plot_theme %+replace%
      theme(legend.position = "right",
            legend.title = element_blank()) +
      labs(x = xlab,
           y = ylab)

    if(text_col_pos == "top"){
      p <- p + theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0))
    } else{
      p <- p + theme(axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1))
    }
  }
  if(type == 2){
    p <- ggplot(.data, aes(x = {{env}}, y = {{resp}}))
    if (colour == TRUE) {
      p <- p +
        stat_summary(aes(colour = {{gen}},
                         group = {{gen}}),
                     fun = mean,
                     geom = "line",
                     na.rm = TRUE)
    } else {
      p <- p +
        stat_summary(aes(group = {{gen}}),
                     fun = mean,
                     geom = "line",
                     colour = "black",
                     na.rm = TRUE)
    }
    p <- p + geom_point(stat = "summary",
                        fun = mean,
                        size = 3,
                        shape = 18) +
      plot_theme %+replace%
      theme(legend.position = "right")
  }

  return(p +
           labs(x = xlab,
                y = ylab))
}
