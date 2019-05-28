#' @title Fast way to create a bar plot
#' @description Create a bar plot based on categorical (one or two) variables
#' and one numeric variable.
#' @param data The data set
#' @param ... A comma-separated list of unquoted variable names. Must be up to two
#' variables.
#' @param resp The response variable
#' @param errorbar Logical argument, set to TRUE. In this case, an error bar is shown.
#' @param stat.erbar The statistic to be shown in the errorbar. Must be one of the 'se'
#' (standard error, default), 'sd' (standard deviation), or 'ci' confidence interval, based
#' on the confidence level
#' @param level The fonfidence level
#' @param invert Logical argument. If \code{TRUE}, the order of the factors entered 
#' in changes in the graph
#' @param xlab The x lab
#' @param ylab The y lab
#' @param col Logical argument. If \code {FALSE}, a gray scale is used.
#' @param palette The color palette to be used. For more details, see 
#' \code{?Scale_colour_brewer}
#' @param width.bar The width of the bars in the graph. Defaults to 0.9
#' possible values [0-1].
#' @param cex.angle The angle of the caption text. Default is 0.
#' @param cex.hjust The horizontal adjustment of the caption text. Defaults to 0.5.
#' Use this argument to adjust the text when the angle of the text is different from 0.}
#' @param legend.position The legend position.
#' @param alpha The alpha for the color in confidence band
#' @param size.shape The size for the shape in plot
#' @param size.line The size for the line in the plot
#' @param cex The size of the text
#' @param fontfam The family of the font text
#' @param na.rm Should 'NA' values be removed to compute the statistics? 
#' Defaults to true
#' @param verbose Logical argument. If FALSE the code will run silently
#' @export
#' @seealso \code{\link{plot_line}, \code{\link{plot_factbars}}
#' 
#' 
plot_factbars = function(.data,
                         ...,
                         resp,
                         errorbar = TRUE,
                         stat.erbar = "se",
                         width.erbar = 0.3,
                         level= .95,
                         invert = FALSE,
                         xlab = NULL,
                         ylab = NULL,
                         col = TRUE,
                         palette = "Spectral",
                         width.bar = 0.9,
                         cex.angle = 0,
                         cex.hjust = 0.5,
                         legend.position = "bottom",
                         cex = 12,
                         fontfam = "sans",
                         na.rm = TRUE,
                         verbose = TRUE){

  cl = match.call()
  datac = select(.data, ..., Y = !!enquo(resp)) %>%
     group_by(...) %>%
    summarise(N = n(),
              mean_var = mean(Y, na.rm = na.rm),
              sd = sd(Y, na.rm = na.rm),
              se = sd / sqrt(n()),
              ci = se * qt(level/2 + .5, n()-1))
  nam = names(select(.data, ...))
  if(length(nam)>1){
    names(datac) = c("x", "y", "N", "mean_var", "sd", "se", "ci")
  } else {
    names(datac) = c("x", "N", "mean_var", "sd", "se", "ci")  
  }
  if (is.null(ylab) == T){
    ylab = cl$resp
  }else {ylab = ylab}
  
  if(invert == FALSE){
  if (is.null(xlab) == T){
    xlab = nam[1]
  }else {xlab = xlab}
  } else{
    if (is.null(xlab) == T){
      xlab = nam[2]
    }else {xlab = xlab}  
    
  }
  
  pd = ggplot2::position_dodge(width.bar)
  if(length(nam)>1){
  if(invert == FALSE){
p = ggplot2::ggplot(data=datac, aes(x=x, y=mean_var, fill=y))+
    geom_bar(aes(fill = y), colour = "black", stat="identity", position=position_dodge(), width = width.bar)+
    scale_fill_brewer(type = "qualitative", palette = palette)
  } else{
p = ggplot2::ggplot(data=datac, aes(x=y, y=mean_var, fill=x))+
    geom_bar(aes(fill = x), colour = "black", stat="identity", position=position_dodge(),width = width.bar)
    scale_fill_brewer(type = "qualitative", palette = palette)
  }
  } else{
        p = ggplot2::ggplot(data=datac, aes(x=x, y=mean_var))+
        geom_bar(stat="identity", position=position_dodge(), width = width.bar)
  }
  
  
 p = p + ggplot2::theme_bw()
 
 if(col == FALSE){
p = p + scale_fill_grey(start = 0, end = .9)
} else{p = p}
   
 if (errorbar == TRUE){
 if(stat.erbar == "ci"){
 p = p + geom_errorbar(aes(ymin=mean_var-ci, ymax=mean_var+ci), width= width.erbar, position=pd)
 }

 if(stat.erbar == "sd"){
   p = p + geom_errorbar(aes(ymin=mean_var-sd, ymax=mean_var+sd), width =  width.erbar, position=pd)
 }
 
 if(stat.erbar == "se"){
   p = p + geom_errorbar(aes(ymin=mean_var-se, ymax=mean_var+se), width =  width.erbar, position=pd)
 }
 }

  p = p + ggplot2::theme(axis.ticks.length = unit(.2, "cm"),
                   axis.text = element_text(size = cex, family = fontfam, colour = "black"),
                   axis.title = element_text(size = cex,  family = fontfam, colour = "black"),
                   axis.text.x = element_text(angle = cex.angle, hjust = cex.hjust, size = cex, colour = "black"),
                   axis.ticks = element_line(colour = "black"),
                   plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
                   legend.title = element_blank(),
                   legend.position = legend.position,
                   legend.text = element_text(size = cex, family = fontfam),
                   panel.border = element_rect(colour = "black", fill=NA, size=1),
                   panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
                   panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())+
    ggplot2::labs(y = ylab, x = xlab)
  if(verbose == TRUE){
  print(datac)
  }
 return(p)
  
}
