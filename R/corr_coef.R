#' Reorder a correlation matrix
#' @description
#' `r badge('stable')`
#'
#' Reorder the correlation matrix according to the correlation coefficient by
#' using  hclust for hierarchical clustering order. This is useful to identify
#' the hidden pattern in the matrix.
#'
#' @param x The correlation matrix
#'
#' @return The ordered correlation matrix
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(metan)
#' cor_mat <- corr_coef(data_ge2, PH, EH, CD, CL, ED, NKR)
#' cor_mat$cor
#' reorder_cormat(cor_mat$cor)
#' }
#'
reorder_cormat <- function(x){
  if(!is.matrix(x) | nrow(x) != ncol(x)){
    stop("object 'x' must be a square matrix.", call. = FALSE)
  }
  hc <- hclust(as.dist((1 - x) / 2))
  x <- x[hc$order, hc$order]
  return(x)
}
NULL



#' Computes Pearson's correlation matrix with p-values
#'
#' Computes Pearson's correlation matrix with p-values
#'
#'
#' @param data The data set.
#' @param ... Variables to use in the correlation. If no variable is informed
#'   all the numeric variables from `data` are used.
#' @param verbose Logical argument. If `verbose = FALSE` the code is run
#'  silently.
#' @return A list with the correlation coefficients and p-values
#' @importFrom dplyr between
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#' \donttest{
#' library(metan)
#'
#' # All numeric variables
#' all <- corr_coef(data_ge2)
#'
#' # Select variables
#' sel <- corr_coef(data_ge2, EP, EL, CD, CL)
#' print(sel)
#' }
#'
corr_coef <- function(data, ..., verbose = TRUE){
  if(missing(...)){
    x <- select_if(data, is.numeric)
    if(has_na(x)){
      x <- remove_rows_na(x)
      has_text_in_num(x)
    }
  }
  if(!missing(...)){
    x <- select(data, ...) %>%
      select_numeric_cols()
    if(has_na(x)){
      x <- remove_rows_na(x)
      has_text_in_num(x)
    }
  }
  if(has_na(data) ==  TRUE){
    x <- remove_rows_na(x, verbose = verbose)
  }
  apply_r <- function(A, FUN, ...) {
    mapply(function(a, B)
      lapply(B, function(x) FUN(a, x, ...)),
      a = A,
      MoreArgs = list(B = A))
  }
  pval <- apply(apply_r(x, cor.test), 1:2, function(x) x[[1]]$p.value)
  r <- apply(apply_r(x, cor), 1:2, function(x) x[[1]])
  return(structure(list(cor = r, pval = pval), class = "corr_coef"))
}
NULL





#' Print an object of class corr_coef
#'
#' Print the `corr_coef` object in two ways. By default, the results
#' are shown in the R console. The results can also be exported to the directory
#' into a *.txt file.
#'
#'
#' @param x An object of class `corr_coef`
#' @param export A logical argument. If `TRUE`, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if `export = TRUE`
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   [tibble::formatting()] for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print corr_coef
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' sel <- corr_coef(data_ge2, EP, EL, CD, CL)
#' print(sel)
#' }
print.corr_coef <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (!class(x) == "corr_coef") {
    stop("The object must be of class 'corr_coef'")
  }
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "corr_coef print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  cat("---------------------------------------------------------------------------\n")
  cat("Pearson's correlation coefficient\n")
  cat("---------------------------------------------------------------------------\n")
  print(x$cor, digits = digits)
  cat("---------------------------------------------------------------------------\n")
  cat("p-values for the correlation coefficients\n")
  cat("---------------------------------------------------------------------------\n")
  print(x$pval, digits = digits)
  cat("\n\n\n")
  if (export == TRUE) {
    sink()
  }
}
NULL






#' Create a correlation heat map
#'
#' Create a correlation heat map for object of class `corr_coef`
#'
#' @param x The data set.
#' @param type The type of heat map to produce. Either `lower` (default) to
#'   produce a lower triangle heat map or `upper` to produce an upper
#'   triangular heat map.
#' @param diag Plot diagonal elements? Defaults to `FALSE`.
#' @param reorder Reorder the correlation matrix to identify the hidden pattern?
#'   Defaults to `TRUE`.
#' @param signif How to show significant correlations. If `"stars"` is used
#'   (default), stars are used showing the significance at 0.05 ("*"), 0.01
#'   ("**") and 0.001 ("***") probability error. If `signif = "pval"`, then
#'   the p-values are shown.
#' @param caption Logical. If `TRUE` (Default) includes a caption with the
#'   significance meaning for stars.
#' @param digits.cor,digits.pval The significant digits to show for correlations
#'   and p-values, respectively.
#' @param col.low,col.mid,col.high The color for the low (-1), mid(0) and high
#'   (1) points in the color key. Defaults to `blue`, `white`, and
#'   `red`, respectively.
#' @param lab.x.position,lab.y.position The position of the x and y axis label.
#'   Defaults to `"bottom"` and `"right"` if `type = "lower"` or
#'   `"top"` and `"left"` if `type = "upper"`.
#' @param legend.position The legend position in the plot.
#' @param legend.title The title of the color key. Defaults to `"Pearson's
#'   Correlation"`.
#' @param size.text.cor The size of the text for correlation values. Defaults to 3.
#' @param size.text.signif The size of the text for significance values (stars or p-values). Defaults to 3.
#' @param size.text.lab The size of the text for labels. Defaults to 10.
#' @param ... Currently not used.
#' @method plot corr_coef
#' @return An object of class `gg, ggplot`
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' # All numeric variables
#' all <- corr_coef(data_ge2)
#' plot(all)
#' plot(all, reorder = FALSE)
#'
#' # Select variables
#' sel <- corr_coef(data_ge2, EP, EL, CD, CL)
#' plot(sel,
#'      type = "upper",
#'      reorder = FALSE,
#'      size.text.lab = 14,
#'      size.text.plot = 5)
#' }

plot.corr_coef <- function(x,
                           type = "lower",
                           diag = FALSE,
                           reorder = TRUE,
                           signif = "stars",
                           caption = TRUE,
                           digits.cor = 2,
                           digits.pval = 3,
                           col.low = "blue",
                           col.mid = "white",
                           col.high = "red",
                           lab.x.position = NULL,
                           lab.y.position = NULL,
                           legend.position = NULL,
                           legend.title = "Pearson's\nCorrelation",
                           size.text.cor = 3,
                           size.text.signif = 3,
                           size.text.lab = 10,
                           ...){
  if(!signif %in% c("stars", "pval")){
    stop("The argument 'signif' must be one of 'stars' or 'pval'.")
  }
  pval <- x$pval
  correl <- x$cor
  if(reorder == TRUE){
    correl <- reorder_cormat(correl)
    pval <- pval[rownames(correl), colnames(correl)]
  }
  if(type == "lower"){
    pval <- make_lower_tri(pval)
    correl <- make_lower_tri(correl)
  }
  if(type == "upper"){
    pval <- make_upper_tri(pval)
    correl <- make_upper_tri(correl)
  }
  if(signif == "stars"){
    bind_data <- expand.grid(dimnames(correl)) %>%
      mutate(cor = as.vector(correl),
             pval = as.vector(pval),
             pval_star <- case_when(
               pval < 0.001 ~ "***",
               between(pval, 0.001, 0.01) ~ "**",
               between(pval, 0.01, 0.05) ~ "*",
               pval == NA ~ "",
               FALSE  ~ ""
             )) %>%
      set_names("v1", "v2", "cor", "pval", "pval_star") %>%
      dplyr::filter(!is.na(pval))
  } else{
    bind_data <-
      expand.grid(dimnames(correl)) %>%
      mutate(cor = as.vector(correl),
             pval_star = as.vector(signif(pval, digits = digits.pval))) %>%
      set_names("v1", "v2", "cor", "pval_star") %>%
      dplyr::filter(!is.na(pval_star))
  }

  if(diag == FALSE){
    bind_data <- dplyr::filter(bind_data, !v1 == v2)
  }
  if(type == "lower"){
    lab.y.position <- ifelse(missing(lab.y.position), "right", lab.y.position)
    lab.x.position <- ifelse(missing(lab.x.position), "bottom", lab.x.position)
  }
  if(type == "upper"){
    lab.y.position <- ifelse(missing(lab.y.position), "left", lab.y.position)
    lab.x.position <- ifelse(missing(lab.x.position), "top", lab.x.position)
  }
  axis.text.x.hjust <- ifelse(lab.x.position == "bottom", 1, 0)
  if(missing(legend.position) & type == "lower"){
    legend.position <- c(0.2, 0.85)
  }
  if(missing(legend.position) & type == "upper"){
    legend.position <- c(0.8, 0.17)
  }
  if(!missing(legend.position)){
    legend.position <- legend.position
  }
  p <-
    ggplot(bind_data, aes(v1, v2, fill = cor)) +
    geom_tile(aes(fill = cor),
              colour = "white") +
    geom_text(aes(label = round(cor, digits.cor)),
              vjust = 0,
              size = size.text.cor) +
    geom_text(aes(label = pval_star),
              vjust = 1.5,
              size = size.text.signif)+
    scale_fill_gradient2(low = col.low,
                         high = col.high,
                         mid = col.mid,
                         midpoint = 0,
                         limit = c(-1, 1),
                         space = "Lab",
                         name = legend.title)+
    theme_bw() +
    xlab(NULL) +
    ylab(NULL) +
    scale_y_discrete(expand = expansion(mult = c(0,0)),
                     position = lab.y.position)+
    scale_x_discrete(expand = expansion(mult = c(0,0)),
                     position = lab.x.position) +
    coord_fixed() +
    theme(axis.ticks = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          legend.background = element_blank(),
          legend.position = legend.position,
          legend.direction = "horizontal",
          axis.text = element_text(color = "black", size = size.text.lab),
          axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = axis.text.x.hjust))+
    guides(fill = guide_colorbar(barwidth = 6,
                                 barheight = 1,
                                 title.position = "top",
                                 title.hjust = 0.5))
  if(signif == "stars" & caption == TRUE){
    p <-
      p + labs(caption = c("* p < 0.05; ** p < 0.01; and *** p < 0.001"))
  }
  suppressWarnings(return(p))
}
