corr_plot = function(.data,
                     ... = NULL,
                     upper = "corr",
                     lower = "scatter",
                     axis.labels = FALSE,
                     diag = TRUE,
                     col.diag = "gray",
                     alpha.diag = 1,
                     col.up.panel = "gray",
                     col.lw.panel = "gray",
                     col.dia.panel = "gray",
                     prob = 0.05,
                     col.sign = "green",
                     alpha.sign = 0.15,
                     lab.position = "tr",
                     progress = NULL,
                     smooth = FALSE,
                     col.smooth = "red",
                     size.smooth = 0.3,
                     confint = TRUE,
                     size.point = 1,
                     shape.point = 19,
                     alpha.point = 0.7,
                     fill.point = NULL,
                     col.point = "black",
                     minsize = 2,
                     maxsize = 3,
                     pan.spacing = 0.15,
                     digits = 2,
                     export = FALSE,
                     file.type = "pdf",
                     file.name = NULL,
                     width = 8,
                     height = 7,
                     resolution = 300){
  if (!lab.position %in% c("tr", "tl", "br", "bl")){
    stop("The argument 'lab.position' must be one of the 'tr', 'tl', 'br', or 'bl'.")
  }

  if(!is.null(upper)){
    if (!upper %in% c("corr", "scatter", NULL)){
      stop("The argument 'upper' must be one of the 'corr', 'scatter' or 'NULL'.")
    }
  }
  if(!is.null(lower)){
    if (!lower %in% c("corr", "scatter", NULL)){
      stop("The argument 'lower' must be one of the 'corr', 'scatter' or 'NULL'.")
    }
  }
  if (missing(...)){
    data = .data[ , unlist(lapply(.data, is.numeric))]
    if (sum(lapply(.data, is.factor)==TRUE)>0){
      message("The factors ", paste0(collapse = " ", names(.data[ , unlist(lapply(.data, is.factor)) ])),
              " where excluded to perform the analysis. Only numeric variables were used. ")
    }
  }
  if (!missing(...)){
    data= dplyr::select(.data, ...)
  }
  w = c(21:25)
  if(is.null(fill.point)==TRUE && any(w == shape.point)){
    stop("If 'shape.point' is a value between 21 and 25, you must provide a color to fill the shape using the argument 'fill.point.'")
  }
  my_custom_cor = function(data, mapping, color = I("black"), sizeRange = c(minsize, maxsize), ...) {
    x = GGally::eval_data_col(data, mapping$x)
    y = GGally::eval_data_col(data, mapping$y)
    ct = cor.test(x,y)
    sig = symnum(
      ct$p.value, corr = FALSE, na = FALSE,
      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
      symbols = c("***", "**", "*", ".", " ")
    )
    r = unname(ct$estimate)
    rt = format(r, digits = digits)[1]
    cex = max(sizeRange)
    percent_of_range = function(percent, range) {
      percent * diff(range) + min(range, na.rm = TRUE)
    }
    GGally::ggally_text(
      label = as.character(rt),
      mapping = aes(),
      xP = 0.5, yP = 0.5,
      size = I(percent_of_range(cex * abs(r), sizeRange)),
      color = color,
      ...
    ) +
      geom_text(aes_string(x = 0.8, y = 0.8),
                label = sig,
                size = I(cex),
                color = color,
                ...) +
      theme_classic() +
      theme(panel.background = ggplot2::element_rect(color = col.up.panel),
            axis.line = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank())
  }
  my_custom_smooth = function(data, mapping, ...) {
    x = GGally::eval_data_col(data, mapping$x)
    y = GGally::eval_data_col(data, mapping$y)
    ct <- cor.test(x,y)
    r <- unname(ct$p.value)
    rt <- format(r, digits = digits)[1]
    tt <- as.character(rt)
    p = ggplot2::ggplot(data = data, mapping = mapping)
    if(is.null(fill.point)==FALSE){
      p = p + geom_point(color = I(col.point), fill = fill.point, shape = shape.point, size = size.point, alpha = alpha.point)
    } else{
      p = p + geom_point(color = I(col.point), shape = shape.point, size = size.point, alpha = alpha.point)
    }
    p = p + theme_classic() +
           theme(panel.background = ggplot2::element_rect(fill = "white", color = col.lw.panel))
    if(smooth == TRUE){
      p  = p + geom_smooth(method = "lm", se = confint, size = size.smooth, color = col.smooth, ...)
    }
    if (r < prob) {
      p = p + theme(
        panel.background = ggplot2::element_rect(fill=ggplot2::alpha(col.sign, alpha.sign)))
    }
    p
  }
  ggally_mysmooth = function(data, mapping, ...){
    ggplot2::ggplot(data = data, mapping=mapping) +
      ggplot2::geom_density(fill = ggplot2::alpha(col.diag, alpha.diag))+
      ggplot2::theme_classic() +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill=ggplot2::alpha('white', 1), color = col.dia.panel))
  }
  if (diag == TRUE){
    diag = list(continuous = ggally_mysmooth)
  } else {diag = NULL}
  if (!is.null(upper)){
    if (upper == "corr"){
      upper = list(continuous = my_custom_cor)
    }
    if (upper == "scatter"){
      upper = list(continuous = my_custom_smooth)
    }
  }
  if (is.null(upper)){
    upper = NULL
  }
  if (!is.null(lower)){
    if (lower == "corr"){
      lower = list(continuous = my_custom_cor)
    }
    if (lower == "scatter"){
      lower = list(continuous = my_custom_smooth)
    }
  }
  if (is.null(lower)){
    lower = NULL
  }
  if(lab.position == "tr"){
    switch = NULL
  }
  if (lab.position == "tl"){
    switch = "y"
  }
  if (lab.position == "br"){
    switch = "x"
  }
  if (lab.position == "bl"){
    switch = "both"
  }
  if (axis.labels == TRUE){
    axis.labels = "show"
  } else {axis.labels = "none"}
  p1 = GGally::ggpairs(data,
                       upper = upper,
                       lower = lower,
                       switch = switch,
                       diag = diag,
                       progress = progress,
                       axisLabels = axis.labels)
  ggplot2::theme_set(ggplot2::theme_gray() +
                       ggplot2::theme(panel.spacing = grid::unit(pan.spacing,"lines")))

  if (export  ==  FALSE) {
    return(p1)
  } else
    if (file.type == "pdf"){
      if (is.null(file.name)){
        pdf("Scatterplot Correlation.pdf",width = width, height = height)
      } else
        pdf(paste0(file.name, ".pdf"), width = width, height = height)
      print(p1)
      dev.off()
    }

  if (file.type == "tiff"){
    if (is.null(file.name)){
      tiff(filename = "Scatterplot Correlation.tiff",width = width, height = height, units = "in", compression = "lzw", res = resolution)
    } else
      tiff(filename = paste0(file.name, ".tiff"), width = width, height = height, units = "in", compression = "lzw", res = resolution)
    print(p1)
    dev.off()
  }
}
