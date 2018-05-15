plot.WAASBYratio = function(x,
                            type,
                            export = FALSE,
                            file.type = "pdf",
                            file.name = NULL,
                            width = 6,
                            height = 5,
                            size.lab = 1,
                            margins = c(5, 4),
                            y.lab = NULL,
                            x.lab = NULL,
                            key.lab = "Genotype ranking",
                            resolution = 300,
                            ...){
  data = x
  if (type == 1){

    if (is.null(x.lab)){
      x.lab = "Number of axes"}
    else x.lab  = x.lab

    if (is.null(y.lab)){
      y.lab = "Genotypes"}
    else y.lab = y.lab


    if (export  ==  F|FALSE) {

      Rowv =  data$hetdata %>% dist %>% hclust %>% as.dendrogram %>%
        dendextend::set("branches_k_color", k = 4) %>% dendextend::set("branches_lwd", 2) %>%
        sort(type = "nodes")
      Colv = data$hetdata %>% t %>% dist %>% hclust %>% as.dendrogram %>%
        dendextend::set("branches_k_color", k = 3) %>% dendextend::set("branches_lwd", 2) %>%
        sort(type = "nodes")
      colfunc = grDevices::colorRampPalette(c("green","yellow", "red"))
      gplots::heatmap.2(data$hetdata, dendrogram = "both",
                        col = colfunc(nrow(subset(data$WAASY, Code = "GEN"))),
                        Rowv = Rowv,
                        Colv = Colv,
                        density.info = ("none"),
                        offsetRow = -0.4,
                        offsetCol = -0.4,
                        adjCol = c(1,0.5),
                        cexCol = size.lab,
                        cexRow = size.lab,
                        trace = "none",
                        key.title = "",
                        key.xlab = key.lab,
                        key.ylab = "",
                        key.par = list(mgp = c(1.5, 0.5, 0),
                                     mar = c(3, 0.2, 0.2, 0.2)),
                        lhei = c(1,6),
                        xlab = x.lab,
                        ylab = y.lab,
                        margins = margins)

    } else

      if (file.type == "pdf"){
        if (is.null(file.name)){
        pdf("Heat Ranks PCA.pdf",width = width, height = height)
        } else
        pdf(paste0(file.name, ".pdf"), width = width, height = height)

        Rowv =  data$hetdata %>% dist %>% hclust %>% as.dendrogram %>%
          dendextend::set("branches_k_color", k = 4) %>% dendextend::set("branches_lwd", 2) %>%
          sort(type = "nodes")
        Colv = data$hetdata %>% t %>% dist %>% hclust %>% as.dendrogram %>%
          dendextend::set("branches_k_color", k = 3) %>% dendextend::set("branches_lwd", 2) %>%
          sort(type = "nodes")
        colfunc = grDevices::colorRampPalette(c("green","yellow", "red"))
        gplots::heatmap.2(data$hetdata, dendrogram = "both",
                          col = colfunc(nrow(subset(data$WAASY, Code = "GEN"))),
                          Rowv = Rowv,
                          Colv = Colv,
                          density.info = ("none"),
                          offsetRow = -0.4,
                          offsetCol = -0.4,
                          adjCol = c(1,0.5),
                          cexCol = size.lab,
                          cexRow = size.lab,
                          trace = "none",
                          key.title = "",
                          key.xlab = key.lab,
                          key.ylab = "",
                          key.par = list(mgp = c(1.5, 0.5, 0),
                                       mar = c(3, 0.2, 0.2, 0.2)),
                          lhei = c(1,6),
                          xlab = x.lab,
                          ylab = y.lab,
                          margins = margins)
        dev.off()
      }

    if (file.type == "tiff"){
      if (is.null(file.name)){
      tiff(filename = "Heat Ranks PCA.tiff",width = width, height = height, units = "in", compression = "lzw", res = resolution)
      } else
      tiff(filename = paste0(file.name, ".tiff"), width = width, height = height, units = "in", compression = "lzw", res = resolution)
        Rowv =  data$hetdata %>% dist %>% hclust %>% as.dendrogram %>%
        dendextend::set("branches_k_color", k = 4) %>% dendextend::set("branches_lwd", 2) %>%
        sort(type = "nodes")
      Colv = data$hetdata %>% t %>% dist %>% hclust %>% as.dendrogram %>%
        dendextend::set("branches_k_color", k = 3) %>% dendextend::set("branches_lwd", 2) %>%
        sort(type = "nodes")
      colfunc = grDevices::colorRampPalette(c("green","yellow", "red"))
      gplots::heatmap.2(data$hetdata, dendrogram = "both",
                        col = colfunc(nrow(subset(data$WAASY, Code = "GEN"))),
                        Rowv = Rowv,
                        Colv = Colv,
                        density.info = ("none"),
                        offsetRow = -0.4,
                        offsetCol = -0.4,
                        adjCol = c(1,0.5),
                        cexCol = size.lab,
                        cexRow = size.lab,
                        trace = "none",
                        key.title = "",
                        key.xlab = key.lab,
                        key.ylab = "",
                        key.par = list(mgp = c(1.5, 0.5, 0),
                                     mar = c(3, 0.2, 0.2, 0.2)),
                        lhei = c(1,6),
                        xlab = x.lab,
                        ylab = y.lab,
                        margins = margins)
      dev.off()
    }

  }

  if (type == 2){

    if (is.null(x.lab)){
      x.lab = "WAASB/GY ratio"}
    else x.lab  = x.lab

    if (is.null(y.lab)){
      y.lab = "Genotypes"}
    else y.lab = y.lab

    if (export  ==  F|FALSE) {
      Rowv = data$hetcomb %>% dist %>% hclust %>% as.dendrogram %>%
        dendextend::set("branches_k_color", k = 4) %>% dendextend::set("branches_lwd", 2)%>%
        dendextend::rotate_DendSer(ser_weight = dist(data$hetcomb))%>%
        sort(type = "nodes")
      colfunc = grDevices::colorRampPalette(c("green","yellow", "red"))
      gplots::heatmap.2(data$hetcomb, dendrogram = "row",
                        col = colfunc(nrow(subset(data$WAASY, Code = "GEN"))),
                        Rowv = Rowv,
                        Colv = F,
                        density.info = ("none"),
                        offsetRow = -0.4,
                        offsetCol = -0.4,
                        adjCol = c(1,0.5),
                        trace = "none",
                        key.title = "",
                        key.xlab = key.lab,
                        key.ylab = "",
                        key.par = list(mgp = c(1.5, 0.5, 0),
                                     mar = c(3, 0.2, 0.2, 0.2)),
                        lhei = c(1,6),
                        xlab = x.lab,
                        ylab = y.lab,
                        cexCol = size.lab,
                        cexRow = size.lab,
                        margins = margins)

    } else

      if (file.type == "pdf"){

        if (is.null(file.name)){
          pdf("Heat map Ranks WAAS-GY.pdf",width = width, height = height)
        } else
          pdf(paste0(file.name, ".pdf"), width = width, height = height)

        Rowv = data$hetcomb %>% dist %>% hclust %>% as.dendrogram %>%
          dendextend::set("branches_k_color", k = 4) %>% dendextend::set("branches_lwd", 2)%>%
          dendextend::rotate_DendSer(ser_weight = dist(data$hetcomb))%>%
          sort(type = "nodes")
        colfunc = grDevices::colorRampPalette(c("green","yellow", "red"))
        gplots::heatmap.2(data$hetcomb, dendrogram = "row",
                          col = colfunc(nrow(subset(data$WAASY, Code = "GEN"))),
                          Rowv = Rowv,
                          Colv = F,
                          density.info = ("none"),
                          offsetRow = -0.4,
                          offsetCol = -0.4,
                          adjCol = c(1,0.5),
                          trace = "none",
                          key.title = "",
                          key.xlab = key.lab,
                          key.ylab = "",
                          key.par = list(mgp = c(1.5, 0.5, 0),
                                       mar = c(3, 0.2, 0.2, 0.2)),
                          lhei = c(1,6),
                          xlab = x.lab,
                          ylab = y.lab,
                          cexCol = size.lab,
                          cexRow = size.lab,
                          margins = margins)
        dev.off()
      }

    if (file.type == "tiff"){
      if (is.null(file.name)){
        tiff(filename = "Heat map Ranks WAAS-GY.tiff",width = width, height = height, units = "in", compression = "lzw", res = resolution)
      } else
        tiff(filename = paste0(file.name, ".tiff"), width = width, height = height, units = "in", compression = "lzw", res = resolution)

      Rowv = data$hetcomb %>% dist %>% hclust %>% as.dendrogram %>%
        dendextend::set("branches_k_color", k = 4) %>% dendextend::set("branches_lwd", 2)%>%
        dendextend::rotate_DendSer(ser_weight = dist(data$hetcomb))%>%
        sort(type = "nodes")
      colfunc = grDevices::colorRampPalette(c("green","yellow", "red"))
      gplots::heatmap.2(data$hetcomb, dendrogram = "row",
                        col = colfunc(nrow(subset(data$WAASY, Code = "GEN"))),
                        Rowv = Rowv,
                        Colv = F,
                        density.info = ("none"),
                        offsetRow = -0.4,
                        offsetCol = -0.4,
                        adjCol = c(1,0.5),
                        trace = "none",
                        key.title = "",
                        key.xlab = key.lab,
                        key.ylab = "",
                        key.par = list(mgp = c(1.5, 0.5, 0),
                                     mar = c(3, 0.2, 0.2, 0.2)),
                        lhei = c(1,6),
                        xlab = x.lab,
                        ylab = y.lab,
                        cexCol = size.lab,
                        cexRow = size.lab,
                        margins = margins)
      dev.off()
    }

  }

}

