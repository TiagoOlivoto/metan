
plot.WAASBYratio = function(data,
                            type,
                            export = FALSE,
                            file.type = "pdf",
                            width = 6,
                            height = 5,
                            labSize = 1,
                            margins = c(5, 4),
                            y.lab = "Genotypes",
                            key.lab = "Genotype ranking",
                            resolution = 300
){
  if(type==1){


    if(export == F|FALSE) {
      Rowv= data$hetdata %>% dist %>% hclust %>% as.dendrogram %>%
        set("branches_k_color", k = 4) %>% set("branches_lwd", 2)%>%
        sort(type = "nodes")
      Colv=data$hetdata %>% t %>% dist %>% hclust %>% as.dendrogram %>%
        set("branches_k_color", k = 3) %>% set("branches_lwd", 2)%>%
        sort(type = "nodes")
      colfunc <- colorRampPalette(c("green","yellow", "red"))
      gplots::heatmap.2(data$hetdata, dendrogram = "both",
                        col=colfunc(nrow(subset(data$WAASY, Code="GEN"))),
                        Rowv = Rowv,
                        Colv = Colv,
                        density.info=("none"),
                        offsetRow = -0.4,
                        offsetCol = -0.4,
                        adjCol = c(1,0.5),
                        cexCol = labSize,
                        cexRow = labSize,
                        trace = "none",
                        key.title = "",
                        key.xlab = key.lab,
                        key.ylab = "",
                        key.par=list(mgp=c(1.5, 0.5, 0),
                                     mar=c(3, 0.2, 0.2, 0.2)),
                        lhei = c(1,6),
                        xlab = "Number of axis",
                        ylab = y.lab,
                        margins = margins)

    } else

      if(file.type=="pdf"){
        pdf("Heat Ranks PCA.pdf",width=width, height=height)
        Rowv= data$hetdata %>% dist %>% hclust %>% as.dendrogram %>%
          set("branches_k_color", k = 4) %>% set("branches_lwd", 2)%>%
          sort(type = "nodes")
        Colv=data$hetdata %>% t %>% dist %>% hclust %>% as.dendrogram %>%
          set("branches_k_color", k = 3) %>% set("branches_lwd", 2)%>%
          sort(type = "nodes")
        colfunc <- colorRampPalette(c("green","yellow", "red"))
        gplots::heatmap.2(data$hetdata, dendrogram = "both",
                          col=colfunc(nrow(subset(data$WAASY, Code="GEN"))),
                          Rowv = Rowv,
                          Colv = Colv,
                          density.info=("none"),
                          offsetRow = -0.4,
                          offsetCol = -0.4,
                          adjCol = c(1,0.5),
                          cexCol = labSize,
                          cexRow = labSize,
                          trace = "none",
                          key.title = "",
                          key.xlab = key.lab,
                          key.ylab = "",
                          key.par=list(mgp=c(1.5, 0.5, 0),
                                       mar=c(3, 0.2, 0.2, 0.2)),
                          lhei = c(1,6),
                          xlab = "Number of axis",
                          ylab = y.lab,
                          margins = margins)
        dev.off()
      }

    if (file.type=="tiff"){
      tiff(filename="Heat Ranks PCA.tiff",width=width, height=height, units = "in", compression = "lzw", res=resolution)
      Rowv= data$hetdata %>% dist %>% hclust %>% as.dendrogram %>%
        set("branches_k_color", k = 4) %>% set("branches_lwd", 2)%>%
        sort(type = "nodes")
      Colv=data$hetdata %>% t %>% dist %>% hclust %>% as.dendrogram %>%
        set("branches_k_color", k = 3) %>% set("branches_lwd", 2)%>%
        sort(type = "nodes")
      colfunc <- colorRampPalette(c("green","yellow", "red"))
      gplots::heatmap.2(data$hetdata, dendrogram = "both",
                        col=colfunc(nrow(subset(data$WAASY, Code="GEN"))),
                        Rowv = Rowv,
                        Colv = Colv,
                        density.info=("none"),
                        offsetRow = -0.4,
                        offsetCol = -0.4,
                        adjCol = c(1,0.5),
                        cexCol = labSize,
                        cexRow = labSize,
                        trace = "none",
                        key.title = "",
                        key.xlab = key.lab,
                        key.ylab = "",
                        key.par=list(mgp=c(1.5, 0.5, 0),
                                     mar=c(3, 0.2, 0.2, 0.2)),
                        lhei = c(1,6),
                        xlab = "Number of axis",
                        ylab = y.lab,
                        margins = margins)
      dev.off()
    }

  }


  if(type==2){


    if(export == F|FALSE) {
      Rowv=data$hetcomb %>% dist %>% hclust %>% as.dendrogram %>%
        set("branches_k_color", k = 4) %>% set("branches_lwd", 2)%>%
        rotate_DendSer(ser_weight = dist(data$hetcomb))%>%
        sort(type = "nodes")
      colfunc <- colorRampPalette(c("green","yellow", "red"))
      gplots::heatmap.2(data$hetcomb, dendrogram = "row",
                        col=colfunc(nrow(subset(data$WAASY, Code="GEN"))),
                        Rowv = Rowv,
                        Colv = F,
                        density.info=("none"),
                        offsetRow = -0.4,
                        offsetCol = -0.4,
                        adjCol = c(1,0.5),
                        trace = "none",
                        key.title = "",
                        key.xlab = key.lab,
                        key.ylab = "",
                        key.par=list(mgp=c(1.5, 0.5, 0),
                                     mar=c(3, 0.2, 0.2, 0.2)),
                        lhei = c(1,6),
                        xlab = "WAASB/GY ratio",
                        ylab = y.lab,
                        cexCol = labSize,
                        cexRow = labSize,
                        margins = margins)

    } else

      if(file.type=="pdf"){
        pdf("Heat map Ranks WAAS-GY.pdf",width=width, height=height)
        Rowv=data$hetcomb %>% dist %>% hclust %>% as.dendrogram %>%
          set("branches_k_color", k = 4) %>% set("branches_lwd", 2)%>%
          rotate_DendSer(ser_weight = dist(data$hetcomb))%>%
          sort(type = "nodes")
        colfunc <- colorRampPalette(c("green","yellow", "red"))
        gplots::heatmap.2(data$hetcomb, dendrogram = "row",
                          col=colfunc(nrow(subset(data$WAASY, Code="GEN"))),
                          Rowv = Rowv,
                          Colv = F,
                          density.info=("none"),
                          offsetRow = -0.4,
                          offsetCol = -0.4,
                          adjCol = c(1,0.5),
                          trace = "none",
                          key.title = "",
                          key.xlab = key.lab,
                          key.ylab = "",
                          key.par=list(mgp=c(1.5, 0.5, 0),
                                       mar=c(3, 0.2, 0.2, 0.2)),
                          lhei = c(1,6),
                          xlab = "WAASB/GY ratio",
                          ylab = y.lab,
                          cexCol = labSize,
                          cexRow = labSize,
                          margins = margins)
        dev.off()
      }

    if (file.type=="tiff"){
      tiff(filename="Heat map Ranks WAAS-GY.tiff",width=width, height=height, units = "in", compression = "lzw", res=resolution)
      Rowv=data$hetcomb %>% dist %>% hclust %>% as.dendrogram %>%
        set("branches_k_color", k = 4) %>% set("branches_lwd", 2)%>%
        rotate_DendSer(ser_weight = dist(data$hetcomb))%>%
        sort(type = "nodes")
      colfunc <- colorRampPalette(c("green","yellow", "red"))
      gplots::heatmap.2(data$hetcomb, dendrogram = "row",
                        col=colfunc(nrow(subset(data$WAASY, Code="GEN"))),
                        Rowv = Rowv,
                        Colv = F,
                        density.info=("none"),
                        offsetRow = -0.4,
                        offsetCol = -0.4,
                        adjCol = c(1,0.5),
                        trace = "none",
                        key.title = "",
                        key.xlab = key.lab,
                        key.ylab = "",
                        key.par=list(mgp=c(1.5, 0.5, 0),
                                     mar=c(3, 0.2, 0.2, 0.2)),
                        lhei = c(1,6),
                        xlab = "WAASB/GY ratio",
                        ylab = y.lab,
                        cexCol = labSize,
                        cexRow = labSize,
                        margins = margins)
      dev.off()
    }

  }

}

