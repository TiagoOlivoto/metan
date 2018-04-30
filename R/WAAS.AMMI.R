WAASB.AMMI=function(data, resp, p.valuePC = 0.05, naxis = NULL, weight.response, weight.WAASB){

  Y = data[paste(resp)]
  data = as.data.frame(data[,1:3])
  data = cbind(data, Y)
  names(data) = c("ENV", "GEN", "REP", "Y")


ENV = as.factor(data$ENV)
GEN = as.factor(data$GEN)
REP = as.factor(data$REP)
Y = as.numeric(data$Y)
Nenv = length(unique(ENV))
Ngen = length(unique(GEN))
minimo = min(Nenv, Ngen)-1


  model = agricolae::AMMI(ENV, GEN, REP, Y)
  PC=model$analysis
  PC = PC %>%
    select(-percent, -acum, everything())
  anova=model$ANOVA
  anova = cbind(Percent = ".", anova)
  anova = anova %>%
    select(-Percent, everything())
  anova = cbind(Accumul = ".", anova)
  anova = anova %>%
    select(-Accumul, everything())
  sum=as.data.frame(anova[nrow(anova),])
  sum$Df=sum(anova$Df)
  sum$`Sum Sq`=sum(anova$`Sum Sq`)
  sum$`Mean Sq`=sum$`Sum Sq`/sum$Df
  rownames(sum)="Total"
  ERRO = anova[nrow(anova),]
  ERRO=rbind(ERRO, sum)
  anova = anova[-nrow(anova),]
  names(PC) =colnames(anova)
  anova2 = suppressWarnings(suppressMessages(rbind(anova, PC)))
  anova=plyr::rbind.fill(anova, PC)
  rownames(anova)=rownames(anova2)
  anova=rbind(anova, ERRO)
  anova$`Pr(>F)`=format(anova$`Pr(>F)`, scipen=0, digits=5, scientific=FALSE)
  anova$Percent=format(anova$Percent, scipen=0, digits=3, scientific=FALSE)
  anova$Accumul=format(anova$Accumul, scipen=0, digits=3, scientific=FALSE)
  MeansGxE=model$means


  ################################# Save de Scores ####################################
  Escores = model$biplot
  Escores = cbind(Code = row.names(Escores), Escores)
  Escores = Escores %>%
    select(type, everything())

  if(is.null(naxis)){
  SigPC1=nrow(PC[which(PC[,5]<p.valuePC),])
  } else{
    SigPC1 = naxis
  }

  if(SigPC1 > minimo){
    stop("The number of axis to be used must be lesser than or equal to ", minimo, " [min(GEN-1;ENV-1)]")
  } else{


  Pesos=model$analysis[1]
  Pesos=as.data.frame(Pesos[c(1:SigPC1), ])
  colnames(Pesos)="Percent"
  WAASB=Escores
  WAASBAbs=Escores
  for (i in 4:ncol(WAASB)){
    WAASB[,i] <- abs(WAASB[i])
  }
  t_WAASB = transpose(WAASB)
  colnames(t_WAASB)  = rownames(WAASB)
  rownames(t_WAASB)  = colnames(WAASB)
    t_WAASB=t_WAASB[-c(1, 2, 3), ]
  t_WAASB=t_WAASB[c(1:SigPC1), ]
    t_WAASB=cbind(t_WAASB, Pesos)
  for (i in 1:ncol(t_WAASB)){
    t_WAASB[,i] <- as.numeric(as.character(t_WAASB[,i]))
  }
  Ponderado=transpose(as.data.frame(sapply(t_WAASB[ ,-ncol(t_WAASB)], weighted.mean,  w = t_WAASB$Percent)))
  rownames(Ponderado)=c("WAASB")
  t_WAASB=subset(t_WAASB, select = -Percent)
  colnames(Ponderado)=colnames(t_WAASB)
  t_WAASB=rbind(t_WAASB, Ponderado)
  t_WAASB2=transpose(t_WAASB)
  colnames(t_WAASB2)  = rownames(t_WAASB)
  rownames(t_WAASB2)  = colnames(t_WAASB)
  WAASBAbs=cbind(WAASBAbs, subset(t_WAASB2, select = WAASB))
  WAASBAbs2=subset(WAASBAbs, type=="ENV")
  WAASBAbs2$PctResp = (WAASBAbs2$Y / max(WAASBAbs2$Y))*100
  WAASBAbs2$PctWAASB = (100-WAASBAbs2$WAASB / min(WAASBAbs2$WAASB))
  WAASBAbs3=subset(WAASBAbs, type=="GEN")
  WAASBAbs3$PctResp = (WAASBAbs3$Y / max(WAASBAbs3$Y))*100
  WAASBAbs3$PctWAASB = (100-WAASBAbs3$WAASB / min(WAASBAbs3$WAASB))
  WAASBAbs=rbind(WAASBAbs3, WAASBAbs2)
  WAASBAbs=setDT(WAASBAbs)[, OrResp:=rank(-Y), by = type][]
  WAASBAbs=setDT(WAASBAbs)[, OrWAASB:=rank(WAASB), by = type][]
  WAASBAbs=setDT(WAASBAbs)[, OrPC1:=rank(abs(PC1)), by = type][]
  WAASBAbs$PesRes=as.vector(weight.response)
  WAASBAbs$PesWAASB=as.vector(weight.WAASB)
  for (i in 1:nrow(WAASBAbs)){
    WAASBAbs$WAASBY[i] = (WAASBAbs$PctResp[i]*WAASBAbs$PesRes[i]+WAASBAbs$PctWAASB[i]*WAASBAbs$PesWAASB[i])/
      sum(WAASBAbs$PesRes[i]+WAASBAbs$PesWAASB[i])
  }
  WAASB=data.table::setDT(WAASBAbs)[, OrWAASBY:=rank(-WAASBY), by = type][]

  MinENV = WAASBAbs2[which(WAASBAbs2[,3] <=min(WAASBAbs2$Y)),]
  MinENV = paste0("Environment ", MinENV$Code , " (", round(MinENV$Y,4), ") ")

  MaxENV = WAASBAbs2[which(WAASBAbs2[,3] >=max(WAASBAbs2$Y)),]
  MaxENV = paste0("Environment ", MaxENV$Code , " (", round(MaxENV$Y,4), ") ")

  MinGEN = WAASBAbs3[which(WAASBAbs3[,3] <=min(WAASBAbs3$Y)),]
  MinGEN = paste0("Genotype ", MinGEN$Code , " (", round(MinGEN$Y,4), ") ")


  MaxGEN = WAASBAbs3[which(WAASBAbs3[,3] >=max(WAASBAbs3$Y)),]
  MaxGEN = paste0("Genotype ", MaxGEN$Code , " (", round(MaxGEN$Y,4), ") ")


  mean = round(mean(MeansGxE$Y),4)

  min = MeansGxE[which(MeansGxE[,3] <=min(MeansGxE$Y)),]
  min = paste0(round(min$Y,4), " (Genotype ", min$GEN , " in ", min$ENV, " )")

  max = MeansGxE[which(MeansGxE[,3] >=max(MeansGxE$Y)),]
  max = paste0(round(max$Y,4), " (Genotype ", max$GEN , " in ", max$ENV, " )")


  PCA = PC[,4:7]

  return(list(model = WAASB, MeansGxE = MeansGxE, PCA = PCA, anova = anova, WgtResponse=weight.response,
              WgtWAASB=weight.WAASB, Ngen=Ngen, Nenv = Nenv, OVmean=mean, Min=min, Max=max, MinENV=MinENV, MaxENV=MaxENV,
              MinGEN=MinGEN, MaxGEN=MaxGEN ))
}
}



