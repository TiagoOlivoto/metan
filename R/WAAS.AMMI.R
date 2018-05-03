WAAS.AMMI = function(data, resp, p.valuePC = 0.05, naxis = NULL, weight.response, weight.WAAS){

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
  WAAS=Escores
  WAASAbs=Escores
  for (i in 4:ncol(WAAS)){
    WAAS[,i] <- abs(WAAS[i])
  }
  t_WAAS = transpose(WAAS)
  colnames(t_WAAS)  = rownames(WAAS)
  rownames(t_WAAS)  = colnames(WAAS)
    t_WAAS=t_WAAS[-c(1, 2, 3), ]
  t_WAAS=t_WAAS[c(1:SigPC1), ]
    t_WAAS=cbind(t_WAAS, Pesos)
  for (i in 1:ncol(t_WAAS)){
    t_WAAS[,i] <- as.numeric(as.character(t_WAAS[,i]))
  }
  Ponderado=transpose(as.data.frame(sapply(t_WAAS[ ,-ncol(t_WAAS)], weighted.mean,  w = t_WAAS$Percent)))
  rownames(Ponderado)=c("WAAS")
  t_WAAS=subset(t_WAAS, select = -Percent)
  colnames(Ponderado)=colnames(t_WAAS)
  t_WAAS=rbind(t_WAAS, Ponderado)
  t_WAAS2=transpose(t_WAAS)
  colnames(t_WAAS2)  = rownames(t_WAAS)
  rownames(t_WAAS2)  = colnames(t_WAAS)
  WAASAbs=cbind(WAASAbs, subset(t_WAAS2, select = WAAS))
  WAASAbs2=subset(WAASAbs, type=="ENV")
  WAASAbs2$PctResp = (WAASAbs2$Y / max(WAASAbs2$Y))*100
  WAASAbs2$PctWAAS = (100-WAASAbs2$WAAS / min(WAASAbs2$WAAS))
  WAASAbs3=subset(WAASAbs, type=="GEN")
  WAASAbs3$PctResp = (WAASAbs3$Y / max(WAASAbs3$Y))*100
  WAASAbs3$PctWAAS = (100-WAASAbs3$WAAS / min(WAASAbs3$WAAS))
  WAASAbs=rbind(WAASAbs3, WAASAbs2)
  WAASAbs=setDT(WAASAbs)[, OrResp:=rank(-Y), by = type][]
  WAASAbs=setDT(WAASAbs)[, OrWAAS:=rank(WAAS), by = type][]
  WAASAbs=setDT(WAASAbs)[, OrPC1:=rank(abs(PC1)), by = type][]
  WAASAbs$PesRes=as.vector(weight.response)
  WAASAbs$PesWAAS=as.vector(weight.WAAS)
  for (i in 1:nrow(WAASAbs)){
    WAASAbs$WAASY[i] = (WAASAbs$PctResp[i]*WAASAbs$PesRes[i]+WAASAbs$PctWAAS[i]*WAASAbs$PesWAAS[i])/
      sum(WAASAbs$PesRes[i]+WAASAbs$PesWAAS[i])
  }
  WAAS=data.table::setDT(WAASAbs)[, OrWAASY:=rank(-WAASY), by = type][]

  MinENV = WAASAbs2[which(WAASAbs2[,3] <=min(WAASAbs2$Y)),]
  MinENV = paste0("Environment ", MinENV$Code , " (", round(MinENV$Y,4), ") ")

  MaxENV = WAASAbs2[which(WAASAbs2[,3] >=max(WAASAbs2$Y)),]
  MaxENV = paste0("Environment ", MaxENV$Code , " (", round(MaxENV$Y,4), ") ")

  MinGEN = WAASAbs3[which(WAASAbs3[,3] <=min(WAASAbs3$Y)),]
  MinGEN = paste0("Genotype ", MinGEN$Code , " (", round(MinGEN$Y,4), ") ")


  MaxGEN = WAASAbs3[which(WAASAbs3[,3] >=max(WAASAbs3$Y)),]
  MaxGEN = paste0("Genotype ", MaxGEN$Code , " (", round(MaxGEN$Y,4), ") ")


  mean = round(mean(MeansGxE$Y),4)

  min = MeansGxE[which(MeansGxE[,3] <=min(MeansGxE$Y)),]
  min = paste0(round(min$Y,4), " (Genotype ", min$GEN , " in ", min$ENV, " )")

  max = MeansGxE[which(MeansGxE[,3] >=max(MeansGxE$Y)),]
  max = paste0(round(max$Y,4), " (Genotype ", max$GEN , " in ", max$ENV, " )")


  PCA = PC[,4:7]

  Details = list(WgtResponse=weight.response, WgtWAAS=weight.WAAS, Ngen=Ngen, Nenv = Nenv,
                 OVmean=mean, Min=min, Max=max, MinENV=MinENV, MaxENV=MaxENV, MinGEN=MinGEN, MaxGEN=MaxGEN)
  Details = do.call(rbind.data.frame, Details)
  names(Details) = "Values"
  rownames(Details) = c("WgtResponse", "WgtWAAS", "Ngen", "Nenv", "OVmean", "Min",
                          "Max", "MinENV", "MaxENV", "MinGEN", "MaxGEN")

  return(list(model = WAAS, MeansGxE = MeansGxE, PCA = PCA,
              anova = anova, Details = Details, object = "WAAS"))
}
}



