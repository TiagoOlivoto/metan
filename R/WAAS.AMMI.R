WAAS.AMMI = function(data,
                     resp,
                     gen,
                     env,
                     rep,
                     p.valuePC = 0.05,
                     naxis = NULL,
                     weight.response = 50,
                     weight.WAAS = 50){

Y = eval(substitute(resp), eval(data))
GEN = factor(eval(substitute(gen), eval(data)))
ENV = factor(eval(substitute(env), eval(data)))
REP = factor(eval(substitute(rep), eval(data)))
data = data.frame(ENV, GEN, REP, Y)
Nenv = length(unique(ENV))
Ngen = length(unique(GEN))
minimo = min(Nenv, Ngen) - 1

if (minimo < 2) {
  cat("\nWarning. The analysis is not possible.")
  cat("\nThe number of environments and number of genotypes must be greater than 2\n")
}

temp = data.frame(matrix(".",length(unique(data$ENV)),13))
actualenv = 0
for (n in c(1:13)) {
  temp[,n] = as.numeric(temp[,n])
}
names(temp) = c("ENV", "Mean", "MSblock", "MSgen", "MSres", "Fcal(Blo)", "Pr>F(Blo)","Fcal(Gen)", "Pr>F(Gen)", "CV(%)", "h2", "AS", "R2")
for (i in 1:length(unique(data$ENV))){
  envnam = levels(data$ENV)[actualenv + 1]
  data2 = subset(data, ENV == paste0(envnam))
  anova = anova(aov(Y ~ GEN + REP, data = data2))
  MSB = anova[2, 3]
  MSG = anova[1, 3]
  MSE = anova[3, 3]
  NR = length(unique(data2$REP))
  CV = sqrt(MSE)/mean(data2$Y)*100
  h2 = (MSG - MSE)/MSG
  AS = sqrt(h2)
  temp[i,1] = paste(envnam)
  temp[i,2] = mean(data2$Y)
  temp[i,3] = MSB
  temp[i,4] = MSG
  temp[i,5] = MSE
  temp[i,6] = anova[2, 4]
  temp[i,7] = anova[2, 5]
  temp[i,8] = anova[1, 4]
  temp[i,9] = anova[1, 5]
  temp[i,10] = CV
  temp[i,11] = h2
  temp[i,12] = AS
  temp[i,13] = 1/(2-AS^2)

  actualenv = actualenv + 1

}
MSEratio = max(temp$MSres) / min(temp$MSres)
individual = list(individual = temp,
                  MSEratio = MSEratio)

model = performs_ammi(ENV, GEN, REP, Y)
anova = model$ANOVA
PC = model$analysis
MeansGxE = model$means[,1:3]
Escores = model$biplot
Escores = cbind(Code = row.names(Escores), Escores)
Escores = Escores %>%
          dplyr::select(type, everything())

EscGEN = subset(Escores, type == "GEN")
names(EscGEN)[2] = "GEN"
names(EscGEN)[3] = "y"
EscENV = subset(Escores, type == "ENV")
names(EscENV)[2] = "ENV"
MeansGxE = suppressMessages(suppressWarnings(dplyr::mutate(MeansGxE,
                                        envPC1 = left_join(MeansGxE, EscENV %>% select(ENV, PC1))$PC1,
                                        genPC1 = left_join(MeansGxE, EscGEN %>% select(GEN, PC1))$PC1,
                                        nominal = left_join(MeansGxE, EscGEN %>% select(GEN, y))$y + genPC1*envPC1)))


  if (is.null(naxis)){
  SigPC1 = nrow(PC[which(PC[,5]<p.valuePC),])
  } else{
    SigPC1 = naxis
  }
  if (SigPC1 > minimo){
    stop("The number of axis to be used must be lesser than or equal to ", minimo, " [min(GEN-1;ENV-1)]")
  } else{
Pesos = model$analysis[1]
Pesos = as.data.frame(Pesos[c(1:SigPC1), ])
colnames(Pesos) = "Percent"
WAAS = Escores
WAASAbs = Escores
  for (i in 4:ncol(WAAS)){
    WAAS[,i] <- abs(WAAS[i])
  }
t_WAAS = data.table::transpose(WAAS)
colnames(t_WAAS)  = rownames(WAAS)
rownames(t_WAAS)  = colnames(WAAS)
t_WAAS = t_WAAS[-c(1, 2, 3), ]
t_WAAS = t_WAAS[c(1:SigPC1), ]
t_WAAS = cbind(t_WAAS, Pesos)
  for (i in 1:ncol(t_WAAS)){
t_WAAS[,i] <- as.numeric(as.character(t_WAAS[,i]))
}
Ponderado = data.table::transpose(as.data.frame(sapply(t_WAAS[ ,-ncol(t_WAAS)], weighted.mean,  w = t_WAAS$Percent)))
rownames(Ponderado) = c("WAAS")
t_WAAS = subset(t_WAAS, select = -Percent)
colnames(Ponderado) = colnames(t_WAAS)
t_WAAS = rbind(t_WAAS, Ponderado)
t_WAAS2 = data.table::transpose(t_WAAS)
colnames(t_WAAS2)  = rownames(t_WAAS)
rownames(t_WAAS2)  = colnames(t_WAAS)
WAASAbs = cbind(WAASAbs, subset(t_WAAS2, select = WAAS))
WAASAbs2 = subset(WAASAbs, type == "ENV")
WAASAbs2$PctResp = (WAASAbs2$Y / max(WAASAbs2$Y))*100
WAASAbs2$PctWAAS = (100-WAASAbs2$WAAS / min(WAASAbs2$WAAS))
WAASAbs3 = subset(WAASAbs, type == "GEN")
WAASAbs3$PctResp = (WAASAbs3$Y / max(WAASAbs3$Y))*100
WAASAbs3$PctWAAS = (100-WAASAbs3$WAAS / min(WAASAbs3$WAAS))
WAASAbs = rbind(WAASAbs3, WAASAbs2)
WAASAbs = data.table::setDT(WAASAbs)[, OrResp:= rank(-Y), by = type][]
WAASAbs = data.table::setDT(WAASAbs)[, OrWAAS:= rank(WAAS), by = type][]
WAASAbs = data.table::setDT(WAASAbs)[, OrPC1:= rank(abs(PC1)), by = type][]
WAASAbs$PesRes = as.vector(weight.response)
WAASAbs$PesWAAS = as.vector(weight.WAAS)
  for (i in 1:nrow(WAASAbs)){
    WAASAbs$WAASY[i] = (WAASAbs$PctResp[i]*WAASAbs$PesRes[i] + WAASAbs$PctWAAS[i]*WAASAbs$PesWAAS[i])/
      sum(WAASAbs$PesRes[i] + WAASAbs$PesWAAS[i])
  }
WAAS = data.table::setDT(WAASAbs)[, OrWAASY:= rank(-WAASY), by = type][]
MinENV = WAASAbs2[which(WAASAbs2[,3] <=  min(WAASAbs2$Y)),]
MinENV = paste0("Environment ", MinENV$Code , " (", round(MinENV$Y,4), ") ")
MaxENV = WAASAbs2[which(WAASAbs2[,3] >=  max(WAASAbs2$Y)),]
MaxENV = paste0("Environment ", MaxENV$Code , " (", round(MaxENV$Y,4), ") ")
MinGEN = WAASAbs3[which(WAASAbs3[,3] <=  min(WAASAbs3$Y)),]
MinGEN = paste0("Genotype ", MinGEN$Code , " (", round(MinGEN$Y,4), ") ")
MaxGEN = WAASAbs3[which(WAASAbs3[,3] >=  max(WAASAbs3$Y)),]
MaxGEN = paste0("Genotype ", MaxGEN$Code , " (", round(MaxGEN$Y,4), ") ")
mean = round(mean(MeansGxE$Y),4)
min = MeansGxE[which(MeansGxE[,3] <=  min(MeansGxE$Y)),]
min = paste0(round(min$Y,4), " (Genotype ", min$GEN , " in ", min$ENV, " )")
max = MeansGxE[which(MeansGxE[,3] >=  max(MeansGxE$Y)),]
max = paste0(round(max$Y,4), " (Genotype ", max$GEN , " in ", max$ENV, " )")
PCA = PC[,4:7]
Details = list(WgtResponse = weight.response,
               WgtWAAS = weight.WAAS,
               Ngen = Ngen,
               Nenv = Nenv,
               OVmean = mean,
               Min = min,
               Max = max,
               MinENV = MinENV,
               MaxENV = MaxENV,
               MinGEN = MinGEN,
               MaxGEN = MaxGEN)
Details = do.call(rbind.data.frame, Details)
names(Details) = "Values"
Details = plyr::mutate(Details,
                         Parameters = c("WgtResponse", "WgtWAAS", "Ngen", "Nenv", "OVmean", "Min",
                                        "Max", "MinENV", "MaxENV", "MinGEN", "MaxGEN"))
Details = Details %>%
          dplyr::select(Parameters, everything())
return(structure(list(individual = individual,
                      model = WAAS,
                      MeansGxE = MeansGxE,
                      PCA = PCA,
                      anova = anova,
                      Details = Details),
                 class = "WAAS.AMMI"))
  }
}



