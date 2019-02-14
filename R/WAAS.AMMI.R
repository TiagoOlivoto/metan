WAAS.AMMI = function(data,
                     resp,
                     gen,
                     env,
                     rep,
                     mresp = NULL,
                     wresp = NULL,
                     prob = 0.05,
                     naxis = NULL,
                     verbose = TRUE){

datain = data
GEN = factor(eval(substitute(gen), eval(datain)))
ENV = factor(eval(substitute(env), eval(datain)))
REP = factor(eval(substitute(rep), eval(datain)))
listres = list()
d = match.call()
nvar = as.numeric(ifelse(length(d$resp)>1, length(d$resp) -1, length(d$resp)))

if (!is.null(naxis)){
  if(length(d$resp) > 1) {
if(length(naxis) != length(d$resp) -1){
    stop("The argument 'naxix' must length of ", nvar,
         ", the same number of variables in object 'resp'.")
}
} else{
    if(length(naxis) != length(d$resp)){
      stop("The argument 'naxix' must length of ", nvar,
           ", the same number of variables in object 'resp'.")
}
}
}

if (is.null(mresp)) {
  mresp = replicate(nvar, 100)
  minresp = 100 - mresp
} else {
  if (length(mresp) != nvar) {
    stop("The length of the numeric vector 'mresp' must be equal the number of variables in argument 'resp'")
  }
  if (sum(mresp == 100) + sum(mresp == 0) != nvar) {
    stop("The values of the numeric vector 'mresp' must be 0 or 100.")
  }
  mresp = mresp
  minresp = 100 - mresp
}

if (is.null(wresp)) {
  PesoResp = replicate(nvar, 50)
  PesoWAASB = 100 - PesoResp
} else{
  if (length(wresp) != nvar) {
    stop("The length of the numeric vector 'wresp' must be equal the number of variables in argument 'resp'")
  }
  if (min(wresp) < 0 | max(wresp) > 100) {
    stop("The range of the numeric vector 'wresp' must be equal between 0 and 100.")
  }
  PesoResp = wresp
  PesoWAASB = 100 - PesoResp
}
vin = 0
for (var in 2:length(d$resp)){
if(length(d$resp)>1){
Y = eval(substitute(resp)[[var]], eval(datain))
} else{
Y = eval(substitute(resp), eval(datain))
}
data = data.frame(ENV, GEN, REP, Y)
Nenv = length(unique(ENV))
Ngen = length(unique(GEN))
minimo = min(Nenv, Ngen) - 1
vin = vin + 1
if (minimo < 2) {
  cat("\nWarning. The analysis is not possible.")
  cat("\nThe number of environments and number of genotypes must be greater than 2\n")
}

temp = data.frame(matrix(0,length(unique(data$ENV)),12))
actualenv = 0
names(temp) = c("ENV", "Mean", "MSblock", "MSgen", "MSres", "Fcal(Blo)", "Pr>F(Blo)","Fcal(Gen)", "Pr>F(Gen)", "CV(%)", "h2", "AS")
for (i in 1:length(unique(data$ENV))){
  envnam = levels(data$ENV)[actualenv + 1]
  data2 = subset(data, ENV == paste0(envnam))
  anova = suppressWarnings(anova(aov(Y ~ GEN + REP, data = data2)))
  h2 = (anova[1, 3] - anova[3, 3])/anova[1, 3]
  temp[i,1] = paste(envnam)
  temp[i,2] = mean(data2$Y)
  temp[i,3] = anova[2, 3]
  temp[i,4] = anova[1, 3]
  temp[i,5] = anova[3, 3]
  temp[i,6] = anova[2, 4]
  temp[i,7] = anova[2, 5]
  temp[i,8] = anova[1, 4]
  temp[i,9] = anova[1, 5]
  temp[i,10] = sqrt(anova[3, 3])/mean(data2$Y)*100
  temp[i,11] = ifelse(h2 < 0, 0, h2)
  temp[i,12] = ifelse(h2 < 0, 0, sqrt(h2))
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
  SigPC1 = nrow(PC[which(PC[,5]<prob),])
  } else{
    if(nvar == 1){
    SigPC1 = naxis
    } else{
    SigPC1 = naxis[vin]
  }
}
  if (SigPC1 > minimo){
    stop("The number of axis to be used must be lesser than or equal to ", minimo, " [min(GEN-1;ENV-1)]")
  } else{
Pesos = model$analysis[1]
Pesos = as.data.frame(Pesos[c(1:SigPC1), ])
colnames(Pesos) = "Percent"
WAAS = Escores
WAASAbs = Escores
WAAS[, 4:ncol(WAAS)] = lapply(WAAS[,4:ncol(WAAS)], abs)

t_WAAS = data.frame(t(WAAS))
colnames(t_WAAS)  = rownames(WAAS)
rownames(t_WAAS)  = colnames(WAAS)
t_WAAS = t_WAAS[-c(1, 2, 3), ]
t_WAAS = t_WAAS[c(1:SigPC1), ]
t_WAAS = cbind(t_WAAS, Pesos)
  for (i in 1:ncol(t_WAAS)){
t_WAAS[,i] <- as.numeric(as.character(t_WAAS[,i]))
  }
Ponderado = t(as.data.frame(sapply(t_WAAS[ ,-ncol(t_WAAS)], weighted.mean,  w = t_WAAS$Percent)))
rownames(Ponderado) = c("WAAS")
t_WAAS = subset(t_WAAS, select = -Percent)
colnames(Ponderado) = colnames(t_WAAS)
t_WAAS = rbind(t_WAAS, Ponderado)
t_WAAS2 = data.frame(t(t_WAAS))
WAASAbs = cbind(WAASAbs, subset(t_WAAS2, select = WAAS))
WAASAbs2 = subset(WAASAbs, type == "ENV")
if (nvar > 1){
WAASAbs2$PctResp = resca(WAASAbs2$Y, new_min = minresp[vin], new_max = mresp[vin])
} else{
WAASAbs2$PctResp = resca(WAASAbs2$Y, new_min = minresp, new_max = mresp)
}
WAASAbs2$PctWAAS = resca(WAASAbs2$WAAS, new_min =  100, new_max = 0)
WAASAbs3 = subset(WAASAbs, type == "GEN")
if (nvar > 1){
WAASAbs3$PctResp = resca(WAASAbs3$Y, new_min = minresp[vin], new_max = mresp[vin])
} else {
WAASAbs3$PctResp = resca(WAASAbs3$Y, new_min = minresp, new_max = mresp)
}
WAASAbs3$PctWAAS = resca(WAASAbs3$WAAS, 100, 0)
WAASAbs = rbind(WAASAbs3, WAASAbs2) %>%
          dplyr::group_by(type) %>%
          dplyr::mutate(OrResp = rank(-Y),
                        OrWAAS = rank(WAAS),
                        OrPC1 = rank(abs(PC1)))

if (nvar > 1) {
WAASAbs$PesRes = as.vector(PesoResp)[vin]
WAASAbs$PesWAAS = as.vector(PesoWAASB)[vin]
} else {
WAASAbs$PesRes = as.vector(65)
WAASAbs$PesWAAS = as.vector(35)
}

WAAS = WAASAbs %>% dplyr::mutate(WAASY = ((PctResp*PesRes) + (PctWAAS * PesWAAS)) / (PesRes + PesWAAS)) %>%
         dplyr::group_by(type) %>%
         dplyr::mutate(OrWAASY = rank(-WAASY))

MinENV = WAASAbs2[head(which(WAASAbs2[,3] <=  min(WAASAbs2$Y)), n = 1),]
MinENV = paste0("Environment ", MinENV$Code , " (", round(MinENV$Y,4), ") ")
MaxENV = WAASAbs2[head(which(WAASAbs2[,3] >=  max(WAASAbs2$Y)), n = 1),]
MaxENV = paste0("Environment ", MaxENV$Code , " (", round(MaxENV$Y,4), ") ")
MinGEN = WAASAbs3[head(which(WAASAbs3[,3] <=  min(WAASAbs3$Y)), n = 1),]
MinGEN = paste0("Genotype ", MinGEN$Code , " (", round(MinGEN$Y,4), ") ")
MaxGEN = WAASAbs3[head(which(WAASAbs3[,3] >=  max(WAASAbs3$Y)), n = 1),]
MaxGEN = paste0("Genotype ", MaxGEN$Code , " (", round(MaxGEN$Y,4), ") ")
mean = round(mean(MeansGxE$Y),4)
min = MeansGxE[head(which(MeansGxE[,3] <=  min(MeansGxE$Y)), n = 1),]
min = paste0(round(min$Y,4), " (Genotype ", min$GEN , " in ", min$ENV, " )")
max = MeansGxE[head(which(MeansGxE[, 3] >= max(MeansGxE$Y)), n = 1),]
max = paste0(round(max$Y,4), " (Genotype ", max$GEN , " in ", max$ENV, " )")
PCA = PC[,4:7]
Details = list(Ngen = Ngen,
               Nenv = Nenv,
               OVmean = mean,
               Min = min,
               Max = max,
               MinENV = MinENV,
               MaxENV = MaxENV,
               MinGEN = MinGEN,
               MaxGEN = MaxGEN,
               SigPC = SigPC1)
Details = do.call(rbind.data.frame, Details)
names(Details) = "Values"
Details = plyr::mutate(Details,
                         Parameters = c("Ngen", "Nenv", "OVmean", "Min",
                                        "Max", "MinENV", "MaxENV", "MinGEN", "MaxGEN", "SigPC"))
Details = Details %>%
          dplyr::select(Parameters, everything())

temp = structure(list(individual = individual,
                      model = WAAS,
                      MeansGxE = MeansGxE,
                      PCA = PCA,
                      anova = anova,
                      Details = Details,
                      residuals = model$residuals,
                      probint = model$probint),
                 class = "WAAS.AMMI")

if(length(d$resp)>1){
listres[[paste(d$resp[var])]] = temp
if (verbose == TRUE){
cat("Evaluating variable", paste(d$resp[var]), round((var - 1)/(length(d$resp) - 1)*100, 1), "%", "\n")
}
} else{
listres[[paste(d$resp)]] = temp
}
  }
}
if (verbose == T){
  if(length(which(unlist(lapply(listres, function(x){
    x[["probint"]]
  })) > prob)) > 0){
  cat("------------------------------------------------------------\n")
  cat("Variables with nonsignificant GxE interaction\n")
  cat(names(which(unlist(lapply(listres, function(x){
    x[["probint"]]
  })) > prob)), "\n")
  cat("------------------------------------------------------------\n")
  }
  cat("Done!\n")
}
invisible(structure(listres, class = "WAAS.AMMI"))

}


