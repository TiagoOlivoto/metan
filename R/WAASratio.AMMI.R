WAASratio.AMMI = function(data,
                           resp,
                           p.valuePC = 0.05,
                           increment = 5,
                           saveWAASY = 50,
                           progbar = TRUE){
PesoWAAS = 100
PesoResp = 0
Y = data[paste(resp)]
data = as.data.frame(data[,1:3])
data = cbind(data, Y)
names(data) = c("ENV", "GEN", "REP", "Y")
Ngen = nrow(plyr::count(data, "GEN"))
Nenv = nrow(plyr::count(data, "ENV"))
ncomb = (100/increment) + 1

test = saveWAASY %% increment == 0
if (test == FALSE){
  stop("The argument 'saveWAASY = ", saveWAASY,"' must be divisible by 'increment' (", increment, "). Please, consider changing the values.")
} else{


CombWAASY = data.frame(type = matrix(".",(Ngen + Nenv),1))
WAASY.Values = list()
model = with(data, agricolae::AMMI(ENV, GEN, REP, Y))
PC = model$analysis
PC = PC %>%
     dplyr::select(-percent, -acum, everything())
anova = model$ANOVA
anova = cbind(Percent = ".", anova)
anova = anova %>%
        dplyr::select(-Percent, everything())
anova = cbind(Accumul = ".", anova)
anova = anova %>%
        dplyr::select(-Accumul, everything())
sum = as.data.frame(anova[nrow(anova),])
sum$Df = sum(anova$Df)
sum$`Sum Sq` = sum(anova$`Sum Sq`)
sum$`Mean Sq` = sum$`Sum Sq`/sum$Df
rownames(sum) = "Total"
ERRO = anova[nrow(anova),]
ERRO = rbind(ERRO, sum)
anova = anova[-nrow(anova),]
names(PC) = colnames(anova)
anova2 = suppressWarnings(suppressMessages(rbind(anova, PC)))
anova = plyr::rbind.fill(anova, PC)
rownames(anova) = rownames(anova2)
anova = rbind(anova, ERRO)
anova$`Pr(>F)` = format(anova$`Pr(>F)`, scipen = 0, digits = 5, scientific = FALSE)
anova$Percent = format(anova$Percent, scipen = 0, digits = 3, scientific = FALSE)
anova$Accumul = format(anova$Accumul, scipen = 0, digits = 3, scientific = FALSE)
MeansGxE = model$means
totalcomb = ncomb*nrow(PC)
initial = 0
if (progbar == TRUE){
pb = winProgressBar(title = "the model is being built, please, wait.",
min = 1, max = totalcomb, width = 570)
}
for (k in 1:ncomb){
Escores = model$biplot
Escores = cbind(Code = row.names(Escores), Escores)
Escores = Escores %>%
          dplyr::select(type, everything())
SigPC1 = nrow(PC[which(PC[,5] < p.valuePC),])
Pesos = model$analysis[1]
Pesos = as.data.frame(Pesos[c(1:SigPC1), ])
colnames(Pesos) = "Percent"
WAAS = Escores
WAASAbs = Escores
for (i in 4:ncol(WAAS)){
WAAS[,i] <- abs(WAAS[i])
}
t_WAAS = data.table::transpose(WAAS)
colnames(t_WAAS) = rownames(WAAS)
rownames(t_WAAS) = colnames(WAAS)
t_WAAS = t_WAAS[-c(1, 2, 3), ]
t_WAAS = t_WAAS[c(1:SigPC1), ]
t_WAAS = cbind(t_WAAS, Pesos)
for (i in 1:ncol(t_WAAS)){
t_WAAS[,i] <- as.numeric(as.character(t_WAAS[,i]))
}
Ponderado = data.table::transpose(as.data.frame(sapply(t_WAAS[ ,-ncol(t_WAAS)], weighted.mean,w = t_WAAS$Percent)))
rownames(Ponderado) = c("WAAS")
t_WAAS = subset(t_WAAS, select = -Percent)
colnames(Ponderado) = colnames(t_WAAS)
t_WAAS = rbind(t_WAAS, Ponderado)
t_WAAS2 = data.table::transpose(t_WAAS)
colnames(t_WAAS2) = rownames(t_WAAS)
rownames(t_WAAS2) = colnames(t_WAAS)
WAASAbs = cbind(WAASAbs, subset(t_WAAS2, select = WAAS))
WAASAbs2 = subset(WAASAbs, type  ==  "ENV")
WAASAbs2$PctResp = (WAASAbs2$Y / max(WAASAbs2$Y))*100
WAASAbs2$PctWAAS = (100  -WAASAbs2$WAAS / min(WAASAbs2$WAAS))
WAASAbs3 = subset(WAASAbs, type  ==  "GEN")
WAASAbs3$PctResp = (WAASAbs3$Y / max(WAASAbs3$Y))*100
WAASAbs3$PctWAAS = (100 - WAASAbs3$WAAS / min(WAASAbs3$WAAS))
WAASAbs = rbind(WAASAbs3, WAASAbs2)
WAASAbs = data.table::setDT(WAASAbs)[, OrResp := rank(-Y), by = type][]
WAASAbs = data.table::setDT(WAASAbs)[, OrWAAS := rank(WAAS), by = type][]
WAASAbs = data.table::setDT(WAASAbs)[, OrPC1 := rank(abs(PC1)), by = type][]
WAASAbs$PesRes = as.vector(PesoResp)
WAASAbs$PesWAAS = as.vector(PesoWAAS)
WAASAbs = plyr::mutate(WAASAbs,
WAASY = (PctResp * PesRes + PctWAAS * PesWAAS)/(PesRes + PesWAAS))
WAASAbsInicial = data.table::setDT(WAASAbs)[, OrWAASY := rank(-WAASY), by = type][]
inicial = as.data.frame(WAASAbsInicial$OrWAAS)
colnames(inicial) = paste0(SigPC1,"PC")
SigPC2 = 1
for (j in 1:nrow(PC)){
Escores = model$biplot
Escores = cbind(Code = row.names(Escores), Escores)
Escores = Escores %>%
          dplyr::select(type, everything())
Pesos = model$analysis[1]
Pesos = as.data.frame(Pesos[c(1:SigPC2), ])
colnames(Pesos) = "Percent"
WAAS = Escores
WAASAbs = Escores
for (i in 4:ncol(WAAS)){
WAAS[,i] <- abs(WAAS[i])
}
t_WAAS = data.table::transpose(WAAS)
colnames(t_WAAS) = rownames(WAAS)
rownames(t_WAAS) = colnames(WAAS)
t_WAAS = t_WAAS[-c(1, 2, 3), ]
t_WAAS = t_WAAS[c(1:SigPC2), ]
t_WAAS = cbind(t_WAAS, Pesos)
for (i in 1:ncol(t_WAAS)){
t_WAAS[,i] <- as.numeric(as.character(t_WAAS[,i]))
}
Ponderado = data.table::transpose(as.data.frame(sapply(t_WAAS[ ,-ncol(t_WAAS)], weighted.mean,w = t_WAAS$Percent)))
rownames(Ponderado) = c("WAAS")
t_WAAS = subset(t_WAAS, select = -Percent)
colnames(Ponderado) = colnames(t_WAAS)
t_WAAS = rbind(t_WAAS, Ponderado)
t_WAAS2 = data.table::transpose(t_WAAS)
colnames(t_WAAS2) = rownames(t_WAAS)
rownames(t_WAAS2) = colnames(t_WAAS)
WAASAbs = cbind(WAASAbs, subset(t_WAAS2, select = WAAS))
WAASAbs2 = subset(WAASAbs, type  ==  "ENV")
WAASAbs2$PctResp = (WAASAbs2$Y / max(WAASAbs2$Y))*100
WAASAbs2$PctWAAS = (100 - WAASAbs2$WAAS / min(WAASAbs2$WAAS))
WAASAbs3 = subset(WAASAbs, type  ==  "GEN")
WAASAbs3$PctResp = (WAASAbs3$Y / max(WAASAbs3$Y))*100
WAASAbs3$PctWAAS = (100 - WAASAbs3$WAAS / min(WAASAbs3$WAAS))
WAASAbs = rbind(WAASAbs3, WAASAbs2)
WAASAbs = data.table::setDT(WAASAbs)[, OrResp := rank(-Y), by = type][]
WAASAbs = data.table::setDT(WAASAbs)[, OrWAAS := rank(WAAS), by = type][]
WAASAbs = data.table::setDT(WAASAbs)[, OrPC1 := rank(abs(PC1)), by = type][]
WAASAbs$PesRes = as.vector(PesoResp)
WAASAbs$PesWAAS = as.vector(PesoWAAS)
for (i in 1:nrow(WAASAbs)){
WAASAbs$WAASY[i] = (WAASAbs$PctResp[i]*WAASAbs$PesRes[i] + WAASAbs$PctWAAS[i]*WAASAbs$PesWAAS[i]) /
sum(WAASAbs$PesRes[i] + WAASAbs$PesWAAS[i])
}
WAASAbs = data.table::setDT(WAASAbs)[, OrWAASY := rank(-WAASY), by = type][]
results = as.data.frame(WAASAbs$OrWAAS)
names(results) = paste0(SigPC2,"PC")
final = cbind(results, inicial)
inicial = final
SigPC2 = SigPC2 + 1
ProcdAtua = j
initial = initial + 1
Sys.sleep(0.1)
if (progbar == TRUE){
setWinProgressBar(pb, initial, title = paste("Obtaining the Ranks considering",ProcdAtua,
 " of ",nrow(PC),"Principal components:",
 "|WAAS:",PesoWAAS,"% ",
 "GY:",PesoResp,"%|" ,
 "-",round(initial/totalcomb*100,2),"% Concluded -"))
}
}
initial = initial
WAAS = WAASAbsInicial
WAAS$type = ifelse(WAAS$type  ==  "GEN", "Genotype", "Environment")
CombWAASY[[sprintf("%.0f/%.0f",PesoWAAS,PesoResp)]] = WAASAbsInicial$OrWAASY
WAASY.Values[[paste("WAAS/GY",PesoWAAS, "/",PesoResp)]] = data.frame( WAAS)
PesoResp = PesoResp + increment
PesoWAAS = PesoWAAS - increment
if (PesoWAAS + increment  ==  saveWAASY){
genotypes = subset(WAAS, type  ==  "Genotype")
genotypes = genotypes[,c("Code", "PesRes", "PesWAAS", "WAASY")]
genotypes = genotypes[order(genotypes$WAASY), ]
genotypes$Code = factor(genotypes$Code, levels = genotypes$Code)
genotypes$Mean = ifelse(genotypes$WAASY < mean(genotypes$WAASY), "below", "above")}
}
Rank = final[,-(SigPC2)]
Names = WAASAbsInicial[,c("type", "Code", "OrResp", "OrPC1", "OrWAAS") ]
Rank = cbind(Names, Rank)
hetdata = as.data.frame(subset(Rank, type  ==  "GEN"))
rownames(hetdata) = hetdata$Code
hetdata = hetdata[,-c(1,2, 3, 4)]
colnames(hetdata)[1:1] = c("WAAS")
hetdata = as.matrix(hetdata)
CombWAASY = CombWAASY[,-1]
CombWAASY = cbind(Names[,c(1,2)], CombWAASY)
hetcomb = as.data.frame(subset(CombWAASY, type  ==  "GEN"))
rownames(hetcomb) = hetcomb$Code
hetcomb = hetcomb[,-c(1,2)]
hetcomb = as.matrix(hetcomb)
CorrRank = Rank[,c("type", "Code", "OrResp","OrWAAS")]
CorrCWA = CombWAASY[,-c(1,2)]
CorcombWAASY = as.data.frame(cbind(CorrRank, CorrCWA))
CorcombWAASY = subset(CorcombWAASY, type  ==  "GEN")
rownames(CorcombWAASY) = CorcombWAASY$Code
CorcombWAASY = CorcombWAASY[,-c(1)]
PC1 = Pesos[1,1]
PC2 = Pesos[2,1]
mean = mean(WAAS$Y)
if (progbar == TRUE){
close(pb)
utils::winDialog(type = "ok", "Procedure suceful! Check the results in R environment")
}
return(structure(list(anova = anova,
                      PC = PC,
                      MeansGxE = MeansGxE,
                      WAAS = WAAS,
                      WAASxGY = WAASY.Values,
                      WAASY = genotypes,
                      hetcomb = hetcomb,
                      hetdata = hetdata,
                      Ranks = Rank),
                      class = "WAASratio.AMMI"))
}
}
