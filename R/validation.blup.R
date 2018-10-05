validation.blup = function(data,
                           resp,
                           gen,
                           env,
                           rep,
                           nboot,
                           nrepval,
                           progbar = TRUE){

RMSPDres  = data.frame(RMSPD = matrix(".",nboot,1))
for (n in c(1,1:ncol(RMSPDres))) {
  RMSPDres[,n] = as.numeric(RMSPDres[,n])
}
Y = eval(substitute(resp), eval(data))
GEN = factor(eval(substitute(gen), eval(data)))
ENV = factor(eval(substitute(env), eval(data)))
REP = factor(eval(substitute(rep), eval(data)))
data = data.frame(cbind(ENV, GEN, REP, Y))
data$ID = as.numeric(rownames(data))
Nbloc = length(unique(REP))
Nenv = length(unique(ENV))

if (nrepval !=  Nbloc - 1){
  stop("The number replications used for validation must be equal to total number of replications -1 (In this case ", (Nbloc-1),").")
} else{

if (progbar  ==  TRUE){
pb = winProgressBar(title = "the model is being built, please, wait.",
                  min = 1, max = nboot, width = 570)
}
for (b in 1:nboot) {
  temp = data.frame(matrix(".",0,ncol(data)))
  for(n in c(4:5)) {
    temp[,n] = as.numeric(temp[,n])
  }
  names(temp) = names(data)
  actualenv = 0
  for (K in 1:Nenv){
    X = sample(1:10000,1)
    set.seed(X)
    X2 = sample(1:Nbloc,nrepval , replace = F)
    names = factor(data$ENV)
    names = levels(names)[actualenv + 1]
    actualenv = actualenv + 1
    temp2 = dplyr::filter(data, ENV  ==  names)
    modeling = temp2 %>%
    dplyr::group_by(GEN) %>%
    dplyr::filter(REP %in% c(X2))
    modeling = as.data.frame(modeling)
    modeling = rbind(temp, modeling)
    temp = modeling
  }
rownames(modeling ) = modeling$ID
testing = suppressWarnings(dplyr::anti_join(data, modeling, by = c("ENV", "GEN", "REP", "Y", "ID")))
testing = testing[order(testing[,1], testing[,2], testing[,3] ),]
ENV = as.factor(modeling$ENV)
GEN = as.factor(modeling$GEN)
BLOCO = as.factor(modeling$REP)
Resp = as.numeric(modeling$Y)
ovmean = mean(Resp)
model = suppressWarnings(suppressMessages(lme4::lmer(Resp~  BLOCO%in%ENV + (1|GEN) + (1|ENV)  + (1|GEN:ENV))))
bups = suppressWarnings(suppressMessages(lme4::ranef(model)))
blupGEN = bups$GEN
blupGEN = mutate(blupGEN,
                 GEN = factor(rownames(blupGEN)))
blupGEN = blupGEN %>%
  dplyr::select(GEN, everything())
colnames(blupGEN) = c("GEN", "BLUPg")
blupGEN = mutate(blupGEN,
               Predicted = BLUPg + ovmean)
blupGEN = blupGEN[order(-blupGEN[,3]),]
blupGEN = mutate(blupGEN,
               Rank = rank(-blupGEN[,3]))
blupGEN = blupGEN %>%
  dplyr::select(Rank, everything())
blupINT = bups$`GEN:ENV`
blups = data.frame(Names = rownames(blupINT))
blups = data.frame(do.call('rbind', strsplit(as.character(blups$Names),':',fixed = TRUE)))
blups = blups %>%
  dplyr::select(-X1, everything())
blups = cbind(blups, blupINT)
names(blups) = c("Code", "GEN", "BLUPge")
blups = blups[gtools::mixedorder(blups[,1]),]
raw = data.frame(ENV, GEN, Resp)
MEDIAS = tapply(raw[, 3], raw[, c(1, 2)], mean)
xx = rownames(MEDIAS)
yy = colnames(MEDIAS)
fila = length(xx)
col = length(yy)
total = fila * col
x = character(length = total)
y = character(length = total)
z = numeric(length = total)
k = 0
for (i in 1:fila) {
  for (j in 1:col) {
    k = k + 1
    x[k] = xx[i]
    y[k] = yy[j]
    z[k] = MEDIAS[i, j]
  }
}
MEDIAS = data.frame(ENV = x, GEN = y, Resp = z)
OUTMED = by(MEDIAS[, 3], MEDIAS[, c(2, 1)], function(x) sum(x,
                                                             na.rm = TRUE))
MGEN = data.frame(type = "GEN", Code = rownames(OUTMED),  Resp = apply(OUTMED, 1, mean))
MENV = data.frame(type = "ENV", Code = colnames(OUTMED), Resp = apply(OUTMED, 2, mean))
selectioNenv = suppressMessages(dplyr::left_join(blups, blupGEN  %>%
                                            dplyr::select(GEN, BLUPg)))
selectioNenv = suppressMessages(mutate(selectioNenv,
                      gge = BLUPge + BLUPg,
                      Predicted = BLUPge + BLUPg +
                      left_join(blups, MENV %>% dplyr::select(Code, Resp))$Resp))
names(selectioNenv) = c("ENV", "GEN", "BLUPge", "BLUPg", "BLUPg+ge", "Predicted")
validation  = mutate(selectioNenv,
                      testing = testing$Y,
                      error = Predicted - testing)
RMSPD = sqrt(sum(validation$error^2)/length(validation$error))
RMSPDres[,1][b] = RMSPD
if (progbar  ==  TRUE){
ProcdAtua = b
setWinProgressBar(pb, b, title = paste("Estimating BLUPs for ",ProcdAtua,
                                     " of ",nboot," total validation datasets",
                                     "-",round(b/nboot*100,3),"% Concluded -"))
 }
}
RMSPDres = dplyr::mutate(RMSPDres,
                        MODEL = "BLUP")
RMSPDres = RMSPDres %>%
          dplyr::select(MODEL, everything())
RMSPDmean = plyr::ddply(RMSPDres, .(MODEL), summarize, mean = mean(RMSPD))
if (progbar  ==  TRUE){
close(pb)
#utils::winDialog(type = "ok", "Validation sucessful! Check the results in R environment")
  }
return(structure(list(RMSPD = RMSPDres,
              RMSPDmean = RMSPDmean),
              class = "validation.blup"))
 }
}

