WAASB = function(data,
                 resp,
                 gen,
                 env,
                 rep,
                 random = "gen",
                 prob = 0.95,
                 weight.response = 50,
                 weight.WAAS = 50){

  Y = eval(substitute(resp), eval(data))
  GEN = factor(eval(substitute(gen), eval(data)))
  ENV = factor(eval(substitute(env), eval(data)))
  REP = factor(eval(substitute(rep), eval(data)))
  data = data.frame(ENV, GEN, REP, Y)
  Nenv = length(unique(ENV))
  Ngen = length(unique(GEN))
  Nbloc = length(unique(REP))
  minimo = min(Nenv, Ngen) - 1
  ovmean = mean(Y)
  PesoWAASB = weight.WAAS
  PesoResp = weight.response
  options(lmerControl = list(check.nobs.vs.rankZ = "warning",
                             check.nobs.vs.nlev = "warning",
                             check.nobs.vs.nRE = "warning",
                             check.nlev.gtreq.5 = "warning",
                             check.nlev.gtr.1 = "warning"))

  if (minimo < 2) {
    cat("\nWarning. The analysis is not possible.")
    cat("\nThe number of environments and number of genotypes must be greater than 2\n")
  }

temp = NULL
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
    if (h2 < 0) {AS = 0} else {AS = sqrt(h2)}
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

    actualenv = actualenv + 1

  }
names(temp) = c("ENV", "Mean", "MSblock", "MSgen", "MSres", "Fcal(Blo)", "Pr>F(Blo)","Fcal(Gen)", "Pr>F(Gen)", "CV(%)", "h2", "AS")
  MSEratio = max(temp$MSres) / min(temp$MSres)
  individual = list(individual = temp,
                    MSEratio = MSEratio)

  if (random  ==  "all"){

    Complete = suppressWarnings(suppressMessages(lme4::lmer(Y ~  REP%in%ENV + (1|ENV) + (1|GEN) + (1|GEN:ENV))))
    reducedG = suppressWarnings(suppressMessages(lme4::lmer(Y ~  REP%in%ENV + (1|ENV) + (1|GEN:ENV))))
    reducedE = suppressWarnings(suppressMessages(lme4::lmer(Y ~  REP%in%ENV + (1|GEN) + (1|GEN:ENV))))
    ReducedGE = suppressWarnings(suppressMessages(lme4::lmer(Y ~  REP%in%ENV + (1|ENV) + (1|GEN))))
    gtest = suppressWarnings(suppressMessages(data.frame(t(anova(Complete, reducedG)))))
    etest = suppressWarnings(suppressMessages(data.frame(t(anova(Complete, reducedE)))))
    geatest = suppressWarnings(suppressMessages(data.frame(t(anova(Complete, ReducedGE)))))
    LRT = cbind(gtest, etest, geatest)
    model = Complete
    summary(model)
    fixed = broom::tidy(model, effects = "fixed", conf.int = TRUE)
    random = broom::tidy(model, effects = "ran_pars")
    random = random[with(random, order(group)), ]
    statistics = broom::glance(model)
    REML =  list(fixed = fixed, random = random, statistics = statistics)
    EV = as.numeric(random[1,3]^2)
    GV = as.numeric(random[2,3]^2)
    GEV = as.numeric(random[3,3]^2)
    RV = as.numeric(random[4,3]^2)
    FV =  GEV + GV + RV
    h2g = GV / (GV + GEV + RV)
    h2mg = GV/(GV + GEV/Nenv + RV/(Nenv * Nbloc))
    GEr2 = GEV / (GV + GEV + RV)
    AccuGen = sqrt(h2mg)
    rge = GEV / (GEV + RV)
    CVg = (sqrt(GV)/ovmean) * 100
    CVr = (sqrt(RV)/ovmean) * 100
    CVratio = CVg/CVr
    PROB = ((1-prob)/2)+prob
    t = qt(PROB, 100)
    Limits = t * sqrt(((1-AccuGen) * GV))
    GEVper = (GEV/FV) * 100
    GVper = (GV/FV) * 100
    RVper = (RV/FV) * 100
    GEV = paste0(round(GEV, 6), " (", round(GEVper,2), "% of fenotypic variance.)")
    GV = paste0(round(GV,6), " (", round(GVper,2), "% of fenotypic variance.)")
    RV = paste0(round(RV,6), " (", round(RVper,2), "% of fenotypic variance.)")
    ESTIMATES = list(GEV = GEV,
                     GV = GV,
                     EV = EV,
                     RV = RV,
                     FV = FV,
                     h2g = h2g,
                     GEr2 = GEr2,
                     h2mg = h2mg,
                     AccuGen = AccuGen,
                     rge = rge,
                     CVg = CVg,
                     CVr = CVr,
                     CVratio = CVratio)
    ESTIMATES = do.call(rbind.data.frame, ESTIMATES)
    names(ESTIMATES) = "Values"
    ESTIMATES = dplyr::mutate(ESTIMATES,
                             Parameters  = c("GEI variance", "Genotypic variance", "Environmental variance", "Residual variance",
                                             "Phenotypic variance", "Heritability", "GEIr2",
                                             "Heribatility of means", "Accuracy", "rge", "CVg", "CVr", "CV ratio") )
    ESTIMATES = ESTIMATES %>%
      dplyr::select(Parameters, everything())
    bups = lme4::ranef(model)
    blupGEN = bups$GEN
    blupINT = bups$`GEN:ENV`
    blups = data.frame(Names = rownames(blupINT))
    blups = data.frame(do.call('rbind', strsplit(as.character(blups$Names),':',fixed = TRUE)))
    blups = blups %>%
      select(-X1, everything())
    blups = cbind(blups, blupINT)
    names(blups) = c("Code", "GEN", "BLUPge")
    blups = blups[gtools::mixedorder(blups[,1]),]
    intmatrix = by(blups[, 3], blups[, c(2, 1)], function(x) sum(x, na.rm = TRUE))
    s = svd(intmatrix)
    U = s$u[,1:minimo]
    LL = diag(s$d[1:minimo])
    V = s$v[,1:minimo]
    Eigenvalue = data.frame(Eigenvalue = s$d[1:minimo]^2)
    Eigenvalue =dplyr::mutate(Eigenvalue, Proportion = s$d[1:minimo]^2/sum(s$d[1:minimo]^2) * 100)
    Eigenvalue =dplyr::mutate(group_by(Eigenvalue), Accumulated = cumsum(Proportion))
    Eigenvalue$PC  = rownames(Eigenvalue)
    Eigenvalue = Eigenvalue %>%
      dplyr::select(PC, everything())
    Eigenvalue = as.data.frame(Eigenvalue)
    SCOREG = U%*%LL^0.5
    SCOREE = V%*%LL^0.5
    Escores = rbind(SCOREG, SCOREE)
    colnames(Escores) = paste("PC", 1:minimo, sep = "")
    raw = data.frame(ENV, GEN, Y)
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
    MEDIAS = data.frame(ENV = x, GEN = y, Y = z)
    OUTMED = by(MEDIAS[, 3], MEDIAS[, c(2, 1)], function(x) sum(x,
                                                                na.rm = TRUE))
    MEscores = Escores[1:Ngen, ]
    NEscores = Escores[(Ngen + 1):(Ngen + Nenv), ]
    MGEN = data.frame(type = "GEN", Code = rownames(OUTMED),  Y = apply(OUTMED, 1, mean),
                      MEscores)
    MENV = data.frame(type = "ENV", Code = colnames(OUTMED), Y = apply(OUTMED, 2, mean),
                      NEscores)
    Escores = rbind(MGEN, MENV)


    names(MENV)[2] = c("ENV")
    names(MGEN)[2] = c("GEN")
    names(MGEN)[3] = c("y")
    MEDIAS = suppressMessages(dplyr::mutate(MEDIAS,
                                            envPC1 = left_join(MEDIAS, MENV %>% select(ENV, PC1))$PC1,
                                            genPC1 = left_join(MEDIAS, MGEN %>% select(GEN, PC1))$PC1,
                                            nominal = left_join(MEDIAS, MGEN %>% select(GEN, y))$y + genPC1*envPC1))
    names(MENV)[2] = c("Code")
    names(MGEN)[2] = c("Code")
    names(MGEN)[3] = c("Y")

    Pesos = data.frame(Percent = Eigenvalue$Proportion)
    WAAS = Escores
    WAASAbs = Escores
    for (i in 4:ncol(WAAS)){
      WAAS[,i] = abs(WAAS[i])
    }
    t_WAAS = data.table::transpose(WAAS)
    colnames(t_WAAS)  = rownames(WAAS)
    rownames(t_WAAS)  = colnames(WAAS)
    t_WAAS = t_WAAS[-c(1, 2, 3), ]
    t_WAAS = t_WAAS[c(1:minimo), ]
    t_WAAS = cbind(t_WAAS, Pesos)
    for (i in 1:ncol(t_WAAS)){
      t_WAAS[,i] = as.numeric(as.character(t_WAAS[,i]))
    }

    Ponderado = data.table::transpose(as.data.frame(sapply(t_WAAS[ ,-ncol(t_WAAS)], weighted.mean,  w = t_WAAS$Percent)))
    rownames(Ponderado) = c("WAASB")
    t_WAAS = subset(t_WAAS, select = -Percent)
    colnames(Ponderado) = colnames(t_WAAS)
    t_WAAS = rbind(t_WAAS, Ponderado)
    t_WAAS2 = data.table::transpose(t_WAAS)
    colnames(t_WAAS2)  = rownames(t_WAAS)
    rownames(t_WAAS2)  = colnames(t_WAAS)
    WAASAbs = cbind(WAASAbs, subset(t_WAAS2, select = WAASB))
    WAASAbs2 = subset(WAASAbs, type == "ENV")
    WAASAbs2$PctResp = (WAASAbs2$Y / max(WAASAbs2$Y)) * 100
    WAASAbs2$PctWAASB = (100-WAASAbs2$WAASB / min(WAASAbs2$WAASB))
    WAASAbs3 = subset(WAASAbs, type == "GEN")
    WAASAbs3$PctResp = (WAASAbs3$Y / max(WAASAbs3$Y)) * 100
    WAASAbs3$PctWAASB = (100-WAASAbs3$WAASB / min(WAASAbs3$WAASB))
    WAASAbs = rbind(WAASAbs3, WAASAbs2)
    WAASAbs = data.table::setDT(WAASAbs)[, OrResp:= rank(-Y), by = type][]
    WAASAbs = data.table::setDT(WAASAbs)[, OrWAASB:= rank(WAASB), by = type][]
    WAASAbs = data.table::setDT(WAASAbs)[, OrPC1:= rank(abs(PC1)), by = type][]
    WAASAbs$PesRes = as.vector(PesoResp)
    WAASAbs$PesWAASB = as.vector(PesoWAASB)
    WAASAbs = dplyr::mutate(WAASAbs,
                           WAASY = (PctResp * PesRes + PctWAASB * PesWAASB)/(PesRes + PesWAASB))
    WAASAbsInicial = data.table::setDT(WAASAbs)[, OrWAASY:= rank(-WAASY), by = type][]
    MinENV = WAASAbs2[head(which(WAASAbs2[,3] <= min(WAASAbs2$Y)), n = 1),]
    MinENV = paste0("Environment ", MinENV$Code , " (", round(MinENV$Y,4), ") ")
    MaxENV = WAASAbs2[head(which(WAASAbs2[,3] >= max(WAASAbs2$Y)), n = 1),]
    MaxENV = paste0("Environment ", MaxENV$Code , " (", round(MaxENV$Y,4), ") ")
    MinGEN = WAASAbs3[head(which(WAASAbs3[,3] <= min(WAASAbs3$Y)),n = 1),]
    MinGEN = paste0("Genotype ", MinGEN$Code , " (", round(MinGEN$Y,4), ") ")
    MaxGEN = WAASAbs3[head(which(WAASAbs3[,3] >= max(WAASAbs3$Y)), n =1),]
    MaxGEN = paste0("Genotype ", MaxGEN$Code , " (", round(MaxGEN$Y,4), ") ")
    mean = round(mean(MEDIAS$Y),4)
    min = MEDIAS[head(which(MEDIAS[,3] <= min(MEDIAS$Y)),n = 1),]
    min = paste0(round(min$Y,4), " (Genotype ", min$GEN , " in ", min$ENV, " )")
    max = MEDIAS[head(which(MEDIAS[,3] >= max(MEDIAS$Y)),n = 1),]
    max = paste0(round(max$Y,4), " (Genotype ", max$GEN , " in ", max$ENV, " )")
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
    Details = dplyr::mutate(Details,
                           Parameters = c("WgtResponse", "WgtWAAS", "Ngen", "Nenv", "OVmean", "Min",
                                          "Max", "MinENV", "MaxENV", "MinGEN", "MaxGEN"))
    Details = Details %>%
      dplyr::select(Parameters, everything())
    blupGEN = cbind(GEN = MGEN$Code, BLUP = blupGEN)
    colnames(blupGEN) = c("GEN", "BLUPg")
    blupGEN =dplyr::mutate(blupGEN,
                     Predicted = BLUPg + ovmean)
    blupGEN = blupGEN[order(-blupGEN[,3]),]
    blupGEN =dplyr::mutate(blupGEN,
                     Rank = rank(-blupGEN[,3]),
                     LL = Predicted - Limits,
                     UL = Predicted + Limits)
    blupGEN = blupGEN %>%
      dplyr::select(Rank, everything())
    selectioNenv = suppressMessages(dplyr::left_join(blups,
                                                     blupGEN %>% select(GEN, BLUPg)))
    selectioNenv = suppressMessages(dplyr::mutate(selectioNenv,
                                                 gge = BLUPge + BLUPg,
                                                 Predicted = BLUPge + BLUPg + left_join(blups, MENV %>% select(Code, Y))$Y,
                                                 LL = Predicted - Limits,
                                                 UL = Predicted + Limits))
    names(selectioNenv) = c("ENV", "GEN", "BLUPge", "BLUPg", "BLUPg+ge", "Predicted", "LL", "UL")
    return(structure(list(individual = individual,
                          WAASB = WAASAbsInicial,
                          BLUPgen = blupGEN,
                          BLUPgge = selectioNenv,
                          PCA = Eigenvalue,
                          MeansGxE = MEDIAS,
                          Details = Details,
                          REML = REML,
                          ESTIMATES = ESTIMATES,
                          LRT = LRT),
                     class = "WAASB"))
  }
  if (random  ==  "gen"){

    Complete = suppressWarnings(suppressMessages(lme4::lmer(Y ~  REP%in%ENV + ENV + (1|GEN)+ (1|GEN:ENV))))
    reducedG = suppressWarnings(suppressMessages(lme4::lmer(Y ~  REP%in%ENV + ENV + (1|GEN:ENV))))
    ReducedGE = suppressWarnings(suppressMessages(lme4::lmer(Y ~  REP%in%ENV + ENV + (1|GEN))))
    gtest = suppressWarnings(suppressMessages(data.frame(t(anova(Complete, reducedG)))))
    geatest = suppressWarnings(suppressMessages(data.frame(t(anova(Complete, ReducedGE)))))
    LRT = cbind(gtest, geatest)
    model = Complete
    summary(model)
    fixed = broom::tidy(model, effects = "fixed", conf.int = TRUE)
    random = broom::tidy(model, effects = "ran_pars")
    random = random[with(random, order(group)), ]
    statistics = broom::glance(model)
    REML =  list(fixed = fixed, random = random, statistics = statistics)
    GV = as.numeric(random[1,3]^2)
    GEV = as.numeric(random[2,3]^2)
    RV = as.numeric(random[3,3]^2)
    FV =  GEV + GV + RV
    h2g = GV / (GV + GEV + RV)
    h2mg = GV/(GV + GEV/Nenv + RV/(Nenv * Nbloc))
    GEr2 = GEV / (GV + GEV + RV)
    AccuGen = sqrt(h2mg)
    rge = GEV / (GEV + RV)
    CVg = (sqrt(GV)/ovmean) * 100
    CVr = (sqrt(RV)/ovmean) * 100
    CVratio = CVg/CVr
    PROB = ((1-prob)/2)+prob
    t = qt(PROB, 100)
    Limits = t * sqrt(((1-AccuGen) * GV))
    GEVper = (GEV/FV) * 100
    GVper = (GV/FV) * 100
    RVper = (RV/FV) * 100
    GEV = paste0(round(GEV, 6), " (", round(GEVper,2), "% of fenotypic variance.)")
    GV = paste0(round(GV,6), " (", round(GVper,2), "% of fenotypic variance.)")
    RV = paste0(round(RV,6), " (", round(RVper,2), "% of fenotypic variance.)")
    ESTIMATES = list(GEV = GEV,
                     GV = GV,
                     RV = RV,
                     FV = FV,
                     h2g = h2g,
                     GEr2 = GEr2,
                     h2mg = h2mg,
                     AccuGen = AccuGen,
                     rge = rge,
                     CVg = CVg,
                     CVr = CVr,
                     CVratio = CVratio)
    ESTIMATES = do.call(rbind.data.frame, ESTIMATES)
    names(ESTIMATES) = "Values"
    ESTIMATES = dplyr::mutate(ESTIMATES,
                             Parameters = c("GEI variance", "Genotypic variance", "Residual variance",
                                            "Phenotypic variance", "Heritability", "GEIr2",
                                            "Heribatility of means", "Accuracy", "rge", "CVg", "CVr", "CV ratio") )
    ESTIMATES = ESTIMATES %>%
      dplyr::select(Parameters, everything())
    bups = lme4::ranef(model)
    blupGEN = bups$GEN
    blupINT = bups$`GEN:ENV`
    blups = data.frame(Names = rownames(blupINT))
    blups = data.frame(do.call('rbind', strsplit(as.character(blups$Names),':',fixed = TRUE)))
    blups = blups %>%
      dplyr::select(-X1, everything())
    blups = cbind(blups, blupINT)
    names(blups) = c("Code", "GEN", "BLUPge")
    blups = blups[gtools::mixedorder(blups[,1]),]
    intmatrix = by(blups[, 3], blups[, c(2, 1)], function(x) sum(x, na.rm = TRUE))
    s = svd(intmatrix)
    U = s$u[,1:minimo]
    LL = diag(s$d[1:minimo])
    V = s$v[,1:minimo]
    Eigenvalue = data.frame(Eigenvalue = s$d[1:minimo]^2)
    Eigenvalue =dplyr::mutate(Eigenvalue, Proportion = s$d[1:minimo]^2/sum(s$d[1:minimo]^2) * 100)
    Eigenvalue =dplyr::mutate(group_by(Eigenvalue), Accumulated = cumsum(Proportion))
    Eigenvalue$PC  = rownames(Eigenvalue)
    Eigenvalue = Eigenvalue %>%
      dplyr::select(PC, everything())
    Eigenvalue = as.data.frame(Eigenvalue)
    SCOREG = U %*% LL^0.5
    SCOREE = V %*% LL^0.5
    Escores = rbind(SCOREG, SCOREE)
    colnames(Escores) = paste("PC", 1:minimo, sep = "")
    raw = data.frame(ENV, GEN, Y)
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
    MEDIAS = data.frame(ENV = x, GEN = y, Y = z)
    OUTMED = by(MEDIAS[, 3], MEDIAS[, c(2, 1)], function(x) sum(x,
                                                                na.rm = TRUE))
    MEscores = Escores[1:Ngen, ]
    NEscores = Escores[(Ngen + 1):(Ngen + Nenv), ]
    MGEN = data.frame(type = "GEN", Code = rownames(OUTMED),  Y = apply(OUTMED, 1, mean),
                      MEscores)
    MENV = data.frame(type = "ENV", Code = colnames(OUTMED), Y = apply(OUTMED, 2, mean),
                      NEscores)
    Escores = rbind(MGEN, MENV)

    names(MENV)[2] = c("ENV")
    names(MGEN)[2] = c("GEN")
    names(MGEN)[3] = c("y")
    MEDIAS = suppressMessages(dplyr::mutate(MEDIAS,
                           envPC1 = left_join(MEDIAS, MENV %>% select(ENV, PC1))$PC1,
                           genPC1 = left_join(MEDIAS, MGEN %>% select(GEN, PC1))$PC1,
                           nominal = left_join(MEDIAS, MGEN %>% select(GEN, y))$y + genPC1*envPC1))
    names(MENV)[2] = c("Code")
    names(MGEN)[2] = c("Code")
    names(MGEN)[3] = c("Y")

    Pesos = data.frame(Percent = Eigenvalue$Proportion)
    WAAS = Escores
    WAASAbs = Escores
    for (i in 4:ncol(WAAS)){
      WAAS[,i] = abs(WAAS[i])
    }
    t_WAAS = data.table::transpose(WAAS)
    colnames(t_WAAS)  = rownames(WAAS)
    rownames(t_WAAS)  = colnames(WAAS)
    t_WAAS = t_WAAS[-c(1, 2, 3), ]
    t_WAAS = t_WAAS[c(1:minimo), ]
    t_WAAS = cbind(t_WAAS, Pesos)
    for (i in 1:ncol(t_WAAS)){
      t_WAAS[,i] = as.numeric(as.character(t_WAAS[,i]))
    }
    Ponderado = data.table::transpose(as.data.frame(sapply(t_WAAS[ ,-ncol(t_WAAS)], weighted.mean,  w = t_WAAS$Percent)))
    rownames(Ponderado) = c("WAASB")
    t_WAAS = subset(t_WAAS, select = -Percent)
    colnames(Ponderado) = colnames(t_WAAS)
    t_WAAS = rbind(t_WAAS, Ponderado)
    t_WAAS2 = data.table::transpose(t_WAAS)
    colnames(t_WAAS2)  = rownames(t_WAAS)
    rownames(t_WAAS2)  = colnames(t_WAAS)
    WAASAbs = cbind(WAASAbs, subset(t_WAAS2, select = WAASB))
    WAASAbs2 = subset(WAASAbs, type == "ENV")
    WAASAbs2$PctResp = (WAASAbs2$Y / max(WAASAbs2$Y)) * 100
    WAASAbs2$PctWAASB = (100-WAASAbs2$WAASB / min(WAASAbs2$WAASB))
    WAASAbs3 = subset(WAASAbs, type == "GEN")
    WAASAbs3$PctResp = (WAASAbs3$Y / max(WAASAbs3$Y)) * 100
    WAASAbs3$PctWAASB = (100-WAASAbs3$WAASB / min(WAASAbs3$WAASB))
    WAASAbs = rbind(WAASAbs3, WAASAbs2)
    WAASAbs = data.table::setDT(WAASAbs)[, OrResp:= rank(-Y), by = type][]
    WAASAbs = data.table::setDT(WAASAbs)[, OrWAASB:= rank(WAASB), by = type][]
    WAASAbs = data.table::setDT(WAASAbs)[, OrPC1:= rank(abs(PC1)), by = type][]
    WAASAbs$PesRes = as.vector(PesoResp)
    WAASAbs$PesWAASB = as.vector(PesoWAASB)
    WAASAbs = dplyr::mutate(WAASAbs,
                           WAASY = (PctResp * PesRes + PctWAASB * PesWAASB)/(PesRes + PesWAASB))
    WAASAbsInicial = data.table::setDT(WAASAbs)[, OrWAASY:= rank(-WAASY), by = type][]
    MinENV = WAASAbs2[head(which(WAASAbs2[,3] <= min(WAASAbs2$Y)), n = 1),]
    MinENV = paste0("Environment ", MinENV$Code , " (", round(MinENV$Y,4), ") ")
    MaxENV = WAASAbs2[head(which(WAASAbs2[,3] >= max(WAASAbs2$Y)), n = 1),]
    MaxENV = paste0("Environment ", MaxENV$Code , " (", round(MaxENV$Y,4), ") ")
    MinGEN = WAASAbs3[head(which(WAASAbs3[,3] <= min(WAASAbs3$Y)),n = 1),]
    MinGEN = paste0("Genotype ", MinGEN$Code , " (", round(MinGEN$Y,4), ") ")
    MaxGEN = WAASAbs3[head(which(WAASAbs3[,3] >= max(WAASAbs3$Y)), n =1),]
    MaxGEN = paste0("Genotype ", MaxGEN$Code , " (", round(MaxGEN$Y,4), ") ")
    mean = round(mean(MEDIAS$Y),4)
    min = MEDIAS[head(which(MEDIAS[,3] <= min(MEDIAS$Y)),n = 1),]
    min = paste0(round(min$Y,4), " (Genotype ", min$GEN , " in ", min$ENV, " )")
    max = MEDIAS[head(which(MEDIAS[,3] >= max(MEDIAS$Y)),n = 1),]
    max = paste0(round(max$Y,4), " (Genotype ", max$GEN , " in ", max$ENV, " )")
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
    Details = dplyr::mutate(Details,
                           Parameters = c("WgtResponse", "WgtWAAS", "Ngen", "Nenv", "OVmean", "Min",
                                          "Max", "MinENV", "MaxENV", "MinGEN", "MaxGEN"))
    Details = Details %>%
      dplyr::select(Parameters, everything())
    blupGEN = cbind(GEN = MGEN$Code, BLUP = blupGEN)
    colnames(blupGEN) = c("GEN", "BLUPg")
    blupGEN = dplyr::mutate(blupGEN,
                     Predicted = BLUPg+ovmean)
    blupGEN = blupGEN[order(-blupGEN[,3]),]


    blupGEN = dplyr::mutate(blupGEN,
                     Rank = rank(-blupGEN[,3]),
                     LL = Predicted - Limits,
                     UL = Predicted + Limits)


    blupGEN = blupGEN %>%
      dplyr::select(Rank, everything())
    selectioNenv = suppressMessages(dplyr::left_join(blups,
                                                     blupGEN %>% select(GEN, BLUPg)))
    selectioNenv = suppressMessages(dplyr::mutate(selectioNenv,
                                                 gge = BLUPge + BLUPg,
                                                 Predicted = BLUPge + BLUPg + left_join(blups, MENV %>% select(Code, Y))$Y,
                                                 LL = Predicted - Limits,
                                                 UL = Predicted + Limits))
    names(selectioNenv) = c("ENV", "GEN", "BLUPge", "BLUPg", "BLUPg+ge", "Predicted", "LL", "UL")

    return(structure(list(individual = individual,
                          model = WAASAbsInicial,
                          BLUPgen = blupGEN,
                          BLUPgge = selectioNenv,
                          PCA = Eigenvalue,
                          MeansGxE = MEDIAS,
                          Details = Details,
                          REML = REML,
                          ESTIMATES = ESTIMATES,
                          LRT = LRT),
                     class = "WAASB"))
  }
}
