
WAASB = function(data, resp, random = "gen", prob = 0.95, weight.response = 50, weight.WAAS = 50){

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
  Nbloc = length(unique(REP))
  minimo = min(Nenv, Ngen)-1
  ovmean = mean(Y)
  PesoWAASB = weight.WAAS                   # Initial Weigth for stability
  PesoResp = weight.response                     # Initial Weigth for productivity



  ## Calculate variance components
  # requires lme4 package
  options(lmerControl=list(check.nobs.vs.rankZ = "warning",
                           check.nobs.vs.nlev = "warning",
                           check.nobs.vs.nRE = "warning",
                           check.nlev.gtreq.5 = "warning",
                           check.nlev.gtr.1 = "warning"))

  if (random == "all"){

    # GENar Model with random effects for variance components
    model = suppressWarnings(suppressMessages(lme4::lmer(Y~  REP%in%ENV + (1|GEN) + (1|ENV)  + (1|GEN:ENV))))

    # Extract variance components
    summary(model)
    fixed = broom::tidy(model, effects = "fixed", conf.int=TRUE)
    random = broom::tidy(model, effects = "ran_pars")
    random = random[with(random, order(group)), ]

    statistics = broom::glance(model)

    REML =  list(fixed = fixed, random = random, statistics = statistics)

    EV = random[1,3]^2
    GV = random[2,3]^2
    GEV = random[3,3]^2
    RV = random[4,3]^2
    FV= GEV + GV + RV
    h2g = GV / (GV + GEV + RV)

    h2mg = GV/(GV + GEV/Nbloc + RV/(Nbloc*nrow(data)))


    GEr2 = GEV / (GV + GEV + RV)

    AccuGen = sqrt(h2mg)

    rge = GEV / (GEV + RV)

    CVg = (sqrt(GV)/ovmean)*100
    CVr = (sqrt(RV)/ovmean)*100
    CVratio = CVg/CVr

    PROB = ((1-prob)/2)+prob

    t=qt(PROB, 100)

    Limits = t*sqrt(((1-AccuGen)*GV))

    GEVper = (GEV/FV)*100
    GVper = (GV/FV)*100
    RVper = (RV/FV)*100

    GEV = paste0(round(GEV, 6), " (", round(GEVper,2), "% of fenotypic variance.)")
    GV = paste0(round(GV,6), " (", round(GVper,2), "% of fenotypic variance.)")
    RV = paste0(round(RV,6), " (", round(RVper,2), "% of fenotypic variance.)")

    ESTIMATES = list(GEV = GEV, GV = GV, EV = EV, RV = RV,
                     FV = FV, h2g = h2g, GEr2 = GEr2, h2mg = h2mg,
                     AccuGen = AccuGen, rge = rge, CVg = CVg, CVr = CVr, CVratio = CVratio)

    ESTIMATES = do.call(rbind.data.frame, ESTIMATES)
    names(ESTIMATES) = "Values"
    rownames(ESTIMATES) = c("GEI variance", "Genotypic variance", "Residual variance",
                    "Phenotypic variance", "Heritability", "GEIr2",
                    "Heribatility of means", "Accuracy", "rge", "CVg", "CVr", "CV ratio")

    ## BLUPS
    # estimate BLUPS
    bups = lme4::ranef(model)
    # extract blup for GEN
    blupGEN = bups$GEN
    blupINT = bups$`GEN:ENV`

    blups = data.frame(Names=rownames(blupINT))
    blups = data.frame(do.call('rbind', strsplit(as.character(blups$Names),':',fixed=TRUE)))

    blups = blups %>%
      select(-X1, everything())
    blups = cbind(blups, blupINT)
    names(blups) = c("Code", "GEN", "BLUPge")
    blups = blups[mixedorder(blups[,1]),]


    intmatrix <- by(blups[, 3], blups[, c(2, 1)], function(x) sum(x, na.rm = TRUE))
    s <- svd(intmatrix)
    U <- s$u[,1:minimo]             # initial guess for singular vector G
    LL <- diag(s$d[1:minimo]) # initial guess for singular value
    V <- s$v[,1:minimo]             # initial guess for singular vector E

    # Eigenvalues and the proportion of variation explained by the principal components.
    Eigenvalue = data.frame(Eigenvalue=s$d[1:minimo]^2)
    Eigenvalue = mutate(Eigenvalue, Proportion=s$d[1:minimo]^2/sum(s$d[1:minimo]^2)*100)
    Eigenvalue = mutate(group_by(Eigenvalue), Accumulated=cumsum(Proportion))
    Eigenvalue$PC =rownames(Eigenvalue)

    Eigenvalue = Eigenvalue %>%
      select(PC, everything())
    Eigenvalue = as.data.frame(Eigenvalue)

    SCOREG <- U %*% LL
    SCOREE <- V %*% LL
    Escores <- rbind(SCOREG, SCOREE)
    colnames(Escores) <- paste("PC", 1:minimo, sep = "")



    raw <- data.frame(ENV, GEN, Y)
    MEDIAS <- tapply(raw[, 3], raw[, c(1, 2)], mean)
    xx <- rownames(MEDIAS)
    yy <- colnames(MEDIAS)


    fila <- length(xx)
    col <- length(yy)
    total <- fila * col
    x <- character(length = total)
    y <- character(length = total)
    z <- numeric(length = total)
    k <- 0
    for (i in 1:fila) {
      for (j in 1:col) {
        k <- k + 1
        x[k] <- xx[i]
        y[k] <- yy[j]
        z[k] <- MEDIAS[i, j]
      }
    }
    MEDIAS <- data.frame(ENV = x, GEN = y, Y = z)

    OUTMED <- by(MEDIAS[, 3], MEDIAS[, c(2, 1)], function(x) sum(x,
                                                                 na.rm = TRUE))
    MEscores <- Escores[1:Ngen, ]
    NEscores <- Escores[(Ngen + 1):(Ngen + Nenv), ]
    MGEN <- data.frame(type = "GEN", Code=rownames(OUTMED),  Y = apply(OUTMED, 1, mean),
                       MEscores)
    MENV <- data.frame(type = "ENV", Code=colnames(OUTMED), Y = apply(OUTMED, 2, mean),
                       NEscores)
    Escores <- rbind(MGEN, MENV)



    Pesos=data.frame(Percent=Eigenvalue$Proportion)


    WAAS=Escores
    WAASAbs=Escores
    for (i in 4:ncol(WAAS)){
      WAAS[,i] <- abs(WAAS[i])
    }
    t_WAAS = transpose(WAAS)
    colnames(t_WAAS)  = rownames(WAAS)
    rownames(t_WAAS)  = colnames(WAAS)
    t_WAAS=t_WAAS[-c(1, 2, 3), ]
    t_WAAS=t_WAAS[c(1:minimo), ]
    t_WAAS=cbind(t_WAAS, Pesos)
    for (i in 1:ncol(t_WAAS)){
      t_WAAS[,i] <- as.numeric(as.character(t_WAAS[,i]))
    }

    Ponderado=transpose(as.data.frame(sapply(t_WAAS[ ,-ncol(t_WAAS)], weighted.mean,  w = t_WAAS$Percent)))



    rownames(Ponderado)=c("WAASB")
    t_WAAS=subset(t_WAAS, select = -Percent)
    colnames(Ponderado)=colnames(t_WAAS)
    t_WAAS=rbind(t_WAAS, Ponderado)
    t_WAAS2=transpose(t_WAAS)
    colnames(t_WAAS2)  = rownames(t_WAAS)
    rownames(t_WAAS2)  = colnames(t_WAAS)
    WAASAbs=cbind(WAASAbs, subset(t_WAAS2, select = WAASB))
    WAASAbs2=subset(WAASAbs, type=="ENV")
    WAASAbs2$PctResp = (WAASAbs2$Y / max(WAASAbs2$Y))*100
    WAASAbs2$PctWAASB = (100-WAASAbs2$WAASB / min(WAASAbs2$WAASB))
    WAASAbs3=subset(WAASAbs, type=="GEN")
    WAASAbs3$PctResp = (WAASAbs3$Y / max(WAASAbs3$Y))*100
    WAASAbs3$PctWAASB = (100-WAASAbs3$WAASB / min(WAASAbs3$WAASB))
    WAASAbs=rbind(WAASAbs3, WAASAbs2)
    WAASAbs=setDT(WAASAbs)[, OrResp:=rank(-Y), by = type][]
    WAASAbs=setDT(WAASAbs)[, OrWAASB:=rank(WAASB), by = type][]
    WAASAbs=setDT(WAASAbs)[, OrPC1:=rank(abs(PC1)), by = type][]
    WAASAbs$PesRes=as.vector(PesoResp)
    WAASAbs$PesWAASB=as.vector(PesoWAASB)

    WAASAbs = plyr::mutate(WAASAbs,
                           WAASY = (PctResp * PesRes + PctWAASB * PesWAASB)/(PesRes + PesWAASB))

    WAASAbsInicial=data.table::setDT(WAASAbs)[, OrWAASY:=rank(-WAASY), by = type][]

    MinENV = WAASAbs2[which(WAASAbs2[,3] <=min(WAASAbs2$Y)),]
    MinENV = paste0("Environment ", MinENV$Code , " (", round(MinENV$Y,4), ") ")

    MaxENV = WAASAbs2[which(WAASAbs2[,3] >=max(WAASAbs2$Y)),]
    MaxENV = paste0("Environment ", MaxENV$Code , " (", round(MaxENV$Y,4), ") ")

    MinGEN = WAASAbs3[which(WAASAbs3[,3] <=min(WAASAbs3$Y)),]
    MinGEN = paste0("Genotype ", MinGEN$Code , " (", round(MinGEN$Y,4), ") ")


    MaxGEN = WAASAbs3[which(WAASAbs3[,3] >=max(WAASAbs3$Y)),]
    MaxGEN = paste0("Genotype ", MaxGEN$Code , " (", round(MaxGEN$Y,4), ") ")


    mean = round(mean(MEDIAS$Y),4)

    min = MEDIAS[which(MEDIAS[,3] <=min(MEDIAS$Y)),]
    min = paste0(round(min$Y,4), " (Genotype ", min$GEN , " in ", min$ENV, " )")

    max = MEDIAS[which(MEDIAS[,3] >=max(MEDIAS$Y)),]
    max = paste0(round(max$Y,4), " (Genotype ", max$GEN , " in ", max$ENV, " )")


    Details = list(WgtResponse=weight.response, WgtWAAS=weight.WAAS, Ngen=Ngen,
                   Nenv = Nenv, OVmean=mean, Min=min, Max=max, MinENV=MinENV, MaxENV=MaxENV,
                   MinGEN=MinGEN, MaxGEN=MaxGEN)

    Details = do.call(rbind.data.frame, Details)
    names(Details) = "Values"
    rownames(Details) = c("WgtResponse", "WgtWAAS", "Ngen", "Nenv", "OVmean", "Min",
                          "Max", "MinENV", "MaxENV", "MinGEN", "MaxGEN")

    ## Creating plots with the BLUPs
    # Create a numeric vector with the BLUP for each GEN
    blupGEN=cbind(GEN=MGEN$Code, BLUP=blupGEN)
    colnames(blupGEN) = c("GEN", "BLUPg")
    blupGEN=mutate(blupGEN,
                   Predicted = BLUPg+ovmean)
    blupGEN=blupGEN[order(-blupGEN[,3]),]
    blupGEN=mutate(blupGEN,
                   Rank = rank(-blupGEN[,3]),
                   LL = Predicted - Limits,
                   UL = Predicted + Limits)

    blupGEN = blupGEN %>%
      select(Rank, everything())


    selectioNenv = suppressMessages(dplyr::left_join(blups,
                                                     blupGEN %>% select(GEN, BLUPg)))

    selectioNenv = suppressMessages(plyr::mutate(selectioNenv,
                                                 gge = BLUPge + BLUPg,
                                                 Predicted = BLUPge + BLUPg + left_join(blups, MENV %>% select(Code, Y))$Y,
                                                 LL = Predicted - Limits,
                                                 UL = Predicted + Limits))

    names(selectioNenv) = c("ENV", "GEN", "BLUPge", "BLUPg", "BLUPg+ge", "Predicted", "LL", "UL")


    return(list(WAASB=WAASAbsInicial, BLUPgen = blupGEN, BLUPgge = selectioNenv, PCA = Eigenvalue,
                MeansGxE = MEDIAS, Details = Details, REML = REML, ESTIMATES = ESTIMATES))
  }

  if (random == "gen"){
    # GENar Model with random effects for variance components
    model = suppressWarnings(suppressMessages(lme4::lmer(Y~  REP%in%ENV + (1|GEN) + ENV + (1|GEN:ENV))))

    # Extract variance components
    summary(model)
    fixed = broom::tidy(model, effects = "fixed", conf.int = TRUE)
    random = broom::tidy(model, effects = "ran_pars")
    random = random[with(random, order(group)), ]

    statistics = broom::glance(model)

    REML =  list(fixed = fixed, random = random, statistics = statistics)

    GV = random[1,3]^2
    GEV = random[2,3]^2
    RV = random[3,3]^2
    FV= GEV + GV + RV
    h2g = GV / (GV + GEV + RV)

    h2mg = GV/(GV + GEV/Nbloc + RV/(Nbloc*nrow(data)))

    GEr2 = GEV / (GV + GEV + RV)

    AccuGen = sqrt(h2mg)

    rge = GEV / (GEV + RV)

    CVg = (sqrt(GV)/ovmean)*100
    CVr = (sqrt(RV)/ovmean)*100
    CVratio = CVg/CVr

    PROB = ((1-prob)/2)+prob

    t=qt(PROB, 100)

    Limits = t*sqrt(((1-AccuGen)*GV))

    GEVper = (GEV/FV)*100
    GVper = (GV/FV)*100
    RVper = (RV/FV)*100

    GEV = paste0(round(GEV, 6), " (", round(GEVper,2), "% of fenotypic variance.)")
    GV = paste0(round(GV,6), " (", round(GVper,2), "% of fenotypic variance.)")
    RV = paste0(round(RV,6), " (", round(RVper,2), "% of fenotypic variance.)")

    ESTIMATES = list(GEV = GEV, GV = GV, RV = RV,
                     FV = FV, h2g = h2g, GEr2 = GEr2, h2mg = h2mg,
                     AccuGen = AccuGen, rge = rge, CVg = CVg, CVr = CVr, CVratio = CVratio)

    ESTIMATES = do.call(rbind.data.frame, ESTIMATES)
    names(ESTIMATES) = "Values"
    rownames(ESTIMATES) = c("GEI variance", "Genotypic variance", "Residual variance",
                            "Phenotypic variance", "Heritability", "GEIr2",
                            "Heribatility of means", "Accuracy", "rge", "CVg", "CVr", "CV ratio")



    ## BLUPS
    # estimate BLUPS
    bups = lme4::ranef(model)
    # extract blup for GEN
    blupGEN = bups$GEN
    blupINT = bups$`GEN:ENV`

    blups = data.frame(Names=rownames(blupINT))
    blups = data.frame(do.call('rbind', strsplit(as.character(blups$Names),':',fixed=TRUE)))

    blups = blups %>%
      select(-X1, everything())
    blups = cbind(blups, blupINT)
    names(blups) = c("Code", "GEN", "BLUPge")
    blups=blups[mixedorder(blups[,1]),]


    intmatrix <- by(blups[, 3], blups[, c(2, 1)], function(x) sum(x, na.rm = TRUE))
    s <- svd(intmatrix)
    U <- s$u[,1:minimo]             # initial guess for singular vector G
    LL <- diag(s$d[1:minimo])    # initial guess for singular value
    V <- s$v[,1:minimo]             # initial guess for singular vector E

    # Eigenvalues and the proportion of variation explained by the principal components.
    Eigenvalue = data.frame(Eigenvalue=s$d[1:minimo]^2)
    Eigenvalue = mutate(Eigenvalue, Proportion=s$d[1:minimo]^2/sum(s$d[1:minimo]^2)*100)
    Eigenvalue = mutate(group_by(Eigenvalue), Accumulated=cumsum(Proportion))
    Eigenvalue$PC =rownames(Eigenvalue)

    Eigenvalue = Eigenvalue %>%
      select(PC, everything())
    Eigenvalue = as.data.frame(Eigenvalue)

    SCOREG <- U %*% LL
    SCOREE <- V %*% LL
    Escores <- rbind(SCOREG, SCOREE)
    colnames(Escores) <- paste("PC", 1:minimo, sep = "")



    raw <- data.frame(ENV, GEN, Y)
    MEDIAS <- tapply(raw[, 3], raw[, c(1, 2)], mean)
    xx <- rownames(MEDIAS)
    yy <- colnames(MEDIAS)


    fila <- length(xx)
    col <- length(yy)
    total <- fila * col
    x <- character(length = total)
    y <- character(length = total)
    z <- numeric(length = total)
    k <- 0
    for (i in 1:fila) {
      for (j in 1:col) {
        k <- k + 1
        x[k] <- xx[i]
        y[k] <- yy[j]
        z[k] <- MEDIAS[i, j]
      }
    }
    MEDIAS <- data.frame(ENV = x, GEN = y, Y = z)

    OUTMED <- by(MEDIAS[, 3], MEDIAS[, c(2, 1)], function(x) sum(x,
                                                                 na.rm = TRUE))
    MEscores <- Escores[1:Ngen, ]
    NEscores <- Escores[(Ngen + 1):(Ngen + Nenv), ]
    MGEN <- data.frame(type = "GEN", Code=rownames(OUTMED),  Y = apply(OUTMED, 1, mean),
                       MEscores)
    MENV <- data.frame(type = "ENV", Code=colnames(OUTMED), Y = apply(OUTMED, 2, mean),
                       NEscores)
    Escores <- rbind(MGEN, MENV)



    Pesos=data.frame(Percent=Eigenvalue$Proportion)


    WAAS=Escores
    WAASAbs=Escores
    for (i in 4:ncol(WAAS)){
      WAAS[,i] <- abs(WAAS[i])
    }
    t_WAAS = transpose(WAAS)
    colnames(t_WAAS)  = rownames(WAAS)
    rownames(t_WAAS)  = colnames(WAAS)
    t_WAAS=t_WAAS[-c(1, 2, 3), ]
    t_WAAS=t_WAAS[c(1:minimo), ]
    t_WAAS=cbind(t_WAAS, Pesos)
    for (i in 1:ncol(t_WAAS)){
      t_WAAS[,i] <- as.numeric(as.character(t_WAAS[,i]))
    }

    Ponderado=transpose(as.data.frame(sapply(t_WAAS[ ,-ncol(t_WAAS)], weighted.mean,  w = t_WAAS$Percent)))



    rownames(Ponderado)=c("WAASB")
    t_WAAS=subset(t_WAAS, select = -Percent)
    colnames(Ponderado)=colnames(t_WAAS)
    t_WAAS=rbind(t_WAAS, Ponderado)
    t_WAAS2=transpose(t_WAAS)
    colnames(t_WAAS2)  = rownames(t_WAAS)
    rownames(t_WAAS2)  = colnames(t_WAAS)
    WAASAbs=cbind(WAASAbs, subset(t_WAAS2, select = WAASB))
    WAASAbs2=subset(WAASAbs, type=="ENV")
    WAASAbs2$PctResp = (WAASAbs2$Y / max(WAASAbs2$Y))*100
    WAASAbs2$PctWAASB = (100-WAASAbs2$WAASB / min(WAASAbs2$WAASB))
    WAASAbs3=subset(WAASAbs, type=="GEN")
    WAASAbs3$PctResp = (WAASAbs3$Y / max(WAASAbs3$Y))*100
    WAASAbs3$PctWAASB = (100-WAASAbs3$WAASB / min(WAASAbs3$WAASB))
    WAASAbs=rbind(WAASAbs3, WAASAbs2)
    WAASAbs=setDT(WAASAbs)[, OrResp:=rank(-Y), by = type][]
    WAASAbs=setDT(WAASAbs)[, OrWAASB:=rank(WAASB), by = type][]
    WAASAbs=setDT(WAASAbs)[, OrPC1:=rank(abs(PC1)), by = type][]
    WAASAbs$PesRes=as.vector(PesoResp)
    WAASAbs$PesWAASB=as.vector(PesoWAASB)
    WAASAbs = plyr::mutate(WAASAbs,
                           WAASY = (PctResp * PesRes + PctWAASB * PesWAASB)/(PesRes + PesWAASB))

    WAASAbsInicial=data.table::setDT(WAASAbs)[, OrWAASY:=rank(-WAASY), by = type][]

    MinENV = WAASAbs2[which(WAASAbs2[,3] <=min(WAASAbs2$Y)),]
    MinENV = paste0("Environment ", MinENV$Code , " (", round(MinENV$Y,4), ") ")

    MaxENV = WAASAbs2[which(WAASAbs2[,3] >=max(WAASAbs2$Y)),]
    MaxENV = paste0("Environment ", MaxENV$Code , " (", round(MaxENV$Y,4), ") ")

    MinGEN = WAASAbs3[which(WAASAbs3[,3] <=min(WAASAbs3$Y)),]
    MinGEN = paste0("Genotype ", MinGEN$Code , " (", round(MinGEN$Y,4), ") ")


    MaxGEN = WAASAbs3[which(WAASAbs3[,3] >=max(WAASAbs3$Y)),]
    MaxGEN = paste0("Genotype ", MaxGEN$Code , " (", round(MaxGEN$Y,4), ") ")


    mean = round(mean(MEDIAS$Y),4)

    min = MEDIAS[which(MEDIAS[,3] <=min(MEDIAS$Y)),]
    min = paste0(round(min$Y,4), " (Genotype ", min$GEN , " in ", min$ENV, " )")

    max = MEDIAS[which(MEDIAS[,3] >=max(MEDIAS$Y)),]
    max = paste0(round(max$Y,4), " (Genotype ", max$GEN , " in ", max$ENV, " )")


    Details = list(WgtResponse=weight.response, WgtWAAS=weight.WAAS, Ngen=Ngen,
                   Nenv = Nenv, OVmean=mean, Min=min, Max=max, MinENV=MinENV, MaxENV=MaxENV,
                   MinGEN=MinGEN, MaxGEN=MaxGEN)

    Details = do.call(rbind.data.frame, Details)
    names(Details) = "Values"
    rownames(Details) = c("WgtResponse", "WgtWAAS", "Ngen", "Nenv", "OVmean", "Min",
                          "Max", "MinENV", "MaxENV", "MinGEN", "MaxGEN")

    ## Creating plots with the BLUPs
    # Create a numeric vector with the BLUP for each GEN
    blupGEN=cbind(GEN=MGEN$Code, BLUP=blupGEN)
    colnames(blupGEN) = c("GEN", "BLUPg")
    blupGEN=mutate(blupGEN,
                   Predicted = BLUPg+ovmean)
    blupGEN=blupGEN[order(-blupGEN[,3]),]
    blupGEN=mutate(blupGEN,
                   Rank = rank(-blupGEN[,3]),
                   LL = Predicted - Limits,
                   UL = Predicted + Limits)

    blupGEN = blupGEN %>%
      select(Rank, everything())


    selectioNenv = suppressMessages(dplyr::left_join(blups,
                                                     blupGEN %>% select(GEN, BLUPg)))

    selectioNenv = suppressMessages(plyr::mutate(selectioNenv,
                                                 gge = BLUPge + BLUPg,
                                                 Predicted = BLUPge + BLUPg + left_join(blups, MENV %>% select(Code, Y))$Y,
                                                 LL = Predicted - Limits,
                                                 UL = Predicted + Limits))

    names(selectioNenv) = c("ENV", "GEN", "BLUPge", "BLUPg", "BLUPg+ge", "Predicted", "LL", "UL")

    return(list(model = WAASAbsInicial, BLUPgen = blupGEN, BLUPgge = selectioNenv, PCA = Eigenvalue,
                MeansGxE = MEDIAS, Details = Details, REML = REML, ESTIMATES = ESTIMATES, object = "WAASB"))

  }

}
