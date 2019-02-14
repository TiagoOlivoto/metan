
WAASB = function(data,
                 resp,
                 gen,
                 env,
                 rep,
                 mresp = NULL,
                 wresp = NULL,
                 random = "gen",
                 prob = 0.05,
                 verbose = TRUE) {


  if (!random %in% c("env", "gen", "all")) {
    stop("The argument 'random' must be one of the 'gen', 'env', or 'all'.")
  }

  datain = data
  GEN = factor(eval(substitute(gen), eval(datain)))
  ENV = factor(eval(substitute(env), eval(datain)))
  REP = factor(eval(substitute(rep), eval(datain)))
  listres = list()
  d = match.call()
  nvar = as.numeric(ifelse(length(d$resp) > 1, length(d$resp) -1, length(d$resp)))

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

  if (random == "env") {
    for (var in 2:length(d$resp)){
      if(length(d$resp)>1){
        Y = eval(substitute(resp)[[var]], eval(datain))
        data = data.frame(ENV, GEN, REP, Y)
      } else{
        Y = eval(substitute(resp), eval(datain))
        data = data.frame(ENV, GEN, REP, Y)
      }
      Nenv = length(unique(ENV))
      Ngen = length(unique(GEN))
      minimo = min(Nenv, Ngen) - 1
      vin = vin + 1
      Nbloc = length(unique(REP))
      ovmean = mean(Y)

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
      MSEratio = max(temp$MSres)/min(temp$MSres)
      individual = list(individual = temp, MSEratio = MSEratio)

      Complete = lmerTest::lmer(data = data, Y ~  GEN + (1|ENV/REP) + (1|GEN:ENV))
      LRT = lmerTest::ranova(Complete, reduce.terms = FALSE)
      rownames(LRT) = c("Complete", "Env/Rep", "Env", "Gen vs Env")
      random = as.data.frame(lme4::VarCorr(Complete))[, c(1, 4)]
      random = random[with(random, order(grp)), ]
      names(random) = c("Group", "Variance")
      fixed = anova(Complete)
      ENVIR = as.numeric(random[1, 2])
      GEV = as.numeric(random[2, 2])
      BEV = as.numeric(random[3, 2])
      RV = as.numeric(random[4, 2])
      FV = ENVIR + GEV + BEV + RV
      ENVper = (ENVIR/FV) * 100
      GEVper = (GEV/FV) * 100
      RVper = (RV/FV) * 100
      BEVper = (BEV/FV) * 100
      GEV = paste0(round(GEV, 6), " (", round(GEVper, 2), "% of phenotypic variance.)")
      ENVIR = paste0(round(ENVIR, 6), " (", round(ENVper, 2), "% of phenotypic variance.)")
      RV = paste0(round(RV, 6), " (", round(RVper, 2), "% of phenotypic variance.)")
      BEV = paste0(round(BEV, 6), " (", round(BEVper, 2), "% of phenotypic variance.)")
      ESTIMATES = list(GEV = GEV, ENVIR = ENVIR, RV = RV, BEV = BEV,  FV = FV)
      ESTIMATES = do.call(rbind.data.frame, ESTIMATES)
      names(ESTIMATES) = "Values"
      ESTIMATES = dplyr::mutate(ESTIMATES, Parameters = c("GEI variance",
                                                          "Environment variance", "Residual variance", "Env/block variance", "Phenotypic variance"))
      ESTIMATES = ESTIMATES %>% dplyr::select(Parameters, everything())
      bups = lme4::ranef(Complete)
      blupINT = bups$`GEN:ENV`
      blups = data.frame(Names = rownames(blupINT))
      blups = data.frame(do.call("rbind", strsplit(as.character(blups$Names),
                                                   ":", fixed = TRUE)))
      blups = blups %>% dplyr::select(-X1, everything())
      blups = cbind(blups, blupINT)
      names(blups) = c("Code", "GEN", "BLUPge")
      blups = blups[gtools::mixedorder(blups[, 1]), ]
      intmatrix = by(blups[, 3], blups[, c(2, 1)], function(x) sum(x,
                                                                   na.rm = TRUE))
      s = svd(intmatrix)
      U = s$u[, 1:minimo]
      LL = diag(s$d[1:minimo])
      V = s$v[, 1:minimo]
      Eigenvalue = data.frame(Eigenvalue = s$d[1:minimo]^2)
      Eigenvalue = dplyr::mutate(Eigenvalue, Proportion = s$d[1:minimo]^2/sum(s$d[1:minimo]^2) *
                                   100)
      Eigenvalue = dplyr::mutate(group_by(Eigenvalue), Accumulated = cumsum(Proportion))
      Eigenvalue$PC = rownames(Eigenvalue)
      Eigenvalue = Eigenvalue %>% dplyr::select(PC, everything())
      Eigenvalue = as.data.frame(Eigenvalue)
      SCOREG = U %*% LL^0.5
      SCOREE = V %*% LL^0.5
      Escores = rbind(SCOREG, SCOREE)
      colnames(Escores) = paste("PC", 1:minimo, sep = "")
      MEDIAS = data.frame(data %>% dplyr::group_by(ENV, GEN) %>% dplyr::summarise(Y = mean(Y)))
      OUTMED = by(MEDIAS[, 3], MEDIAS[, c(2, 1)], function(x) sum(x,
                                                                  na.rm = TRUE))
      MEscores = Escores[1:Ngen, ]
      NEscores = Escores[(Ngen + 1):(Ngen + Nenv), ]
      MGEN = data.frame(type = "GEN", Code = rownames(OUTMED), Y = apply(OUTMED,
                                                                         1, mean), MEscores)
      MENV = data.frame(type = "ENV", Code = colnames(OUTMED), Y = apply(OUTMED,
                                                                         2, mean), NEscores)
      Escores = rbind(MGEN, MENV)
      names(MENV)[2] = c("ENV")
      names(MGEN)[2] = c("GEN")
      names(MGEN)[3] = c("y")
      MEDIAS = suppressMessages(dplyr::mutate(MEDIAS, envPC1 = left_join(MEDIAS,
                                                                         MENV %>% select(ENV, PC1))$PC1,
                                              genPC1 = left_join(MEDIAS, MGEN %>% select(GEN, PC1))$PC1,
                                              nominal = left_join(MEDIAS, MGEN %>% select(GEN, y))$y + genPC1 * envPC1))
      names(MENV)[2] = c("Code")
      names(MGEN)[2] = c("Code")
      names(MGEN)[3] = c("Y")
      Pesos = data.frame(Percent = Eigenvalue$Proportion)
      WAAS = Escores
      WAASAbs = Escores
      WAAS[, 4:ncol(WAAS)] = lapply(WAAS[,4:ncol(WAAS)], abs)
      t_WAAS = data.frame(t(WAAS))
      colnames(t_WAAS) = rownames(WAAS)
      rownames(t_WAAS) = colnames(WAAS)
      t_WAAS = t_WAAS[-c(1, 2, 3), ]
      t_WAAS = t_WAAS[c(1:minimo), ]
      t_WAAS = cbind(t_WAAS, Pesos)
      for (i in 1:ncol(t_WAAS)) {
        t_WAAS[, i] = as.numeric(as.character(t_WAAS[, i]))
      }
      Ponderado = t(as.data.frame(sapply(t_WAAS[ ,-ncol(t_WAAS)], weighted.mean,  w = t_WAAS$Percent)))
      rownames(Ponderado) = c("WAASB")
      t_WAAS = subset(t_WAAS, select = -Percent)
      colnames(Ponderado) = colnames(t_WAAS)
      t_WAAS = rbind(t_WAAS, Ponderado)
      t_WAAS2 = data.frame(t(t_WAAS))
      colnames(t_WAAS2) = rownames(t_WAAS)
      rownames(t_WAAS2) = colnames(t_WAAS)
      WAASAbs = cbind(WAASAbs, subset(t_WAAS2, select = WAASB))
      WAASAbs2 = subset(WAASAbs, type == "ENV")
      if (nvar > 1){
        WAASAbs2$PctResp = resca(WAASAbs2$Y, new_min = minresp[vin], new_max = mresp[vin])
      } else {
        WAASAbs2$PctResp = resca(WAASAbs2$Y, new_min = minresp, new_max = mresp)
      }
      WAASAbs2$PctWAASB = resca(WAASAbs2$WAASB, new_min =  100, new_max = 0)
      WAASAbs3 = subset(WAASAbs, type == "GEN")
      if (nvar > 1){
        WAASAbs3$PctResp = resca(WAASAbs3$Y, new_min = minresp[vin], new_max = mresp[vin])
      } else{
        WAASAbs3$PctResp = resca(WAASAbs3$Y, new_min = minresp, new_max = mresp)
      }
      WAASAbs3$PctWAASB = resca(WAASAbs3$WAASB, new_min =  100, new_max = 0)

      WAASAbs = rbind(WAASAbs3, WAASAbs2) %>%
        dplyr::group_by(type) %>%
        dplyr::mutate(OrResp = rank(-Y),
                      OrWAASB = rank(WAASB),
                      OrPC1 = rank(abs(PC1)))

      if (nvar > 1) {
        WAASAbs$PesRes = as.vector(PesoResp)[vin]
        WAASAbs$PesWAASB = as.vector(PesoWAASB)[vin]
      } else{
        WAASAbs$PesRes = as.vector(PesoResp)
        WAASAbs$PesWAASB = as.vector(PesoWAASB)
      }

      WAASAbs = dplyr::mutate(WAASAbs, WAASBY = (PctResp * PesRes + PctWAASB *
                                                   PesWAASB)/(PesRes + PesWAASB))

      WAASAbsInicial = WAASAbs %>%
        dplyr::group_by(type) %>%
        dplyr::mutate(OrWAASBY = rank(-WAASBY))

      MinENV = WAASAbs2[head(which(WAASAbs2[, 3] <= min(WAASAbs2$Y)),
                             n = 1), ]
      MinENV = paste0("Environment ", MinENV$Code, " (", round(MinENV$Y,
                                                               4), ") ")
      MaxENV = WAASAbs2[head(which(WAASAbs2[, 3] >= max(WAASAbs2$Y)),
                             n = 1), ]
      MaxENV = paste0("Environment ", MaxENV$Code, " (", round(MaxENV$Y,
                                                               4), ") ")
      MinGEN = WAASAbs3[head(which(WAASAbs3[, 3] <= min(WAASAbs3$Y)),
                             n = 1), ]
      MinGEN = paste0("Genotype ", MinGEN$Code, " (", round(MinGEN$Y,
                                                            4), ") ")
      MaxGEN = WAASAbs3[head(which(WAASAbs3[, 3] >= max(WAASAbs3$Y)),
                             n = 1), ]
      MaxGEN = paste0("Genotype ", MaxGEN$Code, " (", round(MaxGEN$Y,
                                                            4), ") ")
      mean = round(mean(MEDIAS$Y), 4)
      min = MEDIAS[head(which(MEDIAS[, 3] <= min(MEDIAS$Y)), n = 1),
                   ]
      min = paste0(round(min$Y, 4), " (Genotype ", min$GEN, " in ", min$ENV,
                   " )")
      max = MEDIAS[head(which(MEDIAS[, 3] >= max(MEDIAS$Y)), n = 1),
                   ]
      max = paste0(round(max$Y, 4), " (Genotype ", max$GEN, " in ", max$ENV,
                   " )")
      Details = list(Ngen = Ngen, Nenv = Nenv, OVmean = mean, Min = min, Max = max,
                     MinENV = MinENV, MaxENV = MaxENV, MinGEN = MinGEN, MaxGEN = MaxGEN)
      Details = do.call(rbind.data.frame, Details)
      names(Details) = "Values"
      Details = dplyr::mutate(Details, Parameters = c("Ngen", "Nenv", "OVmean", "Min", "Max", "MinENV",
                                                      "MaxENV", "MinGEN", "MaxGEN"))
      Details = Details %>% dplyr::select(Parameters, everything())
      Predicted = data %>% mutate(Predicted = predict(Complete))
      residuals = data.frame(fortify.merMod(Complete))
      temp = structure(list(individual = individual, fixed = fixed, random = random, LRT = LRT,
                  model = WAASAbsInicial, BLUPgen = NULL, BLUPgge = Predicted,
                  PCA = Eigenvalue, MeansGxE = MEDIAS, Details = Details,
                  ESTIMATES = ESTIMATES, residuals = residuals),
                  class = "WAASB")
      if(length(d$resp)>1){
        if (verbose == T) {
        cat("Evaluating variable", paste(d$resp[var]), round((var - 1)/(length(d$resp) - 1)*100, 1), "%", "\n")
        }
        listres[[paste(d$resp[var])]] = temp
      } else{
        listres[[paste(d$resp)]] = temp
      }
    }
  } else
    if (random == "gen") {
    for (var in 2:length(d$resp)){
      if(length(d$resp)>1){
        Y = eval(substitute(resp)[[var]], eval(datain))
        data = data.frame(ENV, GEN, REP, Y)
      } else{
        Y = eval(substitute(resp), eval(datain))
        data = data.frame(ENV, GEN, REP, Y)
      }
      Nenv = length(unique(ENV))
      Ngen = length(unique(GEN))
      minimo = min(Nenv, Ngen) - 1
      vin = vin + 1
      Nbloc = length(unique(REP))
      ovmean = mean(Y)

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

      MSEratio = max(temp$MSres)/min(temp$MSres)
      individual = list(individual = temp, MSEratio = MSEratio)

    Complete = suppressWarnings(suppressMessages(lmerTest::lmer(data = data, Y ~ REP %in% ENV + ENV + (1 | GEN) + (1 | GEN:ENV))))
    LRT = lmerTest::ranova(Complete, reduce.terms = FALSE)
    rownames(LRT) = c("Complete", "Genotype", "Gen vs Env")
    random = as.data.frame(lme4::VarCorr(Complete))[, c(1, 4)]
    random = random[with(random, order(grp)), ]
    names(random) = c("Group", "Variance")
    fixed = anova(Complete)
    GV = as.numeric(random[1, 2])
    GEV = as.numeric(random[2, 2])
    RV = as.numeric(random[3, 2])
    FV = GEV + GV + RV
    h2g = GV/FV
    h2mg = GV/(GV + GEV/Nenv + RV/(Nenv * Nbloc))
    GEr2 = GEV/(GV + GEV + RV)
    AccuGen = sqrt(h2mg)
    rge = GEV/(GEV + RV)
    CVg = (sqrt(GV)/ovmean) * 100
    CVr = (sqrt(RV)/ovmean) * 100
    CVratio = CVg/CVr
    PROB = ((1 - (1 - prob))/2) + (1 - prob)
    t = qt(PROB, 100)
    Limits = t * sqrt(((1 - AccuGen) * GV))
    GEVper = (GEV/FV) * 100
    GVper = (GV/FV) * 100
    RVper = (RV/FV) * 100
    GEV = paste0(round(GEV, 6), " (", round(GEVper, 2), "% of phenotypic variance.)")
    GV = paste0(round(GV, 6), " (", round(GVper, 2), "% of phenotypic variance.)")
    RV = paste0(round(RV, 6), " (", round(RVper, 2), "% of phenotypic variance.)")
    ESTIMATES = list(GEV = GEV, GV = GV, RV = RV, FV = FV, h2g = h2g,
                     GEr2 = GEr2, h2mg = h2mg, AccuGen = AccuGen, rge = rge, CVg = CVg,
                     CVr = CVr, CVratio = CVratio)
    ESTIMATES = do.call(rbind.data.frame, ESTIMATES)
    names(ESTIMATES) = "Values"
    ESTIMATES = dplyr::mutate(ESTIMATES, Parameters = c("GEI variance",
                                                        "Genotypic variance", "Residual variance", "Phenotypic variance",
                                                        "Heritability", "GEIr2", "Heribatility of means", "Accuracy",
                                                        "rge", "CVg", "CVr", "CV ratio"))
    ESTIMATES = ESTIMATES %>% dplyr::select(Parameters, everything())
    bups = lme4::ranef(Complete)
    blupGEN = bups$GEN
    blupINT = bups$`GEN:ENV`
    blups = data.frame(Names = rownames(blupINT))
    blups = data.frame(do.call("rbind", strsplit(as.character(blups$Names),
                                                 ":", fixed = TRUE)))
    blups = blups %>% dplyr::select(-X1, everything())
    blups = cbind(blups, blupINT)
    names(blups) = c("Code", "GEN", "BLUPge")
    blups = blups[gtools::mixedorder(blups[, 1]), ]
    intmatrix = by(blups[, 3], blups[, c(2, 1)], function(x) sum(x,
                                                                 na.rm = TRUE))
    s = svd(intmatrix)
    U = s$u[, 1:minimo]
    LL = diag(s$d[1:minimo])
    V = s$v[, 1:minimo]
    Eigenvalue = data.frame(Eigenvalue = s$d[1:minimo]^2)
    Eigenvalue = dplyr::mutate(Eigenvalue, Proportion = s$d[1:minimo]^2/sum(s$d[1:minimo]^2) *
                                 100)
    Eigenvalue = dplyr::mutate(group_by(Eigenvalue), Accumulated = cumsum(Proportion))
    Eigenvalue$PC = rownames(Eigenvalue)
    Eigenvalue = Eigenvalue %>% dplyr::select(PC, everything())
    Eigenvalue = as.data.frame(Eigenvalue)
    SCOREG = U %*% LL^0.5
    SCOREE = V %*% LL^0.5
    Escores = rbind(SCOREG, SCOREE)
    colnames(Escores) = paste("PC", 1:minimo, sep = "")
    raw = data.frame(ENV, GEN, Y)
    MEDIAS = data.frame(data %>% dplyr::group_by(ENV, GEN) %>% dplyr::summarise(Y = mean(Y)))
    OUTMED = by(MEDIAS[, 3], MEDIAS[, c(2, 1)], function(x) sum(x,
                                                                na.rm = TRUE))
    MEscores = Escores[1:Ngen, ]
    NEscores = Escores[(Ngen + 1):(Ngen + Nenv), ]
    MGEN = data.frame(type = "GEN", Code = rownames(OUTMED), Y = apply(OUTMED,
                                                                       1, mean), MEscores)
    MENV = data.frame(type = "ENV", Code = colnames(OUTMED), Y = apply(OUTMED,
                                                                       2, mean), NEscores)
    Escores = rbind(MGEN, MENV)

    names(MENV)[2] = c("ENV")
    names(MGEN)[2] = c("GEN")
    names(MGEN)[3] = c("y")
    MEDIAS = suppressMessages(dplyr::mutate(MEDIAS, envPC1 = left_join(MEDIAS,
                                                                       MENV %>% select(ENV, PC1))$PC1, genPC1 = left_join(MEDIAS,
                                                                                                                          MGEN %>% select(GEN, PC1))$PC1, nominal = left_join(MEDIAS,
                                                                                                                                                                              MGEN %>% select(GEN, y))$y + genPC1 * envPC1))
    names(MENV)[2] = c("Code")
    names(MGEN)[2] = c("Code")
    names(MGEN)[3] = c("Y")

    Pesos = data.frame(Percent = Eigenvalue$Proportion)
    WAAS = Escores
    WAASAbs = Escores
    WAAS[, 4:ncol(WAAS)] = lapply(WAAS[,4:ncol(WAAS)], abs)
    t_WAAS = data.frame(t(WAAS))
    colnames(t_WAAS) = rownames(WAAS)
    rownames(t_WAAS) = colnames(WAAS)
    t_WAAS = t_WAAS[-c(1, 2, 3), ]
    t_WAAS = t_WAAS[c(1:minimo), ]
    t_WAAS = cbind(t_WAAS, Pesos)
    for (i in 1:ncol(t_WAAS)) {
      t_WAAS[, i] = as.numeric(as.character(t_WAAS[, i]))
    }
    Ponderado = t(as.data.frame(sapply(t_WAAS[ ,-ncol(t_WAAS)], weighted.mean,  w = t_WAAS$Percent)))
    rownames(Ponderado) = c("WAASB")
    t_WAAS = subset(t_WAAS, select = -Percent)
    colnames(Ponderado) = colnames(t_WAAS)
    t_WAAS = rbind(t_WAAS, Ponderado)
    t_WAAS2 = data.frame(t(t_WAAS))
    colnames(t_WAAS2) = rownames(t_WAAS)
    rownames(t_WAAS2) = colnames(t_WAAS)
    WAASAbs = cbind(WAASAbs, subset(t_WAAS2, select = WAASB))
    WAASAbs2 = subset(WAASAbs, type == "ENV")
    if (nvar > 1){
      WAASAbs2$PctResp = resca(WAASAbs2$Y, new_min = minresp[vin], new_max = mresp[vin])
    } else {
      WAASAbs2$PctResp = resca(WAASAbs2$Y, new_min = minresp, new_max = mresp)
    }
    WAASAbs2$PctWAASB = resca(WAASAbs2$WAASB, new_min =  100, new_max = 0)
    WAASAbs3 = subset(WAASAbs, type == "GEN")
    if (nvar > 1){
      WAASAbs3$PctResp = resca(WAASAbs3$Y, new_min = minresp[vin], new_max = mresp[vin])
    } else{
      WAASAbs3$PctResp = resca(WAASAbs3$Y, new_min = minresp, new_max = mresp)
    }
    WAASAbs3$PctWAASB = resca(WAASAbs3$WAASB, new_min =  100, new_max = 0)

    WAASAbs = rbind(WAASAbs3, WAASAbs2) %>%
      dplyr::group_by(type) %>%
      dplyr::mutate(OrResp = rank(-Y),
                    OrWAASB = rank(WAASB),
                    OrPC1 = rank(abs(PC1)))

    if (nvar > 1) {
      WAASAbs$PesRes = as.vector(PesoResp)[vin]
      WAASAbs$PesWAASB = as.vector(PesoWAASB)[vin]
    } else{
      WAASAbs$PesRes = as.vector(PesoResp)
      WAASAbs$PesWAASB = as.vector(PesoWAASB)
    }

    WAASAbs = dplyr::mutate(WAASAbs, WAASBY = (PctResp * PesRes + PctWAASB *
                                                 PesWAASB)/(PesRes + PesWAASB))

    WAASAbsInicial = WAASAbs %>%
      dplyr::group_by(type) %>%
      dplyr::mutate(OrWAASBY = rank(-WAASBY))

    MinENV = WAASAbs2[head(which(WAASAbs2[, 3] <= min(WAASAbs2$Y)),
                           n = 1), ]
    MinENV = paste0("Environment ", MinENV$Code, " (", round(MinENV$Y,
                                                             4), ") ")
    MaxENV = WAASAbs2[head(which(WAASAbs2[, 3] >= max(WAASAbs2$Y)),
                           n = 1), ]
    MaxENV = paste0("Environment ", MaxENV$Code, " (", round(MaxENV$Y,
                                                             4), ") ")
    MinGEN = WAASAbs3[head(which(WAASAbs3[, 3] <= min(WAASAbs3$Y)),
                           n = 1), ]
    MinGEN = paste0("Genotype ", MinGEN$Code, " (", round(MinGEN$Y,
                                                          4), ") ")
    MaxGEN = WAASAbs3[head(which(WAASAbs3[, 3] >= max(WAASAbs3$Y)),
                           n = 1), ]
    MaxGEN = paste0("Genotype ", MaxGEN$Code, " (", round(MaxGEN$Y,
                                                          4), ") ")
    mean = round(mean(MEDIAS$Y), 4)
    min = MEDIAS[head(which(MEDIAS[, 3] <= min(MEDIAS$Y)), n = 1),
                 ]
    min = paste0(round(min$Y, 4), " (Genotype ", min$GEN, " in ", min$ENV,
                 " )")
    max = MEDIAS[head(which(MEDIAS[, 3] >= max(MEDIAS$Y)), n = 1),
                 ]
    max = paste0(round(max$Y, 4), " (Genotype ", max$GEN, " in ", max$ENV,
                 " )")
    Details = list(Ngen = Ngen, Nenv = Nenv, OVmean = mean, Min = min, Max = max,
                   MinENV = MinENV, MaxENV = MaxENV, MinGEN = MinGEN, MaxGEN = MaxGEN)
    Details = do.call(rbind.data.frame, Details)
    names(Details) = "Values"
    Details = dplyr::mutate(Details, Parameters = c("Ngen", "Nenv", "OVmean", "Min", "Max", "MinENV",
                                                    "MaxENV", "MinGEN", "MaxGEN"))
    Details = Details %>% dplyr::select(Parameters, everything())
    blupGEN = cbind(GEN = MGEN$Code, BLUP = blupGEN)
    colnames(blupGEN) = c("GEN", "BLUPg")
    blupGEN = dplyr::mutate(blupGEN, Predicted = BLUPg + ovmean)
    blupGEN = blupGEN[order(-blupGEN[, 3]), ]


    blupGEN = dplyr::mutate(blupGEN, Rank = rank(-blupGEN[, 3]), LL = Predicted -
                              Limits, UL = Predicted + Limits)


    blupGEN = blupGEN %>% dplyr::select(Rank, everything())
    selectioNenv = suppressMessages(dplyr::left_join(blups, blupGEN %>%
                                                       select(GEN, BLUPg)))
    selectioNenv = suppressMessages(dplyr::mutate(selectioNenv, gge = BLUPge +
                                                    BLUPg, Predicted = BLUPge + BLUPg + left_join(blups, MENV %>%
                                                                                                    select(Code, Y))$Y, LL = Predicted - Limits, UL = Predicted +
                                                    Limits))
    names(selectioNenv) = c("ENV", "GEN", "BLUPge", "BLUPg", "BLUPg+ge",
                            "Predicted", "LL", "UL")
    residuals = data.frame(fortify.merMod(Complete))
    residuals$reff = selectioNenv$BLUPge
    temp = structure(list(individual = individual, fixed = fixed, random = random, LRT = LRT,
                model = WAASAbsInicial, blupGEN = blupGEN, BLUPgge = selectioNenv,
                PCA = Eigenvalue, MeansGxE = MEDIAS, Details = Details,
                ESTIMATES = ESTIMATES, residuals = residuals),
                class = "WAASB")

  if(length(d$resp)>1){
    if (verbose == T) {
    cat("Evaluating variable", paste(d$resp[var]), round((var - 1)/(length(d$resp) - 1)*100, 1), "%", "\n")
    }
    listres[[paste(d$resp[var])]] = temp
  } else{
    listres[[paste(d$resp)]] = temp
  }
  }
  } else {
    for (var in 2:length(d$resp)){
      if(length(d$resp)>1){
        Y = eval(substitute(resp)[[var]], eval(datain))
        data = data.frame(ENV, GEN, REP, Y)
      } else{
        Y = eval(substitute(resp), eval(datain))
        data = data.frame(ENV, GEN, REP, Y)
      }
      Nenv = length(unique(ENV))
      Ngen = length(unique(GEN))
      minimo = min(Nenv, Ngen) - 1
      vin = vin + 1
      Nbloc = length(unique(REP))
      ovmean = mean(Y)

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

      MSEratio = max(temp$MSres)/min(temp$MSres)
      individual = list(individual = temp, MSEratio = MSEratio)

      Complete = suppressWarnings(suppressMessages(lmerTest::lmer(data = data,
                                                                  Y ~ (1 | GEN) +
                                                                    (1|ENV/REP) +
                                                                    (1|GEN:ENV))))
      LRT = lmerTest::ranova(Complete, reduce.terms = FALSE)
      rownames(LRT) = c("Complete", "Genotype", "Env/Rep", "Environment", "Gen:Env")
      random = as.data.frame(lme4::VarCorr(Complete))[, c(1, 4)]
      random = random[with(random, order(grp)), ]
      names(random) = c("Group", "Variance")
      EV = as.numeric(random[1, 2])
      GV = as.numeric(random[2, 2])
      GEV = as.numeric(random[3, 2])
      BWE = as.numeric(random[4, 2])
      RV = as.numeric(random[5, 2])
      FV = GEV + GV + EV + RV
      h2g = GV/FV
      h2mg = GV/(GV + GEV/Nenv + RV/(Nenv * Nbloc))
      GEr2 = GEV/(GV + GEV + RV)
      AccuGen = sqrt(h2mg)
      rge = GEV/(GEV + RV)
      CVg = (sqrt(GV)/ovmean) * 100
      CVr = (sqrt(RV)/ovmean) * 100
      CVratio = CVg/CVr
      PROB = ((1 - (1 - prob))/2) + (1 - prob)
      t = qt(PROB, 100)
      Limits = t * sqrt(((1 - AccuGen) * GV))
      GEVper = (GEV/FV) * 100
      GVper = (GV/FV) * 100
      RVper = (RV/FV) * 100
      EVper = (EV/FV) * 100
      GEV = paste0(round(GEV, 6), " (", round(GEVper, 2), "% of phenotypic variance.)")
      GV = paste0(round(GV, 6), " (", round(GVper, 2), "% of phenotypic variance.)")
      RV = paste0(round(RV, 6), " (", round(RVper, 2), "% of phenotypic variance.)")
      EV = paste0(round(EV, 6), " (", round(EVper, 2), "% of phenotypic variance.)")
      ESTIMATES = list(GEV = GEV, GV = GV, EV = EV, RV = RV, FV = FV,
                       h2g = h2g, GEr2 = GEr2, h2mg = h2mg, AccuGen = AccuGen, rge = rge,
                       CVg = CVg, CVr = CVr, CVratio = CVratio)
      ESTIMATES = do.call(rbind.data.frame, ESTIMATES)
      names(ESTIMATES) = "Values"
      ESTIMATES = dplyr::mutate(ESTIMATES, Parameters = c("GEI variance",
                                                          "Genotypic variance", "Environmental variance", "Residual variance",
                                                          "Phenotypic variance", "Heritability", "GEIr2", "Heribatility of means",
                                                          "Accuracy", "rge", "CVg", "CVr", "CV ratio"))
      ESTIMATES = ESTIMATES %>% dplyr::select(Parameters, everything())
      bups = lme4::ranef(Complete)
      blupINT = bups$`GEN:ENV`
      blups = data.frame(Names = rownames(blupINT))
      blups = data.frame(do.call("rbind", strsplit(as.character(blups$Names),
                                                   ":", fixed = TRUE)))
      blups = blups %>% select(-X1, everything())
      blups = cbind(blups, blupINT)
      names(blups) = c("Code", "GEN", "BLUPge")
      blups = blups[gtools::mixedorder(blups[, 1]), ]
      intmatrix = by(blups[, 3], blups[, c(2, 1)], function(x) sum(x,
                                                                   na.rm = TRUE))
      s = svd(intmatrix)
      U = s$u[, 1:minimo]
      LL = diag(s$d[1:minimo])
      V = s$v[, 1:minimo]
      Eigenvalue = data.frame(Eigenvalue = s$d[1:minimo]^2)
      Eigenvalue = dplyr::mutate(Eigenvalue, Proportion = s$d[1:minimo]^2/sum(s$d[1:minimo]^2) *
                                   100)
      Eigenvalue = dplyr::mutate(group_by(Eigenvalue), Accumulated = cumsum(Proportion))
      Eigenvalue$PC = rownames(Eigenvalue)
      Eigenvalue = data.frame(Eigenvalue %>% dplyr::select(PC, everything()))
      SCOREG = U %*% LL^0.5
      SCOREE = V %*% LL^0.5
      Escores = rbind(SCOREG, SCOREE)
      colnames(Escores) = paste("PC", 1:minimo, sep = "")
      raw = data.frame(ENV, GEN, Y)
      MEDIAS = data.frame(raw %>% group_by(ENV, GEN) %>% summarize(Y = mean(Y)))
      OUTMED = by(MEDIAS[, 3], MEDIAS[, c(2, 1)], function(x) sum(x,
                                                                  na.rm = TRUE))
      MEscores = Escores[1:Ngen, ]
      NEscores = Escores[(Ngen + 1):(Ngen + Nenv), ]
      MGEN = data.frame(type = "GEN", Code = rownames(OUTMED), Y = apply(OUTMED,
                                                                         1, mean), MEscores)
      MENV = data.frame(type = "ENV", Code = colnames(OUTMED), Y = apply(OUTMED,
                                                                         2, mean), NEscores)
      Escores = rbind(MGEN, MENV)


      names(MENV)[2] = c("ENV")
      names(MGEN)[2] = c("GEN")
      names(MGEN)[3] = c("y")
      MEDIAS = suppressMessages(dplyr::mutate(MEDIAS, envPC1 = left_join(MEDIAS,
                                                                         MENV %>% select(ENV, PC1))$PC1,
                                              genPC1 = left_join(MEDIAS, MGEN %>% select(GEN, PC1))$PC1,
                                              nominal = left_join(MEDIAS, MGEN %>% select(GEN, y))$y + genPC1 * envPC1))
      names(MENV)[2] = c("Code")
      names(MGEN)[2] = c("Code")
      names(MGEN)[3] = c("Y")

      Pesos = data.frame(Percent = Eigenvalue$Proportion)
      WAAS = Escores
      WAASAbs = Escores
      WAAS[, 4:ncol(WAAS)] = lapply(WAAS[,4:ncol(WAAS)], abs)
      t_WAAS = data.frame(t(WAAS))
      colnames(t_WAAS) = rownames(WAAS)
      rownames(t_WAAS) = colnames(WAAS)
      t_WAAS = t_WAAS[-c(1, 2, 3), ]
      t_WAAS = t_WAAS[c(1:minimo), ]
      t_WAAS = cbind(t_WAAS, Pesos)
      for (i in 1:ncol(t_WAAS)) {
        t_WAAS[, i] = as.numeric(as.character(t_WAAS[, i]))
      }
      Ponderado = t(as.data.frame(sapply(t_WAAS[ ,-ncol(t_WAAS)], weighted.mean,  w = t_WAAS$Percent)))
      rownames(Ponderado) = c("WAASB")
      t_WAAS = subset(t_WAAS, select = -Percent)
      colnames(Ponderado) = colnames(t_WAAS)
      t_WAAS = rbind(t_WAAS, Ponderado)
      t_WAAS2 = data.frame(t(t_WAAS))
      colnames(t_WAAS2) = rownames(t_WAAS)
      rownames(t_WAAS2) = colnames(t_WAAS)
      WAASAbs = cbind(WAASAbs, subset(t_WAAS2, select = WAASB))
      WAASAbs2 = subset(WAASAbs, type == "ENV")
      if (nvar > 1){
        WAASAbs2$PctResp = resca(WAASAbs2$Y, new_min = minresp[vin], new_max = mresp[vin])
      } else {
        WAASAbs2$PctResp = resca(WAASAbs2$Y, new_min = minresp, new_max = mresp)
      }
      WAASAbs2$PctWAASB = resca(WAASAbs2$WAASB, new_min =  100, new_max = 0)
      WAASAbs3 = subset(WAASAbs, type == "GEN")
      if (nvar > 1){
        WAASAbs3$PctResp = resca(WAASAbs3$Y, new_min = minresp[vin], new_max = mresp[vin])
      } else{
        WAASAbs3$PctResp = resca(WAASAbs3$Y, new_min = minresp, new_max = mresp)
      }
      WAASAbs3$PctWAASB = resca(WAASAbs3$WAASB, new_min =  100, new_max = 0)

      WAASAbs = rbind(WAASAbs3, WAASAbs2) %>%
        dplyr::group_by(type) %>%
        dplyr::mutate(OrResp = rank(-Y),
                      OrWAASB = rank(WAASB),
                      OrPC1 = rank(abs(PC1)))

      if (nvar > 1) {
        WAASAbs$PesRes = as.vector(PesoResp)[vin]
        WAASAbs$PesWAASB = as.vector(PesoWAASB)[vin]
      } else{
        WAASAbs$PesRes = as.vector(PesoResp)
        WAASAbs$PesWAASB = as.vector(PesoWAASB)
      }

      WAASAbs = dplyr::mutate(WAASAbs, WAASBY = (PctResp * PesRes + PctWAASB *
                                                   PesWAASB)/(PesRes + PesWAASB))

      WAASAbsInicial = WAASAbs %>%
        dplyr::group_by(type) %>%
        dplyr::mutate(OrWAASBY = rank(-WAASBY))

      MinENV = WAASAbs2[head(which(WAASAbs2[, 3] <= min(WAASAbs2$Y)),
                             n = 1), ]
      MinENV = paste0("Environment ", MinENV$Code, " (", round(MinENV$Y,
                                                               4), ") ")
      MaxENV = WAASAbs2[head(which(WAASAbs2[, 3] >= max(WAASAbs2$Y)),
                             n = 1), ]
      MaxENV = paste0("Environment ", MaxENV$Code, " (", round(MaxENV$Y,
                                                               4), ") ")
      MinGEN = WAASAbs3[head(which(WAASAbs3[, 3] <= min(WAASAbs3$Y)),
                             n = 1), ]
      MinGEN = paste0("Genotype ", MinGEN$Code, " (", round(MinGEN$Y,
                                                            4), ") ")
      MaxGEN = WAASAbs3[head(which(WAASAbs3[, 3] >= max(WAASAbs3$Y)),
                             n = 1), ]
      MaxGEN = paste0("Genotype ", MaxGEN$Code, " (", round(MaxGEN$Y,
                                                            4), ") ")
      mean = round(mean(MEDIAS$Y), 4)
      min = MEDIAS[head(which(MEDIAS[, 3] <= min(MEDIAS$Y)), n = 1),
                   ]
      min = paste0(round(min$Y, 4), " (Genotype ", min$GEN, " in ", min$ENV,
                   ")")
      max = MEDIAS[head(which(MEDIAS[, 3] >= max(MEDIAS$Y)), n = 1),
                   ]
      max = paste0(round(max$Y, 4), " (Genotype ", max$GEN, " in ", max$ENV,
                   ")")
      Details = list(Ngen = Ngen, Nenv = Nenv, OVmean = mean, Min = min, Max = max,
                     MinENV = MinENV, MaxENV = MaxENV, MinGEN = MinGEN, MaxGEN = MaxGEN)
      Details = do.call(rbind.data.frame, Details)
      names(Details) = "Values"
      Details = dplyr::mutate(Details, Parameters = c("Ngen", "Nenv", "OVmean", "Min", "Max", "MinENV",
                                                      "MaxENV", "MinGEN", "MaxGEN"))
      Details = Details %>% dplyr::select(Parameters, everything())

      blupGEN = bups$GEN
      blupGEN = cbind(GEN = MGEN$Code, BLUP = blupGEN)
      colnames(blupGEN) = c("GEN", "BLUPg")
      blupGEN = dplyr::mutate(blupGEN, Predicted = BLUPg + ovmean)
      blupGEN = blupGEN[order(-blupGEN[, 3]), ]
      blupGEN = dplyr::mutate(blupGEN, Rank = rank(-blupGEN[, 3]), LL = Predicted -
                                Limits, UL = Predicted + Limits)
      blupGEN = blupGEN %>% dplyr::select(Rank, everything())


      blupENV = bups$ENV
      blupENV = cbind(ENV = MENV$Code, BLUP = blupENV)
      colnames(blupENV) = c("Code", "BLUPe")
      blupENV = dplyr::mutate(blupENV, Predicted = BLUPe + ovmean)
      blupENV = blupENV[order(-blupENV[, 3]), ]
      blupENV = dplyr::mutate(blupENV, Rank = rank(-blupENV[, 3]))
      blupENV = blupENV %>% dplyr::select(Rank, everything())

      selectioNenv = suppressMessages(dplyr::left_join(blups, blupGEN %>%
                                                         select(GEN, BLUPg)))
      selectioNenv = suppressMessages(dplyr::mutate(selectioNenv, BLUPe = left_join(blups, blupENV %>% select(Code, BLUPe))$BLUPe,
                                                    ggee = BLUPge + BLUPg + BLUPe,
                                                    Predicted = ggee + ovmean))
      names(selectioNenv) = c("ENV", "GEN", "BLUPge", "BLUPg", "BLUPe",
                              "BLUPge+g+e", "Predicted")
      residuals = fortify.merMod(Complete)
      residuals$reff = selectioNenv$BLUPge
      temp = structure(list(individual = individual, fixed = NULL, random = random, LRT = LRT,
                  model = WAASAbsInicial, blupGEN = blupGEN, BLUPgge = selectioNenv,
                  PCA = Eigenvalue, MeansGxE = MEDIAS, Details = Details,
                  ESTIMATES = ESTIMATES, residuals = residuals),
                  class = "WAASB")

      if (length(d$resp) > 1) {
        if (verbose == T) {
        cat("Evaluating variable", paste(d$resp[var]), round((var - 1)/(length(d$resp) - 1)*100, 1), "%", "\n")
        }
        listres[[paste(d$resp[var])]] = temp
      } else{
        listres[[paste(d$resp)]] = temp
      }
    }
  }
if (verbose == T){
  if(length(which(unlist(lapply(listres, function(x){
    x[["LRT"]][3, 6]
  })) > prob)) > 0){
    cat("------------------------------------------------------------\n")
    cat("Variables with nonsignificant GxE interaction\n")
    cat(names(which(unlist(lapply(listres, function(x){
      x[["LRT"]][3, 6]
    })) > prob)), "\n")
    cat("------------------------------------------------------------\n")
  }
cat("Done!\n")
}
  invisible(structure(listres, class = "WAASB"))

}
