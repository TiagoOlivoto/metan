performs_ammi = function (ENV, GEN, REP, Y)
{
  resp = paste(deparse(substitute(Y)))
  ENV = as.factor(ENV)
  GEN = as.factor(GEN)
  REP = as.factor(REP)
  nenv = length(unique(ENV))
  ngen = length(unique(GEN))
  minimo = min(ngen, nenv)
  nrep = length(unique(REP))
  model = aov(Y ~ ENV + REP %in% ENV + GEN + ENV:GEN)
  datares = model$model
  datares$factors = paste(datares$ENV, datares$GEN)
  df = fortify(model)
  df = data.frame(fitted = df[, 8], resid = df[, 9], stdres = df[,
                                                                 10])
  residuals = cbind(datares, df)
  mm = anova(model)
  nn = mm[2, ]
  mm[2, ] = mm[3, ]
  mm[3, ] = nn
  row.names(mm)[2] = "REP(ENV)"
  row.names(mm)[3] = "GEN     "
  mm[1, 4] = mm[1, 3]/mm[2, 3]
  mm[1, 5] = 1 - pf(mm[1, 4], mm[1, 1], mm[2, 1])
  anova = mm
  probint = anova[4, 5]
  DFE = df.residual(model)
  MSE = deviance(model)/DFE
  raw = data.frame(ENV, GEN, Y)
  MEANS = tapply(raw[, 3], raw[, c(1, 2)], mean)
  xx = rownames(MEANS)
  yy = colnames(MEANS)
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
      z[k] = MEANS[i, j]
    }
  }
  MEANS = data.frame(ENV = x, GEN = y, Y = z)
  x = MEANS[, 1]
  y = MEANS[, 2]
  z = MEANS[, 3]
  model2 = lm(z ~ x + y)
  for (i in 1:length(z)) {
    if (is.na(z[i]))
      z[i] = predict(model2, data.frame(x = MEANS[i, 1],
                                        y = MEANS[i, 2]))
  }
  MEANS = data.frame(ENV = x, GEN = y, Y = z)
  model1 = lm(Y ~ ENV + GEN, data = MEANS)
  residual = model1$residuals
  MEANS = data.frame(MEANS, RESIDUAL = residual)
  mlabel = names(MEANS)
  names(MEANS) = c(mlabel[1:2], resp, mlabel[4])
  OUTRES = MEANS[order(MEANS[, 1], MEANS[, 2]), ]
  OUTRES2 = by(OUTRES[, 4], OUTRES[, c(2, 1)], function(x) sum(x,
                                                               na.rm = TRUE))
  OUTMED = by(OUTRES[, 3], OUTRES[, c(2, 1)], function(x) sum(x,
                                                              na.rm = TRUE))
  s = svd(OUTRES2)
  U = s$u
  L = s$d
  V = s$v
  L = L[1:minimo]
  SS = (L^2) * nrep
  SUMA = sum(SS)
  percent = (1/SUMA) * SS * 100
  DFAMMI = rep(0, minimo)
  acum = DFAMMI
  MSAMMI = DFAMMI
  F.AMMI = DFAMMI
  PROBF = DFAMMI
  acumula = 0
  for (i in 1:(minimo - 1)) {
    DF = (ngen - 1) + (nenv - 1) - (2 * i - 1)
    if (DF <= 0)
      break
    DFAMMI[i] = DF
    acumula = acumula + percent[i]
    acum[i] = acum[i] + acumula
    MSAMMI[i] = SS[i]/DFAMMI[i]
    if (MSE > 0)
      F.AMMI[i] = round(MSAMMI[i]/MSE, 2)
    else F.AMMI[i] = NA
    if (DFE > 0)
      PROBF[i] = round(1 - pf(F.AMMI[i], DFAMMI[i], DFE),
                       4)
    else PROBF[i] = NA
  }
  percent = round(percent, 1)
  acum = round(acum, 1)
  SS = round(SS, 5)
  MSAMMI = round(MSAMMI, 5)
  SSAMMI = data.frame(percent, acum, Df = DFAMMI, `Sum Sq` = SS,
                      `Mean Sq` = MSAMMI, `F value` = F.AMMI, Pr.F = PROBF)
  nssammi = nrow(SSAMMI)
  SSAMMI = SSAMMI[SSAMMI$Df > 0, ]
  nss = nrow(SSAMMI)
  row.names(SSAMMI) = paste("PC", 1:nss, sep = "")
  LL = sqrt(diag(L))
  SCOREG = U %*% LL
  SCOREE = V %*% LL
  SCORES = rbind(SCOREG, SCOREE)
  colnames(SCORES) = paste("PC", 1:nssammi, sep = "")
  MSCORES = SCORES[1:ngen, ]
  NSCORES = SCORES[(ngen + 1):(ngen + nenv), ]
  MGEN = data.frame(type = "GEN", Y = apply(OUTMED, 1, mean),
                    MSCORES)
  MENV = data.frame(type = "ENV", Y = apply(OUTMED, 2, mean),
                    NSCORES)
  bplot = rbind(MGEN, MENV)
  bplot = bplot[, 1:(nss + 2)]
  mlabel = names(bplot)
  names(bplot) = c(mlabel[1], resp, mlabel[c(-1, -2)])
  if (minimo <= 2) {
    cat("\nWarning. The analysis AMMI is not possible.")
    cat("\nThe number of environments and number of genotypes must be greater than 2\n")
  }
  PC = SSAMMI %>% dplyr::select(-percent, -acum, everything())
  anova = mm
  anova = cbind(Percent = ".", anova)
  anova = anova %>% dplyr::select(-Percent, everything())
  anova = cbind(Accumul = ".", anova)
  anova = anova %>% dplyr::select(-Accumul, everything())
  sum = as.data.frame(anova[nrow(anova), ])
  sum$Df = sum(anova$Df)
  sum$`Sum Sq` = sum(anova$`Sum Sq`)
  sum$`Mean Sq` = sum$`Sum Sq`/sum$Df
  rownames(sum) = "Total"
  ERRO = anova[nrow(anova), ]
  ERRO = rbind(ERRO, sum)
  anova = anova[-nrow(anova), -c(6,7) ]
  names(PC) = paste(c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)", "Percent", "Accumul"))
  anova = rbind_fill(anova, PC, ERRO)
  object = list(ANOVA = anova,
                analysis = PC,
                means = MEANS,
                biplot = bplot,
                residuals = residuals,
                probint = probint)
  class(object) = "WAASB"
  return(object)
}
