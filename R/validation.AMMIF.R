validation.AMMIF = function(data,
                            resp,
                            design = "RCBD",
                            nboot,
                            nrepval,
                            progbar = TRUE){
  Y = data[paste(resp)]
  data = as.data.frame(data[,1:3])
  data = cbind(data, Y)
  names(data) = c("ENV", "GEN", "REP", "Y")
  data$ENV = as.factor(data$ENV)
  data$GEN = as.factor(data$GEN)
  data$REP= as.factor(data$REP)
  data$ID = as.numeric(rownames(data))
  Nenv = length(unique(data$ENV))
  Ngen = length(unique(data$GEN))
  Nbloc = length(unique(data$REP))
  minimo = min(Nenv, Ngen)-1
  naxisvalidation = minimo + 1
  totalboot = naxisvalidation * nboot
  initial = 0
  NAXIS = minimo
  AMMIval =data.frame(matrix(".",nboot,1))
  for(n in c(1,1:ncol(AMMIval))) {
    AMMIval[,n] = as.numeric(AMMIval[,n])
  }

  if (design == "RCBD" | design == "CRD"){
    if (nrepval != Nbloc-1){
      stop("The number replications used for validation must be equal to total number of replications -1 (In this case ", (Nbloc-1),").")
    } else{

if (progbar == TRUE){
      pb = winProgressBar(title = "the model is being built, please, wait.",
                          min = 1, max = totalboot, width = 570)
}

      for (y in 1:naxisvalidation){

        RMSPDres = data.frame(matrix(".",nboot,1))
        for (n in c(1,1:ncol(RMSPDres))) {
          RMSPDres[,n] = as.numeric(RMSPDres[,n])
        }

        for (b in 1:nboot) {

          if (design == "CRD"){
            X = sample(1:10000,1)
            set.seed(X)
            modeling = data %>%
              dplyr::group_by(ENV, GEN) %>%
              dplyr::sample_n(nrepval, replace=F)
            modeling = as.data.frame(modeling[order(modeling$ID),])
            rownames(modeling) = modeling$ID
          }

          if (design == "RCBD"){
            temp = data.frame(matrix(".",0,ncol(data)))
            for (n in c(4:5)) {
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
              temp2 = dplyr::filter(data, ENV == names)
              modeling = temp2 %>%
                dplyr::group_by(GEN) %>%
                dplyr::filter(REP %in% c(X2))
              modeling = as.data.frame(modeling)
              modeling = rbind(temp, modeling)
              temp = modeling
            }
            rownames(modeling) = modeling$ID
          }
          testing = suppressWarnings(dplyr::anti_join(data, modeling, by = c("ENV", "GEN", "REP", "Y", "ID")))
          testing = suppressWarnings(testing[order(testing[,1], testing[,2], testing[,3] ),])
          x1  = factor(testing$ENV)
          z1  = factor(testing$GEN)
          ENV = as.factor(modeling$ENV)
          GEN = as.factor(modeling$GEN)
          BLOCO = as.factor(modeling$REP)
          Resp = as.numeric(modeling$Y)
          ovmean = mean(Resp)
          raw = data.frame(ENV, GEN, Resp)
          MEDIAS = tapply(raw[, 3], raw[, c(1, 2)], mean)
          xx = rownames(MEDIAS)
          yy = colnames(MEDIAS)
          fila = length(xx)
          col = length(yy)
          total = fila * col
          x = numeric(length = total)
          y = numeric(length = total)
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
          modelo1 = lm(Y ~ ENV + GEN, data = MEDIAS)
          residual = modelo1$residuals
          MEDIAS = data.frame(MEDIAS, resOLS = residual)
          intmatrix = t(matrix(MEDIAS$resOLS, Nenv, byrow = T))
          s = svd(intmatrix)
          U = s$u[,1:NAXIS]
          LL = s$d[1:NAXIS]
          V = s$v[,1:NAXIS]
          x1  = model.matrix(~x1 -1)
          z1  = model.matrix(~z1 -1)
          AMMI = ((z1 %*% U) * (x1 %*% V)) %*% LL
          MEDIAS = mutate(MEDIAS,
                          Ypred = Y - resOLS,
                          ResAMMI =  AMMI,
                          YpredAMMI = Ypred + ResAMMI,
                          testing = testing$Y,
                          error = YpredAMMI - testing,
                          errrorAMMI0 = Ypred - testing )
          if (NAXIS == 0){
            RMSPD = sqrt(sum(MEDIAS$errrorAMMI0^2)/length(MEDIAS$errrorAMMI0))
          } else{
            RMSPD = sqrt(sum(MEDIAS$error^2)/length(MEDIAS$error))
          }

          RMSPDres[,1][b] = RMSPD

          if (NAXIS == minimo){
            ACTUAL = "AMMIF"
          }else
            ACTUAL = sprintf("AMMI%.0f",NAXIS)
          initial = initial + 1
          if (progbar == TRUE){
          ProcdAtua = b
          setWinProgressBar(pb, initial, title=paste("|Family = ",ACTUAL,"| Validating ",ProcdAtua,
                                                     " of ",nboot,"validation datasets, considering", NAXIS, "axes",
                                                     "-",round(initial/totalboot*100,2),"% Concluded -"))
          }
        }

        if (NAXIS == minimo){
          AMMIval[["AMMIF"]] = RMSPDres[,1]
        }else

          AMMIval[[sprintf("AMMI%.0f",NAXIS)]]=RMSPDres[,1]
        NAXIS = NAXIS - 1
        initial = initial
      }
      if (progbar == TRUE){
      close(pb)
     # utils::winDialog(type = "ok", "Validation sucessful! Check the results in R environment")
 }
}
    RMSPD = AMMIval[,-c(1)]
    xx = rownames(RMSPD)
    yy = colnames(RMSPD)
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
        z[k] = RMSPD[i, j]
      }
    }
    RMSPD = data.frame(MODEL = y, RMSPD = z)
    RMSPD = RMSPD[mixedorder(RMSPD[,1]),]
    RMSPDmean = plyr::ddply(RMSPD, .(MODEL), summarize, mean = mean(RMSPD))
    RMSPDmean = RMSPDmean[order(RMSPDmean[,2]),]
return(structure(list(RMSPD = RMSPD,
                RMSPDmean = RMSPDmean,
                Estimated = MEDIAS,
                Modeling = modeling,
                Testing = testing),
                class = "validation.AMMIF"))
  } else{
    stop("Incorrect experimental design informed! Plesease inform RCBD for randomized complete block or CRD for completely randomized design.")
  }
}

