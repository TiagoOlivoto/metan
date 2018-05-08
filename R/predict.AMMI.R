predict.AMMI = function(object,
                        resp,
                        naxis,
                        ...){

  Y = object[paste(resp)]
  data = as.data.frame(object[,1:3])
  data = cbind(data, Y)
  names(data) = c("ENV", "GEN", "REP", "Y")
  data$ENV = as.factor(data$ENV)
  data$GEN = as.factor(data$GEN)
  data$REP = as.factor(data$REP)
  Nenv = length(unique(data$ENV))
  Ngen = length(unique(data$GEN))
  Nbloc = length(unique(data$REP))
  minimo = min(Nenv, Ngen) - 1
    if (naxis > minimo){
      stop("The number of axis to be used must be lesser than or equal to min(GEN-1;ENV-1), in this case, ", minimo,".")
    } else{

      if (naxis == 0){
        stop("Invalid argument. The AMMI0 model is calculated automatically. Please, inform naxis > 0")
      } else{

        ENV = as.factor(data$ENV)
        GEN = as.factor(data$GEN)
        BLOCO = as.factor(data$REP)
        Resp = as.numeric(data$Y)
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
        x1  = model.matrix(~factor(MEDIAS$ENV) -1)
        z1  = model.matrix(~factor(MEDIAS$GEN) -1)
        modelo1 = lm(Y ~ ENV + GEN, data = MEDIAS)
        residual = modelo1$residuals
        MEDIAS = data.frame(MEDIAS, resOLS = residual)
        intmatrix = t(matrix(MEDIAS$resOLS, Nenv, byrow = T))
        s = svd(intmatrix)
        U = s$u[,1:naxis]
        LL = s$d[1:naxis]
        V = s$v[,1:naxis]
        AMMI = ((z1 %*% U)*(x1 %*% V)) %*% LL
        Estimated = mutate(MEDIAS,
                        Ypred = Y - resOLS,
                        ResAMMI =  AMMI,
                        YpredAMMI = Ypred + ResAMMI,
                        AMMI0 = Ypred)
        return(Estimated)
      }
    }
  }
