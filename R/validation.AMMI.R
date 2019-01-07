validation.AMMI = function(data,
                           resp,
                           gen,
                           env,
                           rep,
                           design = "RCBD",
                           nboot,
                           nrepval,
                           naxis,
                           progbar = TRUE){

  RMSPDres  = data.frame(RMSPD  = matrix(".",nboot,1))
  for (n in c(1,1:ncol(RMSPDres))) {
    RMSPDres[,n] = as.numeric(RMSPDres[,n])
  }

  Y = eval(substitute(resp), eval(data))
  GEN = factor(eval(substitute(gen), eval(data)))
  ENV = factor(eval(substitute(env), eval(data)))
  REP = factor(eval(substitute(rep), eval(data)))
  REPS = eval(substitute(rep), eval(data))
  data = data.frame(ENV, GEN, REP, Y)
  data = mutate(data, ID = rownames(data))
  Nenv = length(unique(ENV))
  Ngen = length(unique(GEN))
  Nbloc = length(unique(REP))
  minimo = min(Nenv, Ngen) - 1

  if (design == "RCBD" | design == "CRD"){
    if (nrepval > Nbloc - 1){
      stop("The number replications used for validation must be equal to total number of replications -1 (In this case ", (Nbloc-1),").")
    } else{

      if (naxis > minimo){
        stop("The number of axis to be used must be lesser than or equal to ", minimo, " [min(GEN-1;ENV-1)]")
      } else{

        if (progbar == TRUE){
          pb = winProgressBar(title = "the model is being built, please, wait.",
                              min = 1, max = nboot, width = 570)
        }

        for (b in 1:nboot) {

                    if (design == "CRD"){
            X = sample(1:10000,1)
            set.seed(X)
            modeling = data %>%
              dplyr::group_by(ENV, GEN) %>%
              dplyr::sample_n(nrepval, replace = F)
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
              X2 = sample(unique(REPS),nrepval , replace = F)
              names = factor(data$ENV)
              names = levels(names)[actualenv + 1]
              actualenv = actualenv + 1
              temp2 = dplyr::filter(data, ENV  ==   names)
              modeling = temp2 %>%
              dplyr::group_by(GEN) %>%
              dplyr::filter(REP %in% c(X2))
              modeling = as.data.frame(modeling)
              modeling = rbind(temp, modeling)
              temp = modeling
            }
            rownames(modeling) = modeling$ID
          }
          testing = dplyr::anti_join(data, modeling, by = c("ENV", "GEN", "REP", "Y", "ID"))
          testing = testing[order(testing[,1], testing[,2], testing[,3] ),]
          x1  = factor(testing$ENV)
          z1  = factor(testing$GEN)
          MEDIAS = data.frame(modeling %>% dplyr::group_by(ENV, GEN) %>% dplyr::summarise(Y = mean(Y)))
          modelo1 = lm(Y ~ ENV + GEN, data = MEDIAS)
          residual = modelo1$residuals
          intmatrix = t(matrix(residual, Nenv, byrow = T))
          s = svd(intmatrix)
          U = s$u[,1:NAXIS]
          LL = s$d[1:NAXIS]
          V = s$v[,1:NAXIS]
          x1  = model.matrix(~x1 -1)
          z1  = model.matrix(~z1 -1)
          AMMI = ((z1 %*% U) * (x1 %*% V)) %*% LL
          MEDIAS = mutate(MEDIAS,
                          Ypred = Y - residual,
                          ResAMMI =  AMMI,
                          YpredAMMI = Ypred + ResAMMI,
                          testing = testing$Y,
                          error = YpredAMMI - testing,
                          errrorAMMI0 = Ypred - testing )
          if (naxis == 0){
            RMSPD = sqrt(sum(MEDIAS$errrorAMMI0^2)/length(MEDIAS$errrorAMMI0))
          } else{
            RMSPD = sqrt(sum(MEDIAS$error^2)/length(MEDIAS$error))
          }
          RMSPDres[,1][b] = RMSPD
          if (progbar == TRUE){
            ProcdAtua = b
            setWinProgressBar(pb, b, title = paste("Validating ",ProcdAtua,
                                                   " of ",nboot,"validation datasets, considering", naxis, "axes",
                                                   "-",round(b/nboot*100,2),"% Concluded -"))
          }

        }
        if (progbar == TRUE){
          close(pb)
          utils::winDialog(type = "ok", "Validation sucessful! Check the results in R environment")
        }
      }
    }
    RSMEmean = mean(RMSPDres$RMSPD)
return(structure(list(RMSPD = RMSPDres,
                RSMEmean  = RSMEmean,
                Estimated = MEDIAS,
                Modeling = modeling,
                Testing = testing),
                class = "validation.AMMI"))
  }
  else{
    stop("Incorrect experimental design informed! Plesease inform RCBD for randomized complete block or CRD for completely randomized design.")
  }
}

