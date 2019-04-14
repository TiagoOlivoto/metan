anova_ind = function(.data,
                     env,
                     gen,
                     rep,
                     resp,
                     verbose = TRUE){
  datain <- .data
  GEN <- factor(eval(substitute(gen), eval(datain)))
  ENV <- factor(eval(substitute(env), eval(datain)))
  REP <- factor(eval(substitute(rep), eval(datain)))
  listres <- list()
  d <- match.call()
  nvar <- as.numeric(ifelse(length(d$resp) > 1, length(d$resp) - 1, length(d$resp)))
  for (var in 2:length(d$resp)) {
    if (length(d$resp) > 1) {
      Y <- eval(substitute(resp)[[var]], eval(datain))
      varnam = paste(d$resp[var])
    } else {
      Y <- eval(substitute(resp), eval(datain))
      varnam = paste(d$resp)
    }
    data <- data.frame(ENV, GEN, REP, Y)
    names(data) = c("ENV", "GEN", "REP", "mean")
  grouped =  data  %>% split(dplyr::pull(., ENV))
  formula = as.formula(paste0("mean ~ ", substitute(gen), "+", substitute(rep)))
  individual = do.call(rbind,
                       lapply(grouped, function(x){
                         anova = anova(aov(formula, data = x))
                         MSB = anova[2, 3]
                         MSG = anova[1, 3]
                         MSE = anova[3, 3]
                         h2 = (MSG - MSE)/MSG
                         if (h2 < 0) {AS = 0} else {AS = sqrt(h2)}
                         final = data.frame(cbind(MEAN =  mean(x$mean),
                                                  MSB = MSB,
                                                  MSG = MSG,
                                                  MSR = MSE,
                                                  FCB = anova[2, 4],
                                                  PRFB = anova[2, 5],
                                                  FCG = anova[1, 4],
                                                  PRFG = anova[1, 5],
                                                  CV = sqrt(MSE)/mean(x$mean)*100,
                                                  h2 = h2,
                                                  AS = AS))

                       }))
  temp = list(individual = individual,
              MSRratio = max(individual$MSR)/min(individual$MSR))
  if (length(d$resp) > 1) {
    listres[[paste(d$resp[var])]] <- temp
    if (verbose == TRUE) {
      cat("Evaluating variable", paste(d$resp[var]), round((var - 1)/(length(d$resp) -
                                                                        1) * 100, 1), "%", "\n")
    }
  } else {
    listres[[paste(d$resp)]] <- temp
  }
  }
  return(listres)
}

