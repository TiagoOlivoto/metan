anova_ind = function(.data,
                     env,
                     gen,
                     rep,
                     resp){
  group_var <- dplyr::enquo(env)
  genot <- dplyr::enquo(gen)
  replicate <- dplyr::enquo(rep)
  resp <- dplyr::enquo(resp)
  grouped =  .data  %>%
    dplyr::select(!!group_var, !!genot, !!replicate, !!resp) %>%
    dplyr::group_by(!!!group_var, !!!genot, !!!replicate) %>%
    dplyr::summarise(mean = mean(!!resp)) %>%
    dplyr::ungroup() %>%
    split(dplyr::pull(., !!group_var))
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
  return(list(individual = individual,
              MSRratio = max(individual$MSR)/min(individual$MSR)))
}

