ge_stats = function(.data,
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
individual <- data %>% anova_ind(ENV, GEN, REP, mean)
data2 =  data  %>%
  dplyr::group_by(ENV, GEN) %>%
  dplyr::summarise(mean = mean(mean)) %>%
  as.data.frame()
data3 = mutate(data2,
               ge = residuals(lm(mean ~ ENV + GEN, data = data2)),
               gge = residuals(lm(mean ~ ENV, data = data2)))
ge_mean = mmat(data3, GEN, ENV, mean)
ge_effect = mmat(data3, GEN, ENV, ge)
gge_effect = mmat(data3, GEN, ENV, gge)
Mean = apply(ge_mean, 1, mean)
Ecoval = rowSums(ge_effect^2 * length(unique(data$REP)))
Variance = rowSums(apply(ge_mean, 2, function(x) (x - Mean)^2))
GENSS = Variance * length(unique(data$REP))
Ecov.perc = (Ecoval / sum(Ecoval)) * 100
ge_stat = data.frame(Mean = Mean,
               Variance = Variance,
               GENSS = GENSS,
               Ecoval = Ecoval,
               Ecov.perc = Ecov.perc)
model1 <- lm(mean ~ GEN + ENV + ENV/REP + ENV * GEN, data = data)
modav <- anova(model1)
mydf = data.frame(aggregate(mean ~ GEN + ENV, data = data, mean))
myAgg = aggregate(mean ~ GEN, mydf, "c")
iamb = data.frame(aggregate(mean ~ ENV, data = data, mean))
iamb = dplyr::mutate(iamb, IndAmb = mean - mean(mean))
iamb2 = data.frame(aggregate(mean ~ ENV + GEN, data = data, mean))
iamb2 = suppressMessages(dplyr::mutate(iamb2,
                        IndAmb = dplyr::left_join(iamb2, iamb %>% select(ENV, IndAmb))$IndAmb))
matx <- myAgg$mean
myAgg$GEN
meandf = data.frame(GEN = myAgg$GEN, myAgg$mean)
names(meandf) = c("GEN", levels(mydf$ENV))
gradyt = mean(matx)
iij = apply(matx, 2, mean) - gradyt
sqiij = sum((iij)^2)
YiIj = matx %*% iij
bij = YiIj/sqiij
svar = (apply(matx^2, 1, sum)) - (((apply(matx, 1, sum))^2)/ncol(matx))
bYijIj = bij * YiIj
dij = svar - bYijIj
devtab <- data.frame(GEN = meandf$GEN, svar, bij, YiIj, bYijIj, dij)
S2e = modav$"Mean Sq"[5]
rps = length(levels(data$REP))
en = length(levels(data$ENV))
S2di = (dij/(en - 2)) - (S2e/rps)
data2 = data2
model2 <- lm(mean ~ GEN + ENV, data = data2)
amod2 <- anova(model2)
SSL = amod2$"Sum Sq"[2]
SSGxL = amod2$"Sum Sq"[3]
SS.L.GxL = SSL + SSGxL
SSL.Linear = (1/length(levels(data$GEN))) * (colSums(matx) %*% iij)^2/sum(iij^2)
SS.L.GxL.linear = sum(bYijIj) - SSL.Linear
ge = length(levels(mydf$GEN))
Df <- c(en * ge - 1, ge - 1, ge * (en - 1), 1, ge - 1, ge * (en - 2),
        replicate(length(dij), en - 2), en * ge * (rps - 1))
poolerr = modav$"Sum Sq"[5]/rps
SSS <- c(sum(amod2$"Sum Sq"), amod2$"Sum Sq"[1], SSL + SSGxL,
         SSL.Linear, SS.L.GxL.linear, sum(dij), dij, poolerr) * rps
MSSS = (SSS/Df)
FVAL = c(NA, MSSS[2]/MSSS[6], NA, NA, MSSS[5]/MSSS[6], NA,
           MSSS[7:(length(MSSS) - 1)]/MSSS[length(MSSS)], NA)
PLINES = 1 - pf(FVAL[7:(length(MSSS) - 1)], Df[7], Df[length(Df)])
pval = c(NA, 1 - pf(FVAL[2], Df[2], Df[6]), NA, NA, 1 -
             pf(FVAL[5], Df[5], Df[6]), NA, PLINES, NA)
anovadf <- data.frame(Df, `Sum Sq` = SSS, `Mean Sq` = MSSS,
                        `F value` = FVAL, `Pr(>F)` = pval, check.names = FALSE)
rownames(anovadf) <- c("Total", "GEN", "ENV + (GEN x ENV)", "ENV (linear)",
                       " GEN x ENV (linear)", "Pooled deviation",
                         levels(data$GEN), "Pooled error")
ER = list(anova = anovadf,
          regression = data.frame(GEN = devtab$GEN,
                                  bij = devtab$bij,
                                  sdij = S2di))
p = ggplot2::ggplot(iamb2, aes(x = IndAmb, y = mean))+
  ggplot2::geom_point(aes(colour = factor(GEN)), size = 1.5)+
    geom_smooth(aes(colour = factor(GEN)), method = "lm", se = FALSE)+
    ggplot2::theme_bw()+
    ggplot2::labs(x = "Environmental index", y = varnam)+
    ggplot2::theme(axis.ticks.length = unit(.2, "cm"),
                   axis.text = element_text(size = 12, colour = "black"),
                   axis.title = element_text(size = 12, colour = "black"),
                   axis.ticks = element_line(colour = "black"),
                   plot.margin = margin(0.5, 0.5, 0.2, 0.6, "cm"),
                   axis.title.y = element_text(margin = margin(r=16)),
                   legend.title = element_blank(),
                   legend.text = element_text(size=12),
                   panel.border = element_rect(colour = "black", fill=NA, size=1),
                   panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
                   panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())

temp = list(individual = individual,
              ge_mean = ge_mean,
              ge_effect = ge_effect,
              gge_effect = gge_effect,
              ge_stats = ge_stat,
              env_index = iamb2,
              plot = p,
              ER = ER)


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
  return(structure(listres, class = "ge_stats"))
}
