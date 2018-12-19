AMMI_indexes = function(x){

  if (!is(x, "WAAS.AMMI")) {
      stop("The object 'x' must be an object of class \"WAAS.AMMI\"")
  }
model = x
n = sum(model$PCA$`Pr(>F)` <= 0.05, na.rm = TRUE)
meange = model$MeansGxE
effects = residuals(lm(Y ~ ENV + GEN, data = meange))
meange$residual = effects
ge = by(meange[, 7], meange[, c(2, 1)], function(x) sum(x,
                                                     na.rm = TRUE))
ge = array(ge, dim(ge), dimnames(ge))
svdge = svd(ge)
gamma.n = svdge$u[, 1:n]
theta.n = model$PCA$Percent[1:n]/100

# sum of absolute scores
PCA = data.frame(model$model)
SCOR = PCA[PCA[, 1] == "GEN", c(seq(4, n+3))]
SIPC = unname(rowSums(apply(SCOR, 2, FUN = abs)))
mean = PCA[PCA[, 1] == "GEN", c(2:3)]
rS = rank(SIPC)
rY = rank(-mean[2])
ssiSIPC = rS + rY

# Za
Za = rowSums(abs(gamma.n %*% diag(theta.n)))
rZA = rank(Za)
ssiZA = rZA + rY
# averages of the squared eigenvector
EV = rowSums(gamma.n^2/n)
rEV = rank(EV)
ssiEV = rEV + rY

# AMMI stability values
pc = model$anova$`Sum Sq`[5]/model$anova$`Sum Sq`[6]
ASV = sqrt((pc * SCOR[,1])^2+ SCOR[,2]^2)
rASV = rank(ASV)
ssiASV = rASV + rY

return(data.frame(GEN = mean[1],
                  Y = mean[2],
                  rY = rY,
                  ASV = ASV,
                  rASV = rASV,
                  ssiASV = ssiASV,
                  SIPC = SIPC,
                  rSIPC = rS,
                  ssiSIPC = ssiSIPC,
                  EV = EV,
                  rEV = rEV,
                  ssiEV))
}

