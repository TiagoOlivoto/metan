predict.AMMImean = function(object,
                            resp,
                            naxis,
                            ...){
Y = object[paste(resp)]
data = as.data.frame(object[,1:2])
data = cbind(data, Y)
names(data) = c("ENV", "GEN", "Y")
data$ENV = as.factor(data$ENV)
data$GEN = as.factor(data$GEN)
x1  = model.matrix(~factor(data$ENV) -1)
z1  = model.matrix(~factor(data$GEN) -1)
Nenv = ncol(x1)
Ngen = ncol(z1)
min = min(Nenv, Ngen) - 1

if (naxis > min){
  stop("The number of axis to be used must be lesser than or equal to min(GEN-1;ENV-1), in this case, ", min,".")
} else{
model = lm(Y ~ ENV + GEN, data = data)
MODEL = dplyr::mutate(data,
                      ResOLS = model$residuals,
                      PredOLS = Y - ResOLS)
INT = t(matrix(MODEL$ResOLS, Nenv, byrow = T))
s = svd(INT)
U = s$u[,1:naxis]
LL = s$d[1:naxis]
V = s$v[,1:naxis]
AMMI = ((z1 %*% U)*(x1 %*% V)) %*% LL
MODEL = dplyr::mutate(MODEL,
                      ResAMMI = AMMI,
                      PredAMMI = PredOLS + ResAMMI)
return(structure(MODEL), class = "predict.AMMImean")
 }
}


