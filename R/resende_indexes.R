Resende_indexes =  function(x){
  if (!is(x, "WAASB")) {
    stop("The object 'x' must be an object of class \"WAASB\"")
  }

gge = x$BLUPgge
# Helper functions
hmean_fun = function(x){
hmean = length(x)/sum(1/x)
return(hmean)
}
make_matrix = function(data, row, col, value){
  nam = cbind(c(row, col, value))
  data = data.frame(data[(match(c(nam), names(data)))])
  return(data.frame(tapply(data[, 3], data[, c(1, 2)], mean)))
}
# Harmonic mean
GEPRED = make_matrix(gge, "GEN", "ENV", "Predicted")
HMGV = data.frame(HMGV = apply(GEPRED, 1, FUN = hmean_fun),
                  HMGV_order = rank(-apply(GEPRED, 1, FUN = hmean_fun)))
## Relative performance
y = x$MeansGxE
GEMEAN = make_matrix(y, "GEN", "ENV", "Y")
ovmean = mean(y$Y)
mean_env = apply(GEMEAN, 2, FUN = mean)
RPGV = data.frame(RPGV = apply(t(t(GEPRED) / mean_env), 1, mean))
RPGV_data = dplyr::mutate(RPGV,
                          GY_RPGV = RPGV * ovmean,
                          RPGV_order = rank(-GY_RPGV))
## Harmonic mean of Relative performance
HMRPGV = data.frame(HMRPGV = apply(t(t(GEPRED) / mean_env), 1, hmean_fun))
HMRPGV_data = dplyr::mutate(HMRPGV,
                            GY_HMRPGV = HMRPGV * ovmean,
                            HMRPGV_order = rank(-GY_HMRPGV))
final = data.frame(cbind(HMGV, RPGV_data, HMRPGV_data))
rownames(final) = rownames(RPGV)
return(final)
}
