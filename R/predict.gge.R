#' Predict a two-way table based on GGE model
#'
#' Predict the means for a genotype-vs-environment trial based on a Genotype plus
#' Genotype-vs-Environment interction (GGE) model.
#'
#' This function is used to predict the response variable of a two-way table
#' (for examples the yielding of g genotypes in e environments)
#' based on GGE model. This prediction is based on the number of
#' principal components used. For more detais see Yan and Kang (2007).
#'
#' @param object An object of class \code{gge}.
#' @param naxis The the number of principal components to be used in the prediction.
#' Generaly, two axis may be used. In this case, the estimated values will be those shown
#' in the biplot.
#' @param output The type of output. It must be one of the \code{'long'} (default) returning
#' a long-format table with the columns for environment (ENV), genotypes (GEN) and response
#' variable (Y); or \code{'wide'} to return a two-way table with genotypes in the row, environments
#' in the columns, filled by the estimated values.
#' @param ... Additional parameter for the function
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Yan, W., and M.S. Kang. 2003. GGE biplot analysis: a graphical tool for breeders,
#'  geneticists, and agronomists. CRC Press.
#' @method predict gge
#' @export
#' @examples
#'
#' library(metan)
#' mod = gge(data_ge, GEN, ENV, GY)
#' predict(mod)
#'
predict.gge <- function(object, naxis = 2, output = "wide", ...) {
  if (!class(object) == "gge") {
    stop("The object must be of class 'gge'.")
  }
  if (naxis > min(dim(object$coordenv))) {
    stop("The number of principal components cannot be greater than min(g, e), in this case ",
         min(dim(object$coordenv)))
  }
  # SVP
  if (object$svp == "environment" | object$svp == 2) {
    pred <- (object$coordgen[, 1:naxis] * (object$d)) %*%
      t(object$coordenv[, 1:naxis])
  }
  if (object$svp == "genotype" | object$svp == 1) {
    pred <- (object$coordgen[, 1:naxis] %*% t(object$coordenv[,
                                                              1:naxis] * (object$d)))
  }
  if (object$svp == "symmetrical" | object$svp == 3) {
    pred <- (object$coordgen[, 1:naxis] %*% t(object$coordenv[,
                                                              1:naxis]))
  }
  # Scaling
  if (object$scaling == "sd" | object$scaling == 1) {
    pred <- sweep(pred, 2, object$scale_val, FUN = "*")
  }
  # Centering
  if (object$centering == "global" | object$centering == 1) {
    pred <- pred + object$grand_mean
  }
  if (object$centering == "environment" | object$centering ==
      2) {
    pred <- sweep(pred, 2, object$mean_env, FUN = "+")
  }
  if (object$centering == "double" | object$centering == 3) {
    for (i in 1:nrow(pred)) {
      for (j in 1:ncol(pred)) {
        pred[i, j] <- pred[i, j] - object$grand_mean +
          object$mean_env[j] + object$mean_gen[i]
      }
    }
  }
  rownames(pred) <- object$labelgen
  colnames(pred) <- object$labelenv
  if (output == "wide") {
    return(pred)
  }
  if (output == "long") {
    return(pred %>% as.data.frame() %>% rownames_to_column("GEN") %>%
             gather(-GEN, key = "ENV", value = "Y") %>% select(ENV,
                                                               everything()))
  }
}
