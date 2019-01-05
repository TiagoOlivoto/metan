\name{predict.WAAS.AMMI}
\alias{predict.WAAS.AMMI}
\title{Predict the means of a WAAS.AMMI object}
\usage{
\method{predict}{WAAS.AMMI}(object, naxis, digits = 4, ...)
}

\description{
Predict the means of a WAAS.AMMI object considering a specific number of axis.

}
\arguments{

  \item{object}{An object of class WAAS.AMMI}

  \item{naxis}{The the number of axis to be use in the prediction. }

  \item{digits}{The significant digits to be shown.}

  \item{...}{Additional parameter for the function}
  }

\details{
This function is used to predict the response variable of a two-way table (for examples the yielding of the i-th genotype in the j-th environment) based on AMMI model. This prediction is based on the number of multiplicative terms used. If \code{naxis = 0}, only the main effects (AMMI0) are used. In this case, the predicted mean will be the predicted value from OLS estimation. If \code{naxis = 1} the AMMI1 (with one multiplicative term) is used for predicting the response variable. If \code{naxis = min(gen-1;env-1)}, the AMMIF is fitted and the predicted value will be the cell mean, i.e. the mean of R-replicates of the i-th genotype in the j-th environment. The number of axis to be used must be carrefully chosen. Precures based on Postdictive sucess (such as Gollobs's d.f.) or Predictive sucess (such as cross-validation) should be used to do this. This package provide both. \code{\link{WAAS.AMMI}} function compute traditional AMMI analysis showing the number of significant axis. On the other hand, \code{\link{validation.AMMIF}} function provide a cross-validation, estimating the RMSE of all AMMI-family models, based on resampling procedure. }

\examples{
library(METAAB)
model = WAAS.AMMI(WAASBdata, GY, GEN, ENV, REP)
predict = predict.AMMI(model, naxis = 3)
}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

