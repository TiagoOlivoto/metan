\name{predict.AMMI}
\alias{predict.AMMI}
\title{Predict means in AMMI analysis}
\usage{
\method{predict}{AMMI}(object, resp, naxis, ...)
}

\description{
Predict the means of i-genotypes in j-environment considering a specific number of axis.

}
\arguments{

  \item{object}{The data set. Must be in the following order. Environment, Genotype, Block/Rep, and response variable(s). }

  \item{resp}{The response variable, for example \code{resp = "RG"}.}

  \item{naxis}{The the number of axis to be use in the prediction. }

  \item{...}{Additional parameter for the function}
  }

\details{
This function is used to predict the response variable of a two-way table (for examples the yielding of the i-th genotype in the j-th environment) based on AMMI model. This prediction is based on the number of multiplicative terms used. If \code{naxis = 0}, only the main effects (AMMI0) are used. In this case, the predicted mean will be the predicted value from OLS estimation. If \code{naxis = 1} the AMMI1 (with one multiplicative term) is used for predicting the response variable. If \code{naxis = min(gen-1;env-1)}, the AMMIF is fitted and the predicted value will be the cell mean, i.e. the mean of R-replicates of the i-th genotype in the j-th environment. The number of axis to be used must be carrefully chosen. Precures based on Postdictive sucess (such as Gollobs's d.f.) or Predictive sucess (such as cross-validation) should be used to do this. This package provide both. \code{\link{WAAS.AMMI}} function compute traditional AMMI analysis showing the number of significant axis. On the other hand, \code{\link{validation.AMMIF}} function provide a cross-validation, estimating the RMSE of all AMMI-family models, based on resampling procedure. }

\examples{
library(WAASB)
predict = predict.AMMI(WAASBdata, resp = "GY", naxis = 3)
}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\seealso{\code{\link{predict.AMMImean}}}
