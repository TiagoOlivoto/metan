\name{predict.AMMImean}
\alias{predict.AMMImean}
\title{Predict means in AMMI analysis with data from single observation}
\usage{
\method{predict}{AMMImean}(data, resp, naxis)
}

\description{
Predict the means of i-genotypes in j-environment considering a specific number of axis.

}
\arguments{

  \item{data}{The data set. Must be in the following order. Environment, Genotype, and response variable(s). }

    \item{resp}{The response variable, for example \code{resp = "RG"}.}

  \item{naxis}{The the number of axis to be use in the prediction. If \code{NAXIS = 0}, the predicted mean is the predicted value of OLS estimation. If \code{NAXIS} = "min(gen-1;env-1)", the predicted value will be the cell mean, i.e. the mean of R-replicates of the i-th genotype in the j-th environment }
  }


\examples{
library(WAASB)
predmean = predict.AMMImean(datameans, resp = "Resp", naxis = 3)

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\seealso{\code{\link{predict.AMMI}}}
