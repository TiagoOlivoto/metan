\name{validation.blup}
\alias{validation.blup}
\title{Cross-validation for blup prediction }
\usage{
validation.blup(data, resp, gen, env, rep,
                nboot, nrepval, progbar = TRUE)
}
\description{
Cross-validation for blup prediction.

}

\arguments{

  \item{data}{The dataset containing the columns related to Environments, Genotypes, replication/block and response variable(s).}

  \item{resp}{The response variable, for example \code{resp = GY}.}

  \item{gen}{The name of the column that contains the levels of the genotypes.}

  \item{env}{The name of the column that contains the levels of the environments.}

  \item{rep}{The name of the column that contains the levels of the replications/blocks.}

  \item{nboot}{The number of resamples to be used in the cross-validation }

  \item{nrepval}{The number of replicates (r) from total number of replicates (R) to be used in the modeling dataset. Only one replicate is used as validating data each step, so, \code{Nrepval} must be equal \code{R-1}}

  \item{progbar}{A logical argument to define if a progress bar is shown. Default is \code{TRUE}.}
  }


\details{
This function provides a  cross-validation procedure for mixed models using replicate-based data. By default, complete blocks are randomly selected within each evironment. In each iteraction, the original dataset is split up into two datasets: training and validation data. The "training" set has all combinations (genotype x environment) with the number of replications informed in \code{nrepval}. The "validation" set has the remaining replication. The estimated values are compared with the "validation" data and the Root Means Square Prediction Difference is computed. At the end of boots, a list is returned.}

\examples{
library(METAAB)
model = validation.blup(data_ge,
                        resp = GY,
                        gen = GEN,
                        env = ENV,
                        rep = REP,
                        nboot = 100,
                        nrepval = 2)

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\seealso{\code{\link{plot.scores}}, \code{\link{plot.WAASBY}}}
