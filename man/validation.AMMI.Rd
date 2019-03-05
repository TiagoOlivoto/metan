\name{validation.AMMI}
\alias{validation.AMMI}
\title{Cross-validation for estimation of AMMI model}
\usage{
validation.AMMI(data, resp, gen, env, rep, design = "RCBD",
                nboot, nrepval, naxis, progbar = TRUE)
}

\description{
Cross-validation for estimation of AMMI models

}
\arguments{

 \item{data}{The dataset containing the columns related to Environments, Genotypes, replication/block and response variable(s).}

  \item{resp}{The response variable, for example \code{resp = "RG"}.}

  \item{gen}{The name of the column that contains the levels of the genotypes.}

  \item{env}{The name of the column that contains the levels of the environments.}

  \item{rep}{The name of the column that contains the levels of the replications/blocks.}

  \item{design}{The experimental desig to be considered. Default is \code{RCBD} (Randomized complete Block Design). For Completely Randomized Designs inform \code{design = "CRD"}.}

  \item{nboot}{The number of resamples to be used in the cross-validation }

  \item{nrepval}{The number of replicates (r) from total number of replicates (R) to be used in the modeling dataset. Only one replicate is used as validating data each step, so, \code{Nrepval} must be equal \code{R-1} }

\item{naxis}{The number of axis to be considered for estimation of GE effects.}

\item{progbar}{A logical argument to define if a progress bar is shown. Default is \code{TRUE}.}
}


\details{
For each iteration, the original dataset is split into two datasets: modeling and validating data. The dataset "modeling" has all combinations (genotype x environment) with the number of replications informed in \code{nrepval}. The dataset "validating" has one replication. The splitting of the dataset into modeling and validating data depends on the design informed. For Completely Randomized Block Design (default), compltely blocks are selected within environments. The remained block serves validation data. If \code{design = "RCD"} is informed, completely randomly samples are made for each genotype-by-environment combination. The estimated values (depending on NAXIS informed) are compared with the "validating" data. the Root Means Square error is computed. At the end of boots, a list is returned with the following values.}


 \value{
\item{RMSE}{A vector with Nboot-estimates of the root mean squared error estimated with the difference between predicted and validating data.}

\item{RSMEmean}{The mean of RMSE estimates.}

\item{Estimated}{A data frame that contain the values (predicted, observed, validation) of the last loop.}

\item{Modeling}{The dataset used as modeling data in the last loop.}

\item{Testing}{The dataset used as testing data in the last loop.}

}

\examples{
\dontrun{
library(METAAB)
model = validation.AMMI(data_ge,
                        resp = GY,
                        gen = GEN,
                        env = ENV,
                        rep = REP,
                        nboot = 100,
                        nrepval = 2,
                        naxis = 2)
}
}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\seealso{\code{\link{validation.AMMIF}}}
