\name{MTSI}
\alias{MTSI}
\title{Multi-trait stability index}
\usage{
MTSI(.data, index = "WAASB", SI = 15, mineval = 1, verbose = TRUE)
}

\description{
Multi-trait stability index in analysis of multi-environment trials using AMMI or BLUP methods.


}
\arguments{

  \item{.data}{An object of class \code{WAASB} or \code{WAAS.AMMI}.  }

  \item{index}{If \code{index = "WAASB"} (Default) the multitrait index will be computed considering the stability of genotypes only. If \code{index = "WAASBY"} both stability and mean performance are considered. More details can be seen in \code{\link{WAASB}} and \code{\link{WAAS.AMMI}} functions.}

  \item{SI}{An integer [0-100]. The selection intensity in percentage of the total number of genotypes.}

  \item{mineval}{The minimum value so that an eigenvector is retained in the factor analysis.}

  \item{verbose}{If \code{verbose = TRUE} (Default) then some results are shown in the console.}

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\examples{
library(METAAB)
library(dplyr)

# Based on stability only, for both GY and HM, higher is better
MTSI_MODEL = WAASB(data_ge,
                   env = ENV,
                   gen = GEN,
                   rep = REP.
                   resp = c(GY, HM))
MTSI_index = MTSI(MTSI_MODEL)


# Based on mean performance and stability (using pipe operator \%>\%)
# GY: higher is better
# HM: lower is better
MTSI_index2 = data_ge \%>\%
              WAASB(ENV, GEN, REP, c(GY, HM), mresp = c(100, 0)) \%>\%
              MTSI(index = "WAASBY")
}
