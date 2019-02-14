\name{MTSI}
\alias{MTSI}
\title{Multitrait stability index}
\usage{
MTSI(x, index = "WAASB", show = TRUE, SI = 15, mineval = 1)
}

\description{
Multitrait stability index in analysis of multi-environment trials using AMMI or BLUP methods.


}
\arguments{

  \item{x}{An object of class \code{WAASB} or \code{WAAS.AMMI}.  }

  \item{index}{If \code{index = "WAASB"} (Default) the multitrait index will be computed considering the stability of genotypes only. If \code{index = "WAASBY"} both stability and mean performance are considered. More details can be seen in \code{\link{WAASB}} and \code{\link{WAAS.AMMI}} functions.}

  \item{show}{If \code{show = TRUE} (Default) some results are shown in the console.}

    \item{SI}{An integer [0-100]. The selection intensity in percentage of the total number of genotypes.}

  \item{mineval}{The minimum value so that an eigenvector is retained in the factor analysis.}

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\examples{
library(METAAB)
# Based on stability only
MTSI_MODEL = WAASB(data_ge,
                   resp = c(GY, HM),
                   gen = GEN,
                   env = ENV,
                   rep = REP)

MTSI_index = MTSI(MTSI_MODEL)


# Based on mean performance and stability
MTSI_MODEL2 = WAASB(data_ge,
                    resp = c(GY, HM),
                    gen = GEN,
                    env = ENV,
                    rep = REP,
                    mresp = c(100, 100), #Default
                    wresp = c(65, 65))
MTSI_index2 = MTSI(MTSI_MODEL2,
                   index = "WAASBY")

}
