\name{FAI.BLUP}
\alias{FAI.BLUP}
\title{Multitrait index based on factor analysis and ideotype-design}
\usage{
FAI.BLUP(.data, DI, UI, SI = NULL, mineval = 1, verbose = TRUE)
}

\description{
Multitrait index based on factor analysis and ideotype-design proposed by Rocha et al. (2018).


}
\arguments{


  \item{.data}{An object of class \code{WAASB} }

  \item{DI, UI}{A vector of the same length of \code{.data} to construct the desirable (DI) and undesirable (UI) ideotypes. For each element of the vector, allowed values are \code{"max"}, \code{"min"}, \code{"mean"}, or a numeric value.}


  \item{SI}{An integer [0-100]. The selection intensity in percentage of the total number of genotypes.}

  \item{mineval}{The minimum value so that an eigenvector is retained in the factor analysis.}

  \item{verbose}{Logical value. If \code{TRUE} some results are shown in console.}

}


\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\references{
Rocha, J.R.A.S.C.R, J.C. Machado, and P.C.S. Carneiro. 2018. Multitrait index based on factor analysis and ideotype-design: proposal and application on elephant grass breeding for bioenergy. GCB Bioenergy 10:52-60. doi: \href{https://onlinelibrary.wiley.com/doi/full/10.1111/gcbb.12443}{doi:10.1111/gcbb.12443}.
}

\examples{
library(METAAB)
library(dplyr)

mod = WAASB(data_ge,
            resp = c(GY, HM),
            gen = GEN,
            env = ENV,
            rep = REP)

FAI = FAI.BLUP(mod,
               SI = 15,
               DI = c("max", "max"),
               UI = c("min", "min"))

# Or using the pipe operator \%>\%

FAI = data_ge \%>\%
      WAASB(c(GY, HM), GEN, ENV, REP) \%>\%
      FAI.BLUP(DI = c("max", "max"),
               UI = c("min", "min"),
               SI = 15)
plot(FAI)
}
