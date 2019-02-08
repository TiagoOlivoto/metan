\name{plot.FAI.BLUP}
\alias{plot.FAI.BLUP}
\title{Multitrait index based on factor analysis and ideotype-design}
\usage{
\method{plot}{FAI.BLUP}(x, ideotype = 1, SI = 15, radar = TRUE, ...)
}

\description{
Plot the multitrait index based on factor analysis and ideotype-design proposed by Rocha et al. (2018).


}
\arguments{

  \item{x}{An object of class \code{WAASB} }

  \item{ideotype}{The ideotype to be plotted. Default is 1.}

  \item{SI}{An integer [0-100]. The selection intensity in percentage of the total number of genotypes.}

  \item{radar}{Logical argument. If true (default) a radar plot is generated after using \code{coord_polar()}.}


  \item{...}{Other arguments of the function}

}


\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\references{
Rocha, J.R.A.S.C.R, J.C. Machado, and P.C.S. Carneiro. 2018. Multitrait index based on factor analysis and ideotype-design: proposal and application on elephant grass breeding for bioenergy. GCB Bioenergy 10:52-60. doi: \href{https://onlinelibrary.wiley.com/doi/full/10.1111/gcbb.12443}{doi:10.1111/gcbb.12443}.
}
