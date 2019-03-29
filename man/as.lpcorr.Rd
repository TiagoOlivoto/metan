\name{as.lpcor}
\alias{as.lpcor}
\alias{is.lpcor}
\title{Coerce to an object of class lpcor}
\usage{
\method{as}{lpcor}(x, y, ...)

\method{is}{lpcor}(x)

}
\description{
Functions to check if an object is of class \code{lpcor}, or coerce it if possible.
}

\arguments{
 \item{x}{The first object.}

 \item{y}{The second object.}

  \item{...}{Additional matrices to be coerced to a list.}

}


\examples{
library(METAAB)
library(dplyr)
mt_num = mtcars \%>\% select_if(., is.numeric)
lpdata = as.lpcor(cor(mt_num[1:5]),
                  cor(mt_num[1:5]),
                  cor(mt_num[2:6]),
                  cor(mt_num[4:8]))
is.lpcor(lpdata)

}

