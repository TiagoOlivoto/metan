\name{make_sym}
\alias{make_sym}
\title{Make a symetric matrix based on an upper or lower triangular matrix}
\usage{
make_sym(.matrix, make = "upper")
}
\description{
This function help users to easily make a symmetric matrix using a lower or an upper triangular matrix.
}

\arguments{
 \item{.matrix}{The upper or lower triangular matrix.}

 \item{make}{The triangular to built. Default is \code{"upper"}. In this case, a symetric matrix will be built based on the values of a lower triangular values.}


 }

\examples{
m <- matrix(NA,4,4)
m[lower.tri(m)] <- 1:16
m
make_sym(m)
}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


