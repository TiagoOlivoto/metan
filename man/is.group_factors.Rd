\name{is.group_factors}
\alias{is.group_factors}
\alias{as.group_factors}
\title{Coerce to an object of class group_factors}
\usage{
\method{as}{group_factors}(x, verbose = TRUE)

\method{is}{group_factors}(x)

}
\description{
Functions to check if an object is of class \code{group_factors}, or coerce it if possible.
}

\arguments{
 \item{x}{The input data.}

  \item{verbose}{Logical argument. If \code{verbose = FALSE} the code will run silently.}
}

\details{
A dataframe may be easily coerced to be split into named subsets based on each combination of factors existing in the original dataframe. For example, if the original data has two columns, namely ENV (four levels) and HIB (ten levels), and ten numeric columns, then using \code{as.group_factors} will split the data into 40 10-columns subsets, corresponding to each combination of ENV x HIB.
}
\examples{
library(METAAB)
library(dplyr)
spdata = as.group_factors(iris)

spdata2 = as.group_factors(CO2)

is.group_factors(spdata2)
}

