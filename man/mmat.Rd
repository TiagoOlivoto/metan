\name{mmat}
\alias{mmat}
\title{Make a two-way table based on categorical and numerical arguments}
\usage{
mmat(.data, row, col, value)
}
\description{
This function help users to easely make a two-way table from a "long format" data.

}

\arguments{
 \item{.data}{The dataset. Must contains at least two categorical columns.}

 \item{row}{The column of data in which the mean of each level will correspond to \strong{one line} in the output.}

 \item{col}{The column of data in which the mean of each level will correspond to \strong{one column} in the output.}

 \item{value}{The column of data that contains the values to fill the two-way table.}


 }

\examples{
library(METAAB)
library(dplyr)
matrix = data_ge \%>\% mmat(row = GEN, col = ENV, val = GY)
matrix
}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


