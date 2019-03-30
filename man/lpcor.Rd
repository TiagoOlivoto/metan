\name{lpcor}
\alias{lpcor}
\title{Linear and Partial Correlation Coefficients}
\usage{
lpcor(.data, group_var = NULL, n = NULL, verbose = TRUE)
}
\description{
Estimates the linear and partial correlation coefficients using as input a data frame or a correlation matrix.

}

\arguments{
\item{.data}{The data to be analyzed. Must be a symmetric correlation matrix or a dataframe containing continuous variables and one possible grouping variable.}

 \item{group_var}{A possible factor to group the data. If more than two factors are in the dataset, the extra factor column(s) will be excluded with a message by default.}

  \item{n}{If a correlation matrix is provided, then \code{n} is the number of objects used to compute the correlation coefficients.}

   \item{verbose}{If \code{verbose = TRUE} then some results are shown in the console.}
 }

\value{
If a grouping factor is used then a list is returned with the following values.
\item{linear.mat}{The matrix of linear correlation.}
\item{partial.mat}{The matrix of partial correlations.}
\item{results}{Hypotesis testing for each pairwise comparation.}
}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\examples{
library(METAAB)
library(dplyr)
partial1 = lpcor(iris[,1:4])

# Alternatively using the pipe operator \%>\%
partial2 = iris \%>\%
           select(1:4) \%>\%
           lpcor()

partial3 = cor(iris[1:4]) \%>\%
           lpcor(n = nrow(iris),
                 verbose = FALSE)

partial4 = iris \%>\%
           lpcor(group_var = Species)
}


