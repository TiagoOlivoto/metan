\name{partial.corr}
\alias{partial.corr}
\title{Partial correlation coefficient}
\usage{
partial.corr(data, n = NULL, method = "pearson")
}
\description{
Estimates of partial correlation coefficient based on data raw or correlation matrices.

}

\arguments{
 \item{data}{The data, must be a dataframe or a symetric correlation matrix.}

 \item{n}{The number of subjects used if \code{data} is a correlation matrix.}

 \item{method}{If \code{data} is a correlation matrix, then code{method} is the method for computing the correlation; one of the "pearson" (default), "kendall", or "spearman" }

 }

\value{
\item{linear.mat}{The matrix of linear correlation.}
\item{partial.mat}{The matrix of partial correlations.}
\item{results}{Hypotesis testing for each pairwise comparation}
\item{call}{The function call}
}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\examples{
partial.corr(iris[,1:4])
}


