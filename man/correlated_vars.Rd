% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlated_vars.R
\name{correlated_vars}
\alias{correlated_vars}
\title{Generate correlated variables}
\usage{
correlated_vars(
  y,
  min_cor = -1,
  max_cor = 1,
  nvars,
  constant = NULL,
  operation = "*",
  x = NULL
)
}
\arguments{
\item{y}{A vector to generate variables correlated with.}

\item{min_cor}{The minimum desired correlation.}

\item{max_cor}{The maximum desired correlation.}

\item{nvars}{The number of variables.}

\item{constant}{A constant. Use \code{operation} to define which operation is
used.}

\item{operation}{The operation to be applied to the \code{constant} value.}

\item{x}{An optional vector of the same length of \code{y}. If not informed
(default) then a normally distributed variable (mean = 0, sd = 1) will be
used.}
}
\value{
A data frame with the \code{y} variable and the correlated variables.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}

Generate correlated variables using a vector of know values and desired
maximum and minimum correlations
}
\examples{
\donttest{
library(metan)
y <- rnorm(n = 10)
cor_vars <- correlated_vars(y, nvar = 6)
plot(cor_vars)
}


}
\author{
Tiago Olivoto \email{tiagoolivoto@gmail.com}
}
