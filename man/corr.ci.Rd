\name{corr.ci}
\alias{corr.ci}
\title{Nonparametric confidence interval for Pearson's correlation coefficient}
\usage{
corr.ci(data = NA, r = NULL, n = NULL)
}
\description{
Computes the half-width confidence interval for correlation coefficient based on the nonparametri method proposed by Olivoto et al. (2018).
}

\arguments{
 \item{data}{A dataset containing variables only or a symetric correlation matrix.}

  \item{r}{If \code{data} is not available, provide the value for correlation coefficient.}

   \item{n}{The sample size if \code{data} is a correlation matrix or if r is informed.}
 }

 \details{The half-width confidence interval is computed according to the following equation:
 \deqn{CI_w = 0.45304^r \times 2.25152 \times n^{-0.50089}}

 where \eqn{n} is the sample size and \code{r} is the correlation coefficient.

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\references{
Olivoto, T., A.D.C. Lucio, V.Q. Souza, M. Nardino, M.I. Diel, B.G. Sari, D.. K. Krysczun, D. Meira, and C. Meier. 2018. Confidence interval width for Pearson's correlation coefficient: a Gaussian-independent estimator based on sample size and strength of association. Agron. J. 110:1-8. \href{https://dl.sciencesocieties.org/publications/aj/abstracts/0/0/agronj2017.09.0566}{doi:10.2134/agronj2017.09.0566}.

}

\examples{

corr.ci(iris[,1:3])

}

