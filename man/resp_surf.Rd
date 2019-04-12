\name{resp_surf}
\alias{resp_surf}
\title{Response surface model}
\usage{
resp_surf(.data, factor1, factor2, rep,
          resp, prob = 0.05, verbose = TRUE)
}
\description{
Compute a surface model and find the best combination of factor1 and factor2 to obtain the stacionary point.

}

\arguments{

 \item{.data}{The dataset containing the columns related to Environments, factor1, factor2, replication/block and response variable(s).}

  \item{factor1}{The first factor, for example, dose of Nitrogen.}

  \item{factor2}{The second factor, for example, dose of potassium.}

  \item{rep}{The name of the column that contains the levels of the replications/blocks.}

  \item{resp}{The response variable(s).}

  \item{prob}{The probability error.}

  \item{verbose}{If \code{verbose = TRUE} then some results are shown in the console.}
}
\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}



