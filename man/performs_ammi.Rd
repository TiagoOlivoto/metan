\name{performs_ammi}
\alias{performs_ammi}
\title{Performs Additive Main effects and Multiplicative Interaction}
\usage{
performs_ammi(ENV, GEN, REP, Y)
}
\description{
Compute the Additive Main effects and Multiplicative interaction. This is a helper function for other procedures performed in the WAASB package such as \code{\link{WAAS.AMMI}}, \code{\link{validation.AMMI}}, \code{\link{validation.AMMIF}}, and \code{\link{predict.WAAS.AMMI}}.

}

\arguments{
  \item{ENV}{The name of the collum that contains the levels of the environments}

 \item{GEN}{The name of the collum that contains the levels of the genotypes}

\item{REP}{The name of the collum that contains the levels of the replications/blocks}

\item{Y}{The response variable}

}


\value{

\item{ANOVA}{The analysis of variance for the AMMI model.}

\item{analysis}{The principal component analysis}

\item{means}{means of genotype vs environment}

\item{biplot}{escores for genotypes and environments in all the possible axes.}

}


\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\examples{
library(METAAB)
ammi_model = with(data_ge, performs_ammi(ENV, GEN, REP, GY))

}


