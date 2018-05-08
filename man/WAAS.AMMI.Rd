\name{WAAS.AMMI}
\alias{WAAS.AMMI}
\title{Weighted Average of Absolute Scores for AMMI analysis}
\usage{
WAAS.AMMI(data, resp, p.valuePC = 0.05, naxis = NULL, weight.response = 50, weight.WAAS = 50)
}
\description{
Compute the Weighted Average of Absolute Scores for AMMI analysis.

}

\arguments{
  \item{data}{The data frame. Data must be in the following order: Environment, Genotype, Block, Response variable.}

 \item{resp}{The response variable, for example \code{resp = "RG"}.}

  \item{p.valuePC}{The p-value for considering the PC significant. The number of significant Principal Components to be used for calculating the WAAS will be chosen based on this probability. }

  \item{naxis}{The number of axis to be used for computing the WAAS value. Defaul it \code{NULL}. If informed some values, the number of axis informed is used for computing the WAAS.}

  \item{weight.response}{The weight for the response variable.}

  \item{weight.WAAS}{The weight for the WAAS (stability).}

}

\details{
This function compute the weighted average of absolute scores considering the significant principal component axis. The number of axis to be used depends on the probability level informed in \code{p.valuePC}. For example if \code{p.valuePC = 0.05}, all axis that are significant considering this probability level are used. The number of axis can be also informed by declaring \code{naxis = x}. This comand ignores the \code{p.valuePC} comand.
}

\value{
\item{anova}{Joint analysis of variance for the main effects and Principal Component analysis of the interaction effect.}

\item{PC}{Principal Component Analysis.}


\item{MeansGxE}{The means of genotypes in the environments, with observed, predicted and residual values.}

\item{WAAS}{A data frame with the response variable, the scores of all Principal Components, the estimates of Weighted Average of Absolute Scores, and WAASY (the index that consider the weights for stability and productivity in the genotype ranking.}

\item{WgtResponse}{The weight for the response variable declared. This value is used to create a summary of the analysis.}

\item{WgtWAAS}{The weight for the WAAS (stability) declared. This value is used to create a summary of the analysis.}

\item{Ngen}{The number of genotypes included in the analysis.}

\item{Nenv}{The number of environment included in the analysis.}

\item{NSPC}{The number Principal Components that are the p-value lesser than declared in \\code{pvaluePC}. This number of PC is used to estimate de values of WAAS.}

\item{OVmean}{The overall mean.}

\item{MinENV}{Return both the name and mean of the environment with the lower mean across all environments.}

\item{MaxENV}{Return both the name and mean of the environment with the higher mean across all environments.}

\item{MinGEN}{Return both the name and mean of the genotype with the lower mean across all genotypes.}

\item{MaxGEN}{Return both the name and mean of the genotype with the higher mean across all genotypes.}

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\examples{
# considering p-value <= 0.05 to compute the WAAS
library(WAASB)
model = WAAS.AMMI(WAASBdata, resp = "GY", weight.response = 40, weight.WAAS = 60)


# considering declaring the number of axis to be used for computing WAAS
library(WAASB)
model = WAAS.AMMI(WAASBdata, resp = "GY", naxis = 3, weight.response = 40, weight.WAAS = 60)

}


