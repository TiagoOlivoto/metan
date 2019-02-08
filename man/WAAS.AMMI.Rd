\name{WAAS.AMMI}
\alias{WAAS.AMMI}
\title{Weighted Average of Absolute Scores for AMMI analysis}
\usage{
WAAS.AMMI(data, resp, gen, env, rep, maxresp = NULL,
      minresp = NULL, p.valuePC = 0.05, naxis = NULL,
      weight.response = 50, weight.WAAS = 50)
}
\description{
Compute the Weighted Average of Absolute Scores for AMMI analysis.

}

\arguments{
  \item{data}{The dataset containing the columns related to Environments, Genotypes, replication/block and response variable(s).}

 \item{resp}{The response variable, for example \code{resp = GY}.}

  \item{gen}{The name of the column that contains the levels of the genotypes.}

  \item{env}{The name of the column that contains the levels of the environments.}

  \item{rep}{The name of the column that contains the levels of the replications/blocks.}

  \item{maxresp, minresp}{A numeric value or a numeric vector of length equal the number of variables in the argument \code{resp}. \code{maxresp} and \code{minresp} values will be the new maximum and minimum values, respectively, of the original response variable. By default, all variables in \code{resp} are rescaled so that de maximum value is 100 and the minimum value is 0.}

  \item{p.valuePC}{The p-value for considering the PC significant. The number of significant Principal Components to be used for calculating the WAAS will be chosen based on this probability. }

  \item{naxis}{The number of axis to be used for computing the WAAS value. Defaul it \code{NULL}. If informed some values, the number of axis informed is used for computing the WAAS.}

  \item{weight.response}{The weight for the response variable.}

  \item{weight.WAAS}{The weight for the WAAS (stability).}

}

\details{
This function compute the weighted average of absolute scores considering the significant principal component axis. The number of axis to be used depends on the probability level informed in \code{p.valuePC}. For example if \code{p.valuePC = 0.05}, all axis that are significant considering this probability level are used. The number of axis can be also informed by declaring \code{naxis = x}. This comand ignores the \code{p.valuePC} comand.
}

\value{

\item{individual}{A within-environments ANOVA considering a fixed-effect model.}

\item{model}{A data frame with the response variable, the scores of all Principal Components, the estimates of Weighted Average of Absolute Scores, and WAASY (the index that consider the weights for stability and productivity in the genotype ranking.}

\item{MeansGxE}{The means of genotypes in the environments, with observed, predicted and residual values.}

\item{PCA}{Principal Component Analysis.}


\item{anova}{Joint analysis of variance for the main effects and Principal Component analysis of the interaction effect.}

\item{Details}{A list summarizing the results. The following information are showed. \code{WgtResponse}, the weight for the response variable in estimating WAASB, \code{WgtWAAS} the weight for stability, \code{Ngen} the number of genotypes, \code{Nenv} the number of environments, \code{OVmean} the overall mean, \code{Min} the minimum observed (returning the genotype and environment), \code{Max} the maximum observed, \code{Max} the maximum observed, \code{MinENV} the environment with the lower mean, \code{MaxENV} the environment with the larger mean observed, \code{MinGEN} the genotype with the lower mean,  \code{MaxGEN} the genotype with the larger.}

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\examples{
# Considering p-value <= 0.05 to compute the WAAS
library(METAAB)
model = WAAS.AMMI(WAASBdata, GY, GEN, ENV, REP,
                  weight.response = 40, weight.WAAS = 60)


# Declaring the number of axis to be used for computing WAAS

model = WAAS.AMMI(WAASBdata, resp = GY, gen = GEN, env = ENV, rep = REP,
                  naxis = 3, weight.response = 40, weight.WAAS = 60)

}


