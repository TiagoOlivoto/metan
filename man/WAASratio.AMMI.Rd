\name{WAASratio.AMMI}
\alias{WAASratio.AMMI}
\title{Weighted Average of Absolute Scores for AMMI analysis in different scenarios of WAAS/GY ratio}
\usage{
WAASratio.AMMI(data, resp, p.valuePC = 0.05, increment = 5,saveWAASY = 50,
progbar = TRUE)
}

\description{
Compute the Weighted Average of Absolute Scores for AMMI analysis in different combinations of weights for stability and productivity.

}
\arguments{
 \item{data}{The dataset. The imput must be in the following order: Environment, Genotype, replication/block and response variable(s)}

  \item{resp}{The response variable, for example \code{resp = "RG"}.}

  \item{p.valuePC}{The p-value for considering the PC significant. Default is 0.05. The number of significant Principal Components to be used for calculating the WAASB will be chosen based on this probability. }

  \item{increment}{The range of the increment for WAAS/GY ratio. Default is 5. The function compute the WAASY  values starting with a weight o 100 for stability and 0 for response variable. With the default, the first scenario will be a WAAS/GY ratio = 100/0. In the next scenario, the WAASY values are computed based on a WAAS/GY ratio = 95/5.}

  \item{saveWAASY}{Automatically save the WAASY values when the wheight for WAAS (stability) in the WAAS/GY ratio is "saveWAASY". Default is 100. The value of "saveWAASY" must be multiple of "Increment". If this assumption is not valid, an error will be occour.}

  \item{progbar}{A logical argument to define if a progress bar is shown. Default is \code{TRUE}.}

}

\details{
This function is very similar to the \code{WAASBYratio}. The main difference is that here, the WAASBY values are computed considering a traditional AMMI model.
}


\value{
\item{anova}{Joint analysis of variance for the main effects and Principal Component analysis of the interaction effect.}

\item{PC}{Principal Component Analysis.}


\item{MeansGxE}{The means of genotypes in the environments, with observed, predicted and residual values.}

\item{WAAS}{A data frame with the response variable, the scores of all Principal Components, the estimates of Weighted Average of Absolute Scores, and WAASY (the index that consider the weights for stability and productivity in the genotype ranking.}

\item{WAASY}{The values of the WAASY estimated when the wheight for the stability in the loop match with argument "saveWAASY".}


\item{WAASY.values}{All the values of WAASY estimated in the different scenarios of WAAS/GY weighting ratio.}

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\examples{
# Default, with increment of 5 and saving the WAASY values when weight is 50
#library(WAASB)
#model= WAASratio.AMMI(WAASBdata, resp = "GY")

# Incrementing 2-by-2 and saving the WAASY values when weight is 70
#library(WAASB)
#model= WAASratio.AMMI(WAASBdata, resp = "GY", increment = 2, saveWAASY = 70)
}
