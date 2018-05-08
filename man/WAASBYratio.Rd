\name{WAASBYratio}
\alias{WAASBYratio}
\title{Weighted Average of Absolute Scores from SVD of BLUP-interaction effects matrix in different scenarios of WAASB/GY ratio}
\usage{
WAASBYratio(data, resp, increment = 5, saveWAASY = 50, progbar = TRUE)
}

\description{
Compute the Weighted Average of Absolute Scores for AMMI analysis in different combinations of weights for stability and productivity.

}
\arguments{
 \item{data}{The dataset. The imput must be in the following order: Environment, Genotype, replication/block and response variable(s)}

  \item{resp}{The response variable, for example \code{resp = "RG"}.}

  \item{increment}{The range of the increment for WAASB/GY ratio. Default is 5. The function compute the WAASBY  values starting with a weight o 100 for stability and 0 for response variable. With the default, the first scenario will be a WAASB/GY ratio = 100/0. In the next scenario, the WAASBY values are computed based on a WAASB/GY ratio = 95/5.}

  \item{saveWAASY}{Automatically save the WAASY values when the wheight for WAAS (stability) in the WAAS/GY ratio is "saveWAASY". Default is 100. The value of "saveWAASY" must be multiple of "Increment". If this assumption is not valid, an error will be occour.}

   \item{progbar}{A logical argument to define if a progress bar is shown. Default is \code{TRUE}.}
}

\details{
This function considers both stability (weighted average of absolute scores based on SVD of BLUP-interaction effects matrix) and productivitye for genotype ranking. This function provide the option of attributing weights for stability and productive in genotype ranking. This is important depending on the goal of a selection strategy. For example, if a a goal of a breeding program is to select a genotype whith high yielding (independeltly on the stability performance), that genotype with the first rank in an WAASB/GY = 0/100 ratio should be selected. The reciprocal is true. Aiming at selecting a high-stable genotype (independentely on the productivity), that genotype with the first rank in an WAASB/GY = 100/0 ratio should be selected. By defalut, the increment on the WAASB/GY ratio is equal to 5. In other words, twenty one different combinations are computed. Each combination, the genotypes are ranked regarding the WAASY value.
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
library(WAASB)
model= WAASBYratio(WAASBdata, resp = "GY")

# Incrementing 2-by-2 and saving the WAASY values when weight is 70
library(WAASB)
model= WAASBYratio(WAASBdata, resp = "GY", increment = 2, saveWAASY = 70)
}
