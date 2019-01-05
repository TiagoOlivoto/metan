\name{WAASB}
\alias{WAASB}
\title{Weighted Average of Absolute Scores for BLUP-based PCA analysis}
\usage{
WAASB(data, resp, gen, env, rep, random = "gen",
      prob = 0.95, weight.response = 50, weight.WAAS = 50)
}
\description{
Compute the Weighted Average of Absolute Scores for PCA analysis based on the GE-interaction matrix.

}

\arguments{

 \item{data}{The dataset containing the columns related to Environments, Genotypes, replication/block and response variable(s).}

  \item{resp}{The response variable, for example \code{resp = RG}.}

  \item{gen}{The name of the column that contains the levels of the genotypes.}

  \item{env}{The name of the column that contains the levels of the environments.}

  \item{rep}{The name of the column that contains the levels of the replications/blocks.}

  \item{random}{The effect to be considered as random. Default is genotype "gen". A random effects model can be fitted using \code{random = "all"}.}

  \item{prob}{The probability for estimating confidence interval for BLUP's prediction. }

  \item{weight.response}{The weight for the response variable. Defatul is 50.}

  \item{weight.WAAS}{The weight for the WAAS (stability). Default is 50.}

}

\details{
This function compute the weighted average of absolute scores considering all principal component axis from Singular Value Decomposition of BLUP-interaction effects matrix. As the matrix of genotypic effects (free from error and interaction effects) is used, all multiplicative terms can be used without a noise-adding  effetc. The main advantage of this procedure in relation to the \code{WAAS.AMMI} function is that, by using mixed model-based analysis, random effects can be included in the model. In addition, unbalanced datasets can also be modeled.
}



\value{

\item{individual}{A within-environments ANOVA considering a fixed-effect model.}

\item{model}{A data frame with the response variable, the scores of all Principal Components, the estimates of Weighted Average of Absolute Scores, and WAASY (the index that consider the weights for stability and productivity in the genotype ranking.}

\item{BLUPgen}{The estimated BLUPS for genotypes.}

\item{BLUPenv}{The estimated BLUPS for environments, if \code{random = "all"}.}

\item{BLUPge}{The estimated BLUPS of all genotypes in all environments "BLUPij".}

\item{PCA}{The results of Principal Component Analysis with eigenvalues and explained variance of BLUP-interaction matrix.}

\item{MeansGxE}{The phenotypic means of genotypes in the environments, with observed, predicted and OLS residual prediction.}


\item{Details}{A list summarizing the results. The following information are showed. \code{WgtResponse}, the weight for the response variable in estimating WAASB, \code{WgtWAAS} the weight for stability, \code{Ngen} the number of genotypes, \code{Nenv} the number of environments, \code{OVmean} the overall mean, \code{Min} the minimum observed (returning the genotype and environment), \code{Max} the maximum observed, \code{Max} the maximum observed, \code{MinENV} the environment with the lower mean, \code{MaxENV} the environment with the larger mean observed, \code{MinGEN} the genotype with the lower mean,  \code{MaxGEN} the genotype with the larger.}

\item{REML}{A list with the following data frames: \code{fixed}, that contains the ANOVA table with F-tests and p-values for fixed effects (If \code{random = "gen"}); and \code{random} that contains the estimated variance components of the random effects.}

\item{ESTIMATES}{A list with the following values: \code{GEV} the genotype-by-environment variance (and percentage of phenotypic variance); \code{GV} the genotypic variance (and percentage of phenotypic variance); \code{EV} the environmental variance;\code{RV} the residual variance (and percentage of phenotypic variance); \code{FV} the phenotypic variance;  \code{h2g} the heritability of the trait; \code{GEr2} the coefficient of determination of the interaction effects; \code{h2mg} the heritability of the mean; \code{AccuGen} the selective accuracy; \code{rge} the genotype-environment correlation; \code{CVg} the genotypic coefficient of variation; \code{CVr} the residual coefficient of variation; \code{CVratio} the ratio between genotypic and residual coefficient of variation.}


\item{LRT}{The Likelihood Ratio Test for the random effects of the model.}
}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\examples{
# Ggenotypes as random effects and equal weights for both response variable and stability

library(METAAB)
model = WAASB(WAASBdata, GY, GEN, ENV, REP)

# Higher weight for stability
model2 = WAASB(WAASBdata, rep = REP, env = ENV, resp = GY, gen = GEN,
               weight.response = 70, weight.WAAS = 30)

# A random-effect model
model3 = WAASB(WAASBdata, GY, GEN, ENV, REP, random = "all")

}


