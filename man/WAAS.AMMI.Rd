\name{WAAS.AMMI}
\alias{WAAS.AMMI}
\title{Weighted Average of Absolute Scores for AMMI analysis}
\usage{
WAAS.AMMI(data, resp, gen, env, rep, mresp = NULL,
          wresp = NULL, prob = 0.05, naxis = NULL,
          verbose = TRUE)
}
\description{
Compute the Weighted Average of Absolute Scores for AMMI analysis.

}

\arguments{
  \item{data}{The dataset containing the columns related to Environments, Genotypes, replication/block and response variable(s).}

  \item{resp}{The response variable(s). To analyze multiple variables in a single procedure a vector of variables may be used. For example \code{resp = c(var1, var2, var3}.}

  \item{gen}{The name of the column that contains the levels of the genotypes.}

  \item{env}{The name of the column that contains the levels of the environments.}

  \item{rep}{The name of the column that contains the levels of the replications/blocks.}

  \item{mresp}{A numeric vector of the same length of \code{resp}. The \code{mresp} will be the new maximum value after reescaling. By default, all variables in \code{resp} are rescaled so that de maximum value is 100 and the minimum value is 0.}

  \item{wresp}{The weight for the response variable(s) for computing the WAASBY index. Must be a numeric vector of the same length of \code{resp}. Defatul is 50, i.e., equal weights for stability and mean performance.}

  \item{prob}{The p-value for considering a IPCA significant.}

  \item{naxis}{The number of IPCAs to be used for computing the WAAS index. Default is \code{NULL} (Significant IPCAs are used). If values are informed, the number of IPCAS will be used independently on its significance. Note that if two or more variables are included in \code{resp}, then \code{naxis} must be a vector.}


  \item{verbose}{Logical argument. If \code{verbose = FALSE} the code are run silently.}

}

\details{
This function compute the weighted average of absolute scores, estimated as follows:

\deqn{ WAAS_i  = \sum_{k = 1}^{p} |IPCA_{ik} \times EP_k|/ \sum_{k = 1}^{p}EP_k}

where \eqn{WAAS_i} is the weighted average of absolute scores of the \emph{i}th genotype; \eqn{PCA_{ik}} is the score of the \emph{i}th genotype in the \emph{k}th IPCA; and \eqn{EP_k} is the explained variance of the *k*th IPCA for \emph{k = 1,2,..,p}, considering \emph{p} the number of significant PCAs, or a declared number of PCAs. For example if \code{prob = 0.05}, all axis that are significant considering this probability level are used. The number of axis can be also informed by declaring \code{naxis = x}. This comand ignores the \code{p.valuePC} comand.
}

\value{

\item{individual}{A within-environments ANOVA considering a fixed-effect model.}

\item{model}{A data frame with the response variable, the scores of all Principal Components, the estimates of Weighted Average of Absolute Scores, and WAASY (the index that consider the weights for stability and productivity in the genotype ranking.}

\item{MeansGxE}{The means of genotypes in the environments, with observed, predicted and residual values.}

\item{PCA}{Principal Component Analysis.}


\item{anova}{Joint analysis of variance for the main effects and Principal Component analysis of the interaction effect.}

\item{Details}{A list summarizing the results. The following information are showed. \code{WgtResponse}, the weight for the response variable in estimating WAASB, \code{WgtWAAS} the weight for stability, \code{Ngen} the number of genotypes, \code{Nenv} the number of environments, \code{OVmean} the overall mean, \code{Min} the minimum observed (returning the genotype and environment), \code{Max} the maximum observed, \code{Max} the maximum observed, \code{MinENV} the environment with the lower mean, \code{MaxENV} the environment with the larger mean observed, \code{MinGEN} the genotype with the lower mean,  \code{MaxGEN} the genotype with the larger.}

}

\seealso{\code{\link{WAAS.AMMI}}}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\examples{
library(METAAB)

# Considering p-value <= 0.05 to compute the WAAS

model <- WAAS.AMMI(data_ge,
                  resp = GY,
                  gen = GEN,
                  env = ENV,
                  rep = REP)


# Declaring the number of axis to be used for computing WAAS
# and assigning a larger weight for the response variable when
# computing the WAASBY index.

model2 <- WAAS.AMMI(data_ge,
                   resp = GY,
                   gen = GEN,
                   env = ENV,
                   rep = REP,
                   naxis = 3,
                   wresp = 60)

# Analyzing multiple variables (GY and HM) at the same time
# considering that smaller values of HM are better and higher
# values of GY are better, assigning a larger weight for the GY
# and a smaller weight for HM when computing WAASBY index.

model3 <- WAAS.AMMI(data_ge,
                  resp = c(GY, HM),
                  gen = GEN,
                  env = ENV,
                  rep = REP,
                  mresp = c(100, 0),
                  wresp = c(60, 40))

}


