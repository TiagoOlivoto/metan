\name{ge_stats}
\alias{ge_stats}
\title{Parametric statistics for genotype-vs-environment interaction}
\usage{
ge_stats(.data, env, gen, rep, resp, verbose = TRUE)
}
\description{
Compute statistics for genotype vs environment interaction, including, but not limited to within-environment analysis of variance, Ecovalence (Wricke, 1965)  regression-based stability (Eberhart and Eberhart, 1966).
}

\arguments{
 \item{.data}{The dataset containing the columns related to Environments, Genotypes, replication/block and response variable(s).}

  \item{env}{The name of the column that contains the levels of the environments.}

  \item{gen}{The name of the column that contains the levels of the genotypes.}

  \item{rep}{The name of the column that contains the levels of the replications/blocks.}

  \item{resp}{The response variable(s). To analyze multiple variables in a single procedure use, for example, \code{resp = c(var1, var2, var3)}.}

  \item{verbose}{Logical argument. If \code{verbose = FALSE} the code will run silently.}

}

\references{
Eberhart, S.A., and W.A. Russell. 1966. Stability parameters for comparing Varieties. Crop Sci. 6:36-40. \href{https://www.crops.org/publications/cs/abstracts/6/1/CS0060010036}{doi:10.2135/cropsci1966.0011183X000600010011x}.

Wricke, G. 1965. Zur berechnung der okovalenz bei sommerweizen und hafer. Z. Pflanzenzuchtg 52:127-138.

}


\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\examples{
library(METAAB)
library(dplyr)

model = ge_stats(.data = data_ge2,
                 env = ENV,
                 gen = GEN,
                 rep = REP,
                 resp = PH)

# Alternatively, using the forward-pipe operator \%>\%
# More than one trait

model2 = data_ge2 \%>\%
         ge_stats(ENV, GEN, REP, c(PH, EH, EL))

}


