\name{anova_ind}
\alias{anova_ind}
\title{Within-environment analysis of variance}
\usage{
anova_ind(.data, env, gen, rep, resp)
}
\description{
This is a helper function that performs a within-environment analysis of variance and returns values such as Mean Squares, p-values, coefficient of variation, heritability, and accuracy of selection.

}

\arguments{
  \item{.data}{The dataset containing the columns related to Environments, Genotypes, replication/block and response variable(s).}

  \item{env}{The name of the column that contains the levels of the environments. The analysis of variance is computed for each level of this factor.}

  \item{gen}{The name of the column that contains the levels of the genotypes.}

  \item{rep}{The name of the column that contains the levels of the replications/blocks.}

  \item{resp}{The response variable.}


}


\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\examples{
library(METAAB)
# traditional usage approach
data = data_ge
anova1 = anova_ind(data,
                   env = ENV,
                   gen = GEN,
                   rep = REP,
                   resp = GY)

# Using the pipe operator \%>\%
library(dplyr)
data \%>\% anova_ind(ENV, GEN, REP, GY)

}


