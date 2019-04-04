\name{group_factors}
\alias{group_factors}
\title{Easily splitting a dataframe into subsets grouped by one or more factors}
\usage{
group_factors(.data, ..., keep_factors = FALSE)
}
\description{
Find the required (sufficient) sample size for computing a Pearson correlation coefficient with a desired confidence interval (Olivoto et al., 2018).

}

\arguments{
 \item{.data}{The data that will be split. Must contain at least one grouping variable..}

 \item{...}{The name of grouping variable(s). If more than one grouping variable is used, then they must be unquoted and separated by a comma.}

 \item{keep_factors}{If more than two factors are in the dataframe, should the columns of the non-grouping factors be kept?.}
 }

\value{
\item{A list where each element is a named level of the grouping factors.}
}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}



\examples{

library(METAAB)
library(dplyr)

g1 = iris \%>\% group_factors(Species)

g2 = CO2 \%>\% group_factors(Treatment, Type, keep_factors = TRUE)

}
