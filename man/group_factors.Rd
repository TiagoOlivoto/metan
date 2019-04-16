\name{group_factors}
\alias{group_factors}
\title{Split a dataframe into subsets grouped by one or more factors}
\usage{
group_factors(.data, ..., keep_factors = FALSE, verbose = TRUE)
}
\description{
Split a dataframe into subsets grouped by one or more factors

}

\arguments{
 \item{.data}{The data that will be split. Must contain at least one grouping variable..}

 \item{...}{Comma-separated list of unquoted variable names that will be used to group the data.}

 \item{keep_factors}{If more than two factors are in the dataframe, should the columns of the non-grouping factors be kept?.}

 \item{verbose}{Logical argument. If \code{verbose = FALSE} the code will run silently.}

}

\value{
A list where each element is a named level of the grouping factors.
}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}



\examples{

library(METAAB)
library(dplyr)

g1 = iris \%>\% group_factors(Species)

g2 = CO2 \%>\% group_factors(Treatment, Type, keep_factors = TRUE)

}
