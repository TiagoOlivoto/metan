\name{comb_vars}
\alias{comb_vars}
\title{Generate pairwise combinations of variables by applying one function to each pair}
\usage{
comb_vars(.data, order = "second", FUN = "+", verbose = TRUE)
}
\description{
Pairwise combinations of variables that will be the result of a function applied to each combination.

}

\arguments{
 \item{.data}{A matrix of data with, say, p columns.}

  \item{order}{The order on how the results will appear in the output. Default is \code{order = "second"}. In this case, assuming that .data has four columns, namely, \code{V1, V2, V3, V4}, the order of columns in the output will be \code{V1.V2, V1.V3, V2.V3, V1.V4, V2.V4, V3.V4}. If 'order' is set to "first", the result will be then \code{V1.V2, V1.V3, V1.V4, V2.V3, V2.V4, V3.V4}.}

  \item{FUN}{The function that will be applied to each combination. The default is \code{+}, i.e., V1 + V2.}

 \item{verbose}{Logical argument. If \code{verbose = FALSE} the code will run silently.}

}


\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\examples{
\dontrun{
library(METAAB)
data = data.frame(A = rnorm(n = 5, mean = 10, sd = 3),
                  B = rnorm(n = 5, mean = 120, sd = 30),
                  C = rnorm(n = 5, mean = 40, sd = 10),
                  D = rnorm(n = 5, mean = 1100, sd = 200),
                  E = rnorm(n = 5, mean = 2, sd = 1))
comb_vars(data)
comb_vars(data, FUN = "*", order = "first")
}
}



