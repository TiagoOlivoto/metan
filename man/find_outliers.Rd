\name{find_outliers}
\alias{find_outliers}
\title{Find possible outliers in a dataset }
\usage{
find_outliers(.data, var, plots = FALSE,
              coef = 1.5, verbose = TRUE)
}
\description{
Find possible outliers in the dataset.

}

\arguments{
 \item{.data}{The data to be analyzed. Must be a dataframe or an object of class \code{group_factors}.}

 \item{var}{The variable to be analyzed..}


  \item{plots}{If \code{TRUE}, then histograms and boxplots are shown.}

 \item{coef}{The multiplication coefficient. For more details see \code{?boxplot.stat}.}

\item{verbose}{If \code{verbose = TRUE} then some results are shown in the console.}

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\examples{
library(METAAB)
library(dplyr)

find_outliers(data_ge2, var = PH, plots = TRUE)

data_ge2 \%>\%
group_factors(ENV) \%>\%
find_outliers(var = PH)

}