% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colindiag.R
\name{colindiag}
\alias{colindiag}
\title{Collinearity Diagnostics}
\usage{
colindiag(.data, ..., by = NULL, n = NULL)
}
\arguments{
\item{.data}{The data to be analyzed. It must be a symmetric correlation
matrix, or a data frame, possible with grouped data passed from
\code{\link[dplyr:group_by]{dplyr::group_by()}}.}

\item{...}{Variables to use in the correlation. If \code{...} is null then
all the numeric variables from \code{.data} are used. It must be a single
variable name or a comma-separated list of unquoted variables names.}

\item{by}{One variable (factor) to compute the function by. It is a shortcut
to \code{\link[dplyr:group_by]{dplyr::group_by()}}. To compute the statistics by more than
one grouping variable use that function.}

\item{n}{If a correlation matrix is provided, then \code{n} is the number of
objects used to compute the correlation coefficients.}
}
\value{
If \code{.data} is a grouped data passed from \code{\link[dplyr:group_by]{dplyr::group_by()}}
then the results will be returned into a list-column of data frames.
\itemize{
\item \strong{cormat} A symmetric Pearson's coefficient correlation matrix
between the variables
\item \strong{corlist} A hypothesis testing for each of the correlation
coefficients
\item \strong{evalevet} The eigenvalues with associated eigenvectors of the
correlation matrix
\item \strong{indicators} A \code{data.frame} with the following indicators
}
\itemize{
\item \code{VIF} The Variance Inflation Factors, being the diagonal elements of
the inverse of the correlation matrix.
\item \code{cn} The Condition Number of the correlation matrix, given by the
ratio between the largest and smallest eigenvalue.
\item \code{det} The determinant of the correlation matrix.
\item \code{ncorhigh} Number of correlation greather than |0.8|.
\item \code{largest_corr} The largest correlation (in absolute value) observed.
\item \code{smallest_corr} The smallest correlation (in absolute value)
observed.
\item \code{weight_var} The variables with largest eigenvector (largest weight)
in the eigenvalue of smallest value, sorted in decreasing order.
}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

Perform a (multi)collinearity diagnostic of a correlation matrix of predictor
variables using several indicators, as shown by Olivoto et al. (2017).
}
\examples{
\donttest{
# Using the correlation matrix
library(metan)

cor_iris <- cor(iris[,1:4])
n <- nrow(iris)

col_diag <- colindiag(cor_iris, n = n)


# Using a data frame
col_diag_gen <- data_ge2 \%>\%
                group_by(GEN) \%>\%
                colindiag()

# Diagnostic by levels of a factor
# For variables with "N" in variable name
col_diag_gen <- data_ge2 \%>\%
                group_by(GEN) \%>\%
                colindiag(contains("N"))
}
}
\references{
Olivoto, T., V.Q. Souza, M. Nardino, I.R. Carvalho, M. Ferrari, A.J.
Pelegrin, V.J. Szareski, and D. Schmidt. 2017. Multicollinearity in path
analysis: a simple method to reduce its effects. Agron. J. 109:131-142.
\doi{10.2134/agronj2016.04.0196}
}
\author{
Tiago Olivoto \email{tiagoolivoto@gmail.com}
}
