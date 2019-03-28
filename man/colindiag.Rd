\name{colindiag}
\alias{colindiag}
\title{Collinearity Diagnostics}
\usage{
colindiag(x, group_var = NULL, n = NULL, verbose = TRUE)
}
\description{
Perform a (multi)collinearity diagnostic of a correlation matrix of predictor variables based several indicators, as shown by Olivoto et al. (2017).

}

\arguments{
 \item{x}{The data to be analyzed. Must be a symmetric correlation matrix or a dataframe containing the predictor variables and one possible grouping variable.}

 \item{group_var}{A possible factor to group the data. If more than two factors are in the dataset, the extra factor column(s) will be excluded with a message.}

  \item{n}{If a correlation matrix is provided, then \code{n} is the number of objects used to compute the correlation coefficients.}


   \item{verbose}{If \code{verbose = TRUE} then some results are shown in the console.}
 }

\value{

The following values are returned. Please, note that if a grouping variable is used, then the results are returned into a list.

\item{cormat}{A symmetric Pearson's coefficient correlation matrix between the variables}

\item{corlist}{A hypothesis testing for each of the correlation coefficients}

\item{evalevet}{The eigenvalues with associated eigenvectors of the correlation matrix}

\item{VIF}{The Variance Inflation Factors, being the diagonal elements of the inverse of the correlation matrix.}

\item{CN}{The Condition Number of the correlation matrix, given by the ratio between the largest and smallest eigenvalue.}

\item{det}{The determinant of the correlation matrix.}

\item{largest_corr}{The largest correlation (in absolute value) observed.}

\item{smallest_corr}{The smallest correlation (in absolute value) observed.}

\item{weight_var}{The variables with largest eigenvector (largest weight) in the eigenvalue of smallest value, sorted in decreasing order.}
}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\references{
Olivoto, T., V.Q. Souza, M. Nardino, I.R. Carvalho, M. Ferrari, A.J. Pelegrin, V.J. Szareski, and D. Schmidt. 2017. Multicollinearity in path analysis: a simple method to reduce its effects. Agron. J. 109:131-142. doi:10.2134/agronj2016.04.0196. \href{https://dl.sciencesocieties.org/publications/aj/abstracts/109/1/131}{10.2134/agronj2016.04.0196}.

}

\examples{
# Using the correlation matrix
cor_iris = cor(iris[,1:4])
n = nrow(iris)
colindiag(cor_iris, n = n)

# Using the original data, storing into an object
col_diag = colindiag(iris, verbose = FALSE)

# Diagnostic by species, storing into an object
col_diag_spec = colindiag(iris, group_var = Species)

}

