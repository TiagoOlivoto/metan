\name{path_coeff}
\alias{path_coeff}
\title{Path coefficients with minimal multicollinearity}
\usage{
path_coeff(.data, resp, pred = NULL, exclude = FALSE,
           correction = NULL, knumber = 50, brutstepwise = FALSE,
           maxvif = 10, missingval = "pairwise.complete.obs")
}
\description{
Estimates of direct and indirect effects. An algorithm to select a set of predictors with minimal multicollinearity and high explanatory power is implemented.

}

\arguments{
 \item{.data}{The data. Must be a dataframe or an object of class \code{group_factors}.}

 \item{resp}{The dependent variable.}

  \item{pred}{The predictor variables, set to \code{NULL}, i.e., the predictor variables are all the numeric variables in the data except that in \code{resp}.}

  \item{exclude}{Logical argument, set to false. If \code{exclude = TRUE}, then the variables in \code{pred} are deleted from the data, and the analysis will use as predictor those that remained, except that in \code{resp}.}

  \item{correction}{Set to \code{NULL}. A correction value (k) that will be added into the diagonal elements of the X'X matrix aiming at reducing the multicollinearity problems (Olivoto et al., 2017). Also as known as ridge regression.}

  \item{knumber}{When \code{correction = NULL}, a plot showing the values of direct effects in a set of different k values (0-1) is produced. \code{knumber} is the number of k values used in the range of 0 to 1.}

  \item{brutstepwise}{Logical argument, set to \code{FALSE}. If true, then an algorithm will select a subset of variables with minimal multicollinearity and fit a set of possible models. See the \bold{Details} section for more information.}

  \item{maxvif}{The maximum value for the Variance Inflaction Factor (cut point) that will be accepted. See the \bold{Details} section for more information.}


  \item{missingval}{How to deal with missing values. For more information, please see \code{?cor}.}

}

\value{
\item{Corr.x}{A correlation matrix between the predictor variables.}

\item{Corr.y}{A vector of correlations between each predictor variable with the dependent variable.}

\item{Coefficients}{The path coefficients. Direct effects are the diagonal elements, and the indirect effects those in the off-diagonal elements (column)}

\item{Eigen}{Eigenvectors and eigenvalues of the \code{Corr.x.}}

\item{VIF}{The Variance Inflaction Factors.}

\item{plot}{A ggplot2-based graphic showing the direct effects in 21 different k values..}

\item{Predictors}{The predictor variables used in the model.}

\item{CN}{The Condition Number, i.e., the ratio between the highest and lowest eigenvalue.}

\item{Det}{The matrix determinant of the \code{Corr.x.}.}

\item{R2}{The coefficient of determination of the model.}

\item{Residual}{The residual effect of the model.}

\item{Response}{The response variable.}

\item{weightvar}{The order of the predictor variables with the higest weigth (highest eigenvector) in the lowest eigenvalue.}

}

\details{
When \code{brutstepwise = TRUE}, first, the algorithm will select a set of predictors with minimal multicollinearity. The selection is based on the variance inflation factor (VIF). An iterative process is performed until the maximum VIF observed is less than \code{maxvif}. The variables selected in this iterative process are then used in a series of stepwise-based regressions. The first model is fitted and p-1 predictor variables are retained (p is the number of variables selected in the iterative process. The second model adjusts a regression considering p-2 selected variables, and so on until the last model, which considers only two variables. Three objects are created. \code{Summary}, with the process summary, \code{Models}, containing the aforementioned values for all the adjusted models; and \code{Selectedpred}, a vector with the name of the selected variables in the iterative process.
}

\references{
Olivoto, T., V.Q. Souza, M. Nardino, I.R. Carvalho, M. Ferrari, A.J. Pelegrin, V.J. Szareski, and D. Schmidt. 2017. Multicollinearity in path analysis: a simple method to reduce its effects. Agron. J. 109:131-142. doi:10.2134/agronj2016.04.0196. \href{https://dl.sciencesocieties.org/publications/aj/abstracts/109/1/131}{10.2134/agronj2016.04.0196}.

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\examples{
# Using the correlation matrix
library(METAAB)
library(dplyr)

pcoeff = data_ge2 \%>\%
         path_coeff(resp = KW, pred = c("TKW", "NKE"))

pcoeff2 = data_ge2 \%>\%
         path_coeff(resp = KW,
                    pred = c("PH", "EH", "NKE", "TKW"),
                    exclude = TRUE)

pcoeff3 = data_ge2 \%>\%
          path_coeff(resp = KW,
                     brutstepwise = TRUE,
                     maxvif = 5)

pcoeff4 = data_ge2 \%>\%
          group_factors(ENV) \%>\%
          path_coeff(resp = KW)

}


