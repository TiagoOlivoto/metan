\name{can_corr}
\alias{can_corr}
\title{Canonical correlation analysis}
\usage{
can_corr(.data = NULL, FG = NULL, SG = NULL, use = "cor",
        test = "Bartlett", prob = 0.05, center = TRUE,
        verbose = TRUE, collinearity = TRUE)
}
\description{
Performs canonical correlation analysis with collinearity diagnostic, estimation of canonical loads, canonical scores, and hipotesis testing for correlation pairs.

}

\arguments{
 \item{.data}{The data to be analyzed. Must be a dataframe containing the numeric variables that will be used in the estimation of the correlations. The data can also be passed directly  by the arguments \code{FG} and \code{SG}. Alternatively, \code{.data} may be passed from the function \code{group_factors}. In such case, the canonical correlation will be estimated for each level of the grouping variable in that function.}

 \item{FG}{If a dataframe is informed in \code{.data}, then \code{FG} is a comma-separated list of unquoted variable names that will compose the first (smallest) group of the correlation analysis. FG can also be a cordinate for the variables; for example, \code{FG = data[,1:3]}.}

 \item{SG}{Similar than \code{FG} but for the second group of variables.}


\item{use}{The matrix to be used. Must be one of 'cor' for analysis using the correlation matrix (default) or 'cov' for analysis using the covariance matrix.}


 \item{test}{The test of significance of the relationship between the FG and SG. Must be one of the "Bartlett" (default) or "Rao".}

 \item{prob}{The probability of error assumed. Set to 0.05.}

\item{center}{Should the data be centered to compute the scores?.}


\item{verbose}{Logical argument. If \code{TRUE} (default) then the results are shown in the console.}

\item{collinearity}{Logical argument. If \code{TRUE} (default) then a collinearity diagnostic is performed for each group of variables according to Olivoto et al.(2017).}
}


\value{



\item{Matrix}{The correlation (or covariance) matrix of the variables}

\item{MFG, MSG}{The correlation (or covariance) matrix for the variables of the first group or second group, respectively.}

\item{MFG_SG}{The correlation (or covariance) matrix for the variables of the first group with the second group.}

\item{Coef_FG, Coef_SG}{Matrix of the canonical coefficients of the first group or second group, respectively.}

\item{Loads_FG, Loads_SG}{Matrix of the canonical loadings of the first group or second group, respectively.}

\item{SigTest}{A dataframe with the correlation of the canonical pairs and hypothesis testing results.}

\item{collinearity}{A list with the collinearity diagnostic for each group of variables.}
}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\references{
Olivoto, T., V.Q. Souza, M. Nardino, I.R. Carvalho, M. Ferrari, A.J. Pelegrin, V.J. Szareski, and D. Schmidt. 2017. Multicollinearity in path analysis: a simple method to reduce its effects. Agron. J. 109:131-142. doi:10.2134/agronj2016.04.0196. \href{https://dl.sciencesocieties.org/publications/aj/abstracts/109/1/131}{10.2134/agronj2016.04.0196}.

}

\examples{
\dontrun{
library(METAAB)
library(dplyr)

cc1 = can_corr(data_ge2,
               FG = c(PH, EH, EP),
               SG = c(EL, ED, CL, CD, CW, KW, NR))

cc2 = can_corr(FG = datage_2[, 4:6],
               SG = datage_2[, 7:13],
               verbose = FALSE,
               collinearity = FALSE)

# Canonical correlations for each environment
cc3 = data_ge2 \%>\%
      group_factors(ENV, REP) \%>\%
      can_corr(FG = c(PH, EH, EP),
               SG = c(EL, ED, CL, CD, CW, KW, NR))

}
}

