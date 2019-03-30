\name{find_outliers}
\alias{find_outliers}
\title{Find outliers in a dataset }
\usage{
find_outliers(.data, var, group_var = NULL, plots = FALSE, coef = 1.5)
}
\description{
Find possible outliers in the dataset.

}

\arguments{
 \item{.data}{The data to be analyzed. Must be a symmetric correlation matrix or a dataframe containing the predictor variables and one possible grouping variable.}

 \item{var}{The variable to be analyzed..}

 \item{group_var}{A possible factor to group the data. If more than two factors are in the dataset, the extra factor column(s) will be excluded with a message.}

  \item{plots}{If \code{TRUE}, then histograms and boxplots are shown.}

 \item{coef}{The multiplication coefficient. For more details see \code{?boxplot.stat}.}
 }

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

